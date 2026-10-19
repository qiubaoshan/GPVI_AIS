#!/usr/bin/env Rscript
# =============================================================================
# 01_cellchat_compare.R  --  Cell-cell communication inference and comparative analysis.
#
# Input: the annotated Seurat object exported from 03_scrnaseq
#        (<tag>_01_annotated.h5ad, log-normalised, with cell_type).
#
# Design follows the paper: communication is inferred INDEPENDENTLY in each
# genotype and the two objects are merged afterwards for comparison. Fitting one
# joint model and splitting it is not equivalent and must not be substituted.
#
# GROUP NAMING. `condition` carries WT / KO; the merged objects are named with the
# display labels (WT / Gp6-KO) so every legend matches the manuscript. Output file
# names use the short analysis labels, because those are what the run record lists.
#
# Usage
#   Rscript 01_cellchat_compare.R --h5 all_data_merge_log1p_raw.h5 \
#                      --outdir results/CELLCHAT --tag CELLCHAT \
#                      --senders END,PFB,AST \
#                      --pathways NRXN,CALCR,VISFATIN
#
# See README.md section 8 for parameter rationale.
# =============================================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(patchwork)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# 0. Arguments
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if ("--help" %in% args || "-h" %in% args) {
  cat(
    "Usage: Rscript 01_cellchat_compare.R [options]\n\n",
    "  --h5        FILE  annotated h5 (Seurat-compatible, contains cell_type)\n",
    "  --outdir    DIR   output directory\n",
    "  --tag       STR   prefix for output files\n",
    "  --ctrl      STR   control label in `condition`                [WT]\n",
    "  --case      STR   case label in `condition`                   [KO]\n",
    "  --ctrl-name STR   control display label in figures             [WT]\n",
    "  --case-name STR   case display label in figures                [Gp6-KO]\n",
    "  --senders   STR   comma-separated cell_type suffixes to plot  [END,PFB]\n",
    "  --pathways  STR   comma-separated pathways for heatmaps       [NRXN,CALCR,VISFATIN]\n",
    "  --min-cells INT   filterCommunication min.cells               [10]\n",
    "  --nboot     INT   permutations for communication probability  [100]\n",
    "  --workers   INT   future multisession workers                 [10]\n",
    sep = ""
  )
  quit(save = "no", status = 0)
}

get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
split_arg <- function(x) trimws(strsplit(x, ",")[[1]])

h5_file   <- get_arg("--h5", "all_data_merge_log1p_raw.h5")
outdir    <- get_arg("--outdir", "results/CELLCHAT")
tag       <- get_arg("--tag", "CELLCHAT")
ctrl_lab  <- get_arg("--ctrl", "WT")
case_lab  <- get_arg("--case", "KO")
ctrl_name <- get_arg("--ctrl-name", "WT")
case_name <- get_arg("--case-name", "Gp6-KO")
senders   <- split_arg(get_arg("--senders", "END,PFB"))
pathways  <- split_arg(get_arg("--pathways", "NRXN,CALCR,VISFATIN"))
min_cells <- as.integer(get_arg("--min-cells", 10))
nboot     <- as.integer(get_arg("--nboot", 100))
workers   <- as.integer(get_arg("--workers", 10))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
figdir <- file.path(outdir, "figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

#' Resolve a cell population name to its ident index.
#'
#' Idents carry a numeric prefix ("9:AST", "12:PFB"). The numeric index changes
#' whenever clustering is re-run, so resolve by name, never by hard-coded number.
#' @param x  a CellChat object (single or merged)
#' @param label  suffix ("AST") or full label ("9:AST")
ident_index <- function(x, label) {
  lv <- levels(x@idents)
  hit <- which(lv == label)
  if (length(hit) == 0) hit <- grep(paste0(":", label, "$"), lv)
  if (length(hit) == 0) hit <- grep(label, lv, fixed = TRUE)
  if (length(hit) == 0) {
    stop("Population '", label, "' not found. Available idents: ",
         paste(lv, collapse = ", "))
  }
  hit[1]
}

pdf_dev <- function(file, width = 8, height = 7, expr) {
  pdf(file.path(figdir, file), width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  obj <- force(expr)
  # ggplot / patchwork objects are not drawn automatically on a pdf device.
  if (inherits(obj, c("ggplot", "patchwork", "gg"))) print(obj)
  invisible(obj)
}

# -----------------------------------------------------------------------------
# 1. Load and split
# -----------------------------------------------------------------------------
data.input <- dior::read_h5(h5_file, target.object = "seurat")
stopifnot("cell_type" %in% colnames(data.input@meta.data),
          "condition" %in% colnames(data.input@meta.data))

case_obj <- subset(data.input, cells = colnames(data.input)[data.input$condition == case_lab])
ctrl_obj <- subset(data.input, cells = colnames(data.input)[data.input$condition == ctrl_lab])
message("[cellchat] cells: ", ctrl_lab, " = ", ncol(ctrl_obj),
        ", ", case_lab, " = ", ncol(case_obj))

cc_ctrl <- createCellChat(ctrl_obj, meta = ctrl_obj@meta.data, group.by = "cell_type")
cc_case <- createCellChat(case_obj, meta = case_obj@meta.data, group.by = "cell_type")

CellChatDB <- CellChatDB.mouse
cc_ctrl@DB <- CellChatDB
cc_case@DB <- CellChatDB

cc_ctrl <- subsetData(cc_ctrl)
cc_case <- subsetData(cc_case)

future::plan("multisession", workers = workers)
on.exit(future::plan("sequential"), add = TRUE)

# -----------------------------------------------------------------------------
# 2. Inference, per genotype
# -----------------------------------------------------------------------------
run_cellchat <- function(cc, label) {
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  # raw.use = TRUE: fit on the log-normalised matrix directly, no smoothed data.
  cc <- computeCommunProb(cc, raw.use = TRUE, nboot = nboot)
  cc <- filterCommunication(cc, min.cells = min_cells)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")

  for (t in c("functional", "structural")) {
    cc <- computeNetSimilarity(cc, type = t)
    cc <- netEmbedding(cc, type = t)
    cc <- netClustering(cc, type = t)
  }

  saveRDS(cc, file.path(outdir, paste0(tag, "_cellchat_", label, ".rds")))

  # Export the inference tables: reviewers ask for these, figures are not data.
  comm <- subsetCommunication(cc)
  write.csv(comm, file.path(outdir, paste0(tag, "_communication_", label, ".csv")),
            row.names = FALSE)

  n_lr <- nrow(comm)
  message("[cellchat] ", label, ": ", n_lr, " significant L-R interactions")
  cc
}

cc_ctrl <- run_cellchat(cc_ctrl, ctrl_lab)
cc_case <- run_cellchat(cc_case, case_lab)

# -----------------------------------------------------------------------------
# 3. Merge and compare
# -----------------------------------------------------------------------------
object.list <- list(cc_ctrl, cc_case)
# Dataset names are the display labels: CellChat writes them into every legend and
# into `datasets`, so pos.dataset / datasets below must use the same strings.
names(object.list) <- c(ctrl_name, case_name)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

p1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2))
p2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")
ggsave(file.path(figdir, paste0(tag, "_interaction_counts.pdf")), p1 + p2,
       width = 8, height = 4)

for (t in c("functional", "structural")) {
  cellchat <- computeNetSimilarityPairwise(cellchat, type = t)
  cellchat <- netEmbedding(cellchat, type = t)
  cellchat <- netClustering(cellchat, type = t)
}
pdf_dev(paste0(tag, "_embedding_functional.pdf"), 7, 6,
        netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5))
pdf_dev(paste0(tag, "_embedding_structural.pdf"), 7, 6,
        netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5))

# Ranked differences in outgoing / incoming signalling strength
pdf_dev(paste0(tag, "_ranknet.pdf"), 8, 10, {
  rankNet(cellchat, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE)
})

# -----------------------------------------------------------------------------
# 4. Differential ligand-receptor analysis
# -----------------------------------------------------------------------------
cellchat <- identifyOverExpressedGenes(cellchat,
                                       group.dataset = "datasets",
                                       pos.dataset = case_name,
                                       features.name = case_name,
                                       only.pos = FALSE,
                                       thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 0.05)

net <- netMappingDEG(cellchat, features.name = case_name)

net.up <- subsetCommunication(cellchat, net = net, datasets = case_name,
                              ligand.logFC = 0.1, receptor.logFC = 0.1)
net.down <- subsetCommunication(cellchat, net = net, datasets = ctrl_name,
                                ligand.logFC = -0.1, receptor.logFC = -0.1)

write.csv(net.up,   file.path(outdir, paste0(tag, "_L-R_up_in_", case_lab, ".csv")),
          row.names = FALSE)
write.csv(net.down, file.path(outdir, paste0(tag, "_L-R_up_in_", ctrl_lab, ".csv")),
          row.names = FALSE)
message("[cellchat] differential L-R pairs: higher in ", case_name, " = ", nrow(net.up),
        ", higher in ", ctrl_name, " = ", nrow(net.down))

gene.up   <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

# -----------------------------------------------------------------------------
# 5. Bubble plots for the populations discussed in the text
# -----------------------------------------------------------------------------
targets <- seq_along(levels(cellchat@idents))

for (s in senders) {
  idx <- ident_index(cellchat, s)

  pdf_dev(paste0(tag, "_bubble_all_LR_", gsub("[^A-Za-z0-9]", "", s), ".pdf"), 8, 10,
          netVisual_bubble(cellchat, sources.use = idx, targets.use = targets,
                           comparison = c(1, 2), angle.x = 45,
                           color.text = c("#ABB1A5", "#C14F11")))

  pdf_dev(paste0(tag, "_bubble_up_", gsub("[^A-Za-z0-9]", "", s), ".pdf"), 8, 10,
          netVisual_bubble(cellchat, pairLR.use = net.up[, "interaction_name", drop = FALSE],
                           sources.use = idx, targets.use = targets,
                           comparison = c(1, 2), angle.x = 45, remove.isolate = TRUE,
                           color.text = c("#ABB1A5", "#C14F11")))

  pdf_dev(paste0(tag, "_bubble_down_", gsub("[^A-Za-z0-9]", "", s), ".pdf"), 8, 10,
          netVisual_bubble(cellchat, pairLR.use = net.down[, "interaction_name", drop = FALSE],
                           sources.use = idx, targets.use = targets,
                           comparison = c(1, 2), angle.x = 45, remove.isolate = TRUE,
                           color.text = c("#ABB1A5", "#C14F11")))
}

# Differential interaction heatmaps for those senders
for (s in senders) {
  pdf_dev(paste0(tag, "_diffInteraction_", gsub("[^A-Za-z0-9]", "", s), ".pdf"), 10, 5,
          netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight",
                                    sources.use = ident_index(cellchat, s)))
}

# -----------------------------------------------------------------------------
# 6. Pathway-level heatmaps
# -----------------------------------------------------------------------------
for (pw in pathways) {
  if (!pw %in% cellchat@netP[[1]]$pathways) {
    warning("Pathway ", pw, " was not inferred in ", ctrl_lab, "; skipped.")
    next
  }
  ht <- lapply(object.list, function(x) {
    netVisual_heatmap(x, signaling = pw, color.heatmap = "Reds",
                      title.name = paste(pw, "signaling"))
  })
  pdf(file.path(figdir, paste0(tag, "_pathway_", pw, ".pdf")), width = 12, height = 5)
  ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = grid::unit(0.5, "cm"))
  dev.off()

  pdf(file.path(figdir, paste0(tag, "_genes_", pw, ".pdf")), width = 9, height = 7)
  print(plotGeneExpression(cellchat, signaling = pw, split.by = "datasets",
                           colors.ggplot = TRUE, color.use = c("#ABB1A5", "#C14F11")))
  dev.off()
}

# -----------------------------------------------------------------------------
# 7. Save
# -----------------------------------------------------------------------------
saveRDS(cellchat, file.path(outdir, paste0(tag, "_cellchat_merged.rds")))
writeLines(capture.output(sessionInfo()),
           file.path(outdir, paste0(tag, "_sessionInfo.txt")))
message("[cellchat] done. Outputs in ", normalizePath(outdir))
