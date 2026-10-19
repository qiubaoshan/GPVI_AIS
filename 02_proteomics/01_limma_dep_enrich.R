#!/usr/bin/env Rscript
# =============================================================================
# 01_limma_dep_enrich.R  --  Quantitative proteomics: normalisation, differential
#                            abundance (limma), rank-abundance, heatmap and
#                            GO / KEGG enrichment.
#
# Used for both proteomic matrices in the paper:
#   PROT-PLT     isolated platelets, 3 h after reperfusion   (Fig. 1I-1J, S1D-S1E)
#   PROT-PLASMA  plasma, 3 h after reperfusion               (Fig. S5H-S5I)
#
# Process the two matrices separately and never in one normalisation step: plasma
# spans roughly ten orders of magnitude of dynamic range and would dominate any
# shared scaling.
#
# GROUP NAMING. Genotypes are WT (control) and KO (Gp6 knockout). The exported
# matrix columns carry the legacy submission codes NS (WT) and YY (KO); the
# mapping is declared in the sample sheet (`sample_id`, `genotype`) and never
# inferred from column order. Figures print `--control-label` / `--case-label`.
#
# Input: protein-level intensity matrix exported from Proteome Discoverer 2.4.1.15
#        (1% FDR at PSM, peptide and protein level), tab-separated:
#
#          Protein            Gene    NS1    NS2    NS3    YY1    YY2    YY3
#          P12345             Flna    1.2e8  ...
#
# Usage
#   Rscript 01_limma_dep_enrich.R --matrix platelet_protein_intensity.tsv \
#                              --samples samples/PROT-PLT.tsv \
#                              --outdir results/PROT-PLT --tag PROT-PLT
#
# See README.md section 6 for the acquisition parameters and the threshold caveat.
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

# -----------------------------------------------------------------------------
# 0. Arguments
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if ("--help" %in% args || "-h" %in% args) {
  cat(
    "Usage: Rscript 01_limma_dep_enrich.R [options]\n\n",
    "  --matrix    FILE  protein-level intensity matrix (Protein, Gene, samples)\n",
    "  --samples   FILE  optional sample sheet: sample_id, genotype\n",
    "  --outdir    DIR   output directory\n",
    "  --tag       STR   prefix for output files\n",
    "  --control   STR   control level in `genotype`                [WT]\n",
    "  --case      STR   case level in `genotype`                   [KO]\n",
    "  --control-label STR label printed in figures                 [WT]\n",
    "  --case-label    STR label printed in figures                 [Gp6-KO]\n",
    "  --min-frac  NUM   min detection fraction within a group     [0.7]\n",
    "  --lfc       NUM   |log2FC| threshold                        [0.58]\n",
    "  --p         NUM   significance threshold                    [0.05]\n",
    "  --pcol      STR   adj.P.Val or P.Value                      [adj.P.Val]\n",
    "  --top-n     INT   proteins shown in the heatmap             [40]\n",
    "  --already-log2    intensities are already log2-transformed\n",
    sep = ""
  )
  quit(save = "no", status = 0)
}

get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

matrix_file  <- get_arg("--matrix", "protein_intensity.tsv")
samples_file <- get_arg("--samples", NA)
outdir       <- get_arg("--outdir", "results")
tag          <- get_arg("--tag", "PROT-PLT")
ctrl_level   <- get_arg("--control", "WT")
case_level   <- get_arg("--case", "KO")
ctrl_name    <- get_arg("--control-label", "WT")
case_name    <- get_arg("--case-label", "Gp6-KO")
min_frac     <- as.numeric(get_arg("--min-frac", 0.7))
lfc_cut      <- as.numeric(get_arg("--lfc", 0.58))
p_cut        <- as.numeric(get_arg("--p", 0.05))
p_col        <- get_arg("--pcol", "adj.P.Val")
top_n        <- as.integer(get_arg("--top-n", 40))
already_log2 <- "--already-log2" %in% args

if (!p_col %in% c("adj.P.Val", "P.Value")) {
  stop("--pcol must be adj.P.Val or P.Value")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
figdir <- file.path(outdir, "figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
message("[proteomics] tag = ", tag, " | contrast ", case_name, " vs ", ctrl_name,
        " | p column = ", p_col)

save_pdf <- function(p, name, width = 6, height = 5) {
  ggsave(file.path(figdir, paste0(tag, "_", name, ".pdf")), p, width = width, height = height)
}

# -----------------------------------------------------------------------------
# 1. Load
# -----------------------------------------------------------------------------
prot <- read.delim(matrix_file, check.names = FALSE, stringsAsFactors = FALSE)

gene_col <- if ("Gene" %in% colnames(prot)) "Gene" else
            if ("Gene.names" %in% colnames(prot)) "Gene.names" else NA
if (is.na(gene_col)) stop("No Gene / Gene.names column in the matrix.")

if (!is.na(samples_file)) {
  meta <- read.delim(samples_file, stringsAsFactors = FALSE)
  if (!"genotype" %in% colnames(meta) && "group" %in% colnames(meta)) {
    meta$genotype <- meta$group
  }
  stopifnot(all(c("sample_id", "genotype") %in% colnames(meta)))
  absent <- setdiff(meta$sample_id, colnames(prot))
  if (length(absent) > 0) {
    stop("Sample(s) in the sheet not found in the matrix: ",
         paste(absent, collapse = ", "))
  }
  sample_cols <- meta$sample_id
  group <- factor(meta$genotype, levels = c(ctrl_level, case_level))
  if (any(is.na(group))) {
    stop("`genotype` must be ", ctrl_level, " or ", case_level)
  }
} else {
  sample_cols <- grep(paste0("^(", ctrl_level, "|", case_level, ")"),
                      colnames(prot), value = TRUE)
  if (length(sample_cols) == 0) {
    stop("No sample columns found. Name them ", ctrl_level, "1, ", ctrl_level,
         "2, ... and ", case_level, "1, ... or supply --samples.")
  }
  group <- factor(ifelse(grepl(paste0("^", case_level), sample_cols),
                         case_level, ctrl_level),
                  levels = c(ctrl_level, case_level))
}
group_label <- factor(ifelse(group == ctrl_level, ctrl_name, case_name),
                      levels = c(ctrl_name, case_name))
message("[proteomics] ", length(sample_cols), " samples: ",
        paste(table(group), collapse = " vs "))

# -----------------------------------------------------------------------------
# 2. Filter and normalise
# -----------------------------------------------------------------------------
mat  <- as.matrix(prot[, sample_cols, drop = FALSE])
gene <- prot[[gene_col]]

keep <- rowSums(mat > 0) >= min_frac * min(table(group))
message("[proteomics] proteins retained: ", sum(keep), " / ", nrow(mat))
mat <- mat[keep, , drop = FALSE]
gene <- gene[keep]

if (!already_log2) mat <- log2(mat)
mat[!is.finite(mat)] <- NA

# Quantile normalisation within matrix. Platelet and plasma must be normalised
# in separate calls -- never pool the two matrices.
lmat <- normalizeBetweenArrays(mat, method = "quantile")

# -----------------------------------------------------------------------------
# 3. limma moderated t
# -----------------------------------------------------------------------------
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit  <- lmFit(lmat, design)
cont <- makeContrasts(contrasts = paste0(case_level, "-", ctrl_level), levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont))

tt <- topTable(fit2, number = Inf, sort.by = "P")
tt$Protein <- rownames(tt)
tt$Gene    <- gene[match(rownames(tt), rownames(lmat))]
tt <- tt[, c("Protein", "Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

write.csv(tt, file.path(outdir, paste0(tag, "_limma_all_proteins.csv")), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Volcano
# -----------------------------------------------------------------------------
pval <- tt[[p_col]]
tt$change <- ifelse(abs(tt$logFC) > lfc_cut & pval < p_cut,
                    ifelse(tt$logFC > 0, paste0("Higher in ", case_name),
                           paste0("Lower in ", case_name)), "Not significant")
# Proteins with an NA test statistic cannot be called; keep them out of the counts.
tt$change[is.na(tt$change)] <- "Not significant"
tt$change <- factor(tt$change,
                    levels = c(paste0("Lower in ", case_name), "Not significant",
                               paste0("Higher in ", case_name)))
n_dep <- sum(tt$change != "Not significant")
message("[proteomics] differential proteins (", p_col, "): ", n_dep,
        " (", sum(tt$change == paste0("Lower in ", case_name)), " down, ",
        sum(tt$change == paste0("Higher in ", case_name)), " up)")

write.table(
  data.frame(dataset = tag, contrast = paste0(case_name, " vs ", ctrl_name),
             n_detected = nrow(tt), n_dep = n_dep,
             n_higher_in_case = sum(tt$change == paste0("Higher in ", case_name)),
             n_lower_in_case  = sum(tt$change == paste0("Lower in ", case_name)),
             lfc = lfc_cut, p_column = p_col, p = p_cut),
  file.path(outdir, paste0(tag, "_DEP_summary.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

p_volcano <- ggplot(tt, aes(logFC, -log10(pval), colour = change)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_hline(yintercept = -log10(p_cut), linetype = 5) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 5) +
  scale_colour_manual(values = c("navy", "grey", "firebrick3")) +
  labs(x = paste0("log2 fold change (", case_name, " / ", ctrl_name, ")"),
       y = paste0("-log10 ", p_col)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), legend.title = element_blank())
save_pdf(p_volcano, "volcano", width = 5.5, height = 5)

# -----------------------------------------------------------------------------
# 5. Rank-abundance
# -----------------------------------------------------------------------------
# Methods: "normalized protein abundances were log10-transformed and ranked
# according to mean abundance within each group". lmat is already log2, so the
# conversion is a constant factor; ranking is unaffected either way.
ra <- data.frame(
  Protein    = rownames(lmat),
  Gene       = gene,
  mean_log10_ctrl = rowMeans(lmat[, group == ctrl_level, drop = FALSE], na.rm = TRUE) * log10(2),
  mean_log10_case = rowMeans(lmat[, group == case_level, drop = FALSE], na.rm = TRUE) * log10(2)
)
ra <- ra[order(-ra$mean_log10_ctrl), ]
ra$rank <- seq_len(nrow(ra))
write.csv(ra, file.path(outdir, paste0(tag, "_rank_abundance.csv")), row.names = FALSE)

ra_long <- rbind(
  data.frame(rank = ra$rank, abundance = ra$mean_log10_ctrl, genotype = ctrl_name),
  data.frame(rank = ra$rank, abundance = ra$mean_log10_case, genotype = case_name)
)
ra_long$genotype <- factor(ra_long$genotype, levels = c(ctrl_name, case_name))

p_rank <- ggplot(ra_long, aes(rank, abundance, colour = genotype)) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = setNames(c("#4d4d4d", "#c0392b"),
                                        c(ctrl_name, case_name))) +
  labs(x = "Protein rank", y = "Mean abundance (log10)", colour = "Genotype") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), legend.title = element_blank())
save_pdf(p_rank, "rank_abundance", width = 5.5, height = 4.5)

# -----------------------------------------------------------------------------
# 6. Heatmap of selected proteins (row Z-score)
# -----------------------------------------------------------------------------
selected <- head(tt$Protein[order(tt[[p_col]])], top_n)
selected <- unique(selected)
if (length(selected) >= 2) {
  hm <- lmat[selected, , drop = FALSE]
  rownames(hm) <- make.unique(ifelse(is.na(tt$Gene[match(rownames(hm), tt$Protein)]),
                                     rownames(hm),
                                     tt$Gene[match(rownames(hm), tt$Protein)]))
  hm[is.na(hm)] <- rowMeans(hm, na.rm = TRUE)[row(hm)[is.na(hm)]]
  ann_col <- data.frame(Genotype = group_label, row.names = sample_cols)
  pheatmap(hm, scale = "row",
           annotation_col = ann_col,
           annotation_colors = list(Genotype = setNames(
             c("#4d4d4d", "#c0392b"), c(ctrl_name, case_name))),
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           cluster_cols = FALSE, show_rownames = TRUE, fontsize_row = 7,
           filename = file.path(figdir, paste0(tag, "_heatmap.pdf")),
           width = 5, height = 8)
}

# -----------------------------------------------------------------------------
# 7. GO / KEGG enrichment with an explicit rich-factor denominator
# -----------------------------------------------------------------------------
dep_genes <- unique(na.omit(tt$Gene[tt$change != "Not significant"]))
bg_genes  <- unique(na.omit(gene))

if (length(dep_genes) >= 5) {
  ego <- enrichGO(dep_genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                  universe = bg_genes,
                  pAdjustMethod = "BH", pvalueCutoff = p_cut, qvalueCutoff = p_cut)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    # Rich factor = DEPs in term / all detected proteins annotated to the term.
    # clusterProfiler's BgRatio uses the whole annotation universe, which is not
    # what the figure legends describe; recompute against detected proteins only.
    bg_by_term <- vapply(ego@geneSets[ego@result$ID], function(g)
      length(intersect(bg_genes, g)), integer(1))
    ego@result$RichFactor <- ego@result$Count / bg_by_term
    write.csv(ego@result, file.path(outdir, paste0(tag, "_GO_BP.csv")), row.names = FALSE)

    # Rich-factor dot plot, drawn explicitly because enrichplot::dotplot does not
    # accept a custom numeric column on the x-axis.
    rf <- ego@result[order(-ego@result$RichFactor), ]
    rf <- head(rf, 20)
    rf$Description <- factor(rf$Description, levels = rev(rf$Description))
    p_rf <- ggplot(rf, aes(RichFactor, Description, size = Count, colour = -log10(p.adjust))) +
      geom_point() +
      scale_colour_gradient(low = "navy", high = "firebrick3") +
      labs(x = "Rich factor", y = NULL) +
      theme_bw(base_size = 12)
    save_pdf(p_rf, "GO_BP_richfactor", width = 8, height = 7)
  } else {
    warning("No enriched GO term.")
  }

  entrez <- tryCatch(
    bitr(dep_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db),
    error = function(e) NULL)
  if (!is.null(entrez) && nrow(entrez) > 0) {
    ekegg <- enrichKEGG(entrez$ENTREZID, organism = "mmu",
                        pAdjustMethod = "BH", pvalueCutoff = p_cut)
    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      write.csv(ekegg@result, file.path(outdir, paste0(tag, "_KEGG.csv")), row.names = FALSE)
      save_pdf(dotplot(ekegg, showCategory = 20), "KEGG", width = 8, height = 7)
    }
  }
} else {
  warning("Fewer than 5 differential proteins -- enrichment skipped.")
}

writeLines(c(paste0("tag\t", tag),
             paste0("contrast\t", case_name, " vs ", ctrl_name),
             paste0("p_column\t", p_col), paste0("lfc\t", lfc_cut),
             paste0("p\t", p_cut), paste0("min_frac\t", min_frac),
             paste0("samples\t", paste(sample_cols, collapse = ",")),
             capture.output(sessionInfo())),
           file.path(outdir, paste0(tag, "_sessionInfo.txt")))

message("[proteomics] done. Outputs in ", normalizePath(outdir))
