#!/usr/bin/env Rscript
# =============================================================================
# 01_deseq2_deg_enrich.R  --  Bulk RNA-seq: DESeq2 differential expression, GO
#                             over-representation analysis and GSEA.
#
# Used for both bulk datasets in the paper:
#   BULK-BRAIN   ischemic brain tissue, 24 h after reperfusion  (Fig. S5A-S5D)
#   BULK-NEU     sorted neutrophils, basal and 3 h              (Fig. 4A-4B, S5E-S5G)
#
# The two datasets differ only in the counts matrix and the sample sheet.
# Nothing in this script is dataset specific.
#
# GROUP NAMING. Genotypes are handled as WT (control) and KO (Gp6 knockout), the
# same labels the manuscript uses. The sequencing libraries carry the legacy
# submission codes NS (WT) and YY (KO); the mapping is declared once, in the
# sample sheet (`sample_id`, `genotype`), and never re-derived from the names.
# Figures print `--control-label` / `--case-label`, so the panels read "WT" and
# "Gp6-KO" without touching the code.
#
# CONTRAST. Positive log2FC = higher in KO. The manuscript describes genes as
# "down-regulated in KO"; those rows carry negative log2FC.
#
# Usage
#   Rscript 01_deseq2_deg_enrich.R \
#     --counts all_samples.gene_counts.txt \
#     --samples samples/BULK-BRAIN.tsv \
#     --annot Mus_musculus.GRCm39.gene.csv \
#     --outdir results/BULK-BRAIN \
#     --tag BULK-BRAIN
#
# Run `Rscript 01_deseq2_deg_enrich.R --help` for all options.
# See README.md section 5 for parameter rationale and the published DEG counts.
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(ggplot2)
  library(FactoMineR)
  library(factoextra)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# 0. Arguments
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if ("--help" %in% args || "-h" %in% args) {
  cat(
    "Usage: Rscript 01_deseq2_deg_enrich.R [options]\n\n",
    "  --counts   FILE   featureCounts matrix (GeneID + one column per library)\n",
    "  --samples  FILE   sample sheet: sample_id, genotype[, timepoint]\n",
    "  --annot    FILE   gene annotation: gene_id, gene_name, gene_biotype\n",
    "  --outdir   DIR    output directory (created if absent)\n",
    "  --tag      STR    prefix for output files\n",
    "  --control  STR    control level in `genotype`               [WT]\n",
    "  --case     STR    case level in `genotype`                  [KO]\n",
    "  --control-label STR label printed in figures                [WT]\n",
    "  --case-label    STR label printed in figures                [Gp6-KO]\n",
    "  --min-count INT   pre-filter: min count in >= smallest group [10]\n",
    "  --lfc      NUM    |log2FC| threshold for DEG calling      [1]\n",
    "  --padj     NUM    adjusted P threshold for DEG calling    [0.05]\n",
    "  --hcut-k   INT    k for hierarchical clustering QC plot   [4]\n",
    "  --gsea-term STR   GO term drawn as a GSEA running score   [GO:0032755]\n",
    "  --gsea-metric STR stat | lfc -- gene ranking metric       [stat]\n",
    sep = ""
  )
  quit(save = "no", status = 0)
}

get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

counts_file  <- get_arg("--counts",  "all_samples.gene_counts.txt")
samples_file <- get_arg("--samples", "samples/BULK-BRAIN.tsv")
annot_file   <- get_arg("--annot",   "Mus_musculus.GRCm39.gene.csv")
outdir       <- get_arg("--outdir",  "results")
tag          <- get_arg("--tag",     "BULK-BRAIN")
ctrl_level   <- get_arg("--control", "WT")
case_level   <- get_arg("--case",    "KO")
ctrl_name    <- get_arg("--control-label", "WT")
case_name    <- get_arg("--case-label",    "Gp6-KO")
min_count    <- as.integer(get_arg("--min-count", 10))
lfc_cut      <- as.numeric(get_arg("--lfc", 1))
padj_cut     <- as.numeric(get_arg("--padj", 0.05))
hcut_k       <- as.integer(get_arg("--hcut-k", 4))
gsea_term    <- get_arg("--gsea-term", "GO:0032755")
gsea_metric  <- get_arg("--gsea-metric", "stat")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
figdir <- file.path(outdir, "figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

save_pdf <- function(plot_obj, name, width = 8, height = 6) {
  ggsave(file.path(figdir, paste0(tag, "_", name, ".pdf")),
         plot_obj, width = width, height = height)
}

message("[01_deseq2_deg_enrich.R] tag = ", tag, " | contrast ", case_name, " vs ", ctrl_name)

# -----------------------------------------------------------------------------
# 1. Counts and sample sheet
# -----------------------------------------------------------------------------
raw <- read.delim(counts_file, check.names = FALSE)
rownames(raw) <- raw[[1]]

meta <- read.delim(samples_file, stringsAsFactors = FALSE)
if (!"genotype" %in% colnames(meta) && "group" %in% colnames(meta)) {
  meta$genotype <- meta$group          # tolerate the older column name
}
stopifnot(all(c("sample_id", "genotype") %in% colnames(meta)))

missing_samples <- setdiff(meta$sample_id, colnames(raw))
if (length(missing_samples) > 0) {
  stop("Sample(s) in the sample sheet not found in the counts matrix: ",
       paste(missing_samples, collapse = ", "))
}

# Select columns by name, in sample-sheet order. Never by position: column order
# in a featureCounts matrix follows the order the BAMs were passed on the command
# line and is not stable across re-runs.
count_mat <- as.matrix(raw[, meta$sample_id, drop = FALSE])

meta$genotype <- factor(meta$genotype, levels = c(ctrl_level, case_level))
if (any(is.na(meta$genotype))) {
  stop("`genotype` must be ", ctrl_level, " or ", case_level,
       "; found: ", paste(unique(as.character(meta$genotype)), collapse = ", "))
}
# Display copy: identical grouping, manuscript labels. Used for every axis, legend
# and heatmap annotation, so the group variable never leaks into a figure as WT/KO.
meta$genotype_label <- factor(ifelse(meta$genotype == ctrl_level, ctrl_name, case_name),
                              levels = c(ctrl_name, case_name))

message("[01_deseq2_deg_enrich.R] design: ", paste(table(meta$genotype), collapse = " vs "),
        " (", ctrl_name, " vs ", case_name, ")")

# -----------------------------------------------------------------------------
# 2. Pre-filtering
# -----------------------------------------------------------------------------
# A gene must reach `min_count` in at least as many libraries as the smallest
# group has replicates. This removes genes that cannot possibly be tested while
# keeping enough signal for DESeq2's own independent filtering to work.
keep <- rowSums(count_mat >= min_count) >= min(table(meta$genotype))
message("[01_deseq2_deg_enrich.R] genes retained after pre-filtering: ", sum(keep), " / ", nrow(count_mat))
count_mat <- count_mat[keep, ]

# -----------------------------------------------------------------------------
# 3. DESeq2
# -----------------------------------------------------------------------------
col_data <- data.frame(condition = meta$genotype, row.names = meta$sample_id)
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = col_data,
                              design = ~ condition)
dds <- DESeq(dds)                       # Wald test, Benjamini-Hochberg

norm_counts <- counts(dds, normalized = TRUE)

res <- results(dds, contrast = c("condition", case_level, ctrl_level))

annot <- read.delim(annot_file, row.names = 1, stringsAsFactors = FALSE)
res$gene_name    <- annot[rownames(res), "gene_name"]
res$gene_biotype <- annot[rownames(res), "gene_biotype"]

out <- cbind(as.data.frame(res), as.data.frame(norm_counts))
write.csv(out, file.path(outdir, paste0(tag, "_DE_all_genes.csv")))

# -----------------------------------------------------------------------------
# 4. Unsupervised QC: PCA and hierarchical clustering
# -----------------------------------------------------------------------------
res.pca <- PCA(t(norm_counts), graph = FALSE)
p_pca <- fviz_pca_ind(res.pca, mean.point = FALSE, label = "none",
                      habillage = meta$genotype_label, addEllipses = TRUE) +
  labs(colour = "Genotype", fill = "Genotype")
save_pdf(p_pca, "PCA", width = 6, height = 5)

hc.res <- hcut(t(norm_counts), k = hcut_k, stand = TRUE)
p_dend <- fviz_dend(hc.res, rect = TRUE, rect_border = TRUE, horiz = TRUE, cex = 0.5)
save_pdf(p_dend, "hclust", width = 7, height = 4)

# -----------------------------------------------------------------------------
# 5. DEG table, volcano, heatmap
# -----------------------------------------------------------------------------
res_pc <- out[complete.cases(out), ]
res_pc <- res_pc[res_pc$gene_biotype == "protein_coding", ]

res_pc$change <- ifelse(
  res_pc$padj <= padj_cut & abs(res_pc$log2FoldChange) > lfc_cut,
  ifelse(res_pc$log2FoldChange > lfc_cut, "Up-regulated in KO", "Down-regulated in KO"),
  "Not significant"
)
res_pc$change <- factor(res_pc$change,
                        levels = c("Down-regulated in KO", "Not significant",
                                   "Up-regulated in KO"))

n_down <- sum(res_pc$change == "Down-regulated in KO")
n_up   <- sum(res_pc$change == "Up-regulated in KO")
message("[01_deseq2_deg_enrich.R] DEGs: ", n_down + n_up, " (", n_down, " down, ", n_up, " up)")
write.csv(res_pc, file.path(outdir, paste0(tag, "_DE_protein_coding.csv")))

write.table(
  data.frame(dataset = tag, contrast = paste0(case_name, " vs ", ctrl_name),
             n_tested = nrow(res_pc), n_up = n_up, n_down = n_down,
             lfc = lfc_cut, padj = padj_cut,
             up_definition = paste0("log2FC > ", lfc_cut, " & padj <= ", padj_cut)),
  file.path(outdir, paste0(tag, "_DEG_summary.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

p_volcano <- ggplot(res_pc, aes(log2FoldChange, -log10(padj), colour = change)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_hline(yintercept = -log10(padj_cut), linetype = 5) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 5) +
  scale_colour_manual(values = c("navy", "grey", "firebrick3")) +
  labs(x = paste0("log2 fold change (", case_name, " / ", ctrl_name, ")"),
       y = "-log10 adjusted P") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
save_pdf(p_volcano, "volcano", width = 6, height = 5)

sig_genes <- rownames(res_pc)[res_pc$change != "Not significant"]
if (length(sig_genes) >= 2) {
  # Methods: "variance-stabilized expression values were row-wise Z-score normalized
  # for heatmaps". VST is computed blind to the design, then scaled per row by
  # pheatmap; normalised counts would give a different (and slightly compressed)
  # colour scale for low-count genes.
  vsd <- vst(dds, blind = TRUE)
  ann_col <- data.frame(Genotype = meta$genotype_label, row.names = meta$sample_id)
  pheatmap(SummarizedExperiment::assay(vsd)[sig_genes, meta$sample_id, drop = FALSE],
           scale = "row",
           annotation_col = ann_col,
           annotation_colors = list(Genotype = setNames(
             c("#4d4d4d", "#c0392b"), c(ctrl_name, case_name))),
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           treeheight_row = 0, treeheight_col = 0, show_rownames = FALSE,
           filename = file.path(figdir, paste0(tag, "_DEG_heatmap.pdf")),
           width = 5, height = 7)
} else {
  warning("Fewer than 2 DEGs -- heatmap skipped.")
}

# -----------------------------------------------------------------------------
# 6. GO over-representation analysis (up- and down-regulated genes separately)
# -----------------------------------------------------------------------------
run_go <- function(genes, label) {
  genes <- unique(na.omit(genes))
  if (length(genes) < 5) {
    warning("Too few genes for GO enrichment (", label, "); skipped.")
    return(NULL)
  }
  ego <- enrichGO(gene = genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                  ont = "BP", minGSSize = 10,
                  pAdjustMethod = "BH",
                  pvalueCutoff = padj_cut, qvalueCutoff = padj_cut)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    warning("No enriched GO term for ", label, ".")
    return(NULL)
  }
  write.csv(as.data.frame(ego), file.path(outdir, paste0(tag, "_GO_", label, ".csv")),
            row.names = FALSE)
  save_pdf(dotplot(ego, showCategory = 30), paste0("GO_", label), width = 8, height = 8)
  ego
}

# Direction labels are the manuscript's: enrichment of genes lower / higher in KO.
ego_down <- run_go(res_pc$gene_name[res_pc$change == "Down-regulated in KO"], "down_in_KO")
ego_up   <- run_go(res_pc$gene_name[res_pc$change == "Up-regulated in KO"],   "up_in_KO")

# -----------------------------------------------------------------------------
# 7. GSEA on the ranked gene list
# -----------------------------------------------------------------------------
# --gsea-metric stat -> rank by the signed DESeq2 Wald statistic (Methods)
# --gsea-metric lfc  -> rank by log2FoldChange (kept for sensitivity checks)
if (!gsea_metric %in% c("stat", "lfc")) stop("--gsea-metric must be stat or lfc")

res_rank <- res_pc
metric <- if (gsea_metric == "stat") res_rank$stat else res_rank$log2FoldChange
gene_list <- metric
names(gene_list) <- res_rank$gene_name
gene_list <- gene_list[!is.na(gene_list) & !is.na(names(gene_list))]
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- sort(gene_list, decreasing = TRUE)

gsego <- gseGO(geneList = gene_list, OrgDb = org.Mm.eg.db, ont = "BP",
               keyType = "SYMBOL", minGSSize = 50, maxGSSize = 500,
               pAdjustMethod = "BH", pvalueCutoff = padj_cut, verbose = FALSE)

if (!is.null(gsego) && nrow(as.data.frame(gsego)) > 0) {
  gsea_tab <- as.data.frame(gsego)
  gsea_tab$ranking_metric <- gsea_metric
  write.csv(gsea_tab, file.path(outdir, paste0(tag, "_GSEA_GO_BP.csv")),
            row.names = FALSE)
  save_pdf(dotplot(gsego, showCategory = 8, split = ".sign") + facet_grid(~ .sign),
           "GSEA_dotplot", width = 9, height = 6)
  if (gsea_term %in% gsea_tab$ID) {
    p_score <- gseaplot2(gsego, geneSetID = gsea_term, color = "firebrick3",
                         pvalue_table = TRUE)
    ggsave(file.path(figdir, paste0(tag, "_GSEA_", gsub(":", "_", gsea_term), ".pdf")),
           p_score, width = 6, height = 5)
  } else {
    warning("Term ", gsea_term, " not in the GSEA result table.")
  }
} else {
  warning("GSEA returned no enriched term.")
}

# -----------------------------------------------------------------------------
# 8. Run record
# -----------------------------------------------------------------------------
writeLines(c(paste0("tag\t", tag),
             paste0("contrast\t", case_name, " vs ", ctrl_name),
             paste0("gsea_metric\t", gsea_metric),
             paste0("min_count\t", min_count),
             paste0("lfc\t", lfc_cut), paste0("padj\t", padj_cut),
             paste0("samples\t", paste(meta$sample_id, collapse = ",")),
             capture.output(sessionInfo())),
           file.path(outdir, paste0(tag, "_sessionInfo.txt")))

message("[01_deseq2_deg_enrich.R] done. Outputs in ", normalizePath(outdir))
