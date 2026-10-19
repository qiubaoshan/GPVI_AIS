# GPVI_AIS — Omics Data Processing and Reproducibility Guide

Companion documentation for the manuscript:

> **Platelet GPVI licenses neutrophil capillary plugging and cerebral no-reflow after stroke reperfusion**
> Qiu B. et al.

This repository contains the analysis code for every omics dataset in the paper — bulk RNA-seq
(ischemic brain tissue and sorted circulating neutrophils), quantitative proteomics (platelets and plasma),
and single-cell RNA-seq (including cell–cell communication).

**This README documents the scripts; it does not duplicate them.** Each section below states what a script
does, in order, how to invoke it, which parameters matter, what it writes, and how to tell whether the run
succeeded.

---

## Table of contents

1. [Scope and how to use this repository](#1-scope-and-how-to-use-this-repository)
2. [Datasets at a glance](#2-datasets-at-a-glance)
3. [Software environment](#3-software-environment)
4. [Conventions and shared parameters](#4-conventions-and-shared-parameters)
5. [Module 1 — Bulk RNA-seq](#5-module-1--bulk-rna-seq)
6. [Module 2 — Proteomics](#6-module-2--proteomics)
7. [Module 3 — Single-cell RNA-seq](#7-module-3--single-cell-rna-seq)
8. [Module 4 — Cell–cell communication](#8-module-4--cellcell-communication)
9. [Figure-to-script map](#9-figure-to-script-map)
10. [Manuscript alignment notes](#10-manuscript-alignment-notes)
11. [Data and code availability](#13-data-and-code-availability)

---

## 1. Scope and how to use this repository

**Biological question.** Does platelet GPVI drive post-reperfusion cerebral no-reflow, and through which
transcriptional and proteomic programmes?

**Omics strategy.** Four independent molecular layers were profiled in WT versus *Gp6*-deficient mice after
transient middle cerebral artery occlusion (tMCAO):

| Layer | Matrix | Timepoint | Purpose in the paper |
|---|---|---|---|
| Bulk RNA-seq | Ipsilateral ischemic brain | 24 h after reperfusion | Tissue-level inflammatory remodelling |
| Bulk RNA-seq | Sorted neutrophils (CD45⁺CD11b⁺Ly6G⁺) | Basal and 3 h | Direct effect of GPVI on neutrophil programming |
| Proteomics | Isolated platelets | 3 h | Platelet activation / granule remodelling |
| Proteomics | Plasma | 3 h | Systemic inflammatory signature |
| scRNA-seq | Ipsilateral ischemic brain | 3 h | Cellular resolution; endothelial–neutrophil crosstalk |

**How to use it.** Every module follows the same five-part structure in this document:

```
Inputs  ->  Preprocessing  ->  Script and parameters  ->  Outputs / figures  ->  Sanity checks
```

Run the scripts in the order given. All scripts print their own run record (`*_runinfo.txt` or
`*_sessionInfo.txt`); keep those files with the results — they are what makes a run auditable six months later.

**Script inventory.** One script does one job, and its file name states the job:

| Module | Step | Script | What it does | Section |
|---|---|---|---|---|
| `01_bulk_rnaseq` | 00 | `00_align_quant.sh` | FastQC → HISAT2 → sorted BAM → featureCounts gene counts | 5.2 |
| `01_bulk_rnaseq` | 01 | `01_deseq2_deg_enrich.R` | DESeq2 DEG calling, PCA/hclust QC, VST heatmap, GO:BP, GSEA | 5.3 |
| `02_proteomics` | 01 | `01_limma_dep_enrich.R` | log₂ + quantile normalisation, limma DEP calling, rank-abundance, GO/KEGG | 6.4 |
| `03_scrnaseq` | 00 | `00_cellranger_count.sh` | `cellranger count` per library: filtered matrix and raw loom | 7.1 |
| `03_scrnaseq` | 01 | `01_qc_integrate_annotate.py` | QC, normalisation, HVG/PCA/Harmony, clustering, marker-based annotation | 7.2 |
| `04_cellchat` | 01 | `01_cellchat_compare.R` | CellChat inference per genotype and cross-genotype comparison | 8.1 |

**Naming rule.** Every script is `NN_<object>_<operation>.<ext>`, where `NN` is the execution order inside its
module. `00` is reserved for the step that turns raw reads into a count matrix (`00_align_quant.sh`,
`00_cellranger_count.sh`); everything from `01` onward analyses an existing matrix. In the single-cell chain the
outputs repeat the step number (`<tag>_01_normalized.h5ad`, `<tag>_01_annotated.h5ad`, `figures/01_*.pdf`), so a
file name there tells you which script produced it. The `--tag` argument carries the dataset id from Section 2.

**Reading order for a first-time reproducer.** Sections 2–4 (10 min), then only the module you need.
Section 10 records where the code and the Methods text describe the same analysis in different words, and
what still has to be written into the manuscript.

---

## 2. Datasets at a glance

| ID | Dataset | Assay | Groups (n) | Raw unit | Deposited at |
|---|---|---|---|---|---|
| `BULK-BRAIN` | Ischemic brain, 24 h | Bulk RNA-seq, PE150 | WT (4) vs *Gp6*-KO (4) | featureCounts gene counts | GSA **CRA048632** |
| `BULK-NEU` | Sorted neutrophils, basal / 3 h | Bulk RNA-seq, PE150 | WT vs KO, two timepoints | featureCounts gene counts | GSA **CRA048632** |
| `PROT-PLT` | Platelets, 3 h | LC–MS/MS, DDA | WT vs KO | Protein intensity | OMIX **OMIX019903** |
| `PROT-PLASMA` | Plasma, 3 h | LC–MS/MS, DDA | WT vs KO | Protein intensity | OMIX **OMIX019903** |
| `SC-BRAIN` | Ischemic brain, 3 h | 10x Genomics 3′ scRNA-seq | 1 pooled library per genotype (3 mice pooled each) | Cell Ranger filtered matrix | GSA **CRA048632** |

> **Design note on `SC-BRAIN`.** One pooled library was generated per genotype. Genotype-level single-cell
> comparisons are therefore **exploratory and descriptive**, not replicate-level inference. The manuscript
> states this; the code should be read accordingly. Do not attach a biological-replicate interpretation to
> per-genotype single-cell p-values.

Each dataset has a sample sheet declaring which sequencing library belongs to which genotype. The sheets are
in `01_bulk_rnaseq/samples/` and `02_proteomics/samples/`; Section 4 explains the columns.

---

## 3. Software environment

### Command-line

| Tool | Version | Used by |
|---|---|---|
| FastQC / HISAT2 / samtools / featureCounts | 0.11.x / 2.2.x / 1.x / 2.0.x | `01_bulk_rnaseq/00_align_quant.sh` |
| Cell Ranger | 6.x–7.x | `03_scrnaseq/00_cellranger_count.sh` |
| Proteome Discoverer | 2.4.1.15 | proteomic identification (vendor side) |

Reference genomes:

- **Bulk RNA-seq** — *Mus musculus* **GRCm39**, plus `Mus_musculus.GRCm39.gene.csv`
  (`gene_id`, `gene_name`, `gene_biotype`).
- **scRNA-seq** — **mm10 / GENCODE vM23**, restricted to `Genetype == "protein_coding"`.

> These are **different builds**. Each module's outputs stay on their own build. Convert identifiers
> explicitly before intersecting bulk and single-cell gene sets; do not intersect Ensembl IDs across modules.

### R (v4.4.0)

```r
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Mm.eg.db", "AnnotationDbi",
                       "enrichplot", "limma", "CellChat"))
install.packages(c("ggplot2", "pheatmap", "FactoMineR", "factoextra", "ggrepel",
                   "dplyr", "Seurat", "patchwork", "dior"))
```

`dior` / `diopy` are used only as H5 read/write shuttles between Python and R. Any equivalent converter
(`zellkonverter`, `SeuratDisk`) works; the analysis is unaffected.

### Python (3.9–3.11)

```bash
pip install scanpy anndata harmonypy \
            matplotlib pandas numpy seaborn
```

### Record your session

The R scripts write `*_sessionInfo.txt`; `01_qc_integrate_annotate.py` writes `*_01_runinfo.txt`. Keep both with the
results. Bioconductor minor upgrades change `clusterProfiler` term ordering, which is cosmetic but breaks
byte-identical reproduction.

---

## 4. Conventions and shared parameters

### Group naming

Genotypes are referred to by two labels, and the split matters:

| Role | Analysis label | Display label | Where it appears |
|---|---|---|---|
| Control | `WT` | `WT` | `condition`, `genotype`, file names, `--control` |
| *Gp6* knockout | `KO` | `Gp6-KO` | `condition`, `--case`; `Gp6-KO` in figures |

The sequencing libraries carry the legacy submission codes `NS` (WT) and `YY` (KO). Those codes exist **only**
in file names and in the `sample_id` column of the sample sheets. Every script maps them to the analysis and
display labels declared on the command line, so no axis, legend or table ever shows `NS` / `YY`.

### Sample sheets

Tab-separated, one row per library:

| Column | Meaning |
|---|---|
| `sample_id` | Library / column name as delivered (e.g. `NS1`, `YY3`) |
| `genotype` | `WT` or `KO` — the column every script reads |
| `timepoint` | `basal`, `3h`, `24h` (documentation and grouping checks) |
| `tissue` / `matrix` | `brain`, `neutrophil`, `platelet`, `plasma` |

`01_deseq2_deg_enrich.R` and `01_limma_dep_enrich.R` select matrix columns **by name, in sheet order**, and stop if a sheet entry
is missing from the matrix. Positional assignment is not used anywhere, because featureCounts column order
follows the order the BAMs were passed on the command line and is not stable across re-runs.

### Shared parameters

| Item | Value | Where it is set |
|---|---|---|
| Contrast direction | KO vs WT (positive log₂FC = higher in KO) | `--control WT --case KO` |
| Multiple-testing correction | Benjamini–Hochberg | all scripts |
| Bulk DE thresholds | \|log₂FC\| > 1 **and** adj. P < 0.05 | `01_deseq2_deg_enrich.R --lfc 1 --padj 0.05` |
| Bulk GSEA ranking | signed DESeq2 Wald statistic | `01_deseq2_deg_enrich.R --gsea-metric stat` |
| Proteomics DE thresholds | \|log₂FC\| > 0.58 **and** adj. P < 0.05 | `01_limma_dep_enrich.R --lfc 0.58 --pcol adj.P.Val` |
| scRNA DE | Wilcoxon rank-sum, BH | cluster markers in `01_qc_integrate_annotate.py`; the genotype-level DE behind 4D is not yet scripted (10.6) |
| Enrichment cut-offs | P < 0.05 and q < 0.05, BH | all enrichment scripts |
| Seed | 0 | `--seed` in Python scripts; set `set.seed()` if you add resampling in R |

**Sign convention in every figure and table.** Positive log₂FC = higher in *Gp6*-KO. The manuscript reports
genes as "down-regulated in KO", which corresponds to negative log₂FC; the DEG tables label those rows
`Down-regulated in KO` so the direction cannot be misread. Every figure also states the contrast in its axis
title (`log2 fold change (Gp6-KO / WT)`).

### Two rules that keep runs comparable

1. **Never edit a script to change a threshold.** Pass the flag; the value then appears in the run record.
2. **Never rename a library after `cellranger count`,** and never change `--ctrl-label` / `--case-label`
   between scripts. Single-cell barcodes are `<barcode>_<label>`; a changed label breaks the joins silently.

---

## 5. Module 1 — Bulk RNA-seq

Two libraries were sequenced with an identical pipeline and analysed with the same script. They differ only
in the counts matrix and the sample sheet.

### 5.1 Wet-lab and sequencing

| Step | Detail |
|---|---|
| Tissue | Ipsilateral ischemic brain, consistent anatomical boundaries, 24 h after reperfusion |
| Cells | CD45⁺CD11b⁺Ly6G⁺ neutrophils sorted on BD FACSAria III; basal or 3 h |
| RNA — tissue | FastPure Cell/Tissue Total RNA Isolation Kit V2 (Vazyme, RC112) |
| RNA — neutrophils | PicoPure RNA Isolation Kit (Thermo Fisher, KIT0204) |
| QC | Agilent 2100 Bioanalyzer |
| Library | SMARTer Stranded Total RNA-Seq Kit v2 — Pico Input Mammalian (Takara, 634413); ZapR v2 rRNA depletion |
| Sequencing | Illumina NovaSeq 6000, PE150, ~6–9 Gb per sample (Novogene) |

### 5.2 Alignment and quantification — `01_bulk_rnaseq/00_align_quant.sh`

```bash
cd 01_bulk_rnaseq
SAMPLES=samples/BULK-BRAIN.tsv THREADS=16 \
  ./00_align_quant.sh ../../fastq ../../reference/GRCm39 ../../results/BULK-BRAIN
```

| Step | Tool | Key flags |
|---|---|---|
| Read QC | FastQC | — |
| Alignment | HISAT2 | splice-aware, default |
| Sorting | samtools | — |
| Quantification | featureCounts | `-p --countReadPairs -s 2` |

**The `-s 2` flag is load-bearing.** The libraries are reverse-stranded (SMARTer). Using `-s 0` halves the
usable reads and shifts gene counts, which propagates into every downstream figure. If you change library
kits, re-check strandedness before trusting any DEG number.

The script reads sample names from the sample sheet, so the column order of `all_samples.gene_counts.txt`
always matches the sheet consumed by `01_deseq2_deg_enrich.R`.

### 5.3 Differential expression and enrichment — `01_bulk_rnaseq/01_deseq2_deg_enrich.R`

```bash
cd 01_bulk_rnaseq
Rscript 01_deseq2_deg_enrich.R \
  --counts ../results/BULK-BRAIN/all_samples.gene_counts.txt \
  --samples samples/BULK-BRAIN.tsv \
  --annot Mus_musculus.GRCm39.gene.csv \
  --outdir ../results/BULK-BRAIN --tag BULK-BRAIN
```

The neutrophil libraries are reported as two contrasts (basal, 3 h), so they get two sample sheets over the
same counts matrix:

```bash
Rscript 01_deseq2_deg_enrich.R --counts ../results/BULK-NEU/all_samples.gene_counts.txt \
  --samples samples/BULK-NEU-basal.tsv --annot Mus_musculus.GRCm39.gene.csv \
  --outdir ../results/BULK-NEU-basal --tag BULK-NEU-basal

Rscript 01_deseq2_deg_enrich.R --counts ../results/BULK-NEU/all_samples.gene_counts.txt \
  --samples samples/BULK-NEU-3h.tsv --annot Mus_musculus.GRCm39.gene.csv \
  --outdir ../results/BULK-NEU-3h --tag BULK-NEU-3h
```

`01_deseq2_deg_enrich.R` selects matrix columns by the sheet, so both runs read the same matrix and each model sees only its
own timepoint. Never put both timepoints in one sheet: the design is `~ condition`, and timepoint would be
silently absorbed into the genotype contrast.

**What the script does, in order**

1. Loads the counts matrix and joins it to the sample sheet **by sample name**, with a hard stop
   (`stopifnot`) if a sheet sample is missing from the matrix.
2. Renames `genotype` to the analysis factor with levels `WT`, `KO`, plus a display factor carrying the
   figure labels. `condition` is the only grouping variable in the model.
3. Pre-filters genes: count ≥ 10 in at least as many libraries as the smallest group has replicates.
4. Fits `DESeq2` with design `~ condition`; Wald test; BH correction; contrast `KO` vs `WT`.
5. Annotates results with `gene_name` and `gene_biotype`; writes the full table with normalised counts.
6. Unsupervised QC: PCA (FactoMineR) with genotype ellipses, and hierarchical clustering dendrogram.
7. Calls DEGs on **protein-coding genes**: \|log₂FC\| > 1 and adj. P < 0.05; volcano plot; writes
   `<tag>_DEG_summary.tsv` with the counts and the thresholds that produced them.
8. Row Z-scored heatmap of DEGs on variance-stabilised values (`vst`, blind to the design), annotated by
   genotype, no row/column dendrogram.
9. GO:BP over-representation on up- and down-regulated genes **separately** (`enrichGO`).
10. GSEA on the ranked gene list (`gseGO`), plus a running-score plot for one term.
11. Writes `sessionInfo()` and the effective settings to the run record.

**Parameter table**

| Flag | Default | Effect |
|---|---|---|
| `--control` / `--case` | `WT` / `KO` | Levels of `genotype`; also sets contrast direction |
| `--control-label` / `--case-label` | `WT` / `Gp6-KO` | Labels printed in every figure |
| `--min-count` | 10 | Pre-filter stringency |
| `--lfc` / `--padj` | 1 / 0.05 | DEG thresholds |
| `--hcut-k` | 4 | Cuts the QC dendrogram |
| `--gsea-metric` | `stat` | Signed Wald statistic (Methods). `lfc` ranks by log₂FC; kept for sensitivity checks only |
| `--gsea-term` | `GO:0032755` | Term drawn as a running score |

**Decisions worth knowing**

- Results are **not** LFC-shrunk. Volcano axes use raw `log2FoldChange`. Switching to `lfcShrink()` changes
  the DEG count and the volcano shape; do not mix the two within one figure.
- `results()` applies independent filtering by default, so `padj` is `NA` for low-count genes; those rows are
  dropped before plotting. Report this, because it shrinks the tested universe relative to the annotation.
- Up- and down-regulated sets must be enriched separately. A combined list cancels the biology — the paper's
  central observation is directional (broad suppression of inflammatory programmes in KO).
- GSEA requires unique gene symbols; multi-mapping Ensembl IDs otherwise pick an arbitrary symbol.

**Outputs**

| File | Content |
|---|---|
| `<tag>_DE_all_genes.csv` | Full DESeq2 table + normalised counts |
| `<tag>_DE_protein_coding.csv` | Protein-coding subset with `change` labels |
| `<tag>_DEG_summary.tsv` | DEG counts, thresholds, contrast — paste into the figure legend |
| `<tag>_GO_down_in_KO.csv` / `_GO_up_in_KO.csv` | GO:BP over-representation, per direction |
| `<tag>_GSEA_GO_BP.csv` | Ranked GSEA result, including `ranking_metric` |
| `figures/<tag>_{PCA,hclust,volcano,DEG_heatmap,GO_*,GSEA_*}.pdf` | Panels |
| `<tag>_sessionInfo.txt` | Run record |

### 5.4 Sanity checks

| Dataset | Published checkpoint |
|---|---|
| `BULK-BRAIN` | **69 DEGs (67 down, 2 up)** — Figure S5B |
| `BULK-NEU`, basal | **1 gene** differs — Figure S5E |
| `BULK-NEU`, 3 h | **1,334 DEGs (1,150 down, 184 up)** — Figure S5E |
| PCA | Genotypes separate on PC1/PC2 |

If the DEG count is far off, the usual causes are, in order: `-s 0` instead of `-s 2`; the `gene_biotype`
filter skipped; group labels swapped in the sample sheet; a different `min-count` pre-filter than the
deposited run.

---

## 6. Module 2 — Proteomics

The two matrices share sample-preparation chemistry but differ in complexity, dynamic range and therefore in
what normalisation means. **Process them separately; never pool platelet and plasma intensities in one
normalisation step.**

### 6.1 Sample collection

| Matrix | Procedure |
|---|---|
| Platelets | Cardiac puncture into sodium citrate; 300 × g 5 min RT; supernatant re-spun 800 × g 5 min RT; pellet in 100 µL HEPES–Tyrode's + 20 µL sodium citrate |
| Plasma | Anticoagulated retro-orbital blood, 1,200 × g 10 min 4 °C; aliquoted, −80 °C |
| Timepoint | 3 h after reperfusion |
| Protein assay | BCA after clarification at 12,000 × g 10 min 4 °C |

### 6.2 LC–MS/MS acquisition

| Parameter | Value |
|---|---|
| Digestion | Trypsin 1:50 (w/w) overnight; alkylation 11 mM iodoacetamide, 15 min, dark |
| System | EASY-nLC 1200; 90 min gradient at 500 nL/min |
| Instrument | Orbitrap Exploris 480 |
| Acquisition | DDA; full MS 400–1,200 m/z at 60,000 resolution |
| MS/MS | HCD, top 15 precursors, 30,000 resolution |
| Provider | Jingjie PTM BioLab (Hangzhou) |

### 6.3 Database search

Performed in Proteome Discoverer 2.4.1.15; the parameters below must be matched before exporting the matrix
that the R script consumes.

| Parameter | Setting |
|---|---|
| Database | *Mus musculus* reference proteome |
| Enzyme / missed cleavages | Trypsin / 2 |
| Precursor / fragment tolerance | 10 ppm / 0.02 Da |
| Fixed modification | Carbamidomethyl (C) |
| Variable modifications | Oxidation (M), deamidation (N/Q), N-terminal modifications |
| FDR | 1% at PSM, peptide and protein level |

Export a **protein-level intensity matrix**: `Protein`, `Gene`, one column per sample, tab-separated.

### 6.4 Normalisation, statistics and enrichment — `02_proteomics/01_limma_dep_enrich.R`

```bash
cd 02_proteomics
Rscript 01_limma_dep_enrich.R --matrix platelet_protein_intensity.tsv \
  --samples samples/PROT-PLT.tsv --outdir ../results/PROT-PLT --tag PROT-PLT

Rscript 01_limma_dep_enrich.R --matrix plasma_protein_intensity.tsv \
  --samples samples/PROT-PLASMA.tsv --outdir ../results/PROT-PLASMA --tag PROT-PLASMA
```

**What the script does, in order**

1. Reads the matrix; resolves sample columns from `--samples`, mapping `sample_id` to `genotype`.
   Without a sheet it falls back to `WT*` / `KO*` column prefixes.
2. Filters proteins detected in ≥ 70% of the smallest group.
3. `log2` transform (skip with `--already-log2`), then quantile-normalise **within matrix**.
4. `limma`: `lmFit` → contrast `KO - WT` → `eBayes` → `topTable` (all proteins).
5. Volcano plot at \|log₂FC\| > 0.58 and the chosen P column; writes `<tag>_DEP_summary.tsv` with the count
   and the exact column used.
6. Rank-abundance: mean log₁₀ abundance per genotype, ranked by WT abundance, written as a table plus line
   plot.
7. Row Z-scored heatmap of the top differential proteins, annotated by genotype.
8. GO:BP and KEGG enrichment (`clusterProfiler`), with a **rich-factor denominator computed against detected
   proteins**, not the full annotation universe.
9. Writes `sessionInfo()` and the effective settings.

**Parameter table**

| Flag | Default | Effect |
|---|---|---|
| `--control` / `--case` | `WT` / `KO` | Levels of `genotype` |
| `--control-label` / `--case-label` | `WT` / `Gp6-KO` | Labels printed in every figure |
| `--min-frac` | 0.7 | Detection filter |
| `--lfc` / `--p` | 0.58 / 0.05 | Differential-abundance thresholds |
| `--pcol` | `adj.P.Val` | Significance column; `P.Value` is available but must then be stated in the legend |
| `--top-n` | 40 | Proteins in the heatmap |
| `--already-log2` | off | Skip the log transform |
| `--samples` | — | Sample sheet instead of name prefixes |

**Rich factor, defined explicitly.** The figure legends define rich factor as *differential proteins ÷ all
detected proteins annotated to the pathway*. `clusterProfiler`'s native `BgRatio` uses the whole annotation
universe, which is a different denominator. The script recomputes it from the detected background and writes
the value into the result table, so the dot-plot x-axis is reproducible.

**Outputs:** `<tag>_limma_all_proteins.csv`, `<tag>_DEP_summary.tsv`, `<tag>_rank_abundance.csv`,
`<tag>_GO_BP.csv`, `<tag>_KEGG.csv`,
`figures/<tag>_{volcano,rank_abundance,heatmap,GO_BP_richfactor,KEGG}.pdf`, `<tag>_sessionInfo.txt`.

### 6.5 Sanity checks

| Dataset | Published checkpoint |
|---|---|
| `PROT-PLT` | Global reduction in protein abundance after stroke; pathway-level remodelling of activation, GPVI signalling, surface and granule proteins (Fig. 1I/1J, S1D/E) |
| `PROT-PLASMA` | **13 DEPs: 10 down, 3 up** (Fig. S5H); enrichment of leukocyte migration, positive chemotaxis, inflammatory response, myeloid leukocyte migration, platelet activation, cell chemotaxis (Fig. S5I) |

A plasma DEP count far above 13 usually means the missing-value filter was skipped, the threshold was relaxed
against a much larger protein universe, or the significance column was switched to raw P. Check the
`p_column` field in `<tag>_DEP_summary.tsv` first; 13 proteins is a consequence of the adjusted-P filter, and
the S1D/S5H legends should be reworded to say so (Section 10.3).

---

## 7. Module 3 — Single-cell RNA-seq

`00_cellranger_count.sh` produces one count matrix per library. One Python script (`01`) then runs the
single-cell analysis the manuscript reports — QC, integration, clustering, annotation, marker detection — and
writes a full-gene normalised object plus an annotated one, which is what CellChat reads in module 4.

> **Scope note.** The within-population DE and GO:BP / KEGG step behind Figure 4D was run interactively in the
> deposited notebook and has no script here yet. It is listed as item 10.6 in Section 10 rather than being
> implied by this module.

### 7.1 Library preparation and counting — `03_scrnaseq/00_cellranger_count.sh`

```bash
cd 03_scrnaseq
LIBS="NS1 YY1" ./00_cellranger_count.sh ../../reference/refdata-gex-mm10-2020-A ../../results/10x
```

Runs `cellranger count` once per library. No BAM is written: every downstream step in this paper reads the
count matrix or the loom, not alignments.

> **Library names vs genotypes.** Cell Ranger ids stay as sequenced (`NS1`, `YY1`); the genotype labels are
> declared when the Python scripts are called (`--ctrl-label WT --case-label KO`). Barcodes are built as
> `<barcode>_<label>`, so the label declared here must be the one used in `01_qc_integrate_annotate.py`; changing
> it mid-pipeline silently breaks every join into the CellChat export.

### 7.2 QC, integration, clustering, annotation — `03_scrnaseq/01_qc_integrate_annotate.py`

```bash
python 01_qc_integrate_annotate.py \
  --ctrl-loom 3d_NS1_S4_L004_raw.loom \
  --case-loom 3d_YY1_S3_L003_raw.loom \
  --annot mm10_gencode.vM23.Gene_annot.csv \
  --cluster-map cluster2cell.tsv \
  --outdir ../results/SC-BRAIN --tag SC-BRAIN
```

**What the script does, in order**

1. Reads one loom per genotype, tags `condition` (`WT` / `KO`), and builds unique cell IDs
   (`<barcode>_<label>`).
2. Restricts to protein-coding genes, keyed on gene symbol.
3. Cell and gene filters, mitochondrial content QC, violin and scatter diagnostics.
4. `normalize_total(1e4)` + `log1p`.
5. Writes the **full-gene** normalised object before the HVG subset; CellChat needs a log-normalised matrix
   over all genes.
6. Highly variable genes → PCA → Harmony on `condition` → kNN graph → UMAP → Louvain.
7. Marker detection (Wilcoxon) and marker dot plots.
8. Copies the cluster labels and both embeddings back to the full-gene object.
9. Applies `cluster2cell.tsv` to assign `cell_type`, failing loudly if a cluster has no entry.
10. Cell-type composition (ratio-to-expected) and a per-genotype cell summary.
11. Writes the run record, including QC values, integration method and PC count.

**QC parameters — these are the values stated in Methods**

| Parameter | Value | Flag |
|---|---|---|
| Minimum genes per cell | 200 | `--min-genes` |
| Maximum genes per cell | 5,000 | `--max-genes` |
| Maximum mitochondrial reads | 20% | `--max-mt` |
| Minimum cells per gene | 3 | `--min-cells` |

The values actually used are written to `<tag>_01_runinfo.txt`. Report that file's numbers, not this table,
if a run overrides any of them.

**Other parameters**

| Flag | Default | Note |
|---|---|---|
| `--integration` | `harmony` | `none` leaves the PCA space uncorrected; see the caveat below |
| `--hvg` | 3,000 | `cell_ranger` flavour, appropriate for 10x |
| `--n-pcs` | 30 | Dimensionality stated in Methods |
| `--n-neighbors` | 100 | |
| `--resolution` | 1.0 | 15 clusters expected |
| `--cluster-map` | `cluster2cell.tsv` | Cluster → cell-type decisions live in a table, not in code |

> **Caveat on correcting by genotype.** Harmony is run on `condition`, which is the genotype. With one pooled
> library per genotype there is no batch structure to separate from the biological factor, so the correction
> can shrink the very difference the panel is about. `--integration none` runs the graph on uncorrected PCs.
> Whichever is used, state it in Methods — and state it once.

**Outputs:** `<tag>_01_normalized.h5ad`, `<tag>_01_annotated.h5ad`, `<tag>_composition.csv`,
`<tag>_cell_summary.csv`, `<tag>_01_runinfo.txt`, and `figures/01_*.pdf`.

### 7.3 Sanity checks

| Check | Expected |
|---|---|
| Clusters at resolution 1.0 | 15 |
| Lineages present | MG, MP, NEU, AST, OL, NT, END, PFB |
| Cell types vs `cluster2cell.tsv` | No unmapped cluster (the script exits if there is one) |
| EC and neutrophil DEGs (4D) | Coordinated **down-regulation** of leukocyte migration, chemotaxis, inflammatory signalling, cell adhesion in KO — compare against the notebook run until 10.6 is scripted |

---

## 8. Module 4 — Cell–cell communication

### 8.1 Inference and comparison — `04_cellchat/01_cellchat_compare.R`

```bash
cd 04_cellchat
Rscript 01_cellchat_compare.R \
  --h5 ../results/SC-BRAIN/SC-BRAIN_01_annotated.h5 \
  --outdir ../results/CELLCHAT --tag CELLCHAT \
  --senders END,PFB,AST \
  --pathways NRXN,CALCR,VISFATIN
```

**What the script does, in order**

1. Reads the annotated object and splits it by `condition` (`WT` / `KO`).
2. Creates one CellChat object per genotype (`group.by = "cell_type"`) against `CellChatDB.mouse`.
3. Runs the full inference chain **independently per genotype**: over-expressed genes → over-expressed
   interactions → `computeCommunProb(raw.use = TRUE, nboot = 100)` → `filterCommunication(min.cells = 10)`
   → pathway aggregation → network aggregation → centrality → functional and structural similarity/embedding/
   clustering.
4. Exports each genotype's significant L-R table to CSV. Figures are not data; export the table.
5. Merges the two objects, named with the display labels (`WT`, `Gp6-KO`), and computes interaction
   counts/weights, pairwise network similarity, and `rankNet`.
6. Differential analysis across datasets: `identifyOverExpressedGenes(group.dataset = "datasets",
   pos.dataset = "Gp6-KO")` → `netMappingDEG` → `subsetCommunication` up/down at ligand/receptor log₂FC ±0.1.
7. Bubble plots for each sender in `--senders`, plus differential-interaction heatmaps.
8. Pathway heatmaps and ligand/receptor expression violins for each pathway in `--pathways`.
9. Saves the merged object and `sessionInfo()`.

| Flag | Default | Effect |
|---|---|---|
| `--ctrl` / `--case` | `WT` / `KO` | Levels of `condition` |
| `--ctrl-name` / `--case-name` | `WT` / `Gp6-KO` | Dataset names, i.e. every legend in the panel |
| `--min-cells` | 10 | Interactions from populations below this are dropped |
| `--nboot` | 100 | Permutations behind the communication probability |
| `--workers` | 10 | `future` multisession workers |
| `--senders` | `END,PFB` | Resolved **by name**, not by ident number |
| `--pathways` | `NRXN,CALCR,VISFATIN` | Heatmaps and gene-expression violins |

> **Never hard-code ident indices.** Idents carry a numeric prefix (`9:AST`, `12:PFB`) that changes whenever
> clustering is re-run. The script resolves names via `ident_index()` and stops with the list of available
> idents if a name is missing. A silently wrong index produces a plausible but meaningless panel.

> **Never fit one joint model and split it.** Inference must be independent per genotype; merging is a
> comparison step, not an inference step.

**Outputs:** `<tag>_cellchat_{WT,KO}.rds`, `<tag>_cellchat_merged.rds`,
`<tag>_communication_{WT,KO}.csv`, `<tag>_L-R_up_in_{KO,WT}.csv`,
`figures/<tag>_{interaction_counts,embedding_functional,embedding_structural,ranknet,bubble_*,diffInteraction_*,pathway_*,genes_*}.pdf`.

### 8.2 Sanity checks

| Check | Expected |
|---|---|
| Significant L-R interactions per genotype | Reported in the console log; record both |
| Endothelial-origin signalling | Fewer and weaker interactions in KO |
| Endothelial → neutrophil axis | The most affected sender–target pair |
| `NRXN` / `CALCR` / `VISFATIN` present in the inferred pathways | If not, `--pathways` entries are skipped with a warning |

---

## 9. Figure-to-script map

| Figure | Panel | Dataset | Script | Flag / step |
|---|---|---|---|---|
| 1I | Rank-abundance, platelet proteome | `PROT-PLT` | `01_limma_dep_enrich.R` | step 6 |
| 1J | Platelet protein heatmap | `PROT-PLT` | `01_limma_dep_enrich.R` | step 7 |
| S1D | Platelet proteome volcano | `PROT-PLT` | `01_limma_dep_enrich.R` | `--lfc 0.58` |
| S1E | KEGG enrichment, platelet DEPs | `PROT-PLT` | `01_limma_dep_enrich.R` | step 8 |
| 4A | Neutrophil RNA-seq workflow | `BULK-NEU` | `00_align_quant.sh` | — |
| 4B | GO:BP, KO vs WT neutrophils | `BULK-NEU` | `01_deseq2_deg_enrich.R` | steps 9–10 |
| 4C | scRNA UMAP, major populations | `SC-BRAIN` | `01_qc_integrate_annotate.py` | default QC |
| 4D | GO:BP, endothelial and neutrophil DEGs | `SC-BRAIN` | not yet scripted — see 10.6 | — |
| 4E | CellChat endothelial-origin signalling | `SC-BRAIN` | `01_cellchat_compare.R` | `--senders END` |
| S5A–S5C | Brain RNA-seq workflow, volcano, heatmap | `BULK-BRAIN` | `00_align_quant.sh`, `01_deseq2_deg_enrich.R` | — |
| S5D | GO:BP, brain DEGs | `BULK-BRAIN` | `01_deseq2_deg_enrich.R` | step 9 |
| S5E | Neutrophil volcano, basal vs 3 h | `BULK-NEU` | `01_deseq2_deg_enrich.R` | two runs |
| S5F | Neutrophil DEG heatmap | `BULK-NEU` | `01_deseq2_deg_enrich.R` | step 8 |
| S5G | GSEA, adhesion terms | `BULK-NEU` | `01_deseq2_deg_enrich.R` | `--gsea-term` |
| S5H | Plasma proteome volcano | `PROT-PLASMA` | `01_limma_dep_enrich.R` | `--lfc 0.58` |
| S5I | GO:BP, plasma DEPs | `PROT-PLASMA` | `01_limma_dep_enrich.R` | step 8 |

`BULK-BRAIN` and `BULK-NEU` share one script. Which dataset a run belongs to is defined entirely by the
sample sheet and the counts matrix — keep both next to the output directory.

---

## 10. Manuscript alignment notes

The code and the Methods text are aligned on the analysis. Six points are listed below; 10.4 is already closed
by removing the analysis, so five still need one sentence or one script before the repository is made public.

| # | Point | Current state | Action |
|---|---|---|---|
| 10.1 | Integration method | Methods says Seurat reciprocal PCA; `01_qc_integrate_annotate.py` uses Harmony on 30 PCs | Either change the Methods sentence to Scanpy + Harmony + 30 PCs, or rerun with Seurat rPCA. Do not report "Seurat" for a Scanpy run |
| 10.2 | Genome build | Bulk uses GRCm39; single-cell uses mm10 / GENCODE vM23 | Add one clause to Methods. Both are correct for their assay; the risk is only naive ID intersection across modules |
| 10.3 | Proteomics significance column | Methods says limma with BH; the S1D and S5H legends say `P < 0.05`. The script uses `adj.P.Val`, the column consistent with the reported 13 plasma DEPs | Write "adjusted P" in Methods **and** in both legends, or regenerate the panels with raw P — not one of each |
| 10.4 | Undescribed analyses | **Closed by removal.** MELD, astrocyte substates and RNA velocity were in the deposited pipeline but nowhere in the draft | They are no longer part of this package. If either analysis is added to the paper later, restore the script and the Methods paragraph together — a repo that outruns the text invites the question |
| 10.5 | Volcano threshold wording | The S5B / S5E legends describe the horizontal line as `p = 0.05` while DEG calling uses adjusted P; the S5H legend says ±1 although the proteomics threshold is 0.58 | Reword each legend to the threshold actually used. As written, panels and DEG counts appear to disagree |
| 10.6 | 4D is notebook-only | The endothelial / neutrophil DEG and GO:BP–KEGG step of Figure 4D was run interactively in the deposited notebook (endothelial subcluster 1 vs 2, Wilcoxon + `gseapy`) and has no script in this package | Script it before publication, or state in the legend that the panel derives from the notebook. A repository should not omit a panel it claims to reproduce |

Two further conventions are recorded in the scripts but not currently in the manuscript, and are worth one
sentence each because they change reported numbers:

- **GSEA ranking metric** — signed DESeq2 Wald statistic (not log₂FC). `--gsea-metric lfc` exists only for
  sensitivity checks; the default and the reported result use `stat`.
- **Rich factor denominator** — differential proteins ÷ detected proteins annotated to the pathway, not the
  full annotation universe. The value is written into `<tag>_GO_BP.csv` so the axis can be rebuilt.

### Reproducibility gaps that remain

- `pheatmap` output depends on `pheatmap`'s own clustering when `treeheight_*` is 0 but `cluster_*` is not;
  both are set to `FALSE`/0 consistently here — do not "improve" one without the other.
- CellChat's `netEmbedding` uses UMAP with a random seed. Set `set.seed()` before the embedding step if you
  need identical 2D layouts across runs.

---

## 11. Data and code availability

| Resource | Identifier |
|---|---|
| Bulk and single-cell RNA-seq | GSA CRA048632|
| Platelet and plasma proteomics | OMIX OMIX019903|
| Analysis code | https://github.com/qiubaoshan/GPVI_AIS |
| Clinical trial (ESPRIT) | NCT05836753 |

Raw sequencing and proteomic data are deposited. Processed matrices, the `cluster2cell.tsv` annotation table
and figure source data should be deposited alongside the code, so that the checks in Section 11 can be run
without re-processing raw reads.

---

### Document status

Version 1.4 — the package now contains only the analyses the manuscript reports. MELD, astrocyte substates and
RNA velocity were removed together with their three scripts, and `00_cellranger_velocyto.sh` became
`00_cellranger_count.sh` (no velocyto step, no BAM). Six scripts remain: one per analysis that appears in a
figure, Methods paragraph or figure legend.
Version 1.3 named every script after its task and execution order (`NN_<object>_<operation>`); version 1.2
aligned the code with the manuscript's genotype labels (`WT` / *Gp6*-KO) and with the analysis settings stated
in Methods. Derived from the manuscript draft `GPVI_manuscript_Baoshan Qiu2.docx` and the originally deposited
scripts (`pip.R`, `cellchat.R`, `scRNA-seq_pipe.ipynb`), now reorganised into this repository. The five open
items in Section 10 must be settled before the repository is made public.
