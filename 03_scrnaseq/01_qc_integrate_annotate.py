#!/usr/bin/env python
"""
01_qc_integrate_annotate.py -- scRNA-seq: object assembly, QC, integration,
                              clustering and marker-based annotation.

Produces the annotated object that CellChat consumes in module 4, and the UMAP
panels of Figure 4C.

GENOTYPE LABELS. `condition` holds the analysis labels WT / KO. `genotype` holds
the display labels printed in figures (WT / Gp6-KO). Every axis, legend and table
uses `genotype`, so the manuscript wording never has to be patched downstream.
The legacy library codes NS / YY appear only in the loom file names.

QC. One parameter set, the one stated in Methods: >= 200 genes, <= 5,000 genes,
<= 20% mitochondrial reads, genes detected in >= 3 cells. Each value can still be
overridden from the command line, and the values actually used are written to
`<tag>_01_runinfo.txt`.

INTEGRATION. `--integration harmony` (default) embeds the two libraries with
Harmony on 30 PCs, the dimensionality stated in Methods. `--integration none`
leaves the PCA space uncorrected; use it when the genotype difference itself is
the signal under test, because correcting on `condition` regresses out exactly
that difference.

Usage
    python 01_qc_integrate_annotate.py \
        --ctrl-loom 3d_NS1_S4_L004_raw.loom \
        --case-loom 3d_YY1_S3_L003_raw.loom \
        --annot mm10_gencode.vM23.Gene_annot.csv \
        --cluster-map cluster2cell.tsv \
        --outdir results/SC-BRAIN \
        --tag SC-BRAIN
"""

import argparse
import contextlib
import io
import os
import sys

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce

sc.settings.verbosity = 1

# Manuscript QC thresholds (Methods, "Single-cell RNA sequencing").
QC_DEFAULTS = dict(min_genes=200, max_genes=5000, max_mt=20.0, min_cells=3)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ctrl-loom", required=True, help="Cell Ranger loom, WT library")
    p.add_argument("--case-loom", required=True, help="Cell Ranger loom, Gp6-KO library")
    p.add_argument("--annot", required=True, help="GENCODE gene annotation CSV (mm10 vM23)")
    p.add_argument("--cluster-map", default="cluster2cell.tsv",
                   help="cluster -> cell_type table (louvain, cell_type)")
    p.add_argument("--outdir", default="results/SC-BRAIN")
    p.add_argument("--tag", default="SC-BRAIN")
    p.add_argument("--ctrl-label", default="WT", help="analysis label for the control")
    p.add_argument("--case-label", default="KO", help="analysis label for the knockout")
    p.add_argument("--ctrl-name", default="WT", help="display label for the control")
    p.add_argument("--case-name", default="Gp6-KO", help="display label for the knockout")
    p.add_argument("--min-genes", type=int, default=QC_DEFAULTS["min_genes"])
    p.add_argument("--max-genes", type=int, default=QC_DEFAULTS["max_genes"])
    p.add_argument("--max-mt", type=float, default=QC_DEFAULTS["max_mt"])
    p.add_argument("--min-cells", type=int, default=QC_DEFAULTS["min_cells"])
    p.add_argument("--integration", choices=["harmony", "none"], default="harmony")
    p.add_argument("--hvg", type=int, default=3000)
    p.add_argument("--n-pcs", type=int, default=30,
                   help="PCs kept for the graph (Methods: 30)")
    p.add_argument("--n-neighbors", type=int, default=100)
    p.add_argument("--resolution", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args(argv)


def read_loom(path, label):
    """Read one Cell Ranger loom and give it condition-tagged, unique names."""
    a = sc.read_loom(path)
    a.obs["condition"] = label
    if "barcodes" in a.obs:
        a.obs["cell_id"] = [f"{b}_{label}" for b in a.obs["barcodes"]]
    else:
        a.obs["cell_id"] = [f"{b}_{label}" for b in a.obs_names]

    if "gene_ids" in a.var:
        a.var["gene_id"] = a.var["gene_ids"].values
        a.var_names = a.var["gene_id"].astype(str)
    a.obs_names = a.obs["cell_id"].astype(str)
    a.var_names_make_unique()
    return a


def to_protein_coding(adata, annot_path):
    """Restrict to protein-coding genes, keyed on gene symbol."""
    annot = pd.read_csv(annot_path)
    gene_col = "Gene_name" if "Gene_name" in annot.columns else "gene_name"
    id_col = "GeneID" if "GeneID" in annot.columns else "gene_id"
    type_col = "Genetype" if "Genetype" in annot.columns else "gene_type"

    pc = annot[annot[type_col] == "protein_coding"].drop_duplicates(subset=[id_col])
    adata = adata[:, adata.var["gene_id"].isin(pc[id_col])].copy()

    mapper = dict(zip(pc[id_col].astype(str), pc[gene_col].astype(str)))
    adata.var["gene_name"] = [mapper.get(g, g) for g in adata.var["gene_id"].astype(str)]
    adata.var_names = adata.var["gene_name"].astype(str)
    adata.var_names_make_unique()
    return adata


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(os.path.join(args.outdir, "figures"), exist_ok=True)
    figdir = os.path.join(args.outdir, "figures")

    qc = {k: getattr(args, k) for k in QC_DEFAULTS}
    print(f"[01] QC: {qc}")

    # -- 1. assemble ---------------------------------------------------------
    ctrl = read_loom(args.ctrl_loom, args.ctrl_label)
    case = read_loom(args.case_loom, args.case_label)
    adata = ad.concat([ctrl, case], join="inner", index_unique=None, merge="same")
    adata.obs_names_make_unique()
    n_raw = adata.n_obs
    print(f"[01] cells after concatenation: {n_raw}")

    # -- 2. protein-coding genes --------------------------------------------
    adata = to_protein_coding(adata, args.annot)
    print(f"[01] genes after protein-coding filter: {adata.n_vars}")

    # -- 3. cell / gene filters ---------------------------------------------
    sc.pp.filter_cells(adata, min_genes=qc["min_genes"])
    sc.pp.filter_genes(adata, min_cells=qc["min_cells"])

    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                               log1p=False, inplace=True)

    sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                 jitter=0.4, multi_panel=True, show=False)
    plt.savefig(os.path.join(figdir, "01_qc_violin.pdf"), bbox_inches="tight")
    plt.close()

    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", color="condition", show=False)
    plt.savefig(os.path.join(figdir, "01_qc_scatter_mt.pdf"), bbox_inches="tight")
    plt.close()

    adata = adata[adata.obs["n_genes_by_counts"] <= qc["max_genes"], :].copy()
    adata = adata[adata.obs["pct_counts_mt"] <= qc["max_mt"], :].copy()
    print(f"[01] cells after QC: {adata.n_obs} ({adata.n_obs / n_raw:.1%} of raw)")

    # Display labels used by every downstream figure and table.
    adata.obs["genotype"] = pd.Categorical(
        np.where(adata.obs["condition"] == args.case_label, args.case_name, args.ctrl_name),
        categories=[args.ctrl_name, args.case_name])

    # -- 4. normalisation ----------------------------------------------------
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Write the full-gene normalised object before subsetting to HVGs: CellChat
    # needs a log-normalised matrix over the complete gene set, whereas the graph
    # and clustering below use highly variable genes only.
    adata.write_h5ad(os.path.join(args.outdir, f"{args.tag}_01_normalized.h5ad"))

    # -- 5. HVG, PCA, integration, graph, embedding, clustering --------------
    sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=args.hvg)
    work = adata[:, adata.var["highly_variable"]].copy()

    sc.tl.pca(work, svd_solver="arpack", n_comps=args.n_pcs, random_state=args.seed)
    if args.integration == "harmony":
        sce.pp.harmony_integrate(work, "condition", random_state=args.seed)
        rep = "X_pca_harmony"
    else:
        rep = "X_pca"

    sc.pp.neighbors(work, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs, use_rep=rep)
    sc.tl.umap(work, random_state=args.seed)
    sc.tl.louvain(work, resolution=args.resolution, key_added="louvain_1.0",
                  random_state=args.seed)

    adata.obs["louvain_1.0"] = work.obs["louvain_1.0"].values
    adata.obsm["X_umap"] = work.obsm["X_umap"]
    adata.obsm[rep] = work.obsm[rep]
    adata.uns["louvain_1.0"] = work.uns["louvain_1.0"]
    print(f"[01] clusters: {adata.obs['louvain_1.0'].nunique()} (graph on {rep})")

    # -- 6. markers ----------------------------------------------------------
    sc.tl.rank_genes_groups(adata, "louvain_1.0", method="wilcoxon",
                            use_raw=False, tie_correct=True,
                            key_added="wilcoxon_tie_1.0")
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=10, key="wilcoxon_tie_1.0",
                                    groupby="louvain_1.0", standard_scale="var",
                                    min_logfoldchange=1.2, cmap="vlag",
                                    dendrogram=True, show=False)
    plt.savefig(os.path.join(figdir, "01_cluster_markers.pdf"), bbox_inches="tight")
    plt.close()

    markers = ["Pecam1", "Cldn5", "Aqp4", "Slc1a2", "Cspg4", "Pdgfra", "Tmem119",
               "Gpr34", "Cd163", "Itgam", "Cd83", "Cd68", "S100a9", "Ptprc",
               "Gap43", "Mbp", "Mobp", "Dcn", "Col1a2", "Gad1", "Kit"]
    markers = [m for m in markers if m in adata.var_names]
    sc.pl.dotplot(adata, markers, groupby="louvain_1.0", dendrogram=True,
                  swap_axes=True, var_group_rotation=100, show=False)
    plt.savefig(os.path.join(figdir, "01_lineage_markers.pdf"), bbox_inches="tight")
    plt.close()

    # -- 7. annotation -------------------------------------------------------
    if os.path.exists(args.cluster_map):
        cmap = pd.read_csv(args.cluster_map, sep="\t", dtype=str)
        mapping = dict(zip(cmap.iloc[:, 0], cmap.iloc[:, 1]))
        missing = set(adata.obs["louvain_1.0"]) - set(mapping)
        if missing:
            sys.exit(f"[01] clusters without a cell_type in {args.cluster_map}: "
                     f"{sorted(missing)}")
        adata.obs["cell_type"] = [mapping[c] for c in adata.obs["louvain_1.0"]]
    else:
        print(f"[01] WARNING: {args.cluster_map} not found, cell_type = cluster id")
        adata.obs["cell_type"] = adata.obs["louvain_1.0"].astype(str)

    # -- 8. figures ----------------------------------------------------------
    sc.pl.umap(adata, color=["louvain_1.0"], legend_loc="on data",
               legend_fontsize="large", show=False)
    plt.savefig(os.path.join(figdir, "01_umap_clusters.pdf"), bbox_inches="tight")
    plt.close()

    sc.pl.umap(adata, color=["genotype", "cell_type"],
               legend_loc="on data", show=False)
    plt.savefig(os.path.join(figdir, "01_umap_celltype.pdf"), bbox_inches="tight")
    plt.close()

    # -- 9. composition ------------------------------------------------------
    # Ratio-to-expected: percentage of the cellularity-implied share. Bars are
    # therefore not constrained to sum to 100 and the 50% line is meaningful.
    A = pd.crosstab(adata.obs["cell_type"], adata.obs["genotype"])
    B = A / A.sum(axis=0)                      # per-library fraction
    counts = (B.T / B.sum(axis=1)).T * 100     # ratio to expected, %
    counts.to_csv(os.path.join(args.outdir, f"{args.tag}_composition.csv"))

    ax = counts.sort_values(args.case_name).plot.bar(
        stacked=True, color=["#0343df", "#f97306"], width=0.9, ylim=(0, 100),
        figsize=(15, 8), dpi=300)
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
    ax.axhline(y=50, color="k", linestyle="--")
    ax.set_ylabel("Cell proportion (%)", fontsize=16)
    plt.xticks(rotation=45, ha="right")
    plt.savefig(os.path.join(figdir, "01_composition.pdf"),
                dpi=300, bbox_inches="tight")
    plt.close()

    # -- 10. run record ------------------------------------------------------
    summary = (adata.obs.groupby(["genotype", "cell_type"], observed=True)
               .agg(n_cells=("cell_type", "size"),
                    median_genes=("n_genes_by_counts", "median"),
                    median_counts=("total_counts", "median"),
                    median_mt=("pct_counts_mt", "median")))
    summary.to_csv(os.path.join(args.outdir, f"{args.tag}_cell_summary.csv"))

    with open(os.path.join(args.outdir, f"{args.tag}_01_runinfo.txt"), "w") as fh:
        for k, v in qc.items():
            fh.write(f"{k}\t{v}\n")
        fh.write(f"integration\t{args.integration}\n")
        fh.write(f"hvg\t{args.hvg}\nn_pcs\t{args.n_pcs}\n")
        fh.write(f"n_neighbors\t{args.n_neighbors}\nresolution\t{args.resolution}\n")
        fh.write(f"labels\t{args.ctrl_label}={args.ctrl_name}, "
                 f"{args.case_label}={args.case_name}\n")
        fh.write(f"n_cells_raw\t{n_raw}\nn_cells_final\t{adata.n_obs}\n")
        fh.write(f"n_clusters\t{adata.obs['louvain_1.0'].nunique()}\n\n")
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            sc.logging.print_versions()
        fh.write(buf.getvalue())

    adata.write_h5ad(os.path.join(args.outdir, f"{args.tag}_01_annotated.h5ad"))
    print(f"[01] done -> {args.outdir}/{args.tag}_01_annotated.h5ad")


if __name__ == "__main__":
    main()
