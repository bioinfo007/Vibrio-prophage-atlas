#!/usr/bin/env python3
"""
Per-vOTU tail fiber gene copy number vs. host-range class,
controlling for phage taxonomy.

Pipeline:
    1. Parse the CD-HIT .clstr file to map every tail-category CDS
       (n=3,706) to its cluster number (n=513 clusters, -c 0.5).
    2. Use the confidently PhANNs-classified "tail fiber" cluster list
       (n=130 clusters, from phaggns_to_clusters.csv) to flag every member
       gene in those clusters as a tail-fiber gene.
    3. Join back to combined_CDS_476_vOTUs.tsv to recover each gene's
       source vOTU (contig).
    4. Compute per-vOTU tail-fiber gene count (and total tail gene count
       for reference / proportion).
    5. Merge with host-range class and taxonomy; test the tail-fiber
       count against host-range class, controlling for taxonomy
       (same nested-model F-test approach as the category-level analysis).

Inputs:
    combined_CDS_476_vOTUs.tsv   - per-gene annotation table (gene -> contig, category)
    vOTU_host_range_summary.tsv  - per-vOTU CRISPR-based host-range data
    vOTU_taxonomy.csv           - per-vOTU taxonomy
    tail_genes_50_faa.clstr      - CD-HIT cluster membership (raw .clstr output)
    phaggns_to_clusters.csv      - confidently-classified "tail fiber" clusters
                                    (Gene = representative, Cluster = cluster #,
                                     Cluster_Size = # members)

Outputs:
    gene_tailfiber_flags.csv               - every tail CDS with cluster # and tail-fiber flag
    vOTU_tailfiber_counts.csv              - per-vOTU tail-fiber and total-tail gene counts
    summary_tailfiber_by_hostrange.csv     - mean/SEM tail-fiber count per host-range class
    tailfiber_lm_pvalue.csv                - taxonomy-controlled F-test p-value
    tailfiber_hostrange_barplot.png/.pdf/.svg
"""

import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
import warnings

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# 0. FILE PATHS -- edit if needed
# ---------------------------------------------------------------------------
CDS_FILE = "combined_CDS_476_vOTUs.tsv"
HOSTRANGE_FILE = "vOTU_host_range_summary.tsv"
TAXONOMY_FILE = "vOTU_taxonomy.csv"
CLSTR_FILE = "tail_genes_50_faa.clstr"
TAILFIBER_CLUSTERS_FILE = "phaggns_to_clusters.csv"

OUT_GENE_FLAGS = "gene_tailfiber_flags.csv"
OUT_VOTU_COUNTS = "vOTU_tailfiber_counts.csv"
OUT_SUMMARY = "summary_tailfiber_by_hostrange.csv"
OUT_PVALUE = "tailfiber_lm_pvalue.csv"
OUT_PLOT_PNG = "tailfiber_hostrange_barplot.png"
OUT_PLOT_PDF = "tailfiber_hostrange_barplot.pdf"
OUT_PLOT_SVG = "tailfiber_hostrange_barplot.svg"

# ---------------------------------------------------------------------------
# 1. PARSE THE .clstr FILE -> gene_id -> cluster_number
# ---------------------------------------------------------------------------
gene_cluster = {}
current_cluster = None
member_re = re.compile(r">([^.]+)\.\.\.")

with open(CLSTR_FILE) as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">Cluster"):
            current_cluster = int(line.split()[-1])
        else:
            m = member_re.search(line)
            if m:
                gene_id = m.group(1)
                gene_cluster[gene_id] = current_cluster

clstr_df = pd.DataFrame(
    {"gene": list(gene_cluster.keys()), "cluster": list(gene_cluster.values())}
)
print(f"[.clstr] Parsed {clstr_df.shape[0]} genes across "
      f"{clstr_df['cluster'].nunique()} clusters")

# ---------------------------------------------------------------------------
# 2. IDENTIFY TAIL-FIBER CLUSTERS (confidently classified by PhANNs)
# ---------------------------------------------------------------------------
tf_clusters_df = pd.read_csv(TAILFIBER_CLUSTERS_FILE)
tailfiber_cluster_ids = set(tf_clusters_df["Cluster"].astype(int))
print(f"[Tail-fiber clusters] {len(tailfiber_cluster_ids)} clusters "
      f"confidently classified as tail fiber")

clstr_df["is_tail_fiber"] = clstr_df["cluster"].isin(tailfiber_cluster_ids)
print(f"[Gene-level] {clstr_df['is_tail_fiber'].sum()} of {clstr_df.shape[0]} "
      f"tail CDS flagged as tail-fiber genes")

clstr_df.to_csv(OUT_GENE_FLAGS, index=False)
print(f"[Saved] {OUT_GENE_FLAGS}")

# ---------------------------------------------------------------------------
# 3. JOIN BACK TO combined_CDS_476_vOTUs.tsv TO GET SOURCE vOTU (contig)
# ---------------------------------------------------------------------------
cds = pd.read_csv(CDS_FILE, sep="\t", low_memory=False)
tail_cds = cds[cds["category"] == "tail"][["gene", "contig"]].copy()
print(f"\n[CDS] {tail_cds.shape[0]} tail-category CDS in {CDS_FILE}")

flagged = tail_cds.merge(clstr_df[["gene", "cluster", "is_tail_fiber"]],
                          on="gene", how="left")
n_unmatched = flagged["cluster"].isna().sum()
if n_unmatched > 0:
    print(f"  [!] {n_unmatched} tail CDS from {CDS_FILE} not found in the "
          f".clstr file -- check gene ID formatting consistency")
flagged["is_tail_fiber"] = flagged["is_tail_fiber"].fillna(False)

# ---------------------------------------------------------------------------
# 4. PER-vOTU TAIL-FIBER AND TOTAL-TAIL GENE COUNTS
# ---------------------------------------------------------------------------
votu_counts = (
    flagged.groupby("contig")
    .agg(
        total_tail_genes=("gene", "count"),
        tail_fiber_genes=("is_tail_fiber", "sum"),
    )
    .reset_index()
)
votu_counts["pct_tail_fiber"] = (
    votu_counts["tail_fiber_genes"] / votu_counts["total_tail_genes"] * 100
)
votu_counts.to_csv(OUT_VOTU_COUNTS, index=False)
print(f"[Saved] {OUT_VOTU_COUNTS}")

# ---------------------------------------------------------------------------
# 5. MERGE WITH HOST RANGE + TAXONOMY
# ---------------------------------------------------------------------------
hostrange = pd.read_csv(HOSTRANGE_FILE, sep="\t")
taxonomy = pd.read_csv(TAXONOMY_FILE)


def classify_host_range(n):
    if pd.isna(n):
        return np.nan
    if n <= 1:
        return "specialist"
    elif n <= 3:
        return "moderate"
    else:
        return "generalist"


hostrange["host_range_class"] = hostrange["species_count_robust"].apply(classify_host_range)
hostrange_classified = hostrange.dropna(subset=["host_range_class"]).copy()

merged = votu_counts.merge(
    hostrange_classified[["vOTU", "species_count_robust", "host_range_class"]],
    left_on="contig", right_on="vOTU", how="inner",
)
merged = merged.merge(
    taxonomy[["Genome", "Deepest_taxonomy"]],
    left_on="contig", right_on="Genome", how="left",
)
merged["Deepest_taxonomy"] = merged["Deepest_taxonomy"].fillna("Unassigned")

print(f"\n[Merged] {merged.shape[0]} vOTUs with tail-fiber counts + "
      f"host-range class + taxonomy")

# ---------------------------------------------------------------------------
# 6. SUMMARY: mean +/- SEM tail-fiber count per host-range class
# ---------------------------------------------------------------------------
summary = (
    merged.groupby("host_range_class")
    .agg(
        mean_tail_fiber=("tail_fiber_genes", "mean"),
        sem_tail_fiber=("tail_fiber_genes", stats.sem),
        mean_pct_tail_fiber=("pct_tail_fiber", "mean"),
        sem_pct_tail_fiber=("pct_tail_fiber", stats.sem),
        n=("contig", "count"),
    )
    .reindex(["specialist", "moderate", "generalist"])
    .reset_index()
)
summary.to_csv(OUT_SUMMARY, index=False)
print(f"\n[Saved] {OUT_SUMMARY}")
print(summary.round(3).to_string(index=False))

# ---------------------------------------------------------------------------
# 7. TAXONOMY-CONTROLLED LINEAR MODEL: tail_fiber_genes ~ host_range_class + taxonomy
# ---------------------------------------------------------------------------
sub = merged[["contig", "tail_fiber_genes", "host_range_class", "Deepest_taxonomy"]].copy()
tax_counts = sub["Deepest_taxonomy"].value_counts()
valid_taxa = tax_counts[tax_counts >= 3].index
sub = sub[sub["Deepest_taxonomy"].isin(valid_taxa)]

pval_row = {}
if sub["host_range_class"].nunique() >= 2 and sub["Deepest_taxonomy"].nunique() >= 2:
    model_full = smf.ols(
        "tail_fiber_genes ~ C(host_range_class) + C(Deepest_taxonomy)", data=sub
    ).fit()
    model_reduced = smf.ols(
        "tail_fiber_genes ~ C(Deepest_taxonomy)", data=sub
    ).fit()
    anova_res = anova_lm(model_reduced, model_full)
    p_val = anova_res["Pr(>F)"].iloc[1]
    pval_row = {"n_vOTUs": sub.shape[0], "n_taxa": sub["Deepest_taxonomy"].nunique(),
                "p_host_range_controlling_taxonomy": p_val}
    print(f"\n[Model] tail_fiber_genes ~ host_range_class + taxonomy "
          f"(n={sub.shape[0]}, {sub['Deepest_taxonomy'].nunique()} taxa)")
    print(f"  P (host-range effect, controlling for taxonomy) = {p_val:.4g}")
else:
    print("\n[Model] Insufficient groups to fit taxonomy-controlled model")

pd.DataFrame([pval_row]).to_csv(OUT_PVALUE, index=False)
print(f"[Saved] {OUT_PVALUE}")

# ---------------------------------------------------------------------------
# 8. BAR PLOT: mean tail-fiber gene count per vOTU, by host-range class
# ---------------------------------------------------------------------------
order = ["specialist", "moderate", "generalist"]
palette = {"specialist": "#4C72B0", "moderate": "#55A868", "generalist": "#C44E52"}

fig, axes = plt.subplots(1, 2, figsize=(10, 5))

sns.barplot(
    data=merged, x="host_range_class", y="tail_fiber_genes",
    order=order, palette=palette, estimator=np.mean, errorbar="se", ax=axes[0],
)
axes[0].set_xlabel("Host-range class")
axes[0].set_ylabel("Mean tail-fiber genes per vOTU")
axes[0].set_title("Tail-fiber gene copy number")

sns.barplot(
    data=merged, x="host_range_class", y="pct_tail_fiber",
    order=order, palette=palette, estimator=np.mean, errorbar="se", ax=axes[1],
)
axes[1].set_xlabel("Host-range class")
axes[1].set_ylabel("Tail-fiber genes (% of tail genes)")
axes[1].set_title("Tail-fiber proportion of tail module")

plt.tight_layout()
plt.savefig(OUT_PLOT_PNG, dpi=300)
plt.savefig(OUT_PLOT_PDF)
plt.savefig(OUT_PLOT_SVG, format="svg")
print(f"\n[Saved] {OUT_PLOT_PNG}, {OUT_PLOT_PDF}, {OUT_PLOT_SVG}")

print("\nDone.")
