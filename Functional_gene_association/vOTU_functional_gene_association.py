#!/usr/bin/env python3
"""
Overall functional gene-category composition by host-range class,
controlling for phage taxonomy.

Inputs:
    combined_CDS_476_vOTUs.tsv   - per-gene annotation table (Pharokka output)
    vOTU_host_range_summary.tsv  - per-vOTU CRISPR-based host-range data
    vOTU_taxonomy.csv           - per-vOTU taxonomy (vContact3 deepest taxonomy)

Outputs:
    summary_category_by_hostrange.csv          - mean/SEM gene count per category x host-range class
    summary_category_by_hostrange_taxonomy.csv - same, stratified by taxonomy
    category_lm_pvalues.csv                     - per-category linear model p-value for
                                                   host_range_class effect, controlling for taxonomy
    category_hostrange_barplot.png/.pdf         - grouped bar plot
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.formula.api as smf
import warnings

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# 0. FILE PATHS -- edit these if your files are named/located differently
# ---------------------------------------------------------------------------
CDS_FILE = "combined_CDS_476_vOTUs.tsv"
HOSTRANGE_FILE = "vOTU_host_range_summary.tsv"
TAXONOMY_FILE = "vOTU_taxonomy.csv"

OUT_SUMMARY = "summary_category_by_hostrange.csv"
OUT_SUMMARY_TAX = "summary_category_by_hostrange_taxonomy.csv"
OUT_PVALUES = "category_lm_pvalues.csv"
OUT_PLOT_PNG = "category_hostrange_barplot.png"
OUT_PLOT_PDF = "category_hostrange_barplot.pdf"
OUT_PLOT_SVG = "category_hostrange_barplot.svg"

# ---------------------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------------------
cds = pd.read_csv(CDS_FILE, sep="\t", low_memory=False)
hostrange = pd.read_csv(HOSTRANGE_FILE, sep="\t")
taxonomy = pd.read_csv(TAXONOMY_FILE)

print(f"[CDS] {cds.shape[0]} genes, {cds['contig'].nunique()} unique vOTUs")
print(f"[Host range] {hostrange.shape[0]} vOTUs")
print(f"[Taxonomy] {taxonomy.shape[0]} vOTUs")

# Sanity check on required columns
for col in ["contig", "category"]:
    assert col in cds.columns, f"Missing expected column '{col}' in {CDS_FILE}"
for col in ["vOTU", "species_count_robust"]:
    assert col in hostrange.columns, f"Missing expected column '{col}' in {HOSTRANGE_FILE}"
for col in ["Genome", "Deepest_taxonomy"]:
    assert col in taxonomy.columns, f"Missing expected column '{col}' in {TAXONOMY_FILE}"

# ---------------------------------------------------------------------------
# 2. DERIVE HOST-RANGE CLASS  (matches Fig. 5a thresholds)
#    1 species        -> specialist
#    2-3 species       -> moderate
#    >=4 species       -> generalist
# ---------------------------------------------------------------------------
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

# Restrict to vOTUs with a non-zero robust species count (matches paper's n=476 set)
hostrange_classified = hostrange.dropna(subset=["host_range_class"]).copy()
print(f"[Host range] {hostrange_classified.shape[0]} vOTUs with valid host-range class")
print(hostrange_classified["host_range_class"].value_counts())

# ---------------------------------------------------------------------------
# 3. PER-VOTU GENE COUNTS PER FUNCTIONAL CATEGORY
# ---------------------------------------------------------------------------
gene_counts = (
    cds.groupby(["contig", "category"])
    .size()
    .reset_index(name="gene_count")
)

# Pivot to wide: rows = vOTU, columns = category, values = gene_count
wide = gene_counts.pivot_table(
    index="contig", columns="category", values="gene_count", fill_value=0
).reset_index()

# ---------------------------------------------------------------------------
# 4. MERGE HOST RANGE + TAXONOMY
# ---------------------------------------------------------------------------
merged = wide.merge(
    hostrange_classified[["vOTU", "species_count_robust", "host_range_class"]],
    left_on="contig", right_on="vOTU", how="inner",
)
merged = merged.merge(
    taxonomy[["Genome", "Deepest_taxonomy"]],
    left_on="contig", right_on="Genome", how="left",
)

n_missing_tax = merged["Deepest_taxonomy"].isna().sum()
print(f"[Merge] {merged.shape[0]} vOTUs with host-range class; "
      f"{n_missing_tax} missing taxonomy assignment")

merged["Deepest_taxonomy"] = merged["Deepest_taxonomy"].fillna("Unassigned")

category_cols = [c for c in wide.columns if c != "contig"]

# ---------------------------------------------------------------------------
# 5. SUMMARY 1: mean +/- SEM gene count per category, by host-range class
#    (overall picture, no taxonomy stratification)
# ---------------------------------------------------------------------------
long_df = merged.melt(
    id_vars=["contig", "host_range_class", "Deepest_taxonomy"],
    value_vars=category_cols,
    var_name="category",
    value_name="gene_count",
)

summary = (
    long_df.groupby(["category", "host_range_class"])["gene_count"]
    .agg(mean="mean", sem=stats.sem, n="count")
    .reset_index()
)
summary.to_csv(OUT_SUMMARY, index=False)
print(f"\n[Saved] {OUT_SUMMARY}")

# ---------------------------------------------------------------------------
# 6. SUMMARY 2: same, stratified by taxonomy (only taxa with vOTUs spanning
#    >=2 host-range classes, so the comparison is meaningful within-lineage)
# ---------------------------------------------------------------------------
tax_class_counts = (
    long_df.drop_duplicates(subset=["contig"])
    .groupby("Deepest_taxonomy")["host_range_class"]
    .nunique()
)
taxa_with_multiple_classes = tax_class_counts[tax_class_counts >= 2].index.tolist()
print(f"\n[Taxonomy control] {len(taxa_with_multiple_classes)} taxa have vOTUs "
      f"spanning >=2 host-range classes")

summary_tax = (
    long_df[long_df["Deepest_taxonomy"].isin(taxa_with_multiple_classes)]
    .groupby(["Deepest_taxonomy", "category", "host_range_class"])["gene_count"]
    .agg(mean="mean", sem=stats.sem, n="count")
    .reset_index()
)
summary_tax.to_csv(OUT_SUMMARY_TAX, index=False)
print(f"[Saved] {OUT_SUMMARY_TAX}")

# ---------------------------------------------------------------------------
# 7. PER-CATEGORY LINEAR MODEL: gene_count ~ host_range_class + taxonomy
#    Extracts an F-test p-value for the host_range_class term after
#    controlling for taxonomy (a formal version of the "control" step)
# ---------------------------------------------------------------------------
pval_rows = []
for cat in category_cols:
    sub = merged[["contig", cat, "host_range_class", "Deepest_taxonomy"]].copy()
    sub = sub.rename(columns={cat: "gene_count"})

    # Only keep taxa with >=3 vOTUs total, and require >=2 host-range classes
    # represented overall, else the model is unidentifiable/uninformative
    if sub["host_range_class"].nunique() < 2:
        continue

    tax_counts = sub["Deepest_taxonomy"].value_counts()
    valid_taxa = tax_counts[tax_counts >= 3].index
    sub = sub[sub["Deepest_taxonomy"].isin(valid_taxa)]
    if sub["Deepest_taxonomy"].nunique() < 2 or sub.shape[0] < 10:
        continue

    try:
        model_full = smf.ols(
            "gene_count ~ C(host_range_class) + C(Deepest_taxonomy)", data=sub
        ).fit()
        model_reduced = smf.ols(
            "gene_count ~ C(Deepest_taxonomy)", data=sub
        ).fit()
        # F-test comparing full vs. reduced model (isolates host_range_class effect)
        from statsmodels.stats.anova import anova_lm
        anova_res = anova_lm(model_reduced, model_full)
        p_val = anova_res["Pr(>F)"].iloc[1]
        pval_rows.append({"category": cat, "n_vOTUs": sub.shape[0],
                           "n_taxa": sub["Deepest_taxonomy"].nunique(),
                           "p_host_range_controlling_taxonomy": p_val})
    except Exception as e:
        print(f"  [skip] {cat}: {e}")

pval_df = pd.DataFrame(pval_rows).sort_values("p_host_range_controlling_taxonomy")
pval_df.to_csv(OUT_PVALUES, index=False)
print(f"\n[Saved] {OUT_PVALUES}")
print(pval_df.to_string(index=False))

# ---------------------------------------------------------------------------
# 8. GROUPED BAR PLOT (overall picture, all categories, by host-range class)
# ---------------------------------------------------------------------------
order = ["specialist", "moderate", "generalist"]
palette = {"specialist": "#4C72B0", "moderate": "#55A868", "generalist": "#C44E52"}

plt.figure(figsize=(12, 6))
ax = sns.barplot(
    data=long_df,
    x="category", y="gene_count", hue="host_range_class",
    hue_order=order, palette=palette,
    estimator=np.mean, errorbar="se",
)
ax.set_xlabel("Functional gene category")
ax.set_ylabel("Mean gene count per vOTU")
ax.set_title("Gene-category composition by host-range class")
plt.xticks(rotation=40, ha="right")
plt.legend(title="Host-range class")
plt.tight_layout()
plt.savefig(OUT_PLOT_PNG, dpi=300)
plt.savefig(OUT_PLOT_PDF)
plt.savefig(OUT_PLOT_SVG, format="svg")
print(f"\n[Saved] {OUT_PLOT_PNG}, {OUT_PLOT_PDF}, {OUT_PLOT_SVG}")

print("\nDone.")
