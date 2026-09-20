#!/usr/bin/env python3
"""
Focused association analysis: does phage taxonomy predict lysis cassette
composition across vOTUs?

Tests taxonomy at two resolutions:
  1. Deepest_taxonomy (order/genus-level)
  2. Family_subfamily, nested within each Deepest_taxonomy group that
     shows internal variation - to see whether finer taxonomic resolution
     resolves the "unexplained" variation left at the coarser level.

Also breaks down association by individual lysis gene (not just the
complete Holin+Endolysin+Spanin cassette), since different gene types may
track taxonomy differently.

Requires: pandas, scipy, matplotlib
    pip install pandas scipy matplotlib

Usage:
    python lysis_taxonomy_association.py \
        <presence_absence_tsv> <vOTU_taxonomy_csv> [output_prefix]

Inputs:
    <presence_absence_tsv> - vOTU x gene-type 0/1 matrix
                              (build_lysis_matrices.py output)
    <vOTU_taxonomy_csv>   - columns: Genome, Deepest_taxonomy,
                              Family_subfamily

'Complete' lysis cassette is defined as Holin + Endolysin + Rz-like spanin
all present.

Outputs:
    <output_prefix>_merged_data.tsv         - merged per-vOTU table
    <output_prefix>_gene_by_taxonomy.tsv    - % presence of each gene type
                                               (+ complete cassette) per
                                               Deepest_taxonomy group
    <output_prefix>_association_tests.tsv   - chi-square results: each
                                               gene type / complete
                                               cassette vs Deepest_taxonomy,
                                               plus Family_subfamily tests
                                               nested within variable taxa
    <output_prefix>_heatmap.png/.pdf/.svg   - heatmap of % gene presence
                                               by Deepest_taxonomy
    <output_prefix>_subfamily_<taxon>.png/.pdf/.svg
                                             - bar chart of % complete
                                               cassette by Family_subfamily,
                                               for each taxon with internal
                                               variation
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2_contingency


GENE_COLS_DEFAULT = [
    "holin", "endolysin", "Rz-like spanin",
    "lysis determining protein", "amidase",
]


def collapse_rare(ct, min_n=5):
    small = ct.sum(axis=1) < min_n
    if small.sum() > 1:
        other = ct[small].sum()
        ct = ct[~small].copy()
        ct.loc["Other"] = other
    return ct


def run_chi_square(df, predictor, target, min_n=5):
    sub = df[[predictor, target]].dropna()
    if sub[predictor].nunique() < 2 or sub[target].nunique() < 2:
        return None
    ct = pd.crosstab(sub[predictor], sub[target])
    ct = collapse_rare(ct, min_n=min_n)
    if ct.shape[0] < 2 or ct.shape[1] < 2:
        return None
    chi2, p, dof, expected = chi2_contingency(ct)
    n = ct.values.sum()
    k = min(ct.shape)
    cramers_v = (chi2 / (n * (k - 1))) ** 0.5
    return {"n": n, "chi2": chi2, "dof": dof, "p": p, "cramers_v": cramers_v}


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    pa_path, tax_path = sys.argv[1], sys.argv[2]
    out_prefix = sys.argv[3] if len(sys.argv) > 3 else "lysis_taxonomy"

    pa = pd.read_csv(pa_path, sep="\t", index_col=0)
    tax = pd.read_csv(tax_path)

    gene_cols = [c for c in GENE_COLS_DEFAULT if c in pa.columns]

    df = pa.copy()
    df.index.name = None
    df["vOTU"] = df.index
    df = df.merge(
        tax[["Genome", "Deepest_taxonomy", "Family_subfamily"]],
        left_on="vOTU", right_on="Genome", how="left",
    )
    df["complete_cassette"] = (
        (df["holin"] == 1) & (df["endolysin"] == 1) & (df["Rz-like spanin"] == 1)
    ).astype(int)

    targets = gene_cols + ["complete_cassette"]
    df.to_csv(f"{out_prefix}_merged_data.tsv", sep="\t", index=False)

    print(f"Merged table: {df.shape[0]} vOTUs "
          f"({df['Deepest_taxonomy'].isna().sum()} missing taxonomy)")
    print()

    # --- % presence of each gene / complete cassette, per Deepest_taxonomy ---
    tax_counts = df.groupby("Deepest_taxonomy").size().rename("n")
    pct_table = df.groupby("Deepest_taxonomy")[targets].mean() * 100
    pct_table = pct_table.round(1)
    pct_table.insert(0, "n", tax_counts)
    pct_table = pct_table.sort_values("n", ascending=False)
    pct_table.to_csv(f"{out_prefix}_gene_by_taxonomy.tsv", sep="\t")

    print("=== %% presence by Deepest_taxonomy (n>=5 shown) ===")
    print(pct_table[pct_table["n"] >= 5].to_string())
    print()

    # --- Chi-square: each gene / complete cassette vs Deepest_taxonomy ---
    results = []
    for target in targets:
        res = run_chi_square(df, "Deepest_taxonomy", target)
        if res is None:
            continue
        results.append({"level": "Deepest_taxonomy", "stratum": "(all)",
                         "target": target, **res})

    print("=== Chi-square: gene/cassette presence vs Deepest_taxonomy ===")
    for r in results:
        print(f"{r['target']:30s} n={r['n']:4d}  chi2={r['chi2']:7.2f}  "
              f"dof={r['dof']:3d}  p={r['p']:.3g}  Cramer's V={r['cramers_v']:.3f}")
    print()

    # --- Nested Family_subfamily tests, within taxa showing internal
    #     variation in complete_cassette (i.e. neither ~0% nor ~100%) ---
    var_check = df.groupby("Deepest_taxonomy")["complete_cassette"].agg(["mean", "count"])
    variable_taxa = var_check[
        (var_check["mean"] > 0) & (var_check["mean"] < 1) & (var_check["count"] >= 10)
    ].index.tolist()

    print(f"Deepest_taxonomy groups with internal variation (n>=10): {variable_taxa}")
    print()

    for taxon in variable_taxa:
        sub_df = df[df["Deepest_taxonomy"] == taxon]
        res = run_chi_square(sub_df, "Family_subfamily", "complete_cassette")
        if res is None:
            continue
        results.append({"level": "Family_subfamily", "stratum": taxon,
                         "target": "complete_cassette", **res})
        print(f"=== Family_subfamily vs complete_cassette, within {taxon} ===")
        print(f"n={res['n']}, chi2={res['chi2']:.2f}, dof={res['dof']}, "
              f"p={res['p']:.4g}, Cramer's V={res['cramers_v']:.3f}")
        breakdown = sub_df.groupby("Family_subfamily")["complete_cassette"].agg(["mean", "count"])
        breakdown["mean"] = (breakdown["mean"] * 100).round(1)
        print(breakdown.sort_values("count", ascending=False).to_string())
        print()

        # bar chart for this taxon's subfamily breakdown
        bd = breakdown[breakdown["count"] >= 3].sort_values("mean", ascending=False)
        if not bd.empty:
            fig, ax = plt.subplots(figsize=(max(6, len(bd) * 1.3), 5))
            ax.bar(bd.index.astype(str), bd["mean"], color="#3b3b3b", width=0.6)
            for i, (val, n) in enumerate(zip(bd["mean"], bd["count"])):
                ax.text(i, val, f"{val:.0f}%\n(n={n})", ha="center", va="bottom", fontsize=8)
            ax.set_ylabel("% with complete lysis cassette")
            ax.set_title(f"Complete cassette by Family/Subfamily\nwithin {taxon}")
            ax.set_ylim(0, max(bd["mean"]) * 1.2 + 5)
            ax.spines[["top", "right"]].set_visible(False)
            plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
            fig.tight_layout()
            safe = taxon.replace(" ", "_").replace("|", "_")
            fig.savefig(f"{out_prefix}_subfamily_{safe}.png", dpi=300, bbox_inches="tight")
            fig.savefig(f"{out_prefix}_subfamily_{safe}.pdf", bbox_inches="tight")
            fig.savefig(f"{out_prefix}_subfamily_{safe}.svg", bbox_inches="tight")
            plt.close(fig)

    results_df = pd.DataFrame(results).sort_values("p")
    results_df.to_csv(f"{out_prefix}_association_tests.tsv", sep="\t", index=False)

    # --- Heatmap: % gene presence by Deepest_taxonomy (taxa with n>=5) ---
    heat_df = pct_table[pct_table["n"] >= 5][targets]
    heat_df = heat_df.sort_values("complete_cassette", ascending=False)

    fig, ax = plt.subplots(figsize=(8, max(4, len(heat_df) * 0.4)))
    im = ax.imshow(heat_df.values, cmap="Greys", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(heat_df.columns)))
    ax.set_xticklabels(heat_df.columns, rotation=30, ha="right")
    ax.set_yticks(range(len(heat_df.index)))
    ax.set_yticklabels(heat_df.index)
    for i in range(heat_df.shape[0]):
        for j in range(heat_df.shape[1]):
            val = heat_df.values[i, j]
            color = "white" if val > 60 else "black"
            ax.text(j, i, f"{val:.0f}", ha="center", va="center", fontsize=8, color=color)
    ax.set_title("% gene presence by phage taxonomy")
    fig.colorbar(im, ax=ax, label="% of vOTUs", shrink=0.6)
    fig.tight_layout()
    fig.savefig(f"{out_prefix}_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_prefix}_heatmap.pdf", bbox_inches="tight")
    fig.savefig(f"{out_prefix}_heatmap.svg", bbox_inches="tight")
    plt.close(fig)

    print("=== Summary of all chi-square tests (sorted by p-value) ===")
    print(results_df.to_string(index=False))
    print()
    print(f"Wrote merged data -> {out_prefix}_merged_data.tsv")
    print(f"Wrote %% table -> {out_prefix}_gene_by_taxonomy.tsv")
    print(f"Wrote test summary -> {out_prefix}_association_tests.tsv")
    print(f"Wrote heatmap -> {out_prefix}_heatmap.png / .pdf / .svg")
    print(f"Wrote subfamily bar charts -> {out_prefix}_subfamily_<taxon>.png / .pdf / .svg")


if __name__ == "__main__":
    main()
