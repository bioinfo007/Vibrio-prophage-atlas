#!/usr/bin/env python3
"""
01_filter_true_matches.py
==========================
Step 1: Filter spacer matches with controlled identity, mismatch, and length to obtain the true matches

Purpose
-------
Filters raw BLAST-style spacer-vs-prophage hits down to high-confidence, true
spacer-protospacer matches (near-perfect identity, no gaps, near-full-length
alignment), then:
  1. Extracts the CRISPR spacer's host genome (biosample ID) from the query ID
  2. Extracts the targeted prophage's host genome from the subject ID
  3. Flags self-targeting hits (spacer and prophage originate from the same genome)
  4. Maps each targeted prophage sequence to its vOTU cluster
  5. Summarizes vOTU-level host range based on cross-genome (non-self) hits

Input
-----
- spacer_hits_sample_input.tsv
    Tabular BLAST/outfmt6-style output of spacer-vs-prophage searches, with
    columns: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend,
    sstart, send, evalue, bitscore, qlen, slen (no header row).
    qseqid format: <host_biosample>_<contig>_CRISPR_<n>_spacer_<n>
    sseqid format: <host_biosample>_p<prophage_number>

- Cuadovircetes_ani_clusters.tsv
    vOTU cluster membership table (aniclust-style output), with columns:
    votu_rep, members (comma-separated list of member sequence IDs).

Output
------
- spacer_hits_filtered.tsv
    Filtered true-match hits, annotated with host_sample, prophage_host_genome,
    self_hit flag, and vOTU assignment.

- votu_host_range.tsv
    Per-vOTU count of distinct host genomes hit by spacers, excluding self-hits
    (i.e., historical cross-genome infection evidence), sorted by descending
    host range.

Filtering criteria for a "true" spacer-protospacer match
----------------------------------------------------------
- gapopen == 0        (no gaps in the alignment)
- mismatch <= 0        (perfect identity)
- length >= qlen - 1   (alignment spans effectively the full spacer length,
                         allowing at most a 1 bp truncation)

Notes
-----
- host_sample and prophage_host_genome are both stripped to the bare biosample
  ID (no contig/prophage suffix) so they can be directly compared for
  self-targeting detection.
- Update the input file path (spacer_hits_sample_input.tsv) and vOTU cluster
  file path (Cuadovircetes_ani_clusters.tsv) below as needed for your run.

Usage
-----
    python 01_filter_true_matches.py
"""

import pandas as pd
import re

cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
        "qstart","qend","sstart","send","evalue","bitscore","qlen","slen"]

df = pd.read_csv("spacer_hits_sample_input.tsv", sep="\t", names=cols)  #Change the input file path here

# --- 1. Filter to true spacer-protospacer matches ---
hits = df[
    (df.gapopen == 0) &
    (df.mismatch <= 0) &
    (df.length >= df.qlen - 1)
].copy()

# --- 2. Extract TRUE host genome (biosample only, no contig) from qseqid ---
# e.g. SAMN02435881_1_CRISPR_1_spacer_1 -> SAMN02435881
# (the "_1" after the biosample ID is a contig number, not part of the genome ID,
#  and prophage sequence IDs don't carry a contig suffix, so it must be stripped
#  here for host_sample and prophage_host_genome to be comparable)
hits["host_sample"] = hits.qseqid.str.extract(r"^([A-Z]+\d+)_")

# sanity check: make sure every qseqid matched the pattern
n_missing = hits["host_sample"].isna().sum()
if n_missing:
    print(f"WARNING: {n_missing} qseqid values did not match the biosample ID pattern")
    print(hits.loc[hits["host_sample"].isna(), "qseqid"].unique()[:10])

# --- 3. Extract prophage host genome from sseqid ---
# e.g. SAMN13704165_p2 -> SAMN13704165
hits["prophage_host_genome"] = hits.sseqid.str.extract(r"^(.+?)_p\d+$")

# --- 4. Flag self-targeting ---
hits["self_hit"] = hits.host_sample == hits.prophage_host_genome

# --- 5. Map phage sequences to vOTU clusters ---
# adjust parsing to match your actual aniclust output format
votu_map = pd.read_csv("Cuadovircetes_ani_clusters.tsv", sep="\t", names=["votu_rep","members"])
seq2votu = votu_map.assign(member=votu_map.members.str.split(",")).explode("member")
seq2votu = seq2votu[["member","votu_rep"]].rename(columns={"member":"sseqid","votu_rep":"votu"})
hits = hits.merge(seq2votu, on="sseqid", how="left")

# check for unmapped sequences (sanity check)
unmapped = hits[hits.votu.isna()].sseqid.unique()
if len(unmapped):
    print(f"WARNING: {len(unmapped)} phage sequences had no vOTU mapping, e.g.: {unmapped[:5]}")

hits.to_csv("spacer_hits_filtered.tsv", sep="\t", index=False)

# --- 6. vOTU-level host range (cross-genome only) ---
host_range = (
    hits[~hits.self_hit]
    .groupby("votu")["host_sample"]
    .nunique()
    .reset_index(name="n_distinct_hosts")
    .sort_values("n_distinct_hosts", ascending=False)
)
host_range.to_csv("votu_host_range.tsv", sep="\t", index=False)

# --- 7. Self-targeting summary ---
self_targeting = hits[hits.self_hit]
print(f"Self-targeting spacer-prophage pairs: {len(self_targeting)}")
print(f"Distinct genomes showing self-targeting: {self_targeting.host_sample.nunique()}")
