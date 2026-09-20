#!/usr/bin/env python3
"""
02_create_master_table.py
==========================
Step 2: Create master table combining all information

Purpose
-------
Consolidates the filtered spacer-protospacer hits (output of
01_filter_true_matches.py) with per-prophage and per-genome annotations into a
single master table, by joining on prophage sequence ID (sseqid) and host
genome ID (host_sample):
  1. Host species, from the genome-to-organism mapping
  2. CheckV genome quality, per prophage
  3. BACPHLIP lifestyle prediction scores (virulent / temperate), per prophage
  4. Isolation source, per host genome

Input
-----
- spacer_hits_filtered.tsv
    Filtered true-match spacer hits (output of 01_filter_true_matches.py),
    including qseqid, sseqid, host_sample, votu, self_hit, etc.

- Cuadovircetes_phage_quality.tsv
    CheckV quality summary per prophage. Required columns: prophage_id,
    checkv_quality. Joined on sseqid == prophage_id.

- Cuadovircetes_vOTUs_bacphlip_result.tsv
    BACPHLIP lifestyle prediction output per phage sequence. Required columns:
    PahgeID [sic — matches a typo in the source file/pipeline, kept as-is so
    the join doesn't silently fail], Virulent, Temperate. Joined on
    sseqid == PahgeID.

- organism_name.csv
    Genome ID to species/organism name mapping. Required columns: genome_id,
    species. Joined on host_sample == genome_id.

- Isolation_source.txt
    Genome ID to isolation source mapping (comma-delimited despite the .txt
    extension). Required columns: genome_id, isolation_source. Joined on
    host_sample == genome_id.

Output
------
- master_spacer_table_complete.tsv
    spacer_hits_filtered.tsv with four additional columns appended:
    host_species, checkv_quality, virulent_score, temperate_score,
    isolation_source.

Notes
-----
- All joins are left joins via dict-based .map(), so hits with no match in a
  given lookup table (e.g., a prophage missing from the BACPHLIP results) will
  receive NaN in that column rather than being dropped. No explicit warning is
  currently printed for unmapped IDs — check for NaNs in the output if
  completeness matters downstream.
- 'PahgeID' is a typo present in the BACPHLIP output file itself (not
  introduced by this script); the column name here must match it exactly.

Usage
-----
    python 02_create_master_table.py
"""

import pandas as pd

# Read all files
spacer_hits = pd.read_csv('spacer_hits_filtered.tsv', sep='\t')
phage_quality = pd.read_csv('Cuadovircetes_phage_quality.tsv', sep='\t')
bacphlip = pd.read_csv('Cuadovircetes_vOTUs_bacphlip_result.tsv', sep='\t')
organism_map = pd.read_csv('organism_name.csv')
isolation_source = pd.read_csv('Isolation_source.txt', sep=',')

print("Files loaded successfully!")
print(f"spacer_hits: {len(spacer_hits)} rows")
print(f"phage_quality: {len(phage_quality)} rows")
print(f"bacphlip: {len(bacphlip)} rows")

# Map genome IDs to species
genome_to_species = dict(zip(organism_map['genome_id'], organism_map['species']))

# Add species to spacer hits
spacer_hits['host_species'] = spacer_hits['host_sample'].map(genome_to_species)

# Add phage quality (using prophage_id column)
phage_to_quality = dict(zip(phage_quality['prophage_id'], phage_quality['checkv_quality']))
spacer_hits['checkv_quality'] = spacer_hits['sseqid'].map(phage_to_quality)

# Add Bacphlip results (note the typo: PahgeID)
phage_to_virulent = dict(zip(bacphlip['PahgeID'], bacphlip['Virulent']))
phage_to_temperate = dict(zip(bacphlip['PahgeID'], bacphlip['Temperate']))
spacer_hits['virulent_score'] = spacer_hits['sseqid'].map(phage_to_virulent)
spacer_hits['temperate_score'] = spacer_hits['sseqid'].map(phage_to_temperate)

# Add isolation source
source_to_sample = dict(zip(isolation_source['genome_id'], isolation_source['isolation_source']))
spacer_hits['isolation_source'] = spacer_hits['host_sample'].map(source_to_sample)

# Save master table
spacer_hits.to_csv('master_spacer_table_complete.tsv', sep='\t', index=False)

print(f"\nMaster table created with {len(spacer_hits)} rows")
print("\nColumns:", spacer_hits.columns.tolist())
print("\nFirst few rows:")
print(spacer_hits[['qseqid', 'sseqid', 'votu', 'host_sample', 'host_species']].head())
