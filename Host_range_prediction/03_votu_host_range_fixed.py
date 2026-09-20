#!/usr/bin/env python3
"""
Step 3: Calculate host range at vOTU level.

===============================================================================
HOST-RANGE PREDICTION APPROACH
===============================================================================
Input: master_spacer_table_complete.tsv -- CRISPR spacer to
protospacer BLASTN hits (No mismatch tolerated, Full-length spacer
alignment required), annotated with vOTU cluster, host sample, host
species, and self-hit status.

Filtering applied here:
  1. Self-targeting hits excluded (self_hit == True; spacer's own resident
     genome matching its own prophage -- not evidence of cross-strain
     infection).
  2. A species is only called as a true host of a vOTU if it is supported by
     >= MIN_INDEPENDENT_HOSTS (default 2) distinct host_sample genomes bearing
     an independent spacer match.
      Rationale: on this dataset, 27% of raw vOTU-species calls rested on a single host
     sample / single spacer, and apparent host range (raw species_count)
     correlated strongly with spacer_count and host_count (i.e. with how
     thoroughly a phage happened to be sampled in the CRISPR spacer database,
     rho ~0.7-0.8), not purely with biology. Requiring >=2 independent hosts
     removes this sampling-depth confound. This yields species_count_robust,
     the primary host-range metric to use downstream (e.g. in gene-association
     tests); the original single-hit-inclusive species_count is retained
     alongside it for comparison/transparency only, not as the headline metric.
  3. checkv_quality per vOTU is assigned by true majority vote (Counter) across
     all member sequences with a quality call, not an arbitrary pick.
  4. Specialist / generalist buckets are defined on species_count_robust using
     one explicit, consistent threshold (GENERALIST_MIN_SPECIES), rather than
     a threshold that disagreed with its own print label.
===============================================================================
"""

import pandas as pd
from collections import defaultdict, Counter

MIN_INDEPENDENT_HOSTS = 2   # min. distinct host_sample genomes to count a species as a real host
GENERALIST_MIN_SPECIES = 2  # threshold used for the "generalist" bucket, applied consistently below

# Read the master table
df = pd.read_csv('master_spacer_table_complete.tsv', sep='\t')

# Filter out self-hits
df_filtered = df[df['self_hit'] == False]

print("="*70)
print("vOTU-LEVEL HOST RANGE ANALYSIS")
print("="*70)
print(f"Total spacer hits (non-self): {len(df_filtered)}")
print(f"Unique vOTUs: {df_filtered['votu'].nunique()}")
print(f"Unique species: {df_filtered['host_species'].nunique()}")

# --- Per (vOTU, species) support: how many independent spacers / host genomes back each call ---
support = (
    df_filtered.groupby(['votu', 'host_species'])
    .agg(n_spacers=('qseqid', 'nunique'), n_hosts=('host_sample', 'nunique'))
    .reset_index()
)

# --- Collect per-vOTU sets as before (raw) ---
votu_species = defaultdict(set)
votu_phages = defaultdict(set)
votu_spacers = defaultdict(set)
votu_hosts = defaultdict(set)
votu_quality = defaultdict(list)

for _, row in df_filtered.iterrows():
    votu = row['votu']
    if pd.isna(votu):
        continue
    votu_species[votu].add(row['host_species'])
    votu_phages[votu].add(row['sseqid'])
    votu_spacers[votu].add(row['qseqid'])
    votu_hosts[votu].add(row['host_sample'])
    if pd.notna(row['checkv_quality']):
        votu_quality[votu].append(row['checkv_quality'])

# Robust species set per vOTU: only species with >= MIN_INDEPENDENT_HOSTS distinct host samples
robust_species = defaultdict(set)
for _, r in support[support['n_hosts'] >= MIN_INDEPENDENT_HOSTS].iterrows():
    robust_species[r['votu']].add(r['host_species'])

# Build vOTU summary
votu_results = []
for votu in votu_species:
    quality_counts = Counter(votu_quality.get(votu, ['Unknown']))
    quality = quality_counts.most_common(1)[0][0]  # true mode, not arbitrary pick

    robust_set = robust_species.get(votu, set())
    votu_results.append({
        'vOTU': votu,
        'species_count': len(votu_species[votu]),
        'species_list': ';'.join(sorted(votu_species[votu])),
        'species_count_robust': len(robust_set),
        'species_list_robust': ';'.join(sorted(robust_set)),
        'phage_count': len(votu_phages[votu]),
        'spacer_count': len(votu_spacers[votu]),
        'host_count': len(votu_hosts[votu]),
        'checkv_quality': quality,
    })

votu_df = pd.DataFrame(votu_results)
votu_df = votu_df.sort_values('species_count', ascending=False)

print("\nTOP 20 VOTUS BY HOST RANGE (SPECIES LEVEL) -- raw vs. robust:")
print(votu_df[['vOTU', 'species_count', 'species_count_robust', 'phage_count',
                'spacer_count', 'host_count']].head(20).to_string(index=False))

print(f"\nMean raw species_count:    {votu_df['species_count'].mean():.2f}")
print(f"Mean robust species_count: {votu_df['species_count_robust'].mean():.2f}")
print(f"vOTUs where robust count < raw count: "
      f"{(votu_df['species_count_robust'] < votu_df['species_count']).sum()} "
      f"({(votu_df['species_count_robust'] < votu_df['species_count']).mean()*100:.1f}%)")

# 2. Identify specialist vs generalist vOTUs (using ROBUST count, and one consistent threshold)
specialist_votus = votu_df[votu_df['species_count_robust'] <= 1]
generalist_votus = votu_df[votu_df['species_count_robust'] >= GENERALIST_MIN_SPECIES]

print(f"\nSPECIALIST VOTUS (robust species_count <= 1): "
      f"{len(specialist_votus)} ({len(specialist_votus)/len(votu_df)*100:.1f}%)")
print(f"GENERALIST VOTUS (robust species_count >= {GENERALIST_MIN_SPECIES}): "
      f"{len(generalist_votus)} ({len(generalist_votus)/len(votu_df)*100:.1f}%)")

# Save
votu_df.to_csv('vOTU_host_range_summary.tsv', sep='\t', index=False)
support.to_csv('votu_species_support.tsv', sep='\t', index=False)
specialist_votus.to_csv('specialist_vOTUs.tsv', sep='\t', index=False)
generalist_votus.to_csv('generalist_vOTUs.tsv', sep='\t', index=False)

# 3. Spacers driving broad host range (unchanged logic, kept for continuity)
spacer_species = defaultdict(set)
spacer_votus = defaultdict(set)
spacer_phages = defaultdict(set)
for _, row in df_filtered.iterrows():
    if pd.isna(row['votu']):
        continue
    spacer_species[row['qseqid']].add(row['host_species'])
    spacer_votus[row['qseqid']].add(row['votu'])
    spacer_phages[row['qseqid']].add(row['sseqid'])

broad_spacers = [
    {'spacer': s, 'species_count': len(sp), 'species_list': ';'.join(sorted(sp)),
     'vOTU_count': len(spacer_votus[s]), 'phage_count': len(spacer_phages[s])}
    for s, sp in spacer_species.items() if len(sp) > 1
]
broad_spacer_df = pd.DataFrame(broad_spacers)
if len(broad_spacer_df) > 0:
    broad_spacer_df = broad_spacer_df.sort_values('species_count', ascending=False)
    print("\nTOP 20 SPACERS WITH BROADEST HOST RANGE:")
    print(broad_spacer_df.head(20).to_string(index=False))
    broad_spacer_df.to_csv('broad_host_range_spacers.tsv', sep='\t', index=False)
else:
    print("\nNo broad host range spacers found (no spacer appears in multiple species)")

# 4. Species targeted by most vOTUs (unchanged)
species_votu_count = defaultdict(set)
for _, row in df_filtered.iterrows():
    if pd.notna(row['votu']):
        species_votu_count[row['host_species']].add(row['votu'])

species_summary = pd.DataFrame({
    'species': list(species_votu_count.keys()),
    'vOTU_count': [len(v) for v in species_votu_count.values()]
}).sort_values('vOTU_count', ascending=False)

print("\nTOP 10 SPECIES WITH MOST VOTU MATCHES:")
print(species_summary.head(10).to_string(index=False))

# 5. Distribution stats -- report both metrics
print("\n" + "="*70)
print("SUMMARY STATISTICS:")
for col in ['species_count', 'species_count_robust']:
    print(f"  [{col}] mean={votu_df[col].mean():.2f} median={votu_df[col].median():.0f} "
          f"max={votu_df[col].max()} min={votu_df[col].min()}")
print(f"  Total unique species: {len(species_summary)}")

print("\nDISTRIBUTION OF ROBUST HOST RANGE:")
range_counts = votu_df['species_count_robust'].value_counts().sort_index()
for i in range(0, min(11, int(range_counts.index.max()) + 1)):
    if i in range_counts.index:
        print(f"  {i} species: {range_counts[i]} vOTUs")
print("="*70)

# 6. Summary report
with open('host_range_summary_report.txt', 'w') as f:
    f.write("VOTU HOST RANGE ANALYSIS REPORT (raw + confidence-filtered)\n")
    f.write("="*70 + "\n")
    f.write(f"Total vOTUs: {len(votu_df)}\n")
    f.write(f"MIN_INDEPENDENT_HOSTS threshold for 'robust': {MIN_INDEPENDENT_HOSTS}\n")
    f.write(f"Specialist vOTUs (robust species_count<=1): {len(specialist_votus)} "
            f"({len(specialist_votus)/len(votu_df)*100:.1f}%)\n")
    f.write(f"Generalist vOTUs (robust species_count>={GENERALIST_MIN_SPECIES}): {len(generalist_votus)} "
            f"({len(generalist_votus)/len(votu_df)*100:.1f}%)\n")
    f.write(f"Mean raw species_count: {votu_df['species_count'].mean():.2f}\n")
    f.write(f"Mean robust species_count: {votu_df['species_count_robust'].mean():.2f}\n")
    f.write(f"Max robust species_count: {votu_df['species_count_robust'].max()}\n\n")
    f.write("Top 5 broadest host range vOTUs (by robust count):\n")
    f.write(votu_df.sort_values('species_count_robust', ascending=False).head(5).to_string())
