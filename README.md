# Vibrio-prophage-atlas

Code and analysis pipelines for a reference-independent, genus-scale survey of dsDNA prophages across *Vibrio* genomes, characterizing prophage taxonomic diversity, host-range evolution via diversity-generating retroelements (DGRs), lysis-module architecture, and functional cargo.

## Overview

This repository accompanies the study *"Genomic catalogue of Vibrio prophages reveals expanded viral diversity and convergent retroelement-associated host-range breadth"*, in which 3,916 dsDNA Caudoviricetes prophage sequences were recovered from 6,541 public *Vibrio* genome assemblies and clustered into 1,087 viral operational taxonomic units (vOTUs). The scripts here cover the full analysis pipeline, from prophage prediction through host-range reconstruction, DGR/lysis-module characterization, and functional cargo annotation.

Predicted prophage sequences and associated metadata are archived separately on Zenodo (DOI: *10.5281/zenodo.22844791*); this repository contains only the analysis code.

## Repository structure

```
Vibrio-prophage-atlas/
├── Consensus_prediction/          # Prophage prediction and consensus boundary calling
├── phage_diversity_assessment/    # vOTU clustering, taxonomic assignment, alpha/beta diversity analysis
├── Host_range_prediction/         # CRISPR spacer-based host range reconstruction
├── Functional_gene_association/   # Metabolic, resistance, and virulence cargo annotation
├── Lysis_cassette_association/    # Lysis-module architecture and lineage classification
├── LICENSE
└── README.md
```

### Consensus_prediction
Script for consensus prophage derivation from prophage prediction (using tools: VIBRANT and geNomad), and provides consensus phage sequences for CheckV-based trimming of host contamination and quality assessment.

### phage_diversity_assessment
Scripts for alpha/beta diversity analysis of prophage repertoires across host species and isolation sources (Kruskal-Wallis tests, PERMANOVA, betadisper, PCoA, and balanced subsampling robustness checks).

### Host_range_prediction
Scripts for reconstructing historical host associations from CRISPR spacer-protospacer matches:
- `01_filter_true_matches.py` — filters raw spacer-vs-prophage BLAST hits to high-confidence true matches (zero mismatches, zero gaps, near-full-length alignment), maps hits to vOTU clusters, flags self-targeting, and summarizes vOTU-level host range from cross-genome hits.
- `02_create_master_table.py` — merges filtered spacer hits with host species, CheckV quality, BACPHLIP lifestyle predictions, and isolation source metadata into a single master table.
- '03_votu_host_range_prediction.py' — main script which takes input 'master table' and assesses spacer-based host range both at genome-level, and vOTU-level.

### Functional_gene_association
Scripts for identifying and classifying prophage-encoded metabolic, antimicrobial resistance, and virulence-associated genes ("cargo genes").

### Lysis_cassette_association
Scripts for annotating lysis module architecture (e.g., holin/endolysin/spanin gene content and organization) and classifying lysis cassette lineages across vOTUs.

## Requirements

- Python ≥ 3.9
  - pandas
  - (add other packages as needed, e.g. numpy, biopython)
- R ≥ 4.2
  - tidyverse
  - vegan
  - RColorBrewer
  - (add other packages as needed)


## Usage

Each subdirectory contains scripts intended to be run in numbered order where applicable (e.g., `01_...`, `02_...`). Input file paths are currently set as relative paths at the top of each script and should be updated to match your local directory structure. See comments within each script for a full description of required input files, expected columns, and output files produced.

Example:
```bash
cd Host_range_prediction
python 01_filter_true_matches.py
python 02_create_master_table.py
```

## Data availability

- **Prophage sequences and metadata**: Zenodo, DOI *10.5281/zenodo.22844791* and *publication supplementary data*
- **Host genome assemblies**: all publicly available via NCBI GenBank; accession numbers are provided in prophage metadata tables.


## License

This project is licensed under the GNU General Public License v3.0 — see [LICENSE](LICENSE) for details.

## Contact

Aqib Javaid — *ajg4u.aj@gmail.com*

For questions about specific analyses, please open an issue in this repository.
