## =============================================================================
## vOTU alpha/beta diversity analysis, stratified by Vibrio species and
## isolation source
##
## DATA-ONLY SCRIPT: writes tables only, no figures. Use
## votu_diversity_figures.R (separate script) to generate publication-quality
## plots from these outputs.
##
## Inputs (edit paths below):
##   1. Cuadovircetes_ani_clusters.tsv
##        col1 = representative prophage_id (used as vOTU_id)
##        col2 = comma-separated member prophage_ids (representative included)
##        prophage_id format: <BioSample>_p<N>  e.g. SAMN35557946_p2
##   2. METADATA_phage_hosts.csv
##        columns: Genome, BioSample, Assigned Taxonomy, Strain,
##                 Source_category, Isolation_source
##
## Outputs written to OUT_DIR:
##   genome_votu_matrix_long.tsv / _wide.tsv   genome x vOTU count matrix
##   genome_metadata_final.tsv                 cleaned metadata incl. grouped species
##   alpha_diversity_per_genome.tsv            richness/Shannon/Simpson per genome
##   alpha_diversity_kruskal_results.tsv
##   alpha_diversity_summary_by_species.tsv    median/IQR/mean/sd per group
##   alpha_diversity_summary_by_isolation_source.tsv
##   beta_diversity_permanova_results.tsv
##   beta_diversity_betadisper_results.tsv
##   pcoa_coordinates.tsv                      PCoA1/PCoA2 per genome (+ groups)
##   pcoa_variance_explained.tsv               % variance explained per axis
##   subsampling_permanova_iterations.tsv      per-iteration R2/p (supplementary)
##   subsampling_permanova_summary.tsv         summary stats across iterations
##   subsampling_dropped_groups.tsv            groups excluded from the species
##                                              subsampling check (n below floor)
##
## Note on threading: N_CORES (set below) parallelizes the PERMANOVA
## permutation calls and the subsampling iteration loop. This relies on
## fork-based parallelism (mclapply/vegan's parallel= arg) and works on
## Linux/Mac; on Windows it silently falls back to serial execution.
## =============================================================================

## ---- 0. Setup ---------------------------------------------------------------

library(tidyverse)
library(vegan)
library(parallel)

# Multi-threading: used for the permutation-heavy PERMANOVA calls (Section 8)
# and to parallelize across subsampling iterations (Section 11). Leaves 2
# cores free for the OS/other work rather than claiming all of them — bump
# this up manually if you want to use every core on the 28-core workstation.
N_CORES <- max(1, parallel::detectCores() - 2)
message(sprintf("Detected %d cores; using %d for parallelized steps.", parallel::detectCores(), N_CORES))

CLUSTER_FILE   <- "Cuadovircetes_ani_clusters.tsv"
METADATA_FILE  <- "vOTU_related_genomes_supplementary_Data.csv"
OUT_DIR        <- "diversity_outputs"
MIN_SPECIES_N  <- 5     # species with fewer genomes get pooled into "Other Vibrio spp." (main analysis grouping)

dir.create(OUT_DIR, showWarnings = FALSE)

## ---- 1. Parse the vOTU cluster membership file into long format ------------

cluster_raw <- read_tsv(
  CLUSTER_FILE,
  col_names = c("votu_id", "members"),
  col_types = cols(.default = "c")
)

# One row per prophage instance, tagged with its vOTU cluster (representative ID)
prophage_votu <- cluster_raw %>%
  separate_rows(members, sep = ",") %>%
  transmute(
    prophage_id = members,
    votu_id     = votu_id
  )

# Recover genome_id (BioSample) by stripping the trailing _p<N> suffix
prophage_votu <- prophage_votu %>%
  mutate(genome_id = str_remove(prophage_id, "_p[0-9]+$"))

stopifnot(!any(is.na(prophage_votu$genome_id)))

## ---- 2. Build genome x vOTU count matrix (long -> wide) ---------------------

genome_votu_long <- prophage_votu %>%
  count(genome_id, votu_id, name = "count")

write_tsv(genome_votu_long, file.path(OUT_DIR, "genome_votu_matrix_long.tsv"))

genome_votu_wide <- genome_votu_long %>%
  pivot_wider(names_from = votu_id, values_from = count, values_fill = 0)

## ---- 3. Load and clean metadata ---------------------------------------------

meta_raw <- read_csv(METADATA_FILE, col_types = cols(.default = "c"))

# Species grouping: species with n >= MIN_SPECIES_N kept standalone,
# everything below pooled into "Other Vibrio spp."
species_counts <- meta_raw %>% count(`Assigned Taxonomy`, name = "n_genomes")

meta <- meta_raw %>%
  left_join(species_counts, by = "Assigned Taxonomy") %>%
  mutate(
    species          = `Assigned Taxonomy`,
    species_grouped  = if_else(n_genomes >= MIN_SPECIES_N, species, "Other Vibrio spp."),
    isolation_source = Isolation_source,     # kept as-is; "Unknown" retained as its own category
    genome_id        = BioSample
  ) %>%
  select(genome_id, genome_accession = Genome, species, species_grouped, isolation_source)

write_tsv(meta, file.path(OUT_DIR, "genome_metadata_final.tsv"))

## ---- 4. Join matrix with metadata; sanity-check for zero-prophage rows -----

# METADATA_phage_hosts.csv only includes genomes that carry >=1 Caudoviricetes
# prophage (confirmed: 2,925 host genomes / 3,933 Caudoviricetes prophages),
# so every genome_id here SHOULD match at least one row in the cluster file.
# A genome ending up all-zero after the join is therefore NOT expected
# biology — it almost always means a genome_id (BioSample) mismatch between
# the two files (whitespace, version suffix, formatting difference), so this
# check is a data-integrity sanity check, not a real "zero-prophage" result.

all_genome_ids <- unique(meta$genome_id)

genome_votu_full <- genome_votu_wide %>%
  right_join(tibble(genome_id = all_genome_ids), by = "genome_id") %>%
  mutate(across(-genome_id, ~replace_na(.x, 0)))

zero_prophage_ids <- genome_votu_full %>%
  rowwise() %>%
  mutate(total = sum(c_across(-genome_id))) %>%
  ungroup() %>%
  filter(total == 0) %>%
  pull(genome_id)

n_zero_prophage <- length(zero_prophage_ids)

if (n_zero_prophage > 0) {
  warning(sprintf(
    "%d genome(s) matched zero vOTUs after the join — this is unexpected given the input data and likely indicates a genome_id/BioSample mismatch between METADATA_phage_hosts.csv and %s. Inspect these IDs before trusting downstream results:\n%s",
    n_zero_prophage, CLUSTER_FILE, paste(head(zero_prophage_ids, 10), collapse = ", ")
  ))
} else {
  message("OK: every genome in metadata matched >=1 vOTU in the cluster file, as expected.")
}

write_tsv(genome_votu_full, file.path(OUT_DIR, "genome_votu_matrix_wide.tsv"))

## ---- 5. Alpha diversity per genome -------------------------------------------

mat <- genome_votu_full %>% select(-genome_id) %>% as.matrix()
rownames(mat) <- genome_votu_full$genome_id

alpha_df <- tibble(
  genome_id = rownames(mat),
  richness  = specnumber(mat),
  shannon   = diversity(mat, index = "shannon"),
  simpson   = diversity(mat, index = "simpson")
) %>%
  left_join(meta, by = "genome_id")

write_tsv(alpha_df, file.path(OUT_DIR, "alpha_diversity_per_genome.tsv"))

## ---- 6. Alpha diversity group comparisons (Kruskal-Wallis) -----------------

alpha_long <- alpha_df %>%
  pivot_longer(cols = c(richness, shannon, simpson), names_to = "metric", values_to = "value")

kruskal_results <- alpha_long %>%
  group_by(metric) %>%
  summarise(
    species_kw_p = kruskal.test(value ~ species_grouped)$p.value,
    isolation_kw_p = kruskal.test(value ~ isolation_source)$p.value,
    .groups = "drop"
  )

write_tsv(kruskal_results, file.path(OUT_DIR, "alpha_diversity_kruskal_results.tsv"))

## ---- 6b. Alpha diversity summary tables (for manuscript text/tables) -------
##
## Per-group median/IQR/mean/sd/n for each metric, sorted by descending
## median so the highest- and lowest-diversity groups are easy to read off
## directly.

alpha_summary_species <- alpha_long %>%
  group_by(metric, species_grouped) %>%
  summarise(
    n = n(),
    mean = mean(value),
    sd = sd(value),
    median = median(value),
    Q1 = quantile(value, 0.25),
    Q3 = quantile(value, 0.75),
    min = min(value),
    max = max(value),
    .groups = "drop"
  ) %>%
  arrange(metric, desc(median))

alpha_summary_isolation <- alpha_long %>%
  group_by(metric, isolation_source) %>%
  summarise(
    n = n(),
    mean = mean(value),
    sd = sd(value),
    median = median(value),
    Q1 = quantile(value, 0.25),
    Q3 = quantile(value, 0.75),
    min = min(value),
    max = max(value),
    .groups = "drop"
  ) %>%
  arrange(metric, desc(median))

write_tsv(alpha_summary_species, file.path(OUT_DIR, "alpha_diversity_summary_by_species.tsv"))
write_tsv(alpha_summary_isolation, file.path(OUT_DIR, "alpha_diversity_summary_by_isolation_source.tsv"))

## ---- 7. Beta diversity: drop all-zero genomes, build distance matrices -----

nonzero_ids <- genome_votu_full %>%
  rowwise() %>%
  mutate(total = sum(c_across(-genome_id))) %>%
  ungroup() %>%
  filter(total > 0) %>%
  pull(genome_id)

mat_nonzero <- mat[rownames(mat) %in% nonzero_ids, , drop = FALSE]

meta_nonzero <- meta %>%
  filter(genome_id %in% nonzero_ids) %>%
  arrange(match(genome_id, rownames(mat_nonzero)))

stopifnot(identical(meta_nonzero$genome_id, rownames(mat_nonzero)))

# Presence/absence for Jaccard; raw counts for Bray-Curtis
mat_pa <- decostand(mat_nonzero, method = "pa")

jaccard_dist <- vegdist(mat_pa, method = "jaccard", binary = TRUE)
bray_dist    <- vegdist(mat_nonzero, method = "bray")

## ---- 8. PERMANOVA (adonis2) -------------------------------------------------

set.seed(42)

permanova_species <- adonis2(
  jaccard_dist ~ species_grouped,
  data = meta_nonzero, permutations = 999, parallel = N_CORES
)

permanova_isolation <- adonis2(
  jaccard_dist ~ isolation_source,
  data = meta_nonzero, permutations = 999, parallel = N_CORES
)

permanova_combined <- adonis2(
  jaccard_dist ~ species_grouped + isolation_source,
  data = meta_nonzero, permutations = 999, by = "margin", parallel = N_CORES
)

permanova_results <- bind_rows(
  as_tibble(permanova_species, rownames = "term") %>% mutate(model = "species_only", distance = "jaccard"),
  as_tibble(permanova_isolation, rownames = "term") %>% mutate(model = "isolation_only", distance = "jaccard"),
  as_tibble(permanova_combined, rownames = "term") %>% mutate(model = "combined_margin", distance = "jaccard")
)

# Repeat with Bray-Curtis for comparison (only meaningful if counts vary, i.e.
# some genomes carry multiple copies of the same vOTU — otherwise Bray-Curtis
# on 0/1 data collapses to the same information as Jaccard)
permanova_species_bc <- adonis2(bray_dist ~ species_grouped, data = meta_nonzero, permutations = 999, parallel = N_CORES)
permanova_isolation_bc <- adonis2(bray_dist ~ isolation_source, data = meta_nonzero, permutations = 999, parallel = N_CORES)

permanova_results <- bind_rows(
  permanova_results,
  as_tibble(permanova_species_bc, rownames = "term") %>% mutate(model = "species_only", distance = "bray"),
  as_tibble(permanova_isolation_bc, rownames = "term") %>% mutate(model = "isolation_only", distance = "bray")
)

write_tsv(permanova_results, file.path(OUT_DIR, "beta_diversity_permanova_results.tsv"))

## ---- 9. betadisper: check significant PERMANOVA isn't a dispersion artifact

bd_species <- betadisper(jaccard_dist, meta_nonzero$species_grouped)
bd_isolation <- betadisper(jaccard_dist, meta_nonzero$isolation_source)

bd_species_test <- permutest(bd_species, permutations = 999, parallel = N_CORES)
bd_isolation_test <- permutest(bd_isolation, permutations = 999, parallel = N_CORES)

betadisper_results <- tibble(
  grouping = c("species_grouped", "isolation_source"),
  F_statistic = c(bd_species_test$tab$F[1], bd_isolation_test$tab$F[1]),
  p_value = c(bd_species_test$tab$`Pr(>F)`[1], bd_isolation_test$tab$`Pr(>F)`[1])
)

write_tsv(betadisper_results, file.path(OUT_DIR, "beta_diversity_betadisper_results.tsv"))

## ---- 10. Ordination (PCoA) — coordinates only, no plotting -----------------

pcoa <- cmdscale(jaccard_dist, k = 2, eig = TRUE)
pcoa_df <- as_tibble(pcoa$points, .name_repair = "unique") %>%
  rename(PCoA1 = 1, PCoA2 = 2) %>%
  bind_cols(meta_nonzero)

var_explained <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 2)

write_tsv(pcoa_df, file.path(OUT_DIR, "pcoa_coordinates.tsv"))
write_tsv(
  tibble(axis = c("PCoA1", "PCoA2"), pct_variance_explained = var_explained),
  file.path(OUT_DIR, "pcoa_variance_explained.tsv")
)

## ---- 11. SUPPLEMENTARY: subsampling robustness check for PERMANOVA --------
##
## Group sizes are very unequal (species_grouped: 37-1186; isolation_source:
## 82-853). To check whether the PERMANOVA result in Section 8 is being
## driven by the large groups rather than reflecting a pattern consistent
## across group sizes, repeatedly subsample groups down to a common size and
## rerun PERMANOVA each time. A result that stays significant with a similar
## R^2 across iterations is robust to the sample-size imbalance.
##
## SPECIES_SUBSAMPLE_MIN_N sets a fixed floor (20) rather than using the
## dataset's true minimum (5, "Other Vibrio spp."): with 29 groups, sampling
## everyone down to n=5 fits 28 group-mean parameters to only 145 genomes,
## which mechanically inflates PERMANOVA R^2 regardless of any real effect
## (expected null-baseline R^2 ~ (29-1)/(145-1) = 19%, purely from the
## group-count-to-n ratio). At floor 20, groups below the floor are dropped
## from this check only (NOT from the main analysis above), leaving a much
## healthier ratio of model df to total df. Isolation source's natural
## minimum group (82, Environmental (other)) is already comfortably above
## any reasonable floor, so no groups are dropped there.

N_SUBSAMPLE_ITER          <- 100   # number of random subsampling iterations
SUBSAMPLE_PERMANOVA_PERMS <- 199   # permutations per adonis2 call
SUBSAMPLE_SEED            <- 123
SPECIES_SUBSAMPLE_MIN_N   <- 20    # fixed floor for species_grouped; groups below this are dropped from this check

jaccard_mat_full <- as.matrix(jaccard_dist)  # subset this per iteration rather than recomputing vegdist each time

run_subsample_permanova <- function(grouping_var, n_iter, n_perms, seed, n_cores, min_group_n = NULL) {
  set.seed(seed)
  group_sizes <- table(meta_nonzero[[grouping_var]])

  if (is.null(min_group_n)) {
    # No fixed floor requested: use the dataset's natural minimum group size
    target_n <- min(group_sizes)
    eligible_groups <- names(group_sizes)
    dropped_groups <- tibble(group = character(0), n = integer(0))
  } else {
    # Fixed floor: drop groups below it, sample exactly min_group_n from each retained group
    target_n <- min_group_n
    eligible_groups <- names(group_sizes[group_sizes >= min_group_n])
    dropped_groups <- tibble(
      group = names(group_sizes[group_sizes < min_group_n]),
      n = as.integer(group_sizes[group_sizes < min_group_n])
    )
  }

  message(sprintf(
    "Subsampling '%s': %d/%d groups retained (n>=%d), sampling n=%d from each, %d iterations across %d cores",
    grouping_var, length(eligible_groups), length(group_sizes), target_n, target_n, n_iter, n_cores
  ))
  if (nrow(dropped_groups) > 0) {
    message(sprintf(
      "  Dropped from this check (n < %d): %s",
      target_n, paste(sprintf("%s (n=%d)", dropped_groups$group, dropped_groups$n), collapse = ", ")
    ))
  }

  eligible_meta <- meta_nonzero %>% filter(.data[[grouping_var]] %in% eligible_groups)

  # Pre-generate a distinct random seed per iteration so results are
  # reproducible regardless of core count/scheduling order under mclapply
  # (mclapply does NOT guarantee RNG reproducibility across forked workers
  # unless each task sets its own seed explicitly).
  iter_seeds <- sample.int(.Machine$integer.max, n_iter)

  results_list <- mclapply(seq_len(n_iter), function(i) {
    set.seed(iter_seeds[i])

    sampled_ids <- eligible_meta %>%
      group_by(.data[[grouping_var]]) %>%
      slice_sample(n = target_n) %>%
      ungroup() %>%
      pull(genome_id)

    sub_dist <- as.dist(jaccard_mat_full[sampled_ids, sampled_ids])
    sub_meta <- eligible_meta %>% filter(genome_id %in% sampled_ids) %>%
      arrange(match(genome_id, sampled_ids))

    # parallel = 1 here deliberately: iterations are already parallelized
    # across cores via mclapply, so parallelizing permutations *within* each
    # call as well would oversubscribe cores and slow things down overall.
    fit <- adonis2(
      sub_dist ~ sub_meta[[grouping_var]],
      permutations = n_perms, parallel = 1
    )

    tibble(
      iteration = i,
      R2 = fit$R2[1],
      p_value = fit$`Pr(>F)`[1]
    )
  }, mc.cores = n_cores)

  list(
    results = bind_rows(results_list) %>%
      mutate(grouping = grouping_var, n_per_group = target_n, n_groups = length(eligible_groups)),
    dropped = dropped_groups %>% mutate(grouping = grouping_var)
  )
}

subsample_species   <- run_subsample_permanova("species_grouped", N_SUBSAMPLE_ITER, SUBSAMPLE_PERMANOVA_PERMS, SUBSAMPLE_SEED, N_CORES, min_group_n = SPECIES_SUBSAMPLE_MIN_N)
subsample_isolation <- run_subsample_permanova("isolation_source", N_SUBSAMPLE_ITER, SUBSAMPLE_PERMANOVA_PERMS, SUBSAMPLE_SEED, N_CORES, min_group_n = NULL)

subsample_results <- bind_rows(subsample_species$results, subsample_isolation$results)
write_tsv(subsample_results, file.path(OUT_DIR, "subsampling_permanova_iterations.tsv"))

dropped_groups_all <- bind_rows(subsample_species$dropped, subsample_isolation$dropped)
write_tsv(dropped_groups_all, file.path(OUT_DIR, "subsampling_dropped_groups.tsv"))

subsample_summary <- subsample_results %>%
  group_by(grouping, n_per_group, n_groups) %>%
  summarise(
    n_iterations = n(),
    mean_R2 = mean(R2),
    median_R2 = median(R2),
    R2_IQR_low = quantile(R2, 0.25),
    R2_IQR_high = quantile(R2, 0.75),
    pct_iterations_significant = mean(p_value < 0.05) * 100,
    .groups = "drop"
  )

# Compare against the full-data (unbalanced, all-groups) result from Section 8 for reference
full_data_species_R2   <- permanova_species$R2[1]
full_data_isolation_R2 <- permanova_isolation$R2[1]

subsample_summary <- subsample_summary %>%
  mutate(
    full_data_R2 = case_when(
      grouping == "species_grouped" ~ full_data_species_R2,
      grouping == "isolation_source" ~ full_data_isolation_R2
    )
  )

write_tsv(subsample_summary, file.path(OUT_DIR, "subsampling_permanova_summary.tsv"))

message("Subsampling robustness check complete. See subsampling_permanova_summary.tsv:")
print(subsample_summary)

## ---- 12. Summary printout ----------------------------------------------------

message("== DONE ==")
message(sprintf("Genomes in metadata: %d", nrow(meta)))
message(sprintf("Genomes with >=1 vOTU (used for beta diversity): %d", length(nonzero_ids)))
message(sprintf("vOTU clusters (Caudoviricetes): %d", ncol(mat)))
message("All data tables written to: ", OUT_DIR)
message("Run votu_diversity_figures.R against this directory to generate plots.")
