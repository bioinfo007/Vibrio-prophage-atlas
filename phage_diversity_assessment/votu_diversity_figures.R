## =============================================================================
## Publication-quality figures for the vOTU alpha/beta diversity analysis
##
## Reads the data tables produced by votu_diversity_analysis.R (does NOT
## recompute any statistics — pure plotting script) and writes figures with:
##   - larger, more legible points on the PCoA plots
##   - a qualitative color palette built to stay distinguishable even at
##     29 species-level groups (RColorBrewer Dark2 + Set1 + Set3 concatenated)
##   - 68% confidence ellipses on PCoA plots to make group separation
##     visually readable despite point overlap
##   - a faceted "highlight" version of the species PCoA (one panel per
##     species, that species colored, everything else in grey) as a second,
##     often clearer alternative to the single crowded-legend plot — this is
##     the standard fix for high-cardinality (many-group) ordination plots
##
## Run this AFTER votu_diversity_analysis.R has populated diversity_outputs/.
## =============================================================================

## ---- 0. Setup ---------------------------------------------------------------

library(tidyverse)
library(RColorBrewer)

IN_DIR  <- "diversity_outputs"
OUT_DIR <- "diversity_outputs/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Liberation Sans for consistency with prior figures — optional; falls back
# to the default sans font if showtext isn't installed or the font file
# isn't found.
FONT_FAMILY <- "sans"
if (requireNamespace("showtext", quietly = TRUE)) {
  library(showtext)
  font_path <- "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"
  if (file.exists(font_path)) {
    font_add("LiberationSans", regular = font_path)
    showtext_auto()
    FONT_FAMILY <- "LiberationSans"
  } else {
    message("Liberation Sans font file not found at ", font_path, " — using default sans font instead.")
  }
} else {
  message("Package 'showtext' not installed — using default sans font instead. Install with install.packages('showtext') to match prior figure styling.")
}

plot_theme <- theme_bw(base_size = 13, base_family = FONT_FAMILY) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.key.size = unit(0.4, "cm")
  )

## ---- 1. Qualitative color palette generator ---------------------------------
##
## A single evenly-spaced hue palette (e.g. scales::hue_pal) gets visually
## muddy past ~12-15 groups. Concatenating several ColorBrewer qualitative
## sets gives up to 8 (Dark2) + 9 (Set1) + 12 (Set3) = 29 distinguishable
## colors — exactly enough for the 29-group species_grouped variable.

make_qual_palette <- function(n) {
  base_colors <- c(
    brewer.pal(8, "Dark2"),
    brewer.pal(9, "Set1"),
    brewer.pal(12, "Set3")
  )
  base_colors <- unique(base_colors)
  if (n <= length(base_colors)) {
    base_colors[seq_len(n)]
  } else {
    colorRampPalette(base_colors)(n)
  }
}

## ---- 2. Load data -------------------------------------------------------------

alpha_df   <- read_tsv(file.path(IN_DIR, "alpha_diversity_per_genome.tsv"), show_col_types = FALSE)
alpha_long <- alpha_df %>%
  pivot_longer(cols = c(richness, shannon, simpson), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("richness", "shannon", "simpson")))

pcoa_df   <- read_tsv(file.path(IN_DIR, "pcoa_coordinates.tsv"), show_col_types = FALSE)
var_exp   <- read_tsv(file.path(IN_DIR, "pcoa_variance_explained.tsv"), show_col_types = FALSE)
pc1_label <- sprintf("PCoA1 (%.1f%%)", var_exp$pct_variance_explained[var_exp$axis == "PCoA1"])
pc2_label <- sprintf("PCoA2 (%.1f%%)", var_exp$pct_variance_explained[var_exp$axis == "PCoA2"])

subsample_iters <- read_tsv(file.path(IN_DIR, "subsampling_permanova_iterations.tsv"), show_col_types = FALSE)
subsample_summary <- read_tsv(file.path(IN_DIR, "subsampling_permanova_summary.tsv"), show_col_types = FALSE)

## ---- 3. Alpha diversity boxplots ---------------------------------------------
##
## Groups ordered by descending median (per metric's own ordering isn't
## used for x-axis order — richness ordering is used throughout so panels
## stay visually comparable) so the highest/lowest-diversity groups are
## easy to spot without hunting through a crowded, unordered axis.

species_order <- alpha_long %>%
  filter(metric == "richness") %>%
  group_by(species_grouped) %>%
  summarise(med = median(value), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(species_grouped)

isolation_order <- alpha_long %>%
  filter(metric == "richness") %>%
  group_by(isolation_source) %>%
  summarise(med = median(value), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(isolation_source)

alpha_long_species <- alpha_long %>% mutate(species_grouped = factor(species_grouped, levels = species_order))
alpha_long_isolation <- alpha_long %>% mutate(isolation_source = factor(isolation_source, levels = isolation_order))

species_palette <- setNames(make_qual_palette(length(species_order)), species_order)
isolation_palette <- setNames(brewer.pal(length(isolation_order), "Dark2"), isolation_order)

p_alpha_species <- alpha_long_species %>%
  ggplot(aes(x = species_grouped, y = value, fill = species_grouped)) +
  geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.5, linewidth = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = species_palette) +
  facet_wrap(~metric, scales = "free_y") +
  plot_theme +
  labs(x = "Vibrio species (grouped, ordered by median richness)", y = "Alpha diversity value",
       title = "vOTU alpha diversity by Vibrio species")

p_alpha_isolation <- alpha_long_isolation %>%
  ggplot(aes(x = isolation_source, y = value, fill = isolation_source)) +
  geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.5, linewidth = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = isolation_palette) +
  facet_wrap(~metric, scales = "free_y") +
  plot_theme +
  labs(x = "Isolation source (ordered by median richness)", y = "Alpha diversity value",
       title = "vOTU alpha diversity by isolation source")

ggsave(file.path(OUT_DIR, "alpha_diversity_boxplots_species.svg"), p_alpha_species, width = 12, height = 6.5)
ggsave(file.path(OUT_DIR, "alpha_diversity_boxplots_species.pdf"), p_alpha_species, width = 12, height = 6.5)
ggsave(file.path(OUT_DIR, "alpha_diversity_boxplots_isolation_source.svg"), p_alpha_isolation, width = 11, height = 6)
ggsave(file.path(OUT_DIR, "alpha_diversity_boxplots_isolation_source.pdf"), p_alpha_isolation, width = 11, height = 6)

## ---- 4. PCoA plots: single legend version, larger points + ellipses --------

p_pcoa_isolation <- pcoa_df %>%
  ggplot(aes(x = PCoA1, y = PCoA2, color = isolation_source, fill = isolation_source)) +
  geom_point(size = 2.6, alpha = 0.75, shape = 21, color = "white", stroke = 0.15) +
  stat_ellipse(aes(color = isolation_source), type = "norm", level = 0.68, linewidth = 0.8, fill = NA, geom = "path") +
  scale_fill_manual(values = isolation_palette) +
  scale_color_manual(values = isolation_palette) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(
    x = pc1_label, y = pc2_label,
    title = "vOTU community composition (Jaccard) by isolation source",
    subtitle = "Ellipses = 68% normal-distribution confidence region per group",
    fill = "Isolation source", color = "Isolation source"
  )

ggsave(file.path(OUT_DIR, "pcoa_isolation_source.svg"), p_pcoa_isolation, width = 9, height = 6.5)
ggsave(file.path(OUT_DIR, "pcoa_isolation_source.pdf"), p_pcoa_isolation, width = 9, height = 6.5)

p_pcoa_species <- pcoa_df %>%
  mutate(species_grouped = factor(species_grouped, levels = species_order)) %>%
  ggplot(aes(x = PCoA1, y = PCoA2, color = species_grouped, fill = species_grouped)) +
  geom_point(size = 2.2, alpha = 0.7, shape = 21, color = "white", stroke = 0.1) +
  scale_fill_manual(values = species_palette) +
  scale_color_manual(values = species_palette) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  ) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 3)), color = "none") +
  labs(
    x = pc1_label, y = pc2_label,
    title = "vOTU community composition (Jaccard) by Vibrio species",
    fill = "Species"
  )

ggsave(file.path(OUT_DIR, "pcoa_species_legend.svg"), p_pcoa_species, width = 11, height = 7)
ggsave(file.path(OUT_DIR, "pcoa_species_legend.pdf"), p_pcoa_species, width = 11, height = 7)

## ---- 5. PCoA species: faceted "highlight" small-multiples version ----------
##
## 29 groups in one legend inevitably gets visually crowded no matter how
## good the palette is. The standard fix for high-cardinality ordination
## plots: one small panel per group, that group colored, every other point
## shown in light grey as context. Restricted to species with n>=20 (the
## same well-powered set used in the subsampling robustness check) to keep
## the panel count reasonable and avoid a wall of near-empty panels for
## rare species — the full legend version above still shows everything.

highlight_species <- alpha_df %>% count(species_grouped) %>% filter(n >= 20) %>% pull(species_grouped)

pcoa_highlight <- map_dfr(highlight_species, function(sp) {
  pcoa_df %>%
    mutate(
      facet_species = sp,
      is_highlighted = species_grouped == sp
    )
})

p_pcoa_species_facet <- pcoa_highlight %>%
  arrange(is_highlighted) %>%  # plot grey points first, colored points on top
  ggplot(aes(x = PCoA1, y = PCoA2)) +
  geom_point(
    data = . %>% filter(!is_highlighted),
    color = "grey85", size = 1, alpha = 0.5
  ) +
  geom_point(
    data = . %>% filter(is_highlighted),
    aes(color = facet_species), size = 1.6, alpha = 0.85
  ) +
  scale_color_manual(values = species_palette, guide = "none") +
  facet_wrap(~facet_species, ncol = 4) +
  plot_theme +
  theme(
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 8, face = "bold")
  ) +
  labs(
    x = pc1_label, y = pc2_label,
    title = "vOTU community composition by Vibrio species (highlighted, n>=20 species only)",
    subtitle = "Grey = all other genomes; colored = genomes of the labeled species"
  )

ggsave(file.path(OUT_DIR, "pcoa_species_highlight_facets.svg"), p_pcoa_species_facet, width = 13, height = 10)
ggsave(file.path(OUT_DIR, "pcoa_species_highlight_facets.pdf"), p_pcoa_species_facet, width = 13, height = 10)

## ---- 6. Subsampling robustness: R2 distribution plot -----------------------

p_subsample <- subsample_iters %>%
  left_join(subsample_summary %>% select(grouping, full_data_R2), by = "grouping") %>%
  mutate(grouping = recode(grouping,
                            species_grouped = "Species (n>=20 groups only)",
                            isolation_source = "Isolation source")) %>%
  ggplot(aes(x = R2)) +
  geom_histogram(bins = 25, fill = "grey60", color = "white", linewidth = 0.2) +
  geom_vline(aes(xintercept = full_data_R2), color = "firebrick", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~grouping, scales = "free", ncol = 1) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(
    x = expression(R^2~"(subsampled, balanced PERMANOVA)"),
    y = "Count",
    title = "PERMANOVA R\u00b2 across balanced subsampling iterations",
    subtitle = "Dashed red line = R\u00b2 from full (unbalanced) dataset"
  )

ggsave(file.path(OUT_DIR, "subsampling_permanova_R2_distribution.svg"), p_subsample, width = 8, height = 8)
ggsave(file.path(OUT_DIR, "subsampling_permanova_R2_distribution.pdf"), p_subsample, width = 8, height = 8)

## ---- Done ---------------------------------------------------------------------

message("== Figures written to ", OUT_DIR, " ==")
message("  alpha_diversity_boxplots_species.svg/.pdf")
message("  alpha_diversity_boxplots_isolation_source.svg/.pdf")
message("  pcoa_isolation_source.svg/.pdf")
message("  pcoa_species_legend.svg/.pdf            (all 29 groups, single legend)")
message("  pcoa_species_highlight_facets.svg/.pdf  (n>=20 groups, one highlighted panel each — recommended for readability)")
message("  subsampling_permanova_R2_distribution.svg/.pdf")
