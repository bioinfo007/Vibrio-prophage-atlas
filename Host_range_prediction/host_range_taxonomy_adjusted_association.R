#!/usr/bin/env Rscript

# ================================================================
# ROBUST TAXONOMY-ADJUSTED HOST-RANGE ASSOCIATION
#
# Tests:
#   Host-range class ~ Isolation source, adjusted for taxonomy
#
# Primary robust test:
#   Permutation test of the association between isolation source and
#   host-range class within Deepest_taxonomy strata.
#
# This avoids relying on sparse Poisson log-linear models, which can
# suffer from fitted-zero/separation problems with many rare
# taxonomy x source combinations.
#
# Host range:
#   0 -> discarded
#   1 -> Specialist
#   2-3 -> Moderate
#   >3 -> Generalist
#
# Outputs:
#   SVG figures
#   TSV supplementary data
#   observed contingency tables
#   residuals
#   permutation statistics
# ================================================================

options(stringsAsFactors=FALSE, scipen=999)

HOST_FILE <- "vOTU_host_range_summary.tsv"
TAX_FILE <- "phage_taxonomy.csv"
SOURCE_FILE <- "Isolation_source.txt"

OUTDIR <- "host_range_taxonomy_adjusted_results"
N_PERM <- 99999
SEED <- 20260829

dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

# ---------------- packages ----------------
pkgs <- c("ggplot2","dplyr","tidyr","readr","stringr")
missing <- pkgs[!vapply(pkgs, requireNamespace, quietly=TRUE,
                        FUN.VALUE=logical(1))]
if(length(missing)>0){
  stop(
    "Missing packages: ",paste(missing,collapse=", "),
    "\nInstall with:\ninstall.packages(c(",
    paste(sprintf('"%s"',missing),collapse=", "), "))"
  )
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# ================================================================
# HELPERS
# ================================================================

mc_p <- function(extreme, nperm){
  # +1 correction prevents zero P-values
  (extreme + 1) / (nperm + 1)
}

chi_stat <- function(x, y){
  tab <- table(x, y)
  # Remove empty rows/columns within each taxonomy stratum.
  # This prevents NaN statistics for strata containing only one
  # source or one host-range class.
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop=FALSE]
  if(nrow(tab)<2 || ncol(tab)<2) return(0)
  ans <- suppressWarnings(chisq.test(tab, correct=FALSE))
  val <- as.numeric(ans$statistic)
  if(!is.finite(val)) return(0)
  val
}

cramers_v <- function(tab){
  z <- suppressWarnings(chisq.test(tab, correct=FALSE))
  n <- sum(tab)
  k <- min(nrow(tab)-1,ncol(tab)-1)
  if(n<=0 || k<=0) return(NA_real_)
  sqrt(as.numeric(z$statistic)/(n*k))
}

fmtp <- function(x){
  format.pval(as.numeric(x),digits=6,eps=1e-5)
}

# ================================================================
# READ HOST RANGE
# ================================================================

if(!file.exists(HOST_FILE))
  stop("Cannot find ",HOST_FILE)

host0 <- read_tsv(HOST_FILE,show_col_types=FALSE,progress=FALSE)

if(!all(c("vOTU","species_count_robust") %in% names(host0)))
  stop("Host file must contain vOTU and species_count_robust.")

host <- host0 %>%
  transmute(
    vOTU=as.character(vOTU),
    species_count_robust=
      suppressWarnings(as.numeric(species_count_robust))
  ) %>%
  filter(!is.na(vOTU),vOTU!="",!is.na(species_count_robust))

n_total <- nrow(host)
n_zero <- sum(host$species_count_robust==0,na.rm=TRUE)

# ================================================================
# READ TAXONOMY
# ================================================================

if(!file.exists(TAX_FILE))
  stop("Cannot find ",TAX_FILE)

tax0 <- read_csv(TAX_FILE,show_col_types=FALSE,progress=FALSE)

required_tax <- c("vOTU","Deepest_taxonomy","Family_subfamily")
if(!all(required_tax %in% names(tax0)))
  stop("Taxonomy file must contain: ",paste(required_tax,collapse=", "))

tax <- tax0 %>%
  transmute(
    vOTU=as.character(vOTU),
    Deepest_taxonomy=str_squish(as.character(Deepest_taxonomy)),
    Family_subfamily=str_squish(as.character(Family_subfamily))
  ) %>%
  distinct(vOTU,.keep_all=TRUE)

# ================================================================
# READ ISOLATION SOURCE
# ================================================================

if(!file.exists(SOURCE_FILE))
  stop("Cannot find ",SOURCE_FILE)

z <- readLines(SOURCE_FILE,warn=FALSE,encoding="UTF-8")
z <- z[nzchar(str_trim(z))]

if(length(z)<2)
  stop("Isolation_source.txt contains no records.")

source_df <- tibble(raw=z[-1]) %>%
  mutate(
    genome_id=str_trim(sub(",.*$","",raw)),
    isolation_source=str_trim(sub("^[^,]*,","",raw)),
    isolation_source=str_remove(isolation_source,"^\\("),
    isolation_source=str_remove(isolation_source,"\\)$"),
    isolation_source=str_squish(isolation_source)
  ) %>%
  filter(genome_id!="",isolation_source!="") %>%
  distinct(genome_id,.keep_all=TRUE) %>%
  select(genome_id,isolation_source)

# ================================================================
# MERGE
# ================================================================

dat <- host %>%
  mutate(
    genome_id=str_remove(vOTU,"_p[^_]+$")
  ) %>%
  left_join(tax,by="vOTU") %>%
  left_join(source_df,by="genome_id") %>%
  mutate(
    host_range_class=case_when(
      species_count_robust==1 ~ "Specialist",
      species_count_robust %in% c(2,3) ~ "Moderate",
      species_count_robust>3 ~ "Generalist",
      TRUE ~ NA_character_
    )
  )

analysis <- dat %>%
  filter(
    species_count_robust>0,
    !is.na(host_range_class),
    !is.na(isolation_source),isolation_source!="",
    !is.na(Deepest_taxonomy),Deepest_taxonomy!=""
  ) %>%
  mutate(
    host_range_class=factor(
      host_range_class,
      levels=c("Specialist","Moderate","Generalist")
    ),
    isolation_source=factor(isolation_source),
    taxonomy=factor(Deepest_taxonomy)
  )

if(nrow(analysis)<10)
  stop("Too few complete observations: ",nrow(analysis))

# ================================================================
# SAVE COMPLETE SUPPLEMENTARY DATA
# ================================================================

supp <- analysis %>%
  transmute(
    vOTU,genome_id,species_count_robust,
    host_range_class,
    isolation_source,
    Deepest_taxonomy,
    Family_subfamily
  ) %>%
  arrange(isolation_source,Deepest_taxonomy,host_range_class,vOTU)

write_tsv(
  supp,
  file.path(OUTDIR,"taxonomy_adjusted_supplementary.tsv")
)

# ================================================================
# DIAGNOSTICS
# ================================================================

diagnostics <- data.frame(
  metric=c(
    "Total host-range records",
    "Zero robust species count discarded",
    "Non-zero records",
    "Complete cases",
    "Isolation-source groups",
    "Deepest-taxonomy groups"
  ),
  value=c(
    n_total,
    n_zero,
    sum(host$species_count_robust>0,na.rm=TRUE),
    nrow(analysis),
    nlevels(analysis$isolation_source),
    nlevels(analysis$taxonomy)
  ),
  stringsAsFactors=FALSE
)

write_tsv(
  diagnostics,
  file.path(OUTDIR,"taxonomy_adjusted_diagnostics.tsv")
)

# ================================================================
# OBSERVED SOURCE x HOST-RANGE ASSOCIATION
# ================================================================

tab_source <- table(
  analysis$isolation_source,
  analysis$host_range_class
)

tab_source <- tab_source[
  rowSums(tab_source)>0,
  colSums(tab_source)>0,
  drop=FALSE
]

chi_source <- suppressWarnings(
  chisq.test(tab_source,correct=FALSE)
)

V_source <- cramers_v(tab_source)

# ================================================================
# OBSERVED TAXONOMY x HOST-RANGE ASSOCIATION
# ================================================================

tab_tax <- table(
  analysis$taxonomy,
  analysis$host_range_class
)

tab_tax <- tab_tax[
  rowSums(tab_tax)>0,
  colSums(tab_tax)>0,
  drop=FALSE
]

chi_tax <- suppressWarnings(
  chisq.test(tab_tax,correct=FALSE)
)

V_tax <- cramers_v(tab_tax)

# ================================================================
# TAXONOMY-ADJUSTED TEST
#
# Null hypothesis:
#   After conditioning on taxonomy, isolation source is unrelated
#   to host-range class.
#
# Permutation scheme:
#   Shuffle isolation-source labels WITHIN each taxonomy group.
#
# This preserves:
#   - taxonomy composition
#   - number of observations per taxonomy
#   - number of observations assigned to each isolation source
#     within each taxonomy
#
# It destroys:
#   - source x host-range association within taxonomy
#
# Statistic:
#   Sum of Pearson chi-square statistics across taxonomy strata.
# ================================================================

observed_stratified_stat <- 0

strata <- split(
  seq_len(nrow(analysis)),
  analysis$taxonomy
)

strata <- strata[
  vapply(
    strata,
    function(ii){
      length(unique(analysis$isolation_source[ii]))>=2 &&
        length(unique(analysis$host_range_class[ii]))>=2
    },
    logical(1)
  )
]

if(length(strata)==0){
  stop(
    "No taxonomy strata contain variation in both isolation source ",
    "and host-range class; taxonomy-adjusted permutation test cannot ",
    "be performed."
  )
}

for(ii in strata){
  observed_stratified_stat <-
    observed_stratified_stat +
    chi_stat(
      analysis$isolation_source[ii],
      analysis$host_range_class[ii]
    )
}

# ------------------------------------------------
# Permutation
# ------------------------------------------------

set.seed(SEED)

perm_stats <- numeric(N_PERM)

source_values <- as.character(analysis$isolation_source)

for(b in seq_len(N_PERM)){

  perm_source <- source_values

  for(ii in strata){
    perm_source[ii] <- sample(
      source_values[ii],
      length(ii),
      replace=FALSE
    )
  }

  s <- 0

  for(ii in strata){
    s <- s +
      chi_stat(
        perm_source[ii],
        analysis$host_range_class[ii]
      )
  }

  perm_stats[b] <- s
}

extreme <- sum(
  perm_stats >= observed_stratified_stat
)

if(!is.finite(observed_stratified_stat) || !is.finite(extreme)){
  stop("Adjusted permutation statistic became non-finite. Check taxonomy/source strata.")
}
p_adjusted <- mc_p(extreme,N_PERM)

# ================================================================
# EFFECT SIZE FOR ADJUSTED TEST
# ================================================================

# Convert summed stratified chi-square to a standardized effect size.
# This is a descriptive partial association measure, not identical
# to ordinary Cramer's V.

N <- nrow(analysis)

K_source <- nlevels(droplevels(analysis$isolation_source))
K_hr <- nlevels(droplevels(analysis$host_range_class))

adjusted_effect <- sqrt(
  observed_stratified_stat /
  (N * min(K_source-1,K_hr-1))
)

# ================================================================
# SOURCE RESIDUALS
# ================================================================

source_obs <- as.data.frame(as.table(tab_source))
names(source_obs) <- c(
  "isolation_source",
  "host_range_class",
  "observed_count"
)

source_exp <- as.data.frame(as.table(chi_source$expected))
names(source_exp) <- c(
  "isolation_source",
  "host_range_class",
  "expected_count"
)

source_res <- as.data.frame(as.table(chi_source$stdres))
names(source_res) <- c(
  "isolation_source",
  "host_range_class",
  "standardized_residual"
)

source_res <- source_res %>%
  left_join(source_obs,
            by=c("isolation_source","host_range_class")) %>%
  left_join(source_exp,
            by=c("isolation_source","host_range_class")) %>%
  mutate(
    approximate_p=
      2*pnorm(-abs(standardized_residual)),
    BH_adjusted_p=
      p.adjust(approximate_p,method="BH"),
    interpretation=case_when(
      standardized_residual>=2 ~ "Enriched",
      standardized_residual<=-2 ~ "Depleted",
      TRUE ~ "No strong deviation"
    )
  ) %>%
  arrange(desc(abs(standardized_residual)))

write_tsv(
  source_res,
  file.path(
    OUTDIR,
    "isolation_source_host_range_residuals.tsv"
  )
)

# ================================================================
# TAXONOMY RESIDUALS
# ================================================================

tax_obs <- as.data.frame(as.table(tab_tax))
names(tax_obs) <- c(
  "taxonomy",
  "host_range_class",
  "observed_count"
)

tax_exp <- as.data.frame(as.table(chi_tax$expected))
names(tax_exp) <- c(
  "taxonomy",
  "host_range_class",
  "expected_count"
)

tax_res <- as.data.frame(as.table(chi_tax$stdres))
names(tax_res) <- c(
  "taxonomy",
  "host_range_class",
  "standardized_residual"
)

tax_res <- tax_res %>%
  left_join(tax_obs,
            by=c("taxonomy","host_range_class")) %>%
  left_join(tax_exp,
            by=c("taxonomy","host_range_class")) %>%
  mutate(
    approximate_p=
      2*pnorm(-abs(standardized_residual)),
    BH_adjusted_p=
      p.adjust(approximate_p,method="BH"),
    interpretation=case_when(
      standardized_residual>=2 ~ "Enriched",
      standardized_residual<=-2 ~ "Depleted",
      TRUE ~ "No strong deviation"
    )
  ) %>%
  arrange(desc(abs(standardized_residual)))

write_tsv(
  tax_res,
  file.path(
    OUTDIR,
    "taxonomy_host_range_residuals.tsv"
  )
)

# ================================================================
# PERMUTATION SUMMARY
# ================================================================

perm_summary <- data.frame(
  statistic="Taxonomy-stratified source x host-range association",
  observed_statistic=observed_stratified_stat,
  permutations=N_PERM,
  extreme_permutations=extreme,
  permutation_p=p_adjusted,
  adjusted_effect_size=adjusted_effect,
  taxonomy_strata_tested=length(strata),
  stringsAsFactors=FALSE
)

write_tsv(
  perm_summary,
  file.path(
    OUTDIR,
    "taxonomy_adjusted_permutation_test.tsv"
  )
)

# ================================================================
# MAIN STATISTICS
# ================================================================

stats <- data.frame(
  analysis=c(
    "Unadjusted isolation source",
    "Unadjusted phage taxonomy",
    "Taxonomy-adjusted isolation source"
  ),
  test=c(
    "Pearson chi-square",
    "Pearson chi-square",
    "Within-taxonomy permutation"
  ),
  statistic=as.numeric(c(
    unname(chi_source$statistic),
    unname(chi_tax$statistic),
    observed_stratified_stat
  )),
  df=as.numeric(c(
    unname(chi_source$parameter),
    unname(chi_tax$parameter),
    NA_real_
  )),
  p_value=as.numeric(c(
    chi_source$p.value,
    chi_tax$p.value,
    p_adjusted
  )),
  effect_size=as.numeric(c(
    V_source,
    V_tax,
    adjusted_effect
  )),
  N=as.numeric(c(
    nrow(analysis),
    nrow(analysis),
    nrow(analysis)
  )),
  stringsAsFactors=FALSE
)

write_tsv(
  stats,
  file.path(
    OUTDIR,
    "taxonomy_adjusted_statistics.tsv"
  )
)

# ================================================================
# FIGURE 1: SOURCE x HOST RANGE
# ================================================================

plot_source <- source_obs %>%
  group_by(isolation_source) %>%
  mutate(
    total=sum(observed_count),
    percent=100*observed_count/total
  ) %>%
  ungroup() %>%
  mutate(
    host_range_class=factor(
      host_range_class,
      levels=c("Specialist","Moderate","Generalist")
    )
  )

source_order <- plot_source %>%
  group_by(isolation_source) %>%
  summarise(n=sum(observed_count),.groups="drop") %>%
  arrange(desc(n)) %>%
  mutate(label=paste0(isolation_source,"\n(n=",n,")"))

plot_source <- plot_source %>%
  left_join(source_order,by="isolation_source") %>%
  mutate(x=factor(label,levels=source_order$label))

p1 <- ggplot(
  plot_source,
  aes(x=x,y=percent,fill=host_range_class)
) +
  geom_col(width=.78,colour="white",linewidth=.25) +
  geom_text(
    aes(
      label=ifelse(
        observed_count>0,
        paste0(
          observed_count,
          " (",
          sprintf("%.1f",percent),
          "%)"
        ),
        ""
      )
    ),
    position=position_stack(vjust=.5),
    size=3
  ) +
  scale_y_continuous(
    limits=c(0,100),
    breaks=seq(0,100,20),
    expand=c(0,0),
    labels=function(x) paste0(x,"%")
  ) +
  labs(
    x="Isolation source",
    y="Host-range composition (%)",
    fill="Host-range class"
  ) +
  theme_classic(base_size=12) +
  theme(
    axis.text.x=element_text(angle=40,hjust=1,vjust=1),
    axis.title=element_text(face="bold"),
    legend.title=element_text(face="bold")
  )

ggsave(
  file.path(
    OUTDIR,
    "complete_case_isolation_source_host_range.svg"
  ),
  p1,
  width=max(
    10,
    min(20,4+.45*nlevels(analysis$isolation_source))
  ),
  height=7,
  units="in",
  device="svg",
  bg="white"
)

# ================================================================
# FIGURE 2: TAXONOMY x HOST RANGE
# ================================================================

plot_tax <- tax_obs %>%
  group_by(taxonomy) %>%
  mutate(
    total=sum(observed_count),
    percent=100*observed_count/total
  ) %>%
  ungroup() %>%
  mutate(
    host_range_class=factor(
      host_range_class,
      levels=c("Specialist","Moderate","Generalist")
    )
  )

tax_order <- plot_tax %>%
  group_by(taxonomy) %>%
  summarise(n=sum(observed_count),.groups="drop") %>%
  arrange(desc(n)) %>%
  mutate(label=paste0(taxonomy,"\n(n=",n,")"))

plot_tax <- plot_tax %>%
  left_join(tax_order,by="taxonomy") %>%
  mutate(x=factor(label,levels=tax_order$label))

p2 <- ggplot(
  plot_tax,
  aes(x=x,y=percent,fill=host_range_class)
) +
  geom_col(width=.78,colour="white",linewidth=.25) +
  geom_text(
    aes(
      label=ifelse(
        observed_count>0,
        paste0(
          observed_count,
          " (",
          sprintf("%.1f",percent),
          "%)"
        ),
        ""
      )
    ),
    position=position_stack(vjust=.5),
    size=3
  ) +
  scale_y_continuous(
    limits=c(0,100),
    breaks=seq(0,100,20),
    expand=c(0,0),
    labels=function(x) paste0(x,"%")
  ) +
  labs(
    x="Phage taxonomy",
    y="Host-range composition (%)",
    fill="Host-range class"
  ) +
  theme_classic(base_size=12) +
  theme(
    axis.text.x=element_text(angle=40,hjust=1,vjust=1),
    axis.title=element_text(face="bold"),
    legend.title=element_text(face="bold")
  )

ggsave(
  file.path(
    OUTDIR,
    "taxonomy_host_range_complete_case.svg"
  ),
  p2,
  width=max(
    10,
    min(22,4+.45*nlevels(analysis$taxonomy))
  ),
  height=7,
  units="in",
  device="svg",
  bg="white"
)

# ================================================================
# FIGURE 3: PERMUTATION NULL DISTRIBUTION
# ================================================================

perm_df <- data.frame(statistic=perm_stats[is.finite(perm_stats)])

p3 <- ggplot(
  perm_df,
  aes(x=statistic)
) +
  geom_histogram(
    bins=70,
    linewidth=.2
  ) +
  geom_vline(
    xintercept=observed_stratified_stat,
    linewidth=1
  ) +
  labs(
    x="Taxonomy-stratified association statistic",
    y="Permutation frequency"
  ) +
  theme_classic(base_size=12) +
  theme(
    axis.title=element_text(face="bold")
  )

ggsave(
  file.path(
    OUTDIR,
    "taxonomy_adjusted_permutation_null.svg"
  ),
  p3,
  width=8,
  height=6,
  units="in",
  device="svg",
  bg="white"
)

# ================================================================
# FINAL REPORT
# ================================================================

cat("\n")
cat("============================================================\n")
cat(" TAXONOMY-ADJUSTED HOST-RANGE ASSOCIATION\n")
cat("============================================================\n\n")

cat("Total host-range records: ",n_total,"\n",sep="")
cat("Zero robust species count discarded: ",n_zero,"\n",sep="")
cat("Complete cases: ",nrow(analysis),"\n",sep="")
cat("Isolation-source groups: ",nlevels(analysis$isolation_source),"\n",sep="")
cat("Taxonomy groups: ",nlevels(analysis$taxonomy),"\n",sep="")
cat("Taxonomy strata used in adjusted permutation: ",
    length(strata),"\n\n",sep="")

cat("HOST-RANGE DISTRIBUTION\n")
cat("-----------------------\n")
print(table(analysis$host_range_class))

cat("\nUNADJUSTED ISOLATION SOURCE\n")
cat("---------------------------\n")
cat("Pearson chi-square = ",
    sprintf("%.4f",as.numeric(chi_source$statistic)),"\n",sep="")
cat("df = ",as.numeric(chi_source$parameter),"\n",sep="")
cat("P = ",fmtp(chi_source$p.value),"\n",sep="")
cat("Cramer's V = ",sprintf("%.4f",V_source),"\n",sep="")

cat("\nUNADJUSTED TAXONOMY\n")
cat("-------------------\n")
cat("Pearson chi-square = ",
    sprintf("%.4f",as.numeric(chi_tax$statistic)),"\n",sep="")
cat("df = ",as.numeric(chi_tax$parameter),"\n",sep="")
cat("P = ",fmtp(chi_tax$p.value),"\n",sep="")
cat("Cramer's V = ",sprintf("%.4f",V_tax),"\n",sep="")

cat("\nPRIMARY TEST: ISOLATION SOURCE AFTER TAXONOMY\n")
cat("----------------------------------------------\n")
cat("Observed stratified statistic = ",
    sprintf("%.6f",observed_stratified_stat),"\n",sep="")
cat("Permutation replicates = ",N_PERM,"\n",sep="")
cat("Extreme permutations = ",extreme,"\n",sep="")
cat("Permutation P = ",fmtp(p_adjusted),"\n",sep="")
cat("Adjusted effect size = ",
    sprintf("%.4f",adjusted_effect),"\n",sep="")

if(p_adjusted<0.05){
  cat("\nRESULT: Isolation source retains a significant association\n")
  cat("with host range after conditioning on phage taxonomy.\n")
} else {
  cat("\nRESULT: Isolation source does not retain a significant\n")
  cat("association with host range after conditioning on phage taxonomy.\n")
}

cat("\nTOP UNADJUSTED TAXONOMY RESIDUALS\n")
cat("---------------------------------\n")
print(
  tax_res %>%
    select(
      taxonomy,
      host_range_class,
      observed_count,
      expected_count,
      standardized_residual,
      BH_adjusted_p,
      interpretation
    ) %>%
    head(20),
  row.names=FALSE
)

cat("\nOUTPUT DIRECTORY\n")
cat("----------------\n")
cat(normalizePath(OUTDIR,mustWork=FALSE),"\n")

cat("\n============================================================\n")
cat("Analysis completed successfully.\n")
cat("============================================================\n")
