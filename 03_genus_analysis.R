# ============================================================
# 03_genus_analysis.R
# ============================================================
#
# Purpose:
#   Describe WPFG conservation within genera using classified
#   species.
#
# Important:
#   Genus-level information is DESCRIPTIVE ONLY.
#   It is not used to infer species-level WPFG classifications.
#
# Outputs:
#   - genus_summary
#   - historical_genus_assignments
#   - genus_fuzzy_7
#   - genus_fuzzy_10
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(dplyr)
library(tidyr)
library(readxl)
library(writexl)


# ============================================================
# 2. CHECK SETTINGS
# ============================================================

required_objects <- c(
  "analysis_file",
  "genus_file"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1)
  )
]

if (length(missing_objects) > 0) {
  
  stop(
    "Missing run_all.R settings: ",
    paste(
      missing_objects,
      collapse = ", "
    )
  )
}


# ============================================================
# 3. READ ANALYSIS
# ============================================================

taxon_master <- read_excel(
  analysis_file,
  sheet = "taxon_master"
)

master_records <- read_excel(
  analysis_file,
  sheet = "master_records"
)

fuzzy_7_long <- read_excel(
  analysis_file,
  sheet = "fuzzy_7_long"
)

fuzzy_10_long <- read_excel(
  analysis_file,
  sheet = "fuzzy_10_long"
)

classification_7 <- read_excel(
  analysis_file,
  sheet = "classification_7"
)

classification_10 <- read_excel(
  analysis_file,
  sheet = "classification_10"
)


# ============================================================
# 4. TYPES
# ============================================================

taxon_master <- taxon_master %>%
  mutate(
    taxon_id =
      as.character(taxon_id)
  )

master_records <- master_records %>%
  mutate(
    taxon_id =
      as.character(taxon_id)
  )

fuzzy_7_long <- fuzzy_7_long %>%
  mutate(
    taxon_id =
      as.character(taxon_id)
  )

fuzzy_10_long <- fuzzy_10_long %>%
  mutate(
    taxon_id =
      as.character(taxon_id)
  )

classification_7 <- classification_7 %>%
  mutate(
    taxon_id =
      as.character(taxon_id)
  )

classification_10 <- classification_10 %>%
  mutate(
    taxon_id =
      as.character(taxon_id)
  )


# ============================================================
# 5. VALID WPFG LEVELS
# ============================================================

WPFG_7_levels <- c(
  "Tdr",
  "Tda",
  "ATl",
  "ATe",
  "ARp",
  "ARf",
  "S"
)

WPFG_10_levels <- c(
  "Tdr",
  "Tda",
  "ATw",
  "ATl",
  "ATe",
  "ARp",
  "ARf",
  "Sr",
  "Sk",
  "Se"
)


# ============================================================
# 6. SPECIES -> GENUS LOOKUP
# ============================================================

species_taxa <- taxon_master %>%
  select(
    taxon_id,
    genus
  ) %>%
  filter(
    !is.na(genus),
    genus != ""
  )


# ============================================================
# 7. SPECIES-DERIVED GENUS FUZZY COMPOSITION
# ============================================================
#
# Each species contributes equally to the genus-level fuzzy
# distribution.
#
# ============================================================

genus_fuzzy_7 <- fuzzy_7_long %>%
  
  filter(
    terminal_WPFG %in%
      WPFG_7_levels
  ) %>%
  
  left_join(
    species_taxa,
    by = "taxon_id"
  ) %>%
  
  filter(
    !is.na(genus),
    genus != ""
  ) %>%
  
  group_by(
    genus,
    terminal_WPFG
  ) %>%
  
  summarise(
    genus_membership =
      mean(
        membership,
        na.rm = TRUE
      ),
    .groups = "drop"
  )


genus_fuzzy_10 <- fuzzy_10_long %>%
  
  filter(
    terminal_WPFG %in%
      WPFG_10_levels
  ) %>%
  
  left_join(
    species_taxa,
    by = "taxon_id"
  ) %>%
  
  filter(
    !is.na(genus),
    genus != ""
  ) %>%
  
  group_by(
    genus,
    terminal_WPFG
  ) %>%
  
  summarise(
    genus_membership =
      mean(
        membership,
        na.rm = TRUE
      ),
    .groups = "drop"
  )


# ============================================================
# 8. WIDE GENUS FUZZY TABLES
# ============================================================

genus_fuzzy_7_wide <- genus_fuzzy_7 %>%
  
  pivot_wider(
    names_from =
      terminal_WPFG,
    values_from =
      genus_membership,
    values_fill =
      0,
    names_prefix =
      "P_genus_7_"
  )


genus_fuzzy_10_wide <- genus_fuzzy_10 %>%
  
  pivot_wider(
    names_from =
      terminal_WPFG,
    values_from =
      genus_membership,
    values_fill =
      0,
    names_prefix =
      "P_genus_10_"
  )


# ============================================================
# 9. GENUS-LEVEL DOMINANT WPFG
# ============================================================

genus_summary_7 <- genus_fuzzy_7 %>%
  
  group_by(
    genus
  ) %>%
  
  arrange(
    desc(genus_membership),
    terminal_WPFG
  ) %>%
  
  mutate(
    rank =
      row_number()
  ) %>%
  
  summarise(
    
    dominant_WPFG_7 =
      terminal_WPFG[
        rank == 1
      ][1],
    
    dominant_membership_7 =
      genus_membership[
        rank == 1
      ][1],
    
    second_WPFG_7 =
      terminal_WPFG[
        rank == 2
      ][1],
    
    second_membership_7 =
      genus_membership[
        rank == 2
      ][1],
    
    margin_7 =
      dominant_membership_7 -
      second_membership_7,
    
    n_WPFG_7 =
      sum(
        genus_membership > 0
      ),
    
    .groups = "drop"
  )


genus_summary_10 <- genus_fuzzy_10 %>%
  
  group_by(
    genus
  ) %>%
  
  arrange(
    desc(genus_membership),
    terminal_WPFG
  ) %>%
  
  mutate(
    rank =
      row_number()
  ) %>%
  
  summarise(
    
    dominant_WPFG_10 =
      terminal_WPFG[
        rank == 1
      ][1],
    
    dominant_membership_10 =
      genus_membership[
        rank == 1
      ][1],
    
    second_WPFG_10 =
      terminal_WPFG[
        rank == 2
      ][1],
    
    second_membership_10 =
      genus_membership[
        rank == 2
      ][1],
    
    margin_10 =
      dominant_membership_10 -
      second_membership_10,
    
    n_WPFG_10 =
      sum(
        genus_membership > 0
      ),
    
    .groups = "drop"
  )


# ============================================================
# 10. SPECIES BEST-GUESS CONSERVATION
# ============================================================
#
# This provides the most intuitive genus-level evidence for
# users:
#
#   "What proportion of classified species in this genus have
#    the genus-dominant WPFG as their species-level best guess?"
#
# This is NOT a species-level recommendation.
#
# ============================================================

species_best_7 <- classification_7 %>%
  
  transmute(
    taxon_id,
    species_best_guess_7 =
      most_probable_WPFG
  ) %>%
  
  left_join(
    species_taxa,
    by = "taxon_id"
  ) %>%
  
  filter(
    !is.na(genus),
    genus != ""
  )


species_best_10 <- classification_10 %>%
  
  transmute(
    taxon_id,
    species_best_guess_10 =
      most_probable_WPFG
  ) %>%
  
  left_join(
    species_taxa,
    by = "taxon_id"
  ) %>%
  
  filter(
    !is.na(genus),
    genus != ""
  )


# ============================================================
# 11. NUMBER OF CLASSIFIED SPECIES PER GENUS
# ============================================================

genus_species_count_7 <- species_best_7 %>%
  
  group_by(
    genus
  ) %>%
  
  summarise(
    n_species_7 =
      n_distinct(
        taxon_id
      ),
    .groups = "drop"
  )


genus_species_count_10 <- species_best_10 %>%
  
  group_by(
    genus
  ) %>%
  
  summarise(
    n_species_10 =
      n_distinct(
        taxon_id
      ),
    .groups = "drop"
  )


# ============================================================
# 12. GENUS-DOMINANT BEST GUESS PROPORTION: 7 LEVEL
# ============================================================

genus_best_7 <- genus_summary_7 %>%
  
  select(
    genus,
    dominant_WPFG_7
  ) %>%
  
  left_join(
    species_best_7,
    by = "genus"
  ) %>%
  
  group_by(
    genus,
    dominant_WPFG_7
  ) %>%
  
  summarise(
    
    n_species_7 =
      n_distinct(
        taxon_id
      ),
    
    n_species_best_guess_dominant_7 =
      sum(
        species_best_guess_7 ==
          dominant_WPFG_7,
        na.rm = TRUE
      ),
    
    proportion_species_dominant_7 =
      if_else(
        n_species_7 > 0,
        n_species_best_guess_dominant_7 /
          n_species_7,
        NA_real_
      ),
    
    .groups = "drop"
  )


# ============================================================
# 13. GENUS-DOMINANT BEST GUESS PROPORTION: 10 LEVEL
# ============================================================

genus_best_10 <- genus_summary_10 %>%
  
  select(
    genus,
    dominant_WPFG_10
  ) %>%
  
  left_join(
    species_best_10,
    by = "genus"
  ) %>%
  
  group_by(
    genus,
    dominant_WPFG_10
  ) %>%
  
  summarise(
    
    n_species_10 =
      n_distinct(
        taxon_id
      ),
    
    n_species_best_guess_dominant_10 =
      sum(
        species_best_guess_10 ==
          dominant_WPFG_10,
        na.rm = TRUE
      ),
    
    proportion_species_dominant_10 =
      if_else(
        n_species_10 > 0,
        n_species_best_guess_dominant_10 /
          n_species_10,
        NA_real_
      ),
    
    .groups = "drop"
  )


# ============================================================
# 14. HISTORICAL GENUS-LEVEL ASSIGNMENTS
# ============================================================
#
# These are descriptive only.
# They are NOT used to infer species-level WPFG.
#
# ============================================================

historical_genus_assignments <-
  master_records %>%
  
  filter(
    !WPFG_evidence_eligible
  ) %>%
  
  filter(
    !is.na(genus),
    genus != "",
    !is.na(WPFG),
    WPFG != ""
  ) %>%
  
  select(
    
    taxon_id,
    genus,
    original_name,
    taxon_rank,
    WPFG,
    source,
    chrono,
    source_system
    
  ) %>%
  
  distinct()


historical_genus_summary <-
  historical_genus_assignments %>%
  
  group_by(
    genus
  ) %>%
  
  summarise(
    
    reported_genus_WPFG =
      paste(
        sort(
          unique(
            WPFG
          )
        ),
        collapse = " | "
      ),
    
    n_genus_assignments =
      n(),
    
    n_genus_sources =
      n_distinct(
        source
      ),
    
    .groups = "drop"
  )


# ============================================================
# 15. COMBINE GENUS SUMMARY
# ============================================================

genus_summary <-
  
  full_join(
    genus_species_count_7,
    genus_summary_7,
    by = "genus"
  ) %>%
  
  full_join(
    genus_best_7 %>%
      select(
        genus,
        n_species_best_guess_dominant_7,
        proportion_species_dominant_7
      ),
    by = "genus"
  ) %>%
  
  full_join(
    genus_species_count_10,
    by = "genus"
  ) %>%
  
  full_join(
    genus_summary_10,
    by = "genus"
  ) %>%
  
  full_join(
    genus_best_10 %>%
      select(
        genus,
        n_species_best_guess_dominant_10,
        proportion_species_dominant_10
      ),
    by = "genus"
  ) %>%
  
  left_join(
    genus_fuzzy_7_wide,
    by = "genus"
  ) %>%
  
  left_join(
    genus_fuzzy_10_wide,
    by = "genus"
  ) %>%
  
  left_join(
    historical_genus_summary,
    by = "genus"
  ) %>%
  
  arrange(
    genus
  )


# ============================================================
# 16. ORDER FINAL GENUS TABLE
# ============================================================

genus_summary <- genus_summary %>%
  
  select(
    
    # --------------------------------------------------------
    # Genus
    # --------------------------------------------------------
    
    genus,
    
    
    # --------------------------------------------------------
    # Practical 7-level evidence
    # --------------------------------------------------------
    
    n_species_7,
    
    dominant_WPFG_7,
    
    proportion_species_dominant_7,
    
    n_species_best_guess_dominant_7,
    
    
    # --------------------------------------------------------
    # Detailed 7-level evidence
    # --------------------------------------------------------
    
    dominant_membership_7,
    
    second_WPFG_7,
    
    second_membership_7,
    
    margin_7,
    
    n_WPFG_7,
    
    
    # --------------------------------------------------------
    # Practical 10-level evidence
    # --------------------------------------------------------
    
    n_species_10,
    
    dominant_WPFG_10,
    
    proportion_species_dominant_10,
    
    n_species_best_guess_dominant_10,
    
    
    # --------------------------------------------------------
    # Detailed 10-level evidence
    # --------------------------------------------------------
    
    dominant_membership_10,
    
    second_WPFG_10,
    
    second_membership_10,
    
    margin_10,
    
    n_WPFG_10,
    
    
    # --------------------------------------------------------
    # Historical genus assignments
    # --------------------------------------------------------
    
    reported_genus_WPFG,
    
    n_genus_assignments,
    
    n_genus_sources,
    
    
    # --------------------------------------------------------
    # Full genus fuzzy distributions
    # --------------------------------------------------------
    
    starts_with(
      "P_genus_"
    )
  )


# ============================================================
# 17. EXPORT
# ============================================================

write_xlsx(
  
  list(
    
    genus_summary =
      genus_summary,
    
    historical_genus_assignments =
      historical_genus_assignments,
    
    genus_fuzzy_7 =
      genus_fuzzy_7,
    
    genus_fuzzy_10 =
      genus_fuzzy_10
    
  ),
  
  genus_file
)


# ============================================================
# 18. QA
# ============================================================

cat(
  "\n============================================\n",
  "GENUS ANALYSIS COMPLETE\n",
  "============================================\n\n"
)

cat(
  "Genera analysed: ",
  nrow(
    genus_summary
  ),
  "\n"
)

cat(
  "7-level genus summaries: ",
  sum(
    !is.na(
      genus_summary$dominant_WPFG_7
    )
  ),
  "\n"
)

cat(
  "10-level genus summaries: ",
  sum(
    !is.na(
      genus_summary$dominant_WPFG_10
    )
  ),
  "\n"
)

cat(
  "\nOutput:\n",
  genus_file,
  "\n"
)