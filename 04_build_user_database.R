# ============================================================
# 04_build_user_database.R
# ============================================================
#
# Purpose:
#   Build the final user-facing WPFG database.
#
# Inputs:
#   - WPFG classification analysis
#   - genus-level analysis
#
# Output:
#   WPFG_user_database_VerX.xlsx
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
# 2. CHECK RELEASE SETTINGS
# ============================================================

required_objects <- c(
  "analysis_file",
  "genus_file",
  "user_file",
  "version_label",
  "input_file",
  "apc_data_type",
  "apc_version"
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
# 3. READ MAIN ANALYSIS
# ============================================================

taxon_master <- read_excel(
  analysis_file,
  sheet = "taxon_master"
)

master_records <- read_excel(
  analysis_file,
  sheet = "master_records"
)

evidence_records <- read_excel(
  analysis_file,
  sheet = "evidence_records"
)

species_source_7 <- read_excel(
  analysis_file,
  sheet = "species_source_7"
)

species_source_10 <- read_excel(
  analysis_file,
  sheet = "species_source_10"
)

classification_7 <- read_excel(
  analysis_file,
  sheet = "classification_7"
)

classification_10 <- read_excel(
  analysis_file,
  sheet = "classification_10"
)

fuzzy_7_long <- read_excel(
  analysis_file,
  sheet = "fuzzy_7_long"
)

fuzzy_10_long <- read_excel(
  analysis_file,
  sheet = "fuzzy_10_long"
)


# ============================================================
# 4. READ GENUS ANALYSIS
# ============================================================

genus_summary <- read_excel(
  genus_file,
  sheet = "genus_summary"
)

historical_genus_assignments <-
  read_excel(
    genus_file,
    sheet = "historical_genus_assignments"
  )


# ============================================================
# 5. STANDARDISE IDS
# ============================================================

taxon_master <- taxon_master %>%
  mutate(
    taxon_id =
      as.character(
        taxon_id
      )
  )

master_records <- master_records %>%
  mutate(
    taxon_id =
      as.character(
        taxon_id
      )
  )

classification_7 <- classification_7 %>%
  mutate(
    taxon_id =
      as.character(
        taxon_id
      )
  )

classification_10 <- classification_10 %>%
  mutate(
    taxon_id =
      as.character(
        taxon_id
      )
  )

fuzzy_7_long <- fuzzy_7_long %>%
  mutate(
    taxon_id =
      as.character(
        taxon_id
      )
  )

fuzzy_10_long <- fuzzy_10_long %>%
  mutate(
    taxon_id =
      as.character(
        taxon_id
      )
  )


# ============================================================
# 6. VALID WPFG LEVELS
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
# 7. VALIDATION GATES
# ============================================================

invalid_7 <- setdiff(
  unique(
    fuzzy_7_long$terminal_WPFG
  ),
  WPFG_7_levels
)

if (length(invalid_7) > 0) {
  
  stop(
    "7-level fuzzy output contains invalid levels: ",
    paste(
      invalid_7,
      collapse = ", "
    )
  )
}


invalid_10 <- setdiff(
  unique(
    fuzzy_10_long$terminal_WPFG
  ),
  WPFG_10_levels
)

if (length(invalid_10) > 0) {
  
  stop(
    "10-level fuzzy output contains invalid levels: ",
    paste(
      invalid_10,
      collapse = ", "
    )
  )
}


# ============================================================
# 8. USER TAXON IDENTITY
# ============================================================
#
# Family and genus are retained as APC taxonomic information.
# They play no role in species-level WPFG inference.
#
# ============================================================

user_identity <- taxon_master %>%
  
  transmute(
    
    taxon_id,
    
    original_names =
      original_name,
    
    recommended_name,
    
    accepted_name,
    
    suggested_name,
    
    scientific_name,
    
    family,
    
    genus,
    
    taxon_rank,
    
    taxonomic_status
  )


# ============================================================
# 9. FUZZY 7
# ============================================================

fuzzy_7_wide <- fuzzy_7_long %>%
  
  filter(
    terminal_WPFG %in%
      WPFG_7_levels
  ) %>%
  
  select(
    taxon_id,
    terminal_WPFG,
    membership
  ) %>%
  
  pivot_wider(
    
    names_from =
      terminal_WPFG,
    
    values_from =
      membership,
    
    values_fill =
      0,
    
    names_prefix =
      "P_7_"
  )


# ============================================================
# 10. FUZZY 10
# ============================================================

fuzzy_10_wide <- fuzzy_10_long %>%
  
  filter(
    terminal_WPFG %in%
      WPFG_10_levels
  ) %>%
  
  select(
    taxon_id,
    terminal_WPFG,
    membership
  ) %>%
  
  pivot_wider(
    
    names_from =
      terminal_WPFG,
    
    values_from =
      membership,
    
    values_fill =
      0,
    
    names_prefix =
      "P_10_"
  )


# ============================================================
# 11. 7-LEVEL USER RESULTS
# ============================================================

user_7 <- classification_7 %>%
  
  transmute(
    
    taxon_id,
    
    best_guess_7 =
      most_probable_WPFG,
    
    support_7 =
      support_best,
    
    second_guess_7 =
      second_WPFG,
    
    support_second_7 =
      support_second,
    
    margin_7 =
      margin_best_vs_second,
    
    n_sources_7 =
      n_source_assignments,
    
    effective_evidence_7 =
      effective_evidence
  )


# ============================================================
# 12. 10-LEVEL USER RESULTS
# ============================================================

user_10 <- classification_10 %>%
  
  transmute(
    
    taxon_id,
    
    best_guess_10 =
      most_probable_WPFG,
    
    support_10 =
      support_best,
    
    second_guess_10 =
      second_WPFG,
    
    support_second_10 =
      support_second,
    
    margin_10 =
      margin_best_vs_second,
    
    n_sources_10 =
      n_source_assignments,
    
    effective_evidence_10 =
      effective_evidence
  )


# ============================================================
# 13. MAIN DECISION DATABASE
# ============================================================

WPFG_decision_database <-
  user_identity %>%
  
  left_join(
    user_7,
    by = "taxon_id"
  ) %>%
  
  left_join(
    user_10,
    by = "taxon_id"
  ) %>%
  
  left_join(
    fuzzy_7_wide,
    by = "taxon_id"
  ) %>%
  
  left_join(
    fuzzy_10_wide,
    by = "taxon_id"
  )


WPFG_decision_database <-
  WPFG_decision_database %>%
  
  select(
    
    # --------------------------------------------------------
    # Names / taxonomy
    # --------------------------------------------------------
    
    taxon_id,
    
    original_names,
    
    recommended_name,
    
    accepted_name,
    
    suggested_name,
    
    scientific_name,
    
    family,
    
    genus,
    
    taxon_rank,
    
    taxonomic_status,
    
    
    # --------------------------------------------------------
    # Seven-level evidence
    # --------------------------------------------------------
    
    best_guess_7,
    
    support_7,
    
    second_guess_7,
    
    support_second_7,
    
    margin_7,
    
    n_sources_7,
    
    effective_evidence_7,
    
    
    # --------------------------------------------------------
    # Ten-level evidence
    # --------------------------------------------------------
    
    best_guess_10,
    
    support_10,
    
    second_guess_10,
    
    support_second_10,
    
    margin_10,
    
    n_sources_10,
    
    effective_evidence_10,
    
    
    # --------------------------------------------------------
    # Fuzzy membership
    # --------------------------------------------------------
    
    starts_with(
      "P_"
    )
  )


# ============================================================
# 14. NAME LOOKUP
# ============================================================
#
# The name_lookup sheet is the intended entry point for users.
#
# original_name includes:
#   - names occurring in the historical source datasets
#   - recommended_name entries added below
#
# Once a candidate name has been located, the accompanying
# APC taxonomic fields should be used to confirm identity.
#
# taxon_id then provides the link to WPFG_decision_database.
#
# ============================================================

name_lookup <- master_records %>%
  
  select(
    taxon_id,
    original_name
  ) %>%
  
  filter(
    !is.na(original_name),
    original_name != ""
  ) %>%
  
  distinct()


recommended_name_lookup <-
  WPFG_decision_database %>%
  
  transmute(
    
    taxon_id,
    
    original_name =
      recommended_name
  ) %>%
  
  filter(
    !is.na(original_name),
    original_name != ""
  )


name_lookup <-
  bind_rows(
    name_lookup,
    recommended_name_lookup
  ) %>%
  
  distinct(
    original_name,
    taxon_id
  ) %>%
  
  left_join(
    
    WPFG_decision_database %>%
      
      select(
        
        taxon_id,
        
        recommended_name,
        
        accepted_name,
        
        suggested_name,
        
        scientific_name,
        
        family,
        
        genus,
        
        taxon_rank,
        
        taxonomic_status,
        
        best_guess_7,
        
        support_7,
        
        second_guess_7,
        
        support_second_7,
        
        margin_7,
        
        n_sources_7,
        
        effective_evidence_7,
        
        best_guess_10,
        
        support_10,
        
        second_guess_10,
        
        support_second_10,
        
        margin_10,
        
        n_sources_10,
        
        effective_evidence_10
      ),
    
    by =
      "taxon_id"
  ) %>%
  
  arrange(
    original_name,
    recommended_name
  )


# ============================================================
# 15. SPECIES EVIDENCE
# ============================================================

species_evidence <-
  evidence_records


# ============================================================
# 16. README: RELEASE METADATA
# ============================================================

readme_version <- tibble(
  
  section = c(
    "Database version",
    "Source file",
    "APC data type",
    "APC version"
  ),
  
  column = c(
    "",
    "",
    "",
    ""
  ),
  
  description = c(
    version_label,
    basename(input_file),
    apc_data_type,
    apc_version
  )
)


# ============================================================
# 17. README: HOW TO USE THIS DATABASE
# ============================================================

readme_usage <- tibble(
  
  section = c(
    
    "How to use this database - overview",
    
    "Step 1 - Find the taxon",
    
    "Step 1 - Search field",
    
    "Step 1 - Confirm the taxon",
    
    "Step 1 - Taxonomic caution",
    
    "Step 2 - Link to WPFG evidence",
    
    "Step 3 - Choose classification resolution",
    
    "Step 4 - Examine support",
    
    "Step 4 - Examine alternatives",
    
    "Step 4 - Examine evidence depth",
    
    "Step 5 - Examine fuzzy memberships",
    
    "Step 6 - Species-level evidence unavailable",
    
    "Step 6 - Interpret genus evidence",
    
    "Step 6 - Genus-level caution",
    
    "Step 7 - Make the final decision",
    
    "Step 7 - Record the decision"
  ),
  
  information = c(
    
    "For a new species, first identify the appropriate taxonomic record in the name_lookup worksheet, use taxon_id to locate the corresponding record in WPFG_decision_database, then assess the WPFG recommendation together with its supporting evidence. If direct species-level evidence is unavailable, genus_summary can provide descriptive decision-support information.",
    
    "Start with the name_lookup worksheet rather than searching WPFG_decision_database directly.",
    
    "Search the original_name field for the taxon name. This lookup field contains names occurring in the historical source datasets and also includes the database recommended_name, allowing historical and currently recommended names to resolve to the same taxon_id where appropriate.",
    
    "After finding a candidate record, confirm that it represents the intended taxon using recommended_name, accepted_name, suggested_name, scientific_name, family, genus, taxon_rank and taxonomic_status.",
    
    "Family and genus are provided to help confirm taxonomic identity. They should not by themselves be used to assume that a species has the same WPFG as another member of the family or genus.",
    
    "Once the correct record has been identified in name_lookup, note its taxon_id. Use that taxon_id to locate the corresponding row in WPFG_decision_database. taxon_id is the internal link between taxonomic identity and WPFG evidence.",
    
    "WPFG_decision_database provides separate best_guess_7 and best_guess_10 recommendations. The seven-level and ten-level systems are alternative classification resolutions; neither is inherently preferred. Select the resolution appropriate to the study.",
    
    "Interpret best_guess_7 or best_guess_10 together with support_7 or support_10. Support is the proportion of dependence-adjusted historical evidence supporting the best guess.",
    
    "second_guess_7 or second_guess_10 identifies the strongest alternative classification. support_second and margin show how strongly the best guess is separated from that alternative. Small margins indicate greater ambiguity.",
    
    "n_sources gives the number of source lists contributing species-level evidence, while effective_evidence gives the sum of dependence-adjusted assignment weights. A high support value based on one source is not equivalent to the same support based on several sources or greater effective evidence.",
    
    "Fields beginning P_7_ and P_10_ show the full fuzzy distribution of dependence-adjusted evidence among WPFGs. They are particularly useful where evidence is divided among plausible alternatives and a single best guess would conceal meaningful disagreement.",
    
    "If the species has no direct species-level WPFG recommendation, consult the genus_summary worksheet.",
    
    "For each classification resolution, genus_summary reports the dominant WPFG among classified species in the genus, the proportion of classified species whose species-level best guess is that dominant WPFG, and the number of classified species contributing to that estimate. Both the proportion and its evidence base should be considered.",
    
    "Genus-level information is descriptive decision-support evidence, not an automatic species-level WPFG assignment. Users should combine it with relevant biological knowledge and the requirements of their study.",
    
    "The database provides evidence to support WPFG assignment; the final classification remains a user decision.",
    
    "Where the database is used to classify new taxa, retain enough information to reproduce the decision, including the taxonomic name used, matched taxon_id, classification resolution, selected WPFG and database version."
  )
)


# ============================================================
# 18. README: GENERAL INFORMATION
# ============================================================

readme_intro <- tibble(
  
  section = c(
    
    "Purpose",
    
    "7-level system",
    
    "10-level system",
    
    "Support",
    
    "Sources",
    
    "Effective evidence",
    
    "Second guess and margin",
    
    "Fuzzy memberships",
    
    "Taxonomic names",
    
    "Genus-level information",
    
    "Important note"
  ),
  
  information = c(
    
    "This workbook provides taxonomic name harmonisation and evidence-based WPFG classifications derived from multiple historical source lists.",
    
    "The seven-level system comprises Tdr, Tda, ATl, ATe, ARp, ARf and S. Historical Sr, Sk and Se records are represented as S, and ATw is represented as ATe.",
    
    "The ten-level system comprises Tdr, Tda, ATw, ATl, ATe, ARp, ARf, Sr, Sk and Se. It uses only sources from the expanded classification period.",
    
    "Support is the proportion of dependence-adjusted historical evidence supporting the best guess. Higher values indicate stronger agreement among the available evidence.",
    
    "The number of sources is the number of distinct source lists contributing species-level evidence. A result based on one source should be interpreted differently from a result supported by several sources.",
    
    "Effective evidence is the total dependence-adjusted evidence contributing to a classification and is calculated as the sum of assignment-chain weights. Repeated unchanged classifications receive reduced weight because successive lists may have inherited earlier classifications rather than independently reassessing the taxon.",
    
    "The second guess is the strongest alternative WPFG. The margin is the difference between support for the best and second guesses. Small margins indicate greater ambiguity between the leading classifications.",
    
    "Fuzzy membership fields show the full distribution of dependence-adjusted evidence across WPFGs. They sum to 1 within each classification system.",
    
    "Names were harmonised using APCalign. accepted_name, suggested_name, scientific_name, family and genus are retained as taxonomic information from the APC alignment. Original source names are retained for traceability.",
    
    "The genus_summary sheet describes WPFG composition among classified species within each genus. The proportion_species_dominant fields show the proportion of classified species whose best guess matches the dominant WPFG for the genus, while n_species gives the evidence base for that proportion. Genus-level information is descriptive and does not automatically assign a WPFG to an unclassified species.",
    
    "The seven-level and ten-level systems are alternative resolutions. Neither is inherently preferred. Users should select the system appropriate to their study."
  )
)


# ============================================================
# 19. README: MAIN DATABASE COLUMNS
# ============================================================

readme_columns <- tibble(
  
  column = c(
    
    "taxon_id",
    
    "original_names",
    
    "recommended_name",
    
    "accepted_name",
    
    "suggested_name",
    
    "scientific_name",
    
    "family",
    
    "genus",
    
    "taxon_rank",
    
    "taxonomic_status",
    
    "best_guess_7",
    
    "support_7",
    
    "second_guess_7",
    
    "support_second_7",
    
    "margin_7",
    
    "n_sources_7",
    
    "effective_evidence_7",
    
    "best_guess_10",
    
    "support_10",
    
    "second_guess_10",
    
    "support_second_10",
    
    "margin_10",
    
    "n_sources_10",
    
    "effective_evidence_10"
  ),
  
  description = c(
    
    "Stable internal identifier for the resolved taxon. Use this field to link a record in name_lookup to the corresponding record in WPFG_decision_database.",
    
    "All original names from source lists that resolved to this taxon. Multiple names are concatenated with ' | '.",
    
    "Preferred name for general use.",
    
    "APC accepted name, where available.",
    
    "APC suggested name, retained even where it is not accepted.",
    
    "Scientific name returned by APCalign.",
    
    "Family from the APC taxonomic hierarchy. Where APC does not provide an unambiguous family assignment, the field is left missing.",
    
    "Genus from the APC taxonomic hierarchy.",
    
    "Taxonomic rank returned by APC.",
    
    "Taxonomic status returned by APC.",
    
    "WPFG with greatest weighted support under the seven-level system.",
    
    "Proportion of dependence-adjusted evidence supporting best_guess_7.",
    
    "Second-most strongly supported seven-level WPFG.",
    
    "Support for second_guess_7.",
    
    "Difference between support_7 and support_second_7.",
    
    "Number of distinct sources contributing a seven-level classification.",
    
    "Total dependence-adjusted evidence contributing to the seven-level classification; calculated as the sum of assignment-chain weights.",
    
    "WPFG with greatest weighted support under the ten-level system.",
    
    "Proportion of dependence-adjusted evidence supporting best_guess_10.",
    
    "Second-most strongly supported ten-level WPFG.",
    
    "Support for second_guess_10.",
    
    "Difference between support_10 and support_second_10.",
    
    "Number of distinct sources contributing a ten-level classification.",
    
    "Total dependence-adjusted evidence contributing to the ten-level classification; calculated as the sum of assignment-chain weights."
  )
)


# ============================================================
# 20. README: NAME LOOKUP FIELDS
# ============================================================

name_lookup_readme <- tibble(
  
  column = c(
    
    "original_name",
    
    "taxon_id",
    
    "recommended_name",
    
    "accepted_name",
    
    "suggested_name",
    
    "scientific_name",
    
    "family",
    
    "genus",
    
    "taxon_rank",
    
    "taxonomic_status"
  ),
  
  description = c(
    
    "Primary field to search when looking up a taxon. Contains names occurring in the historical source datasets and recommended_name entries added by the database workflow.",
    
    "Stable internal identifier linking the lookup record to the corresponding taxon in WPFG_decision_database.",
    
    "Preferred name for general use. Use to help confirm that the lookup record represents the intended taxon.",
    
    "APC accepted name, where available. Use to help confirm taxonomic identity.",
    
    "Alternative or suggested name returned by APCalign.",
    
    "Scientific name returned by APCalign.",
    
    "Family from the APC taxonomic hierarchy; useful for confirming identity but not for automatically assigning a species-level WPFG.",
    
    "Genus from the APC taxonomic hierarchy; useful for confirming identity but not for automatically assigning a species-level WPFG.",
    
    "Taxonomic rank returned by APC.",
    
    "Taxonomic status returned by APC."
  )
)


# ============================================================
# 21. README: FUZZY FIELDS
# ============================================================

fuzzy_readme <- tibble(
  
  column = c(
    
    "P_7_Tdr",
    "P_7_Tda",
    "P_7_ATl",
    "P_7_ATe",
    "P_7_ARp",
    "P_7_ARf",
    "P_7_S",
    
    "P_10_Tdr",
    "P_10_Tda",
    "P_10_ATw",
    "P_10_ATl",
    "P_10_ATe",
    "P_10_ARp",
    "P_10_ARf",
    "P_10_Sr",
    "P_10_Sk",
    "P_10_Se"
  ),
  
  description = c(
    
    "Fuzzy membership for Tdr under the seven-level system.",
    
    "Fuzzy membership for Tda under the seven-level system.",
    
    "Fuzzy membership for ATl under the seven-level system.",
    
    "Fuzzy membership for ATe under the seven-level system.",
    
    "Fuzzy membership for ARp under the seven-level system.",
    
    "Fuzzy membership for ARf under the seven-level system.",
    
    "Fuzzy membership for S under the seven-level system.",
    
    "Fuzzy membership for Tdr under the ten-level system.",
    
    "Fuzzy membership for Tda under the ten-level system.",
    
    "Fuzzy membership for ATw under the ten-level system.",
    
    "Fuzzy membership for ATl under the ten-level system.",
    
    "Fuzzy membership for ATe under the ten-level system.",
    
    "Fuzzy membership for ARp under the ten-level system.",
    
    "Fuzzy membership for ARf under the ten-level system.",
    
    "Fuzzy membership for Sr under the ten-level system.",
    
    "Fuzzy membership for Sk under the ten-level system.",
    
    "Fuzzy membership for Se under the ten-level system."
  )
)


# ============================================================
# 22. README: GENUS SUMMARY FIELDS
# ============================================================

genus_readme <- tibble(
  
  column = c(
    
    "n_species_7",
    
    "dominant_WPFG_7",
    
    "proportion_species_dominant_7",
    
    "n_species_best_guess_dominant_7",
    
    "dominant_membership_7",
    
    "second_WPFG_7",
    
    "second_membership_7",
    
    "margin_7",
    
    "n_WPFG_7",
    
    "n_species_10",
    
    "dominant_WPFG_10",
    
    "proportion_species_dominant_10",
    
    "n_species_best_guess_dominant_10",
    
    "dominant_membership_10",
    
    "second_WPFG_10",
    
    "second_membership_10",
    
    "margin_10",
    
    "n_WPFG_10",
    
    "reported_genus_WPFG",
    
    "n_genus_assignments",
    
    "n_genus_sources"
  ),
  
  description = c(
    
    "Number of species in the genus with a seven-level species classification.",
    
    "WPFG receiving the greatest average fuzzy membership across classified species under the seven-level system.",
    
    "Proportion of classified species in the genus whose seven-level species best guess is the genus-dominant WPFG.",
    
    "Number of classified species whose seven-level species best guess is the genus-dominant WPFG.",
    
    "Average fuzzy membership of the dominant seven-level WPFG across classified species in the genus.",
    
    "Second-most strongly supported seven-level WPFG within the genus.",
    
    "Average fuzzy membership of the second-ranked seven-level WPFG across classified species in the genus.",
    
    "Difference between dominant and second-ranked seven-level genus membership.",
    
    "Number of seven-level WPFGs with non-zero genus-level membership.",
    
    "Number of species in the genus with a ten-level species classification.",
    
    "WPFG receiving the greatest average fuzzy membership across classified species under the ten-level system.",
    
    "Proportion of classified species in the genus whose ten-level species best guess is the genus-dominant WPFG.",
    
    "Number of classified species whose ten-level species best guess is the genus-dominant WPFG.",
    
    "Average fuzzy membership of the dominant ten-level WPFG across classified species in the genus.",
    
    "Second-most strongly supported ten-level WPFG within the genus.",
    
    "Average fuzzy membership of the second-ranked ten-level WPFG across classified species in the genus.",
    
    "Difference between dominant and second-ranked ten-level genus membership.",
    
    "Number of ten-level WPFGs with non-zero genus-level membership.",
    
    "WPFGs explicitly reported in historical genus-level assignments.",
    
    "Number of historical genus-level WPFG assignments.",
    
    "Number of distinct sources contributing historical genus-level assignments."
  )
)


# ============================================================
# 23. ASSEMBLE README
# ============================================================

readme <- bind_rows(
  
  # ----------------------------------------------------------
  # Release metadata
  # ----------------------------------------------------------
  
  readme_version,
  
  
  # ----------------------------------------------------------
  # Detailed usage instructions
  # ----------------------------------------------------------
  
  tibble(
    section = "",
    column = "",
    description = ""
  ),
  
  readme_usage %>%
    transmute(
      section,
      column = "",
      description =
        information
    ),
  
  
  # ----------------------------------------------------------
  # General information
  # ----------------------------------------------------------
  
  tibble(
    section = "",
    column = "",
    description = ""
  ),
  
  readme_intro %>%
    transmute(
      section,
      column = "",
      description =
        information
    ),
  
  
  # ----------------------------------------------------------
  # name_lookup documentation
  # ----------------------------------------------------------
  
  tibble(
    section = "",
    column = "",
    description = ""
  ),
  
  name_lookup_readme %>%
    mutate(
      section =
        "name_lookup"
    ) %>%
    select(
      section,
      column,
      description
    ),
  
  
  # ----------------------------------------------------------
  # Main database documentation
  # ----------------------------------------------------------
  
  tibble(
    section = "",
    column = "",
    description = ""
  ),
  
  readme_columns %>%
    mutate(
      section =
        "WPFG_decision_database"
    ) %>%
    select(
      section,
      column,
      description
    ),
  
  
  # ----------------------------------------------------------
  # Fuzzy membership documentation
  # ----------------------------------------------------------
  
  fuzzy_readme %>%
    mutate(
      section =
        "Fuzzy membership fields"
    ) %>%
    select(
      section,
      column,
      description
    ),
  
  
  # ----------------------------------------------------------
  # Genus summary documentation
  # ----------------------------------------------------------
  
  genus_readme %>%
    mutate(
      section =
        "genus_summary"
    ) %>%
    select(
      section,
      column,
      description
    )
)


# ============================================================
# 24. FINAL EXPORT
# ============================================================

write_xlsx(
  
  list(
    
    README =
      readme,
    
    WPFG_decision_database =
      WPFG_decision_database,
    
    name_lookup =
      name_lookup,
    
    genus_summary =
      genus_summary,
    
    species_evidence =
      species_evidence,
    
    species_source_7 =
      species_source_7,
    
    species_source_10 =
      species_source_10
    
  ),
  
  user_file
)


# ============================================================
# 25. QA
# ============================================================

cat(
  "\n============================================\n",
  "WPFG USER DATABASE COMPLETE\n",
  "============================================\n\n"
)

cat(
  "Version: ",
  version_label,
  "\n"
)

cat(
  "APC: ",
  apc_data_type,
  " ",
  apc_version,
  "\n"
)

cat(
  "Resolved taxa: ",
  nrow(
    WPFG_decision_database
  ),
  "\n"
)

cat(
  "7-level best guesses: ",
  sum(
    !is.na(
      WPFG_decision_database$best_guess_7
    )
  ),
  "\n"
)

cat(
  "10-level best guesses: ",
  sum(
    !is.na(
      WPFG_decision_database$best_guess_10
    )
  ),
  "\n"
)

cat(
  "Genera in genus summary: ",
  nrow(
    genus_summary
  ),
  "\n"
)

cat(
  "Name lookup rows: ",
  nrow(
    name_lookup
  ),
  "\n"
)

cat(
  "\nOutput:\n",
  user_file,
  "\n"
)