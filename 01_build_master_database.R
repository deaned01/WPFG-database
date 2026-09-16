# ============================================================
# 01_build_master_database.R
# ============================================================
#
# Purpose:
#   Build the taxonomically harmonised master database.
#
# Important:
#   All taxonomic ranks are retained.
#   Only species/infraspecific records are flagged as
#   WPFG evidence eligible.
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(tidyverse)
library(writexl)
library(APCalign)


# ============================================================
# 2. CHECK RELEASE SETTINGS
# ============================================================

required_objects <- c(
  "input_file",
  "master_file",
  "version_label",
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
# 3. READ RAW INPUT
# ============================================================

dat <- read_csv(
  input_file,
  show_col_types = FALSE
) %>%
  mutate(
    record_id = row_number()
  )


required_input <- c(
  "original_name",
  "WPFG",
  "source",
  "chrono"
)

missing_input <- setdiff(
  required_input,
  names(dat)
)

if (length(missing_input) > 0) {
  
  stop(
    "Input is missing required columns: ",
    paste(
      missing_input,
      collapse = ", "
    )
  )
}


raw_records <- dat


# ============================================================
# 4. APC TAXONOMIC LOOKUP
# ============================================================

resources <- load_taxonomic_resources(
  stable_or_current_data =
    apc_data_type,
  version =
    apc_version
)


# ------------------------------------------------------------
# APCalign name resolution
# ------------------------------------------------------------

apc_lookup <- create_taxonomic_update_lookup(
  dat$original_name,
  resources = resources
) %>%
  select(
    original_name,
    aligned_name,
    accepted_name,
    suggested_name,
    taxon_rank,
    taxonomic_dataset,
    taxonomic_status,
    scientific_name,
    aligned_reason,
    update_reason,
    number_of_collapsed_taxa
  ) %>%
  distinct()


# ------------------------------------------------------------
# APC taxonomic hierarchy
#
# The main APC table contains the authoritative family and
# genus associated with the resolved scientific name.
# ------------------------------------------------------------

apc_taxon_lookup <- resources$APC %>%
  select(
    scientific_name,
    family,
    genus
  ) %>%
  filter(
    !is.na(scientific_name),
    scientific_name != ""
  ) %>%
  distinct(
    scientific_name,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# Join family and genus to APCalign result
# ------------------------------------------------------------

apc_lookup <- apc_lookup %>%
  left_join(
    apc_taxon_lookup,
    by = "scientific_name"
  )

# ============================================================
# 5. JOIN APC RESULTS
# ============================================================

resolved <- dat %>%
  left_join(
    apc_lookup,
    by = "original_name"
  )


# ============================================================
# 6. SOURCE CLASSIFICATION SYSTEM
# ============================================================

resolved <- resolved %>%
  mutate(
    
    chrono =
      as.numeric(
        chrono
      ),
    
    source_system =
      case_when(
        
        chrono %in% c(1, 2) ~
          7L,
        
        chrono > 2 ~
          10L,
        
        TRUE ~
          NA_integer_
      )
  )


source_system_check <- resolved %>%
  distinct(
    source,
    chrono,
    source_system
  ) %>%
  count(source) %>%
  filter(n > 1)


if (nrow(source_system_check) > 0) {
  
  stop(
    "A source has multiple chrono/source_system values."
  )
}


# ============================================================
# 7. HELPER FUNCTIONS
# ============================================================

first_non_missing <- function(x) {
  
  x <- x[
    !is.na(x) &
      trimws(
        as.character(x)
      ) != ""
  ]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  as.character(x[1])
}


collapse_unique <- function(x) {
  
  x <- x[
    !is.na(x) &
      trimws(
        as.character(x)
      ) != ""
  ]
  
  x <- unique(
    as.character(x)
  )
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(
    sort(x),
    collapse = " | "
  )
}


# ============================================================
# 8. SPECIES / INFRASPECIFIC EVIDENCE FLAG
# ============================================================

resolved <- resolved %>%
  mutate(
    
    taxon_rank_clean =
      str_to_lower(
        str_trim(
          coalesce(
            taxon_rank,
            ""
          )
        )
      ),
    
    WPFG_evidence_eligible =
      
      str_detect(
        taxon_rank_clean,
        "species"
      ) |
      
      str_detect(
        taxon_rank_clean,
        "subspecies"
      ) |
      
      str_detect(
        taxon_rank_clean,
        "variety"
      ) |
      
      str_detect(
        taxon_rank_clean,
        "forma"
      ) |
      
      str_detect(
        taxon_rank_clean,
        "\\bform\\b"
      ) |
      
      str_detect(
        taxon_rank_clean,
        "infraspecific"
      )
  )


rank_summary <- resolved %>%
  count(
    taxon_rank,
    WPFG_evidence_eligible,
    sort = TRUE
  )


# ============================================================
# 9. TAXON IDENTITY
# ============================================================

resolved <- resolved %>%
  mutate(
    
    taxon_key =
      case_when(
        
        !is.na(accepted_name) &
          accepted_name != "" ~
          accepted_name,
        
        !is.na(aligned_name) &
          aligned_name != "" ~
          aligned_name,
        
        !is.na(scientific_name) &
          scientific_name != "" ~
          scientific_name,
        
        TRUE ~
          original_name
      )
  )


taxon_lookup <- resolved %>%
  distinct(
    taxon_key
  ) %>%
  arrange(
    taxon_key
  ) %>%
  mutate(
    taxon_id =
      sprintf(
        "T%05d",
        row_number()
      )
  )


resolved <- resolved %>%
  left_join(
    taxon_lookup,
    by = "taxon_key"
  )


# ============================================================
# 10. TAXON MASTER
# ============================================================

taxon_master <- resolved %>%
  group_by(
    taxon_id,
    taxon_key
  ) %>%
  summarise(
    
    accepted_name =
      first_non_missing(
        accepted_name
      ),
    
    scientific_name =
      collapse_unique(
        scientific_name
      ),
    
    suggested_name =
      collapse_unique(
        suggested_name
      ),
    
    aligned_name =
      collapse_unique(
        aligned_name
      ),
    
    family =
      collapse_unique(
        family
      ),
    
    genus =
      collapse_unique(
        genus
      ),
    
    original_name =
      collapse_unique(
        original_name
      ),
    
    taxon_rank =
      collapse_unique(
        taxon_rank
      ),
    
    taxonomic_status =
      collapse_unique(
        taxonomic_status
      ),
    
    taxonomic_dataset =
      collapse_unique(
        taxonomic_dataset
      ),
    
    aligned_reason =
      collapse_unique(
        aligned_reason
      ),
    
    update_reason =
      collapse_unique(
        update_reason
      ),
    
    number_of_collapsed_taxa =
      suppressWarnings(
        max(
          number_of_collapsed_taxa,
          na.rm = TRUE
        )
      ),
    
    n_total_records =
      n(),
    
    n_sources_total =
      n_distinct(
        source
      ),
    
    n_species_level_records =
      sum(
        WPFG_evidence_eligible,
        na.rm = TRUE
      ),
    
    n_species_level_sources =
      n_distinct(
        source[
          WPFG_evidence_eligible
        ]
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    
    number_of_collapsed_taxa =
      ifelse(
        is.infinite(
          number_of_collapsed_taxa
        ),
        NA,
        number_of_collapsed_taxa
      )
  )


# ============================================================
# 11. RECOMMENDED NAME
# ============================================================

taxon_master <- taxon_master %>%
  mutate(
    
    recommended_name =
      case_when(
        
        !is.na(accepted_name) &
          accepted_name != "" ~
          accepted_name,
        
        !is.na(suggested_name) &
          suggested_name != "" ~
          suggested_name,
        
        !is.na(scientific_name) &
          scientific_name != "" ~
          scientific_name,
        
        TRUE ~
          original_name
      )
  ) %>%
  relocate(
    
    taxon_id,
    original_name,
    recommended_name,
    accepted_name,
    suggested_name,
    scientific_name,
    family,
    genus
  )


# ============================================================
# 12. WPFG HIERARCHY
# ============================================================

group_priority <- tibble(
  
  WPFG = c(
    "Sk",
    "Sr",
    "Se",
    "S",
    "ARf",
    "ARp",
    "ATe",
    "ATl",
    "ATw",
    "Tda",
    "Tdr"
  ),
  
  priority = c(
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11
  )
)


# ============================================================
# 13. WITHIN-SOURCE RESOLUTION
# ============================================================

resolve_WPFG <- function(
    df,
    priority_table
) {
  
  df <- df %>%
    filter(
      !is.na(WPFG),
      WPFG != ""
    )
  
  
  if (nrow(df) == 0) {
    
    return(
      tibble(
        selected_WPFG = NA_character_,
        selection_method = "No WPFG",
        conflict = FALSE,
        n_source_records = 0
      )
    )
  }
  
  
  counts <- df %>%
    count(
      WPFG,
      name = "n"
    )
  
  
  if (nrow(counts) == 1) {
    
    return(
      tibble(
        selected_WPFG =
          counts$WPFG[1],
        selection_method =
          "Single group",
        conflict =
          FALSE,
        n_source_records =
          nrow(df)
      )
    )
  }
  
  
  max_n <- max(
    counts$n
  )
  
  
  majority_groups <- counts %>%
    filter(
      n == max_n
    )
  
  
  if (nrow(majority_groups) == 1) {
    
    return(
      tibble(
        selected_WPFG =
          majority_groups$WPFG[1],
        selection_method =
          "Majority",
        conflict =
          TRUE,
        n_source_records =
          nrow(df)
      )
    )
  }
  
  
  ranked <- counts %>%
    left_join(
      priority_table,
      by = "WPFG"
    ) %>%
    mutate(
      priority =
        coalesce(
          priority,
          999
        )
    ) %>%
    arrange(
      priority,
      WPFG
    )
  
  
  tibble(
    selected_WPFG =
      ranked$WPFG[1],
    
    selection_method =
      "Hierarchy",
    
    conflict =
      TRUE,
    
    n_source_records =
      nrow(df)
  )
}


# ============================================================
# 14. BUILD SOURCE ASSIGNMENTS
# ============================================================

build_source_assignments <- function(
    data,
    wpfg_column
) {
  
  if (!wpfg_column %in% names(data)) {
    
    stop(
      "Column '",
      wpfg_column,
      "' not found."
    )
  }
  
  
  source_assignments <- data %>%
    
    group_by(
      taxon_id,
      source
    ) %>%
    
    nest() %>%
    
    mutate(
      
      resolution =
        map(
          
          data,
          
          ~ resolve_WPFG(
            
            tibble(
              WPFG =
                .x[[wpfg_column]]
            ),
            
            group_priority
          )
        )
    ) %>%
    
    select(
      taxon_id,
      source,
      resolution
    ) %>%
    
    unnest(
      resolution
    ) %>%
    
    ungroup()
  
  
  source_details <- data %>%
    group_by(
      taxon_id,
      source
    ) %>%
    summarise(
      
      chrono =
        first_non_missing(
          chrono
        ),
      
      source_system =
        first_non_missing(
          source_system
        ),
      
      original_names =
        collapse_unique(
          original_name
        ),
      
      aligned_names =
        collapse_unique(
          aligned_name
        ),
      
      suggested_names =
        collapse_unique(
          suggested_name
        ),
      
      source_WPFG =
        collapse_unique(
          .data[[wpfg_column]]
        ),
      
      taxon_ranks =
        collapse_unique(
          taxon_rank
        ),
      
      .groups = "drop"
    )
  
  
  source_assignments %>%
    left_join(
      source_details,
      by = c(
        "taxon_id",
        "source"
      )
    ) %>%
    select(
      taxon_id,
      source,
      chrono,
      source_system,
      original_names,
      aligned_names,
      suggested_names,
      taxon_ranks,
      source_WPFG,
      selected_WPFG,
      selection_method,
      conflict,
      n_source_records
    )
}


source_assignments <-
  build_source_assignments(
    resolved,
    "WPFG"
  )


duplicate_source_keys <-
  source_assignments %>%
  count(
    taxon_id,
    source
  ) %>%
  filter(
    n > 1
  )


if (nrow(duplicate_source_keys) > 0) {
  
  stop(
    "Duplicate taxon_id x source records remain."
  )
}


conflict_report <-
  source_assignments %>%
  filter(
    conflict
  ) %>%
  arrange(
    source,
    taxon_id
  )


# ============================================================
# 15. MASTER DATABASE
# ============================================================

combined_WPFG <- source_assignments %>%
  group_by(
    taxon_id
  ) %>%
  summarise(
    
    combined_WPFG =
      collapse_unique(
        selected_WPFG
      ),
    
    sources =
      collapse_unique(
        source
      ),
    
    .groups = "drop"
  )


master_database <-
  taxon_master %>%
  left_join(
    combined_WPFG,
    by = "taxon_id"
  )


wide_WPFG <- source_assignments %>%
  select(
    taxon_id,
    source,
    selected_WPFG
  ) %>%
  pivot_wider(
    names_from =
      source,
    values_from =
      selected_WPFG,
    names_prefix =
      "WPFG_"
  )


master_wide <-
  master_database %>%
  left_join(
    wide_WPFG,
    by = "taxon_id"
  )


# ============================================================
# 16. METADATA
# ============================================================

source_metadata <- resolved %>%
  distinct(
    source,
    chrono,
    source_system
  ) %>%
  arrange(
    chrono,
    source
  )


database_metadata <- tibble(
  
  field = c(
    "database_version",
    "source_file",
    "APC_data_type",
    "APC_version"
  ),
  
  value = c(
    version_label,
    basename(input_file),
    apc_data_type,
    apc_version
  )
)


# ============================================================
# 17. EXPORT
# ============================================================

write_xlsx(
  
  list(
    
    database_metadata =
      database_metadata,
    
    raw_records =
      raw_records,
    
    resolved_records =
      resolved,
    
    taxon_master =
      taxon_master,
    
    source_metadata =
      source_metadata,
    
    rank_summary =
      rank_summary,
    
    source_assignments =
      source_assignments,
    
    conflict_report =
      conflict_report,
    
    master_database =
      master_database,
    
    master_wide =
      master_wide
    
  ),
  
  master_file
)


cat(
  "\nScript 1 complete.\n",
  "Taxa: ",
  n_distinct(
    resolved$taxon_id
  ),
  "\n"
)