# ============================================================
# 02_wpfg_analysis.R
# ============================================================

library(tidyverse)
library(readxl)
library(writexl)
library(tibble)


# ============================================================
# CHECK SETTINGS
# ============================================================

required_objects <- c(
  "master_file",
  "analysis_file",
  "version_label"
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
# SETTINGS
# ============================================================

alpha_penalty <- 0.2


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


WPFG_children <- list(
  S = c(
    "Sr",
    "Sk",
    "Se"
  )
)


# ============================================================
# READ MASTER
# ============================================================

resolved <- read_excel(
  master_file,
  sheet = "resolved_records"
) %>%
  mutate(
    
    taxon_id =
      as.character(taxon_id),
    
    source =
      as.character(source),
    
    original_name =
      as.character(original_name),
    
    WPFG =
      as.character(WPFG),
    
    taxon_rank =
      as.character(taxon_rank),
    
    chrono =
      as.numeric(chrono)
  )


taxon_master <- read_excel(
  master_file,
  sheet = "taxon_master"
) %>%
  mutate(
    taxon_id =
      as.character(taxon_id)
  )


# ============================================================
# SOURCE METADATA
# ============================================================

source_metadata <- resolved %>%
  distinct(
    source,
    chrono,
    source_system
  )


source_check <- source_metadata %>%
  count(source) %>%
  filter(
    n > 1
  )


if (nrow(source_check) > 0) {
  
  stop(
    "A source has multiple chrono/source_system values."
  )
}


# ============================================================
# SPECIES-LEVEL EVIDENCE ONLY
# ============================================================

evidence_records <- resolved %>%
  filter(
    WPFG_evidence_eligible
  )


# ============================================================
# HELPERS
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
# WPFG HIERARCHY
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
# SOURCE-LEVEL RESOLUTION
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
        selected_WPFG =
          NA_character_,
        
        selection_method =
          "No WPFG",
        
        conflict =
          FALSE,
        
        n_source_records =
          0
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
  
  
  max_n <-
    max(
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


# ============================================================
# 7-LEVEL EVIDENCE
# ============================================================

evidence_7 <- evidence_records %>%
  mutate(
    
    WPFG_analysis =
      case_when(
        
        WPFG %in% c(
          "Sr",
          "Sk",
          "Se",
          "S"
        ) ~
          "S",
        
        WPFG %in% c(
          "ATw",
          "ATe"
        ) ~
          "ATe",
        
        TRUE ~
          WPFG
      )
  )


invalid_evidence_7 <- setdiff(
  unique(
    na.omit(
      evidence_7$WPFG_analysis
    )
  ),
  WPFG_7_levels
)


if (length(invalid_evidence_7) > 0) {
  
  stop(
    "Invalid 7-level evidence values: ",
    paste(
      invalid_evidence_7,
      collapse = ", "
    )
  )
}


species_source_assignments_7 <-
  build_source_assignments(
    evidence_7,
    "WPFG_analysis"
  )


# ============================================================
# 10-LEVEL EVIDENCE
# ============================================================

evidence_10 <- evidence_records %>%
  filter(
    source_system == 10
  ) %>%
  mutate(
    WPFG_analysis =
      WPFG
  )


invalid_evidence_10 <- setdiff(
  unique(
    na.omit(
      evidence_10$WPFG_analysis
    )
  ),
  c(
    WPFG_10_levels,
    "S"
  )
)


if (length(invalid_evidence_10) > 0) {
  
  stop(
    "Invalid 10-level evidence values: ",
    paste(
      invalid_evidence_10,
      collapse = ", "
    )
  )
}


species_source_assignments_10 <-
  build_source_assignments(
    evidence_10,
    "WPFG_analysis"
  )


# ============================================================
# HIERARCHICAL RELATIONSHIPS
# ============================================================

get_terminal_set <- function(
    wpfg,
    system
) {
  
  if (
    is.na(wpfg) ||
    wpfg == ""
  ) {
    return(character(0))
  }
  
  
  if (system == 7) {
    return(wpfg)
  }
  
  
  if (wpfg == "S") {
    
    return(
      WPFG_children$S
    )
  }
  
  
  wpfg
}


classification_relation <- function(
    previous,
    current,
    system
) {
  
  if (
    is.na(previous) ||
    previous == "" ||
    is.na(current) ||
    current == ""
  ) {
    return("initial")
  }
  
  
  if (previous == current) {
    return("repeat")
  }
  
  
  previous_set <- get_terminal_set(
    previous,
    system
  )
  
  
  current_set <- get_terminal_set(
    current,
    system
  )
  
  
  if (
    length(
      intersect(
        previous_set,
        current_set
      )
    ) > 0
  ) {
    
    return(
      "refinement"
    )
  }
  
  
  "change"
}


# ============================================================
# CHAINS
# ============================================================

build_chain <- function(
    data,
    wpfg_column,
    system,
    alpha = 0.2
) {
  
  data %>%
    
    select(
      taxon_id,
      source,
      chrono,
      all_of(wpfg_column)
    ) %>%
    
    rename(
      WPFG =
        all_of(
          wpfg_column
        )
    ) %>%
    
    mutate(
      chrono =
        as.numeric(
          chrono
        )
    ) %>%
    
    filter(
      !is.na(WPFG),
      WPFG != ""
    ) %>%
    
    arrange(
      taxon_id,
      chrono,
      source
    ) %>%
    
    group_by(
      taxon_id
    ) %>%
    
    mutate(
      
      previous_WPFG =
        lag(WPFG),
      
      relation =
        map2_chr(
          previous_WPFG,
          WPFG,
          
          ~ classification_relation(
            .x,
            .y,
            system
          )
        ),
      
      weight =
        case_when(
          
          relation == "initial" ~
            1.0,
          
          relation == "repeat" ~
            alpha,
          
          relation == "refinement" ~
            1.0,
          
          relation == "change" ~
            1.0,
          
          TRUE ~
            NA_real_
        )
    ) %>%
    
    ungroup()
}


chain_7 <- build_chain(
  species_source_assignments_7,
  "selected_WPFG",
  system = 7,
  alpha = alpha_penalty
)


chain_10 <- build_chain(
  species_source_assignments_10,
  "selected_WPFG",
  system = 10,
  alpha = alpha_penalty
)


# ============================================================
# FUZZY EVIDENCE
# ============================================================

fuzzy_evidence_7 <- chain_7 %>%
  transmute(
    
    taxon_id,
    source,
    chrono,
    
    observed_WPFG =
      WPFG,
    
    weight,
    
    terminal_WPFG =
      WPFG,
    
    terminal_share =
      1
  )


fuzzy_evidence_10 <- chain_10 %>%
  
  mutate(
    
    terminal_WPFG =
      map(
        
        WPFG,
        
        ~ get_terminal_set(
          .x,
          system = 10
        )
      )
  ) %>%
  
  unnest_longer(
    terminal_WPFG
  ) %>%
  
  group_by(
    
    taxon_id,
    source,
    chrono,
    observed_WPFG =
      WPFG,
    weight
    
  ) %>%
  
  mutate(
    
    terminal_share =
      1 / n()
  ) %>%
  
  ungroup()


calculate_membership <- function(
    evidence,
    levels
) {
  
  evidence %>%
    
    mutate(
      
      terminal_evidence =
        weight *
        terminal_share
    ) %>%
    
    group_by(
      taxon_id,
      terminal_WPFG
    ) %>%
    
    summarise(
      
      evidence_weight =
        sum(
          terminal_evidence,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    ) %>%
    
    complete(
      
      taxon_id,
      
      terminal_WPFG =
        levels,
      
      fill = list(
        evidence_weight =
          0
      )
    ) %>%
    
    group_by(
      taxon_id
    ) %>%
    
    mutate(
      
      effective_evidence =
        sum(
          evidence_weight
        ),
      
      membership =
        if_else(
          effective_evidence > 0,
          
          evidence_weight /
            effective_evidence,
          
          0
        )
    ) %>%
    
    ungroup()
}


fuzzy_7_long <- calculate_membership(
  fuzzy_evidence_7,
  WPFG_7_levels
)


fuzzy_10_long <- calculate_membership(
  fuzzy_evidence_10,
  WPFG_10_levels
)


# ============================================================
# CLASSIFICATION SUPPORT
# ============================================================

calculate_taxon_support <- function(
    fuzzy_long
) {
  
  rankings <- fuzzy_long %>%
    arrange(
      taxon_id,
      desc(membership),
      terminal_WPFG
    ) %>%
    group_by(
      taxon_id
    ) %>%
    mutate(
      rank =
        row_number()
    ) %>%
    ungroup()
  
  
  best <- rankings %>%
    filter(
      rank == 1
    ) %>%
    transmute(
      
      taxon_id,
      
      most_probable_WPFG =
        terminal_WPFG,
      
      support_best =
        membership
    )
  
  
  second <- rankings %>%
    filter(
      rank == 2
    ) %>%
    transmute(
      
      taxon_id,
      
      second_WPFG =
        terminal_WPFG,
      
      support_second =
        membership
    )
  
  
  best %>%
    left_join(
      second,
      by = "taxon_id"
    ) %>%
    mutate(
      
      support_second =
        replace_na(
          support_second,
          0
        ),
      
      margin_best_vs_second =
        support_best -
        support_second
    )
}


# ============================================================
# CONCENTRATION
# ============================================================

calculate_concentration <- function(
    fuzzy_long,
    levels
) {
  
  full <- fuzzy_long %>%
    
    group_by(
      taxon_id
    ) %>%
    
    summarise(
      
      entropy_full =
        -sum(
          membership[
            membership > 0
          ] *
            log(
              membership[
                membership > 0
              ]
            )
        ),
      
      .groups = "drop"
    ) %>%
    
    mutate(
      
      concentration_full =
        1 -
        entropy_full /
        log(
          length(levels)
        )
    )
  
  
  active <- fuzzy_long %>%
    
    filter(
      membership > 0
    ) %>%
    
    group_by(
      taxon_id
    ) %>%
    
    summarise(
      
      n_active_levels =
        n(),
      
      entropy_active =
        -sum(
          membership *
            log(
              membership
            )
        ),
      
      .groups = "drop"
    ) %>%
    
    mutate(
      
      concentration_active =
        ifelse(
          
          n_active_levels <= 1,
          
          1,
          
          1 -
            entropy_active /
            log(
              n_active_levels
            )
        )
    )
  
  
  full %>%
    left_join(
      active,
      by = "taxon_id"
    )
}


# ============================================================
# EVIDENCE DEPTH
# ============================================================

calculate_evidence_depth <- function(
    chain
) {
  
  source_counts <- chain %>%
    distinct(
      taxon_id,
      source
    ) %>%
    count(
      taxon_id,
      name =
        "n_source_assignments"
    )
  
  
  chain_totals <- chain %>%
    group_by(
      taxon_id
    ) %>%
    summarise(
      
      effective_evidence =
        sum(
          weight,
          na.rm = TRUE
        ),
      
      n_genuine_changes =
        sum(
          relation ==
            "change",
          na.rm = TRUE
        ),
      
      n_refinements =
        sum(
          relation ==
            "refinement",
          na.rm = TRUE
        ),
      
      n_repeats =
        sum(
          relation ==
            "repeat",
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  
  source_counts %>%
    left_join(
      chain_totals,
      by = "taxon_id"
    )
}


assemble_classification <- function(
    fuzzy_long,
    chain,
    levels,
    system_name
) {
  
  support <-
    calculate_taxon_support(
      fuzzy_long
    )
  
  concentration <-
    calculate_concentration(
      fuzzy_long,
      levels
    )
  
  evidence <-
    calculate_evidence_depth(
      chain
    )
  
  
  support %>%
    
    left_join(
      concentration,
      by = "taxon_id"
    ) %>%
    
    left_join(
      evidence,
      by = "taxon_id"
    ) %>%
    
    mutate(
      
      # Analytical only.
      # Not exported as a user-facing confidence score.
      
      evidence_strength =
        1 -
        exp(
          -effective_evidence
        ),
      
      classification_confidence_index =
        support_best *
        concentration_active *
        evidence_strength,
      
      classification_system =
        system_name
    )
}


classification_7 <- assemble_classification(
  fuzzy_7_long,
  chain_7,
  WPFG_7_levels,
  "7-level"
)


classification_10 <- assemble_classification(
  fuzzy_10_long,
  chain_10,
  WPFG_10_levels,
  "10-level"
)


# ============================================================
# WIDE FUZZY MATRICES
# ============================================================

fuzzy_7_wide <- fuzzy_7_long %>%
  
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


fuzzy_10_wide <- fuzzy_10_long %>%
  
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
# MASTER 7 / MASTER 10
# ============================================================

master_7 <- taxon_master %>%
  
  left_join(
    classification_7,
    by = "taxon_id"
  ) %>%
  
  left_join(
    fuzzy_7_wide,
    by = "taxon_id"
  )


master_10 <- taxon_master %>%
  
  left_join(
    classification_10,
    by = "taxon_id"
  ) %>%
  
  left_join(
    fuzzy_10_wide,
    by = "taxon_id"
  )


# ============================================================
# CONFUSION
# ============================================================

confusion_events_7 <- chain_7 %>%
  
  filter(
    
    relation ==
      "change",
    
    previous_WPFG %in%
      WPFG_7_levels,
    
    WPFG %in%
      WPFG_7_levels
  ) %>%
  
  transmute(
    
    taxon_id,
    source,
    chrono,
    
    from_WPFG =
      previous_WPFG,
    
    to_WPFG =
      WPFG,
    
    weight
  )


confusion_events_10 <- chain_10 %>%
  
  filter(
    
    relation ==
      "change",
    
    previous_WPFG %in%
      WPFG_10_levels,
    
    WPFG %in%
      WPFG_10_levels
  ) %>%
  
  transmute(
    
    taxon_id,
    source,
    chrono,
    
    from_WPFG =
      previous_WPFG,
    
    to_WPFG =
      WPFG,
    
    weight
  )


# ============================================================
# DIRECTIONAL CONFUSION
# ============================================================

build_directional_confusion <- function(
    events,
    levels
) {
  
  events %>%
    
    count(
      
      from_WPFG,
      to_WPFG,
      
      wt =
        weight,
      
      name =
        "weighted_confusion"
    ) %>%
    
    complete(
      
      from_WPFG =
        levels,
      
      to_WPFG =
        levels,
      
      fill =
        list(
          weighted_confusion =
            0
        )
    ) %>%
    
    group_by(
      from_WPFG
    ) %>%
    
    mutate(
      
      row_total =
        sum(
          weighted_confusion
        ),
      
      row_probability =
        if_else(
          
          row_total == 0,
          
          0,
          
          weighted_confusion /
            row_total
        )
    ) %>%
    
    ungroup()
}


directional_confusion_7 <-
  build_directional_confusion(
    confusion_events_7,
    WPFG_7_levels
  )


directional_confusion_10 <-
  build_directional_confusion(
    confusion_events_10,
    WPFG_10_levels
  )


# ============================================================
# POOLED CONFUSION
# ============================================================

build_pooled_confusion <- function(
    events
) {
  
  if (nrow(events) == 0) {
    
    return(
      tibble(
        
        Level_Min =
          character(),
        
        Level_Max =
          character(),
        
        weighted_cross_classification =
          numeric(),
        
        n_events =
          integer(),
        
        n_taxa =
          integer(),
        
        proportion_of_all_conflicts =
          numeric()
      )
    )
  }
  
  
  pooled <- events %>%
    
    mutate(
      
      Level_Min =
        pmin(
          from_WPFG,
          to_WPFG
        ),
      
      Level_Max =
        pmax(
          from_WPFG,
          to_WPFG
        )
    ) %>%
    
    group_by(
      
      Level_Min,
      Level_Max
      
    ) %>%
    
    summarise(
      
      weighted_cross_classification =
        sum(
          weight,
          na.rm = TRUE
        ),
      
      n_events =
        n(),
      
      n_taxa =
        n_distinct(
          taxon_id
        ),
      
      .groups = "drop"
    )
  
  
  total_conflicts <-
    sum(
      pooled$weighted_cross_classification,
      na.rm = TRUE
    )
  
  
  pooled %>%
    
    mutate(
      
      proportion_of_all_conflicts =
        
        if (
          total_conflicts == 0
        ) {
          
          0
          
        } else {
          
          weighted_cross_classification /
            total_conflicts
          
        }
    ) %>%
    
    arrange(
      desc(
        weighted_cross_classification
      )
    )
}


pooled_confusion_7 <-
  build_pooled_confusion(
    confusion_events_7
  )


pooled_confusion_10 <-
  build_pooled_confusion(
    confusion_events_10
  )


# ============================================================
# SYMMETRIC MATRICES
# ============================================================

make_symmetric_matrix <- function(
    pooled,
    levels
) {
  
  result <- matrix(
    
    0,
    
    nrow =
      length(levels),
    
    ncol =
      length(levels),
    
    dimnames =
      list(
        levels,
        levels
      )
  )
  
  
  pooled <- pooled %>%
    filter(
      
      Level_Min %in%
        levels,
      
      Level_Max %in%
        levels
    )
  
  
  if (nrow(pooled) > 0) {
    
    for (
      i in seq_len(
        nrow(pooled)
      )
    ) {
      
      a <-
        pooled$Level_Min[i]
      
      b <-
        pooled$Level_Max[i]
      
      v <-
        pooled$proportion_of_all_conflicts[i]
      
      result[a, b] <-
        v
      
      result[b, a] <-
        v
    }
  }
  
  
  as.data.frame(
    result
  ) %>%
    
    rownames_to_column(
      "WPFG"
    )
}


conf_mat_7 <-
  make_symmetric_matrix(
    pooled_confusion_7,
    WPFG_7_levels
  )


conf_mat_10 <-
  make_symmetric_matrix(
    pooled_confusion_10,
    WPFG_10_levels
  )


# ============================================================
# REFINEMENT EVENTS
# ============================================================

refinement_events_10 <- chain_10 %>%
  
  filter(
    relation ==
      "refinement"
  ) %>%
  
  transmute(
    
    taxon_id,
    source,
    chrono,
    
    from_WPFG =
      previous_WPFG,
    
    to_WPFG =
      WPFG,
    
    weight
  )


# ============================================================
# QA
# ============================================================

evidence_summary <- resolved %>%
  
  group_by(
    source
  ) %>%
  
  summarise(
    
    n_total_records =
      n(),
    
    n_species_level_records =
      sum(
        WPFG_evidence_eligible,
        na.rm = TRUE
      ),
    
    n_higher_rank_records =
      sum(
        !WPFG_evidence_eligible,
        na.rm = TRUE
      ),
    
    proportion_species_level =
      mean(
        WPFG_evidence_eligible,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  left_join(
    source_metadata,
    by = "source"
  )


system_summary <- tibble(
  
  system =
    c(
      "7-level",
      "10-level"
    ),
  
  n_taxa =
    c(
      n_distinct(
        chain_7$taxon_id
      ),
      
      n_distinct(
        chain_10$taxon_id
      )
    ),
  
  n_source_assignments =
    c(
      nrow(chain_7),
      nrow(chain_10)
    ),
  
  n_genuine_conflicts =
    c(
      nrow(confusion_events_7),
      nrow(confusion_events_10)
    ),
  
  n_refinements =
    c(
      sum(
        chain_7$relation ==
          "refinement"
      ),
      
      sum(
        chain_10$relation ==
          "refinement"
      )
    )
)


# ============================================================
# VALIDATION GATES
# ============================================================

invalid_fuzzy_7 <- setdiff(
  
  unique(
    na.omit(
      fuzzy_7_long$terminal_WPFG
    )
  ),
  
  WPFG_7_levels
)


if (length(invalid_fuzzy_7) > 0) {
  
  stop(
    "7-level fuzzy output failed validation."
  )
}


invalid_fuzzy_10 <- setdiff(
  
  unique(
    na.omit(
      fuzzy_10_long$terminal_WPFG
    )
  ),
  
  WPFG_10_levels
)


if (length(invalid_fuzzy_10) > 0) {
  
  stop(
    "10-level fuzzy output failed validation."
  )
}


# ============================================================
# EXPORT
# ============================================================

write_xlsx(
  
  list(
    
    taxon_master =
      taxon_master,
    
    master_records =
      resolved,
    
    evidence_records =
      evidence_records,
    
    evidence_summary =
      evidence_summary,
    
    species_source_7 =
      species_source_assignments_7,
    
    species_source_10 =
      species_source_assignments_10,
    
    chain_7 =
      chain_7,
    
    chain_10 =
      chain_10,
    
    fuzzy_7_long =
      fuzzy_7_long,
    
    fuzzy_10_long =
      fuzzy_10_long,
    
    classification_7 =
      classification_7,
    
    classification_10 =
      classification_10,
    
    master_7 =
      master_7,
    
    master_10 =
      master_10,
    
    confusion_events_7 =
      confusion_events_7,
    
    confusion_events_10 =
      confusion_events_10,
    
    directional_confusion_7 =
      directional_confusion_7,
    
    directional_confusion_10 =
      directional_confusion_10,
    
    pooled_confusion_7 =
      pooled_confusion_7,
    
    pooled_confusion_10 =
      pooled_confusion_10,
    
    conf_mat_7 =
      conf_mat_7,
    
    conf_mat_10 =
      conf_mat_10,
    
    refinement_events_10 =
      refinement_events_10,
    
    system_summary =
      system_summary
    
  ),
  
  analysis_file
)


cat(
  "\nScript 2 complete.\n"
)