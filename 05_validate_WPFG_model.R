# ============================================================
# 05_validate_WPFG_model.R
# ============================================================
#
# PURPOSE
# -------
# Audit and validate the WPFG classification model.
#
# This script is NOT part of the routine database release
# pipeline. It is run separately when the model or source
# database is being evaluated.
#
# INPUT:
#   WPFG_classification_analysis_VerX.xlsx
#
# OUTPUT:
#   WPFG_validation_VerX.xlsx
#
# TESTS:
#   1. Structural validation
#   2. Synthetic weighting validation
#   3. Real historical assignment-chain diagnostics
#   4. S / hierarchical relationship validation
#   5. Metric distributions
#   6. Relationships among evidence metrics
#   7. Sensitivity to alpha_penalty
#   8. Identification of ambiguous / high-agreement taxa
#
# IMPORTANT:
# This script does not modify the classification model.
# It reports diagnostics that should be interpreted alongside
# biological knowledge and manual inspection of selected taxa.
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(tidyverse)
library(readxl)
library(writexl)


# ============================================================
# 2. RELEASE SETTINGS
# ============================================================
#
# Run this script after run_all.R, so the current release
# objects are available in the R session.
#
# Alternatively, set db_version manually here.
#
# ============================================================

if (!exists("db_version")) {
  db_version <- 1
}

if (!exists("version_label")) {
  version_label <- paste0(
    "Ver",
    db_version
  )
}

if (!exists("analysis_file")) {
  analysis_file <- file.path(
    "data_out",
    paste0(
      "WPFG_classification_analysis_",
      version_label,
      ".xlsx"
    )
  )
}

validation_file <- file.path(
  "data_out",
  paste0(
    "WPFG_validation_",
    version_label,
    ".xlsx"
  )
)


# ============================================================
# 3. MODEL SETTINGS
# ============================================================

alpha_penalty <- 0.2

alpha_values <- c(
  0,
  0.1,
  0.2,
  0.3,
  0.5,
  1
)


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
# 4. READ ANALYSIS OUTPUT
# ============================================================

if (!file.exists(analysis_file)) {
  
  stop(
    "Analysis workbook not found:\n",
    analysis_file
  )
}


chain_7 <- read_excel(
  analysis_file,
  sheet = "chain_7"
)

chain_10 <- read_excel(
  analysis_file,
  sheet = "chain_10"
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

species_source_7 <- read_excel(
  analysis_file,
  sheet = "species_source_7"
)

species_source_10 <- read_excel(
  analysis_file,
  sheet = "species_source_10"
)

taxon_master <- read_excel(
  analysis_file,
  sheet = "taxon_master"
)


# ============================================================
# 5. STANDARDISE IDS
# ============================================================

chain_7 <- chain_7 %>%
  mutate(
    taxon_id =
      as.character(taxon_id)
  )

chain_10 <- chain_10 %>%
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

species_source_7 <- species_source_7 %>%
  mutate(
    taxon_id =
      as.character(taxon_id)
  )

species_source_10 <- species_source_10 %>%
  mutate(
    taxon_id =
      as.character(taxon_id)
  )

taxon_master <- taxon_master %>%
  mutate(
    taxon_id =
      as.character(taxon_id)
  )


# ============================================================
# 6. README / TEST GUIDE
# ============================================================

readme <- tibble(
  
  section = c(
    
    "Purpose",
    "How to use this workbook",
    "Structural validation",
    "Synthetic weighting tests",
    "Historical chain diagnostics",
    "Seven- vs ten-level comparison",
    "S / hierarchical validation",
    "Metric summaries",
    "Metric correlations",
    "Alpha sensitivity",
    "Ambiguous taxa",
    "High-agreement taxa",
    "Single-source taxa",
    "Many-source taxa",
    "Interpretation of support",
    "Interpretation of effective evidence",
    "Interpretation of margin",
    "Interpretation of the confidence index",
    "What constitutes good validation",
    "What this workbook does not establish"
  ),
  
  information = c(
    
    paste0(
      "Validation diagnostics for WPFG classification model ",
      version_label,
      ". The tests assess whether the computational implementation behaves as intended and whether the selected classifications are robust to key modelling assumptions."
    ),
    
    "Start with structural_validation and weighting_tests. Then inspect alpha sensitivity, metric summaries and selected real taxa. Ambiguous and alpha-sensitive taxa should be manually reviewed for biological plausibility.",
    
    "These tests check that only valid WPFG categories occur, fuzzy memberships sum to one, and each taxon contributes at most one source-level assignment per source. All reported problems should be resolved before interpreting model results.",
    
    "Synthetic chains test the mechanics of the dependence weighting. For alpha = 0.2, an initial assignment receives weight 1.0, an unchanged repeat receives 0.2, and a genuine change or hierarchical refinement receives 1.0.",
    
    "Real chains show the chronological evidence history for individual taxa, including the relation between successive classifications and the weights applied. These are useful for spot-checking whether the model's interpretation agrees with the actual source history.",
    
    "The two systems are alternative resolutions. Similar results are not inherently expected or required. The purpose is to confirm that each system is operating on its intended evidence rather than to establish that one resolution is superior.",
    
    "In the ten-level system, S is an unspecified parent category covering Sr, Sk and Se. S-to-subclass and subclass-to-S transitions should be treated as refinements/coarsenings, not as genuine conflicts. S evidence is distributed equally across Sr, Sk and Se in the fuzzy classification.",
    
    "Metric summaries describe the empirical distribution of support, margin, effective evidence and concentration across all classified taxa. They are descriptive rather than predefined acceptance criteria.",
    
    "Correlations show whether the composite confidence index is providing information distinct from support, margin, concentration and effective evidence. A strong correlation does not indicate a computational error; it indicates that metrics may be measuring related properties.",
    
    "The model is recalculated across alpha values from 0 to 1. Alpha = 0 treats unchanged repeats as providing no additional evidence; alpha = 1 treats them as fully independent. Stability of the winning WPFG across this range indicates robustness to the dependence assumption.",
    
    "Ambiguous taxa have weak separation between the leading and second-ranked WPFGs. They are the most useful cases for manual biological review and for judging whether fuzzy membership provides an adequate representation of uncertainty.",
    
    "High-agreement taxa have strong support for one WPFG. These are useful for confirming that the model behaves sensibly in uncomplicated cases.",
    
    "Single-source taxa illustrate an important distinction: support can equal 1.0 even when evidence depth is minimal. Such cases should not be interpreted as equivalent to strong multi-source evidence.",
    
    "Many-source taxa show how the model behaves where repeated classification provides substantially more evidence and where genuine disagreement is more likely to emerge.",
    
    "Support is the proportion of dependence-adjusted evidence supporting the leading WPFG. High support indicates agreement among the available evidence, but does not by itself indicate how much evidence exists.",
    
    "Effective evidence is the amount of dependence-adjusted historical information contributing to a classification. It distinguishes, for example, a single-source classification from one supported across several source lists.",
    
    "Margin is the difference between the best and second-best WPFG memberships. A small margin means that the best guess is only weakly separated from its principal alternative.",
    
    "The composite confidence index was examined during model development but is not treated as the primary user-facing measure. In this dataset it is strongly correlated with support, margin and concentration, while effective evidence shows a negative relationship because taxa with more historical coverage also have more opportunities for genuine disagreement.",
    
    "A well-behaved model should pass the structural checks, reproduce the expected synthetic weighting patterns, correctly identify S-related refinements, and produce classifications that are robust to plausible alpha values. Real ambiguous cases should also be biologically interpretable when inspected manually.",
    
    "These tests validate computational behaviour and robustness of the classification framework. They do not establish that any individual WPFG assignment is biologically correct, nor do they calibrate support as a probability of correctness."
  )
)


readme_metadata <- tibble(
  
  field = c(
    "Database version",
    "Analysis workbook",
    "Alpha used for primary analysis",
    "Alpha values used for sensitivity analysis"
  ),
  
  value = c(
    version_label,
    basename(analysis_file),
    alpha_penalty,
    paste(
      alpha_values,
      collapse = ", "
    )
  )
)


# ============================================================
# 7. STRUCTURAL VALIDATION
# ============================================================

invalid_fuzzy_7 <- setdiff(
  unique(
    na.omit(
      fuzzy_7_long$terminal_WPFG
    )
  ),
  WPFG_7_levels
)


invalid_fuzzy_10 <- setdiff(
  unique(
    na.omit(
      fuzzy_10_long$terminal_WPFG
    )
  ),
  WPFG_10_levels
)


membership_sum_7 <- fuzzy_7_long %>%
  group_by(
    taxon_id
  ) %>%
  summarise(
    membership_sum =
      sum(
        membership,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  mutate(
    valid =
      abs(
        membership_sum - 1
      ) < 1e-10
  )


membership_sum_10 <- fuzzy_10_long %>%
  group_by(
    taxon_id
  ) %>%
  summarise(
    membership_sum =
      sum(
        membership,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  mutate(
    valid =
      abs(
        membership_sum - 1
      ) < 1e-10
  )


duplicate_7 <- species_source_7 %>%
  count(
    taxon_id,
    source
  ) %>%
  filter(
    n > 1
  )


duplicate_10 <- species_source_10 %>%
  count(
    taxon_id,
    source
  ) %>%
  filter(
    n > 1
  )


structural_validation <- tibble(
  
  check = c(
    
    "Invalid 7-level fuzzy categories",
    
    "Invalid 10-level fuzzy categories",
    
    "7-level fuzzy memberships sum to 1",
    
    "10-level fuzzy memberships sum to 1",
    
    "Duplicate 7-level taxon x source assignments",
    
    "Duplicate 10-level taxon x source assignments"
  ),
  
  n_problem = c(
    
    length(invalid_fuzzy_7),
    
    length(invalid_fuzzy_10),
    
    sum(
      !membership_sum_7$valid
    ),
    
    sum(
      !membership_sum_10$valid
    ),
    
    nrow(duplicate_7),
    
    nrow(duplicate_10)
  ),
  
  status = if_else(
    
    c(
      
      length(invalid_fuzzy_7),
      
      length(invalid_fuzzy_10),
      
      sum(
        !membership_sum_7$valid
      ),
      
      sum(
        !membership_sum_10$valid
      ),
      
      nrow(duplicate_7),
      
      nrow(duplicate_10)
    ) == 0,
    
    "PASS",
    
    "CHECK"
  )
)


# ============================================================
# 8. SYNTHETIC WEIGHTING TESTS
# ============================================================

simulate_chain <- function(
    assignments,
    alpha = 0.2
) {
  
  assignments <-
    as.character(
      assignments
    )
  
  previous <-
    lag(assignments)
  
  relation <- case_when(
    
    is.na(previous) ~
      "initial",
    
    assignments == previous ~
      "repeat",
    
    TRUE ~
      "change"
  )
  
  weight <- case_when(
    
    relation == "initial" ~
      1,
    
    relation == "repeat" ~
      alpha,
    
    relation == "change" ~
      1
  )
  
  tibble(
    assignment =
      assignments,
    
    previous =
      previous,
    
    relation =
      relation,
    
    weight =
      weight
  )
}


test_chains <- list(
  
  single =
    c("A"),
  
  repeat_two =
    c("A", "A"),
  
  repeat_three =
    c("A", "A", "A"),
  
  simple_change =
    c("A", "B"),
  
  repeat_then_change =
    c("A", "A", "B"),
  
  change_then_repeat =
    c("A", "B", "B"),
  
  alternating =
    c("A", "B", "A")
)


weighting_tests <- map_dfr(
  
  names(test_chains),
  
  function(name) {
    
    result <-
      simulate_chain(
        test_chains[[name]],
        alpha = alpha_penalty
      )
    
    tibble(
      
      test =
        name,
      
      chain =
        paste(
          test_chains[[name]],
          collapse = " -> "
        ),
      
      weights =
        paste(
          result$weight,
          collapse = " + "
        ),
      
      effective_evidence =
        sum(
          result$weight
        )
    )
  }
)


expected_weighting <- tibble(
  
  test = c(
    "single",
    "repeat_two",
    "repeat_three",
    "simple_change",
    "repeat_then_change",
    "change_then_repeat",
    "alternating"
  ),
  
  expected_effective_evidence = c(
    1,
    1 + alpha_penalty,
    1 + 2 * alpha_penalty,
    2,
    2 + alpha_penalty,
    2 + alpha_penalty,
    3
  )
)


weighting_tests <- weighting_tests %>%
  left_join(
    expected_weighting,
    by = "test"
  ) %>%
  mutate(
    
    matches_expected =
      abs(
        effective_evidence -
          expected_effective_evidence
      ) < 1e-10
  )


# ============================================================
# 9. REAL HISTORICAL CHAIN DIAGNOSTICS
# ============================================================

real_chain_summary_7 <- chain_7 %>%
  
  group_by(
    taxon_id
  ) %>%
  
  summarise(
    
    chain =
      paste(
        WPFG,
        collapse = " -> "
      ),
    
    relations =
      paste(
        relation,
        collapse = " -> "
      ),
    
    weights =
      paste(
        weight,
        collapse = " + "
      ),
    
    n_sources =
      n(),
    
    effective_evidence =
      sum(
        weight,
        na.rm = TRUE
      ),
    
    n_changes =
      sum(
        relation == "change"
      ),
    
    n_refinements =
      sum(
        relation == "refinement"
      ),
    
    .groups = "drop"
  )


real_chain_summary_10 <- chain_10 %>%
  
  group_by(
    taxon_id
  ) %>%
  
  summarise(
    
    chain =
      paste(
        WPFG,
        collapse = " -> "
      ),
    
    relations =
      paste(
        relation,
        collapse = " -> "
      ),
    
    weights =
      paste(
        weight,
        collapse = " + "
      ),
    
    n_sources =
      n(),
    
    effective_evidence =
      sum(
        weight,
        na.rm = TRUE
      ),
    
    n_changes =
      sum(
        relation == "change"
      ),
    
    n_refinements =
      sum(
        relation == "refinement"
      ),
    
    .groups = "drop"
  )


# ============================================================
# 10. S / HIERARCHICAL VALIDATION
# ============================================================

s_transition_check <- chain_10 %>%
  
  filter(
    previous_WPFG == "S" |
      WPFG == "S"
  ) %>%
  
  select(
    taxon_id,
    source,
    chrono,
    previous_WPFG,
    WPFG,
    relation,
    weight
  )


# Expected S relationships:
#
# S <-> Sr
# S <-> Sk
# S <-> Se
#
# should be refinements.
#
# All other S-related transitions should be inspected.


s_expected <- s_transition_check %>%
  
  filter(
    (
      previous_WPFG == "S" &
        WPFG %in% c(
          "Sr",
          "Sk",
          "Se"
        )
    ) |
      (
        WPFG == "S" &
          previous_WPFG %in% c(
            "Sr",
            "Sk",
            "Se"
          )
      )
  )


s_expected_check <- s_expected %>%
  summarise(
    n_expected_transitions =
      n(),
    
    n_correctly_classified_as_refinement =
      sum(
        relation == "refinement"
      ),
    
    n_incorrect =
      sum(
        relation != "refinement"
      )
  )


s_classification_10 <-
  classification_10 %>%
  
  semi_join(
    s_transition_check,
    by = "taxon_id"
  ) %>%
  
  arrange(
    desc(effective_evidence)
  )


# ============================================================
# 11. METRIC SUMMARIES
# ============================================================

metric_summary_7 <- classification_7 %>%
  summarise(
    
    n_taxa =
      n(),
    
    support_min =
      min(
        support_best,
        na.rm = TRUE
      ),
    
    support_median =
      median(
        support_best,
        na.rm = TRUE
      ),
    
    support_mean =
      mean(
        support_best,
        na.rm = TRUE
      ),
    
    support_max =
      max(
        support_best,
        na.rm = TRUE
      ),
    
    margin_median =
      median(
        margin_best_vs_second,
        na.rm = TRUE
      ),
    
    effective_evidence_median =
      median(
        effective_evidence,
        na.rm = TRUE
      ),
    
    concentration_median =
      median(
        concentration_active,
        na.rm = TRUE
      ),
    
    confidence_index_median =
      median(
        classification_confidence_index,
        na.rm = TRUE
      )
  )


metric_summary_10 <- classification_10 %>%
  summarise(
    
    n_taxa =
      n(),
    
    support_min =
      min(
        support_best,
        na.rm = TRUE
      ),
    
    support_median =
      median(
        support_best,
        na.rm = TRUE
      ),
    
    support_mean =
      mean(
        support_best,
        na.rm = TRUE
      ),
    
    support_max =
      max(
        support_best,
        na.rm = TRUE
      ),
    
    margin_median =
      median(
        margin_best_vs_second,
        na.rm = TRUE
      ),
    
    effective_evidence_median =
      median(
        effective_evidence,
        na.rm = TRUE
      ),
    
    concentration_median =
      median(
        concentration_active,
        na.rm = TRUE
      ),
    
    confidence_index_median =
      median(
        classification_confidence_index,
        na.rm = TRUE
      )
  )


# ============================================================
# 12. METRIC CORRELATIONS
# ============================================================

metric_correlations_7 <- classification_7 %>%
  summarise(
    
    support_vs_confidence =
      cor(
        support_best,
        classification_confidence_index,
        use = "complete.obs"
      ),
    
    evidence_vs_confidence =
      cor(
        effective_evidence,
        classification_confidence_index,
        use = "complete.obs"
      ),
    
    margin_vs_confidence =
      cor(
        margin_best_vs_second,
        classification_confidence_index,
        use = "complete.obs"
      ),
    
    concentration_vs_confidence =
      cor(
        concentration_active,
        classification_confidence_index,
        use = "complete.obs"
      )
  )


metric_correlations_10 <- classification_10 %>%
  summarise(
    
    support_vs_confidence =
      cor(
        support_best,
        classification_confidence_index,
        use = "complete.obs"
      ),
    
    evidence_vs_confidence =
      cor(
        effective_evidence,
        classification_confidence_index,
        use = "complete.obs"
      ),
    
    margin_vs_confidence =
      cor(
        margin_best_vs_second,
        classification_confidence_index,
        use = "complete.obs"
      ),
    
    concentration_vs_confidence =
      cor(
        concentration_active,
        classification_confidence_index,
        use = "complete.obs"
      )
  )


# ============================================================
# 13. AMBIGUOUS TAXA
# ============================================================
#
# These are diagnostic samples, not formal "bad" cases.
# They are useful for manual review.
#
# ============================================================

ambiguous_7 <- classification_7 %>%
  filter(
    !is.na(
      margin_best_vs_second
    )
  ) %>%
  arrange(
    margin_best_vs_second
  ) %>%
  slice_head(
    n = 50
  )


ambiguous_10 <- classification_10 %>%
  filter(
    !is.na(
      margin_best_vs_second
    )
  ) %>%
  arrange(
    margin_best_vs_second
  ) %>%
  slice_head(
    n = 50
  )


# ============================================================
# 14. HIGH AGREEMENT TAXA
# ============================================================

high_agreement_7 <- classification_7 %>%
  filter(
    support_best >= 0.95
  ) %>%
  arrange(
    desc(
      effective_evidence
    )
  )


high_agreement_10 <- classification_10 %>%
  filter(
    support_best >= 0.95
  ) %>%
  arrange(
    desc(
      effective_evidence
    )
  )


# ============================================================
# 15. SINGLE-SOURCE TAXA
# ============================================================

single_source_7 <- classification_7 %>%
  filter(
    n_source_assignments == 1
  )


single_source_10 <- classification_10 %>%
  filter(
    n_source_assignments == 1
  )


# ============================================================
# 16. MANY-SOURCE TAXA
# ============================================================

many_source_7 <- classification_7 %>%
  filter(
    n_source_assignments >= 4
  ) %>%
  arrange(
    desc(
      effective_evidence
    )
  )


many_source_10 <- classification_10 %>%
  filter(
    n_source_assignments >= 4
  ) %>%
  arrange(
    desc(
      effective_evidence
    )
  )


# ============================================================
# 17. ALPHA SENSITIVITY
# ============================================================

calculate_alpha_result <- function(
    chain,
    levels,
    alpha,
    system
) {
  
  chain_alpha <- chain %>%
    mutate(
      
      alpha_weight =
        case_when(
          
          relation == "initial" ~
            1.0,
          
          relation == "repeat" ~
            alpha,
          
          relation == "change" ~
            1.0,
          
          relation == "refinement" ~
            1.0,
          
          TRUE ~
            NA_real_
        )
    )
  
  
  if (system == 7) {
    
    fuzzy_alpha <- chain_alpha %>%
      
      transmute(
        taxon_id,
        source,
        chrono,
        WPFG,
        alpha_weight,
        terminal_WPFG =
          WPFG,
        terminal_share =
          1
      )
    
  } else {
    
    fuzzy_alpha <- chain_alpha %>%
      
      mutate(
        
        terminal_WPFG =
          map(
            WPFG,
            ~ if (
              .x == "S"
            ) {
              
              c(
                "Sr",
                "Sk",
                "Se"
              )
              
            } else {
              
              .x
              
            }
          )
      ) %>%
      
      unnest_longer(
        terminal_WPFG
      ) %>%
      
      group_by(
        
        taxon_id,
        source,
        chrono,
        WPFG,
        alpha_weight
        
      ) %>%
      
      mutate(
        
        terminal_share =
          1 / n()
        
      ) %>%
      
      ungroup()
  }
  
  
  membership <- fuzzy_alpha %>%
    
    mutate(
      
      terminal_evidence =
        alpha_weight *
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
        evidence_weight = 0
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
        ifelse(
          effective_evidence > 0,
          evidence_weight /
            effective_evidence,
          0
        )
    ) %>%
    
    arrange(
      taxon_id,
      desc(membership),
      terminal_WPFG
    ) %>%
    
    mutate(
      rank =
        row_number()
    ) %>%
    
    ungroup()
  
  
  membership %>%
    
    filter(
      rank == 1
    ) %>%
    
    transmute(
      
      taxon_id,
      
      best_WPFG =
        terminal_WPFG,
      
      support =
        membership,
      
      effective_evidence =
        effective_evidence,
      
      alpha =
        alpha
    )
}


alpha_results_7 <- map_dfr(
  
  alpha_values,
  
  ~ calculate_alpha_result(
    chain =
      chain_7,
    
    levels =
      WPFG_7_levels,
    
    alpha =
      .x,
    
    system =
      7
  )
)


alpha_results_10 <- map_dfr(
  
  alpha_values,
  
  ~ calculate_alpha_result(
    chain =
      chain_10,
    
    levels =
      WPFG_10_levels,
    
    alpha =
      .x,
    
    system =
      10
  )
)


# ============================================================
# 18. ALPHA STABILITY
# ============================================================

alpha_stability_7 <- alpha_results_7 %>%
  
  group_by(
    taxon_id
  ) %>%
  
  summarise(
    
    n_distinct_best_guesses =
      n_distinct(
        best_WPFG
      ),
    
    proportion_same_as_alpha_02_across_alpha =
      mean(
        best_WPFG ==
          best_WPFG[
            alpha == alpha_penalty
          ][1]
      ),
    
    min_support =
      min(
        support
      ),
    
    max_support =
      max(
        support
      ),
    
    .groups = "drop"
  )


alpha_stability_10 <- alpha_results_10 %>%
  
  group_by(
    taxon_id
  ) %>%
  
  summarise(
    
    n_distinct_best_guesses =
      n_distinct(
        best_WPFG
      ),
    
    proportion_same_as_alpha_02_across_alpha =
      mean(
        best_WPFG ==
          best_WPFG[
            alpha == alpha_penalty
          ][1]
      ),
    
    min_support =
      min(
        support
      ),
    
    max_support =
      max(
        support
      ),
    
    .groups = "drop"
  )


overall_alpha_stability <- tibble(
  
  system =
    c(
      "7-level",
      "10-level"
    ),
  
  proportion_classifications_stable_all_alpha =
    c(
      
      mean(
        alpha_stability_7$
          n_distinct_best_guesses ==
          1
      ),
      
      mean(
        alpha_stability_10$
          n_distinct_best_guesses ==
          1
      )
    ),
  
  proportion_same_as_alpha_02_across_alpha =
    c(
      
      mean(
        alpha_stability_7$
          proportion_same_as_alpha_02_across_alpha ==
          1
      ),
      
      mean(
        alpha_stability_10$
          proportion_same_as_alpha_02_across_alpha ==
          1
      )
    )
)


# ============================================================
# 19. ALPHA-SENSITIVE TAXA
# ============================================================

alpha_sensitive_7 <- alpha_stability_7 %>%
  
  filter(
    n_distinct_best_guesses > 1
  ) %>%
  
  arrange(
    desc(
      n_distinct_best_guesses
    ),
    min_support
  )


alpha_sensitive_10 <- alpha_stability_10 %>%
  
  filter(
    n_distinct_best_guesses > 1
  ) %>%
  
  arrange(
    desc(
      n_distinct_best_guesses
    ),
    min_support
  )


# ============================================================
# 20. VALIDATION SUMMARY
# ============================================================

validation_summary <- tibble(
  
  test = c(
    
    "Structural validation",
    
    "Synthetic weighting",
    
    "S / hierarchical transitions",
    
    "Alpha sensitivity, 7-level",
    
    "Alpha sensitivity, 10-level"
  ),
  
  result = c(
    
    ifelse(
      all(
        structural_validation$status ==
          "PASS"
      ),
      "PASS",
      "CHECK"
    ),
    
    ifelse(
      all(
        weighting_tests$matches_expected
      ),
      "PASS",
      "CHECK"
    ),
    
    ifelse(
      nrow(s_expected_check) == 0 ||
        s_expected_check$n_incorrect == 0,
      "PASS",
      "CHECK"
    ),
    
    paste0(
      round(
        overall_alpha_stability$
          proportion_classifications_stable_all_alpha[1] *
          100,
        1
      ),
      "% stable across tested alpha values"
    ),
    
    paste0(
      round(
        overall_alpha_stability$
          proportion_classifications_stable_all_alpha[2] *
          100,
        1
      ),
      "% stable across tested alpha values"
    )
  )
)


# ============================================================
# 21. EXPORT VALIDATION WORKBOOK
# ============================================================

write_xlsx(
  
  list(
    
    README =
      bind_rows(
        
        readme_metadata %>%
          transmute(
            section =
              "Release metadata",
            
            column =
              field,
            
            description =
              value
          ),
        
        tibble(
          section = "",
          column = "",
          description = ""
        ),
        
        readme %>%
          transmute(
            section,
            column = "",
            description = information
          ),
        
        tibble(
          section = "",
          column = "",
          description = ""
        ),
        
        validation_summary %>%
          transmute(
            section =
              "Validation summary",
            
            column =
              test,
            
            description =
              result
          )
        
      ),
    
    structural_validation =
      structural_validation,
    
    weighting_tests =
      weighting_tests,
    
    weighting_expected =
      expected_weighting,
    
    real_chain_summary_7 =
      real_chain_summary_7,
    
    real_chain_summary_10 =
      real_chain_summary_10,
    
    s_transition_check =
      s_transition_check,
    
    s_expected_check =
      s_expected_check,
    
    s_classification_10 =
      s_classification_10,
    
    metric_summary_7 =
      metric_summary_7,
    
    metric_summary_10 =
      metric_summary_10,
    
    metric_correlations_7 =
      metric_correlations_7,
    
    metric_correlations_10 =
      metric_correlations_10,
    
    high_agreement_7 =
      high_agreement_7,
    
    high_agreement_10 =
      high_agreement_10,
    
    ambiguous_7 =
      ambiguous_7,
    
    ambiguous_10 =
      ambiguous_10,
    
    single_source_7 =
      single_source_7,
    
    single_source_10 =
      single_source_10,
    
    many_source_7 =
      many_source_7,
    
    many_source_10 =
      many_source_10,
    
    alpha_results_7 =
      alpha_results_7,
    
    alpha_results_10 =
      alpha_results_10,
    
    alpha_stability_7 =
      alpha_stability_7,
    
    alpha_stability_10 =
      alpha_stability_10,
    
    overall_alpha_stability =
      overall_alpha_stability,
    
    alpha_sensitive_7 =
      alpha_sensitive_7,
    
    alpha_sensitive_10 =
      alpha_sensitive_10
    
  ),
  
  validation_file
)


# ============================================================
# 22. CONSOLE SUMMARY
# ============================================================

cat(
  "\n============================================\n",
  "WPFG VALIDATION COMPLETE\n",
  "============================================\n\n"
)


cat(
  "Version: ",
  version_label,
  "\n\n"
)


print(
  validation_summary
)


cat(
  "\nValidation workbook:\n",
  validation_file,
  "\n"
)