# ============================================================
# run_all.R
# ============================================================
#
# Rebuild the complete WPFG database.
#
# Change ONLY:
#   db_version
#   apc_data_type
#   apc_version
#
# and provide the matching versioned input CSV.
#
# ============================================================
renv::snapshot()
# renv::restore() # to restore environment
# ============================================================
# 1. RELEASE SETTINGS
# ============================================================

db_version <- 1

apc_data_type <- "stable"

apc_version <- "2026-03-25"


# ============================================================
# 2. VERSION LABEL
# ============================================================

version_label <- paste0(
  "Ver",
  db_version
)


# ============================================================
# 3. DIRECTORIES
# ============================================================

data_in <- "data_in"

data_out <- "data_out"

if (!dir.exists(data_out)) {
  dir.create(
    data_out,
    recursive = TRUE
  )
}


# ============================================================
# 4. VERSIONED FILE PATHS
# ============================================================

input_file <- file.path(
  data_in,
  paste0(
    "wpfg_lists_",
    version_label,
    ".csv"
  )
)

master_file <- file.path(
  data_out,
  paste0(
    "plant_functional_group_master_",
    version_label,
    ".xlsx"
  )
)

analysis_file <- file.path(
  data_out,
  paste0(
    "WPFG_classification_analysis_",
    version_label,
    ".xlsx"
  )
)

genus_file <- file.path(
  data_out,
  paste0(
    "WPFG_genus_analysis_",
    version_label,
    ".xlsx"
  )
)

user_file <- file.path(
  data_out,
  paste0(
    "WPFG_user_database_",
    version_label,
    ".xlsx"
  )
)


# ============================================================
# 5. CHECK INPUT
# ============================================================

if (!file.exists(input_file)) {
  
  stop(
    "\nInput file not found:\n",
    input_file,
    "\n\nExpected filename:\n",
    paste0(
      "wpfg_lists_",
      version_label,
      ".csv"
    ),
    "\n"
  )
}


# ============================================================
# 6. MAKE RELEASE SETTINGS AVAILABLE
# ============================================================

assign(
  "db_version",
  db_version,
  envir = .GlobalEnv
)

assign(
  "version_label",
  version_label,
  envir = .GlobalEnv
)

assign(
  "apc_data_type",
  apc_data_type,
  envir = .GlobalEnv
)

assign(
  "apc_version",
  apc_version,
  envir = .GlobalEnv
)

assign(
  "input_file",
  input_file,
  envir = .GlobalEnv
)

assign(
  "master_file",
  master_file,
  envir = .GlobalEnv
)

assign(
  "analysis_file",
  analysis_file,
  envir = .GlobalEnv
)

assign(
  "genus_file",
  genus_file,
  envir = .GlobalEnv
)

assign(
  "user_file",
  user_file,
  envir = .GlobalEnv
)


# ============================================================
# 7. START
# ============================================================

cat(
  "\n============================================\n",
  "WPFG DATABASE BUILD ",
  version_label,
  "\n",
  "============================================\n\n"
)

cat(
  "Source file:\n",
  input_file,
  "\n\n"
)

cat(
  "APC data type: ",
  apc_data_type,
  "\n"
)

cat(
  "APC version: ",
  apc_version,
  "\n\n"
)


# ============================================================
# 8. RUN PIPELINE
# ============================================================

cat(
  "1/4  Building master database...\n"
)

source(
  "01_build_master_database.R"
)


cat(
  "\n2/4  Running WPFG analysis...\n"
)

source(
  "02_wpfg_analysis.R"
)


cat(
  "\n3/4  Running genus analysis...\n"
)

source(
  "03_genus_analysis.R"
)


cat(
  "\n4/4  Building user database...\n"
)

source(
  "04_build_user_database.R"
)


# ============================================================
# 9. COMPLETE
# ============================================================

cat(
  "\n============================================\n",
  "BUILD COMPLETE - ",
  version_label,
  "\n",
  "============================================\n\n"
)

cat(
  "Final database:\n",
  user_file,
  "\n\n"
)

cat(
  "Supporting files:\n",
  master_file,
  "\n",
  analysis_file,
  "\n",
  genus_file,
  "\n"
)