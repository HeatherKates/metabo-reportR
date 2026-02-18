#' SECIM Metabolomics Pipeline - Master Execution Script
#' 
#' This script orchestrates the two-stage metabolomics analysis workflow:
#' Stage 1: Data processing, normalization, and statistical analysis
#' Stage 2: Interactive report generation via Quarto
#' 
#' @usage Rscript Run_Pipeline.R [path/to/config.yaml]
#' @author SECIM - University of Florida
#' @date 2026

# ===========================
# Load Required Libraries
# ===========================
suppressPackageStartupMessages({
  library(yaml)        # Configuration parsing
  library(qs)          # Fast serialization
  library(dplyr)       # Data manipulation
  library(quarto)      # Report rendering
})

# ===========================
# LOAD CONFIGURATION
# ===========================

# Check for config file argument
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  config_file <- args[1]
} else {
  # Default config file location
  config_file <- "config/study_params.yaml"
}

# Verify config file exists
if (!file.exists(config_file)) {
  stop(sprintf("ERROR: Configuration file not found: %s\n\nUsage: Rscript Run_Pipeline.R [path/to/config.yaml]\n", config_file))
}

cat("\n========================================\n")
cat("SECIM Metabolomics Pipeline\n")
cat("========================================\n\n")

cat(sprintf("Loading configuration: %s\n", config_file))

# Parse YAML configuration
config <- yaml::read_yaml(config_file)

# Extract nested parameters into flat structure
client <- config$project$client
PI_name <- config$project$PI_name

Input <- config$input$file
mzmine_version <- config$input$mzmine_version

samples_to_drop_pre_norm <- config$filtering$samples_to_drop_pre_norm
samples_to_drop_post_norm <- config$filtering$samples_to_drop_post_norm

test_type <- config$statistics$test_type
contrast_var <- config$statistics$contrast_var
ReferenceLevel <- config$statistics$reference_level
paired <- config$statistics$paired
num_meta <- config$statistics$num_meta
subset <- config$statistics$subset

# Convert formula strings to R formulas
anova_formula <- if (!is.null(config$statistics$anova_formula)) {
  as.formula(config$statistics$anova_formula)
} else {
  NULL
}

lm_model <- if (!is.null(config$statistics$lm_model)) {
  as.formula(config$statistics$lm_model)
} else {
  NULL
}

SECIM_column <- config$annotation$SECIM_column
metid_DB_file <- config$annotation$metid_DB_file

batch_correct <- config$processing$batch_correct
rowNorm <- config$processing$rowNorm
transNorm <- config$processing$transNorm
scaleNorm <- config$processing$scaleNorm

cat("✓ Configuration loaded\n\n")

# Source the main processing function
source("R/Generate_Reporting_Inputs.R")

# ===========================
# VALIDATION CHECKS
# ===========================

cat("Validating configuration...\n")

# Check input file exists
if (!file.exists(Input)) {
  stop(sprintf("ERROR: Input file not found: %s\n", Input))
}

# Check MZmine version
if (!mzmine_version %in% c(2, 3)) {
  stop("ERROR: mzmine_version must be 2 or 3\n")
}

# Check test type
valid_tests <- c("t.test", "anova", "lm", "lme", "nostats")
if (!test_type %in% valid_tests) {
  stop(sprintf("ERROR: test_type must be one of: %s\n", paste(valid_tests, collapse = ", ")))
}

# Check SECIM column
valid_columns <- c("ACE-PFP", "Evospere-PFP")
if (!SECIM_column %in% valid_columns) {
  stop(sprintf("ERROR: SECIM_column must be one of: %s\n", paste(valid_columns, collapse = ", ")))
}

# Check metabolite database
metid_DB_path <- file.path("reference", metid_DB_file)
if (!file.exists(metid_DB_path)) {
  stop(sprintf("ERROR: Metabolite database not found: %s\n", metid_DB_path))
}

# Validate formulas for advanced test types
if (test_type == "anova" && is.null(anova_formula)) {
  stop("ERROR: anova_formula is required when test_type = 'anova'\n")
}

if (test_type %in% c("lm", "lme") && is.null(lm_model)) {
  stop("ERROR: lm_model is required when test_type = 'lm' or 'lme'\n")
}

cat("✓ Configuration validated\n\n")

# ===========================
# STAGE 1: DATA PROCESSING
# ===========================

cat("========================================\n")
cat("Stage 1: Processing and Statistical Analysis\n")
cat("========================================\n\n")

cat(sprintf("Client: %s\n", client))
cat(sprintf("PI: %s\n", PI_name))
cat(sprintf("Input: %s\n", Input))
cat(sprintf("MZmine Version: %d\n", mzmine_version))
cat(sprintf("Test Type: %s\n", test_type))
cat(sprintf("Contrast Variable: %s\n", contrast_var))
cat(sprintf("SECIM Column: %s\n", SECIM_column))
cat(sprintf("Metabolite Database: %s\n", metid_DB_file))
cat("\n")

# Start timer
start_time <- Sys.time()

# Run the main processing function
tryCatch({
  
  ReportInput <- Generate_Report_Inputs(
    client = client,
    samples_to_drop_pre_norm = samples_to_drop_pre_norm,
    samples_to_drop_post_norm = samples_to_drop_post_norm,
    mzmine_version = mzmine_version,
    ReferenceLevel = ReferenceLevel,
    Input = Input,
    contrast_var = contrast_var,
    num_meta = num_meta,
    SECIM_column = SECIM_column,
    anova_formula = anova_formula,
    lm_model = lm_model,
    test_type = test_type,
    subset = subset,
    metid_DB_file = metid_DB_file,
    paired = paired,
    batch_correct = batch_correct,
    rowNorm = rowNorm,
    transNorm = transNorm,
    scaleNorm = scaleNorm
  )
  
}, error = function(e) {
  cat("\n❌ ERROR during data processing:\n")
  cat(conditionMessage(e), "\n")
  stop("Pipeline halted due to processing error.")
})

# Calculate processing time
end_time <- Sys.time()
processing_time <- round(difftime(end_time, start_time, units = "mins"), 2)

cat("\n✓ Stage 1 complete\n")
cat(sprintf("Processing time: %.2f minutes\n\n", processing_time))

# ===========================
# SAVE RESULTS
# ===========================

cat("Saving results to disk...\n")

# Create Results directory if it doesn't exist
results_dir <- file.path("Results", client)
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Save serialized results using qs (fast compression)
output_file <- file.path(results_dir, paste0(client, "_ReportData.qs"))
qs::qsave(ReportInput, output_file, preset = "balanced")

cat(sprintf("✓ Results saved: %s\n", output_file))
cat(sprintf("File size: %.2f MB\n\n", file.size(output_file) / 1024^2))

# ===========================
# STAGE 2: REPORT GENERATION
# ===========================

cat("========================================\n")
cat("Stage 2: Generating Interactive Report\n")
cat("========================================\n\n")

# Prepare report parameters
report_params <- list(
  client = client,
  PI_name = PI_name,
  test_type = test_type,
  contrast_var = contrast_var,
  input_file = Input,
  data_file = output_file,
  processing_time = as.numeric(processing_time)
)

# Render Quarto report
report_template <- "Templates/Report_Template.qmd"

if (!file.exists(report_template)) {
  warning(sprintf("Report template not found: %s\nSkipping report generation.", report_template))
} else {
  
  tryCatch({
    
    quarto::quarto_render(
      input = report_template,
      output_format = "html",
      output_file = paste0(client, "_Report.html"),
      execute_params = report_params,
      execute_dir = results_dir
    )
    
    cat(sprintf("\n✓ Report generated: %s/%s_Report.html\n", results_dir, client))
    
  }, error = function(e) {
    cat("\n⚠ WARNING: Report generation failed\n")
    cat(conditionMessage(e), "\n")
    cat("Processed data is still available in the .qs file.\n")
  })
}

# ===========================
# CLEANUP
# ===========================

cat("\nCleaning up temporary files...\n")

# Remove MetaboAnalystR intermediate files
temp_files <- c(
  "data_orig.qs", "data_prefilter_iqr.csv", "data_proc.qs", 
  "prenorm.qs", "preproc.qs", "row_norm.qs", "complete_norm.qs"
)

for (file in temp_files) {
  if (file.exists(file)) {
    file.remove(file)
  }
}

# Remove temporary CSV files
metab_csv_files <- list.files(pattern = ".*metab.*\\.csv$", full.names = TRUE)
if (length(metab_csv_files) > 0) {
  file.remove(metab_csv_files)
}

cat("✓ Cleanup complete\n")

# ===========================
# PIPELINE SUMMARY
# ===========================

cat("\n========================================\n")
cat("PIPELINE COMPLETE\n")
cat("========================================\n\n")

cat(sprintf("Total runtime: %.2f minutes\n", as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
cat(sprintf("Results directory: %s/\n", results_dir))
cat(sprintf("Data file: %s\n", basename(output_file)))
cat(sprintf("Report file: %s_Report.html\n", client))
cat("\nTo re-render the report without reprocessing data:\n")
cat(sprintf("  quarto::quarto_render('Templates/Report_Template.qmd', execute_params = list(data_file = '%s'))\n", output_file))
cat("\n")
