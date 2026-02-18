#' Generate Metabolomics Report Inputs
#'
#' Master function to process both ion modes (Positive/Negative), perform statistical
#' analysis, and prepare data objects for downstream reporting via Quarto templates.
#'
#' @param client Character. Client/project identifier for naming outputs.
#' @param samples_to_drop_pre_norm Character vector. Sample names to exclude before normalization.
#' @param samples_to_drop_post_norm Character vector. Sample names to exclude after normalization.
#' @param mzmine_version Integer (2 or 3). MZmine version used for peak detection.
#' @param ReferenceLevel Character. Factor level to use as statistical reference (optional).
#' @param Input Character. Path to input Excel file (relative to project root).
#' @param contrast_var Character. Metadata column name for primary contrast (e.g., "Class").
#' @param num_meta Integer. Number of metadata columns in the dataset.
#' @param SECIM_column Character vector. SECIM column names for annotation.
#' @param anova_formula Formula. Model formula for ANOVA (if test_type="anova").
#' @param lm_model Formula. Linear mixed model formula (if test_type="lm"/"lme").
#' @param test_type Character. Statistical test: "t.test", "anova", "lm", "lme", or "nostats".
#' @param subset List of lists. Pairwise contrasts for subset analyses.
#' @param metid_DB_file Character. Metabolite ID database file (e.g., "kegg_ms1_database0.0.3.rda").
#' @param paired Logical. Whether samples are paired.
#' @param batch_correct Logical. Apply batch correction.
#' @param rowNorm Character. Row normalization method (default: "SumNorm").
#' @param transNorm Character. Transformation method (default: "LogNorm").
#' @param scaleNorm Character. Scaling method (default: "ParetoNorm").
#'
#' @return List. Named list containing combined results, mode-specific outputs, and metadata.
#'
#' @details
#' Normalization options:
#' - rowNorm: "QuantileNorm", "CompNorm", "SumNorm", "MedianNorm", "SpecNorm"
#' - transNorm: "LogNorm", "CrNorm" (Cubic Root)
#' - scaleNorm: "MeanCenter", "AutoNorm", "ParetoNorm", "RangeNorm"
#'
#' @export

Generate_Report_Inputs <- function(
  client,
  samples_to_drop_pre_norm = NULL,
  samples_to_drop_post_norm = NULL,
  mzmine_version,
  ReferenceLevel = NULL,
  Input,
  contrast_var,
  num_meta,
  SECIM_column,
  anova_formula = NULL,
  lm_model = NULL,
  test_type,
  subset = NULL,
  metid_DB_file,
  paired = FALSE,
  batch_correct = FALSE,
  rowNorm = "SumNorm",
  transNorm = "LogNorm",
  scaleNorm = "ParetoNorm"
) {
  
  # ===========================
  # Load Required Libraries
  # ===========================
  suppressPackageStartupMessages({
    library(dplyr)
    library(qs)
    library(lme4)
    library(impute)
    library(data.table)
    library(foreach)
    library(parallel)
    library(emmeans)
    library(broom)
    library(ggplot2)
    library(gridExtra)
    library(grid)
    library(stringr)
    library(nlme)
    library(readxl)
    library(omu)
    library(metid)
    library(MetaboAnalystR)
  })
  
  # Source internal scripts (relative paths)
  source("R/SECIM_Metabolomics.R")
  source("R/Norm_Plots.R")
  source("R/metid_SECIM-main/R/annotate_metabolites_mass_dataset.R")
  source("R/metid_SECIM-main/R/mzIdentify_mass_dataset.R")
  source("R/metid_SECIM-main/R/convert_mzmine2mass_dataset.R")
  source("R/SanityCheck.HRK.R")
  
  # ===========================
  # Internal Helper Function: Process Ion Mode Data
  # ===========================
  process_mode_data <- function(mode, Input, mzmine_version, samples_to_drop_pre_norm, ReferenceLevel) {
    
    cat(sprintf("\n========== Processing %s Mode ==========\n", mode))
    
    # Read metadata and peakdata
    metadata <- read_excel(Input, sheet = "Sample.data")
    sheet_name <- ifelse(mode == "Neg", "Peaktable.neg", "Peaktable.pos")
    peakdata <- read_excel(Input, sheet = sheet_name, col_names = FALSE)
    
    # Extract and clean header row
    header_row <- as.character(peakdata[1, ])
    dup_cols <- header_row[duplicated(header_row)]
    
    # Remove duplicated columns
    peakdata <- peakdata[-1, ]
    peakdata <- peakdata[, !duplicated(header_row)]
    colnames(peakdata) <- header_row[!duplicated(header_row)]
    
    if (length(dup_cols) > 0) {
      cat("Removed duplicated columns:", paste(dup_cols, collapse = ", "), "\n")
    }
    
    # ===========================
    # Standardize Peak Table Columns (MZmine Version-Specific)
    # ===========================
    if (mzmine_version == 3) {
      colnames(peakdata) <- c("id", "rt", "mz", "internal_ID", "spectral_ID", colnames(peakdata)[6:ncol(peakdata)])
    } else {
      colnames(peakdata) <- c("id", "mz", "rt", "compound", colnames(peakdata)[5:ncol(peakdata)])
      peakdata <- peakdata %>% dplyr::relocate("rt", .before = "mz")
    }
    
    # ===========================
    # Process Confidence Levels (MZmine 3 Only)
    # ===========================
    if ("internal_ID" %in% colnames(peakdata)) {
      peakdata$Confidence <- ""
      
      # Confidence level 1: internal_ID present, spectral_ID absent
      peakdata$Confidence <- ifelse(is.na(peakdata$spectral_ID) & !is.na(peakdata$internal_ID), 1, peakdata$Confidence)
      
      # Confidence level 2: spectral_ID present, internal_ID absent
      peakdata$Confidence <- ifelse(!is.na(peakdata$spectral_ID) & is.na(peakdata$internal_ID), 2, peakdata$Confidence)
      peakdata$spectral_ID <- ifelse(!is.na(peakdata$spectral_ID) & is.na(peakdata$internal_ID), 
                                     toupper(peakdata$spectral_ID), peakdata$spectral_ID)
      
      # Clear spectral_ID if both are present
      peakdata$spectral_ID <- ifelse(!is.na(peakdata$spectral_ID) & !is.na(peakdata$internal_ID), "", peakdata$spectral_ID)
      
      # Combine into single compound column
      peakdata$compound <- ifelse(!is.na(peakdata$internal_ID), peakdata$internal_ID, peakdata$spectral_ID)
      
      # Save confidence levels for post-processing
      pre_stats_conf <- peakdata %>% dplyr::select(id, Confidence)
      
      # Remove temporary columns
      peakdata <- peakdata %>% 
        dplyr::select(-internal_ID, -spectral_ID, -Confidence) %>%
        dplyr::relocate(compound, .after = mz)
      
    } else {
      pre_stats_conf <- NULL
    }
    
    # ===========================
    # Clean Column Names
    # ===========================
    # Remove SECIM identifiers (Q###_###_###_)
    colnames(peakdata) <- gsub("^Q[^_]+_[^_]+_[^_]+_(.*)", "\\1", colnames(peakdata))
    
    # Handle duplicate column names
    if (any(duplicated(colnames(peakdata)))) {
      duplicate_cols <- colnames(peakdata)[duplicated(colnames(peakdata))]
      colnames(peakdata)[5:ncol(peakdata)] <- make.names(colnames(peakdata)[5:ncol(peakdata)], unique = TRUE)
      warning("Duplicate column names found and made unique: ", paste(duplicate_cols, collapse = ", "))
    }
    
    # Replace special characters (MetaboAnalyst compatibility)
    colnames(peakdata) <- gsub("-", "_", colnames(peakdata))
    metadata$Sample.Name <- gsub("-", "_", metadata$Sample.Name)
    
    # ===========================
    # Filter Adducts and Clean Compound Names (MZmine 2)
    # ===========================
    if (mzmine_version == 2) {
      peakdata <- peakdata %>%
        dplyr::filter(!grepl("adduct|Complex", compound, ignore.case = TRUE))
      
      # Remove numeric suffixes from compound names (except LYSO* compounds)
      exclude_rows <- grepl("^LYSO", peakdata$compound)
      peakdata$compound[!exclude_rows] <- gsub(":\\s*\\d+\\.?\\d*", "", peakdata$compound[!exclude_rows])
    }
    
    # ===========================
    # Process Metadata
    # ===========================
    # Prepend reference level for sorting (legacy compatibility)
    if (!is.null(ReferenceLevel)) {
      metadata$Class <- gsub(ReferenceLevel, paste0("Zref_", ReferenceLevel), metadata$Class)
    }
    
    # Drop pre-normalization samples
    if (length(samples_to_drop_pre_norm) > 0) {
      samples_to_drop_pre_norm <- gsub("-", "_", samples_to_drop_pre_norm)
      metadata <- metadata %>% dplyr::filter(!Sample.Name %in% samples_to_drop_pre_norm)
    }
    
    # ===========================
    # Map Sample Names to Peak Table Columns
    # ===========================
    problematic_samples <- c()
    name_mapping <- sapply(metadata$Sample.Name, function(sample_name) {
      
      # Pattern 1: Bracketed sample names [SampleName]
      if (any(grepl("\\[", colnames(peakdata)))) {
        matched_colname <- grep(paste0("\\[", sample_name, "\\]"), colnames(peakdata), value = TRUE)
        
      # Pattern 2: Sample names starting with numbers
      } else if (any(grepl("^[0-9]+_", metadata$Sample.Name))) {
        matched_colname <- grep(paste0("^", sample_name, "(_|$)"), colnames(peakdata), value = TRUE)
        
      # Pattern 3: Standard numeric prefix
      } else {
        matched_colname <- grep(paste0("^[0-9]+_", sample_name, "_"), colnames(peakdata), value = TRUE)
      }
      
      if (length(matched_colname) == 1) {
        return(matched_colname)
      } else {
        warning("Multiple or no matches found for sample name: ", sample_name)
        problematic_samples <<- c(problematic_samples, sample_name)
        return(NA)
      }
    })
    
    name_mapping <- na.omit(name_mapping)
    metadata <- metadata %>% dplyr::filter(!Sample.Name %in% problematic_samples)
    
    # Subset and rename peak table columns
    columns_to_keep <- c(colnames(peakdata)[1:4], name_mapping)
    peakdata <- peakdata[, colnames(peakdata) %in% columns_to_keep, drop = FALSE]
    
    for (i in 5:ncol(peakdata)) {
      colname <- colnames(peakdata)[i]
      if (colname %in% name_mapping) {
        colnames(peakdata)[i] <- names(name_mapping)[name_mapping == colname]
      }
    }
    
    # ===========================
    # Construct Final Data Object
    # ===========================
    peaks <- peakdata[, metadata$Sample.Name]
    data <- data.frame(t(rbind(
      c(colnames(metadata), peakdata$id),
      merge(metadata, data.frame(t(peaks)), by.x = "Sample.Name", by.y = 0)
    )))
    
    colnames(data) <- data[1, ]
    data <- data[-1, ]
    data <- data %>% dplyr::rename("id" = "Sample.Name")
    rownames(data) <- 1:nrow(data)
    
    cat(sprintf("Processed %d peaks across %d samples\n", nrow(peakdata), ncol(peaks)))
    
    return(list(
      data = data,
      peakdata = peakdata,
      metadata = metadata,
      pre_stats_conf = pre_stats_conf
    ))
  }
  
  # ===========================
  # Process Negative Mode
  # ===========================
  neg_data <- process_mode_data("Neg", Input, mzmine_version, samples_to_drop_pre_norm, ReferenceLevel)
  
  neg.output <- SECIM_Metabolomics(
    dataset = neg_data$data,
    peakdata = neg_data$peakdata,
    num_meta = num_meta,
    original_data = neg_data$data,
    contrast_var = contrast_var,
    subset = subset,
    anova_formula = anova_formula,
    SECIM_column = SECIM_column,
    lm_model = lm_model,
    test_type = test_type,
    emmeans_var = contrast_var,
    mode = "Neg",
    metid_DB_file = metid_DB_file,
    client = client,
    metadata = neg_data$metadata,
    paired = paired,
    batch_correct = batch_correct,
    samples_to_drop_post_norm = samples_to_drop_post_norm,
    rowNorm = rowNorm,
    transNorm = transNorm,
    scaleNorm = scaleNorm
  )
  
  # ===========================
  # Process Positive Mode
  # ===========================
  pos_data <- process_mode_data("Pos", Input, mzmine_version, samples_to_drop_pre_norm, ReferenceLevel)
  
  pos.output <- SECIM_Metabolomics(
    dataset = pos_data$data,
    peakdata = pos_data$peakdata,
    num_meta = num_meta,
    original_data = pos_data$data,
    contrast_var = contrast_var,
    subset = subset,
    anova_formula = anova_formula,
    SECIM_column = SECIM_column,
    lm_model = lm_model,
    test_type = test_type,
    emmeans_var = contrast_var,
    mode = "Pos",
    metid_DB_file = metid_DB_file,
    client = client,
    metadata = pos_data$metadata,
    paired = paired,
    batch_correct = batch_correct,
    samples_to_drop_post_norm = samples_to_drop_post_norm,
    rowNorm = rowNorm,
    transNorm = transNorm,
    scaleNorm = scaleNorm
  )
  
  # ===========================
  # Post-Normalization Sample Filtering
  # ===========================
  if (!is.null(samples_to_drop_post_norm)) {
    samples_to_drop_post_norm <- gsub("-", "_", samples_to_drop_post_norm)
    neg_data$metadata <- neg_data$metadata %>% dplyr::filter(!Sample.Name %in% samples_to_drop_post_norm)
  }
  
  # ===========================
  # Combine and Deduplicate Results
  # ===========================
  combined_results <- rbind(neg.output[[1]], pos.output[[1]])
  combined_results$temp.lc.Metabolite <- tolower(combined_results$compound)
  
  # Identify high-confidence KEGG IDs (Level 1)
  HP.KEGG.list <- combined_results %>%
    dplyr::filter(Level == "1") %>%
    dplyr::select(KEGG) %>%
    unlist() %>%
    na.omit()
  
  combined_results <- combined_results %>%
    dplyr::mutate(status = case_when(KEGG %in% HP.KEGG.list ~ "HP"))
  
  orig.combined_results <- combined_results  # Preserve for ANOVA post-processing
  
  # ===========================
  # Deduplicate by Test Type
  # ===========================
  if (test_type == "t.test") {
    # High-confidence KEGG IDs: prioritize Level 1, then p-value
    temp.HP.Keggs <- combined_results %>%
      dplyr::filter(status == "HP") %>%
      dplyr::arrange(KEGG, Level, p.value) %>%
      dplyr::distinct(KEGG, .keep_all = TRUE)
    
    # Lower-confidence KEGG IDs: prioritize mz.match.score, then p-value
    temp.notHP.Keggs <- combined_results %>%
      dplyr::filter(is.na(status), !is.na(KEGG)) %>%
      dplyr::arrange(KEGG, desc(mz.match.score), p.value) %>%
      dplyr::distinct(KEGG, .keep_all = TRUE)
    
    # Process compound names without KEGG IDs
    HP.names.list <- combined_results %>%
      dplyr::filter(Level == "1") %>%
      dplyr::select(compound) %>%
      unlist() %>%
      na.omit()
    
    combined_results <- combined_results %>%
      dplyr::mutate(status = case_when(compound %in% HP.names.list ~ "HP"))
    
    temp.HP.names <- combined_results %>%
      dplyr::filter(status == "HP", is.na(KEGG)) %>%
      dplyr::arrange(temp.lc.Metabolite, Level, p.value) %>%
      dplyr::distinct(temp.lc.Metabolite, .keep_all = TRUE)
    
    temp.notHP.names <- combined_results %>%
      dplyr::filter(is.na(status)) %>%
      dplyr::arrange(temp.lc.Metabolite, desc(mz.match.score), p.value) %>%
      dplyr::distinct(temp.lc.Metabolite, .keep_all = TRUE)
  }
  
  if (test_type == "nostats") {
    # Sort by log2FC instead of p-value
    temp.HP.Keggs <- combined_results %>%
      dplyr::filter(status == "HP") %>%
      dplyr::arrange(KEGG, Level, mz.match.score) %>%
      dplyr::distinct(KEGG, .keep_all = TRUE)
    
    temp.notHP.Keggs <- combined_results %>%
      dplyr::filter(is.na(status), !is.na(KEGG)) %>%
      dplyr::arrange(KEGG, desc(mz.match.score), desc(log2FC)) %>%
      dplyr::distinct(KEGG, .keep_all = TRUE)
    
    HP.names.list <- combined_results %>%
      dplyr::filter(Level == "1") %>%
      dplyr::select(compound) %>%
      unlist() %>%
      na.omit()
    
    combined_results <- combined_results %>%
      dplyr::mutate(status = case_when(compound %in% HP.names.list ~ "HP"))
    
    temp.HP.names <- combined_results %>%
      dplyr::filter(status == "HP", is.na(KEGG)) %>%
      dplyr::arrange(temp.lc.Metabolite, Level, desc(log2FC)) %>%
      dplyr::distinct(temp.lc.Metabolite, .keep_all = TRUE)
    
    temp.notHP.names <- combined_results %>%
      dplyr::filter(is.na(status)) %>%
      dplyr::arrange(temp.lc.Metabolite, desc(mz.match.score), desc(log2FC)) %>%
      dplyr::distinct(temp.lc.Metabolite, .keep_all = TRUE)
  }
  
  if (test_type %in% c("lm", "lme")) {
    temp.HP.Keggs <- combined_results %>%
      dplyr::filter(status == "HP") %>%
      dplyr::arrange(KEGG, Level, p.value) %>%
      dplyr::distinct(KEGG, .keep_all = TRUE)
    
    temp.notHP.Keggs <- combined_results %>%
      dplyr::filter(is.na(status), !is.na(KEGG)) %>%
      dplyr::arrange(KEGG, desc(mz.match.score), p.value) %>%
      dplyr::distinct(KEGG, .keep_all = TRUE)
    
    temp.HP.names <- combined_results %>%
      dplyr::filter(status == "HP", is.na(KEGG)) %>%
      dplyr::arrange(temp.lc.Metabolite, Level, p.value) %>%
      dplyr::distinct(temp.lc.Metabolite, .keep_all = TRUE)
    
    temp.notHP.names <- combined_results %>%
      dplyr::filter(is.na(status)) %>%
      dplyr::arrange(temp.lc.Metabolite, desc(mz.match.score), p.value) %>%
      dplyr::distinct(temp.lc.Metabolite, .keep_all = TRUE)
  }
  
    if (test_type == "anova") {
    # For ANOVA: use adj.p.value instead of p.value
    temp.HP.Keggs <- combined_results %>%
      dplyr::filter(status == "HP") %>%
      dplyr::arrange(KEGG, Level, adj.p.value) %>%
      dplyr::distinct(KEGG, .keep_all = TRUE)
    
    temp.notHP.Keggs <- combined_results %>%
      dplyr::filter(is.na(status), !is.na(KEGG)) %>%
      dplyr::arrange(KEGG, desc(mz.match.score), adj.p.value) %>%
      dplyr::distinct(KEGG, .keep_all = TRUE)
    
    temp.HP.names <- combined_results %>%
      dplyr::filter(status == "HP", is.na(KEGG)) %>%
      dplyr::arrange(temp.lc.Metabolite, Level, adj.p.value) %>%
      dplyr::distinct(temp.lc.Metabolite, .keep_all = TRUE)
    
    temp.notHP.names <- combined_results %>%
      dplyr::filter(is.na(status)) %>%
      dplyr::arrange(temp.lc.Metabolite, desc(mz.match.score), adj.p.value) %>%
      dplyr::distinct(temp.lc.Metabolite, .keep_all = TRUE)
  }
  
  # ===========================
  # Reassemble Deduplicated Results
  # ===========================
  combined_results <- rbind(temp.HP.Keggs, temp.notHP.Keggs, temp.HP.names, temp.notHP.names) %>%
    dplyr::select(-status, -temp.lc.Metabolite) %>%
    dplyr::distinct()
  
  combined_results <- combined_results %>%
    dplyr::mutate(Confidence = case_when(Level == 3 ~ "LOW CONFIDENCE ID", .default = ""))
  
  combined_results$row.mode <- paste(combined_results$id, combined_results$mode, sep = "_")
  
  # Filter original results to retained peaks
  orig.combined_results$row.mode <- paste(orig.combined_results$id, orig.combined_results$mode, sep = "_")
  retained_peaks <- combined_results$row.mode
  orig.combined_results <- orig.combined_results %>% dplyr::filter(row.mode %in% retained_peaks)
  combined_results <- orig.combined_results
  
  # ===========================
  # Merge Raw and Normalized Intensities
  # ===========================
  # Raw intensities
  combined_proc <- rbind(neg.output[[3]], pos.output[[3]]) %>%
    dplyr::select(c(id, mode, neg_data$metadata$Sample.Name))
  combined_proc$row.mode <- paste(combined_proc$id, combined_proc$mode, sep = "_")
  combined_results <- merge(combined_results, combined_proc, by = "row.mode")
  
  # Normalized intensities
  combined_norm <- rbind(neg.output[[4]], pos.output[[4]]) %>%
    dplyr::select(c(id, mode, neg_data$metadata$Sample.Name))
  combined_norm$row.mode <- paste(combined_norm$id, combined_norm$mode, sep = "_")
  
  # Prefix normalized columns with "norm."
  colnames(combined_norm) <- gsub(
    paste0("(", paste(neg_data$metadata$Sample.Name, collapse = "|"), ")"),
    "norm.\\1",
    colnames(combined_norm)
  )
  combined_results <- merge(combined_results, combined_norm, by = "row.mode")
  
  # ===========================
  # Assemble Final Output List
  # ===========================
  Client_Data_Download <- list(
    report_results = combined_results
  )
  
  # Append mode-specific outputs
  Client_Data_Download <- c(Client_Data_Download, pos.output, neg.output)
  
  # Remove empty placeholders (t-test studies)
  Client_Data_Download <- Client_Data_Download[!names(Client_Data_Download) %in% c("PosEmpty", "NegEmpty")]
  
  # ===========================
  # Clean Redundant Columns
  # ===========================
  for (i in 1:length(Client_Data_Download)) {
    tryCatch({
      # Remove duplicate mode columns
      if ("mode.x" %in% colnames(Client_Data_Download[[i]])) {
        Client_Data_Download[[i]] <- Client_Data_Download[[i]] %>% dplyr::select(-c(mode.x, mode.y))
      }
      
      # Remove duplicate row name columns
      if ("Row.names.y" %in% colnames(Client_Data_Download[[i]])) {
        Client_Data_Download[[i]] <- Client_Data_Download[[i]] %>% dplyr::select(-c(Row.names.x, Row.names.y))
      }
      
      # Remove status column
      if ("status" %in% colnames(Client_Data_Download[[i]])) {
        Client_Data_Download[[i]] <- Client_Data_Download[[i]] %>% dplyr::select(-status)
      }
      
      # Remove duplicate ID columns
      if ("id.y" %in% colnames(Client_Data_Download[[i]])) {
        Client_Data_Download[[i]] <- Client_Data_Download[[i]] %>% dplyr::select(-c(id.x, id.y))
      }
      
      # Remove temporary lowercase metabolite column
      if ("temp.lc.Metabolite" %in% colnames(Client_Data_Download[[i]])) {
        Client_Data_Download[[i]] <- Client_Data_Download[[i]] %>% dplyr::select(-temp.lc.Metabolite)
      }
      
      # Standardize ID column name to "row.ID"
      if ("row.ID" %in% colnames(Client_Data_Download[[i]]) & "id" %in% colnames(Client_Data_Download[[i]])) {
        Client_Data_Download[[i]] <- Client_Data_Download[[i]] %>% dplyr::select(-id)
      } else if ("id" %in% colnames(Client_Data_Download[[i]])) {
        Client_Data_Download[[i]] <- Client_Data_Download[[i]] %>% dplyr::rename(row.ID = id)
      }
      
      # Add Sample.name to metadata
      if (i == which(names(Client_Data_Download) == "metadata")) {
        Client_Data_Download[[i]]["Sample.name"] <- rownames(Client_Data_Download[[i]])
        Client_Data_Download[[i]] <- Client_Data_Download[[i]] %>% relocate(Sample.name)
      }
      
      # Remove empty columns
      empty_columns <- colSums(is.na(Client_Data_Download[[i]]) | Client_Data_Download[[i]] == "") == nrow(Client_Data_Download[[i]])
      Client_Data_Download[[i]] <- Client_Data_Download[[i]][, !empty_columns]
      
    }, error = function(e) {
      message(sprintf("Skipping cleanup for element %d: %s", i, conditionMessage(e)))
    })
  }
  
  # ===========================
  # Apply Pre-Stats Confidence Levels (MZmine 3)
  # ===========================
  if (!is.null(pos_data$pre_stats_conf) && !is.null(neg_data$pre_stats_conf)) {
    pos_data$pre_stats_conf$Mode <- "Pos"
    neg_data$pre_stats_conf$Mode <- "Neg"
    
    Pre_stats_conf <- rbind(pos_data$pre_stats_conf, neg_data$pre_stats_conf)
    Pre_stats_conf$id <- as.character(Pre_stats_conf$id)
    
    for (n in 1:length(Client_Data_Download)) {
      for (i in 1:nrow(Pre_stats_conf)) {
        row_value <- Pre_stats_conf$id[i]
        mode_value <- Pre_stats_conf$Mode[i]
        
        matching_row <- Client_Data_Download[[n]]$id == row_value & 
                        Client_Data_Download[[n]]$mode == mode_value
        
        if (any(matching_row) && Pre_stats_conf$Confidence[i] != "") {
          Client_Data_Download[[n]]$Level[matching_row] <- Pre_stats_conf$Confidence[i]
        }
      }
    }
  }
  
  cat("\n========== Pipeline Complete ==========\n")
  cat(sprintf("Combined results: %d metabolites\n", nrow(combined_results)))
  cat(sprintf("Output contains %d data elements\n", length(Client_Data_Download)))
  
  return(Client_Data_Download)
}

# ===========================
# NOTE: Temporary file cleanup removed
# ===========================
# The original script included hardcoded absolute paths for cleanup:
# /blue/timgarrett/hkates/*.qs, *.csv
# These should be handled by the master Run_Pipeline.R script
# using relative paths in a temporary Results/ subdirectory.