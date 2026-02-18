#' Perform Metabolite Set Enrichment Analysis (MSEA) Using SMPDB
#'
#' Conducts over-representation analysis (ORA) for significantly changed metabolites
#' using the Small Molecule Pathway Database (SMPDB) via MetaboAnalystR.
#'
#' @param results_df Data frame. Statistical results containing compound names, KEGG IDs, and p-values.
#' @param contrast Character. Name of the contrast to analyze (e.g., "Treatment - Control").
#' @param p_threshold Numeric. P-value threshold for significance (default: 0.05).
#' @param p_type Character. P-value column name ("p.value" or "adj.p.value").
#' @param output_dir Character. Directory to save output files (default: "Results/MSEA").
#' @param plot_filename Character. Filename for output plot (default: "MSEA_plot.png").
#' @param dpi Numeric. Resolution for output plot (default: 300).
#'
#' @return List containing:
#'   - `result_table`: Data frame with enrichment results
#'   - `plot_path`: Path to saved plot
#'   - `mSet`: MetaboAnalystR mSet object
#'
#' @details
#' This function:
#' 1. Extracts high-confidence significantly changed metabolites
#' 2. Removes common adducts from compound names
#' 3. Maps compound names to SMPDB pathways
#' 4. Calculates enrichment scores (hypergeometric test)
#' 5. Generates a bar plot of enriched pathways
#'
#' Intermediate files are created in a temporary directory and cleaned up automatically.
#'
#' @export

perform_msea_smpdb <- function(results_df,
                                contrast,
                                p_threshold = 0.05,
                                p_type = "adj.p.value",
                                output_dir = "Results/MSEA",
                                plot_filename = "MSEA_plot.png",
                                dpi = 300) {
  
  # ===========================
  # Validate Inputs
  # ===========================
  
  if (!requireNamespace("MetaboAnalystR", quietly = TRUE)) {
    stop("Package 'MetaboAnalystR' is required for MSEA analysis.")
  }
  
  if (!contrast %in% results_df$contrast) {
    stop(sprintf("Contrast '%s' not found in results_df", contrast))
  }
  
  if (!p_type %in% colnames(results_df)) {
    stop(sprintf("Column '%s' not found in results_df", p_type))
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ===========================
  # Extract Significant Metabolites
  # ===========================
  
  cat(sprintf("Extracting significant metabolites for contrast: %s\n", contrast))
  
  cmpd_vec <- results_df %>%
    filter(contrast == !!contrast) %>%
    filter(!!sym(p_type) < p_threshold) %>%
    filter(ID_confidence == "High") %>%
    dplyr::select(compound) %>%
    unlist()
  
  if (length(cmpd_vec) == 0) {
    warning(sprintf("No significant high-confidence metabolites found for contrast: %s", contrast))
    return(NULL)
  }
  
  cat(sprintf("Found %d significant metabolites\n", length(cmpd_vec)))
  
  # ===========================
  # Remove Adducts
  # ===========================
  
  # Define adducts (ordered longest to shortest)
  adducts <- c(
    "_M\\+HCOO", "_M\\+CH3COO", "_M\\+NH4", "_M\\+Na", "_M-H", "_M\\+K", "_M\\+H",
    "\\+HCOO", "\\+CH3COO", "\\+NH4", "\\+Na", "\\+K", "\\+H", "-H", "\\+Cl", "-Cl"
  )
  
  remove_adducts <- function(compound, adduct_list) {
    for (adduct in adduct_list) {
      pattern <- paste0("([+_-])?", adduct, "$")
      if (grepl(pattern, compound)) {
        compound <- sub(pattern, "", compound)
        break
      }
    }
    return(compound)
  }
  
  cmpd_vec_cleaned <- sapply(cmpd_vec, remove_adducts, adduct_list = adducts)
  
  cat(sprintf("Cleaned compound names (removed adducts)\n"))
  
  # ===========================
  # Perform MSEA
  # ===========================
  
  # Create temporary directory for intermediate files
  temp_dir <- tempfile("msea_")
  dir.create(temp_dir)
  original_wd <- getwd()
  setwd(temp_dir)
  
  tryCatch({
    
    # Initialize MetaboAnalyst data object
    mSet <- InitDataObjects("conc", "msetora", FALSE)
    
    # Setup and map metabolite names
    mSet <- Setup.MapData(mSet, cmpd_vec_cleaned)
    mSet <- CrossReferencing(mSet, "name")
    mSet <- CreateMappingResultTable(mSet)
    
    # Set pathway library (SMPDB)
    mSet <- SetMetabolomeFilter(mSet, FALSE)
    mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway", 2)
    
    # Calculate enrichment scores
    mSet <- CalculateHyperScore(mSet)
    
    # Generate plot
    plot_path <- file.path(output_dir, plot_filename)
    mSet <- PlotORA(mSet, plot_path, dpi = dpi, width = NA)
    
    # Read results table
    if (file.exists("msea_ora_result.csv")) {
      result_table <- read.csv("msea_ora_result.csv")
      colnames(result_table)[1] <- "SMPDB_Pathway"
    } else {
      warning("MSEA results file not found")
      result_table <- NULL
    }
    
    cat(sprintf("✓ MSEA complete. Plot saved: %s\n", plot_path))
    
  }, error = function(e) {
    cat(sprintf("❌ MSEA failed: %s\n", conditionMessage(e)))
    result_table <- NULL
    plot_path <- NULL
    mSet <- NULL
  }, finally = {
    # Restore working directory and clean up temp files
    setwd(original_wd)
    unlink(temp_dir, recursive = TRUE)
  })
  
  # ===========================
  # Return Results
  # ===========================
  
  return(list(
    result_table = result_table,
    plot_path = plot_path,
    mSet = mSet,
    n_input_compounds = length(cmpd_vec_cleaned),
    contrast = contrast
  ))
}


#' Wrapper for MSEA in Reporting Pipeline
#'
#' Convenience function that calls perform_msea_smpdb() for each contrast
#' in the analysis and returns results suitable for report generation.
#'
#' @param Client_Data_Download List. Pipeline output containing report_results.
#' @param contrasts Character vector. Contrasts to analyze.
#' @param p_type Character. P-value type to use.
#' @param client Character. Client name for output filenames.
#' @param output_dir Character. Base output directory.
#'
#' @return List of MSEA results (one per contrast)
#'
#' @export

run_msea_for_contrasts <- function(Client_Data_Download,
                                    contrasts,
                                    p_type = "adj.p.value",
                                    client = "Client",
                                    output_dir = "Results/MSEA") {
  
  msea_results <- list()
  
  for (i in seq_along(contrasts)) {
    cat(sprintf("\n========== MSEA for Contrast %d/%d ==========\n", i, length(contrasts)))
    
    plot_filename <- sprintf("%s_MSEA_contrast_%d.png", client, i)
    
    msea_results[[i]] <- perform_msea_smpdb(
      results_df = Client_Data_Download[["report_results"]],
      contrast = contrasts[i],
      p_threshold = 0.05,
      p_type = p_type,
      output_dir = output_dir,
      plot_filename = plot_filename
    )
    
    names(msea_results)[i] <- contrasts[i]
  }
  
  return(msea_results)
}