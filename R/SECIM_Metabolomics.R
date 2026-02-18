#' SECIM Metabolomics Statistical Analysis Pipeline
#'
#' Core function to perform normalization, statistical testing, metabolite annotation,
#' and group-wise fold-change calculations for both ion modes.
#'
#' @param dataset Data frame. Combined metadata and peak intensity data.
#' @param peakdata Data frame. Peak table from MZmine.
#' @param num_meta Integer. Number of metadata rows in dataset.
#' @param original_data Data frame. Copy of dataset for reference.
#' @param contrast_var Character. Metadata column for primary statistical contrast.
#' @param anova_formula Formula. Model formula for ANOVA (if test_type="anova").
#' @param lm_model Formula. Linear model formula (if test_type="lm"/"lme").
#' @param test_type Character. Statistical test: "t.test", "anova", "lm", "lme", or "nostats".
#' @param subset List of lists. Pairwise contrasts for subset analyses.
#' @param SECIM_column Character. SECIM column type for annotation ("ACE-PFP" or "Evospere-PFP").
#' @param emmeans_var Character. Variable for emmeans contrasts.
#' @param mode Character. Ion mode: "Pos" or "Neg".
#' @param metid_DB_file Character. Metabolite ID database file (relative to reference/).
#' @param client Character. Client/project identifier for file naming.
#' @param metadata Data frame. Sample metadata.
#' @param paired Logical. Whether samples are paired (for t-test).
#' @param batch_correct Character or NULL. Batch correction method: "ComBat", "limma", or NULL.
#' @param rowNorm Character. Row normalization method (default: "SumNorm").
#' @param transNorm Character. Transformation method (default: "LogNorm").
#' @param scaleNorm Character. Scaling method (default: "ParetoNorm").
#' @param samples_to_drop_post_norm Character vector. Sample names to exclude after normalization.
#'
#' @return Named list containing:
#'   - Statistical results (t-test/ANOVA/emmeans)
#'   - Model fit results (for ANOVA/lm/lme)
#'   - Processed data
#'   - Normalized data
#'   - Normalization plots (ggplot grobs)
#'   - Metadata (if mode="Neg")
#'
#' @details
#' This function performs the following steps:
#' 1. Data import and sanity checks via MetaboAnalystR
#' 2. Missing value imputation (KNN)
#' 3. Blank-feature filtering (IQR)
#' 4. Normalization (row/transform/scale)
#' 5. Optional batch correction (ComBat/limma)
#' 6. Statistical testing (t-test, ANOVA, linear models)
#' 7. Metabolite annotation (SECIM libraries + metid)
#' 8. Group means and fold-change calculations
#'
#' @export

SECIM_Metabolomics <- function(
  dataset,
  peakdata,
  num_meta,
  original_data,
  contrast_var,
  anova_formula,
  lm_model,
  test_type,
  subset,
  SECIM_column,
  emmeans_var,
  mode,
  metid_DB_file,
  client,
  metadata,
  paired = FALSE,
  batch_correct = NULL,
  rowNorm = "SumNorm",
  transNorm = "LogNorm",
  scaleNorm = "ParetoNorm",
  samples_to_drop_post_norm = NULL
) {
  
  # ===========================
  # Setup Temporary Directory
  # ===========================
  temp_dir <- file.path("Results", "temp", paste0(client, "_", mode))
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
  }
  
  # ===========================
  # Load Metabolite Database
  # ===========================
  before <- ls()
  metid_DB_path <- file.path("reference", metid_DB_file)
  
  if (!file.exists(metid_DB_path)) {
    stop(sprintf("Metabolite database not found: %s", metid_DB_path))
  }
  
  load(metid_DB_path)
  metid_DB <- setdiff(ls(), before)[2]
  cat(sprintf("Loaded metabolite database: %s\n", metid_DB))
  
  # ===========================
  # Prepare Data for MetaboAnalystR
  # ===========================
  dataset[is.na(dataset)] <- 0
  
  # Extract metadata rows
  md <- data.frame(t(data.frame(dataset[1:num_meta, ])))
  colnames(md) <- md[1, ]
  md <- md %>% dplyr::slice(-1)
  rownames(md) <- gsub("^X", "", rownames(md))
  
  # ===========================
  # MetaboAnalystR Workflow
  # ===========================
  csv_path <- file.path(temp_dir, paste0(client, "_", mode, "_metab.in.csv"))
  
  if (num_meta == 1) {
    dataset[dataset == 0] <- NA
    write.csv(file = csv_path, dataset, row.names = FALSE)
    
    mSet <- InitDataObjects("pktable", "stat", FALSE)
    mSet <- Read.TextData(mSet, csv_path, "colu", "disc")
    mSet <- SanityCheckDataHRK(mSet)
    
    # Extract ordered metadata
    md <- data.frame(mSet[["dataSet"]][["meta.info"]])
    rownames(md) <- mSet[["dataSet"]][["url.smp.nms"]]
    
  } else {
    dataset[dataset == 0] <- NA
    dataset <- dataset[-c(1:num_meta), ]
    write.csv(file = csv_path, dataset, row.names = FALSE)
    
    meta_csv_path <- file.path(temp_dir, paste0(client, "_", mode, "_metab.meta.in.csv"))
    write.csv(file = meta_csv_path, md, row.names = TRUE)
    
    mSet <- InitDataObjects("pktable", "mf", FALSE)
    mSet <- SetDesignType(mSet, "multi")
    mSet <- Read.TextDataTs(mSet, csv_path, "colmf")
    mSet <- ReadMetaData(mSet, meta_csv_path)
    mSet <- SanityCheckDataHRK(mSet)
    
    # Extract ordered metadata
    md <- data.frame(mSet[["dataSet"]][["meta.info"]])
  }
  
  # ===========================
  # Data Preprocessing
  # ===========================
  cat("Removing features with >50% missing values...\n")
  mSet <- RemoveMissingPercent(mSet, percent = 0.5)
  
  cat("Imputing missing values (KNN)...\n")
  mSet <- ImputeMissingVar(mSet, method = "knn_var")
  SanityCheckDataHRK(mSet)
  
  cat("Filtering features by IQR...\n")
  mSet <- FilterVariable(mSet, var.filter = "iqr", "F", 10)
  
  # ===========================
  # Normalization
  # ===========================
  cat(sprintf("Normalizing: %s / %s / %s\n", rowNorm, transNorm, scaleNorm))
  mSet <- PreparePrenormData(mSet)
  mSet <- Normalization(mSet, rowNorm = rowNorm, transNorm = transNorm, 
                        scaleNorm = scaleNorm, ratio = FALSE, ratioNum = 20)
  
  # Read MetaboAnalystR output files (created in working directory)
  proc.data <- qs::qread("data_proc.qs")
  norm.data <- qs::qread("complete_norm.qs")
  
  # ===========================
  # Generate Normalization Plots
  # ===========================
  plots <- Norm_Plots(proc.data = proc.data, norm.data = norm.data)
  
  # ===========================
  # Prepare Final Normalized Data
  # ===========================
  data.final <- data.frame(t(norm.data))
  colnames(data.final) <- gsub("^X", "", colnames(data.final))
  data.final$id <- rownames(data.final)
  data.final <- data.final %>% dplyr::relocate(id)
  
  # Processed data
  data.proc <- data.frame(t(proc.data))
  colnames(data.proc) <- gsub("^X", "", colnames(data.proc))
  data.proc$id <- rownames(data.proc)
  data.proc <- data.proc %>% dplyr::relocate(id)
  
  # Add metadata rows back
  if (num_meta == 1) {
    data.final <- rbind(dataset[1:num_meta, ], data.final)
    data.proc <- rbind(dataset[1:num_meta, ], data.proc)
  } else {
    tmd <- data.frame(t(md))
    tmd$id <- ""
    tmd <- tmd %>% dplyr::relocate(id)
    colnames(tmd) <- gsub("^X", "", colnames(tmd))
    
    data.final <- rbind(tmd, data.final)
    data.proc <- rbind(tmd, data.proc)
  }
  
  # ===========================
  # Batch Correction (Optional)
  # ===========================
  if (!is.null(batch_correct)) {
    if (batch_correct == "ComBat") {
      cat("Applying ComBat batch correction...\n")
      library(sva)
      
      batch_info <- as.numeric(data.final[num_meta, -1])
      data_matrix <- as.matrix(data.final[-(1:num_meta), -1])
      rownames(data_matrix) <- data.final[-(1:num_meta), 1]
      data_matrix <- apply(data_matrix, 2, as.numeric)
      
      original_feature_names <- rownames(data_matrix)
      original_sample_names <- colnames(data_matrix)
      
      data_corrected <- ComBat(dat = data_matrix, batch = batch_info, 
                               mod = NULL, par.prior = TRUE, prior.plots = FALSE)
      rownames(data_corrected) <- original_feature_names
      colnames(data_corrected) <- original_sample_names
      
      data_final_corrected <- as.data.frame(data_corrected)
      data_final_corrected$id <- rownames(data.final)[(num_meta + 1):nrow(data.final)]
      data_final_corrected <- data_final_corrected %>% dplyr::relocate(id)
      data_final_corrected <- rbind(data.final[1:num_meta, ], data_final_corrected)
      
      data.final <- data_final_corrected
    }
    
    if (batch_correct == "limma") {
      cat("Applying limma batch correction...\n")
      library(limma)
      
      batch_info <- as.numeric(data.final[num_meta, -1])
      data_matrix <- as.matrix(data.final[-(1:num_meta), -1])
      rownames(data_matrix) <- data.final[-(1:num_meta), 1]
      data_matrix <- apply(data_matrix, 2, as.numeric)
      
      original_feature_names <- rownames(data_matrix)
      original_sample_names <- colnames(data_matrix)
      
      data_matrix_corrected <- removeBatchEffect(data_matrix, batch = batch_info)
      rownames(data_matrix_corrected) <- original_feature_names
      
      data_corrected_df <- as.data.frame(data_matrix_corrected)
      data_corrected_df$id <- original_feature_names
      data_corrected_df <- dplyr::relocate(data_corrected_df, id)
      
      data_final_corrected <- rbind(data.final[1:num_meta, ], data_corrected_df)
      data.final <- data_final_corrected
    }
  }
  
  # ===========================
  # Post-Normalization Sample Filtering
  # ===========================
  if (!is.null(samples_to_drop_post_norm)) {
    samples_to_drop_post_norm <- gsub("-", "_", samples_to_drop_post_norm)
    data.final.total <- data.final  # Save full dataset for PCA/download
    data.final <- data.final %>% dplyr::select(!samples_to_drop_post_norm)
  }
  
  # ===========================
  # Statistical Testing
  # ===========================
  
  ## T-TEST
  if (test_type == "t.test") {
    cat("Performing t-tests...\n")
    
    ttest <- foreach(i = (num_meta + 1):nrow(data.final), .packages = c("dplyr", "stats", "broom")) %do% {
      temp <- data.frame(t(data.final[c(1:num_meta, i), ]))
      
      if (num_meta == 1) {
        colnames(temp) <- c(contrast_var, "Metabolite")
      } else {
        colnames(temp) <- c(colnames(temp[1, ][1:num_meta]), "Metabolite")
      }
      
      temp <- temp[-1, ]
      temp$Metabolite <- as.numeric(temp$Metabolite)
      temp$rowID <- data.final$id[i]
      
      exp1 <- expr(Metabolite ~ !!ensym(contrast_var))
      
      if (is.null(subset)) {
        if (paired) {
          temp <- temp[order(temp$rowID, temp$Subject), ]
          group1 <- temp$Metabolite[temp$Class == levels(as.factor(temp[[contrast_var]]))[1]]
          group2 <- temp$Metabolite[temp$Class == levels(as.factor(temp[[contrast_var]]))[2]]
          
          if (length(group1) == length(group2)) {
            ttest.res <- tidy(t.test(group1, group2, paired = TRUE))
          } else {
            stop("Groups do not have the same number of observations for paired t-test.")
          }
          
          ttest.res$contrast <- paste(levels(as.factor(temp[[contrast_var]]))[1], "-", 
                                      levels(as.factor(temp[[contrast_var]]))[2])
        } else {
          ttest.res <- tidy(t.test(formula = eval(exp1), data = temp))
          ttest.res$contrast <- paste(levels(as.factor(temp[[contrast_var]]))[1], "-", 
                                      levels(as.factor(temp[[contrast_var]]))[2])
        }
        ttest.res
      } else {
        ttest.res <- list()
        for (n in 1:length(subset)) {
          if (paired) {
            temp_subset <- temp %>% dplyr::filter(Class %in% subset[[n]])
            temp_subset <- temp_subset[order(temp_subset$rowID, temp_subset$Subject), ]
            
            group1 <- temp_subset$Metabolite[temp_subset$Class == subset[[n]][1]]
            group2 <- temp_subset$Metabolite[temp_subset$Class == subset[[n]][2]]
            
            if (length(group1) == length(group2)) {
              ttest.res[[n]] <- tidy(t.test(group1, group2, paired = TRUE))
            } else {
              stop("Groups do not have the same number of observations for paired t-test.")
            }
            
            ttest.res[[n]]$contrast <- paste(subset[[n]][1], "-", subset[[n]][2])
          } else {
            temp_subset <- temp %>% dplyr::filter(Class %in% subset[[n]])
            ttest.res[[n]] <- tidy(t.test(formula = eval(exp1), data = temp_subset))
            ttest.res[[n]]$contrast <- paste(subset[[n]][1], "-", subset[[n]][2])
          }
        }
        ttest.res <- do.call("rbind", ttest.res)
        ttest.res$id <- temp$rowID[1]
        ttest.res
      }
    }
    
    ttest.results <- do.call("rbind", ttest)
    if (is.null(subset)) {
      ttest.results$id <- data.final$id[(num_meta + 1):nrow(data.final)]
    }
    ttest.results <- dplyr::relocate(ttest.results, id)
  }
  
  ## ANOVA / LINEAR MODELS
  if (test_type %in% c("anova", "lm", "lme")) {
    cat(sprintf("Performing %s analysis...\n", toupper(test_type)))
    
    if (test_type == "anova") {
      fit_emmeans <- foreach(i = (num_meta + 1):nrow(data.final), 
                             .packages = c("dplyr", "emmeans", "stats")) %do% {
        tempdf <- data.frame(t(data.final[c(1:num_meta, i), ]))
        
        if (num_meta == 1) {
          colnames(tempdf) <- c(tempdf[1, 1], "id")
        } else {
          colnames(tempdf) <- c(colnames(tempdf)[1:num_meta], "id")
        }
        
        tempdf <- tempdf[-1, ]
        tempdf$id <- as.numeric(tempdf$id)
        tempdf$Class <- as.factor(as.character(tempdf$Class))
        
        if ("ID" %in% colnames(tempdf)) {
          tempdf$ID <- as.factor(tempdf$ID)
        }
        
        fit <- do.call(aov, args = list(anova_formula, tempdf))
        emmeans_obj <- tidy(pairs(emmeans(fit, emmeans_var, data = tempdf)))
        
        list(fit = fit, emmeans = emmeans_obj)
      }
      
      fit <- lapply(fit_emmeans, function(x) x$fit)
      emmeans <- lapply(fit_emmeans, function(x) x$emmeans)
    }
    
    if (test_type == "lm") {
      fit_emmeans <- foreach(i = (num_meta + 1):nrow(data.final), 
                             .packages = c("dplyr", "stats", "emmeans")) %do% {
        temp <- data.frame(t(data.final[c(1:num_meta, i), ]))
        colnames(temp) <- c(temp[1, ][1:num_meta], "Metabolite")
        temp <- temp[-1, ]
        temp$Metabolite <- as.numeric(temp$Metabolite)
        
        fit <- lm(lm_model, data = temp)
        emmeans_obj <- emmeans(fit, specs = emmeans_var)
        pairwise_emmeans_obj <- tidy(pairs(emmeans_obj))
        
        list(fit = fit, emmeans = pairwise_emmeans_obj)
      }
      
      fit <- lapply(fit_emmeans, function(x) x$fit)
      emmeans <- lapply(fit_emmeans, function(x) x$emmeans)
    }
    
    if (test_type == "lme") {
      fit_emmeans <- foreach(i = (num_meta + 1):nrow(data.final), 
                             .packages = c("dplyr", "stats", "lmerTest", "emmeans", "broom")) %do% {
        temp <- data.frame(t(data.final[c(1:num_meta, i), ]))
        colnames(temp) <- c(colnames(temp)[1:num_meta], "id")
        temp <- temp[-1, ]
        temp$id <- as.numeric(temp$id)
        
        temp$Class <- as.factor(temp$Class)
        temp$Subject <- as.factor(temp$Subject)
        temp$Batch <- as.factor(temp$Batch)
        
        fit <- lmerTest::lmer(lm_model, data = temp)
        emmeans_obj <- emmeans(fit, specs = emmeans_var)
        pairwise_emmeans_obj <- tidy(pairs(emmeans_obj))
        
        list(fit = fit, emmeans = pairwise_emmeans_obj)
      }
      
      fit <- lapply(fit_emmeans, function(x) x$fit)
      emmeans <- lapply(fit_emmeans, function(x) x$emmeans)
    }
    
    # Tidy model fit results
    if (test_type %in% c("anova", "lm")) {
      fit.results <- lapply(fit, tidy)
    }
    
    if (test_type == "lme") {
      library(broom.mixed)
      fit.results <- lapply(fit, function(x) tidy(x))
    }
    
    names(fit.results) <- data.final$id[(num_meta + 1):nrow(data.final)]
    fit.results <- do.call("rbind", fit.results)
    fit.results$id <- rownames(fit.results)
    fit.results$id <- gsub("\\.[0-9\\+]", "", fit.results$id)
    
    names(emmeans) <- data.final$id[(num_meta + 1):nrow(data.final)]
    emmeans.results <- do.call("rbind", emmeans)
    emmeans.results$id <- rownames(emmeans.results)
    emmeans.results$id <- gsub("\\.[0-9\\+]", "", emmeans.results$id)
  }
  
  # ===========================
  # Calculate Group Means and Fold Changes
  # ===========================
  cat("Calculating group means and fold changes...\n")
  
    # Define contrasts and groups
  if (test_type == "t.test") {
    contrast_vec <- gsub(" ", "", levels(as.factor(ttest.results$contrast)))
    contrast_vec <- sapply(contrast_vec, function(x) gsub("[()]", "", x))
  }
  
  if (test_type %in% c("anova", "lm", "lme")) {
    contrast_vec <- gsub(" ", "", levels(as.factor(emmeans.results$contrast)))
    contrast_vec <- sapply(contrast_vec, function(x) gsub("[()]", "", x))
  }
  
  if (test_type == "nostats") {
    contrast_vec <- combn(levels(as.factor(metadata$Class)), 2, 
                          FUN = function(x) paste0(x[1], "-", x[2]), simplify = FALSE)
  }
  
  # Extract unique group names from contrasts
  group_vec <- unique(unlist(sapply(as.list(contrast_vec), function(x) str_split(x, "-"))))
  
  # Map samples to groups
  group_samples <- list()
  for (i in 1:length(group_vec)) {
    group_samples[[i]] <- rownames(md %>% dplyr::filter(!!as.symbol(contrast_var) == group_vec[[i]]))
  }
  
  # Calculate group means from processed data
  group_means <- list()
  for (i in 1:length(group_vec)) {
    group_means[[i]] <- data.frame(t(proc.data)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(!!group_vec[[i]] := mean(c_across(matches(group_samples[[i]]))))
  }
  
  # Assemble means matrix
  means <- data.frame(matrix(ncol = 0, nrow = nrow(data.frame(t(proc.data)))))
  for (i in 1:length(group_vec)) {
    means <- cbind(means, group_means[[i]][, ncol(group_means[[i]])])
  }
  rownames(means) <- colnames(proc.data)
  means$id <- as.numeric(rownames(means))
  
  # Calculate per-contrast fold changes
  contrast_fold_changes <- list()
  for (i in 1:length(contrast_vec)) {
    contrast_fold_changes[[i]] <- means %>%
      dplyr::rowwise() %>%
      dplyr::mutate(!!contrast_vec[[i]] := 
                      !!as.symbol(str_split(contrast_vec[[i]], "-")[[1]][1]) / 
                      !!as.symbol(str_split(contrast_vec[[i]], "-")[[1]][2])) %>%
      dplyr::select(c(!!contrast_vec[[i]], id))
    
    log2FC <- log2(contrast_fold_changes[[i]][contrast_vec[[i]]])
    contrast_fold_changes[[i]] <- cbind(contrast_fold_changes[[i]], log2FC)
    colnames(contrast_fold_changes[[i]]) <- c("FC", "id", "log2FC")
    contrast_fold_changes[[i]]["contrast"] <- contrast_vec[[i]]
    contrast_fold_changes[[i]] <- contrast_fold_changes[[i]][, c(2, 4, 1, 3)]
  }
  
  names(contrast_fold_changes) <- contrast_vec
  all_fold_changes <- do.call("rbind", contrast_fold_changes)
  means_FC <- merge(means, all_fold_changes, by = "id")
  means_FC$id <- as.character(means_FC$id)
  means_FC$contrast <- gsub("-", " - ", means_FC$contrast)
  
  # ===========================
  # Metabolite Annotation
  # ===========================
  cat("Annotating metabolites...\n")
  
  # Prepare peakdata for metid (expects old MZmine column names)
  peakdataformetid <- peakdata
  peakdataformetid <- peakdataformetid[, c(1, 3, 2, 4:ncol(peakdataformetid))]
  peakdataformetid[-c(1, 4)] <- lapply(peakdataformetid[-c(1, 4)], function(x) as.numeric(as.character(x)))
  
  colnames(peakdataformetid) <- c("row ID", "row m/z", "row retention time", 
                                  colnames(peakdata)[4:ncol(peakdata)])
  
  # Convert to metid mass_dataset object
  metid <- convet_mzmine2mass_dataset(
    x = peakdataformetid %>% dplyr::select(!compound),
    rt_unit = "minute"
  )
  
  # Perform MS1 annotation
  if (mode == "Pos") {
    metid <- annotate_metabolites_mass_dataset(
      object = metid,
      ms1.match.ppm = 5,
      rt.match.tol = 10001,
      polarity = "positive",
      database = get(metid_DB),
      column = "rp_custom",
      threads = 4,
      candidate.num = 1
    )
  } else if (mode == "Neg") {
    metid <- annotate_metabolites_mass_dataset(
      object = metid,
      ms1.match.ppm = 10,
      rt.match.tol = 10001,
      polarity = "negative",
      database = get(metid_DB),
      column = "rp_custom",
      threads = 4,
      candidate.num = 1
    )
  }
  
  metid.result <- merge(metid@variable_info, metid@annotation_table, by = "variable_id", all.x = TRUE)
  
  # Format retention time and mass
  metid.result <- metid.result %>%
    dplyr::mutate(
      rt = rt / 60,
      rt = format(round(rt, digits = 2), nsmall = 2),
      mz = format(round(mz, digits = 4), nsmall = 4)
    )
  
  # Get KEGG hierarchy
  metid.result <- metid.result %>% dplyr::rename(KEGG = KEGG.ID)
  metid.result <- as.data.frame(assign_hierarchy(
    count_data = metid.result,
    keep_unknowns = TRUE,
    identifier = "KEGG"
  ))
  
  # ===========================
  # SECIM Library Annotation
  # ===========================
  # Load appropriate library based on mode and column type
  library_files <- list(
    Pos = list(
      "Evospere-PFP" = "library_positive1_MK_QE2.csv",
      "ACE-PFP" = "Positive_ACE-PFPs.csv"
    ),
    Neg = list(
      "Evospere-PFP" = "library_negative1_MK02.csv",
      "ACE-PFP" = "Negative_ACE-PFP.csv"
    )
  )
  
  library_file <- library_files[[mode]][[SECIM_column]]
  if (is.null(library_file)) {
    stop(sprintf("Invalid SECIM_column: %s for mode %s", SECIM_column, mode))
  }
  
  library_path <- file.path("reference", library_file)
  if (!file.exists(library_path)) {
    stop(sprintf("SECIM library not found: %s", library_path))
  }
  
  KEGG.compound <- read.csv(library_path) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      KEGG = na_if(KEGG, "null"),
      KEGG = na_if(KEGG, "None"),
      KEGG = na_if(KEGG, "")
    )
  
  # Keep only first KEGG ID per compound
  KEGG.compound <- KEGG.compound %>%
    dplyr::group_by(name) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  peakdata.KEGG <- merge(peakdata, KEGG.compound, by.x = "compound", by.y = "name", all.x = TRUE)
  peakdata.KEGG <- as.data.frame(assign_hierarchy(
    count_data = peakdata.KEGG,
    keep_unknowns = TRUE,
    identifier = "KEGG"
  ))
  
  peakdata.KEGG <- peakdata.KEGG[, setdiff(colnames(peakdata.KEGG), metadata$Sample.Name)]
  
  # ===========================
  # Merge SECIM and metid Annotations
  # ===========================
  # Use SECIM annotation where available, metid for the rest
  rowid.SECIM.na <- dplyr::filter(data.frame(peakdata.KEGG), is.na(peakdata.KEGG$compound)) %>%
    dplyr::select(id) %>%
    unlist()
  
  metid.result <- data.frame(metid.result) %>% dplyr::filter(variable_id %in% rowid.SECIM.na)
  peakdata.KEGG <- data.frame(peakdata.KEGG) %>% dplyr::filter(!id %in% rowid.SECIM.na)
  
  # Standardize columns for binding
  peakdata.KEGG$Level <- 1  # SECIM annotations have confidence level 1
  metid.result <- metid.result %>% dplyr::rename("id" = "variable_id", "compound" = "Compound.name")
  
  metid.result$id <- as.double(metid.result$id)
  metid.result$mz <- as.double(metid.result$mz)
  metid.result$rt <- as.double(metid.result$rt)
  
  peakdata.KEGG$id <- as.double(peakdata.KEGG$id)
  peakdata.KEGG$mz <- as.double(peakdata.KEGG$mz)
  peakdata.KEGG$rt <- as.double(peakdata.KEGG$rt)
  
  peak_annotation <- dplyr::bind_rows(peakdata.KEGG, metid.result)
  
  # Use mz_rt as compound name if no annotation available
  peak_annotation <- peak_annotation %>%
    dplyr::mutate(
      compound = case_when(is.na(compound) ~ paste(mz, rt, sep = "_"), .default = compound),
      KEGG = case_when(KEGG == "" ~ NA, .default = KEGG)
    )
  
  peak_annotation$mode <- mode
  
  # ===========================
  # Assemble Final Output List
  # ===========================
  outputs_list <- list()
  
  if (test_type == "t.test") {
    outputs_list[[1]] <- merge(ttest.results, peak_annotation, by = "id", all.x = TRUE)
    outputs_list[[1]]$contrast <- gsub("[()]", "", outputs_list[[1]]$contrast)
    outputs_list[[1]] <- outputs_list[[1]] %>% inner_join(means_FC, by = c("id", "contrast"))
    outputs_list[[1]] <- outputs_list[[1]] %>% dplyr::relocate(contrast, compound)
    outputs_list[[1]]$adj.p.value <- p.adjust(outputs_list[[1]]$p.value, method = "fdr")
    outputs_list[[1]] <- outputs_list[[1]] %>% dplyr::relocate(adj.p.value, .after = p.value)
    outputs_list[[2]] <- "Empty for t-test"
  }
  
  if (test_type %in% c("lm", "anova", "lme")) {
    outputs_list[[1]] <- merge(emmeans.results, peak_annotation, by = "id", all.x = TRUE)
    outputs_list[[1]]$contrast <- gsub("[()]", "", outputs_list[[1]]$contrast)
    
    if (length(unique(emmeans.results$contrast)) > 1) {
      outputs_list[[1]] <- outputs_list[[1]] %>% inner_join(means_FC, by = c("id", "contrast"))
      outputs_list[[1]] <- outputs_list[[1]] %>% dplyr::relocate(contrast, compound)
    } else {
      outputs_list[[1]]$adj.p.value <- p.adjust(outputs_list[[1]]$p.value, method = "fdr")
      outputs_list[[1]] <- outputs_list[[1]] %>% inner_join(means_FC, by = c("id", "contrast"))
      outputs_list[[1]] <- outputs_list[[1]] %>% dplyr::relocate(contrast, compound)
    }
    
    outputs_list[[2]] <- merge(fit.results, peak_annotation, by = "id", all.x = TRUE)
    outputs_list[[2]] <- outputs_list[[2]] %>% dplyr::relocate(compound)
  }
  
  if (test_type == "nostats") {
    outputs_list[[1]] <- merge(means_FC, peak_annotation, by = "id", all.x = TRUE)
    outputs_list[[1]]$contrast <- gsub("[()]", "", outputs_list[[1]]$contrast)
    outputs_list[[1]] <- outputs_list[[1]] %>% dplyr::relocate(contrast, compound)
    outputs_list[[2]] <- "Empty for nostats"
  }
  
  # Add processed and normalized data
  outputs_list[[3]] <- merge(data.proc, peak_annotation, by = "id", all.x = TRUE)
  outputs_list[[3]] <- outputs_list[[3]][c(
    which(outputs_list[[3]]$id == "Class"),
    setdiff(seq_len(nrow(outputs_list[[3]])), which(outputs_list[[3]]$id == "Class"))
  ), ]
  outputs_list[[3]] <- outputs_list[[3]][c(
    which(is.na(outputs_list[[3]]$id)),
    setdiff(seq_len(nrow(outputs_list[[3]])), which(is.na(outputs_list[[3]]$id)))
  ), ]
  outputs_list[[3]] <- outputs_list[[3]] %>% dplyr::relocate(compound)
  
  outputs_list[[4]] <- merge(data.final, peak_annotation, by = "id", all.x = TRUE)
  outputs_list[[4]] <- outputs_list[[4]][c(
    which(outputs_list[[4]]$id == "Class"),
    setdiff(seq_len(nrow(outputs_list[[4]])), which(outputs_list[[4]]$id == "Class"))
  ), ]
  outputs_list[[4]] <- outputs_list[[4]][c(
    which(is.na(outputs_list[[4]]$id)),
    setdiff(seq_len(nrow(outputs_list[[4]])), which(is.na(outputs_list[[4]]$id)))
  ), ]
  outputs_list[[4]] <- outputs_list[[4]] %>% dplyr::relocate(compound)
  
  # Add total normalized data if samples were dropped post-normalization
  if (!is.null(samples_to_drop_post_norm)) {
    idx <- if (mode == "Neg") 8 else 7
    outputs_list[[idx]] <- merge(data.final.total, peak_annotation, by = "id", all.x = TRUE)
    outputs_list[[idx]] <- outputs_list[[idx]][c(
      which(outputs_list[[idx]]$id == "Class"),
      setdiff(seq_len(nrow(outputs_list[[idx]])), which(outputs_list[[idx]]$id == "Class"))
    ), ]
    outputs_list[[idx]] <- outputs_list[[idx]][c(
      which(is.na(outputs_list[[idx]]$id)),
      setdiff(seq_len(nrow(outputs_list[[idx]])), which(is.na(outputs_list[[idx]]$id)))
    ), ]
    outputs_list[[idx]] <- outputs_list[[idx]] %>% dplyr::relocate(compound)
  }
  
  # Add normalization plots
  outputs_list[[5]] <- arrangeGrob(
    plots[["plot_env"]][["p1"]], plots[["plot_env"]][["p3"]],
    plots[["plot_env"]][["p2"]], plots[["plot_env"]][["p4"]],
    ncol = 2, top = textGrob("Feature View")
  )
  
  outputs_list[[6]] <- arrangeGrob(
    plots[["plot_env"]][["p5"]], plots[["plot_env"]][["p7"]],
    plots[["plot_env"]][["p6"]], plots[["plot_env"]][["p8"]],
    ncol = 2, top = textGrob("Sample View")
  )
  
  # Add metadata (Negative mode only)
  if (mode == "Neg") {
    if (!is.null(samples_to_drop_post_norm)) {
      outputs_list[[7]] <- md %>% dplyr::filter(!rownames(.) %in% samples_to_drop_post_norm)
    } else {
      outputs_list[[7]] <- md
    }
  }
  
  # ===========================
  # Name Output Elements
  # ===========================
  if (mode == "Neg") {
    if (test_type == "t.test") {
      names(outputs_list) <- c(
        paste0(mode, ".ttest.metab"),
        paste0(mode, "Empty"),
        paste0(mode, ".processed.data"),
        paste0(mode, ".normalized.data"),
        paste0(mode, ".FeatureView"),
        paste0(mode, ".SampleView"),
        "metadata"
      )
    }
    if (test_type %in% c("anova", "lm", "lme")) {
      names(outputs_list) <- c(
        paste0(mode, ".emmeans.results.metab"),
        paste0(mode, ".fit.results.metab"),
        paste0(mode, ".processed.data"),
        paste0(mode, ".normalized.data"),
        paste0(mode, ".FeatureView"),
        paste0(mode, ".SampleView"),
        "metadata"
      )
    }
    if (test_type == "nostats") {
      names(outputs_list) <- c(
        paste0(mode, ".FCanlaysis.metab"),
        paste0(mode, "Empty"),
        paste0(mode, ".processed.data"),
        paste0(mode, ".normalized.data"),
        paste0(mode, ".FeatureView"),
        paste0(mode, ".SampleView"),
        "metadata"
      )
    }
    
    if (!is.null(samples_to_drop_post_norm)) {
      names(outputs_list)[8] <- paste0(mode, ".total.normalized.data")
    }
    
  } else {  # Positive mode
    if (test_type == "t.test") {
      names(outputs_list) <- c(
        paste0(mode, ".ttest.metab"),
        paste0(mode, "Empty"),
        paste0(mode, ".processed.data"),
        paste0(mode, ".normalized.data"),
        paste0(mode, ".FeatureView"),
        paste0(mode, ".SampleView")
      )
    }
    if (test_type %in% c("anova", "lm", "lme")) {
      names(outputs_list) <- c(
        paste0(mode, ".emmeans.results.metab"),
        paste0(mode, ".fit.results.metab"),
        paste0(mode, ".processed.data"),
        paste0(mode, ".normalized.data"),
        paste0(mode, ".FeatureView"),
        paste0(mode, ".SampleView")
      )
    }
    if (test_type == "nostats") {
      names(outputs_list) <- c(
        paste0(mode, ".FCanlaysis.metab"),
        paste0(mode, "Empty"),
        paste0(mode, ".processed.data"),
        paste0(mode, ".normalized.data"),
        paste0(mode, ".FeatureView"),
        paste0(mode, ".SampleView")
      )
    }
    
    if (!is.null(samples_to_drop_post_norm)) {
      names(outputs_list)[7] <- paste0(mode, ".total.normalized.data")
    }
  }
  
  cat(sprintf("\n========== %s Mode Analysis Complete ==========\n", mode))
  
  return(outputs_list)
}