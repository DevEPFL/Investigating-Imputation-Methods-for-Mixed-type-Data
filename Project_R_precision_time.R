# ---------------------------
# 0. Set Global Seed
# ---------------------------
set.seed(123)  

# Installer 'pacman' s'il n'est pas présent
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

# Charger (et installer si besoin) tous les packages requis
pacman::p_load(
  # Imputation des données
  VIM, mice, missForest, missMDA,
  
  # Manipulation des données
  tidyr, dplyr, purrr, readr,
  
  # Visualisation
  ggplot2, scales, patchwork,
  
  # Autres utilitaires
  readxl, clustMixType, DescTools,
)

# ---------------------------
# 1. Load Libraries & Functions
# ---------------------------
library(VIM); library(mice); library(missForest); library(missMDA); library(tidyr) 
library(dplyr); library(ggplot2); library(clustMixType); library(scales) 
library(readxl); library(purrr); library(patchwork); library(readr)
library(DescTools)
source("functions_clustering_test.R")

# ---------------------------
# 2. Core Functions
# ---------------------------

add_missing <- function(df, prop) {
  total_cells <- nrow(df) * ncol(df)
  na_count <- pmax(1, pmin(round(total_cells * prop), total_cells - 1))
  
  # Create mask matrix
  mask <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df))
  na_indices <- sample.int(total_cells, size = na_count)
  mask[na_indices] <- TRUE
  
  # Critical fix: Ensure no column becomes fully NA
  col_na_counts <- colSums(mask)
  problematic_cols <- which(col_na_counts == nrow(df))
  
  # Fix columns that would become fully NA
  for (col in problematic_cols) {
    # Randomly preserve 1 value in the column
    row_to_keep <- sample(1:nrow(df), 1)
    mask[row_to_keep, col] <- FALSE
  }
  
  # Ensure no row becomes fully NA
  row_na_counts <- rowSums(mask)
  problematic_rows <- which(row_na_counts == ncol(df))
  
  # Fix rows that would become fully NA
  for (row in problematic_rows) {
    # Randomly preserve 1 value in the row
    col_to_keep <- sample(1:ncol(df), 1)
    mask[row, col_to_keep] <- FALSE
  }
  
  df_na <- df
  df_na[mask] <- NA
  
  list(data = df_na, mask = mask)
}


impute_data <- function(df_missing) {
  # Ensure reproducibility for stochastic methods
  set.seed(123)
  
  # Create method vector for MICE matching column types
  mice_methods <- rep("", ncol(df_missing))
  message("1")
  
  # Assign methods based on column types
  mice_methods[sapply(df_missing, is.numeric)] <- "pmm"
  message("2")
  
  mice_methods[sapply(df_missing, is.factor)] <- "polyreg"
  message("3")
  
  # Handle binary factors with logistic regression
  binary_cols <- sapply(df_missing, function(x) is.factor(x) && nlevels(x) == 2)
  mice_methods[binary_cols] <- "logreg"
  message("4")
  
  # Ensure at least one valid method exists
  if(all(mice_methods == "")) stop("No valid MICE methods assigned")
  message("mice_methods non vide")
  message("5")
  
  ## Create a copy for imputation using droplevels()
  df_missing_impute <- df_missing  # Copie du dataframe original
  for(j in seq_along(df_missing_impute)) {
    if(is.factor(df_missing_impute[[j]])) {
      df_missing_impute[[j]] <- droplevels(df_missing_impute[[j]])
    }
  }
  ## -------------------------------------------------------------------
  
  # ## --- Vérification des variances des variables numériques ---
  # num_vars <- df_missing_impute[, sapply(df_missing_impute, is.numeric)]
  # variance_values <- apply(num_vars, 2, function(x) var(x, na.rm = TRUE))
  # message("Variances des variables numériques:")
  # print(variance_values)
  # ## -------------------------------------------------------------------
  
  ## --- Vérification des distributions et associations des variables catégoriques ---
  
  # Calcul de la matrice des associations avec Cramér's V 
  fact_vars_names <- names(df_missing_impute)[sapply(df_missing_impute, is.factor)]
  cramer_matrix <- matrix(NA, nrow = length(fact_vars_names), ncol = length(fact_vars_names),
                          dimnames = list(fact_vars_names, fact_vars_names))
  for(i in seq_along(fact_vars_names)) {
    for(j in seq_along(fact_vars_names)) {
      if(i == j) {
        cramer_matrix[i, j] <- 1
      } else {
        tbl <- table(df_missing_impute[[fact_vars_names[i]]], df_missing_impute[[fact_vars_names[j]]])
        cramer_matrix[i, j] <- DescTools::CramerV(tbl)
      }
    }
  }

  # Remplacer les NaN par 0, si nécessaire
  cramer_matrix[is.nan(cramer_matrix)] <- 0
  
  
  # --- Création d'une predictorMatrix adaptée ---
  # Quickpred pour générer une matrice initiale
  predMatrix <- quickpred(df_missing_impute, mincor = 0.3)
  
  # Pour les colonnes numériques : supprimer comme prédicteurs 
  # ceux présentant une corrélation trop élevée (threshold = 0.9)
  num_cols <- sapply(df_missing_impute, is.numeric)
  if(sum(num_cols) > 1) {
    # Calculer la matrice de corrélation sur les colonnes numériques
    cor_matrix <- abs(cor(df_missing_impute[, num_cols], use = "pairwise.complete.obs"))
    diag(cor_matrix) <- 0  # Ignorer la diagonale
    threshold <- 0.8
    high_cor_pairs <- which(cor_matrix > threshold, arr.ind = TRUE)
    if(nrow(high_cor_pairs) > 0) {
      for(i in seq_len(nrow(high_cor_pairs))) {
        # Récupérer le nom des variables en corrélation
        var1 <- colnames(cor_matrix)[high_cor_pairs[i, 1]]
        var2 <- colnames(cor_matrix)[high_cor_pairs[i, 2]]
        # Pour réduire la collinéarité, on retire var2 comme prédicteur de var1.
        if(var1 %in% colnames(predMatrix) && var2 %in% colnames(predMatrix)) {
          predMatrix[var1, var2] <- 0
        }
      }
    }
  }
  # ---
  
  # Suppression des prédicteurs trop fortement associés pour les facteurs
  threshold_cat <- 0.8  
  for(i in seq_along(fact_vars_names)) {
    for(j in seq_along(fact_vars_names)) {
      if(i != j && cramer_matrix[i, j] > threshold_cat) {
        if(fact_vars_names[i] %in% colnames(predMatrix) && fact_vars_names[j] %in% colnames(predMatrix)) {
          predMatrix[fact_vars_names[i], fact_vars_names[j]] <- 0
        }
      }
    }
  }
  ## -------------------------------------------------------------------
  
  imputations <- list(
    KNN = {
      start_time <- Sys.time()
      tryCatch({
        
        message("Running KNN on ", nrow(df_missing), " rows, ", ncol(df_missing), " columns")
        knn_data <- as.data.frame(df_missing)
        knn_data <- knn_data %>% mutate(across(where(is.character), as.factor))
        
        # Lancer kNN sans sous-sélection supplémentaire
        imp_df <- VIM::kNN(knn_data, imp_var = FALSE)
        
        message("6: KNN imputation successful")
        message("KNN completed. Dimensions: ", nrow(imp_df), " x ", ncol(imp_df))
        
        end_time <- Sys.time()
        list(imputed_data = imp_df, time = as.numeric(end_time - start_time, units = "secs"))

      }, error = function(e) {
        message("KNN failed: ", conditionMessage(e))
        # Créer un data frame vide de remplacement
        dummy_df <- as.data.frame(matrix(NA, nrow = nrow(df_missing), ncol = ncol(df_missing)))
        colnames(dummy_df) <- colnames(df_missing)
        for (col in names(df_missing)) {
          if (is.factor(df_missing[[col]])) {
            dummy_df[[col]] <- factor(dummy_df[[col]], levels = levels(df_missing[[col]]))
          }
        }
        return(dummy_df)
      })
    },
    
    MICE = {
      start_time <- Sys.time()
      message("7: Avant MICE imputation")
      View(df_missing_impute)
      
      # Vérifier que toutes les colonnes ont une méthode définie
      if(any(mice_methods == "")) {
        stop("Erreur MICE: Certaines colonnes n'ont pas de méthode d'imputation")
      }
  
      # Test rapide pour voir si mice() fonctionne avant d'appeler complete()
      mice_test <- tryCatch({
        mice(df_missing_impute, m = 1, maxit = 5,
             method = mice_methods,
             predictorMatrix = predMatrix,
             ridge = 0.01,
             printFlag = FALSE)
      }, error = function(e) {
        message("Erreur MICE: ", conditionMessage(e))
        NULL
      })
      
      if(is.null(mice_test)) {
        stop("MICE a échoué avant l'imputation")
      }
      
      message("7: MICE a réussi, maintenant complete()")
      
      # Exécuter l'imputation complète
      imputations_out <- mice::complete(mice_test)
      
      message("7: Après MICE imputation")
      end_time <- Sys.time()
      list(imputed_data = imputations_out, time = as.numeric(end_time - start_time, units = "secs"))
    },
    
    
    MissForest = {
      start_time <- Sys.time()
      message("8: MissForest start")
      
      # # --- Debug: Print column types before any conversion ---
      # message("\nColumn types BEFORE conversion:")
      # print(sapply(df_missing, class))
      
      # Convert characters to factors
      df_missing <- df_missing %>%
        mutate(across(where(is.character), as.factor))
      
      # # --- Debug: Print column types after conversion ---
      # message("\nColumn types AFTER character-to-factor conversion:")
      # print(sapply(df_missing, class))
      
      # Identify invalid columns (e.g., logical, Date, POSIXct)
      invalid_cols <- sapply(df_missing, function(x) {
        !(is.numeric(x) | is.factor(x))
      })
      
      if (any(invalid_cols)) {
        stop(
          "MissForest requires numeric/factor columns. Problematic columns:\n",
          paste(names(df_missing)[invalid_cols], " (type: ", sapply(df_missing[invalid_cols], class), ")", collapse = "\n")
        )
      }
      
      # --- Debug: Force data.frame format (not tibble) ---
      df_missing <- as.data.frame(df_missing)
      
      # Run missForest
      mf <- missForest::missForest(df_missing)
      mf_imp <- mf$ximp
      
      # Reconstruct factors with original levels
      factor_cols <- sapply(df_missing, is.factor)
      mf_imp[factor_cols] <- lapply(names(df_missing)[factor_cols], function(col) {
        factor(
          as.character(mf_imp[[col]]),
          levels = levels(df_missing[[col]])
        )
      })
      
      message("8: MissForest end")
      end_time <- Sys.time()
      list(imputed_data = mf_imp, time = as.numeric(end_time - start_time, units = "secs"))
    },
    
    FAMD = {
      start_time <- Sys.time()
      message("9: FAMD start")
      #View(df_missing)
      # message("FAMD Debug - Data Structure:")
      # str(df_missing)
      # 
      # message("\nColumn Types Summary:")
      # print(sapply(df_missing, class))
      # 
      # message("\nMissing Value Count per Column:")
      # print(colSums(is.na(df_missing)))
      
      # message("\nFirst 6 Rows:")
      # print(head(df_missing))
      
      
      ncp <- missMDA::estim_ncpFAMD(df_missing, ncp.max = 5)$ncp
      
      # message("9: FAMD after ncp choice")
      
      # Effectuer l'imputation FAMD et stocker le résultat
      famd_result <- missMDA::imputeFAMD(df_missing, ncp = max(ncp, 1))$completeObs
      
      # message("9: FAMD midway")
      
      # Forcer la correspondance des colonnes avec le jeu de données original
      famd_result <- famd_result[, colnames(df_missing)]
      
      message("9: FAMD end")
      end_time <- Sys.time()
      list(imputed_data = famd_result, time = as.numeric(end_time - start_time, units = "secs"))
    },
    
    
    KProto = {
      start_time <- Sys.time()
      message("10, kProto start")
      tryCatch({
        imp_df <- kproto_kPOD(df_missing, k = 5, na.rm = "imp_onestep")$data
        
        # Ensure we have correct dimensions
        if (!identical(dim(imp_df), dim(df_missing))) {
          stop("KProto returned incorrect dimensions")
        }
        
        message("10, kProto end")
        end_time <- Sys.time()
        list(imputed_data = imp_df, time = as.numeric(end_time - start_time, units = "secs"))
      }, error = function(e) {
        message("KProto failed: ", conditionMessage(e))
        
        # Create dummy dataframe with correct dimensions
        dummy_df <- as.data.frame(matrix(NA, 
                                         nrow = nrow(df_missing), 
                                         ncol = ncol(df_missing)))
        colnames(dummy_df) <- colnames(df_missing)
        
        # Preserve factor types and levels
        for (col in colnames(df_missing)) {
          if (is.factor(df_missing[[col]])) {
            dummy_df[[col]] <- factor(NA, 
                                      levels = levels(df_missing[[col]]))
          } else if (is.numeric(df_missing[[col]])) {
            dummy_df[[col]] <- as.numeric(dummy_df[[col]])
          }
        }
        return(dummy_df)
      })
    })

  View(imputations)
  return(imputations)
}

evaluate <- function(orig, imputed, mask) {
  # Ensure inputs are data frames and mask is logical
  stopifnot(is.data.frame(orig), is.data.frame(imputed), is.logical(mask))
  
  # Identify numeric and factor columns
  numeric_cols <- sapply(orig, is.numeric)
  factor_cols <- sapply(orig, is.factor)
  
  # Initialize metrics
  metrics <- list()
  
  # Calculate NRMSE for numeric columns
  if (any(numeric_cols)) {
    orig_numeric <- orig[, numeric_cols, drop = FALSE]
    imputed_numeric <- imputed[, numeric_cols, drop = FALSE]
    numeric_mask <- mask[, numeric_cols, drop = FALSE]
    
    rmse <- sqrt(mean((orig_numeric[numeric_mask] - imputed_numeric[numeric_mask])^2))
    data_range <- diff(range(orig_numeric, na.rm = TRUE))
    
    metrics$NRMSE <- ifelse(data_range != 0, rmse/data_range, NA_real_)
  }
  
  # Calculate accuracy for factor columns
  if (any(factor_cols)) {
    orig_factor <- orig[, factor_cols, drop = FALSE]
    imputed_factor <- imputed[, factor_cols, drop = FALSE]
    factor_mask <- mask[, factor_cols, drop = FALSE]
    
    matches <- orig_factor[factor_mask] == imputed_factor[factor_mask]
    metrics$Accuracy <- mean(matches, na.rm = TRUE)
  }
  
  return(metrics)
}

# ---------------------------
# 3. Run Full Pipeline
# ---------------------------
run_pipeline <- function(data_loader, props = c(0.05, 0.10, 0.25, 0.50)) {
  # Load data with validation
  df <- tryCatch({
    data <- data_loader()
    validate_data(data)
    data
  }, error = function(e) {
    message("Data loading failed: ", conditionMessage(e))
    return(NULL)
  })
  
  if (is.null(df)) return(NULL)
  
  results <- list()
  original_dims <- dim(df)
  
  for (prop in props) {
    cat("\n=== Processing missingness:", scales::percent(prop), "===\n")
    
    # Introduce missing values with dimension protection
    missing_data <- tryCatch({
      md <- add_missing(df, prop)
      validate_missing_data(md$data, original_dims)
      md
    }, error = function(e) {
      message("Missing generation failed: ", conditionMessage(e))
      return(NULL)
    })
    # View(missing_data$data) 
    if (is.null(missing_data)) next
    
    # Run imputations with critical checks
    imputations <- tryCatch({
      message("Avant imps")
      Avant_imps <- missing_data$data
      # View(Avant_imps)
      imps <- impute_data(missing_data$data)
      message("Après imps")
      # View(imps)
      validate_imputations(imps, missing_data$data)
      imps
      
    }, error = function(e) {
      message("Erreur en dessous de imps")
      message("Imputation failed: ", conditionMessage(e))
      return(NULL)
    })
    
    # Evaluate and store results
    results[[as.character(prop)]] <- process_imputations(
      imputations, df, missing_data$mask
    )
  }
  
  # Return compiled results or NULL
  compile_results(results)
}

# Helper functions --------------------------------------------------------

# Shared validation function
validate_dataset <- function(df) {
  if (!any(sapply(df, is.numeric))) stop("No numeric columns")
  if (!any(sapply(df, is.factor))) stop("No factor columns")
  if (nrow(df) < 50) warning("Small sample (", nrow(df), " rows)")
}

validate_data <- function(data) {
  if (!any(sapply(data, is.numeric))) stop("No numeric columns")
  if (!any(sapply(data, is.factor))) stop("No factor columns")
  if (ncol(data) < 2) stop("Need ≥2 columns")
}

validate_missing_data <- function(md, original_dims) {
  if (!identical(dim(md), original_dims)) {
    stop(sprintf("Dimension mismatch: %dx%d vs original %dx%d",
                 nrow(md), ncol(md), original_dims[1], original_dims[2]))
  }
  if (any(colSums(is.na(md)) == nrow(md))) {
    stop("Columns with 100% NAs detected")
  }
}

validate_imputations <- function(imps, original) {
  lapply(names(imps), function(method_name) {
    imp <- imps[[method_name]]$imputed_data 
    if (!identical(dim(imp), dim(original))) {
      # Print detailed debug info
      message("\n=== Dimension Mismatch in ", method_name, " ===")
      message("Original dimensions: ", nrow(original), " rows x ", ncol(original), " cols")
      message("Imputed dimensions: ", nrow(imp), " rows x ", ncol(imp), " cols")
      
      # Column differences
      missing_cols <- setdiff(colnames(original), colnames(imp))
      extra_cols <- setdiff(colnames(imp), colnames(original))
      
      if (length(missing_cols) > 0) 
        message("Missing columns: ", paste(missing_cols, collapse = ", "))
      if (length(extra_cols) > 0) 
        message("Extra columns: ", paste(extra_cols, collapse = ", "))
      
      stop("Dimension validation failed for ", method_name)
    }
  })
}

process_imputations <- function(imps, original, mask) {
  lapply(imps, function(imp) {
    tryCatch({
      # Use imp$imputed_data for evaluation
      metrics <- evaluate(original, imp$imputed_data, mask)
      # add the runtime 
      metrics$Time <- imp$time
      metrics
    }, error = function(e) {
      message("Evaluation failed: ", conditionMessage(e))
      list(NRMSE = NA_real_, Accuracy = NA_real_, Time = NA_real_)
    })
  })
}

compile_results <- function(results) {
  if (length(results) == 0) return(NULL)
  
  # Create list to store processed results
  compiled <- list()
  
  for (prop in names(results)) {
    prop_results <- results[[prop]]
    
    # Convert each missingness proportion's results to a data frame with proper structure
    df_prop <- bind_rows(
      lapply(names(prop_results), function(method) {
        data.frame(
          Method = method,
          NRMSE = if (!is.null(prop_results[[method]]$NRMSE)) prop_results[[method]]$NRMSE else NA_real_,
          Accuracy = if (!is.null(prop_results[[method]]$Accuracy)) prop_results[[method]]$Accuracy else NA_real_,
          Time = if (!is.null(prop_results[[method]]$Time)) prop_results[[method]]$Time else NA_real_
        )
      })
    )
    df_prop$Proportion <- as.numeric(prop)
    compiled[[prop]] <- df_prop
  }
  
  bind_rows(compiled) %>%
    mutate(
      Proportion = as.numeric(Proportion),
      Method = factor(Method, levels = c("KNN", "MICE", "MissForest", "FAMD", "KProto"))
    )
}


# ---------------------------
# 5. Dataset-specific wrapper 
# ---------------------------

# Specific for WDBC (Wisconsin Diagnostic Breast Cancer) dataset
#***************************************************************
load_breast_cancer <- function(n_rows = 50) {
  # Define vector of 32 columns names
  column_names <- c(
    "ID", "Diagnosis",
    # Define 10 column names with prefix "Mean_"   
    paste0("Mean_", c("Radius", "Texture", "Perimeter", "Area", 
                      "Smoothness", "Compactness", "Concavity", 
                      "Concave_Points", "Symmetry", "Fractal_Dimension")),
    # Define 10 column names with prefix "SE_"  for Standard Error  
    paste0("SE_", c("Radius", "Texture", "Perimeter", "Area", 
                    "Smoothness", "Compactness", "Concavity", 
                    "Concave_Points", "Symmetry", "Fractal_Dimension")),
    # Define 10 column names with prefix "worst_"
    paste0("Worst_", c("Radius", "Texture", "Perimeter", "Area", 
                       "Smoothness", "Compactness", "Concavity", 
                       "Concave_Points", "Symmetry", "Fractal_Dimension"))
  )
  
  
  # Dataset URL
  url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
  
  # Load data with column_names header (there are no header in loaded csv file) and put the dataframe in the variable df
  df <- read.csv(url, header = FALSE, col.names = column_names, nrows = n_rows)
  
  
  df_modif <- df %>%
    
    select(-ID) %>%
    mutate(
      Diagnosis = factor(
        Diagnosis, 
        levels = c("B", "M"),
        labels = c("Benign", "Malignant")
      )
    )
  
  
  return(df_modif)
}


# Specific for UCI credit cards data
#***************************************************************
load_credit_default <- function(n_rows = 1000) {
  url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/00350/default%20of%20credit%20card%20clients.xls"
  
  tryCatch({
    # Télécharger le fichier Excel en mode binaire dans un fichier temporaire
    temp_file <- tempfile(fileext = ".xls")
    download.file(url, temp_file, mode = "wb", quiet = TRUE)
    
    # Lecture du fichier Excel depuis la feuille "Data", en sautant la première ligne
    df_credit_modif <- readxl::read_excel(temp_file, sheet = "Data", skip = 1, n_max = n_rows) %>%
      rename(default = "default payment next month") %>%
      select(-ID) %>%
      mutate(
        default = factor(default, levels = c(0, 1), labels = c("No", "Yes")),
        SEX = factor(SEX, levels = c(1, 2), labels = c("Male", "Female")),
        EDUCATION = factor(EDUCATION,
                           levels = c(0, 1, 2, 3, 4, 5, 6),
                           labels = c("Unknown", "Graduate", "University", "High School", "Other", "Other", "Other")),
        MARRIAGE = factor(MARRIAGE,
                          levels = c(0, 1, 2, 3),
                          labels = c("Unknown", "Married", "Single", "Other"))
      ) %>%
      mutate(across(matches("^PAY_[0-9]$"), 
                    ~ factor(.,
                             levels = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                             labels = c("No Consumption", "Paid in Full", "On Time",
                                        "Delay 1M", "Delay 2M", "Delay 3M", "Delay 4M",
                                        "Delay 5M", "Delay 6M", "Delay 7M", "Delay 8M",
                                        "Delay 9M")))) %>%
      mutate(across(c(LIMIT_BAL, AGE, BILL_AMT1:BILL_AMT6, PAY_AMT1:PAY_AMT6), as.numeric))
    
    unlink(temp_file)
    View(df_credit_modif)
    return(df_credit_modif)
  }, error = function(e) {
    message("Error: ", conditionMessage(e))
    if(exists("temp_file")) unlink(temp_file)
    return(NULL)
  })
}

# Specific for UCI retail data
#***************************************************************
load_online_retail <- function(n_rows = 2000) {
  #Install package readxl if not already available
  if (!requireNamespace("readxl", quietly = TRUE)) {
    install.packages("readxl")
  }
  
  url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/00352/Online%20Retail.xlsx"
  temp_file <- tempfile(fileext = ".xlsx")
  download.file(url, temp_file, mode = "wb", quiet = TRUE)
  
  df_retail <- readxl::read_excel(temp_file, col_names = TRUE, n_max = n_rows)
  
  View(df_retail)
  
  df_retail_modif <- readxl::read_excel(temp_file, col_names = TRUE, n_max = n_rows) %>% 
    # Now remove unwanted columns - KEEP ALL OTHER COLUMNS
    select(-InvoiceDate, -InvoiceNo, -StockCode, -CustomerID, -Description) %>%
    # Convert and rename remaining columns with safeguards
    mutate(
      country = factor(Country, levels = unique(Country)),  # Preserve all levels
      quantity = as.numeric(Quantity),
      unit_price = as.numeric(UnitPrice)
    ) %>%
    select(-Country) %>%
    # Remove columns that become all NA
    select(where(~!all(is.na(.x))))  # Keep only columns with some non-NA values
  
  
  # Suppress of temp file in the path temp-file 
  unlink(temp_file)
  
  # ******************************************
  
  View(df_retail_modif)
  
  # validation
  if(ncol(df_retail_modif) < 3) stop("Too many columns dropped during processing")
  if(!any(sapply(df_retail_modif, is.numeric))) stop("No numeric columns")
  if(!any(sapply(df_retail_modif, is.factor))) stop("No factor columns")
  
  return(df_retail_modif)
}

# Specific for UCI obesity data
#***************************************************************
load_obesity_data <- function(n_rows = 1000) {
  # Try local file first, then online
  online_url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/00544/ObesityDataSet_raw_and_data_sinthetic.csv"
  local_file <- "ObesityDataSet_raw_and_data_sinthetic.csv"
  
  tryCatch({
    if (file.exists(local_file)) {
      df_obesity <- read.csv(local_file, stringsAsFactors = FALSE, nrows = n_rows)
    } else {
      df_obesity <- read.csv(online_url, stringsAsFactors = FALSE, nrows = n_rows)
    }
    
    
    df_obesity_modif <- df_obesity %>%
      mutate(
        across(any_of(c("Gender", "family_history_with_overweight", "FAVC", "CAEC", 
                        "SMOKE", "SCC", "CALC", "MTRANS", "NObeyesdad")), as.factor),
        across(any_of(c("Age", "Height", "Weight", "FCVC", "NCP", "CH2O", "FAF", "TUE")), as.numeric),
        obesity_level = factor(NObeyesdad) 
      ) %>%
      select(-NObeyesdad) 
    
    # Validation
    validate_dataset(df_obesity_modif)
    
    message("Obestity data number of levels : ")
    print(sapply(df_obesity_modif, function(x) if (is.factor(x)) nlevels(x) else NA))

    return(df_obesity_modif)
  }, error = function(e) {
    message("Obesity load failed: ", conditionMessage(e))
    return(NULL)
  })
}



# Specific for UCI Steel industry  data
#***************************************************************
load_steel_energy <- function(n_rows = 4500) {
  # URL alternatives
  online_url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/00616/Steel_industry_data.csv"
  local_file <- "Steel_industry_data.csv"
  
  tryCatch({
    if (!requireNamespace("janitor", quietly = TRUE)) install.packages("janitor")
    if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr") 
    
    if (file.exists(local_file)) {
      df_steel <- read.csv(local_file, stringsAsFactors = FALSE, nrows = n_rows) 
    } else {
      df_steel <- read.csv(online_url, stringsAsFactors = FALSE, nrows = n_rows)
    }
    
    View(df_steel)
    
    df_steel_modif <- df_steel %>%
      { 
        # Conversion des noms de colonnes en minuscules sans ajouter d'underscores supplémentaires
        names(.) <- tolower(names(.))  
        .
      } %>%
      select(-any_of("date")) %>%
      mutate(
        day_of_week = factor(day_of_week, 
                             levels = c("Monday","Tuesday","Wednesday",
                                        "Thursday","Friday","Saturday","Sunday"),
                             ordered = FALSE),
        weekstatus = factor(weekstatus,levels = c("Weekday", "Weekend"), ordered = FALSE),
        load_type = factor(load_type, levels = c("Light_Load", "Maximum_Load", "Medium_Load"), ordered = FALSE),
        across(c(usage_kwh:nsm), as.numeric)
      )
    
    
    validate_dataset(df_steel_modif)
    
    View(df_steel_modif)
    
    return(df_steel_modif)
  }, error = function(e) {
    message("Steel data error: ", conditionMessage(e))
    return(NULL)
  })
}

# Specific for UCI german Credit  data
#***************************************************************
load_german_credit <- function(n_rows = Inf) {
  online_url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/statlog/german/german.data"
  
  column_names <- c(
    "status", "duration", "credit_history", "purpose", "amount",
    "savings", "employment", "installment_rate", "personal_status",
    "other_debtors", "residence", "property", "age", "other_installments",
    "housing", "existing_credits", "job", "dependents", "telephone", 
    "foreign_worker", "credit_risk"
  )
  
  # Lecture du fichier initial (df_german)
  df_german <- read.csv(online_url, stringsAsFactors = FALSE, nrows = n_rows)
  
  # Lecture du fichier avec les colonnes et conversion (df_german_modif)
  df_german_modif <- read.csv(online_url, header = FALSE, sep = " ", 
                              col.names = column_names, nrows = n_rows) %>%  
    mutate(
      # Conversion des variables catégorielles avec niveaux et labels définis
      status = factor(status, 
                      levels = c("A11", "A12", "A13", "A14"),
                      labels = c("<0 DM", "0-200 DM", ">200 DM", "none")),
      credit_history = factor(credit_history,
                              levels = c("A30", "A31", "A32", "A33", "A34"),
                              labels = c("no credits", "all paid", "existing paid", 
                                         "delay previously", "critical accounts")),
      purpose = factor(purpose,
                       levels = paste0("A4", 0:9),
                       labels = c("car (new)", "car (used)", "furniture", 
                                  "radio/TV", "appliances", "repairs", 
                                  "education", "vacation", "retraining", 
                                  "business")),
      savings = factor(savings,
                       levels = c("A61", "A62", "A63", "A64", "A65"),
                       labels = c("<100 DM", "100-500 DM", "500-1000 DM",
                                  ">=1000 DM", "unknown/no account")),
      employment = factor(employment,
                          levels = c("A71", "A72", "A73", "A74", "A75"),
                          labels = c("unemployed", "<1 year", "1-4 years", 
                                     "4-7 years", ">=7 years")),
      personal_status = factor(personal_status,
                               levels = c("A91", "A92", "A93", "A94"),
                               labels = c("male-divorced", "female-divorced",
                                          "male-single", "male-married")),
      other_debtors = factor(other_debtors,
                             levels = c("A101", "A102", "A103"),
                             labels = c("none", "co-applicant", "guarantor")),
      credit_risk = factor(credit_risk,
                           levels = c(1, 2),
                           labels = c("Good", "Bad"))
    ) %>%
    # Conversion des variables numériques
    mutate(across(c(duration, amount, installment_rate, residence, age, 
                    existing_credits, dependents), as.numeric)) %>%
    # Conversion des variables qui doivent être catégorielles (initialement en character)
    mutate(
      property = as.factor(property),
      other_installments = as.factor(other_installments),
      housing = as.factor(housing),
      job = as.factor(job),
      telephone = as.factor(telephone),
      foreign_worker = as.factor(foreign_worker)
    )
  
  # Vérifications finales
  if(!any(sapply(df_german_modif, is.numeric))) stop("No numeric columns")
  if(!any(sapply(df_german_modif, is.factor))) stop("No factor columns")
  
  # message("German Credit Debug - Data Structure:")
  # str(df_german_modif)
  # 
  # message("\nColumn Types Summary:")
  # print(sapply(df_german_modif, class))
  # 
  # message("\nFirst 6 Rows:")
  # print(head(df_german_modif))

  message("German Credit number of levels : ")
  print(sapply(df_german_modif, function(x) if (is.factor(x)) nlevels(x) else NA))

  return(df_german_modif)
}



# --------------------------
# 6. execution
# --------------------------

# Define dataset list
datasets <- list(
  #"Breast Cancer" = load_breast_cancer#,
  #"Credit Default" = load_credit_default#,
  #"Online Retail" = load_online_retail#,
  #"Obesity" = load_obesity_data#,
  "Steel Energy" = load_steel_energy#,
  #"German Credit" = load_german_credit
)

# Modified all_results generation with robust error handling
all_results <- imap_dfr(datasets, ~{
  tryCatch({
    # Print name of dataset in console
    message("Processing dataset: ", .y)
    
    # Run pipeline and add dataset identifier
    results <- run_pipeline(.x)
    
    # Validate results before processing
    if (is.null(results) || nrow(results) == 0) {
      message("No valid results for ", .y)
      return(NULL)
    }
    
    # Add dataset column with type safety
    mutate(results, Dataset = factor(.y), .before = 1)
    
  }, error = function(e) {
    # Enhanced error logging
    error_msg <- paste("Skipped", .y, ":", conditionMessage(e))
    message(error_msg)
    write(paste(Sys.time(), error_msg), "errors.log", append = TRUE)
    NULL
  })
}) %>% 
  # Ensure proper typing of results columns
  mutate(
    across(c(NRMSE, Accuracy), as.numeric),
    Dataset = factor(Dataset),
    Method = factor(Method, levels = c("KNN", "MICE", "MissForest", "FAMD", "KProto"))
  )

# Calculate average metrics across datasets
avg_results <- all_results %>%
  group_by(Proportion, Method) %>%
  summarise(
    Avg_NRMSE = mean(NRMSE, na.rm = TRUE),
    Avg_Accuracy = mean(Accuracy, na.rm = TRUE),
    Avg_Time = mean(Time, na.rm = TRUE),  
    SD_NRMSE = sd(NRMSE, na.rm = TRUE),
    SD_Accuracy = sd(Accuracy, na.rm = TRUE),
    SD_Time = sd(Time, na.rm = TRUE),     
    N_Datasets = n_distinct(Dataset),
    .groups = "drop"
  ) %>%
  filter(N_Datasets >= 1) %>%
  mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA_real_, .)))

# Save both raw and averaged results
write_csv(all_results, "full_results.csv")
write_csv(avg_results, "average_results.csv")

# Structure of all_results:
# | Proportion | Method     | NRMSE | Accuracy | Dataset       |
# |------------|------------|-------|----------|---------------|
# | 0.05       | KNN        | 0.012 | 0.92     | Breast Cancer |
# | ...        | ...        | ...   | ...      | ...           |

# Generate comparison plots
# Generate individual plots for each dataset
plot_results <- function(results) {
  ggplot(
    data = results %>%
      pivot_longer(
        cols = c(NRMSE, Accuracy),
        names_to = "Metric",
        values_to = "Value"
      ),
    mapping = aes(
      x = Proportion,
      y = Value,
      color = Method,
      group = Method
    )
  ) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    facet_wrap(~Metric, scales = "free_y", ncol = 1) +
    scale_x_continuous(
      labels = scales::percent_format(accuracy = 1),
      breaks = c(0.05, 0.10, 0.25, 0.50)
    ) +
    labs(
      x = "Missing Data Proportion",
      y = "Metric Value",
      title = "Imputation Method Comparison",
      color = "Imputation Method"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      panel.spacing = unit(1, "lines"),
      strip.text = element_text(face = "bold")
    )
}

# Time comparison per dataset
plot_time_results <- function(results) {
  ggplot(results, aes(x = Proportion, y = Time, color = Method)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_x_continuous(labels = scales::percent) +
    labs(
      title = "Imputation Method Runtimes",
      x = "Missing Data Proportion",
      y = "Time (seconds)",
      color = "Method"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Average time across datasets
plot_avg_time <- function(avg_results) {
  ggplot(avg_results, aes(x = Proportion, y = Avg_Time, color = Method)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_x_continuous(labels = scales::percent) +
    labs(
      title = "Average Runtime Across Datasets",
      x = "Missing Data Proportion",
      y = "Time (seconds)"
    ) +
    theme_minimal()
}

# Generate Individual Plots per Dataset
walk(unique(all_results$Dataset), function(dataset) {
  dataset_results <- filter(all_results, Dataset == dataset)
  
  p <- plot_results(dataset_results) + 
    ggtitle(paste("Imputation Comparison -", dataset))
  
  ggsave(paste0("plots/", gsub(" ", "_", dataset), ".png"), 
         plot = p, width = 8, height = 6)
})

# 1. Combined faceted plot
combined_plot <- avg_results %>%
  pivot_longer(
    cols = c(Avg_NRMSE, Avg_Accuracy),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(Metric = gsub("Avg_", "", Metric)) %>%
  ggplot(aes(x = Proportion, y = Value, color = Method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  facet_wrap(~Metric, scales = "free_y", ncol = 1) +
  scale_x_continuous(labels = scales::percent) +
  labs(title = "Average Performance Across Datasets",
       x = "Missing Data Proportion",
       y = "Metric Value") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# 2. Separate plots
nrmse_plot <- ggplot(avg_results, aes(x = Proportion, y = Avg_NRMSE, color = Method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_x_continuous(labels = scales::percent) +
  labs(title = "Average Normalized RMSE",
       y = "NRMSE") +
  theme_minimal(base_size = 14)

accuracy_plot <- ggplot(avg_results, aes(x = Proportion, y = Avg_Accuracy, color = Method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_x_continuous(labels = scales::percent) +
  labs(title = "Average Classification Accuracy",
       y = "Accuracy") +
  theme_minimal(base_size = 14)

# Save all versions
ggsave("average_combined.png", combined_plot, width = 8, height = 10)
ggsave("average_nrmse.png", nrmse_plot, width = 8, height = 5)
ggsave("average_accuracy.png", accuracy_plot, width = 8, height = 5)

# Combined side-by-side view
side_by_side <- nrmse_plot + accuracy_plot + plot_layout(guides = "collect")
ggsave("average_side_by_side.png", side_by_side, width = 14, height = 6)


walk(unique(all_results$Dataset), function(dataset) {
  dataset_results <- filter(all_results, Dataset == dataset)
  
  p_time <- plot_time_results(dataset_results) + 
    ggtitle(paste("Runtime Comparison -", dataset))
  
  ggsave(paste0("plots/time_", gsub(" ", "_", dataset), ".png"), 
         plot = p_time, width = 8, height = 6)
})

# Generate average time plot
avg_time_plot <- plot_avg_time(avg_results)
ggsave("average_runtime.png", avg_time_plot, width = 8, height = 6)

# Faceted time comparison
faceted_time_plot <- ggplot(all_results, aes(x = Proportion, y = Time, color = Method)) +
  geom_line() +
  geom_point() +
  facet_wrap(~Dataset, scales = "free_y") +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Missing Proportion", y = "Time (seconds)") +
  theme_bw()

ggsave("faceted_runtimes.png", faceted_time_plot, width = 12, height = 8)
