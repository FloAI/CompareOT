# ==============================================================================
# Multi-Omics Imputation & Wasserstein Distance Comparison
# Scenarios: Complete (CC), Incomplete (IncC), Imputed (ImpC)
# Methods: OTrecod, kNN, MICE, missForest (MF), missMDA (MDA)
# ==============================================================================

suppressPackageStartupMessages({
  library(MASS)
  library(VIM)
  library(mice)
  library(OTrecod)
  library(missForest)
  library(missMDA)
  library(dplyr)
})

set.seed(42)

# ------------------------
# 1. Load Raw Datasets
# ------------------------
data(ncds_14)
data(ncds_5)

# ------------------------
# 2. Evaluation Helper (Full Column Wasserstein)
# FIX 1: sort levels numerically, not lexicographically
# ------------------------
discrete_wass <- function(vec1, vec2) {
  vec1 <- trimws(as.character(vec1))
  vec2 <- trimws(as.character(vec2))

  tab1 <- table(vec1) / length(vec1)
  tab2 <- table(vec2) / length(vec2)

  # Numeric sort prevents "10" < "2" corrupting the CDF
  all_levels <- sort(unique(c(names(tab1), names(tab2))), method = "radix")

  p1 <- tab1[all_levels]; p1[is.na(p1)] <- 0
  p2 <- tab2[all_levels]; p2[is.na(p2)] <- 0

  sum(abs(cumsum(p1) - cumsum(p2)))
}

# ------------------------
# 3. Core Simulation Function
# ------------------------
run_scenario <- function(scenario_name) {
  cat(sprintf("\n--- Running Scenario: %s ---\n", scenario_name))

  if (scenario_name == "CC") {
    db1 <- ncds_14[complete.cases(ncds_14), ]
    db2 <- ncds_5[complete.cases(ncds_5), ]
    merged_tab <- merge_dbs(db1, db2, row_ID1 = 1, row_ID2 = 1, NAME_Y = "GO90", NAME_Z = "RG91", ordinal_DB1 = 3, ordinal_DB2 = 4, seed_choice = 3023)
  } else if (scenario_name == "IncC") {
    merged_tab <- merge_dbs(ncds_14, ncds_5, row_ID1 = 1, row_ID2 = 1, NAME_Y = "GO90", NAME_Z = "RG91", ordinal_DB1 = 3, ordinal_DB2 = 4, seed_choice = 3023)
  } else if (scenario_name == "ImpC") {
    merged_tab <- merge_dbs(ncds_14, ncds_5, row_ID1 = 1, row_ID2 = 1, NAME_Y = "GO90", NAME_Z = "RG91", ordinal_DB1 = 3, ordinal_DB2 = 4, impute = "MICE", R_MICE = 2, seed_choice = 3023)
  }

  df_base <- merged_tab$DB_READY[, -4]

  n_Y <- sum(is.na(df_base$Y))
  n_Z <- sum(is.na(df_base$Z))

  # --- Method 1: OTrecod ---
  df_ot <- df_base
  outj1 <- try(OT_joint(df_ot, nominal = c(1:4), ordinal = 5:6, dist.choice = "E", which.DB = "BOTH"), silent = TRUE)
  if (!inherits(outj1, "try-error")) {
    df_ot$Y[is.na(df_ot$Y)] <- outj1$DATA2_OT[, "OTpred"]
    df_ot$Z[is.na(df_ot$Z)] <- outj1$DATA1_OT[, "OTpred"]
  }
  df_ot <- df_ot[order(as.numeric(rownames(df_ot))), ]

  # --- Method 2: kNN ---
  capture.output(df_knn <- try(kNN(df_base, k = 5, imp_var = FALSE), silent = TRUE))
  if (inherits(df_knn, "try-error")) df_knn <- df_base
  df_knn <- df_knn[order(as.numeric(rownames(df_knn))), ]

  # --- Method 3: MICE ---
  capture.output({
    init <- mice(df_base, maxit = 0)
    meth <- init$method
    predM <- init$predictorMatrix
    meth[c("Y", "Z")] <- "rf"
    imputed_mice <- try(mice(df_base, method = meth, predictorMatrix = predM, m = 1, maxit = 5, printFlag = FALSE), silent = TRUE)
    if (!inherits(imputed_mice, "try-error")) {
      df_mice <- complete(imputed_mice)
    } else {
      df_mice <- df_base
    }
  })
  df_mice <- df_mice[order(as.numeric(rownames(df_mice))), ]

  # --- Method 4: missForest ---
  df_impute_rf <- df_base %>% mutate_if(is.character, as.factor)
  capture.output(mf_res <- try(missForest(df_impute_rf, maxiter = 5, ntree = 50), silent = TRUE))
  df_mf <- if (!inherits(mf_res, "try-error")) mf_res$ximp else df_base
  # FIX 3: reorder rows to match other methods
  df_mf <- df_mf[order(as.numeric(rownames(df_mf))), ]

  # --- Method 5: missMDA ---
  capture.output(famd_res <- try(imputeFAMD(df_impute_rf, ncp = 2), silent = TRUE))
  if (!inherits(famd_res, "try-error")) {
    df_mda <- famd_res$completeObs
    # FIX 2: safe factor -> integer -> character, avoids regex fragility
    df_mda$Y <- as.character(as.integer(df_mda$Y))
    df_mda$Z <- as.character(as.integer(df_mda$Z))
  } else {
    df_mda <- df_base
  }
  # FIX 3: reorder rows to match other methods
  df_mda <- df_mda[order(as.numeric(rownames(df_mda))), ]

  # --- Calculate Wasserstein Distances (Full Column) ---
  res <- data.frame(
    Scenario = scenario_name,
    OT_kNN_Y = discrete_wass(df_ot$Y, df_knn$Y),
    OT_kNN_Z = discrete_wass(df_ot$Z, df_knn$Z),
    OT_MICE_Y = discrete_wass(df_ot$Y, df_mice$Y),
    OT_MICE_Z = discrete_wass(df_ot$Z, df_mice$Z),
    OT_MF_Y = discrete_wass(df_ot$Y, df_mf$Y),
    OT_MF_Z = discrete_wass(df_ot$Z, df_mf$Z),
    OT_MDA_Y = discrete_wass(df_ot$Y, df_mda$Y),
    OT_MDA_Z = discrete_wass(df_ot$Z, df_mda$Z),

    kNN_MICE_Y = discrete_wass(df_knn$Y, df_mice$Y),
    kNN_MICE_Z = discrete_wass(df_knn$Z, df_mice$Z),
    kNN_MF_Y = discrete_wass(df_knn$Y, df_mf$Y),
    kNN_MF_Z = discrete_wass(df_knn$Z, df_mf$Z),
    kNN_MDA_Y = discrete_wass(df_knn$Y, df_mda$Y),
    kNN_MDA_Z = discrete_wass(df_knn$Z, df_mda$Z),

    MICE_MF_Y = discrete_wass(df_mice$Y, df_mf$Y),
    MICE_MF_Z = discrete_wass(df_mice$Z, df_mf$Z),
    MICE_MDA_Y = discrete_wass(df_mice$Y, df_mda$Y),
    MICE_MDA_Z = discrete_wass(df_mice$Z, df_mda$Z),

    MF_MDA_Y = discrete_wass(df_mf$Y, df_mda$Y),
    MF_MDA_Z = discrete_wass(df_mf$Z, df_mda$Z),
    n_Y = n_Y,
    n_Z = n_Z
  )

  return(res)
}

# ------------------------
# 4. Execute All Scenarios
# ------------------------
scenarios <- c("CC", "IncC", "ImpC")
results_list <- lapply(scenarios, run_scenario)
final_results <- bind_rows(results_list)
rownames(final_results) <- final_results$Scenario

# ------------------------
# 5. Transposed LaTeX Table Output
# ------------------------
comparisons <- c("OT_kNN", "OT_MICE", "OT_MF", "OT_MDA", "kNN_MICE", "kNN_MF", "kNN_MDA", "MICE_MF", "MICE_MDA", "MF_MDA")
comp_labels <- c("OT-kNN", "OT-MICE", "OT-MF", "OT-MDA", "kNN-MICE", "kNN-MF", "kNN-MDA", "MICE-MF", "MICE-MDA", "MF-MDA")

cat("\n\n")
cat("\\begin{table}[h!]\n")
cat("\\caption{Discrete 1D Wasserstein distances among OTrecod, kNN, MICE, missForest (MF), and missMDA (MDA) for $Y_A$ and $Y_B$}\n")
cat("\\centering\n")
cat("\\begin{tabular}{l c c c c c c}\n")
cat("\\toprule\n")
cat(" & \\multicolumn{2}{c}{CC} & \\multicolumn{2}{c}{IncC} & \\multicolumn{2}{c}{ImpC} \\\\\n")
cat("\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}\n")
cat("Comparison & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ \\\\\n")
cat("\\midrule\n")

for (i in 1:length(comparisons)) {
  comp <- comparisons[i]
  label <- comp_labels[i]
  col_Y <- paste0(comp, "_Y")
  col_Z <- paste0(comp, "_Z")
  cat(sprintf("%s & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n",
              label,
              final_results[1, col_Y], final_results[1, col_Z],
              final_results[2, col_Y], final_results[2, col_Z],
              final_results[3, col_Y], final_results[3, col_Z]))
}

cat("\\midrule\n")
cat(sprintf("$n$ missing & %d & %d & %d & %d & %d & %d \\\\\n",
            final_results[1, "n_Y"], final_results[1, "n_Z"],
            final_results[2, "n_Y"], final_results[2, "n_Z"],
            final_results[3, "n_Y"], final_results[3, "n_Z"]))
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\label{tab:wass_dist_transposed}\n")
cat("\\end{table}\n")