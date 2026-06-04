# ==============================================================================
# Multi-Omics Imputation & Concordance Comparison
# Scenarios: Complete, Incomplete, Imputed
# Methods: OTrecod, kNN, MICE, missForest (MF), missMDA (MDA)
# ==============================================================================

suppressPackageStartupMessages({
  library(OTrecod)
  library(VIM)
  library(mice)
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
# 2. Evaluation Helpers
# ------------------------
calc_concordance <- function(v1, v2, mask) {
  if (sum(mask) == 0) return(NA_real_)
  mean(as.character(v1[mask]) == as.character(v2[mask]), na.rm = TRUE) * 100
}

calc_multi_concordance <- function(df_list, var_name, mask) {
  if (sum(mask) == 0) return(NA_real_)
  extract_preds <- lapply(df_list, function(df) as.character(df[[var_name]][mask]))
  pred_matrix <- do.call(cbind, extract_preds)
  all_agree <- apply(pred_matrix, 1, function(row) length(unique(row)) == 1)
  mean(all_agree, na.rm = TRUE) * 100
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

  # FIX: record which original row indices are missing BEFORE any reordering,
  # then after sorting we reconstruct masks aligned to the sorted frame
  orig_order  <- as.numeric(rownames(df_base))
  sorted_order <- sort(orig_order)
  # Boolean masks in sorted-row space
  missing_Y_rows <- orig_order[is.na(df_base$Y)]
  missing_Z_rows <- orig_order[is.na(df_base$Z)]
  mask_Y <- sorted_order %in% missing_Y_rows
  mask_Z <- sorted_order %in% missing_Z_rows

  n_Y <- sum(mask_Y)
  n_Z <- sum(mask_Z)

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
    meth[c("Y", "Z")] <- "rf"
    imputed_mice <- try(mice(df_base, method = meth, m = 1, maxit = 5, printFlag = FALSE), silent = TRUE)
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
  df_mf <- df_mf[order(as.numeric(rownames(df_mf))), ]

  # --- Method 5: missMDA ---
  capture.output(famd_res <- try(imputeFAMD(df_impute_rf, ncp = 2), silent = TRUE))
  if (!inherits(famd_res, "try-error")) {
    df_mda <- famd_res$completeObs
    # FIX: safe factor -> integer -> character instead of fragile regex
    df_mda$Y <- as.character(as.integer(df_mda$Y))
    df_mda$Z <- as.character(as.integer(df_mda$Z))
  } else {
    df_mda <- df_base
  }
  df_mda <- df_mda[order(as.numeric(rownames(df_mda))), ]

  # --- Calculate Concordance (masks now correctly aligned to sorted rows) ---
  list_of_dfs <- list(OT = df_ot, kNN = df_knn, MICE = df_mice, MF = df_mf, MDA = df_mda)

  res <- data.frame(
    Scenario   = scenario_name,
    OT_kNN_Y   = calc_concordance(df_ot$Y, df_knn$Y,   mask_Y),
    OT_kNN_Z   = calc_concordance(df_ot$Z, df_knn$Z,   mask_Z),
    OT_MICE_Y  = calc_concordance(df_ot$Y, df_mice$Y,  mask_Y),
    OT_MICE_Z  = calc_concordance(df_ot$Z, df_mice$Z,  mask_Z),
    OT_MF_Y    = calc_concordance(df_ot$Y, df_mf$Y,    mask_Y),
    OT_MF_Z    = calc_concordance(df_ot$Z, df_mf$Z,    mask_Z),
    OT_MDA_Y   = calc_concordance(df_ot$Y, df_mda$Y,   mask_Y),
    OT_MDA_Z   = calc_concordance(df_ot$Z, df_mda$Z,   mask_Z),
    kNN_MICE_Y = calc_concordance(df_knn$Y, df_mice$Y, mask_Y),
    kNN_MICE_Z = calc_concordance(df_knn$Z, df_mice$Z, mask_Z),
    kNN_MF_Y   = calc_concordance(df_knn$Y, df_mf$Y,   mask_Y),
    kNN_MF_Z   = calc_concordance(df_knn$Z, df_mf$Z,   mask_Z),
    kNN_MDA_Y  = calc_concordance(df_knn$Y, df_mda$Y,  mask_Y),
    kNN_MDA_Z  = calc_concordance(df_knn$Z, df_mda$Z,  mask_Z),
    MICE_MF_Y  = calc_concordance(df_mice$Y, df_mf$Y,  mask_Y),
    MICE_MF_Z  = calc_concordance(df_mice$Z, df_mf$Z,  mask_Z),
    MICE_MDA_Y = calc_concordance(df_mice$Y, df_mda$Y, mask_Y),
    MICE_MDA_Z = calc_concordance(df_mice$Z, df_mda$Z, mask_Z),
    MF_MDA_Y   = calc_concordance(df_mf$Y, df_mda$Y,   mask_Y),
    MF_MDA_Z   = calc_concordance(df_mf$Z, df_mda$Z,   mask_Z),
    ALL_Y      = calc_multi_concordance(list_of_dfs, "Y", mask_Y),
    ALL_Z      = calc_multi_concordance(list_of_dfs, "Z", mask_Z),
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

# ------------------------
# 5. LaTeX Table Output
# ------------------------
cat("\n\n")
cat("\\begin{table}[h!]\n")
cat("\\caption{Concordance percentages among OTrecod, kNN, MICE, missForest (MF), and missMDA (MDA) for $Y_A$ and $Y_B$}\n")
cat("\\centering\n")
cat("\\resizebox{\\textwidth}{!}{\n")
cat("\\begin{tabular}{l *{11}{c c} c c}\n")
cat("\\toprule\n")
cat(" & \\multicolumn{2}{c}{OT-kNN} & \\multicolumn{2}{c}{OT-MICE} & \\multicolumn{2}{c}{OT-MF} & \\multicolumn{2}{c}{OT-MDA} & \\multicolumn{2}{c}{kNN-MICE} & \\multicolumn{2}{c}{kNN-MF} & \\multicolumn{2}{c}{kNN-MDA} & \\multicolumn{2}{c}{MICE-MF} & \\multicolumn{2}{c}{MICE-MDA} & \\multicolumn{2}{c}{MF-MDA} & \\multicolumn{2}{c}{All 5} & \\multicolumn{2}{c}{$n$} \\\\\n")
cat("\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9} \\cmidrule(lr){10-11} \\cmidrule(lr){12-13} \\cmidrule(lr){14-15} \\cmidrule(lr){16-17} \\cmidrule(lr){18-19} \\cmidrule(lr){20-21} \\cmidrule(lr){22-23} \\cmidrule(lr){24-25}\n")
cat("Cases & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ \\\\\n")
cat("\\midrule\n")

for (i in 1:nrow(final_results)) {
  r <- final_results[i, ]
  # CC row will have NA concordances (no missing values); print as "--"
  fmt_val <- function(x) if (is.na(x)) "--" else sprintf("%.2f", x)
  cat(sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %d & %d \\\\\n",
              r$Scenario,
              fmt_val(r$OT_kNN_Y),   fmt_val(r$OT_kNN_Z),
              fmt_val(r$OT_MICE_Y),  fmt_val(r$OT_MICE_Z),
              fmt_val(r$OT_MF_Y),    fmt_val(r$OT_MF_Z),
              fmt_val(r$OT_MDA_Y),   fmt_val(r$OT_MDA_Z),
              fmt_val(r$kNN_MICE_Y), fmt_val(r$kNN_MICE_Z),
              fmt_val(r$kNN_MF_Y),   fmt_val(r$kNN_MF_Z),
              fmt_val(r$kNN_MDA_Y),  fmt_val(r$kNN_MDA_Z),
              fmt_val(r$MICE_MF_Y),  fmt_val(r$MICE_MF_Z),
              fmt_val(r$MICE_MDA_Y), fmt_val(r$MICE_MDA_Z),
              fmt_val(r$MF_MDA_Y),   fmt_val(r$MF_MDA_Z),
              fmt_val(r$ALL_Y),      fmt_val(r$ALL_Z),
              r$n_Y, r$n_Z))
}

cat("\\bottomrule\n")
cat("\\end{tabular}\n")# ==============================================================================
# Multi-Omics Imputation & Concordance Comparison
# Scenarios: Complete, Incomplete, Imputed
# Methods: OTrecod, kNN, MICE, missForest (MF), missMDA (MDA)
# ==============================================================================

suppressPackageStartupMessages({
  library(OTrecod)
  library(VIM)
  library(mice)
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
# 2. Evaluation Helpers
# ------------------------
calc_concordance <- function(v1, v2, mask) {
  if (sum(mask) == 0) return(NA_real_)
  mean(as.character(v1[mask]) == as.character(v2[mask]), na.rm = TRUE) * 100
}

calc_multi_concordance <- function(df_list, var_name, mask) {
  if (sum(mask) == 0) return(NA_real_)
  extract_preds <- lapply(df_list, function(df) as.character(df[[var_name]][mask]))
  pred_matrix <- do.call(cbind, extract_preds)
  all_agree <- apply(pred_matrix, 1, function(row) length(unique(row)) == 1)
  mean(all_agree, na.rm = TRUE) * 100
}

# ------------------------
# 3. Core Simulation Function
# ------------------------
run_scenario <- function(scenario_name) {
  cat(sprintf("\n--- Running Scenario: %s ---\n", scenario_name))

  if (scenario_name == "CC") {
    db1 <- ncds_14[complete.cases(ncds_14), ]
    db2 <- ncds_5[complete.cases(ncds_5), ]
    merged_tab <- merge_dbs(db1, db2, row_ID1 = 1, row_ID2 = 1, NAME_Y = "GO90", NAME_Z = "RG91",
                            ordinal_DB1 = 3, ordinal_DB2 = 4, seed_choice = 3023)
  } else if (scenario_name == "IncC") {
    merged_tab <- merge_dbs(ncds_14, ncds_5, row_ID1 = 1, row_ID2 = 1, NAME_Y = "GO90", NAME_Z = "RG91",
                            ordinal_DB1 = 3, ordinal_DB2 = 4, seed_choice = 3023)
  } else if (scenario_name == "ImpC") {
    merged_tab <- merge_dbs(ncds_14, ncds_5, row_ID1 = 1, row_ID2 = 1, NAME_Y = "GO90", NAME_Z = "RG91",
                            ordinal_DB1 = 3, ordinal_DB2 = 4, impute = "MICE", R_MICE = 2, seed_choice = 3023)
  }

  # Sort df_base FIRST so masks and all imputed frames share the same row order
  df_base <- merged_tab$DB_READY[, -4]
  df_base <- df_base[order(as.numeric(rownames(df_base))), ]

  # Masks defined on sorted df_base — valid for all methods below
  mask_Y <- is.na(df_base$Y)
  mask_Z <- is.na(df_base$Z)
  n_Y    <- sum(mask_Y)
  n_Z    <- sum(mask_Z)

  # --- Method 1: OTrecod ---
  df_ot <- df_base
  outj1 <- try(OT_joint(df_ot, nominal = c(1:4), ordinal = 5:6,
                         dist.choice = "E", which.DB = "BOTH"), silent = TRUE)
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
    meth[c("Y", "Z")] <- "rf"
    imputed_mice <- try(mice(df_base, method = meth, m = 1, maxit = 5, printFlag = FALSE), silent = TRUE)
    df_mice <- if (!inherits(imputed_mice, "try-error")) complete(imputed_mice) else df_base
  })
  df_mice <- df_mice[order(as.numeric(rownames(df_mice))), ]

  # --- Method 4: missForest ---
  df_impute_rf <- df_base %>% mutate_if(is.character, as.factor)
  capture.output(mf_res <- try(missForest(df_impute_rf, maxiter = 5, ntree = 50), silent = TRUE))
  df_mf <- if (!inherits(mf_res, "try-error")) mf_res$ximp else df_base
  df_mf <- df_mf[order(as.numeric(rownames(df_mf))), ]

  # --- Method 5: missMDA ---
  capture.output(famd_res <- try(imputeFAMD(df_impute_rf, ncp = 2), silent = TRUE))
  if (!inherits(famd_res, "try-error")) {
    df_mda <- famd_res$completeObs
    df_mda$Y <- as.character(as.integer(df_mda$Y))
    df_mda$Z <- as.character(as.integer(df_mda$Z))
  } else {
    df_mda <- df_base
  }
  df_mda <- df_mda[order(as.numeric(rownames(df_mda))), ]

  # --- Concordance over originally-missing positions only ---
  list_of_dfs <- list(OT = df_ot, kNN = df_knn, MICE = df_mice, MF = df_mf, MDA = df_mda)

  res <- data.frame(
    Scenario   = scenario_name,
    OT_kNN_Y   = calc_concordance(df_ot$Y,    df_knn$Y,   mask_Y),
    OT_kNN_Z   = calc_concordance(df_ot$Z,    df_knn$Z,   mask_Z),
    OT_MICE_Y  = calc_concordance(df_ot$Y,    df_mice$Y,  mask_Y),
    OT_MICE_Z  = calc_concordance(df_ot$Z,    df_mice$Z,  mask_Z),
    OT_MF_Y    = calc_concordance(df_ot$Y,    df_mf$Y,    mask_Y),
    OT_MF_Z    = calc_concordance(df_ot$Z,    df_mf$Z,    mask_Z),
    OT_MDA_Y   = calc_concordance(df_ot$Y,    df_mda$Y,   mask_Y),
    OT_MDA_Z   = calc_concordance(df_ot$Z,    df_mda$Z,   mask_Z),
    kNN_MICE_Y = calc_concordance(df_knn$Y,   df_mice$Y,  mask_Y),
    kNN_MICE_Z = calc_concordance(df_knn$Z,   df_mice$Z,  mask_Z),
    kNN_MF_Y   = calc_concordance(df_knn$Y,   df_mf$Y,    mask_Y),
    kNN_MF_Z   = calc_concordance(df_knn$Z,   df_mf$Z,    mask_Z),
    kNN_MDA_Y  = calc_concordance(df_knn$Y,   df_mda$Y,   mask_Y),
    kNN_MDA_Z  = calc_concordance(df_knn$Z,   df_mda$Z,   mask_Z),
    MICE_MF_Y  = calc_concordance(df_mice$Y,  df_mf$Y,    mask_Y),
    MICE_MF_Z  = calc_concordance(df_mice$Z,  df_mf$Z,    mask_Z),
    MICE_MDA_Y = calc_concordance(df_mice$Y,  df_mda$Y,   mask_Y),
    MICE_MDA_Z = calc_concordance(df_mice$Z,  df_mda$Z,   mask_Z),
    MF_MDA_Y   = calc_concordance(df_mf$Y,    df_mda$Y,   mask_Y),
    MF_MDA_Z   = calc_concordance(df_mf$Z,    df_mda$Z,   mask_Z),
    ALL_Y      = calc_multi_concordance(list_of_dfs, "Y", mask_Y),
    ALL_Z      = calc_multi_concordance(list_of_dfs, "Z", mask_Z),
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

# ------------------------
# 5. LaTeX Table Output
# ------------------------
fmt_val <- function(x) if (is.na(x)) "--" else sprintf("%.2f", x)

cat("\n\n")
cat("\\begin{table}[h!]\n")
cat("\\caption{Concordance percentages among OTrecod, kNN, MICE, missForest (MF), and missMDA (MDA) for $Y_A$ and $Y_B$}\n")
cat("\\centering\n")
cat("\\resizebox{\\textwidth}{!}{\n")
cat("\\begin{tabular}{l *{11}{c c} c c}\n")
cat("\\toprule\n")
cat(" & \\multicolumn{2}{c}{OT-kNN} & \\multicolumn{2}{c}{OT-MICE} & \\multicolumn{2}{c}{OT-MF} & \\multicolumn{2}{c}{OT-MDA} & \\multicolumn{2}{c}{kNN-MICE} & \\multicolumn{2}{c}{kNN-MF} & \\multicolumn{2}{c}{kNN-MDA} & \\multicolumn{2}{c}{MICE-MF} & \\multicolumn{2}{c}{MICE-MDA} & \\multicolumn{2}{c}{MF-MDA} & \\multicolumn{2}{c}{All 5} & \\multicolumn{2}{c}{$n$} \\\\\n")
cat("\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9} \\cmidrule(lr){10-11} \\cmidrule(lr){12-13} \\cmidrule(lr){14-15} \\cmidrule(lr){16-17} \\cmidrule(lr){18-19} \\cmidrule(lr){20-21} \\cmidrule(lr){22-23} \\cmidrule(lr){24-25}\n")
cat("Cases & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ & $Y_A$ & $Y_B$ \\\\\n")
cat("\\midrule\n")

for (i in 1:nrow(final_results)) {
  r <- final_results[i, ]
  cat(sprintf(
    "%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %d & %d \\\\\n",
    r$Scenario,
    fmt_val(r$OT_kNN_Y),   fmt_val(r$OT_kNN_Z),
    fmt_val(r$OT_MICE_Y),  fmt_val(r$OT_MICE_Z),
    fmt_val(r$OT_MF_Y),    fmt_val(r$OT_MF_Z),
    fmt_val(r$OT_MDA_Y),   fmt_val(r$OT_MDA_Z),
    fmt_val(r$kNN_MICE_Y), fmt_val(r$kNN_MICE_Z),
    fmt_val(r$kNN_MF_Y),   fmt_val(r$kNN_MF_Z),
    fmt_val(r$kNN_MDA_Y),  fmt_val(r$kNN_MDA_Z),
    fmt_val(r$MICE_MF_Y),  fmt_val(r$MICE_MF_Z),
    fmt_val(r$MICE_MDA_Y), fmt_val(r$MICE_MDA_Z),
    fmt_val(r$MF_MDA_Y),   fmt_val(r$MF_MDA_Z),
    fmt_val(r$ALL_Y),      fmt_val(r$ALL_Z),
    r$n_Y, r$n_Z
  ))
}

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("}\n")
cat("\\label{tab:agreement_split}\n")
cat("\\end{table}\n")
cat("}\n")
cat("\\label{tab:agreement_split}\n")
cat("\\end{table}\n")