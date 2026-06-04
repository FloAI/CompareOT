suppressPackageStartupMessages({
  require(MASS)
  require(VIM)
  require(mice)
  require(OTrecod)
  require(missForest)
  require(missMDA)
  require(dplyr)
})

# ------------------------
# 1. Helper function for metrics
# ------------------------
calc_metrics <- function(actual_YA, actual_YB, imp_YA, imp_YB, mask_A, mask_B) {

  act1 <- as.character(actual_YA[mask_A])
  imp1 <- as.character(imp_YA[mask_A])
  act2 <- as.character(actual_YB[mask_B])
  imp2 <- as.character(imp_YB[mask_B])

  act_all <- c(act1, act2)
  imp_all <- c(imp1, imp2)

  # Precision (Accuracy)
  precision <- sum(act_all == imp_all, na.rm=TRUE) / length(act_all)

  # Cohen's Kappa
  all_levels <- unique(c(act_all, imp_all))
  act_f <- factor(act_all, levels = all_levels)
  imp_f <- factor(imp_all, levels = all_levels)
  tab <- table(imp_f, act_f)
  p_o <- sum(diag(tab)) / sum(tab)
  p_e <- sum(rowSums(tab) * colSums(tab)) / (sum(tab)^2)
  kappa_val <- ifelse(p_e == 1, 0, (p_o - p_e) / (1 - p_e))

  # MAE (Using raw numeric levels)
  mae_val <- mean(abs(as.numeric(act_all) - as.numeric(imp_all)), na.rm = TRUE)

  return(c(Precision = precision, Kappa = kappa_val, MAE = mae_val))
}

get_mode <- function(v) {
  uniqv <- unique(na.omit(v))
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# ------------------------
# 2. Asymmetry Simulation Function
# ------------------------
simulate_asymmetry_rep <- function(n, p_miss, R2, K_A, total_K) {
  p <- 5

  # Calculate K_B based on the fixed sum
  K_B <- total_K - K_A
  if (K_B < 2) stop("K_B must be at least 2 for classification.")

  # LINEAR CASE: Independent continuous covariates
  X <- matrix(rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
  X_df <- as.data.frame(X)
  colnames(X_df) <- paste0("X", 1:p)

  # Linear outcome generation
  beta <- c(20, 5, 10, 2, -8, 0.3)
  Z_hat <- beta[1] + X %*% beta[2:6]

  # R2 = 1 handling (zero noise)
  if (R2 == 1) {
    sigma2 <- 0
  } else {
    var_zhat <- var(as.vector(Z_hat))
    sigma2 <- var_zhat * (1 - R2) / R2
  }

  Z <- Z_hat + rnorm(n, 0, sqrt(sigma2))

  # ASYMMETRIC MODALITIES:
  YA_true <- cut(Z, breaks = quantile(Z, probs = seq(0, 1, length.out = K_A + 1)), labels = 1:K_A, include.lowest = TRUE)
  YB_true <- cut(Z, breaks = quantile(Z, probs = seq(0, 1, length.out = K_B + 1)), labels = 1:K_B, include.lowest = TRUE)

  # Introduce Missingness based on p_miss
  missing_mask <- runif(n) < p_miss
  DB <- as.factor(ifelse(missing_mask, "A", "B"))

  YA <- YA_true
  YB <- YB_true
  YA[DB == "B"] <- NA
  YB[DB == "A"] <- NA

  df <- data.frame(DB = DB, YA = YA, YB = YB)
  df <- cbind(df, X_df)

  mask_A <- is.na(df$YA)
  mask_B <- is.na(df$YB)

  # --- 0. Baseline (Mode) ---
  df_base <- df
  df_base$YA[is.na(df_base$YA)] <- get_mode(df$YA)
  df_base$YB[is.na(df_base$YB)] <- get_mode(df$YB)
  res_base <- calc_metrics(YA_true, YB_true, df_base$YA, df_base$YB, mask_A, mask_B)

  # --- 1. missForest ---
  capture.output(mf_res <- missForest(df[, -1], maxiter = 5, ntree = 50))
  res_mf <- calc_metrics(YA_true, YB_true, mf_res$ximp$YA, mf_res$ximp$YB, mask_A, mask_B)

  # --- 2. kNN ---
  capture.output(df_knn <- kNN(df[, -1], k = 5, imp_var = FALSE))
  res_knn <- calc_metrics(YA_true, YB_true, df_knn$YA, df_knn$YB, mask_A, mask_B)

  # --- 3. MICE ---
  capture.output({
    init <- mice(df[, -1], maxit=0)
    meth <- init$method
    meth[c("YA", "YB")] <- "rf"
    imputed_mice <- mice(df[, -1], method=meth, m=1, maxit=5, printFlag=FALSE)
    df_mice <- complete(imputed_mice)
  })
  res_mice <- calc_metrics(YA_true, YB_true, df_mice$YA, df_mice$YB, mask_A, mask_B)

  # --- 4. OTrecod ---
  capture.output(R_OUTC1 <- try(OT_joint(df,  prox.X=0.10, convert.num=8, convert.class=1,
        nominal=c(1,4:7), ordinal=2:3, dist.choice="H",
        maxrelax=0.1, which.DB="BOTH"), silent = TRUE))
  df_ot <- df
  if(!inherits(R_OUTC1, "try-error")) {
    df_ot$YB[is.na(df_ot$YB)] <- R_OUTC1$DATA1_OT[, "OTpred"]
    df_ot$YA[is.na(df_ot$YA)] <- R_OUTC1$DATA2_OT[, "OTpred"]
  }
  res_ot <- calc_metrics(YA_true, YB_true, df_ot$YA, df_ot$YB, mask_A, mask_B)

  # --- 5. missMDA (FAMD) FIXED ---
  capture.output({
    famd_res <- try(imputeFAMD(df[, -1], ncp = 2), silent = TRUE)
  })

  if(!inherits(famd_res, "try-error")) {
    res_famd <- calc_metrics(YA_true, YB_true, famd_res$completeObs$YA, famd_res$completeObs$YB, mask_A, mask_B)
  } else {
    # If SVD fails due to block missingness, gracefully assign NAs instead of crashing the run
    res_famd <- c(Precision = NA, Kappa = NA, MAE = NA)
  }

  methods <- c("Baseline", "missForest", "kNN", "MICE", "OTrecod", "missMDA")
  return(data.frame(
    Method = methods,
    Precision = c(res_base["Precision"], res_mf["Precision"], res_knn["Precision"], res_mice["Precision"], res_ot["Precision"], res_famd["Precision"]),
    Kappa = c(res_base["Kappa"], res_mf["Kappa"], res_knn["Kappa"], res_mice["Kappa"], res_ot["Kappa"], res_famd["Kappa"]),
    MAE = c(res_base["MAE"], res_mf["MAE"], res_knn["MAE"], res_mice["MAE"], res_ot["MAE"], res_famd["MAE"])
  ))
}

# ------------------------
# 3. Execution Wrapper
# ------------------------
simulate_asymmetry_scenario <- function(n_list, p_miss_list, R2_list, K_A_list, total_K, n_sim, seed_list) {
  results_list <- list()

  for (current_seed in seed_list) {
    set.seed(current_seed)

    for (n in n_list) {
      for (p_miss in p_miss_list) {
        for (R2 in R2_list) {
          for (K_A in K_A_list) {

            K_B <- total_K - K_A
            diff_K <- abs(K_A - K_B)

            for (rep_i in 1:n_sim) {
              sim_res <- tryCatch({
                simulate_asymmetry_rep(n, p_miss, R2, K_A, total_K)
              }, error = function(e) {
                message(sprintf("Error rep %d (Seed=%d, K_A=%d): %s", rep_i, current_seed, K_A, e$message))
                return(NULL)
              })

              if (!is.null(sim_res)) {
                sim_res$n <- n
                sim_res$perc_missing <- p_miss
                sim_res$R2 <- R2
                sim_res$K_A <- K_A
                sim_res$K_B <- K_B
                sim_res$Difference_K <- diff_K
                sim_res$Total_K <- total_K
                sim_res$Seed <- current_seed
                sim_res$Repetition <- rep_i
                results_list[[length(results_list) + 1]] <- sim_res
              }
            }
            cat(sprintf("Completed: n=%d, K_A=%d, K_B=%d (Seed=%d, Runs=%d)\n", n, K_A, K_B, current_seed, n_sim))
          }
        }
      }
    }
  }

  # Aggregate results
  final_data <- do.call(rbind, results_list)

  summary_results <- final_data %>%
    group_by(n, perc_missing, R2, Total_K, K_A, K_B, Difference_K, Method) %>%
    summarize(
      total_runs = n(),
      mean_precision = mean(Precision, na.rm = TRUE),
      se_precision = sd(Precision, na.rm = TRUE) / sqrt(n()),
      mean_kappa = mean(Kappa, na.rm = TRUE),
      se_kappa = sd(Kappa, na.rm = TRUE) / sqrt(n()),
      mean_mae = mean(MAE, na.rm = TRUE),
      se_mae = sd(MAE, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    ) %>%
    mutate(
      margin_err_precision = 1.96 * se_precision,
      margin_err_kappa = 1.96 * se_kappa,
      margin_err_mae = 1.96 * se_mae
    ) %>%
    select(-starts_with("se_"))

  write.csv(final_data, "asymmetry_raw_results_multiseed2.csv", row.names = FALSE)
  write.csv(summary_results, "asymmetry_summary_results_multiseed2.csv", row.names = FALSE)
  cat("\nSIMULATION COMPLETE AND SAVED TO CSV!\n")

  return(summary_results)
}

# ------------------------
# 4. Run the simulation
# ------------------------
final_summary_asymmetry <- simulate_asymmetry_scenario(
  n_list = c(500),
  p_miss_list = c(0.5),
  R2_list = c(1),
  total_K = 7,
  K_A_list = c(2, 3, 4, 5),
  n_sim = 6,
  seed_list = c(42, 101, 555, 999, 2026)
)

print(final_summary_asymmetry)