suppressPackageStartupMessages({
  library(MASS)
  library(missMDA)
  library(dplyr)
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

# ------------------------
# 2. Asymmetry Simulation Function (missMDA Only)
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
  if (R2 >= 1) {
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

  # Ensure categorical variables are strictly factors
  df_impute <- df[, -1] %>% mutate_if(is.character, as.factor)

  # --- missMDA (FAMD) ---
  capture.output(famd_res <- try(imputeFAMD(df_impute, ncp = 2), silent = TRUE))

  if(!inherits(famd_res, "try-error")) {

    # CRITICAL FIX: Strip the "YA." and "YB." prefixes missMDA adds to factor levels
    # By replacing all non-numeric characters with "", "YA.1" perfectly becomes "1"
    imp_YA_clean <- gsub("[^0-9]", "", as.character(famd_res$completeObs$YA))
    imp_YB_clean <- gsub("[^0-9]", "", as.character(famd_res$completeObs$YB))

    res_famd <- calc_metrics(YA_true, YB_true, imp_YA_clean, imp_YB_clean, mask_A, mask_B)
  } else {
    res_famd <- c(Precision = NA, Kappa = NA, MAE = NA)
  }

  return(data.frame(
    Method = "missMDA",
    Precision = res_famd["Precision"],
    Kappa = res_famd["Kappa"],
    MAE = res_famd["MAE"],
    stringsAsFactors = FALSE,
    row.names = NULL
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
              }, error = function(e) { NULL })

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
            cat(sprintf("Completed: n=%d, K_A=%d, K_B=%d (Seed=%d)\n", n, K_A, K_B, current_seed))
          }
        }
      }
    }
  }

  # Aggregate results
  final_data <- bind_rows(results_list)

  summary_results <- final_data %>%
    group_by(n, perc_missing, R2, Total_K, K_A, K_B, Difference_K, Method) %>%
    summarize(
      total_successful_runs = sum(!is.na(Precision)),
      mean_precision = mean(Precision, na.rm = TRUE),
      err_precision = 1.96 * (sd(Precision, na.rm = TRUE) / sqrt(sum(!is.na(Precision)))),
      mean_kappa = mean(Kappa, na.rm = TRUE),
      err_kappa = 1.96 * (sd(Kappa, na.rm = TRUE) / sqrt(sum(!is.na(Kappa)))),
      mean_mae = mean(MAE, na.rm = TRUE),
      err_mae = 1.96 * (sd(MAE, na.rm = TRUE) / sqrt(sum(!is.na(MAE)))),
      .groups = 'drop'
    )

  write.csv(final_data, "../asymmetry_raw_results_missMDA.csv", row.names = FALSE)
  write.csv(summary_results, "../asymmetry_summary_results_missMDA.csv", row.names = FALSE)
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