suppressPackageStartupMessages({
  require(MASS)
  require(norm)
  require(VIM)
  require(ggplot2)
  require(naniar)
  require(mice)
  require(OTrecod)
  require(missForest)
  require(caret)
  require(dplyr)
  require(missMDA)
})

source('https://raw.githubusercontent.com/R-miss-tastic/website/master/static/how-to/generate/amputation.R')
set.seed(42)

# ------------------------
# Define simulation parameters (From Code 1 Structure)
# ------------------------
num_samples_list <- c(100, 200, 300, 400, 500, 1000)
perc_missing_list <- c(0.5)
target_R2_list <- c(0.9)
n_reps <- 30  # repetitions per combination

# ------------------------
# Helper function for metrics (Precision, Kappa, MAE)
# ------------------------
calc_metrics <- function(actual_df, imputed_df, m_mask_yb1, m_mask_yb2) {

  # Extract missing subsets
  act1 <- as.character(actual_df$Yb1[m_mask_yb1])
  imp1 <- as.character(imputed_df$Yb1[m_mask_yb1])

  act2 <- as.character(actual_df$Yb2[m_mask_yb2])
  imp2 <- as.character(imputed_df$Yb2[m_mask_yb2])

  act_all <- c(act1, act2)
  imp_all <- c(imp1, imp2)

  # 1. Precision (Accuracy)
  precision <- sum(act_all == imp_all, na.rm=TRUE) / length(act_all)

  # 2. Cohen's Kappa
  all_levels <- unique(c(act_all, imp_all))
  act_f <- factor(act_all, levels = all_levels)
  imp_f <- factor(imp_all, levels = all_levels)

  # Calculate expected frequency and observed accuracy for Kappa
  tab <- table(imp_f, act_f)
  p_o <- sum(diag(tab)) / sum(tab)
  p_e <- sum(rowSums(tab) * colSums(tab)) / (sum(tab)^2)
  kappa_val <- ifelse(p_e == 1, 0, (p_o - p_e) / (1 - p_e))

  # 3. MAE (requires numeric mapping for ordinals)
  map_yb1 <- c("[Q1]"=1, "[Q2]"=2, "[Q3]"=3, "[Q4]"=4)
  map_yb2 <- c("[Q1]"=1, "[Q2]"=2, "[Q3]"=3)

  act_num <- c(map_yb1[act1], map_yb2[act2])
  imp_num <- c(map_yb1[imp1], map_yb2[imp2])
  mae_val <- mean(abs(act_num - imp_num), na.rm = TRUE)

  return(c(Precision = precision, Kappa = kappa_val, MAE = mae_val))
}

# ------------------------
# Core Simulation Function (Strictly faithful to Code 2 data generation)
# ------------------------
run_sim_combination <- function(num_samples, perc_missing, target_R2) {

  # Generate covariates using generic variable names X1-X5
  X1 <- rbinom(num_samples, 1, 0.5)  # binary
  X2_raw <- t(rmultinom(num_samples, 1, prob = c(0.3, 0.4, 0.3)))
  X2 <- apply(X2_raw, 1, function(x) c("A","B","Placebo")[which(x==1)])

  X3_raw <- t(rmultinom(num_samples, 1, prob = c(0.1,0.2,0.3,0.4)))
  X3 <- apply(X3_raw, 1, function(x) paste0("L", which(x==1)))

  X4 <- rbinom(num_samples, 1, 0.5)  # binary
  X5 <- runif(num_samples, 20, 80)   # continuous

  # Define the true signal
  X1_num <- X1
  X2_num <- ifelse(X2=="A",1,ifelse(X2=="B",2,0))
  X3_num <- as.numeric(factor(X3))
  X4_num <- X4
  X5_num <- X5

  signal <- 20 + 5*X1_num + 10*X2_num + 2*X3_num - 8*X4_num + 0.3*X5_num

  # Compute sigma for target R²
  sigma <- sqrt(var(signal) * (1/target_R2 - 1))

  # Generate continuous outcomes
  Yb1_cont <- signal + sigma*rnorm(num_samples)
  Yb2_cont <- signal + sigma*rnorm(num_samples)

  # Categorize outcomes strictly as Code 2
  categorize_Yb1 <- cut(Yb1_cont,
                  breaks = quantile(Yb1_cont, probs = seq(0,1,.25)),
                  include.lowest = TRUE, labels = c("[Q1]","[Q2]","[Q3]","[Q4]"))
  categorize_Yb2 <- cut(Yb2_cont,
                  breaks = quantile(Yb2_cont, probs = seq(0,1,1/3)),
                  include.lowest = TRUE, labels = c("[Q1]","[Q2]","[Q3]"))

  Yb1 <- sapply(Yb1_cont, categorize_Yb1)
  Yb2 <- sapply(Yb2_cont, categorize_Yb2)
  df_true <- data.frame(Yb1, Yb2)

  df <- data.frame(X1, X2, X3, X4, X5)

  # Introduce missingness strictly as Code 2
  missing_mask <- runif(num_samples) < perc_missing
  DB <- ifelse(missing_mask, "A", "B")
  Yb1[DB=="B"] <- NA
  Yb2[DB=="A"] <- NA

  df <- cbind(DB, Yb1, Yb2, df)
  df <- df[order(df$DB), ] # This sort changes row indices

  df$DB <- as.factor(df$DB)
  df$Yb1 <- as.factor(df$Yb1)
  df$Yb2 <- as.factor(df$Yb2)
  df$X1 <- as.factor(df$X1)
  df$X2 <- as.factor(df$X2)
  df$X3 <- as.factor(df$X3)
  df$X4 <- as.factor(df$X4)
  df$X5 <- as.numeric(df$X5)

  # Save mask configuration based on sorting behavior in Code 2
  df2 <- df[, c("Yb1","Yb2")]
  df2 <- df2[order(as.numeric(rownames(df2))), ]
  m_mask_yb1 <- is.na(df2$Yb1)
  m_mask_yb2 <- is.na(df2$Yb2)

  # Helper function to align results back to df_true order for all methods
  extract_and_sort <- function(imputed_full_df) {
    temp <- imputed_full_df[, c("Yb1", "Yb2")]
    temp <- temp[order(as.numeric(rownames(temp))), ]
    return(temp)
  }

  # --- 1. OTrecod (Strictly executing Code 2 implementation) ---
  capture.output(
    R_OUTC1 <- OT_joint(df,
                        prox.X = 0.10,
                        convert.num = 8, convert.class = 1,
                        nominal = c(1,4:7), ordinal = 2:3,
                        dist.choice = "H",
                        maxrelax = 0.1,
                        which.DB = "BOTH")
  )

  y_imputed <- R_OUTC1$DATA2_OT[, "OTpred"]
  z_imputed <- R_OUTC1$DATA1_OT[, "OTpred"]

  df_ot_temp <- df
  df_ot_temp$Yb2[is.na(df_ot_temp$Yb2)] <- z_imputed
  df_ot_temp$Yb1[is.na(df_ot_temp$Yb1)] <- y_imputed

  imputated_df_ot <- extract_and_sort(df_ot_temp)
  res_ot <- calc_metrics(df_true, imputated_df_ot, m_mask_yb1, m_mask_yb2)

  # --- 2. missForest ---
  capture.output(mf_res <- missForest(df, maxiter = 5, ntree = 50))
  res_mf <- calc_metrics(df_true, extract_and_sort(mf_res$ximp), m_mask_yb1, m_mask_yb2)

  # --- 3. kNN ---
  capture.output(df_knn <- kNN(df, k = 5, imp_var = FALSE))
  res_knn <- calc_metrics(df_true, extract_and_sort(df_knn), m_mask_yb1, m_mask_yb2)

  # --- 4. MICE ---
  capture.output({
    init <- mice(df, maxit=0)
    meth <- init$method
    meth[c("Yb1", "Yb2")] <- "rf"
    imputed_mice <- mice(df, method=meth, m=1, maxit=5, printFlag=FALSE)
    df_mice <- complete(imputed_mice)
  })
  res_mice <- calc_metrics(df_true, extract_and_sort(df_mice), m_mask_yb1, m_mask_yb2)

  # --- 5. missMDA (FAMD) ---
  capture.output({
    mda_res <- imputeFAMD(df, ncp = 2)
    df_mda <- mda_res$completeObs
  })
  res_mda <- calc_metrics(df_true, extract_and_sort(df_mda), m_mask_yb1, m_mask_yb2)

  # Combine results
  methods <- c("missForest", "kNN", "MICE", "OTrecod", "missMDA")
  return(data.frame(
    Method = methods,
    Precision = c(res_mf["Precision"], res_knn["Precision"], res_mice["Precision"], res_ot["Precision"], res_mda["Precision"]),
    Kappa = c(res_mf["Kappa"], res_knn["Kappa"], res_mice["Kappa"], res_ot["Kappa"], res_mda["Kappa"]),
    MAE = c(res_mf["MAE"], res_knn["MAE"], res_mice["MAE"], res_ot["MAE"], res_mda["MAE"])
  ))
}

# ------------------------
# Run simulation for all combinations (From Code 1 Structure)
# ------------------------
results_list <- list()

for (n_samp in num_samples_list) {
  for (p_miss in perc_missing_list) {
    for (R2 in target_R2_list) {

      for (rep_i in 1:n_reps) {
        sim_res <- tryCatch({
          run_sim_combination(n_samp, p_miss, R2)
        }, error = function(e) {
          message(sprintf("Error in rep %d (n=%d, miss=%.1f, R2=%.1f): %s", rep_i, n_samp, p_miss, R2, e$message))
          return(NULL)
        })

        if (!is.null(sim_res)) {
          sim_res$num_samples <- n_samp
          sim_res$perc_missing <- p_miss
          sim_res$target_R2 <- R2
          results_list[[length(results_list) + 1]] <- sim_res
        }
      }
      cat(sprintf("Completed: Samples=%d, Missing=%.1f, R2=%.1f\n", n_samp, p_miss, R2))
    }
  }
}

# ------------------------
# Aggregate and Summarize Results
# ------------------------
final_data <- do.call(rbind, results_list)

summary_results <- final_data %>%
  group_by(num_samples, perc_missing, target_R2, Method) %>%
  summarize(
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

print(summary_results)

write.csv(final_data, "../simulation_raw_results_n.csv", row.names = FALSE)
write.csv(summary_results, "../simulation_summary_results_n.csv", row.names = FALSE)
