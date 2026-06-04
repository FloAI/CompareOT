suppressPackageStartupMessages({
  require(MASS)
  require(VIM)
  require(mice)
  require(missForest)
  require(missMDA)
  require(dplyr)
})

# ------------------------
# 1. Helper function for metrics
# ------------------------
calc_metrics <- function(actual_df, imputed_df, m_mask_yb1, m_mask_yb2) {
  act1 <- as.character(actual_df$Yb1[m_mask_yb1])
  imp1 <- as.character(imputed_df$Yb1[m_mask_yb1])
  act2 <- as.character(actual_df$Yb2[m_mask_yb2])
  imp2 <- as.character(imputed_df$Yb2[m_mask_yb2])

  act_all <- c(act1, act2)
  imp_all <- c(imp1, imp2)

  precision <- sum(act_all == imp_all, na.rm=TRUE) / length(act_all)

  all_levels <- unique(c(act_all, imp_all))
  act_f <- factor(act_all, levels = all_levels)
  imp_f <- factor(imp_all, levels = all_levels)
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
# 2. Data generation function
# ------------------------
generate_dataset <- function(num_samples, perc_missing, target_R2) {
  X1 <- rbinom(num_samples, 1, 0.5)
  X2_raw <- t(rmultinom(num_samples, 1, prob = c(0.3, 0.4, 0.3)))
  X2 <- apply(X2_raw, 1, function(x) c("A","B","Placebo")[which(x==1)])
  X3_raw <- t(rmultinom(num_samples, 1, prob = c(0.1,0.2,0.3,0.4)))
  X3 <- apply(X3_raw, 1, function(x) paste0("L", which(x==1)))
  X4 <- rbinom(num_samples, 1, 0.5)
  X5 <- runif(num_samples, 20, 80)

  X1_num <- X1
  X2_num <- ifelse(X2=="A",1,ifelse(X2=="B",2,0))
  X3_num <- as.numeric(factor(X3))

  signal <- 20 + 5*X1_num + 10*X2_num + 2*X3_num - 8*X4 + 0.3*X5
  sigma <- sqrt(var(signal) * (1/target_R2 - 1))

  Yb1_cont <- signal + sigma*rnorm(num_samples)
  Yb2_cont <- signal + sigma*rnorm(num_samples)

  categorize_Yb1 <- cut(Yb1_cont,
                  breaks = quantile(Yb1_cont, probs = seq(0,1,.25)),
                  include.lowest = TRUE, labels = c("[Q1]","[Q2]","[Q3]","[Q4]"))
  categorize_Yb2 <- cut(Yb2_cont,
                  breaks = quantile(Yb2_cont, probs = seq(0,1,1/3)),
                  include.lowest = TRUE, labels = c("[Q1]","[Q2]","[Q3]"))

  df_true <- data.frame(Yb1 = sapply(Yb1_cont, categorize_Yb1),
                        Yb2 = sapply(Yb2_cont, categorize_Yb2))

  df <- data.frame(X1, X2, X3, X4, X5)

  missing_mask <- runif(num_samples) < perc_missing
  DB <- ifelse(missing_mask, "A", "B")
  Yb1 <- df_true$Yb1; Yb2 <- df_true$Yb2
  Yb1[DB=="B"] <- NA; Yb2[DB=="A"] <- NA

  df <- cbind(DB, Yb1, Yb2, df)
  df <- df[order(df$DB), ]

  df$DB  <- as.factor(df$DB)
  df$Yb1 <- as.factor(df$Yb1)
  df$Yb2 <- as.factor(df$Yb2)
  df$X1  <- as.factor(df$X1)
  df$X2  <- as.factor(df$X2)
  df$X3  <- as.factor(df$X3)
  df$X4  <- as.factor(df$X4)
  df$X5  <- as.numeric(df$X5)

  df_true <- df_true[order(as.numeric(rownames(df_true))), ]
  m_mask_yb1 <- is.na(df$Yb1)
  m_mask_yb2 <- is.na(df$Yb2)

  return(list(df = df, df_true = df_true, m_mask_yb1 = m_mask_yb1, m_mask_yb2 = m_mask_yb2))
}

# ------------------------
# 3. Sensitivity Analysis
# ------------------------
run_sensitivity_analysis <- function(n_sim = 10) {
  results_list <- list()

  # Fixed data parameters (middle of your simulation grid)
  num_samples <- 100
  perc_missing <- 0.1
  target_R2    <- 0.9

  for (rep_i in 1:n_sim) {

    dat        <- generate_dataset(num_samples, perc_missing, target_R2)
    df         <- dat$df
    df_true    <- dat$df_true
    m_mask_yb1 <- dat$m_mask_yb1
    m_mask_yb2 <- dat$m_mask_yb2

    # -----------------------------------------------------
    # TEST 1: missMDA (FAMD) - Tuning 'ncp'
    # -----------------------------------------------------
    for (ncp_val in c(1, 2, 3, 4, 5)) {
      capture.output(famd_res <- try(imputeFAMD(df[, -1], ncp = ncp_val), silent = TRUE))
      if (!inherits(famd_res, "try-error")) {
        imp_df <- df
        imp_df$Yb1 <- factor(famd_res$completeObs$Yb1, levels = levels(df$Yb1))
        imp_df$Yb2 <- factor(famd_res$completeObs$Yb2, levels = levels(df$Yb2))
        res <- calc_metrics(df_true, imp_df, m_mask_yb1, m_mask_yb2)
        results_list[[length(results_list)+1]] <- data.frame(
          Rep=rep_i, Method="missMDA", Param="ncp", Value=ncp_val,
          Precision=res[1], Kappa=res[2], MAE=res[3]
        )
      }
    }

    # -----------------------------------------------------
    # TEST 2: kNN - Tuning 'k'
    # -----------------------------------------------------
    for (k_val in c(3, 5, 10, 15)) {
      capture.output(knn_res <- kNN(df, k = k_val, imp_var = FALSE))
      res <- calc_metrics(df_true, knn_res, m_mask_yb1, m_mask_yb2)
      results_list[[length(results_list)+1]] <- data.frame(
        Rep=rep_i, Method="kNN", Param="k", Value=k_val,
        Precision=res[1], Kappa=res[2], MAE=res[3]
      )
    }

    # -----------------------------------------------------
    # TEST 3: missForest - Tuning 'ntree'
    # -----------------------------------------------------
    for (ntree_val in c(10, 50, 100)) {
      capture.output(mf_res <- missForest(df, maxiter = 5, ntree = ntree_val))
      res <- calc_metrics(df_true, mf_res$ximp, m_mask_yb1, m_mask_yb2)
      results_list[[length(results_list)+1]] <- data.frame(
        Rep=rep_i, Method="missForest", Param="ntree", Value=ntree_val,
        Precision=res[1], Kappa=res[2], MAE=res[3]
      )
    }

    # -----------------------------------------------------
    # TEST 4: MICE - Tuning 'm' (number of imputations)
    # -----------------------------------------------------
    for (m_val in c(1, 3, 5, 10)) {
      capture.output({
        init <- mice(df, maxit = 0)
        meth <- init$method
        meth[c("Yb1", "Yb2")] <- "rf"
        imputed_mice <- mice(df, method = meth, m = m_val, maxit = 5, printFlag = FALSE)
        df_mice <- complete(imputed_mice)
      })
      res <- calc_metrics(df_true, df_mice, m_mask_yb1, m_mask_yb2)
      results_list[[length(results_list)+1]] <- data.frame(
        Rep=rep_i, Method="MICE", Param="m", Value=m_val,
        Precision=res[1], Kappa=res[2], MAE=res[3]
      )
    }

    cat(sprintf("Completed Hyperparameter Repetition %d of %d\n", rep_i, n_sim))
  }

  # Aggregate
  final_data <- do.call(rbind, results_list)

  summary_results <- final_data %>%
    group_by(Method, Param, Value) %>%
    summarize(
      mean_precision = mean(Precision, na.rm = TRUE),
      mean_kappa     = mean(Kappa,     na.rm = TRUE),
      mean_mae       = mean(MAE,       na.rm = TRUE),
      .groups = 'drop'
    )

  write.csv(summary_results, "../hyperparameter_sensitivity_summary_100_2.csv", row.names = FALSE)
  return(summary_results)
}

# Run
sensitivity_results <- run_sensitivity_analysis(n_sim = 10)
print(sensitivity_results)
