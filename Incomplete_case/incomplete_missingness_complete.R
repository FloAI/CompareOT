suppressPackageStartupMessages({
  library(MASS)
  library(VIM)
  library(mice)
  library(OTrecod)
  library(missForest)
  library(missMDA)
  library(dplyr)
})

# Source the ampute function
source('https://raw.githubusercontent.com/R-miss-tastic/website/master/static/how-to/generate/amputation.R')

# --- Configuration ---
seed_list <- c(42, 123, 456, 789, 999)
n_reps <- 1
num_samples_list <- c(100, 200, 300, 400, 500, 1000)
perc_missing_list <- seq(0.1, 0.9, by = 0.4)
mechanism_list <- c("MCAR", "MAR", "MNAR")
target_R2 <- 1.0
OT_MAX_RETRIES <- 10  # Max attempts to get a valid OTrecod result per rep

# ------------------------
# Helper function for metrics
# ------------------------
calc_metrics <- function(actual, imputed, mask) {
  act <- as.character(actual[mask])
  imp <- as.character(imputed[mask])

  precision <- sum(act == imp, na.rm=TRUE) / length(act)

  all_levels <- unique(c(act, imp))
  act_f <- factor(act, levels = all_levels)
  imp_f <- factor(imp, levels = all_levels)
  tab <- table(imp_f, act_f)
  p_o <- sum(diag(tab)) / sum(tab)
  p_e <- sum(rowSums(tab) * colSums(tab)) / (sum(tab)^2)
  kappa_val <- ifelse(p_e == 1, 0, (p_o - p_e) / (1 - p_e))

  act_num <- as.numeric(factor(act, levels = sort(unique(act))))
  imp_num <- as.numeric(factor(imp, levels = sort(unique(act))))
  mae_val <- mean(abs(act_num - imp_num), na.rm = TRUE)

  return(c(Precision = precision, Kappa = kappa_val, MAE = mae_val))
}

# ------------------------
# Data generation function (separated for retry use)
# ------------------------
generate_data <- function(num_samples, perc_missing, mech) {

  X1 <- rbinom(num_samples, 1, 0.5)
  X2_raw <- t(rmultinom(num_samples, 1, prob = c(0.3,0.4,0.3)))
  X2 <- apply(X2_raw, 1, function(x) c("A","B","Placebo")[which(x==1)])
  X3_raw <- t(rmultinom(num_samples, 1, prob = c(0.1,0.2,0.3,0.4)))
  X3 <- apply(X3_raw, 1, function(x) paste0("L", which(x==1)))
  X4 <- rbinom(num_samples, 1, 0.5)
  X5 <- runif(num_samples, 20, 80)

  X1_num <- X1
  X2_num <- ifelse(X2=="A",1,ifelse(X2=="B",2,0))
  X3_num <- as.numeric(factor(X3))
  X4_num <- X4
  X5_num <- X5

  signal <- 20 + 5*X1_num + 10*X2_num + 2*X3_num - 8*X4_num + 0.3*X5_num

  if(target_R2 >= 1.0) {
    sigma <- 0
  } else {
    sigma <- sqrt(var(signal) * (1/target_R2 - 1))
  }

  Yb1_cont <- signal + sigma*rnorm(num_samples)

  categorize_Yb1 <- cut(Yb1_cont,
                  breaks = quantile(Yb1_cont, probs = seq(0,1,.25)),
                  include.lowest = TRUE, labels = c("[Q1]","[Q2]","[Q3]","[Q4]"))
  categorize_Yb2 <- cut(Yb2_cont,
                  breaks = quantile(Yb2_cont, probs = seq(0,1,1/3)),
                  include.lowest = TRUE, labels = c("[Q1]","[Q2]","[Q3]"))

  Yb1_true <- factor(sapply(Yb1_cont, categorize_Yb1), ordered = TRUE, levels = c("[0-25]", "[25-50]", "[50-75]", "[75+]"))
  Yb2_true <- factor(sapply(Yb1_cont, categorize_Yb2), ordered = TRUE, levels = c("A_Low", "B_Medium", "C_High"))

  df_X <- data.frame(X1=factor(X1), X2=factor(X2), X3=factor(X3), X4=factor(X4), X5=X5)
  mar <- produce_NA(df_X, mechanism=mech, perc.missing = perc_missing, by.patterns=FALSE)
  X_mar <- mar$data.incomp

  DB <- as.factor(ifelse(runif(num_samples) < 0.5, "A", "B"))
  Yb1 <- Yb1_true
  Yb2 <- Yb2_true
  Yb1[DB=="B"] <- NA
  Yb2[DB=="A"] <- NA

  df_working <- data.frame(DB=DB, Yb1=Yb1, Yb2=Yb2)
  df_working <- cbind(df_working, X_mar)

  mask_Yb1 <- is.na(df_working$Yb1)
  mask_Yb2 <- is.na(df_working$Yb2)

  df_impute <- df_working[, -1] %>% mutate_if(is.character, as.factor)

  return(list(
    df_working = df_working,
    df_impute  = df_impute,
    Yb1_true   = Yb1_true,
    Yb2_true   = Yb2_true,
    mask_Yb1   = mask_Yb1,
    mask_Yb2   = mask_Yb2
  ))
}

# ------------------------
# Core Simulation Function
# ------------------------
run_incomplete_case <- function(num_samples, perc_missing, mech) {

  # Generate data once for all methods except OTrecod
  dat <- generate_data(num_samples, perc_missing, mech)
  df_working <- dat$df_working
  df_impute  <- dat$df_impute
  Yb1_true   <- dat$Yb1_true
  Yb2_true   <- dat$Yb2_true
  mask_Yb1   <- dat$mask_Yb1
  mask_Yb2   <- dat$mask_Yb2

  iter_results <- list()
  na_metrics <- c(Precision = NA, Kappa = NA, MAE = NA)

  # --- 1. missForest ---
  capture.output(mf_res <- try(missForest(df_impute, maxiter = 5, ntree = 50), silent = TRUE))
  if(!inherits(mf_res, "try-error")) {
    res_mf_y1 <- calc_metrics(Yb1_true, mf_res$ximp$Yb1, mask_Yb1)
    res_mf_y2 <- calc_metrics(Yb2_true, mf_res$ximp$Yb2, mask_Yb2)
    iter_results[["missForest"]] <- (res_mf_y1 + res_mf_y2) / 2
  } else {
    iter_results[["missForest"]] <- na_metrics
  }

  # --- 2. missMDA (FAMD) ---
  capture.output(famd_res <- try(imputeFAMD(df_impute, ncp = 2), silent = TRUE))
  if(!inherits(famd_res, "try-error")) {
    res_famd_y1 <- calc_metrics(Yb1_true, famd_res$completeObs$Yb1, mask_Yb1)
    res_famd_y2 <- calc_metrics(Yb2_true, famd_res$completeObs$Yb2, mask_Yb2)
    iter_results[["missMDA"]] <- (res_famd_y1 + res_famd_y2) / 2
  } else {
    iter_results[["missMDA"]] <- na_metrics
  }

  # --- 3. kNN ---
  capture.output(df_knn <- try(kNN(df_impute, k = 5, imp_var = FALSE), silent = TRUE))
  if(!inherits(df_knn, "try-error")) {
    res_knn_y1 <- calc_metrics(Yb1_true, df_knn$Yb1, mask_Yb1)
    res_knn_y2 <- calc_metrics(Yb2_true, df_knn$Yb2, mask_Yb2)
    iter_results[["kNN"]] <- (res_knn_y1 + res_knn_y2) / 2
  } else {
    iter_results[["kNN"]] <- na_metrics
  }

  # --- 4. MICE ---
  capture.output(df_mice <- try({
    init_m <- mice(df_impute, maxit=0)
    meth_m <- init_m$method
    meth_m[c("Yb1", "Yb2")] <- "rf"
    complete(mice(df_impute, method=meth_m, m=1, maxit=5, printFlag=FALSE))
  }, silent = TRUE))
  if(!inherits(df_mice, "try-error")) {
    res_mice_y1 <- calc_metrics(Yb1_true, df_mice$Yb1, mask_Yb1)
    res_mice_y2 <- calc_metrics(Yb2_true, df_mice$Yb2, mask_Yb2)
    iter_results[["MICE"]] <- (res_mice_y1 + res_mice_y2) / 2
  } else {
    iter_results[["MICE"]] <- na_metrics
  }

  # --- 5. OTrecod (with retry loop on fresh data each attempt) ---
  ot_success <- FALSE
  for (ot_try in 1:OT_MAX_RETRIES) {
    # Regenerate fresh data on each retry so OTrecod gets a new draw
    dat_ot <- generate_data(num_samples, perc_missing, mech)

    capture.output(R_OUTC1 <- try(
      OT_joint(
        dat_ot$df_working, prox.X=0.10, convert.num=8, convert.class=1,
        nominal=c(1,4:7), ordinal=2:3, dist.choice="H",
        maxrelax=0.1, which.DB="BOTH"
      ), silent = TRUE))

    if(!inherits(R_OUTC1, "try-error")) {
      df_ot <- dat_ot$df_working
      df_ot$Yb2[dat_ot$mask_Yb2] <- R_OUTC1$DATA1_OT[, "OTpred"]
      df_ot$Yb1[dat_ot$mask_Yb1] <- R_OUTC1$DATA2_OT[, "OTpred"]
      res_ot_y1 <- calc_metrics(dat_ot$Yb1_true, df_ot$Yb1, dat_ot$mask_Yb1)
      res_ot_y2 <- calc_metrics(dat_ot$Yb2_true, df_ot$Yb2, dat_ot$mask_Yb2)
      iter_results[["OTrecod"]] <- (res_ot_y1 + res_ot_y2) / 2
      ot_success <- TRUE
      if(ot_try > 1) cat(sprintf("  OTrecod succeeded on retry %d\n", ot_try))
      break
    }
  }

  if(!ot_success) {
    cat(sprintf("  OTrecod failed after %d retries — recording NA\n", OT_MAX_RETRIES))
    iter_results[["OTrecod"]] <- na_metrics
  }

  # Compile and return
  out_df <- data.frame(
    Method    = names(iter_results),
    Precision = sapply(iter_results, function(x) x["Precision"]),
    Kappa     = sapply(iter_results, function(x) x["Kappa"]),
    MAE       = sapply(iter_results, function(x) x["MAE"]),
    stringsAsFactors = FALSE
  )

  return(out_df)
}

# ------------------------
# Execute Simulation Wrapper
# ------------------------
results_list <- list()

for (seed in seed_list) {
  set.seed(seed)

  for (n in num_samples_list) {
    for (p_miss in perc_missing_list) {
      for (mech in mechanism_list) {
        cat(sprintf("Running seed=%d, n=%d, missing=%.1f, mech=%s...\n", seed, n, p_miss, mech))

        for (rep_i in 1:n_reps) {
          sim_data <- tryCatch({
            run_incomplete_case(n, p_miss, mech)
          }, error = function(e){
            cat(sprintf("  Hard error in rep %d: %s\n", rep_i, e$message))
            NULL
          })

          if(!is.null(sim_data)) {
            sim_data$Seed <- seed
            sim_data$n <- n
            sim_data$perc_missing <- p_miss
            sim_data$mechanism <- mech
            results_list[[length(results_list) + 1]] <- sim_data
          }
        }
      }
    }
  }
}

final_raw_data <- bind_rows(results_list)

summary_results <- final_raw_data %>%
  group_by(n, perc_missing, mechanism, Method) %>%
  summarize(
    total_successful_runs = sum(!is.na(Precision)),
    mean_precision = mean(Precision, na.rm = TRUE),
    err_precision  = 1.96 * (sd(Precision, na.rm = TRUE) / sqrt(sum(!is.na(Precision)))),
    mean_kappa     = mean(Kappa, na.rm = TRUE),
    err_kappa      = 1.96 * (sd(Kappa, na.rm = TRUE) / sqrt(sum(!is.na(Kappa)))),
    mean_mae       = mean(MAE, na.rm = TRUE),
    err_mae        = 1.96 * (sd(MAE, na.rm = TRUE) / sqrt(sum(!is.na(MAE)))),
    .groups = 'drop'
  )

print(summary_results)
write.csv(summary_results, "../incomplete_case_all_4.csv", row.names = FALSE)
