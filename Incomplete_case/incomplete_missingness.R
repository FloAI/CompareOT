suppressPackageStartupMessages({
  library(MASS)
  library(VIM)
  library(missForest)
  library(missMDA)
  library(dplyr)
})

# Source the ampute function
source('https://raw.githubusercontent.com/R-miss-tastic/website/master/static/how-to/generate/amputation.R')

# --- Configuration ---
seed_list <- c(42, 123, 456, 789, 999) # 5 seeds
n_reps <- 6                            # 6 runs per seed
num_samples_list <- c(100, 200, 300, 400, 500, 1000)
perc_missing_list <- seq(0.1, 0.9, by = 0.4)
mechanism_list <- c("MCAR", "MAR", "MNAR") # Added mechanisms
target_R2 <- 1.0                           # Set to 1

# ------------------------
# Helper function for Precision
# ------------------------
calc_precision <- function(actual, imputed, mask) {
  act <- as.character(actual[mask])
  imp <- as.character(imputed[mask])
  return(sum(act == imp, na.rm=TRUE) / length(act))
}

# ------------------------
# Core Simulation Function
# ------------------------
run_incomplete_case <- function(num_samples, perc_missing, mech) {

  # 1. Generate covariates X1-X5
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

  # Handle R2 = 1 (Zero noise)
  if(target_R2 >= 1.0) {
    sigma <- 0
  } else {
    sigma <- sqrt(var(signal) * (1/target_R2 - 1))
  }

  Yb1_cont <- signal + sigma*rnorm(num_samples)

  categorize_Yb1 <- function(v) if(v<25) "[0-25]" else if(v<50) "[25-50]" else if(v<75) "[50-75]" else "[75+]"
  categorize_Yb2 <- function(v) if(v<30) "A_Low" else if(v<60) "B_Medium" else "C_High"

  Yb1_true <- factor(sapply(Yb1_cont, categorize_Yb1), ordered = TRUE, levels = c("[0-25]", "[25-50]", "[50-75]", "[75+]"))
  Yb2_true <- factor(sapply(Yb1_cont, categorize_Yb2), ordered = TRUE, levels = c("A_Low", "B_Medium", "C_High"))

  # 2. Introduce Missingness in X based on mechanism
  df_X <- data.frame(X1=factor(X1), X2=factor(X2), X3=factor(X3), X4=factor(X4), X5=X5)
  mar <- produce_NA(df_X, mechanism=mech, perc.missing = perc_missing, by.patterns=FALSE)
  X_mar <- mar$data.incomp # Left as incomplete

  # 3. Introduce Missingness in Y
  DB <- as.factor(ifelse(runif(num_samples) < 0.5, "A", "B"))
  Yb1 <- Yb1_true
  Yb2 <- Yb2_true
  Yb1[DB=="B"] <- NA
  Yb2[DB=="A"] <- NA

  df_working <- data.frame(DB=DB, Yb1=Yb1, Yb2=Yb2)
  df_working <- cbind(df_working, X_mar)

  mask_Yb1 <- is.na(df_working$Yb1)
  mask_Yb2 <- is.na(df_working$Yb2)

  iter_results <- list()

  # Ensure all characters are factors for FAMD and missForest
  df_impute <- df_working[, -1] %>% mutate_if(is.character, as.factor)

  # --- 1. missForest ---
  capture.output(mf_res <- try(missForest(df_impute, maxiter = 5, ntree = 50), silent = TRUE))
  if(!inherits(mf_res, "try-error")) {
    prec_mf_y1 <- calc_precision(Yb1_true, mf_res$ximp$Yb1, mask_Yb1)
    prec_mf_y2 <- calc_precision(Yb2_true, mf_res$ximp$Yb2, mask_Yb2)
    iter_results[["missForest"]] <- mean(c(prec_mf_y1, prec_mf_y2), na.rm=TRUE)
  } else {
    iter_results[["missForest"]] <- NA
  }

  # --- 2. missMDA (FAMD) ---
  capture.output(famd_res <- try(imputeFAMD(df_impute, ncp = 2), silent = TRUE))
  if(!inherits(famd_res, "try-error")) {
    prec_famd_y1 <- calc_precision(Yb1_true, famd_res$completeObs$Yb1, mask_Yb1)
    prec_famd_y2 <- calc_precision(Yb2_true, famd_res$completeObs$Yb2, mask_Yb2)
    iter_results[["missMDA"]] <- mean(c(prec_famd_y1, prec_famd_y2), na.rm=TRUE)
  } else {
    iter_results[["missMDA"]] <- NA
  }

  # Compile and return
  return(data.frame(
    Method = names(iter_results),
    Precision = unlist(iter_results),
    stringsAsFactors = FALSE
  ))
}

# ------------------------
# Execute Simulation Wrapper
# ------------------------
results_list <- list()

for (seed in seed_list) {
  set.seed(seed) # Set the new seed for this batch

  for (n in num_samples_list) {
    for (p_miss in perc_missing_list) {
      for (mech in mechanism_list) {
        cat(sprintf("Running seed=%d, n=%d, missing=%.1f, mech=%s...\n", seed, n, p_miss, mech))

        for (rep_i in 1:n_reps) { # 6 runs per seed
          sim_data <- tryCatch({
            run_incomplete_case(n, p_miss, mech)
          }, error = function(e){ NULL })

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

# Aggregate results across all seeds, runs, and mechanisms
summary_results <- final_raw_data %>%
  group_by(n, perc_missing, mechanism, Method) %>%
  summarize(
    total_successful_runs = sum(!is.na(Precision)),
    mean_precision = mean(Precision, na.rm = TRUE),
    err_precision = 1.96 * (sd(Precision, na.rm = TRUE) / sqrt(sum(!is.na(Precision)))),
    .groups = 'drop'
  )

print(summary_results)
write.csv(summary_results, "../incomplete_case_missforest_famd_all_mechs_R2_1.csv", row.names = FALSE)