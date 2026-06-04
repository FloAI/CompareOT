suppressPackageStartupMessages({
  library(MASS)
  library(VIM)
  library(mice)
  library(OTrecod)
  library(missForest)
  library(missMDA)
  library(dplyr)
  library(randomForest)
  library(pROC)
})

# ================================
# NCDS Dataset Setup & Imputation
# ================================

# Set seed for reproducibility
set.seed(42)

# ------------------------
# 1. Load & Prepare Datasets
# ------------------------
data(ncds_14)
data(ncds_5)

# Ensure datasets are complete (no missing values initially)
ncds_14_clean <- ncds_14[complete.cases(ncds_14), ]
ncds_5_clean  <- ncds_5[complete.cases(ncds_5), ]

# Merge the two datasets
merged_tab <- merge_dbs(
  ncds_14_clean, ncds_5_clean,
  row_ID1 = 1, row_ID2 = 1,
  NAME_Y = "GO90", NAME_Z = "RG91",
  ordinal_DB1 = 3, ordinal_DB2 = 4,
  seed_choice = 3023
)

# Extract the ready database (removing redundant identifier in col 4 if applicable)
merged_fin <- merged_tab$DB_READY[, -4]
merged_finish <- merged_tab$DB_READY[, -4]

# Ensure Y and Z are factors for classification models
merged_finish$Y <- as.factor(merged_finish$Y)
merged_finish$Z <- as.factor(merged_finish$Z)

# ------------------------
# 2. Imputation Methods
# ------------------------
imputed_datasets <- list()

# --- A. Optimal Transport (OTrecod) ---
cat("\nRunning OTrecod...\n")
capture.output(outj1 <- try(OT_joint(
  merged_fin,
  nominal = c(1:4),    # Specify nominal columns
  ordinal = 5:6,       # Specify ordinal columns
  dist.choice = "E",   # Euclidean distance
  which.DB = "BOTH"
), silent = TRUE))

df_OT <- merged_fin
df_OT$Z[is.na(df_OT$Z)] <- outj1$DATA1_OT[, "OTpred"]
df_OT$Y[is.na(df_OT$Y)] <- outj1$DATA2_OT[, "OTpred"]
df_OT$Y <- as.factor(df_OT$Y)
df_OT$Z <- as.factor(df_OT$Z)
imputed_datasets[["OTrecod"]] <- df_OT

# --- B. kNN ---
cat("Running kNN...\n")
capture.output(df_kNN_raw <- kNN(merged_finish, k = 5, imp_var = FALSE))
imputed_datasets[["kNN"]] <- df_kNN_raw

# --- C. MICE ---
cat("Running MICE...\n")
init <- mice(merged_finish, maxit = 0)
meth <- init$method
predM <- init$predictorMatrix
meth[c("Y", "Z")] <- "rf"
capture.output(imputed_mice <- mice(merged_finish, method = meth, predictorMatrix = predM, m = 1, maxit = 5))
df_MICE <- complete(imputed_mice)
imputed_datasets[["MICE"]] <- df_MICE

# --- D. missForest ---
cat("Running missForest...\n")
capture.output(mf_res <- missForest(merged_finish, maxiter = 5, ntree = 50))
imputed_datasets[["missForest"]] <- mf_res$ximp

# --- E. missMDA (FAMD) ---
cat("Running missMDA (FAMD)...\n")
# FAMD requires categorical variables to be strictly factors
df_for_famd <- merged_finish %>% mutate_if(is.character, as.factor)
capture.output(famd_res <- imputeFAMD(df_for_famd, ncp = 2))
df_missMDA <- famd_res$completeObs
df_missMDA$Y <- as.factor(df_missMDA$Y)
df_missMDA$Z <- as.factor(df_missMDA$Z)
imputed_datasets[["missMDA"]] <- df_missMDA

# ------------------------
# 3. Compute 1D Wasserstein Distance
# ------------------------
discrete_wass <- function(vec1, vec2) {
  # Strip whitespace and force character to ensure perfect alignment
  vec1 <- trimws(as.character(vec1))
  vec2 <- trimws(as.character(vec2))

  tab1 <- table(vec1) / length(vec1)
  tab2 <- table(vec2) / length(vec2)

  # Crucial: Sort the levels logically so cumsum() calculates distance correctly
  all_levels <- sort(unique(c(names(tab1), names(tab2))))

  p1 <- tab1[all_levels]; p1[is.na(p1)] <- 0
  p2 <- tab2[all_levels]; p2[is.na(p2)] <- 0

  sum(abs(cumsum(p1) - cumsum(p2)))
}

wass_results <- data.frame(Method = names(imputed_datasets), Wass_Y = NA, Wass_Z = NA)

for(i in 1:length(imputed_datasets)) {
  method_name <- names(imputed_datasets)[i]
  wass_results$Wass_Y[i] <- discrete_wass(imputed_datasets[["OTrecod"]]$Y, imputed_datasets[[method_name]]$Y)
  wass_results$Wass_Z[i] <- discrete_wass(imputed_datasets[["OTrecod"]]$Z, imputed_datasets[[method_name]]$Z)
}

cat("\n--- Wasserstein Distances (Compared to OTrecod) ---\n")
print(wass_results)

# ------------------------
# 4. Downstream Classification Evaluation
# ------------------------
evaluate_classification <- function(df, target_var) {

  # Force target to be an UNORDERED factor to prevent Ops.ordered error
  df[[target_var]] <- factor(as.character(df[[target_var]]))

  # Train/Test Split (70/30)
  train_idx <- sample(seq_len(nrow(df)), size = 0.7 * nrow(df))
  train_data <- df[train_idx, ]
  test_data <- df[-train_idx, ]

  # Train Random Forest
  formula_str <- paste(target_var, "~ .")
  rf_model <- randomForest(as.formula(formula_str), data = train_data, ntree = 100)

  # Predictions
  preds_class <- predict(rf_model, test_data)
  preds_prob <- predict(rf_model, test_data, type = "prob")

  # Force both to standard characters for safe comparison
  actual_char <- as.character(test_data[[target_var]])
  preds_char <- as.character(preds_class)

  # Accuracy
  accuracy <- sum(preds_char == actual_char) / length(actual_char)

  # Multiclass ROC AUC (handles multi-level factors automatically)
  auc_val <- multiclass.roc(test_data[[target_var]], preds_prob)$auc

  return(c(Accuracy = accuracy, ROC_AUC = as.numeric(auc_val)))
}

cat("\n--- Downstream Classification Performance (Predicting Y and Z) ---\n")
class_results_list <- list()

for(method in names(imputed_datasets)) {
  df_model <- imputed_datasets[[method]]

  # Drop DB column so the model doesn't just learn the missingness blocks
  if("DB" %in% colnames(df_model)) {
    df_model <- df_model %>% dplyr::select(-DB)
  }

  # Strip ALL ordered factors from the dataframe before modeling
  df_model <- df_model %>% mutate_if(is.ordered, ~ factor(as.character(.)))

  # Predict Y
  metrics_Y <- evaluate_classification(df_model, target_var = "Y")

  # Predict Z
  metrics_Z <- evaluate_classification(df_model, target_var = "Z")

  # Store both
  class_results_list[[method]] <- data.frame(
    Method = method,
    Target = c("Y", "Z"),
    Accuracy = c(metrics_Y["Accuracy"], metrics_Z["Accuracy"]),
    ROC_AUC = c(metrics_Y["ROC_AUC"], metrics_Z["ROC_AUC"])
  )
}

final_class_results <- bind_rows(class_results_list)
print(final_class_results)

# ------------------------
# 5. Save Results to CSV
# ------------------------
write.csv(wass_results, "../wasserstein_distances_ncds.csv", row.names = FALSE)
write.csv(final_class_results, "../classification_performance_ncds.csv", row.names = FALSE)

cat("\nFiles successfully saved: 'wasserstein_distances_ncds.csv' and 'classification_performance_ncds.csv'.\n")