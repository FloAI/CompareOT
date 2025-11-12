# ================================
#ncds dataset
# ================================

# Set seed for reproducibility
set.seed(42)

# ------------------------
# Load example datasets
# ------------------------
data(ncds_14)
data(ncds_5)

# Ensure datasets are complete (no missing values)
ncds_14_clean <- ncds_14[complete.cases(ncds_14), ]
ncds_5_clean  <- ncds_5[complete.cases(ncds_5), ]

# ------------------------
# Merge the two datasets
# ------------------------
merged_tab <- merge_dbs(
  ncds_14_clean, ncds_5_clean,
  row_ID1 = 1, row_ID2 = 1,
  NAME_Y = "GO90", NAME_Z = "RG91",
  ordinal_DB1 = 3, ordinal_DB2 = 4,
  seed_choice = 3023
)

# Keep relevant columns for further processing
merged_fin   <- merged_tab$DB_READY[, -4]
merged_finish <- merged_tab$DB_READY[, -4]

# ------------------------
# Perform Optimal Transport recoding
# ------------------------
outj1 <- OT_joint(
  merged_fin,
  nominal = c(1:4),    # Specify nominal columns
  ordinal = 5:6,       # Specify ordinal columns
  dist.choice = "E",   # Euclidean distance
  which.DB = "BOTH"
)

# Verify OT output stability
verif_outj1 <- verif_OT(
  outj1,
  ordinal = FALSE,
  stab.prob = TRUE,
  min.neigb = 5
)

# Extract imputed values and fill missing entries
y_imputed <- outj1$DATA2_OT[, "OTpred"]
z_imputed <- outj1$DATA1_OT[, "OTpred"]

merged_fin$Z[is.na(merged_fin$Z)] <- z_imputed
merged_fin$Y[is.na(merged_fin$Y)] <- y_imputed

impu_df <- merged_fin[, c("Y", "Z")]
imp_df <- impu_df[order(as.numeric(rownames(impu_df))), ]

# ------------------------
# kNN recoding
# ------------------------
imputed_df <- kNN(merged_finish, k = 5, imp_var = FALSE)
imputated_df <- imputed_df[, c("Y", "Z")]

# ------------------------
# MICE recoding
# ------------------------

# Initialize MICE
init <- mice(merged_finish, maxit = 0)
meth <- init$method
predM <- init$predictorMatrix

# Use random forest for Y and Z
meth[c("Y", "Z")] <- "rf"

# Perform MICE
imputed <- mice(merged_finish, method = meth, predictorMatrix = predM, m = 5)
imputation <- complete(imputed)

imputated_df1 <- imputation[, c("Y", "Z")]
imputated_df1 <- imputated_df1[order(as.numeric(rownames(imputated_df1))), ]

# ------------------------
# Compute Discrete 1D Wasserstein Distance
# ------------------------

discrete_wass <- function(vec1, vec2) {
  # Probability distributions
  tab1 <- table(vec1) / length(vec1)
  tab2 <- table(vec2) / length(vec2)

  # Align levels
  all_levels <- union(names(tab1), names(tab2))
  p1 <- tab1[all_levels]; p1[is.na(p1)] <- 0
  p2 <- tab2[all_levels]; p2[is.na(p2)] <- 0

  # L1 distance between cumulative distributions
  sum(abs(cumsum(p1) - cumsum(p2)))
}

# Prepare vectors for each imputation method
Y_OT   <- impu_df$Y
Y_kNN  <- imputated_df$Y
Y_MICE <- imputated_df1$Y

Z_OT   <- impu_df$Z
Z_kNN  <- imputated_df$Z
Z_MICE <- imputated_df1$Z

# Compute Wasserstein distances for Y
wass_Y_OT_kNN  <- discrete_wass(Y_OT, Y_kNN)
wass_Y_OT_MICE <- discrete_wass(Y_OT, Y_MICE)
wass_Y_kNN_MICE <- discrete_wass(Y_kNN, Y_MICE)

wass_Y <- c(OT_kNN = wass_Y_OT_kNN, OT_MICE = wass_Y_OT_MICE, kNN_MICE = wass_Y_kNN_MICE)

# Compute Wasserstein distances for Z
wass_Z_OT_kNN  <- discrete_wass(Z_OT, Z_kNN)
wass_Z_OT_MICE <- discrete_wass(Z_OT, Z_MICE)
wass_Z_kNN_MICE <- discrete_wass(Z_kNN, Z_MICE)

wass_Z <- c(OT_kNN = wass_Z_OT_kNN, OT_MICE = wass_Z_OT_MICE, kNN_MICE = wass_Z_kNN_MICE)

# ------------------------
# Display Results
# ------------------------
wass_Y
wass_Z

