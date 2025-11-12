set.seed(42)

# Load example data
data(ncds_14)
data(ncds_5)

# Assuming ncds_14 and ncds_5 are loaded as data frames

# Merge the two datasets
merged_tab <- merge_dbs(ncds_14, ncds_5,
                        row_ID1 = 1, row_ID2 = 1,
                        NAME_Y = "GO90", NAME_Z = "RG91",
                        ordinal_DB1 = 3, ordinal_DB2 = 4,
                        impute = "MICE", R_MICE = 2,
                        seed_choice = 3023)


merged_fin = merged_tab$DB_READY[, -4]
merged_finish = merged_tab$DB_READY[, -4]
merged_finished =merged_tab$DB_READY[, -4]

outj1 <- OT_joint(
  merged_fin,
  nominal = c(1:4),
  ordinal = 5:6,
  dist.choice = "E",
  which.DB = "BOTH"
)

# Verify the output
verif_outj1 <- verif_OT(
  outj1,
  ordinal = FALSE,
  stab.prob = TRUE,
  min.neigb = 5
)

# Extract imputed data correctly from outj1
y_imputed <- outj1$DATA2_OT[, "OTpred"]
z_imputed <- outj1$DATA1_OT[, "OTpred"]
na_indices <- which(is.na(merged_fin$Z))
merged_fin$Z[na_indices] <- z_imputed
na_indices <- which(is.na(merged_fin$Y))
merged_fin$Y[na_indices] <- y_imputed

impu_df <- merged_fin[, c("Y", "Z")]
imp_df <- impu_df[order(as.numeric(rownames(impu_df))),]



imputed_df <- kNN(merged_finish, k = 5, imp_var = FALSE)
imputated_df<- imputed_df[, c("Y", "Z")]

init <- mice(merged_finished, maxit=0)
init
meth <- init$method
predM <- init$predictorMatrix

meth[c("Y")]="rf"
meth[c("Z")]="rf"

meth
imputed <- mice(merged_finished, method=meth, predictorMatrix=predM, m=5)
imputation <- complete(imputed)
imputated_df1<- imputation[, c("Y", "Z")]
imputated_df1 <- imputated_df1[order(as.numeric(rownames(imputated_df1))),]

# Load necessary library
library(dplyr)

# Function to compute discrete 1D Wasserstein distance for categorical variables
discrete_wass <- function(vec1, vec2) {
  # Compute probability distributions
  tab1 <- table(vec1) / length(vec1)
  tab2 <- table(vec2) / length(vec2)

  # Ensure same levels
  all_levels <- union(names(tab1), names(tab2))
  p1 <- tab1[all_levels]; p1[is.na(p1)] <- 0
  p2 <- tab2[all_levels]; p2[is.na(p2)] <- 0

  # Compute discrete Wasserstein: L1 distance between cumulative distributions
  sum(abs(cumsum(p1) - cumsum(p2)))
}

# Prepare vectors for each method and variable
Y_OT   <- impu_df$Y
Y_kNN  <- imputated_df$Y
Y_MICE <- imputated_df1$Y

Z_OT   <- impu_df$Z
Z_kNN  <- imputated_df$Z
Z_MICE <- imputated_df1$Z

# Compute pairwise distances for Y
wass_Y_OT_kNN  <- discrete_wass(Y_OT, Y_kNN)
wass_Y_OT_MICE <- discrete_wass(Y_OT, Y_MICE)
wass_Y_kNN_MICE <- discrete_wass(Y_kNN, Y_MICE)

wass_Y <- c(OT_kNN = wass_Y_OT_kNN, OT_MICE = wass_Y_OT_MICE, kNN_MICE = wass_Y_kNN_MICE)

# Compute pairwise distances for Z
wass_Z_OT_kNN  <- discrete_wass(Z_OT, Z_kNN)
wass_Z_OT_MICE <- discrete_wass(Z_OT, Z_MICE)
wass_Z_kNN_MICE <- discrete_wass(Z_kNN, Z_MICE)

wass_Z <- c(OT_kNN = wass_Z_OT_kNN, OT_MICE = wass_Z_OT_MICE, kNN_MICE = wass_Z_kNN_MICE)

# Display results
wass_Y
wass_Z


