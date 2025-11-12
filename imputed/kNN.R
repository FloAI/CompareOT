suppressPackageStartupMessages({
  require(MASS)
  require(norm)
  require(VIM)
  require(ggplot2)
  require(naniar)
})
library(mice)
library(OTrecod)
source('https://raw.githubusercontent.com/R-miss-tastic/website/master/static/how-to/generate/amputation.R')
set.seed(42)

num_samples_list <- c(500)
perc_missing_list <- c(0.5)
n_reps <- 30

calculate_precision_for_combinations <- function(num_samples, perc_missing) {

  # ------------------------
  # Generate covariates X1-X5
  # ------------------------
  X1 <- rbinom(num_samples, 1, 0.5)                   # binary
  X2_raw <- t(rmultinom(num_samples, 1, prob = c(0.3,0.4,0.3)))
  X2 <- apply(X2_raw, 1, function(x) c("A","B","Placebo")[which(x==1)])
  
  X3_raw <- t(rmultinom(num_samples, 1, prob = c(0.1,0.2,0.3,0.4)))
  X3 <- apply(X3_raw, 1, function(x) paste0("L", which(x==1)))
  
  X4 <- rbinom(num_samples, 1, 0.5)                   # binary
  X5 <- runif(num_samples, 20, 80)                   # continuous

  # ------------------------
  # Define continuous outcomes
  # ------------------------
  X1_num <- X1
  X2_num <- ifelse(X2=="A",1,ifelse(X2=="B",2,0))
  X3_num <- as.numeric(factor(X3))
  X4_num <- X4
  X5_num <- X5

  Yb1_cont <- 20 + 5*X1_num + 10*X2_num + 2*X3_num - 8*X4_num + 0.3*X5_num + rnorm(num_samples,0,5)
  Yb2_cont <- 15 + 4*X1_num + 8*X2_num + 1.5*X3_num - 6*X4_num + 0.25*X5_num + rnorm(num_samples,0,4)

  # Categorize outcomes
  categorize_Yb1 <- function(v) if(v<25) "[0-25]" else if(v<50) "[25-50]" else if(v<75) "[50-75]" else "[75+]"
  categorize_Yb2 <- function(v) if(v<30) "A_Low" else if(v<60) "B_Medium" else "C_High"

  Yb1 <- sapply(Yb1_cont, categorize_Yb1)
  Yb2 <- sapply(Yb2_cont, categorize_Yb2)
  df_true <- data.frame(Yb1,Yb2)

  # ------------------------
  # Prepare data frame for OTrecod
  # ------------------------
  df <- data.frame(X1=factor(X1), X2=factor(X2), X3=factor(X3), X4=factor(X4), X5=X5)

  # Introduce MAR missingness
  mar <- produce_NA(df, mechanism="MNAR", perc.missing = perc_missing, by.patterns=FALSE)
  X.mar <- mar$data.incomp

  # Impute missing covariates with mice
  init <- mice(X.mar, maxit=0)
  meth <- init$method
  predM <- init$predictorMatrix
  X.mar <- mice(X.mar, method="rf", predictorMatrix=predM, m=5)
  X.mar <- complete(X.mar)

  # Introduce missing values in Yb1/Yb2 using DB A/B
  missing_mask <- rank(runif(num_samples)) <= num_samples * 0.5
  DB <- ifelse(missing_mask,"A","B")
  Yb1[DB=="B"] <- NA
  Yb2[DB=="A"] <- NA

  df <- cbind(DB,Yb1,Yb2,X.mar)
  df <- df[order(df$DB),]
  df$DB <- as.factor(df$DB)
  df$Yb1 <- as.factor(df$Yb1)
  df$Yb2 <- as.factor(df$Yb2)

  # ------------------------
  # kNN recoding
  # ------------------------
  
  imputed_df <- kNN(df, k = 5, imp_var = FALSE)
  imputated_df<- imputed_df[, c("Yb1", "Yb2")]
  imputated_df <- imputated_df[order(as.numeric(rownames(imputated_df))),]

  imputated_df <- df[,c("Yb1","Yb2")]
  imputated_df <- imputated_df[order(as.numeric(rownames(imputated_df))),]

  # ------------------------
  # Compute precision
  # ------------------------
  m_mask <- is.na(df_true)
  tp <- sum(actual <- df_true[m_mask]==imputated_df[m_mask], na.rm=TRUE)
  fp <- sum(!is.na(imputated_df[m_mask]) & df_true[m_mask]!=imputated_df[m_mask])
  precision <- tp / (tp+fp)
  return(precision)
}

# ------------------------
# Run simulation for 30 repetitions
# ------------------------
results <- data.frame(num_samples=integer(), perc_missing=numeric(), mean_precision=numeric(), margin_error=numeric())

for (num_samples in num_samples_list) {
  for (perc_missing in perc_missing_list) {
    precisions <- numeric(n_reps)
    for (i in 1:n_reps) {
      precisions[i] <- tryCatch({
        calculate_precision_for_combinations(num_samples, perc_missing)
      }, error=function(e){NA})
    }

    precisions_clean <- na.omit(precisions)
    if(length(precisions_clean)>1){
      mean_prec <- mean(precisions_clean)
      se <- sd(precisions_clean)/sqrt(length(precisions_clean))
      margin_err <- 1.96*se
    } else {
      mean_prec <- NA
      margin_err <- NA
    }

    results <- rbind(results,data.frame(num_samples=num_samples, perc_missing=perc_missing,
                                        mean_precision=mean_prec, margin_error=margin_err))
  }
}

print(results)

