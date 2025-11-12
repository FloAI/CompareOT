source('https://raw.githubusercontent.com/R-miss-tastic/website/master/static/how-to/generate/amputation.R')
set.seed(42)

num_samples_list <- c(100, 200, 300, 400, 500, 1000)
perc_missing_list <- c(0.1, 0.5, 0.9)
short_timeout_iterations <- c(31)

# Function to calculate precision with generic X variables
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

  Yb1_cont <- 20 + 5*X1_num + 10*X2_num + 2*X3_num - 8*X4_num + 0.3*X5_num + rnorm(num_samples, 0, 5)
  
  # Categorize outcomes
  categorize_Yb1 <- function(value) {
    if (value < 25) return("[0-25]")
    else if (value < 50) return("[25-50]")
    else if (value < 75) return("[50-75]")
    else return("[75+]")
  }
  categorize_Yb2 <- function(value) {
    if (value < 30) return("Low")
    else if (value < 60) return("Medium")
    else return("High")
  }

  Yb1 <- sapply(Yb1_cont, categorize_Yb1)
  Yb2 <- sapply(Yb1_cont, categorize_Yb2)
  df_true <- data.frame(Yb1, Yb2)

  # ------------------------
  # Prepare data frame for OTrecod
  # ------------------------
  df <- data.frame(X1=factor(X1), X2=factor(X2), X3=factor(X3), X4=factor(X4), X5=X5)
  
  # Introduce missingness
  missing_mask <- rank(runif(num_samples)) <= num_samples * perc_missing
  DB <- ifelse(missing_mask,"A","B")
  Yb1[DB=="B"] <- NA
  Yb2[DB=="A"] <- NA
  
  df <- cbind(DB, Yb1, Yb2, df)
  df <- df[order(df$DB), ]
  df$DB <- as.factor(df$DB)
  df$Yb1 <- as.factor(df$Yb1)
  df$Yb2 <- as.factor(df$Yb2)

  # ------------------------
  # mice recoding
  # ------------------------
  init <- mice(df, maxit=0)
  init
  meth <- init$method
  predM <- init$predictorMatrix

  meth[c("Yb1")]="rf"
  meth[c("Yb2")]="rf"

  meth
  imputed <- mice(df, method=meth, predictorMatrix=predM, m=5)
  imputation <- complete(imputed)
  imputated_df<- imputation[, c("Yb1", "Yb2")]
  imputated_df <- imputated_df[order(as.numeric(rownames(imputated_df))),]
  # ------------------------
  # Compute precision
  # ------------------------
  df2 <- df[, c("Yb1","Yb2")]
  m_mask <- is.na(df2)
  
  calculate_precision <- function(actual, imputed, m_mask) {
    true_positives <- sum(actual[m_mask]==imputed[m_mask], na.rm=TRUE)
    false_positives <- sum(!is.na(imputed[m_mask]) & actual[m_mask]!=imputed[m_mask])
    precision <- true_positives / (true_positives + false_positives)
    return(precision)
  }

  precision <- calculate_precision(df_true, imputated_df, m_mask)
  return(precision)
}

# ------------------------
# Run simulation over combinations
# ------------------------
results <- data.frame(num_samples=integer(),
                      perc_missing=numeric(),
                      mean_precision=numeric(),
                      margin_error=numeric())

for (num_samples in num_samples_list) {
  for (perc_missing in perc_missing_list) {
    precisions <- numeric()
    
    for (i in 1:10) {
      precision <- tryCatch({
        calculate_precision_for_combinations(num_samples, perc_missing)
      }, error=function(e){NA})
      precisions <- c(precisions, precision)
    }
    
    precisions_clean <- na.omit(precisions)
    
    if (length(precisions_clean) > 1) {
      mean_prec <- mean(precisions_clean)
      se <- sd(precisions_clean)/sqrt(length(precisions_clean))
      margin_err <- 1.96 * se
    } else {
      mean_prec <- NA
      margin_err <- NA
    }
    
    results <- rbind(results, data.frame(num_samples=num_samples,
                                         perc_missing=perc_missing,
                                         mean_precision=mean_prec,
                                         margin_error=margin_err))
  }
}

print(results)

