
library(pROC)

# Prepare data

df <- read.table("example_individual_snpScore_4regression.dataframe")

X <- as.matrix(df[, grep("rs", colnames(df))])
#X <- as.matrix(df[, 3:ncol(df)])
y <- as.factor(df$caseOrControl)


set.seed(123)
n_folds <- 10
folds <- sample(rep(1:n_folds, length.out = length(y)))

aucs <- c()

for (i in 1:n_folds) {
  # Split data
  train_idx <- which(folds != i)
  test_idx <- which(folds == i)
  
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test <- X[test_idx, ]
  y_test <- y[test_idx]
  
  # Prepare data frames
  train_data <- data.frame(y = y_train, X_train)
  test_data <- data.frame(y = y_test, X_test)
  
  # Fit logistic regression
  fit <- glm(y ~ ., data = train_data, family = binomial)
  
  # Predict probabilities on test set
  probs <- predict(fit, newdata = test_data, type = "response")
  
  # Compute AUROC
  library(pROC)
  roc_obj <- roc(y_test, probs)
  aucs[i] <- auc(roc_obj)
}

cat("Mean AUROC across 10 folds:", mean(aucs), "\n")

# Finalize model on the entire dataset
full_data <- data.frame(y = y, X)
final_model <- glm(y ~ ., data = full_data, family = binomial)

# Extract coefficients
coefs <- coef(final_model)
print(coefs)
