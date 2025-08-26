library(glmnet)
library(pROC)

# Prepare data

df <- read.table("example_individual_snpScore_4regression.dataframe")

#X <- as.matrix(df[, grep("rs", colnames(df))])
X <- as.matrix(df[, 3:ncol(df)])
y <- as.factor(df$caseOrControl)

## Fit LASSO (L1 penalty) logistic regression
#cv_fit <- cv.glmnet(X, y, family = "binomial", alpha = 1)### alpha = 0.5, elastic net


set.seed(123)  # for reproducibility

# X: matrix of predictors (samples x features)
# y: binary response vector (0/1)

n_folds <- 10
folds <- sample(rep(1:n_folds, length.out = length(y)))

aucs <- c()

for (i in 1:n_folds) {
  # Split data
  train_idx <- which(folds != i)
  test_idx  <- which(folds == i)

  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test  <- X[test_idx, ]
  y_test  <- y[test_idx]

  # Train model on training set
  cv_fit <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 1)### alpha = 1, lasso regression

  # Predict on test set
  probs <- predict(cv_fit, newx = X_test, s = "lambda.min", type = "response")

  # Compute AUROC
  roc_obj <- roc(y_test, as.numeric(probs))
  aucs[i] <- auc(roc_obj)
}

# Mean AUROC across folds
mean_auc <- mean(aucs)
cat("Mean AUROC across 10 folds:", mean_auc, "\n")

## Coefficients at optimal lambda
#coef(cv_fit, s = "lambda.min")

# Final model on the full dataset
final_cv_fit <- cv.glmnet(X, y, family = "binomial", alpha = 1)

# Extract best lambda
best_lambda <- final_cv_fit$lambda.min

# Extract coefficients at best lambda
coefs <- coef(final_cv_fit, s = best_lambda)

# Convert to matrix and filter non-zero coefficients (excluding intercept)
nonzero_coefs <- as.matrix(coefs)[which(coefs != 0), , drop = FALSE]

# Optional: Remove intercept if not interested
nonzero_coefs <- nonzero_coefs[rownames(nonzero_coefs) != "(Intercept)", , drop = FALSE]

# Print selected SNPs and their coefficients
print(nonzero_coefs)