get_mse <- function(l) {
  l$performance_metrics
}

helmert_coding <- function(n) {
  v <- c()
  for (i in seq(from = 1, to = n-1)) {
    v_col <- c()
    if (i == 1){
      v_col <- c(-(n-1)/n, rep(1/n, n-1))
    } else {
      v_col <- rep(0, i - 1)
      v_col <- c(v_col, -(n - i) / (n - i + 1))
      v_col <- c(v_col, rep(1/(n - i + 1), n - i))
    }
    v <- c(v, v_col)
  }
  matrix(v, nrow = n)
}

sequential_coding <- function(n) {
  v <- c()
  for (i in seq(from = 1, to = n-1)) {
    v_col <- c()
    v_col <- rep(0, i)
    v_col <- c(v_col, rep(1, n - i))
    v <- c(v, v_col)
  }
  matrix(v, nrow = n)
}

reference_mat_entry <- function(strategy, cols_tobeCoded, col_index) {
  col <- cols_tobeCoded[col_index]
  if ((strategy == 1) || (strategy == 2)) {
    seq(length(levels(train_df[, col])))
  } else {
    permn(seq(length(levels(train_df[, col]))))
  }
}

ref_encoding <- function(train_df, 
                         cols_tobeCoded, 
                         other_predictors = c(), 
                         y_train, 
                         strategy_group, 
                         ref_group, 
                         intercept, 
                         fit_family) {
  

  strategy_ref <- c(contr.treatment, contr.sum, helmert_coding, sequential_coding)
  coding_strategy <- c()
  for (i in strategy_group) {
    coding_strategy <- c(coding_strategy, strategy_ref[i])
  }
  
  x_train <- c()
  for (col_index in seq(length(cols_tobeCoded))) {
    col <- cols_tobeCoded[col_index]
    formula <- paste("~", col, sep = "")
    mylist <- list()
    if ((strategy_group[col_index] == 1) || (strategy_group[col_index] == 2)) {
      mylist[[col]] <- coding_strategy[[col_index]](length(levels(train_df[,col])))
      encoded_col <- model.matrix(as.formula(formula), contrasts.arg = mylist, data = train_df)
      x_train <- cbind(x_train, encoded_col[, -ref_group[[col_index]]])
    } else {
      mylist[[col]] <- coding_strategy[[col_index]](length(levels(train_df[,col])))
      train_df_reorder <- train_df
      train_df_reorder[, col] <- factor(train_df[, col], levels(train_df[, col])[unlist(ref_group[[col_index]])])
      encoded_col <- model.matrix(as.formula(formula), contrasts.arg = mylist, data = train_df_reorder)
      x_train <- cbind(x_train, encoded_col[, -1])
    }
  }
  
  for (j in other_predictors) {
    x_train <- cbind(x_train, train_df[,j])
  }
  
  cv <- cv.glmnet(x = x_train, 
                  y = y_train, 
                  alpha = 1, 
                  nlambda = 1000)
  model <- glmnet(x_train, 
                  y_train, 
                  alpha = 1, 
                  lambda = cv$lambda.min, 
                  intercept = intercept, 
                  family = fit_family)
  if (fit_family == "binomial" | fit_family == "multinomial") {
    pre <- predict(model, x_train, s = cv$lambda.min, type = "class")
    return(list("performance_metrics" = auc(y_train, pre),
               "coef" = model))
  } else {
    pre <- predict(model, x_train, s = cv$lambda.min)
    return(list("performance_metrics" = mean ((y_train - pre)^2), 
                "coef" = model))
  }
  
}


ref_encoding_wrapper <- function(train_df, cols_tobeCoded, other_predictors, strategy_group, y_train, intercept, fit_family) {
  reference_mat_list <- list()
  for (col_index in seq(length(cols_tobeCoded))) {
    col <- cols_tobeCoded[col_index]
    reference_mat_list <- c(reference_mat_list, 
                            list(reference_mat_entry(strategy_group[col_index], cols_tobeCoded, col_index)))
  }
  reference_mat <- expand.grid(reference_mat_list)
  
  mse_ref <- future_apply(reference_mat, 1, ref_encoding, 
                          train_df = train_df,
                          cols_tobeCoded = cols_tobeCoded,
                          other_predictors = other_predictors,
                          y_train = y_train,
                          strategy_group = strategy_group,
                          intercept = intercept,
                          fit_family = fit_family
  )
  if (fit_family == "gaussian") {
    mse_result <- future_sapply(mse_ref, get_mse)
    min_mse_index <- which.min(mse_result)
    best_result <- mse_ref[[min_mse_index]]
    best_reference_group <- reference_mat[min_mse_index, ]
    list("best_reference_group" = unlist(best_reference_group), 
         "performance_metrics" = best_result$performance_metrics,
         "coef" = best_result$coef)
  }
  else {
    auc_result <- future_sapply(mse_ref, get_mse)
    max_auc_index <- which.max(auc_result)
    best_result <- mse_ref[[max_auc_index]]
    best_reference_group <- reference_mat[max_auc_index, ]
    list("best_reference_group" = unlist(best_reference_group), 
         "performance_metrics" = best_result$performance_metrics,
         "coef" = best_result$coef)
  }
}

training <- function(train_df, 
                     cols_tobeCoded, 
                     other_predictors = NULL, 
                     y_col, 
                     seed = 99999, 
                     standardize = FALSE,
                     intercept = TRUE,
                     fit_family = "gaussian") {

# Data Preprocessing (prepare for y_train & other_predictors, standardize (or not) other_predictors)  
  set.seed(seed)
  ## input: dataframe, 
  ##        columns to be coded, 
  ##        columns no need for coding (optional = NULL),
  ##        column to be predicted, 
  ##        seed (optional = 99999),
  ##        standardize (optional = FALSE)
  ##        intercept (optional = TRUE)
  ## output: 
  print("data preprocessing")
  for (col in cols_tobeCoded) {
    if (!(col %in% colnames(train_df))) {
      stop(paste("predictor variable", col, "not found in the training dataset"))
    }
  }
  
  for (col in other_predictors) {
    if (!(col %in% colnames(train_df))) {
      stop(paste("predictor variable", col, "not found in the training dataset"))
    }
  }
  
  if (!(y_col[1] %in% colnames(train_df))) {
    stop(paste("outcome variable", y_col[1], "not found in the training dataset"))
  }
  
  y_train <- data.matrix(train_df[y_col[1]])
  
  ## determine other predictors (if default, no other predictors)
  if (is.null(other_predictors)) {
    other_predictors <- c()
  }

  
  if (standardize) {
    for (j in other_predictors) {
      if (is.factor(train_df[, j])) {
        train_df[, j] <- as.numeric(train_df[, j])
        print(paste("Warning: Changing factor variable", j, "to numeric"))
      }
      train_df[, j] <- scale(train_df[, j])
    }
  }
  
  print("choosing for coding strategy and reference group")
  ## coding strategy matrix
  num_strategy <- rep(4, length(cols_tobeCoded))
  cum_prod_strategy <- cumprod(num_strategy)
  strategy_matrix <- matrix(0, nrow = cum_prod_strategy[length(cols_tobeCoded)], ncol = length(cols_tobeCoded))
  
  tmp <- cum_prod_strategy[length(cols_tobeCoded)]
  for (i in seq(length(num_strategy))) {
    if (i == 1) {
      strategy_matrix[, i] <- rep(seq(num_strategy[i]), rep(tmp / num_strategy[i], num_strategy[i]))
    } else {
      strategy_matrix[, i] <- rep(rep(seq(num_strategy[i]), rep(tmp / num_strategy[i], num_strategy[i])), cum_prod_strategy[i-1])
    }
    tmp <- tmp / num_strategy[i]
  }
  
  ## choose the best coding strategy and reference group combination
  result_enco <- future_apply(strategy_matrix, 1, ref_encoding_wrapper,
                              train_df = train_df,
                              cols_tobeCoded = cols_tobeCoded,
                              other_predictors = other_predictors,
                              y_train = y_train,
                              intercept = intercept,
                              fit_family = fit_family)
  
  print("preparing for output")
  
  if (fit_family == "gaussian") {
    mse_result <- future_sapply(result_enco, get_mse)
    min_mse_index <- which.min(mse_result)
    best_ref <- result_enco[[min_mse_index]]
    best_coding_strategy_index <- strategy_matrix[min_mse_index, ]
    coding_strategy <- c("Dummy", "Contrast", "Helmert", "Sequential")
    best_coding_strategy <- coding_strategy[best_coding_strategy_index]
  } else {
    auc_result <- future_sapply(result_enco, get_mse)
    max_auc_index <- which.max(auc_result)
    best_ref <- result_enco[[max_auc_index]]
    best_coding_strategy_index <- strategy_matrix[max_auc_index, ]
    coding_strategy <- c("Dummy", "Contrast", "Helmert", "Sequential")
    best_coding_strategy <- coding_strategy[best_coding_strategy_index]
  }
  
  
  reference_group <- list()
  count <- 1
  for (col_index in seq(length(cols_tobeCoded))) {
    col <- cols_tobeCoded[col_index]
    if (best_coding_strategy_index[col_index] <= 2) {
      reference_group <- c(reference_group, 
                           list(levels(train_df[, col])[count]))
      count <- count + 1
    } else {
      tmp <- best_ref[["best_reference_group"]][seq(count, count - 1 + length(levels(train_df[, col])))]
      reference_group <- c(reference_group, 
                           list(levels(train_df[, col])[tmp]))
      count <- count + length(levels(train_df[, col]))
    }
  }
  names(reference_group) <- cols_tobeCoded
  
  if (fit_family == "gaussian") {
    recommended_model <- list(coding_strategy = best_coding_strategy,
                              reference_group = reference_group,
                              mse = best_ref[["performance_metrics"]],
                              coef = best_ref[["coef"]])
    
    recommended_model_print <- list(coding_strategy = best_coding_strategy,
                                    reference_group = reference_group,
                                    mse = best_ref[["performance_metrics"]])
  } else {
    recommended_model <- list(coding_strategy = best_coding_strategy,
                              reference_group = reference_group,
                              auc = best_ref[["performance_metrics"]],
                              coef = best_ref[["coef"]])
    
    recommended_model_print <- list(coding_strategy = best_coding_strategy,
                                    reference_group = reference_group,
                                    auc = best_ref[["performance_metrics"]])
  }
  
  
  ## output
  result_input <- list("cols_Coded" = cols_tobeCoded, 
                       "other_predictors" = other_predictors, 
                       "y_col" = y_col,
                       "fit_family" = fit_family)
  result <- list("Input" = result_input, 
                 "Recommended_Model" = recommended_model_print)
  
  print(result)
  result$Recommended_Model$coef <- recommended_model$coef
  invisible(result)
}

ref_encoding_test <- function(train_df, 
                              cols_tobeCoded, 
                              other_predictors = c(), 
                              y_train, 
                              strategy_group, 
                              ref_group, 
                              coef, 
                              intercept,
                              fit_family) {
  
  strategy_ref <- c(contr.treatment, contr.sum, helmert_coding, sequential_coding)
  coding_strategy <- c()
  for (i in strategy_group) {
    coding_strategy <- c(coding_strategy, strategy_ref[i])
  }
  
  x_train <- c()
  for (col_index in seq(length(cols_tobeCoded))) {
    col <- cols_tobeCoded[col_index]
    formula <- paste("~", col, sep = "")
    mylist <- list()
    if ((strategy_group[col_index] == 1) || (strategy_group[col_index] == 2)) {
      mylist[[col]] <- coding_strategy[[col_index]](length(levels(train_df[,col])))
      encoded_col <- model.matrix(as.formula(formula), contrasts.arg = mylist, data = train_df)
      x_train <- cbind(x_train, encoded_col[, -ref_group[[col_index]]])
    } else {
      mylist[[col]] <- coding_strategy[[col_index]](length(levels(train_df[,col])))
      train_df_reorder <- train_df
      train_df_reorder[, col] <- factor(train_df[, col], levels(train_df[, col])[unlist(ref_group[[col_index]])])
      encoded_col <- model.matrix(as.formula(formula), contrasts.arg = mylist, data = train_df_reorder)
      x_train <- cbind(x_train, encoded_col[, -1])
    }
  }

  
  for (j in other_predictors) {
    x_train <- cbind(x_train, train_df[,j])
  }
  
  if (fit_family == "gaussian") {
    pre <- predict.glmnet(coef, x_train)
    return(mean((y_train - pre)^2))
  } else {
    x_train <- as.matrix(x_train)
    pre <- predict(coef, x_train, type = "class")
    return(auc(y_train, pre))
  }
}

testing <- function(test, 
                    result_from_training, 
                    ref_groups = NULL, 
                    strategy = NULL, 
                    standardize = FALSE) {
  
  y_test <- data.matrix(test[result_from_training$Input$y_col[1]])
  
  if (standardize) {
    for (j in other_predictors) {
      test[, j] <- scale(test[, j])
    }
  }
  
  coding_strategy <- c()
  if (is.null(strategy)) {
    coding_strategy <- result_from_training$Recommended_Model$coding_strategy
  } else {
    coding_strategy <- strategy
  }
  
  if (is.null(coding_strategy)) {
    stop("No Coding Strategy Chosen")
  }
  
  ## convert coding strategy into code
  strategy <- c()
  for (i in coding_strategy) {
    if (i == "Dummy") {
      strategy <-  c(strategy, 1)
    } else if (i == "Contrast") {
      strategy <-  c(strategy, 2)
    } else if (i == "Helmert") {
      strategy <- c(strategy, 3)
    } else if (i == "Sequential") {
      strategy <- c(strategy, 4)
    }else {
      stop("Coding Strategy not Suitable")
      break
    }
  }
  
  cols_coded <- result_from_training$Input$cols_Coded
  ref_group_index <- list()
  for (i in seq(length(cols_coded))) {
    col <- cols_coded[i]
    if (strategy[i] == 1 | strategy[i] == 2) {
      index <- which(levels(test[, col]) == result_from_training$Recommended_Model$reference_group[[i]])
    } else {
      index <- match(levels(test[, col]), result_from_training$Recommended_Model$reference_group[[i]])
    }
    if (length(index) == 0) {
      stop("Reference Group not Found")
    }
    ref_group_index <- c(ref_group_index, list(index))
  }
  
  performance_metrics <- ref_encoding_test(train_df = test_df, 
                                           cols_tobeCoded = cols_coded, 
                                           other_predictors = result_from_training$Input$other_predictors, 
                                           y_train = y_test, 
                                           strategy_group = strategy, 
                                           ref_group = ref_group_index,
                                           coef = result_from_training$Recommended_Model$coef,
                                           fit_family = result_from_training$Input$fit_family)
  
  if (result_from_training$Input$fit_family == "gaussian") {
    list(mse = performance_metrics)
  } else {
    list(auc = performance_metrics)
  }
}