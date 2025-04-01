library(caret)
library(randomForest)
library(FSelectorRcpp)
library(ggplot2)

OptiSpecPickerRegr <- function(X, y=NULL, patience = 5, decision_threshold = NULL) {
  feature_scores <- FSelectorRcpp::information_gain(y ~ ., data = X, equal = TRUE)
  ranked_features <- order(feature_scores$importance, decreasing = TRUE)
  feature_names <- names(X)[ranked_features]
  
  results <- data.frame(
    Iteration = integer(), 
    Features = character(), 
    RMSE = double(), 
    Training_RMSE = double(),
    Validation_RMSE = double(),
    Computation_Time = double(), 
    Removed_Variables = character()
  )
  
  best_score <- Inf
  second_best_score <- Inf
  no_improvement_count <- 0
  stop_iteration <- NULL
  useful_features <- c()
  removed_features <- c()
  best_subset <- NULL
  second_best_subset <- NULL
  best_iteration <- 0
  second_best_iteration <- 0
  
  for (i in 1:length(ranked_features)) {
    start_time <- Sys.time()  # Start timer
    feature_to_add <- feature_names[i]
    blend <- c(useful_features, feature_to_add)
    subset_X <- X[, blend, drop = FALSE]
    
    model <- train(subset_X, y, method = "rf", 
                   trControl = trainControl(method = "cv", number = 5, 
                                            savePredictions = "final"),
                   metric = "RMSE")
    
    # Calculate RMSE for training and validation
    training_rmse <- sqrt(mean((model$pred$obs - model$pred$pred)^2))
    validation_rmse <- min(model$results$RMSE)
    
    score <- validation_rmse
    computation_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    # Check for improvement
    if (score < best_score) {
      second_best_score <- best_score
      second_best_subset <- best_subset
      second_best_iteration <- best_iteration
      
      best_score <- score
      best_iteration <- i
      best_subset <- blend
      no_improvement_count <- 0
      useful_features <- blend
      removed_variables_iteration <- ""
    } else {
      no_improvement_count <- no_improvement_count + 1
      removed_features <- c(removed_features, feature_to_add)
      removed_variables_iteration <- feature_to_add
    }
    
    results <- rbind(results, data.frame(
      Iteration = i, 
      Features = paste(blend, collapse = ", "),
      RMSE = score,
      Training_RMSE = training_rmse,
      Validation_RMSE = validation_rmse,
      Computation_Time = computation_time,
      Removed_Variables = removed_variables_iteration
    ))
    
    # Print selected features and their RMSE for each CV fold
    cat("Iteration:", i, "\n")
    cat("Selected features:", blend, "\n")
    cat("Training RMSE:", training_rmse, "\n")
    cat("Validation RMSE:", validation_rmse, "\n")
    cat("\n")
    
    if (no_improvement_count >= patience && is.null(stop_iteration)) {
      stop_iteration <- i - patience
    }
  }
  
  # Check if the second best RMSE is within the decision threshold, if provided
  if (!is.null(second_best_subset) && !is.null(decision_threshold)) {
    rmse_difference <- second_best_score - best_score
    if (rmse_difference <= decision_threshold) {
      best_score <- second_best_score
      best_subset <- second_best_subset
      best_iteration <- second_best_iteration
    }
  }
  
  # Plot RMSE and Computation Time over iterations
  max_rmse <- max(results$RMSE, na.rm = TRUE)
  
  p <- ggplot(results, aes(x = Iteration)) +
    geom_line(aes(y = RMSE), color = "blue") +
    geom_point(aes(y = RMSE), color = "blue") +
    geom_line(aes(y = Computation_Time * (max_rmse / max(Computation_Time, na.rm = TRUE))), color = "orange") +
    geom_point(aes(y = Computation_Time * (max_rmse / max(Computation_Time, na.rm = TRUE))), color = "orange") +
    ggtitle("RMSE and Computation Time over Iterations") +
    xlab("Iteration") +
    ylab("RMSEcv") +
    scale_y_continuous(name = "RMSE",
                       limits = c(0, max_rmse * 1.1),
                       sec.axis = sec_axis(~ . * (max(results$Computation_Time, na.rm = TRUE) / max_rmse), name = "Computation Time (s)")) +
    geom_vline(xintercept = best_iteration, linetype = "dashed", color = "blue") +
    annotate("text", x = best_iteration, y = max_rmse, 
             label = paste("Optimal subset \nat iteration", best_iteration), 
             hjust = -0.1, vjust = 1, color = "blue")
  
  if (!is.null(stop_iteration)) {
    p <- p + geom_vline(xintercept = stop_iteration, linetype = "dashed", color = "red") +
      annotate("text", x = stop_iteration, y = max_rmse, 
               label = paste("Patience exceeded\nat iteration", stop_iteration), 
               hjust = -0.1, vjust = 1, color = "red")
  }
  
  print(p)
  final_data <- cbind(X[, best_subset, drop = FALSE], y = y)
  results$Best_Subset <- paste(best_subset, collapse = ", ")
  
  return(list(results = results, final_data = final_data))
}
