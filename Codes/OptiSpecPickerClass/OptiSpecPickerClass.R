library(caret)
library(randomForest)
library(FSelectorRcpp)
library(ggplot2)

#' Perform stepwise feature selection for classification tasks
#'
#' @param X A data frame of predictor variables.
#' @param y A factor vector of the response variable.
#' @param patience An integer indicating the number of iterations to wait for improvement before stopping (default is 5).
#' @param decision_threshold A numeric value to decide if the second-best accuracy is close enough to the best accuracy. If NULL, this comparison is not made (default is NULL).
#'
#' @return A list containing:
#' \item{results}{A data frame with details of each iteration, including selected features, accuracy, training accuracy, validation accuracy, computation time, and removed variables.}
#' \item{final_data}{A data frame containing the final selected features and the response variable.}
#'
#' @examples
#' data(iris)
#' result <- OptiSpecPickerClass(iris[, -5], iris$Species)
#' print(result$results)

OptiSpecPickerClass <- function(X=NULL, y=NULL, patience = 5, decision_threshold = NULL) {
  feature_scores <- FSelectorRcpp::information_gain(y ~ ., data = X, equal = TRUE)
  ranked_features <- order(feature_scores$importance, decreasing = TRUE)
  feature_names <- names(X)[ranked_features]
  
  results <- data.frame(
    Iteration = integer(), 
    Features = character(), 
    Accuracy = double(), 
    Training_Accuracy = double(),
    Validation_Accuracy = double(),
    Computation_Time = double(), 
    Removed_Variables = character()
  )
  
  best_score <- -Inf
  second_best_score <- -Inf
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
    
    model <- caret::train(subset_X, y, method = "rf", 
                   trControl = trainControl(method = "cv", number = 5, 
                                            savePredictions = "final"),
                   metric = "Accuracy")
    
    # Calculate accuracy for training and validation
    training_accuracy <- sum(model$pred$obs == model$pred$pred) / nrow(model$pred)
    validation_accuracy <- max(model$results$Accuracy)
    
    score <- validation_accuracy
    computation_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    # Check for improvement
    if (score > best_score) {
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
      Accuracy = score,
      Training_Accuracy = training_accuracy,
      Validation_Accuracy = validation_accuracy,
      Computation_Time = computation_time,
      Removed_Variables = removed_variables_iteration
    ))
    
    # Print selected features and their accuracy for each CV fold
    cat("Iteration:", i, "\n")
    cat("Selected features:", blend, "\n")
    cat("Training Accuracy:", training_accuracy, "\n")
    cat("Validation Accuracy:", validation_accuracy, "\n")
    cat("\n")
    
    if (no_improvement_count >= patience && is.null(stop_iteration)) {
      stop_iteration <- i - patience
    }
  }
  
  # Check if the second best accuracy is within the decision threshold, if provided
  if (!is.null(second_best_subset) && !is.null(decision_threshold)) {
    accuracy_difference <- second_best_score - best_score
    if (accuracy_difference <= decision_threshold) {
      best_score <- second_best_score
      best_subset <- second_best_subset
      best_iteration <- second_best_iteration
    }
  }
  
  # Plot Accuracy and Computation Time over iterations
  max_accuracy <- max(results$Accuracy, na.rm = TRUE)
  
  p <- ggplot(results, aes(x = Iteration)) +
    geom_line(aes(y = Accuracy), color = "blue") +
    geom_point(aes(y = Accuracy), color = "blue") +
    geom_line(aes(y = Computation_Time * (max_accuracy / max(Computation_Time, na.rm = TRUE))), color = "orange") +
    geom_point(aes(y = Computation_Time * (max_accuracy / max(Computation_Time, na.rm = TRUE))), color = "orange") +
    ggtitle("Accuracy and Computation Time over Iterations") +
    xlab("Iteration") +
    ylab("Accuracy") +
    scale_y_continuous(name = "Accuracy",
                       limits = c(0, max_accuracy * 1.1),
                       sec.axis = sec_axis(~ . * (max(results$Computation_Time, na.rm = TRUE) / max_accuracy), name = "Computation Time (s)")) +
    geom_vline(xintercept = best_iteration, linetype = "dashed", color = "blue") +
    annotate("text", x = best_iteration, y = max_accuracy, 
             label = paste("Optimal subset \nat iteration", best_iteration), 
             hjust = -0.1, vjust = 1, color = "black")
  
  if (!is.null(stop_iteration)) {
    p <- p + geom_vline(xintercept = stop_iteration, linetype = "dashed", color = "red") +
      annotate("text", x = stop_iteration, y = max_accuracy, 
               label = paste("Patience exceeded\nat iteration", stop_iteration), 
               hjust = -0.1, vjust = 1, color = "red")
    # Save the plot as a TIFF at 600 DPI
    ggsave("./Documentation/Results/comp_acc_opticlass.tiff", plot = p, 
           device = "tiff", 
           dpi = 600, 
           width = 8, 
           height = 6, 
           units = "in")
  }
  
  print(p)
  final_data <- cbind(X[, best_subset, drop = FALSE], y = y)
  results$Best_Subset <- paste(best_subset, collapse = ", ")
  
  return(list(results = results, final_data = final_data))
}
