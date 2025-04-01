# - Benchmark example usage a Classification Problem -
#                   No Feature Selection
#---------------------------------------------------------
# Usage with ith resampled EnMap data in classifying Parthenium weed
# --------------------------------------------------------
# Library
library(ggplot2)
library(ggpubr)
library(grid)
library(ggExtra)
library(RColorBrewer)
library(grDevices)
library(GGally)
library(gridExtra)
library(corrplot)
library(psych)
#---------------------------------------------------------------
# Classification using Entire dataset (benchmark)
#---------------------------------------------------------------
# No Feature selection / Entire dataset task
# Load data
enmap_spectra <- read.csv("./Data/Classification-Species/Parthenium_weed_enmap.csv", header=TRUE, as.is=TRUE, check.names=FALSE)
head(enmap_spectra)
# Rename columns with valid names
colnames(enmap_spectra) <- make.names(colnames(enmap_spectra))
# Verify the new column names
colnames(enmap_spectra)

enmap_class_task <- mlr::makeClassifTask(data = enmap_spectra[,c(2:243) ], target = "Class", id="Entire Dataset")
enmap_class_task 

set.seed(153)
n = getTaskSize(enmap_class_task)
train.set = sample(n, size = round(2/3 * n))
test.set = setdiff(seq_len(n), train.set)

# -------------------Algorithm Settings----------------------------
#--------------------------------------------------------------
as.data.frame(mlr::listLearners())
# -------------------Generic Settings---------------------------
# Set Cross-validation strategy
rdesc = makeResampleDesc("CV", stratify = F, iters =10)
# ------------------- Random Forest Parameters----------------------------
# Do not other algorithms after running this section. Jump to parameter tuning. 
# Learner selection and parameter values
# getLearnerParamSet("classif.randomForest")
lrn <- makeLearner("classif.randomForest", predict.type = "response")
# set parameter space
# getParamSet("regr.randomForest")
params <- makeParamSet(makeIntegerParam("ntree",lower = 100L, upper = 1000L),
                       makeIntegerParam("mtry",lower = 1L,upper = 10L))
# Grid-search strategy
ctrl <- makeTuneControlGrid()

# -------------------Computing settings--------------------
# set parallel back-end
require(parallel)
require(parallelMap) 
parallelStartSocket(cpus = detectCores())
#---------------------------------------------------------------
# -------------------Parameter tuning---------------------------
set.seed(123)
system.time(tune.model <- tuneParams(learner = lrn, 
                                     task = enmap_class_task, 
                                     resampling = rdesc, 
                                     measures = list(mlr::acc), 
                                     par.set = params, 
                                     control = ctrl, 
                                     show.info = T))
# Print the optimal hyper-parameters
print(tune.model$x)
# Print CV metrics of hyper-parameters
print(tune.model$y)

# Set the optimal hyper-parameters
lrn_tune.model <- setHyperPars(lrn, par.vals = tune.model$x)
# -------------------Computing settings-----------------------------------------
# set parallel back-end
# require(parallel)
# require(parallelMap) 
parallelStartSocket(cpus = detectCores())
#-------------------------------------------------------------------------------
# -------------------Model Train and Evaluation---------------
# Train  the model 
# Returns an object of class Wrapped Model which encapsulates the fitted model, i.e., the output of the underlying R learning method.
model.split <- mlr::train(learner = lrn_tune.model, task = enmap_class_task, subset = train.set)
print(model.split)

# Extract variable importances (if the algorithm supports it)
getFeatureImportance(model.split)
# Extract feature importance
importance_df <- getFeatureImportance(model.split)$res

# Plot feature importance using ggplot2
var_imp <- ggplot(importance_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip the coordinates for better readability
  xlab("Features") +
  ylab("Importance") +
  ggtitle("") +
  theme_minimal()

# Save the plot as a TIFF at 600 DPI
ggsave("./Documentation/Results/var_imp_optispecpicker.tiff", plot = var_imp, 
       device = "tiff", 
       dpi = 600, 
       width = 8, 
       height = 6, 
       units = "in")

# Predict model
# Predicting the target values for new observations
model.split.pred <- predict(model.split, task = enmap_class_task, subset = test.set)
#-------------------------------------------------------------------------------
# Load necessary libraries
library(mlr)
library(caret)  # For confusion matrix plotting
library(reshape2)

# Assuming model.split.pred contains the predicted values
# Extract true labels and predicted labels
true_labels <- getPredictionTruth(model.split.pred)
predicted_labels <- getPredictionResponse(model.split.pred)

# Create confusion matrix
conf_matrix <- confusionMatrix(as.factor(predicted_labels), as.factor(true_labels))

# Assuming you have already created the confusion matrix 'conf_matrix'
# Convert the confusion matrix table to a dataframe for plotting
conf_matrix_df <- as.data.frame(conf_matrix$table)

# Rename the columns for better understanding
colnames(conf_matrix_df) <- c("True", "Predicted", "Samples")
#-------------------------------------------------------------------------------
# Plot confusion matrix as a heatmap using ggplot2
#-----------------------------------------
cm_plot <- ggplot(conf_matrix_df, aes(x = True, y = Predicted)) +
  geom_tile(aes(fill = Samples), color = "white") +
  scale_fill_gradient(low = "lightblue", high = "red") +
  geom_text(aes(label = Samples), vjust = 1) +
  labs(title = "OptiSpecPicker subset", x = "Reference Class", y = "Predicted Class") + theme_bw()+ theme(plot.title=element_text(size=20, face="bold"), axis.text.x=element_text(size=14, color = "black"), axis.text.y=element_text(size=14, color = "black"), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),legend.title = element_text(size=14, face="bold"), legend.text = element_text(size=14), legend.justification=c(0.5,0.5), text = element_text(family = "serif", face="bold"))

# Save the plot as a TIFF at 600 DPI
ggsave("./Documentation/Results/conf_matrix_entire.tiff", plot = cm_plot, 
       device = "tiff", 
       dpi = 600, 
       width = 8, 
       height = 6, 
       units = "in")

# For two class confusion matrix
# # Plot the confusion matrix
# fourfoldplot(conf_matrix$table, color = c("lightblue", "pink"), 
#              main = "Confusion Matrix", conf.level = 0, margin = 1)
#-------------------------------------------------------------------------------
# Calculate and display overall accuracy
overall_accuracy <- conf_matrix$overall['Accuracy']
print(paste("Overall Accuracy: ", round(overall_accuracy, 2)))

# Producer's Accuracy (Sensitivity): 
# True Positives / (True Positives + False Negatives)
producers_accuracy <- conf_matrix$byClass[, 'Sensitivity']
# Print Producer's Accuracy for each class
print("Producer's Accuracy for each class:")
print(round(producers_accuracy, 2))

# User's Accuracy (Positive Predictive Value): 
# True Positives / (True Positives + False Positives)
users_accuracy <- conf_matrix$byClass[, 'Pos Pred Value']
# Print User's Accuracy for each class
print("User's Accuracy for each class:")
print(round(users_accuracy, 2))

# F1-Score for each class
f1_scores <- conf_matrix$byClass[, 'F1']
print("F1-Score for each class:")
print(round(f1_scores, 2))

macro_f1 <- mean(f1_scores, na.rm = TRUE)
print(paste("Macro F1-Score: ", round(macro_f1, 2)))

# Kappa statistic
kappa_stat <- conf_matrix$overall['Kappa']
print(paste("Kappa Statistic: ", round(kappa_stat, 2)))

# Balanced Accuracy
balanced_accuracy <- conf_matrix$byClass[, 'Balanced Accuracy']
print("Balanced Accuracy for each class:")
print(round(balanced_accuracy, 2))

# Error rate
error_rate <- 1 - overall_accuracy
print(paste("Error Rate: ", round(error_rate, 2)))

# Mean balanced accuracy
mean_balanced_accuracy <- mean(balanced_accuracy, na.rm = TRUE)
print(paste("Mean Balanced Accuracy: ", round(mean_balanced_accuracy, 2)))

# Specificity for each class
specificity <- conf_matrix$byClass[, 'Specificity']
print("Specificity for each class:")
print(round(specificity, 2))

