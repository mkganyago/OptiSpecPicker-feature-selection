# - Benchmark example usage a Classification Problem - 
#         Recursive Feature Elimination
#--------------------------------------------------------------------------
# Usage with resampled EnMap data in classifying Parthenium weed
# -------------------------------------------------------------------------
# Load the data -----------------------------------------------------------
# data preprocessed in source("species_discrimination_example.R") script
enmap_spectra <- read.csv("./Data/Classification-Species/Parthenium_weed_enmap.csv", header=TRUE, as.is=TRUE, check.names=FALSE)

head(enmap_spectra)
head(enmap_spectra)

# Rename columns with valid names
colnames(enmap_spectra) <- make.names(colnames(enmap_spectra))
# Verify the new column names
colnames(enmap_spectra)

enmap_class_task <- mlr::makeClassifTask(data = enmap_spectra[,c(2:243) ], target = "Class", id="Entire Dataset")
enmap_class_task 
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

# #----------------------------------------------------------
# # Recursive Feature Elimination - Random Forest
# #----------------------------------------------------------
# # library(mlr)
# Specify the search strategy
ctrl.feat = makeFeatSelControlSequential(method = "sfs", alpha = 0.000001)

# Select features
rdesc.feat = makeResampleDesc("CV", iters = 10)

set.seed(186)
sfeats.rfe.class = selectFeatures(
  learner = "classif.randomForest", task = enmap_class_task, resampling = rdesc.feat, control = ctrl.feat,
  show.info = FALSE)

sfeats.rfe.class

analyzeFeatSelResult(sfeats.rfe.class)

#-------------------------------------------------------------------------------
# Required Libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Plot the selected features 
enmap_wavelengths <- as.numeric(names(dataset[,-243]))
enmap_spectra <- dataset
names(enmap_spectra)[names(enmap_spectra) == "y"] <- "Class"

# Average the spectra by class and per band
# Assuming enmap_spectra has 242 spectral band columns and one 'Class' column
# Also the spectral bands are in columns 2 to 242
average_spectra_by_class <- enmap_spectra %>%
  group_by(Class) %>%
  summarise(across(2:242, ~mean(.x, na.rm = TRUE)))  # Average reflectance per band by class

# Convert the wide data to long format for easier plotting
average_spectra_long <- average_spectra_by_class %>%
  pivot_longer(cols = 2:242, names_to = "Band", values_to = "Reflectance")

# Convert band names ("423.03", ..., "2438.6") to numeric wavelengths
average_spectra_long$Wavelength <- enmap_wavelengths[as.numeric(gsub("423.03", "", average_spectra_long$Band))]

# Remove rows with missing values
average_spectra_long <- average_spectra_long %>% filter(!is.na(Reflectance))

# Extract the selected bands from Best_Subset
best_subset <- c(563.95, 614.86, 740.4, 1179, 1441.3, 2275.3)
selected_bands <- best_subset

# Find the indices of the selected bands closest to EnMap wavelengths
selected_band_indices <- sapply(selected_bands, function(band) which.min(abs(enmap_wavelengths - band)))

# Plot the average spectra for each class with highlighted selected bands
aver_spectra_rf_rfe <- ggplot(average_spectra_long, aes(x = as.numeric(Band), y = Reflectance, color = Class, group = Class)) +  # Add group = Class
  geom_line(size = 1) +  # Plot average spectral curve for each class
  geom_vline(xintercept = enmap_wavelengths[selected_band_indices], linetype = "dashed", color = "red", size = 1) +  # Highlight selected bands
  geom_point(data = filter(average_spectra_long, Band %in% enmap_wavelengths[selected_band_indices]),
             aes(x = Wavelength, y = Reflectance), color = "red", size = 3) +  # Highlight selected bands as points
  xlab("Wavelength (nm)") +
  ylab("Reflectance (%)") +
  ggtitle("") +
  theme_minimal() +
  scale_color_manual(values = c("AT" = "blue", "GS" = "green", "OPS" = "orange", "PH" = "purple"))

print(aver_spectra_rf_rfe)
# Save the plot as a TIFF at 600 DPI
ggsave("./Documentation/Results/aver_spectra_rf_rfe_class.tiff", plot = aver_spectra_rf_rfe, 
       device = "tiff", 
       dpi = 600, 
       width = 8, 
       height = 6, 
       units = "in")

#-------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Test the subset's predictive capability with RF
# Create a new classif task with selected variables
enmap.subset.rfe <- subset(enmap_spectra, select = c(X563.95, X614.86, X740.4, X1179, X1441.3, X2275.3))
head(enmap.subset.rfe) 

# bind features with class
class = enmap_spectra$Class
enmap.subset.rfe <- cbind(enmap.subset.rfe, class)
head(enmap.subset.rfe)

enmap.subset.rfe.task <- makeClassifTask(id="rfe-rf-class", data = enmap.subset.rfe, target = "class")
enmap.subset.rfe.task

set.seed(153)
n = getTaskSize(enmap.subset.rfe.task)
train.set = sample(n, size = round(2/3 * n))
test.set = setdiff(seq_len(n), train.set)

# -------------------Algorithm Settings----------------------------------------
#------------------------------------------------------------------------------
as.data.frame(mlr::listLearners())
# -------------------Generic Settings------------------------------------------
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
set.seed(143)
system.time(tune.model <- tuneParams(learner = lrn, 
                                     task = enmap.subset.rfe.task, 
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
model.split <- mlr::train(learner = lrn_tune.model, task = enmap.subset.rfe.task, subset = train.set)
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
ggsave("./Documentation/Results/var_imp_rfe_class.tiff", plot = var_imp, 
       device = "tiff", 
       dpi = 600, 
       width = 8, 
       height = 6, 
       units = "in")

# Predict model
# Predicting the target values for new observations
model.split.pred <- predict(model.split, task = enmap.subset.rfe.task, subset = test.set)
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
  labs(title = "", x = "Reference Class", y = "Predicted Class") + theme_bw()+ theme(plot.title=element_text(size=20, face="bold"), axis.text.x=element_text(size=14, color = "black"), axis.text.y=element_text(size=14, color = "black"), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),legend.title = element_text(size=14, face="bold"), legend.text = element_text(size=14), legend.justification=c(0.5,0.5), text = element_text(family = "serif", face="bold"))

# Save the plot as a TIFF at 600 DPI
ggsave("./Documentation/Results/conf_matrix_rfe_class.tiff", plot = cm_plot, 
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

