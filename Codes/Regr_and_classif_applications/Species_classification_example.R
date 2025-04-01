# Classification Problem Example usage
#-------------------------------------------------------------------------------
# Application of Resampled EnMap Hyperspectral Data in Species Classification
# ------------------------------------------------------------------------------
# Source the script ------------------------------------------------------------
source("./Codes/OptiSpecPickerClass/OptiSpecPickerClass.R")
# ------------------------------------------------------------------------------
# Load the data ----------------------------------------------------------------
hyperdata <- read.csv("./Data/Classification-Species/parthenium_weed_data.csv", header=TRUE, as.is=TRUE, check.names=FALSE)
head(hyperdata)
str(hyperdata)

library(hsdar)
# convert the data.frame to matrix
spectra_m<-as.matrix(hyperdata[,2:2152])
str(spectra_m)

# Get the wavelengths from names of original data
wavelength<-as.numeric(names(hyperdata[,2:2152]))
str(wavelength)

# Create Spec lib
field_data<-speclib(spectra_m, wavelength)
str(field_data)

species<-hyperdata$species
SI(field_data) <- as.character(species)
str(field_data)
names(SI(field_data)) <- "Species"

# Plot
par(mfrow = c(1,1))
plot(field_data,col= "red", xlim= c(400,2500), legend = list(x = "topleft"))

# Mask 
hsdar::mask(field_data) <- c(350,399,1350,1460,1790,1960,2300, 2500)

par(mfrow=c(1,1))
plot(field_data, FUN = 1)

# Resampling to EnMap
resampled_enmap_data <- spectralResampling(field_data, "EnMAP", response_function = F)

# Print resampled wv
wavelength(resampled_enmap_data)
# Plot resampled spectra
par(mfrow=c(1,1))
plot(resampled_enmap_data)

#check the structure for Resampled EnMap data
str(resampled_enmap_data)
resampled_enmap_res<-spectra(resampled_enmap_data)

#check structure again
str(resampled_enmap_res)

# Extract the SI (Species names) from speclib and Wavelengths
class_names <- SI(resampled_enmap_data)
column_bandnames <- round(wavelength(resampled_enmap_data), 2)

# Rename the data
resampled_enmap_res <- as.data.frame(cbind(SI,resampled_enmap_res))
head(resampled_enmap_res)
length(resampled_enmap_res)

names(resampled_enmap_res)[2:243] <- as.character(column_bandnames)

head(resampled_enmap_res)

write.csv(resampled_enmap_res, "/Users/mahlatse/Documents/OptiSpecPicker/Data/Classification-Species/Parthenium_weed_enmap.csv")

# Prepare dataset for OptispecPicker and classification
y = resampled_enmap_res$Species
X = resampled_enmap_res[,-1]
dataset <- cbind(X, y)
head(dataset)

# Run the @OptiSpecPickerRegr for feature selection
#--------------OptiSpecPickerRegr Parameters------------------------------------ 
#' @param X A data frame of explanatory variables.
#' @param y A factor vector of the response variable.
#' @param patience An integer indicating the number of iterations to wait for improvement before stopping (default is 5). When working with highly dimensional data such as hyperspectral data, a higher patience can be used.
#' @param decision_threshold A numeric value to decide if the second-best accuracy throughout all iterations is close enough to the best accuracy. It is used to determine the optimal subset. If the difference between best accuracy/RMSE and second-best accuracy(classification)/RMSE (regression) is negligible, fewer features are prioritised. If NULL, this comparison is not made (default is NULL).
#' 
#-------------------------------------------------------------------------------
# OptiSpecPicker Feature Selection ---------------------------------------------
# Reproducibility
set.seed(432)
optimal_features <- OptiSpecPickerClass(dataset[,c(1:242,243)], y, decision_threshold = 0.10)

str(optimal_features)
head(optimal_features$results$Best_Subset)

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
best_subset <- optimal_features$results$Best_Subset[1]
selected_bands <- as.numeric(unlist(strsplit(best_subset, ", ")))

# Find the indices of the selected bands closest to EnMap wavelengths
selected_band_indices <- sapply(selected_bands, function(band) which.min(abs(enmap_wavelengths - band)))

# Plot the average spectra for each class with highlighted selected bands
aver_spectra_optispecpicker <- ggplot(average_spectra_long, aes(x = as.numeric(Band), y = Reflectance, color = Class, group = Class)) +  # Add group = Class
  geom_line(size = 1) +  # Plot average spectral curve for each class
  geom_vline(xintercept = enmap_wavelengths[selected_band_indices], linetype = "dashed", color = "red", size = 1) +  # Highlight selected bands
  geom_point(data = filter(average_spectra_long, Band %in% enmap_wavelengths[selected_band_indices]),
             aes(x = Wavelength, y = Reflectance), color = "red", size = 3) +  # Highlight selected bands as points
  xlab("Wavelength (nm)") +
  ylab("Reflectance (%)") +
  ggtitle("") +
  theme_minimal() +
  scale_color_manual(values = c("AT" = "blue", "GS" = "green", "OPS" = "orange", "PH" = "purple"))

# Save the plot as a TIFF at 600 DPI
ggsave("./Documentation/Results/aver_spectra_optispecpicker.tiff", plot = aver_spectra_optispecpicker, 
       device = "tiff", 
       dpi = 600, 
       width = 8, 
       height = 6, 
       units = "in")

#-------------------------------------------------------------------------------
library(skimr)
skim(optimal_features$results)

print(as.data.frame(optimal_features$results))
write.csv(optimal_features$results, "optimal_featuresclass.csv")

# Print the structure 
str(optimal_features)
# Print the final data 
print(head(as.data.frame(optimal_features$final_data)))
# Export the results to csv.
write.csv(optimal_features$results, "./Documentation/Results/optimal_dataset_class.csv")
# Print the results
print(as.data.frame(optimal_features$results))
# Export the results to csv.
write.csv(optimal_features$results, "./Documentation/Results/optispecpickerclass_results.csv")

# ------------------------------------------------------------------------------
# Stability Test ---------------------------------------------------------------
# load required libraries
#-------------------------------------------------------------------------------
# ------------------------Stability Test----------------------------------------
# Run the function 50 times and store the results
num_runs <- 50
results_list_class <- vector("list", num_runs)

for (run in 1:num_runs) {
  set.seed(run) # Set seed for reproducibility
  results_list_class[[run]] <- OptiSpecPickerClass(dataset[,c(1:242,243)], y, decision_threshold = 0.10)
}

# Optionally, save the results to an RDS file for later use
saveRDS( results_list_class, file = "optispec_results_class.rds")
#-------------------------------------------------------------------------------
# Compute the Frequency of selection from 50 runs ------------------------------
# Libraries
library(dplyr)
library(ggplot2)

# Extract the best subsets from each result
best_subsets <- lapply(results_list_class, function(result) {
  strsplit(result$results$Best_Subset[1], ", ")[[1]]
})

# Flatten the list of best subsets into a single vector
all_selected_variables <- unlist(best_subsets)

# Count the frequency of each variable
variable_frequency <- table(all_selected_variables)
variable_frequency_df <- as.data.frame(variable_frequency)

# Rename columns for clarity
colnames(variable_frequency_df) <- c("Variable", "Frequency")

# Calculate the 25th and 75th percentiles of the frequencies
first_quartile <- quantile(seq(c(1:50)), 0.25)
third_quartile <- quantile(seq(c(1:50)), 0.75)

# Plot the frequency of selected variables [plot frequency over the spectra]
ggplot(variable_frequency_df, aes(x = reorder(Variable, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_hline(yintercept = first_quartile, linetype = "dashed", color = "red") +
  geom_hline(yintercept = third_quartile, linetype = "dashed", color = "green") +
  coord_flip() +
  xlab("Variable") +
  ylab("Frequency") +
  ggtitle("Frequency of Selected Variables over 50 Runs") +
  annotate("text", x = nrow(variable_frequency_df) + 1, y = first_quartile, 
           label = paste("25th percentile:", round(first_quartile, 2)), 
           hjust = -0.1, vjust = -0.5, color = "red") +
  annotate("text", x = nrow(variable_frequency_df) + 1, y = third_quartile, 
           label = paste("75th percentile:", round(third_quartile, 2)), 
           hjust = -0.1, vjust = -0.5, color = "green") +
  theme_minimal()
#include IQR lines on the plot 25 Lower (@13) and upper: 75 (38)

min_frequency_threshold <- 0.20 * num_runs
filtered_df <- variable_frequency_df[variable_frequency_df$Frequency >= min_frequency_threshold, ]

library(ggplot2)

enmap_freq_50runs <- ggplot(filtered_df, aes(x = reorder(Variable, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_hline(yintercept = first_quartile, linetype = "dashed", color = "red") +
  geom_hline(yintercept = third_quartile, linetype = "dashed", color = "green") +
  coord_flip() +
  xlab("Variable") +
  ylab("Frequency") +
  ggtitle("Frequency of Selected Variables (≥ 20% Selection) over 50 Runs") +
  
  # Annotations with smaller text and adjusted positioning
  annotate("text", x = nrow(filtered_df) + 1, y = first_quartile, 
           label = paste("25th percentile:", round(first_quartile, 2)), 
           hjust = -0.1, vjust = -0.5, size = 3, color = "red") +
  annotate("text", x = nrow(filtered_df) + 1, y = third_quartile, 
           label = paste("75th percentile:", round(third_quartile, 2)), 
           hjust = -0.1, vjust = -0.5, size = 3, color = "green") +
  
  # Adjust y-axis labels for better readability
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),     # Reduce y-axis label size
    axis.text.x = element_text(size = 10),    # Adjust x-axis label size
    plot.margin = unit(c(1, 1, 1, 1), "cm"),  # Increase plot margins using unit()
    axis.title.x = element_text(size = 12),   # Adjust x-axis title size
    axis.title.y = element_text(size = 12),   # Adjust y-axis title size
    plot.title = element_text(size = 14, hjust = 0.5)  # Center plot title
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))  # Add space above bars

# Save the plot as a TIFF at 600 DPI
ggsave("./Documentation/Results/enmap_freq_50runs.tiff", plot = enmap_freq_50runs, 
       device = "tiff", 
       dpi = 600, 
       width = 8, 
       height = 6, 
       units = "in")

#-------------------------------------------------------------------------------
# Frequency 
library(plotly)
library(dplyr)

# Extract the best subsets from each result (unique bands per run)
best_subsets_class <- lapply(results_list_class, function(result) {
  unique(strsplit(result$results$Best_Subset[1], ", ")[[1]])
})

# Flatten the list of best subsets into a single vector of selected bands
all_selected_bands <- unlist(best_subsets_class)

# Count the frequency of each selected band (ensuring each band is counted only once per run)
selected_band_frequency <- table(as.numeric(all_selected_bands))
selected_band_frequency_df <- as.data.frame(selected_band_frequency)

# Rename columns for clarity
colnames(selected_band_frequency_df) <- c("Band", "Frequency")

# Convert Band column to numeric for proper ordering on the plot
selected_band_frequency_df$Band <- as.numeric(as.character(selected_band_frequency_df$Band))

# Calculate average spectra by class
avg_spectra_by_class <- enmap_spectra %>%
  group_by(Class) %>%
  summarise(across(2:242, mean, na.rm = TRUE))

# Melt the dataframe for easier plotting (if using tidyr package)
avg_spectra_melt <- avg_spectra_by_class %>%
  pivot_longer(cols = 2:242, names_to = "Band", values_to = "Reflectance")

# Convert Band to numeric for proper ordering on the plot
avg_spectra_melt$Band <- as.numeric(gsub("V", "", avg_spectra_melt$Band))

# Scale the frequency so it fits nicely on the same plot as reflectance
max_reflectance <- max(avg_spectra_melt$Reflectance, na.rm = F)
max_frequency <- max(selected_band_frequency_df$Frequency)

# Scale frequency to match reflectance values
selected_band_frequency_df$Scaled_Frequency <- selected_band_frequency_df$Frequency / max_frequency * max_reflectance

# Scale the frequency of selection based on reflectance range
selected_band_frequency_df$Scaled_Frequency <- selected_band_frequency_df$Frequency * max(avg_spectra_melt$Reflectance) / max(selected_band_frequency_df$Frequency)

# Create a plot using plotly
p <- plot_ly() %>%
  # Add the averaged spectra for each class
  add_trace(data = avg_spectra_melt, 
            x = ~Band, y = ~Reflectance, 
            color = ~Class, type = 'scatter', mode = 'lines',
            name = ~Class) %>%
  # Add bar plot for frequency of selected bands on the secondary axis
  add_trace(data = selected_band_frequency_df, 
            x = ~Band, y = ~Frequency, 
            type = 'bar',name = "",
            yaxis = "y2", marker = list(color = 'red', opacity = 0.6)) %>%
  # Set titles and axis labels
  layout(
    title = "",
    xaxis = list(title = "Wavelength (nm)"),
    yaxis = list(title = "Reflectance"),
    yaxis2 = list(title = "Frequency of Selection", overlaying = "y", side = "right"),
    legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.1)
  )

# Display the plot
p
#-------------------------------------------------------------------------------
# Classification using Optimal Subset ------------------------------------------
library(mlr)
# Enmap OptiSpecPicker
head(optimal_features$final_data)
# View the current column names
colnames(optimal_features$final_data)

# Rename columns with valid names
colnames(optimal_features$final_data) <- make.names(colnames(optimal_features$final_data))

# Verify the new column names
colnames(optimal_features$final_data)

enmap_class_task <- mlr::makeClassifTask(data = optimal_features$final_data, target = "y", id="Optimal subset")
enmap_class_task 

set.seed(153)
n = getTaskSize(enmap_class_task)
train.set = sample(n, size = round(2/3 * n))
test.set = setdiff(seq_len(n), train.set)
#-------------------------------------------------------------------------------
# -------------------Algorithm Settings-----------------------------------------
#-------------------------------------------------------------------------------
as.data.frame(mlr::listLearners())
# -------------------Generic Settings-------------------------------------------
# Set Cross-validation strategy
rdesc = makeResampleDesc("CV", stratify = F, iters =10)
# ------------------- Random Forest Parameters----------------------------------
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

# -------------------Computing settings-----------------------------------------
# set parallel back-end
require(parallel)
require(parallelMap) 
parallelStartSocket(cpus = detectCores())
#-------------------------------------------------------------------------------
# -------------------Parameter tuning-------------------------------------------
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
# -------------------Model Train and Evaluation---------------------------------
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
