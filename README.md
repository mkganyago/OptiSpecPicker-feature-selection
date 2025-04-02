# OptiSpecPicker: Wine-Making Inspired Feature Selection for Remote Sensing  
**Code and Data Repository**  

This repository accompanies the manuscript *"OptiSpecPicker: A Wine-Making Inspired Algorithm for Spectral Feature Space Optimization in Remote Sensing Classification and Regression"* submitted to the *ISPRS Journal of Photogrammetry and Remote Sensing*. It provides code, data, and examples to reproduce the results and apply the algorithm to new datasets.  

---
## 📂 Repository Structure  
```plaintext
├── Code/  
│   ├── benchmark/                  # Benchmark scripts (RF-RFE, full dataset comparisons)  
│   ├── OptiSpecPickerClass/        # Functions for classification tasks  
│   ├── OptiSpecPickerRegr/         # Functions for regression tasks  
│   └── Regr_and_classif_applications/  
│       ├── Regression_LCC_prediction_example.R  # Example: Leaf Chlorophyll Content (LCC) retrieval  
│       └── Species_classification_example.R     # Example: Invasive species discrimination  
├── Data/  
│   ├── EnMAP_resampled/            # Resampled EnMAP hyperspectral data (Parthenium hysterophorus use case)  
│   └── Bothaville_SA/              # Experimental Sentinel-2 & field data (LCC retrieval use case)  
├── LICENSE                         # MIT License  
└── README.md                       # This file  
 
```markdown
# OptiSpecPicker: Wine-Making Inspired Feature Selection for Remote Sensing  
**Code and Data Repository**  

This repository accompanies the manuscript *"OptiSpecPicker: A Wine-Making Inspired Algorithm for Feature Space Optimization in Remote Sensing Classification and Regression"* submitted to the *ISPRS Journal of Photogrammetry and Remote Sensing*. It provides code, data, and examples to reproduce the results and apply the algorithm to new datasets.  

---
## 📂 Repository Structure  
```plaintext
├── Code/  
│   ├── benchmark/                  # Benchmark scripts (RF-RFE, full dataset comparisons)  
│   ├── OptiSpecPickerClass/        # Functions for classification tasks  
│   ├── OptiSpecPickerRegr/         # Functions for regression tasks  
│   └── Regr_and_classif_applications/  
│       ├── Regression_LCC_prediction_example.R  # Example: Leaf Chlorophyll Content (LCC) retrieval  
│       └── Species_classification_example.R     # Example: Invasive species discrimination  
├── Data/  
│   ├── EnMAP_resampled/            # Resampled EnMAP hyperspectral data (Parthenium hysterophorus study)  
│   └── Bothaville_SA/              # Experimental Sentinel-2 & field data (LCC retrieval study)  
├── LICENSE                         # MIT License  
└── README.md                       # This file  
```

---
## 🚀 Quick Start  

### **Dependencies**  
- **R** (≥ v4.1.2)  
- R Packages:  
  ```R
  install.packages(c("randomForest", "caret, "mlr", "dplyr", "ggplot2", "FSelectorRcpp"))
  ```

### **Installation**  
1. Clone the repository:  
   ```bash
   git clone https://github.com/mkganyago/OptiSpecPicker-feature-selection.git  
   ```
2. Navigate to the `Regr_and_classif_applications/` folder to run example scripts.  

---

## 🛠️ Usage  

### **1. Classification Example**  
**Script**: `Species_classification_example.R`  
**Data**: `Data/EnMAP_resampled/`  
**Task**: Discriminate *Parthenium hysterophorus* using resampled EnMAP hyperspectral bands.  
```R
source("../OptiSpecPickerClass/OptiSpecPicker.R")  
data <- read.csv("Data/EnMAP_resampled/enmap_data.csv")  
results <- OptiSpecPickerClass(X = data[, -ncol(data)], Y = data$Class, u = 10, m = 0.01)  
```

### **2. Regression Example**  
**Script**: `Regression_LCC_prediction_example.R`  
**Data**: `Data/Bothaville_SA/`  
**Task**: Retrieve Leaf Chlorophyll Content (LCC) from Sentinel-2 spectral features.  
```R
source("../OptiSpecPickerRegr/OptiSpecPicker.R")  
data <- read.csv("Data/Bothaville_SA/sentinel2_lcc_data.csv")  
results <- OptiSpecPickerRegr(X = data[, -ncol(data)], Y = data$LCC, u = 5, m = 0.05)  
```

### **Key Parameters**  
- `u`: Patience parameter (stops iterations if no improvement in `u` steps).  
- `m`: Decision threshold for a trade-off between accuracy and dimensionality.  

---

## 📊 Data Description  

### **1. EnMAP Resampled Data**  
- **Source**: Field spectroscopy data resampled to EnMAP sensor specifications.  
- **Format**: CSV with 242 spectral bands (350–2500 nm) and class labels (*Parthenium hysterophorus*, grasses, acacia trees).  
- **Use Case**: Hyperspectral classification of invasive species.  

### **2. Bothaville Experimental Data**  
- **Source**: Sentinel-2 MSI bands + vegetation indices paired with field-measured LCC.  
- **Format**: CSV with 22 features (10 spectral bands + 12 VIs) and LCC values (μg/m²).  
- **Use Case**: Multispectral regression for crop health monitoring.  

---

## 🔍 Benchmarking  
The `benchmark/` folder includes scripts to compare OptiSpecPicker against:  
1. **Recursive Feature Elimination (RF-RFE)**  
2. **Full Dataset Models**  
Run `benchmark/classification_benchmark.R` or `benchmark/regression_benchmark.R` to reproduce results from the manuscript.  

---

## 🤝 Contributing  
Contributions are welcome! Open an issue or submit a pull request for:  
- Bug fixes  
- Performance optimizations  
- Extended applications

---

## 📜 License  
MIT License. See [LICENSE](LICENSE) for details.  

---

## 📄 Citation  
If you use this code or data, please cite:  
```bibtex
@Article{Kganyago2024OptiSpecPicker,  
  title = {OptiSpecPicker: A Wine-Making Inspired Algorithm for Feature Space Optimization in Remote Sensing Classification and Regression},  
  author = {Mahlatse Kganyago},  
  journal = {Submitted to a Journal},  
  year = {2025}  
}  
```

---

## 📧 Contact  
For questions or collaborations, contact:  
**Mahlatse Kganyago**  
Email: [mahlatsek@uj.ac.za](mailto:mahlatsek@uj.ac.za)  
Affiliation: University of Johannesburg, South Africa  
``` 
