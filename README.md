BINF702 Project
Date: April 15, 2024

**Overview**
This project analyzes the SRBCT dataset, which consists of gene expression data from 83 patients across four different cancer types:

Ewing sarcoma (EWS) (29 samples)
Burkitt lymphoma (BL) (11 samples)
Neuroblastoma (NB) (18 samples)
Rhabdomyosarcoma (RMS) (25 samples)
The dataset contains 2,308 gene expression measurements per patient. The analysis includes exploratory data analysis, hypothesis testing, clustering, classification, and principal component analysis (PCA).

**Requirements**
To run this project, you need the following dependencies in R:
install.packages(c("plsgenomics", "outliers", "rpart", "rpart.plot"))
Data Description
The project uses the SRBCT dataset from the plsgenomics package. Key variables include:

SRBCT$X: Gene expression data (83 samples × 2,308 genes).
SRBCT$Y: Cancer type labels (1 = EWS, 2 = BL, 3 = NB, 4 = RMS).
SRBCT$gene.names: Gene identifiers and descriptions.

**Analysis Steps**
1. Data Preprocessing
Checked for missing values and duplicate gene identifiers.
Created unique gene IDs.
Sorted patients based on cancer type.
Split the dataset into training (65 samples) and testing (18 samples) sets.
2. Hypothesis Testing & Data Assessment
Normality check: Used the Shapiro-Wilk test to assess whether gene expression values follow a normal distribution.
Outlier detection: Used Grubbs' test to identify potential outliers in gene expression data.
3. Clustering Analysis
Hierarchical clustering (Dendrogram): Used Euclidean distance and complete linkage to group patients by gene expression similarity.
Identified six clusters using hierarchical clustering.
4. Classification (CART Model)
Applied Classification and Regression Trees (CART) to classify cancer types based on gene expression.
Evaluated model performance using a confusion matrix and error rate calculation.
5. Principal Component Analysis (PCA)
Variance analysis:
First 10 principal components explain 81% of variance.
First 2 principal components explain 63% of variance.
Bootstrap confidence intervals estimated for eigenvalues.
Biplot visualization of PCA results.

**How to Run the Code**
Open RStudio or an R environment.

**Load the required libraries:**

library(plsgenomics)
library(outliers)
library(rpart)
library(rpart.plot)
Run the script section by section, following the analysis flow.

**Results & Insights**
EWS and RMS samples showed lower normality and more outliers, suggesting variability in gene expression.
BL samples had the most normally distributed gene expressions.
PCA analysis identified key genes that differentiate BL from other cancer types.
CART model achieved an error rate of 0.036, indicating strong classification performance.
