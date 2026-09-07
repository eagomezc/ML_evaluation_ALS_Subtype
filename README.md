# Machine learning evaluation of gene expression-based ALS subtypes across brain and blood tissues

## Overview:
This repository contains the main **R scripts** used to evaluate the gene expression-based ALS subtypes across brain and blood tissues. 

Transcriptomic data from two **post-mortem motor cortex** ALS samples were obtained from the [**London Neurodegenerative Diseases Brain Bank**](https://www.kcl.ac.uk/neuroscience/facilities/brain-bank) (112 ALS samples) and the [**NYGC ALS Consortium**](https://www.nygenome.org/science-technology/collaborative-research-programs/neurodegenerative-disease-research/als-consortium/) (257 samples, GEO accession code [GSE153960](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153960)). An additional peripheral blood RNA-seq dataset from 96 ALS patients (MQND Biobank dataset) was obtained from GEO (Accession code [GSE234297](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234297)). 

To confirm the presence of ALS subtypes in motor cortex samples, **unsupervised clustering analysis** and **enriched pathway analysis** were performed. Afterwards, three **machine learning strategies** were used to develop ALS subtype classifiers, which were evaluated in an independent cohort. **Multi-class linear discriminant analysis (LDA)** models were then used for ALS subtype classification in blood samples. Evaluation of the LDA models was performed indirectly using a **three-way differential gene expression analysis**. Finally, the association between disease stage and ALS subtype classification of pre-mortem blood samples was evaluated using **Spearman correlation** and **Beta regression analysis**. 

## System Requirements: 

### Hardware requirements: 

All the scripts used for this work were run on a MacBook Air (Apple M4 chip, RAM: 24GB, CPU: 10 cores). 

High-memory and long jobs (including unsupervised clustering and parameter tuning for machine learning models) were performed in the high-performance computing cluster of KCL (CREATE; more information [here](https://docs.er.kcl.ac.uk/)). Running time in the cluster depends of memory availability. 

### System requirements: 

All the scripts were created and used using the **R Programming Language** on a Mac environment. 

**R version**: 4.4.3

**RStudio version**: 2026.01.1+403

To install R and RStudio, follow the installation instructions [here](https://www.stats.bris.ac.uk/R/) and [here](https://www.rstudio.com/products/rstudio/download/). Typical installation times do not exceed a few minutes.

### Required R packages (libraries): 

The required packages to run all scripts should be installed automatically when you run each script; however, in case you need to install them manually, you can do it as follows (showing only the main packages; other packages such as ggplot2, dplyr, etc, should be installed as well):

```
# Main packages for Unsupervised Clustering:
if (!require('NMF')) install.packages('NMF'); library('NMF') # The non-smooth negative factorisation algorithm

# Main packages for Enrichment Analysis:
if (!require('clusterProfiler')) BiocManager::install('clusterProfiler', update = FALSE); library('clusterProfiler') # GO Enrichment
if (!require('org.Hs.eg.db')) BiocManager::install("org.Hs.eg.db", character.only = TRUE); library("org.Hs.eg.db") # H. sapiends database
if (!require('ReactomePA')) BiocManager::install("ReactomePA", character.only = TRUE); library("ReactomePA") # Reactome Enrichment

# Main packages for Supervised Machine Learning:
if (!require('randomForest')) install.packages('randomForest'); library('randomForest') # Random forest
if (!require('caret')) install.packages('caret'); library('caret') # Machine learning framework
if (!require('glmnet')) install.packages('glmnet'); library('glmnet') # To run Elastic Net algorithm
if (!require('xgboost')) install.packages('xgboost'); library('xgboost') # XGBoost

# Main packages for Linear Discriminant Analysis (LDA):
if (!require('caret')) install.packages('caret'); library('caret') # Machine learning framework
if (!require('MASS')) install.packages('MASS'); library('MASS')

# Main packages for three-way differential gene expression analysis:

if (!require('volcano3D')) install.packages('volcano3D'); library('volcano3D') # 3D Volcano Plot
if (!require('htmlwidgets')) install.packages('htmlwidgets'); library('htmlwidgets') # Manage 3D Plots interactions
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if (!require('DESeq2')) BiocManager::install('DESeq2', update = FALSE); library('DESeq2') # For DGE statistics

# Main packages for association between disease stage and ALS subtype classification:
if (!require('betareg')) install.packages('betareg'); library('betareg') # Beta regression models
if (!require('scales')) install.packages('scales'); library('scales')

```

## Content:

The repository contains two folders. The scripts used for the analysis and the expected outputs. 

### [1_Scripts](https://github.com/eagomezc/ML_evaluation_ALS_Subtype/tree/main/1_Scripts)

This folder contains the scripts for clustering, machine learning, blood evaluation and disease status analysis. 

The scripts are: 

**1_unsupervised_clustering.R:** This script takes VST-normalized read counts and performs unsupervised clustering using the Non-negative matrix factorisation algorithm.

**2_pathway_enrichment_analysis.R:** This script takes a list of relevant genes for cluster identification (for each subtype) and performs pathway enrichment analysis to identify which molecular signals are associated with each subtype.

**3_supervised_machine_learning.R:** This script takes VST-normalized read counts and builds machine learning models able to classify ALS samples into different subtypes. Three strategies are used: random forest, elastic net and gradient tree boosting. 

**4_linear_discriminant_analysis.R:** This script takes candidate genes from ALS subtype classifiers and VST-normalized blood read counts for the multi-classification of ALS blood samples into three different subtypes using LDA.

**5_three_way_DGE.R:** This script performs a three-way differential gene expression analysis of blood samples, identifies upregulated and downregulated genes and, based on these genes, performs enrichment analysis as an indirect way to evaluate the classification ability of the LDA models. 

**6_disease_status_and_prediction.R:** This script evaluates the association between disease status in ALS blood samples (Collection time) and the prediction ability of LDA models. 

### [2_Expected_Output](https://github.com/eagomezc/ML_evaluation_ALS_Subtype/tree/main/2_Expected_Output)

This folder contains, separated by subfolders, the different expected outputs that can be obtained after running the scripts. Each subfolder has the name of the specific script that generates it, in addition to the number of the script, to make more clear what file is the result of which script.

The subfolders are:

**1_unsupervised_clustering:** The outputs are clustering plots, metrics of clustering performance, sample classification in each cluster and relevant genes for sample classification.

**2_pathway_enrichment_analysis:** The outputs are enriched pathway plots and tables for Gene Ontology terms and Reactome signals. 

**3_supervised_machine_learning:** The outputs are an accuracy plot, confusion matrices, parameter tuning, feature selection plots and the machine learning models for each methodology. 

**4_linear_discriminant_analysis:** The outputs are LD projection plots for training and evaluation datasets, and tables with the probability of classification for each sample into one of the ALS subtypes. 

**5_three_way_DGE:** The outputs are upregulated and downregulated genes for each of the ALS subtypes, volcano 3D plots (top and interactive), and enriched pathway plots and tables based on differentially expressed genes for each subtype.

**6_disease_status_and_prediction:** The outputs are correlation plots between collection point and probability of prediction, and regression model results (plots, coefficients, etc.) highlighting the association between disease status and the ability of the machine learning models to classify blood samples. 

## Publication: 

If you use these codes for your research please cite these paper:
