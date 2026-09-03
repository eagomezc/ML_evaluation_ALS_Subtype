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


```

## Content:

### a_Scripts

### b_Scripts


