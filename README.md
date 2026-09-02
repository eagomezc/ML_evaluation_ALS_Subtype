# Machine learning evaluation of gene expression-based ALS subtypes across brain and blood tissues

## Overview:
This repository contains the main **R scripts** used to evaluate the gene expression-based ALS subtypes across brain and blood tissues. 

Transcriptomic data from two **post-mortem motor cortex** ALS samples were obtained from the [**London Neurodegenerative Diseases Brain Bank**](https://www.kcl.ac.uk/neuroscience/facilities/brain-bank) (112 ALS samples) and the [**NYGC ALS Consortium**](https://www.nygenome.org/science-technology/collaborative-research-programs/neurodegenerative-disease-research/als-consortium/) (257 samples, GEO accession code [GSE153960](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153960)). An extra, peripheral blood RNA-seq data from 96 ALS patients (MQND Biobank dataset) was obtained from GEO (Accession code [GSE234297](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234297)). 

To confirm the presence of ALS subtypes in motor cortex samples **unsupervised clustering analysis** and **enriched pathway analysis** were performed. After that, three **machine learning strategies** were used to create ALS subtype classifiers that were evaluated in an independent cohort. **Multi-class linear discriminant analysis (LDA)** models were then used for ALS subtype classification in blood samples. Evaluation of the LDA models were performed indirectly using a **three-way diferential gene expression analysis**. Finally, the association between disease stage and ALS subtype classification of pre-mortem blood samples was evaluated using **Spearman correlation** and **Beta regression analysis**. 


