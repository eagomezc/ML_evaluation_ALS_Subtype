#--------- ALS Subtype evaluation: Unsupervised Clustering -------------------#

# This script is to classify the samples from KCL BrainBank and ALS Consortium into
# 3 different subtypes. 

# -------------------------> LIBRARY LOAD: 

# To only use when you're working in the cluster:
 .libPaths("/scratch/users/k2584930/software/R/4.3/")

if (!require('NMF')) install.packages('NMF'); library('NMF')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('data.table')) install.packages('data.table'); library('data.table')

# Set directory to be in the source file location:
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
op <- par(no.readonly = TRUE) # Image display

# Set directory to be in the source file location (Cluster):
this_dir <- function(directory)
setwd( file.path(getwd(), directory) )

# Functions for NMF:

nmf_extract_feature <- function(nmf_final, rawdata=NULL, manual.num=0, method="default", math="mad", FScutoff=0.9){
  feature.score <- featureScore(nmf_final)
  predict.feature <- predict(nmf_final, what="features", prob=T)
  
  data.feature <- data.frame(Gene=names(feature.score),
                             featureScore=feature.score,
                             Group=predict.feature$predict,
                             prob=predict.feature$prob,
                             stringsAsFactors=FALSE)
  if(method=="total"){
    print("return all featureScores")
  }else{
    if(method=="default"){
      # extracted features for each group
      if (manual.num==0){
        extract.feature <- extractFeatures(nmf_final)
      }else if (manual.num>0 && manual.num<=length(featureNames(nmf_final))){
        extract.feature <- extractFeatures(nmf_final, manual.num)
      }else{
        stop("wrong number of (manual num) features ")
      }
      
      data.feature <- extract.feature %>%
        lapply(., function(x) data.feature[x, ]) %>%
        rbindlist %>%
        as.data.frame.matrix
    }else if(method=="rank"){
      if(is.null(rawdata)){
        stop("error: need to provide original expression data if method is 'rank' ")
      }
      data.feature <- cbind(data.feature, math=apply(rawdata, 1, math)) %>%
        filter(featureScore>=FScutoff) %>%
        arrange(Group, dplyr::desc(math), dplyr::desc(prob))
      if(manual.num>0){
        data.feature <- group_by(data.feature, Group) %>%
          top_n(manual.num)
      }
    }
  }
  return(data.feature)
}
nmf_extract_group <- function(nmf_final, type="consensus", matchConseOrder=F){
  data <- NULL
  if(type=="consensus"){
    predict.consensus <- predict(nmf_final, what="consensus")
    silhouette.consensus <- silhouette(nmf_final, what="consensus")
    # It turns out the factor levels is the NMF_assigned_groups from consensus matrix
    # that matches the original sampleNames(nmf_final) order
    # The attributes(a.predict.consensus)$iOrd is the idx order for it to match the
    # order of the samples in consensusmap(nmf_final). It is just for displaying
    # Therefore, the merged data frame sampleNames(nmf_final) + a.predict.consensus is the final
    # consensus nmf_finalults.
    data <- data.frame(Sample_ID=sampleNames(nmf_final),
                       nmf_subtypes = predict.consensus,
                       sil_width = signif(silhouette.consensus[, "sil_width"], 3))
    # If we want to display as we see in consensusmap, we just need to reoder everything.
    # Now re-order data to match consensusmap sample order
    if(matchConseOrder){
      sample.order <- attributes(predict.consensus)$iOrd
      data <- data[sample.order, ]
    }
  }else if(type=="samples"){
    predict.samples <- predict(nmf_final, what="samples", prob=T)
    silhouette.samples <- silhouette(nmf_final, what="samples")
    data <- data.frame(Sample_ID=names(predict.samples$predict),
                       nmf_subtypes = predict.samples$predict,
                       sil_width = signif(silhouette.samples[, "sil_width"], 3),
                       prob = signif(predict.samples$prob, 3))
  }else{
    stop(paste("Wrong type:", type, "Possible options are: 'consensus', 'samples' "))
  }
  return(data)
}
nmf_plot <- function(nmf_final, type="consensus", subsetRow=TRUE, save.image=F, hclustfun="average", silorder=F, add_original_name=T){
  if(save.image)
    pdf(nmf_final.pdf, width=18, height=15)
  
  if(type=="result"){
    print(plot(nmf_final))
  }else{
    si <- silhouette(nmf_final, what=type)
    
    if(type=="features"){
      if(silorder){
        basismap(nmf_final, Rowv = si, subsetRow=subsetRow)
      }else{
        basismap(nmf_final, subsetRow = subsetRow)
      }
    }else if(type=="samples"){
      if(silorder){
        coefmap(nmf_final, Colv = si)
      }else{
        coefmap(nmf_final)
      }
    }else if(type=="consensus"){
      if(add_original_name){
        colnames(nmf_final@consensus) <- sampleNames(nmf_final)
        rownames(nmf_final@consensus) <- sampleNames(nmf_final)
      }
      consensusmap(nmf_final, hclustfun=hclustfun)
    }
  }
  if(save.image)
    dev.off()
}

# -------------------------> DATA LOAD: 

# Open the tsv file the meta data information (diagnose, classification, tissue, 
# sample origin, etc.):

# Brain Bank Meta:
bb_meta <- read.delim("bb_meta_cluster.tsv")
bb_meta <- bb_meta[, c(1,6)]
bb_meta$platform <- "HiSeq_BB"

# ALS Consortium Meta:
es_meta <- read.csv("Cortex_Phenotype.csv")
es_meta <- es_meta[, c(1:5, 9, 10, 12:14)] # Keep just the relevant columns.
colnames(es_meta)[1] <- "Sample"

# Get relevant table:
# We want a table that only contains cases. In addition to that, we want a table
# with only motor cortex samples:
es_rel <- es_meta[es_meta$Subtype != "Control" & es_meta$tissue %like% "Motor", ]
es_rel <- es_rel[order(es_rel$tissue), ]

# Get only relevant table:
es_rel <- es_rel[, c(1,7,4)]
colnames(es_rel) <- c("Sample", "Cluster", "platform")

# Open the tsv file with the VST batch corrected counts for all the samples:
all_genes <- read.delim("Ens_vst_matrix_no_sex_top5000.tsv", row.names = 1)

# -----------------------> DATA MANIPULATION:

# Get a list of the different datasets (All samples will be es_rel):
# So I found that I can separate the samples by tissues, all of them belonging to
# motor or frontal cortex. So I can run the analysis with all the samples, as 
# Eschima did, or separated by tissue (which will get one tissue per sample); this
# coul be interesting because, for example, if Im able to classify subtypes using
# frontal cortex samples with models created with motor cortex, can be a sign of 
# how strong the models are.

# All Samples:
all_samples <- rbind(bb_meta, es_rel)
  
# -----------------------> NMF CLUSTERING: 

# Non-negative matrix factorization: an unsupervised clustering analysis that takes a
# matrix and breaks it down into smaller tables with no negative values, which allows for 
# data reduction, pattern identification and data clustering. 
  
# -----> 1. SELECT K setting with the highest cophenetic correlation coefficient:
# In theory, this has been done already! Eschima paper identified 3 cluster (including
# the one with the TE). I will do it only once, or maybe run to plot it once, but 
# not more. In theory, 3 clusters should be enough: 
  
# Cophenetic Correlation Coefficient:
# What it is: Measures the stability of sample clustering across multiple NMF runs.
# Interpretation: Values close to 1 → high clustering stability.
# Drop in value → instability; ideal K is just before the drop.
# Use it to pick the best K.
  
# Dispersion:
# What it is: Measures how dispersed (variable) the consensus matrix is.
# Interpretation: Lower values → more stable clustering (less variation).
# Flat or minimal dispersion → ideal clustering.
# Helps support cophenetic, not standalone.
  
# Silhouette (Consensus Silhouette Width):
# What it is: Measures how well samples are grouped within their clusters vs. other clusters.
# Range: from -1 to 1.
# Interpretation:
# Near 1 → samples are well clustered (good separation).
# Near 0 → overlapping clusters.
# Negative → likely misclassification.
# Useful to validate biological separability.
  
estim.r <- nmf(all_genes, 2:8, nrun=100, "nsNMF", .opt=paste0("vp", 20), seed=123211, maxIter=1000)
  
# Get the table with the measures:
estim_measures <- estim.r$measures
estim_measures <- estim_measures[, c(1, 13:15)]
  
# Create empty plot list:
estim_list <- list()
  
for (m in 2:4) {
  
  meas <- colnames(estim_measures)[m]
    
   measurements <- 
     ggplot(estim_measures, aes(rank, !!sym(meas))) +
     geom_line() +
     geom_point(size = 3) +
     xlab("Rank") +
     ylab(meas) +
     scale_x_continuous(breaks = unique(estim_measures$rank)) +  
     theme(axis.title = element_text(size = 25, hjust = 0.5),
           plot.title = element_text(size = 28, hjust = 0.5),
           axis.text.x  =  element_text(size = 22, colour = "black"),
           axis.text.y  = element_text(size = 22, colour = "black"),
           legend.title = element_text(size = 22),
           legend.text = element_text(size = 20),
           legend.position = "right",
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           axis.line = element_line(colour = "black")) 
    
  # Save every plot:
  estim_list[[m]] <- measurements # From the previous empty list, save the plots.
  
  }
  
# Since the for start with, we need to delete the empty 3 elements in the list:
estim_list <- estim_list[lengths(estim_list) > 0]
  
png(file = "../../output/nmf_measures.png", width = 1500, height = 800)
  
# Save all the plots in one grid: 
gridExtra::grid.arrange(grobs = estim_list, ncol = 3)

dev.off()

  
# ------> 2. RUN NMF with selected k settings (3 in the all one):
# I checked just in case, but we will need 3 clusters if we want to do the 
# evaluation of the models. For that reason, even I may run the measurement plots
# for all the models, I still going to use 3 K for all of them.
  
final_nmf <- nmf(all_genes, 3, nrun=100, "nsNMF", .opt=paste0("vp", 20), seed=123211, maxIter=1000)
  
# ------> 3. EXTRACT information from the models:

# Extract important features for every cluster:
rel_genes <- nmf_extract_feature(final_nmf, method="default", math="mad")
  
write.table(rel_genes, 
            file = "../../output/Cluster_features.tsv", 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
  
# All the genes by cluster: 
all_genes <- nmf_extract_feature(final_nmf, method="total", math="mad")
  
write.table(all_genes, 
            file = "../../output/Cluster_all_features.tsv", 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
  
# Extract groups. Samples within clusters:
# Get the samples by group after consensus of the analysis:
# IMPORTANT!!!!! Consensus re-lable the groups in a way that can be different from 
# the important features clustering (rel_genes). Meaning that what is a important 
# gene from cluster 1, in the consensus, is refering to cluster 2.
# The tranformation is as follow:
  
# Sample ---> Consensus:
#   1             2
#   2             3
#   3             1
  
con_group <- nmf_extract_group(final_nmf, type="consensus", matchConseOrder = TRUE)
  
write.table(con_group, 
            file = "../../output/Cluster_Consensus.tsv", 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
  
# Get the samples by group after nmf analysis (no dendogram):
sample_group <- nmf_extract_group(final_nmf, type="samples", matchConseOrder = TRUE)
  
write.table(sample_group, 
            file = "../../output/Cluster_Samples.tsv", 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
  
# Together table:
# Put together sample and consensus clusters:
tog_group <- merge(con_group, sample_group, by="Sample_ID", suffixes=c("_con", "_sample"))
tog_group <- tog_group[, c(1,4,2,3,5,6)]
  
write.table(tog_group, 
            file = "../../output/Cluster_consensus_samples.tsv", 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
  
# Plot results:
  
# Plot the genes and how they generate the clusters:
png(file = "../../output/nmf_features.png", width = 1500, height = 800)
  
nmf_plot(final_nmf, type="features", silorder=T)
  
dev.off()
  
# Plot the genes with consensus and samples grouping (as it can be seen grouping is
# almost the same):
png(file = "../../output/nmf_consensus.png", width = 1500, height = 800)
  
nmf_plot(final_nmf, type="consensus", silorder=T)
  
dev.off()
  
# Plot heatmap with samples and features: 
  
png(file = "../../output/nmf_samples.png", width = 1500, height = 800)
  
nmf_plot(final_nmf, type="samples", silorder=T)
  
dev.off()
  
  
  
