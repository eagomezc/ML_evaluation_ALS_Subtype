#--------- ALS Subtype Evaluation: Linear Discriminant Analysis ---------------#

# After creating the machine learning models (and improved them), now we want to 
# explore if we can use motor cortex based models to classify blood samples. 

#-----------------------------> LIBRARY LOAD:

if (!require('randomForest')) install.packages('randomForest'); library('randomForest')
if (!require('caret')) install.packages('caret'); library('caret')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('MASS')) install.packages('MASS'); library('MASS')
if (!require('tidyr')) install.packages('tidyr'); library('tidyr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')

set.seed(42) # To get same results even with the random part.
options(digits = 3) # To get only part of the decimals. 

# Set directory to be in the source file location (Rstudio):
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Function to calculate weighted sensitivity and specificity for multiclass classification:

weighted <- function(conf_table) {
  
  # Create an empty table to save the value: 
  conf_metrics <- data.frame()
  
  # Call every subtype:
  for (n in unique(conf_table$Prediction)) {
    
    # Calculate the Sensitivity (TP / (TP + FN)) and add the weight (Number of samples associated with subtype/All samples):
    sen <- (conf_table[conf_table$Prediction == n & conf_table$Reference == n, 3]/
              sum(conf_table[conf_table$Reference == n, 3]))*(sum(conf_table[conf_table$Reference == n, 3])/sum(conf_table$Freq))
    
    # Calculate the Specificity (TN / (TN + FP)) and add the weight (Number of samples associated with subtype/All samples):
    spec <- (sum(conf_table[conf_table$Prediction != n & conf_table$Reference != n, 3])/
               sum(conf_table[conf_table$Reference != n, 3]))*(sum(conf_table[conf_table$Reference == n, 3])/sum(conf_table$Freq))
    
    # Add the results to a temporary table:
    metrics <- data.frame(subtype = n,
                          sen = sen,
                          spec = spec)
    
    # Add the results to the final table: 
    conf_metrics <- rbind(conf_metrics, metrics)
  }
  
  # Create a vector with the sensitivity and specificity results: 
  final_result <- c(sum(conf_metrics$sen), sum(conf_metrics$spec))  
  
  return <- final_result
  
}

# -------------------------> DATA LOAD: 

# Create the string list to open the different subtypes: 
cluster <- c("Neu", "OxA", "SNs")

# Classification file:
classification <- read.delim("New_clusters_all_samples.tsv")

# -----------------------> LINEAR DISCRIMINATION ANALYSIS:

  
  # Open the training dataset: 
  lda_training <- read.delim("Ens_vst_combat_AllGenes_Train.tsv", row.names = 1)
  lda_training <- as.data.frame(t(lda_training))
  
  # Get the relevant features for cluster selection:
  filenames <- list.files("models", pattern="*AllGenes_Blood.R", full.names=TRUE) 
  models <- lapply(filenames, readRDS)
  
  # Open each model and get the final feature table:
  candidates <- c()
  final_candidates <- data.frame()
  
  for (n in 1:length(models)) {
    # Open each model for each cluster:
    rf_model <- models[[n]]
    
    # Get the features for each model:
    features <- rownames(rf_model$importance)
    
    # Added to the final vector:
    candidates <- c(candidates, features)
    
    # Create a table with features for each model and add info of the cluster that they are relevant to:
    detailed_table <- data.frame(candidate = features,
                                 cluster = cluster[[n]])
    
    # Save all the features in a single table:
    final_candidates <- rbind(final_candidates, detailed_table)
    
  }
  
  # Save the list of features for table modification of blood: 
  feature_table <- data.frame(candidate = candidates)
  
  write.table(feature_table, 
              file = "../../../output/candidate_genes_Blood.tsv", 
              sep = "\t", 
              quote = FALSE,
              row.names = FALSE)
  
  # Subset the table to get only the training dataset and relevant genes:
  sub_training <- lda_training[grepl("^C", row.names(lda_training)), colnames(lda_training) %in% candidates]
  
  # Get the right cluster classification:
  sub_classification <- classification[grepl("^C", classification$Sample), colnames(classification) == "AllGenes" |
                                         colnames(classification) == "Sample"]
  colnames(sub_classification)[2] <- "Cluster" # Rename the cluster column
  
  # Get the cluster classification as a factor:
  sub_classification$Cluster <- factor(sub_classification$Cluster)
  
  # LDA Model:
  # Approach used in supervised ML to solve multiclass classification problems. 
  # Separation through data dimensionality reduction. Bayesian classification (
  # the probability of an event happening due to another one). LDA works by 
  # identifying a linear combination of features that separates or characteries
  # two or more dimension so that can be more easily classified.
  
  # Assumptions:
  # Normal distribution (vst dosent assure full normal distribution, but slightly
  # normalize the data).
  # Data is linearly separable.
  # Each class has the same covariance matrix. 
  
  lda_model <- train(sub_training, sub_classification$Cluster, method = "lda", trControl = 
                       trainControl("boot", number = 70))
  
  # Get confusion matrix:
  conf_lda <- as.data.frame(confusionMatrix(lda_model, "none")$table)
  
  # Calculate weighted sensitivity and specificity:
  # The result is a vector where sensitivy is safe as the firs value and 
  # specificity as the second one:
  conf_values <- weighted(conf_lda)
  
  # Plot model (Just for visualization): 
  sub_training$Cluster <- sub_classification$Cluster
  
  # Define colors:
  my_cols <- c("royalblue2", "green4", "firebrick3")
  cluster_colors <- my_cols[as.numeric(sub_training$Cluster)]
  
  # Get LDA results
  lda_result <- lda(Cluster ~ ., data = sub_training)
  lda_plot <- predict(lda_result, sub_training)
  
  # Create data frame for ggplot
  lda_df <- data.frame(LD1 = lda_plot$x[, 1], LD2 = lda_plot$x[, 2],
                       Cluster = sub_training$Cluster)
  
  write.table(lda_df, 
              file = "../../../output/LDA_table.tsv", 
              sep = "\t", 
              quote = FALSE,
              row.names = TRUE)
  
  # LDA Plot training:
  Lda_train_plot <- 
    ggplot(lda_df, aes(LD1, LD2, color = Cluster)) +
    geom_point(size = 5.5) +
    scale_color_manual(values = c("royalblue2", "green4","firebrick3")) +
    theme(axis.title = element_text(size = 60, hjust = 0.5),
          axis.text.x  =  element_text(size = 55, colour = "black"),
          axis.text.y  = element_text(size = 55, colour = "black"),
          legend.title = element_text(size = 50),
          legend.text = element_text(size = 45),
          legend.position = "right",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_line(colour = "black", linewidth = 2),
          axis.ticks.length = unit(0.6, "cm"),
          axis.line = element_line(colour = "black", linewidth = 2),
          plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 
  
  # Save the plot:
  png(file = "../../../output/LDA_plot.png", width = 1200, height = 1000)
  
  print(Lda_train_plot)
  
  dev.off()
  
  # ----------------> LINEAR DISCRIMINATION MODEL EVALUATION:
    
    # Open the evaluation table:
    lda_eval <- read.delim("Ens_vst_combat_AllGenes_Grima.tsv", row.names = 1)
    lda_eval <- as.data.frame(t(lda_eval))
    
    # Subset to only get relevant genes: 
    sub_eval <- lda_eval[, colnames(lda_eval) %in% candidates]
    sub_eval <- sub_eval[, order(colnames(sub_eval))]
    
    # Make the prediction of the new dataset:
    eval_pred <- predict(lda_model$finalModel, sub_eval)
    
    # Classification by prediction: 
    class_pred <- as.data.frame(predict(lda_model$finalModel, sub_eval))
    class_pred <- class_pred[, 1, drop = FALSE]
    colnames(class_pred)[1] <- "prediction"
    
    # Calculate posterior probability:
    posterior <- as.data.frame(eval_pred$posterior)
    
    # Get the predicted subtypes for every sample: 
    lda_assign <- posterior %>%
      mutate(Sample = rownames(.)) %>%
      pivot_longer(-Sample, names_to = "Cluster", values_to = "Prob") %>%
      group_by(Sample) %>%
      slice_max(Prob, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    # Save output table:
    write.table(lda_assign, 
                file = "../../../output/LDA_cluster_probability_blood.tsv",  
                sep = "\t", 
                quote = FALSE,
                row.names = FALSE)
    
    # Keep samples assigned with >= 0.85 probability:
    high_pos <- posterior[apply(posterior, 1, max) >= 0.85, ]
    
    # Do it for each cluster:
    high_neu <- lda_assign[lda_assign$Cluster == "Neu" & lda_assign$Prob >= 0.85, ]
    high_oxa <- lda_assign[lda_assign$Cluster == "OxA" & lda_assign$Prob >= 0.85, ]
    high_sns <- lda_assign[lda_assign$Cluster == "SNs" & lda_assign$Prob >= 0.85, ]
    
    # LDA Table: 
    lda_table <- data.frame(features = length(candidates),
                            accuracy = (lda_model$results$Accuracy)*100,
                            sensitivity = conf_values[1],
                            specificity = conf_values[2],
                            porcentage_high_prob = (nrow(high_pos)/nrow(posterior))*100,
                            porcentage_high_prob_Neu = (nrow(high_neu)/nrow(lda_assign[lda_assign$Cluster == "Neu", ]))*100,
                            porcentage_high_prob_Oxa = (nrow(high_oxa)/nrow(lda_assign[lda_assign$Cluster == "OxA", ]))*100,
                            porcentage_high_prob_SNs = (nrow(high_sns)/nrow(lda_assign[lda_assign$Cluster == "SNs", ]))*100)
    
    # Create evaluation plot:
    ld_data <- eval_pred$x
    
    # Create an extra column with classification: 
    max_cols <- apply(posterior, 1, function(x) colnames(posterior)[which.max(x)])
    max_cols <- data.frame(row.names = names(max_cols),
                           Prediction = max_cols)
    
    # Add the values to the table:
    ld_data <- merge(ld_data, max_cols, by = "row.names")
    ld_data$Prediction <- factor(ld_data$Prediction)
    
    write.table(ld_data, 
                file = "../../../output/Eval_plot_in_Blood.tsv",
                sep = "\t", 
                quote = FALSE,
                row.names = FALSE)
    
    # Create the plot:
    eval_plot <- 
      ggplot(ld_data, aes(LD1, LD2, color = Prediction)) +
      geom_point(size = 5.5) +
      scale_color_manual(values = c("royalblue2", "green4","firebrick3")) +
      theme(axis.title = element_text(size = 60, hjust = 0.5),
            axis.text.x  =  element_text(size = 55, colour = "black"),
            axis.text.y  = element_text(size = 55, colour = "black"),
            legend.title = element_text(size = 50),
            legend.text = element_text(size = 45),
            legend.position = "right",
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_line(colour = "black", linewidth = 2),
            axis.ticks.length = unit(0.6, "cm"),
            axis.line = element_line(colour = "black", linewidth = 2),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 
    
    # Save the plot:
    png(file = "../../../output/Eval_plot_in_Blood.png", width = 1200, height = 1000)
    
    print(eval_plot)
    
    dev.off()
  
# -----------------------> FINAL RESULTS:
  
# Get final table with right format:

# Round percentages to exact numbers:
lda_table[, c(2)] <- round(lda_table[, c(2)], digits = 0)

# Add 2 decimals to sensitivity and specificity:
lda_table[, c(4:8)] <- round(lda_table[, c(4:8)], digits = 2)

# Output the table: 

write.table(lda_table, 
            file = "../../../output/LDA_blood_evaluation_results.tsv", 
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE)




      
      
    
    
    
    
    
  
  
  

