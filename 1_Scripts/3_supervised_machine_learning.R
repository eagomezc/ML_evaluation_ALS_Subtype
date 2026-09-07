#---------- ALS Subtype Evaluation: Supervised Machine Learning --------------#

# This script takes transcriptomic data to create machine learning models able 
# to classify ALS samples into different subtypes:

#-----------------------------> LIBRARY LOAD:

# To only use when you're working in the cluster:
#.libPaths("/scratch/users/k2584930/software/R/4.3/")

if (!require('randomForest')) install.packages('randomForest'); library('randomForest')
if (!require('caret')) install.packages('caret'); library('caret')
if (!require('glmnet')) install.packages('glmnet'); library('glmnet')
if (!require('xgboost')) install.packages('xgboost'); library('xgboost')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('tidyr')) install.packages('tidyr'); library('tidyr')


set.seed(42) # To get same results even with the random part.
options(digits = 3) # To get only part of the decimals. 

# Set directory to be in the source file location (Rstudio):
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Set directory to be in the source file location (Cluster):
#this_dir <- function(directory)
 #setwd( file.path(getwd(), directory) )

# -------------------------> DATA LOAD: 

# Upload VST counts with gene symbols and cluster information:
train <- read.delim("training_all_genes_cluster.tsv", row.names = 1)

# Upload available meta data to do a ML model with them. This can be modified 
# for a dataset with proper clinical data: 
bb_meta <- read.delim("bb_meta_cluster.tsv", row.names = 1)

# Since from the exploratory analysis the data has already be normalized, no 
# preprocessing steps are neccesary.


# -----------------------> MACHINE LEARNING MODELS:

# Create final tables with the info of all the models created:
final_table <- data.frame() 
  
for (j in unique(train$Cluster)) {
    
  # Get the right format of One vs Rest for each of the clusters:
  train_data <- train
  train_data[train_data$Cluster != j, ]$Cluster <- "Uncluster"
    
  # Add the classification variable to the data frame (One and Rest):
  
  # Getting the explanatory (x) and response (y) variable. By explanatory, 
  # it means all the data that can explain why a patient is or not in an 
  # specific cluster (genes) and the response variable is if the patients belongs
  # to one cluster or not. You then create a formula where the response variable
  # is explain in terms of the explanatory variables:
    
  # Cluster ~ everything else 
    
  # Convert in factors the binary operator:
  train_data$Cluster <- factor(train_data$Cluster)
  
  # Make sure that column names don't represent a problem to the ML models:
  names(train_data) <- make.names(names(train_data))
     
  #--------> RANDOM forest (randomForest R): 
     
  # In Random Forests the idea is to decorrelate the several trees which are 
  # generated on the different bootstrapped samples from training Data and then 
  # reduce the variance in the trees by averaging them.
     
  # Averaging the trees also improve the performance of decision trees on Test 
  # Set and eventually avoid overfitting.
     
  # The idea is to build lots of trees in such a way to make the correlation 
  # between the trees smaller.
     
  # BEST MTRY:
  # mtry is the number of variables available for splitting at each tree node. 
  # Random Forest creates several trees, each one using different variables to 
  # create the best version of it. With mtry we can define how many variables 
  # the data is split to create the different trees.
  # More: https://stats.stackexchange.com/questions/102867/random-forest-mtry-question
     
  # In this case we create a loop to define which mtry is the best one for our models. 
     
  oob_error <- double(ncol(train_data) -1) # Define the number of variables.
                                              # - 1 because the last column is cluster.
     
  # Loop to identify the best mtry:
  for (mtry in 1:(ncol(train_data) - 1)) {
       
  # Run the rf model:
  rf_train_data <- randomForest(Cluster ~ ., data = train_data, mtry = mtry, ntree = 1000)
       
  # Save the mtry value and it's associated accuracy
  oob_error[mtry] <- 100 - ((rf_train_data$err.rate[1000])*100)
  
  }
     
  # Define the best mtry according to the best prediction value. 
  final_mtry <- which.max(oob_error)
     
  # Run the model again with the right mtry value. 
  rf_train_final <- randomForest(Cluster ~ ., data = train_data, mtry = final_mtry, 
                                          importance = TRUE, ntree = 1000)
     
  # Save the model: 
  saveRDS(rf_train_final, file = paste("../../output/models/rf_", j, ".R", sep = ""),  
          ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
     
  # Get the confusion matrix of the model, sensitivity and specificity: 
  confusion_rf <- as.data.frame(rf_train_final$confusion)
  confusion_rf[is.na(confusion_rf)] <- 0
     
  # Calculates sensitivity and specificity.
  specificity_rf <- confusion_rf[2, 2]/(confusion_rf[2, 2] + confusion_rf[2, 1])
  sensitivity_rf <- confusion_rf[1, 1]/(confusion_rf[1, 1] + confusion_rf[1, 2])
     
  # Final table for random forest:
  no_oob_error_table <- data.frame(machine_learning = "Random Forest",
                                    cluster = j, 
                                    features = ncol(train_data)-1,
                                    percentage_accuracy = 100 - ((rf_train_final$err.rate[1000])*100),
                                    sensitivity = sensitivity_rf,
                                    specificity = specificity_rf,
                                    TP_per = confusion_rf[1, 1]/sum(confusion_rf[,-3])*100,
                                    FP_per = confusion_rf[2, 1]/sum(confusion_rf[,-3])*100,
                                    TN_per = confusion_rf[2, 2]/sum(confusion_rf[,-3])*100,
                                    FN_per = confusion_rf[1, 2]/sum(confusion_rf[,-3])*100,
                                    stringsAsFactors = FALSE)
     
  final_table <- rbind(final_table, no_oob_error_table)
     
  # Number of trees plot: 
  # Creates a plot of the ntree tunning parameter:
     
  tree_table <- data.frame(Ensemble = c(1:1000),
                           OBB = rf_train_final$err.rate[, 1],
                           AvgAcc = 100 - ((rf_train_final$err.rate[, 1])*100))
     
  # Tree plot:
  ntree_plot <- 
    ggplot(data = tree_table, aes(x = Ensemble, y = AvgAcc)) + 
    geom_line(linetype = "solid", size = 1.2, color = "navyblue") +
    scale_x_continuous(name = "Trees") +
    scale_y_continuous(name = "Average Percent Accuracy", limits=c(min(tree_table$AvgAcc)-10, max(tree_table$AvgAcc)+10)) +
    theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
          axis.title.y = element_text(size = 30, colour = "black"),
          axis.title.x = element_text(size = 30, colour = "black"),
          axis.text.x  =  element_text(size = 30, colour = "black"), # Put color to the labels
          axis.text.y  = element_text(size = 30, colour = "black"), # Put color to the labels
          axis.line = element_line(colour = 'black', linewidth = 1.5), # Color and thickness of axis
          axis.ticks = element_line(colour = "black", linewidth = 1.5), # Color and thickness of every axis sep. 
          legend.position = ("none"),
          panel.background = element_rect(fill = "white"),
          axis.ticks.length = unit(0.4, "cm"))
     
  write.table(tree_table, 
              file = paste("../../output/tuning_plot/rf_trees_", j, ".tsv", sep = ""), 
              sep = "\t", 
              quote = FALSE,
              row.names = FALSE)
     
  png(file = paste("../../output/3_machine_learning/tuning_plot/rf_trees_", j, ".png", sep = ""),
         width = 1200, height = 1000)
     
  print(ntree_plot)
     
  dev.off()
     
  # Generate the IMPORTANCE FIGURE: 
  # This figure identifies the most relevant features for the classification in RF.
     
  rf_imp_list <- as.data.frame(importance(rf_train_final, type = 1))
  rf_imp_list$features <- rownames(rf_imp_list)
     
  # Order the data in increasing order based on the Mean Decrease Accuracy and 
  # creater a factor column (this with the purpose of get the features in decreasing
  # order in the figure): 
     
  rf_imp_list$features <- factor(rf_imp_list$features,
                                 levels = rf_imp_list$features[order(rf_imp_list$MeanDecreaseAccuracy)])
     
  # Generates the IMPORTANCE plot:
  importance_plot <- 
    ggplot(data = rf_imp_list, mapping = aes(x = features, y = MeanDecreaseAccuracy)) +
    geom_point(size = 8) +
    scale_y_continuous(name = "Mean Decrease Accuracy") +
    labs(x = "Features") +
    coord_flip() +
    theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
          axis.title = element_text(size = 30),
          axis.text.x  = element_text(size = 20, hjust = 0.5, colour = "black"),
          axis.text.y  = element_text(size = 20, hjust = 1, colour = "black"),
          aspect.ratio = 2/1, 
          panel.background = element_rect(fill = "snow"),
          panel.grid.major.y = element_line(linewidth = 0.5, linetype = 2, colour = "black")) 
     
  # Generates the IMPORTANCE plot but with only the top 50 genes: 
  rf_imp_top <- head(rf_imp_list, 50)
     
  impor_top_plot <- 
    ggplot(data = rf_imp_top, mapping = aes(x = features, y = MeanDecreaseAccuracy)) +
    geom_point(size = 8) +
    scale_y_continuous(name = "Mean Decrease Accuracy") +
    labs(x = "Features") +
    coord_flip() +
    theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
          axis.title = element_text(size = 30),
          axis.text.x  = element_text(size = 20, hjust = 0.5, colour = "black"),
          axis.text.y  = element_text(size = 20, hjust = 1, colour = "black"),
          aspect.ratio = 2/1, 
          panel.background = element_rect(fill = "snow"),
          panel.grid.major.y = element_line(linewidth = 0.5, linetype = 2, colour = "black")) 
     
   write.table(rf_imp_list, 
               file = paste("../../output/importance_plot/rf_imp_", 
                             j, ".tsv", sep = ""), 
               sep = "\t", 
               quote = FALSE,
               row.names = FALSE)
     
   png(file = paste("../../output/importance_plot/rf_imp_", j, ".png", sep = ""),
         width = 900, height = 1100)
     
   print(importance_plot)
     
   dev.off()
     
   png(file = paste("../../output/importance_plot/rf_imp_top_50_", j, ".png", sep = ""),
         width = 900, height = 1100)
     
   print(impor_top_plot)
     
   dev.off()
     
   # Generate the CONFUSION Plot: 
     
   # Get the real classification:
   confusion_rf$InitClass <- row.names(confusion_rf)
   # Delete class error:
   confusion_rf$class.error <- NULL
   #Get the table in the right format: 
   confusion_rf <- confusion_rf %>% pivot_longer(cols = colnames(confusion_rf)[-3],
                                      names_to = "PredClass",
                                      values_to = "Overall")
   # Get the percentages: 
   confusion_rf$Percentage <- 0 # Create a numeric column
   # Get the TP and FN:
   confusion_rf[c(1:2), ]$Percentage <- round(confusion_rf[c(1:2), ]$Overall/
                                                  sum(confusion_rf[c(1:2), ]$Overall)*100, digits = 0) 
     
   # Get the TN and FP:
   confusion_rf[c(3:4), ]$Percentage <- round(confusion_rf[c(3:4), ]$Overall/
                                                  sum(confusion_rf[c(3:4), ]$Overall)*100, digits = 0) 
     
   # Generate the plot:
   confusio_plot <- 
       ggplot(data = confusion_rf, mapping = aes(x = InitClass, y = Percentage, fill = PredClass)) +
       geom_bar(stat = "identity", colour = "black") +
       geom_text(data = confusion_rf, aes(x = InitClass, y = Percentage, label = paste(Percentage,"%",sep="")), 
                 size = 12, position = "stack", colour = c("white", "chartreuse3", "blue3", "black"), 
                 alpha = c(1, 0, 0, 1), vjust = 1.2) +
       scale_fill_manual(values = alpha(c("blue3", "chartreuse3"), 0.7), labels = c(j, "Rest")) +
       scale_y_continuous(name = "Class Predictions (%)", breaks = seq(from = 0, to = 100, by = 20),
                          expand = c(0, 1)) +
       scale_x_discrete(labels = c(j, "Rest")) + 
       labs(x = "Classes") +
       coord_cartesian(ylim = c(1, 100)) +
       theme(plot.title = element_text(size = 40, hjust = 0.5),
             axis.title = element_text(size = 40),
             axis.text.x  = element_text(size = 38, hjust = 0.5, colour = "black"),
             axis.text.y  = element_text(size = 40, hjust = 0.5, colour = "black"),
             legend.title = element_blank(),
             legend.key.size = unit(1.3, "cm"), 
             legend.text  = element_text(size = 40),
             axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
             axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
             panel.background = element_rect(fill = "white"),
             axis.ticks.length = unit(0.4, "cm"))
   
   write.table(confusion_rf, 
                 file = paste("../../output/confusion_plot/rf_conf_", j, ".tsv", sep = ""), 
                 sep = "\t", 
                 quote = FALSE,
                 row.names = FALSE)
     
   png(file = paste("../../output/confusion_plot/rf_conf_", j, ".png", sep = ""),
         width = 900, height = 1100)
     
   print(confusio_plot)
     
   dev.off()
     
  
     #-----> ELASTIC Net Regression (caret R):
     
     # Get the explanatory variables as a matrix:
     explanatory <- data.matrix(train_data)[, -(ncol(train_data))]
     
     # LASSO Analysis:
     
     # x: matrix of predictor variables
     # y: the response or outcome variable, which is a binary variable.
     # family: the response type. Use "binomial" for a binary outcome variable
     # alpha: the elasticnet mixing parameter. Allowed values include:
     #   "1": for lasso regression
     #   "0": for ridge regression
     #    A value between 0 and 1 (say 0.3) for elastic net regression.
     # lamba: a numeric value defining the amount of shrinkage. Should be specify by analyst.
     
     # The LASSO method has some limitations:
     # In small-n-large-p dataset the LASSO selects at most n variables before 
     # it saturates. If there are grouped variables (highly correlated between each other) 
     # LASSO tends to select one variable from each group ignoring the others. 
     # Elastic Net overcomes LASSO limitations using a combination of LASSO and 
     # Ridge Regression methods (where no variable is ever eliminated, just reduced).
     
     # alpha: the elasticnet mixing parameter. Allowed values include:
     # A value between 0 and 1 for elastic net regression. The model uses bootstraping 
     # to find the best alpha value. 
     
     # lamba: a numeric value defining the amount of shrinkage. The model uses 
     # bootstraping to find the best alpha value. 
     
     # Define the best alpha and lamba value. The best lambda for your data, can be 
     # defined as the lambda that minimize the cross-validation prediction error rate.    
     
     # Fit predictive models over different tuning parameters:
     model_net <- train(explanatory, train_data$Cluster, method = "glmnet", 
                        trControl = trainControl("boot", number = 70))
     
     # NOTE: As what I understand from the documentation, the function train, applies
     # boostrapping for tunning parameters, but also for calculation of accuracy and 
     # other values, when the best tunning parameters has been identified (this to
     # reduce the variance). So, in this case, model_net, has already the best model
     # wit the applied boostrapping. 
     
     # Save the model: 
     saveRDS(model_net, 
             file = paste("../../output/models/net_", j, ".R", sep = ""),  
             ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
     
     # Get the confusion matrix of the model:
     # Since the best model is created with boostraping (sampling with replacement)
     # the confusion matrix generates more results (for each resampling bassically),
     # that's the reason you dont see 112 sample results, but thousands of them.
     # However, that confusion matrix works as a normal one and you can get stats 
     # from it.
     conf_net <- as.data.frame(confusionMatrix(model_net, "none")$table)
     
     # Final Elastic net model table:
     net_table <- data.frame(machine_learning = "Elastic net",
                             cluster = j, 
                             features = ncol(train_data)-1,
                             percentage_accuracy = max(model_net$results$Accuracy)*100,
                             sensitivity = conf_net[1, 3]/(conf_net[1, 3] + conf_net[2, 3]),
                             specificity = conf_net[4, 3]/(conf_net[4, 3] + conf_net[3, 3]),
                             TP_per = conf_net[1, 3]/(sum(conf_net$Freq))*100,
                             FP_per = conf_net[3, 3]/(sum(conf_net$Freq))*100,
                             TN_per = conf_net[4, 3]/(sum(conf_net$Freq))*100,
                             FN_per = conf_net[2, 3]/(sum(conf_net$Freq))*100,
                             stringsAsFactors = FALSE)
     
     final_table <- rbind(final_table, net_table)
     
     # Parameter tuning figure: 
     # Make a figure of the tunning for alpha and lambda. 
     
    scaleFUN <- function(x) sprintf("%.2f", x)
     
    boot_net_plot <- 
       ggplot(model_net, highlight = TRUE) + 
       scale_x_continuous(name = expression(paste("Alpha (", alpha, ")", sep = ""))) +
       scale_y_continuous(name = "Average Accuracy", labels= scaleFUN) +
       scale_color_manual(values = c("darkorchid3", "orangered1", "chartreuse3")) +
       scale_shape_manual(values=c(16, 16, 16)) +
       labs(color = expression(paste("Lambda (", lambda, ")", sep = ""))) +
       guides(shape = 'none') +
       theme(plot.title = element_text(size = 40, colour = "black", hjust = 0.5),
             axis.title.y = element_text(size = 40, colour = "black"),
             axis.title.x = element_text(size = 40, colour = "black"),
             axis.text.x  =  element_text(size = 40, colour = "black"), # Put color to the labels
             axis.text.y  = element_text(size = 40, colour = "black"), # Put color to the labels
             axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
             axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
             legend.title = element_text(size = 40),
             legend.text  = element_text(size = 35),
             legend.key = element_rect(fill = NA),
             legend.key.size = unit(1.3, "cm"), 
             panel.background = element_rect(fill = "white"),
             axis.ticks.length = unit(0.4, "cm"))
     
    write.table(model_net$results, 
                file = paste("../../output/tuning_plot/net_boot_", j, ".tsv", sep = ""), 
                sep = "\t", 
                quote = FALSE,
                row.names = FALSE)
    
    png(file = paste("../../output/tuning_plot/net_boot_", j, ".png", sep = ""),
        width = 1200, height = 1000)
    
    print(boot_net_plot)
    
    dev.off()
     
     # Generate the IMPORTANCE FIGURE: 
     # This figure identifies the most relevant features for the classification in Elastic Net.
     
     # We use the function VarImp to generates the coefficient associated with the
     # Elastic Net model (the best one, since train function is designed to tune and
     # generate the best model - I TRIED BY TUNING THE PARAMETERS MANNUALY AND I GET
     # SAME RESULTS, see LASSSO_Francesco.R script). Values higher than 0 indicates
     # contribution to the model. The highest the value, the more contribution. 
     
     # There is an option to scale the coefficient from 0 to 100 and also displays
     # the coeficient itself. 
     net_imp_list <- merge(varImp(model_net)$importance, 
                           varImp(model_net, scale = FALSE)$importance, 
                           by = "row.names")
     
     # Get the table in the right format: 
     colnames(net_imp_list) <- c("features", "Scale", "Coefficient")
     net_imp_list <- net_imp_list %>% pivot_longer(cols = colnames(net_imp_list)[-1],
                                                   names_to = "type",
                                                   values_to = "values")
     
     # Get the features in the right order:
     net_imp_list$features <- factor(net_imp_list$features,
                                     levels = unique(net_imp_list$features[order(net_imp_list$type,
                                                                          net_imp_list$values)]))
     
     # Generates the IMPORTANCE plot:
     net_imp_plot <- 
       ggplot(data = net_imp_list, mapping = aes(x = features, y = values)) +
       geom_point(size = 8) +
       labs(x = "Features") +
       coord_flip() +
       facet_wrap("type", scales = "free_x", ncol = 2, nrow = 1, strip.position = "bottom") +
       theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
             axis.title = element_text(size = 30),
             axis.title.x = element_blank(), 
             axis.text.x  = element_text(size = 20, hjust = 0.5, colour = "black"),
             axis.text.y  = element_text(size = 20, hjust = 1, colour = "black"),
             aspect.ratio = 2/1, 
             panel.background = element_rect(fill = "snow"),
             strip.background = element_blank(),
             strip.placement = "outside",
             strip.text = element_text(size = 25, hjust = 0.5, colour = "black"),
             panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black")) 
     
     # Generates importance plot top 50: 
     net_imp_top <- head(net_imp_list[order(net_imp_list$type, net_imp_list$values, 
                                            decreasing = TRUE), ], 50)
     
     net_imp_top_plot <- 
       ggplot(data = net_imp_top, mapping = aes(x = features, y = values)) +
       geom_point(size = 8) +
       labs(x = "Features") +
       coord_flip() +
       facet_wrap("type", scales = "free_x", ncol = 2, nrow = 1, strip.position = "bottom") +
       theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
             axis.title = element_text(size = 30),
             axis.title.x = element_blank(), 
             axis.text.x  = element_text(size = 20, hjust = 0.5, colour = "black"),
             axis.text.y  = element_text(size = 20, hjust = 1, colour = "black"),
             aspect.ratio = 2/1, 
             panel.background = element_rect(fill = "snow"),
             strip.background = element_blank(),
             strip.placement = "outside",
             strip.text = element_text(size = 25, hjust = 0.5, colour = "black"),
             panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black")) 
     
     write.table(net_imp_list, 
                 file = paste("../../output/importance_plot/net_imp_", j, ".tsv", sep = ""), 
                 sep = "\t", 
                 quote = FALSE,
                 row.names = FALSE)
     
     png(file = paste("../../output/importance_plot/net_imp_", j, ".png", sep = ""),
         width = 1000, height = 1100)
     
     print(net_imp_plot)
     
     dev.off()
     
     png(file = paste("../../output/importance_plot/net_imp_top_50_", j, ".png", sep = ""),
         width = 1000, height = 1100)
     
     print(net_imp_top_plot)
     
     dev.off()
     
     # Generate the CONFUSION FIGURE: 
     # Get the table in the right order:
     conf_net <- conf_net[, c(2,1,3)]
     # Get colnames in the right 
     colnames(conf_net) <- c("InitClass", "PredClass", "Overall")
     # Get the percentages: 
     conf_net$Percentage <- 0 # Create a numeric column
     # Get the TP and FN:
     conf_net[c(1:2), ]$Percentage <- round(conf_net[c(1:2), ]$Overall/
                                                  sum(conf_net[c(1:2), ]$Overall)*100, digits = 0) 
     
     # Get the TN and FP:
     conf_net[c(3:4), ]$Percentage <- round(conf_net[c(3:4), ]$Overall/
                                                  sum(conf_net[c(3:4), ]$Overall)*100, digits = 0) 
     
     # Generate the plot:
     net_confusio_plot <- 
       ggplot(data = conf_net, mapping = aes(x = InitClass, y = Percentage, fill = PredClass)) +
       geom_bar(stat = "identity", colour = "black") +
       geom_text(data = conf_net, aes(x = InitClass, y = Percentage, label = paste(Percentage,"%",sep="")), 
                 size = 12, position = "stack", colour = c("white", "chartreuse3", "blue3", "black"), 
                 alpha = c(1, 0, 0, 1), vjust = 1.2) +
       scale_fill_manual(values = alpha(c("blue3", "chartreuse3"), 0.7), labels = c(j, "Rest")) +
       scale_y_continuous(name = "Class Predictions (%)", breaks = seq(from = 0, to = 100, by = 20),
                          expand = c(0, 1)) +
       scale_x_discrete(labels = c(j, "Rest")) + 
       labs(x = "Classes") +
       coord_cartesian(ylim = c(1, 100)) +
       theme(plot.title = element_text(size = 40, hjust = 0.5),
             axis.title = element_text(size = 40),
             axis.text.x  = element_text(size = 38, hjust = 0.5, colour = "black"),
             axis.text.y  = element_text(size = 40, hjust = 0.5, colour = "black"),
             legend.title = element_blank(),
             legend.key.size = unit(1.3, "cm"), 
             legend.text  = element_text(size = 40),
             axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
             axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
             panel.background = element_rect(fill = "white"),
             axis.ticks.length = unit(0.4, "cm"))
     
     write.table(conf_net, 
                 file = paste("../../output/confusion_plot/net_conf_", j, ".tsv", sep = ""), 
                 sep = "\t", 
                 quote = FALSE,
                 row.names = FALSE)
     
     png(file = paste("../../output/confusion_plot/net_conf_", j, ".png", sep = ""),
         width = 900, height = 1100)
     
     print(net_confusio_plot)
     
     dev.off()
     
     #-----> EXTREME gradient boosting (XGBoost R):
     
     # XGBoost also works with decision trees, but different than random forest,
     # where the final model is the consensus of weak decision trees, XGBoost builds
     # a tree after another, improving the previous one. It's based on the bias-variance
     # tradeoff where you keep a balance between a simple but good a predicting model.
     
     # Parameters tunning:
     # One big issue with XGBoost is that request a lot of parameter tunning (a lot of),
     # which is very annoying to do. The parameters are the following:
     
     # eta: scale the contribution of each tree by a factor of 0 < eta < 1. Used to
     #      prevent overfitting by making the boosting process more conservative. Lower values
     #      for eta implies larger value for nrounds (trees) and more robust models to overfitting. 
     #      Default 0.3. [0,1]
     # gamma: minimum lost reduction required to make a further partition on a leaf node
     #       of the tree. The larger gamma is, the more conservative the model will be. 
     #     Default 0. [0, inf] For large features [0, 1]
     # max_depth: Maximum deph of a tree. Increasing this value will make the model more
     #     complex and more likely to overfit. 0 indicates no limit on depth. 
     #     Default 6. [0, inf]. For large features [0, 6]
     # min_child_weight: Minimum sum of instance weight needed in a child. If the 
     #                   tree new partion weight is lower than this value, so there 
     #                   is no partition. The larger, the more conservative. 
     #                   Default 1. [0, inf]. [5-20] Prevent overfitting.
     # subsample: Use portion of the data to create the model. Keep the default of 1.
     #            Since you dont want to use less data that the low data we already have.
     # colsample_bytree: Basically the mtry from randomforest. 
     #                   Default: 1. [0-1]. Large features [0.2, 0.6].
     # lambda: L2 regularization term on weights. Increasing this value will make the
     #         model more conservative. 
     #        Default: 1. [0, inf]. For large features [1-10]. 
     # alpha: L1 regularization term on weights. Increasing this value will make the 
     #        model more conservative. 
     #        Default: 0. [0,inf]. [0.5,10] to lead some variable to zero. 
     
     # Define the parameter grid where you live all the possible parameter options:
     xgb_grid <- expand.grid(nrounds = 10, subsample = 1, eta = c(0.1, 0.3, 0.5), 
                             gamma = c(0, 1, 1.5), max_depth = c(2, 4, 6), 
                             min_child_weight = c(1, 10, 20),
                             colsample_bytree = c(0.2, 0.6, 1))
                             
     # Fit predictive models over different tuning parameters:
     # Build the multiple models:
     xgb_model <- train(explanatory, train_data$Cluster, method = "xgbTree", 
                        trControl = trainControl("boot", number = 70, allowParallel = TRUE),
                        tuneGrid = xgb_grid)
     
     # Get the confusion matrix of the model:
     # Since the best model is created with boostraping (sampling with replacement)
     # the confusion matrix generates more results (for each resampling bassically),
     # that's the reason you dont see 112 sample results, but thousands of them.
     # However, that confusion matrix works as a normal one and you can get stats 
     # from it. 
     conf_xbg <- as.data.frame(confusionMatrix(xgb_model, "none")$table)
     
     # Save the model: 
     saveRDS(xgb_model, 
             file = paste("../../output/models/xgb_", j, ".R", sep = ""),  
             ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
     
     # Final XGB table:
     xgb_table <- data.frame(machine_learning = "Gradient boost trees",
                             cluster = j, 
                             features = ncol(train_data)-1,
                             percentage_accuracy = max(xgb_model$results$Accuracy)*100,
                             sensitivity = conf_xbg[1, 3]/(conf_xbg[1, 3] + conf_xbg[2, 3]),
                             specificity = conf_xbg[4, 3]/(conf_xbg[4, 3] + conf_xbg[3, 3]),
                             TP_per = conf_xbg[1, 3]/(sum(conf_xbg$Freq))*100,
                             FP_per = conf_xbg[3, 3]/(sum(conf_xbg$Freq))*100,
                             TN_per = conf_xbg[4, 3]/(sum(conf_xbg$Freq))*100,
                             FN_per = conf_xbg[2, 3]/(sum(conf_xbg$Freq))*100,
                             stringsAsFactors = FALSE)
     
     final_table <- rbind(final_table, xgb_table)
     
     # Parameter tuning figure: 
     # Make a figure of the tunning for alpha and lambda. 
     
     # Get a table with top 10 models:
     xgb_results <- xgb_model$results
     colnames(xgb_results)[c(1:5)] <- c("Eta", "MaxDepth", "Gamma", "ColsampleByTree",
                                        "MinChildWeight")
     xgb_results_10 <- xgb_results[order(xgb_results$Accuracy, decreasing = TRUE), ]
     xgb_results_10 <- head(xgb_results_10, 10)
     
     # Generate a list with tunning parameters that changes:
     n_unique <- sapply(xgb_results_10[, c(1:5)], function(x) length(unique(x))) # Exclude constant variables
     n_unique <- sort(n_unique, decreasing = TRUE) # Sort for more diverse tunning values
     names_un <- names(n_unique) # Get the tunning names
     
     # To have Accuracy with low digits:
     scaleFUN <- function(x) sprintf("%.2f", x)
     
     boot_xgb_plot <- 
       ggplot(xgb_results, aes(x = .data[[names_un[1]]], y = Accuracy, 
                                     color = as.factor(.data[[names_un[2]]]))) +
         geom_point() +
         facet_grid(as.formula(paste(names(n_unique)[3], "~", names(n_unique)[4])),
                    scales = "free_y", labeller = label_both) + 
       scale_x_continuous(name = names(n_unique[1])) +
       scale_y_continuous(name = "Average Accuracy", labels= scaleFUN) +
       scale_color_manual(values = c("darkorchid3", "orangered1", "chartreuse3")) +
       scale_shape_manual(values=c(16, 16, 16)) +
       labs(color = names(n_unique[2])) +
       guides(shape = 'none') +
       theme(plot.title = element_text(size = 40, colour = "black", hjust = 0.5),
             axis.title.y = element_text(size = 30, colour = "black"),
             axis.title.x = element_text(size = 30, colour = "black"),
             axis.text.x  =  element_text(size = 20, colour = "black"), # Put color to the labels
             axis.text.y  = element_text(size = 20, colour = "black"), # Put color to the labels
             axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
             axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
             legend.title = element_text(size = 25),
             legend.text  = element_text(size = 20),
             legend.key = element_rect(fill = NA),
             legend.key.size = unit(1.3, "cm"), 
             axis.ticks.length = unit(0.4, "cm"),
             panel.background = element_rect(fill = "snow"),
             strip.background = element_rect(fill = "lightblue"),
             strip.text = element_text(size = 15, hjust = 0.5, colour = "black"))
     
     write.table(xgb_results, 
                 file = paste("../../output/tuning_plot/xgb_boot_", j, ".tsv", sep = ""), 
                 sep = "\t", 
                 quote = FALSE,
                 row.names = FALSE)
     
     png(file = paste("../../output/tuning_plot/xgb_boot_", j, ".png", sep = ""),
         width = 1200, height = 1000)
     
     print(boot_xgb_plot)
     
     dev.off()
     
     # Generate the IMPORTANCE FIGURE: 
     # So, for whatever reason, I cant use VarImp in a xgb model because the features
     # got lost in the model. Which is incredibly annoying. However, the boostrapping
     # has been already being done, so I just need to generate the best model based
     # on the bestTune parameters. 
     
     # Obtain the best model after the boostrapping:
     # Since the importance feature of xgboost dosn't  work in the "caret" object, 
     # I will create a regular model using the best parameters for the boostraping 
     # and see if it works.
     
     # Create a list with the parameters of the best model:
     parameters = list(max_depth = xgb_model$bestTune$max_depth,
                       eta = xgb_model$bestTune$eta,
                       gamma = xgb_model$bestTune$gamma,
                       colsample_bytree = xgb_model$bestTune$colsample_bytree,
                       min_child_weight = xgb_model$bestTune$min_child_weight,
                       subsample = xgb_model$bestTune$subsample)
     
     best_xgb <- xgboost(explanatory, train_data$Cluster, params = parameters, verbose = FALSE,
                         nrounds = 1000, eval_metric = "error")
     
     # Generate the importance list:
     # I used the xgb.importance feature from the xgboost package. 
     # It generates a table with the following columns:
     # Feature, Gain (contribution of each feature to the model. Is the improvement
     # in accuracy brougth by a feature to the branches it is on), Cover (number of
     # observations related to this feature - # of observation to be classified and 
     # affected by this feature) and Frequency (% representing the relative
     # number of times a feature have been used in trees).
     
     xgb_imp_list <- xgb.importance(best_xgb$feature_names, best_xgb)
     xgb_imp_list$Gain <- xgb_imp_list$Gain*100
     
     # Get the right order for the plot:
     xgb_imp_list$Feature <- factor(xgb_imp_list$Feature,
                            levels = unique(xgb_imp_list$Feature[order(xgb_imp_list$Gain)]))
     
     # Generates the IMPORTANCE plot:
     xgb_imp_plot <- 
       ggplot(data = xgb_imp_list, mapping = aes(x = Feature, y = Gain)) +
       geom_point(size = 8) +
       labs(x = "Features") +
       coord_flip() + 
         theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
               axis.title = element_text(size = 30),
               axis.text.x  = element_text(size = 20, hjust = 0.5, colour = "black"),
               axis.text.y  = element_text(size = 20, hjust = 1, colour = "black"),
               aspect.ratio = 2/1, 
               panel.background = element_rect(fill = "snow"),
               panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black")) 
     
     # Generates the importance plot with the TOP 50 genes:
     xgb_imp_top <- head(xgb_imp_list, 50)
     
     xgb_imp_top_plot <- 
       ggplot(data = xgb_imp_top, mapping = aes(x = Feature, y = Gain)) +
       geom_point(size = 8) +
       labs(x = "Features") +
       coord_flip() + 
       theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
             axis.title = element_text(size = 30),
             axis.text.x  = element_text(size = 20, hjust = 0.5, colour = "black"),
             axis.text.y  = element_text(size = 20, hjust = 1, colour = "black"),
             aspect.ratio = 2/1, 
             panel.background = element_rect(fill = "snow"),
             panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black")) 
     
     write.table(xgb_imp_list, 
                 file = paste("../../output/importance_plot/xgb_imp_", j, ".tsv", sep = ""), 
                 sep = "\t", 
                 quote = FALSE,
                 row.names = FALSE)
     
     png(file = paste("../../output/importance_plot/xgb_imp_", j, ".png", sep = ""),
         width = 1000, height = 1100)
     
     print(xgb_imp_plot)
     
     dev.off()
     
     png(file = paste("../../output/importance_plot/xgb_imp_top_50_", j, ".png", sep = ""),
         width = 1000, height = 1100)
     
     print(xgb_imp_top_plot)
     
     dev.off()
     
     # Generate the CONFUSION FIGURE: 
     # Get the table in the right order:
     conf_xbg <- conf_xbg[, c(2,1,3)]
     # Get colnames in the right 
     colnames(conf_xbg) <- c("InitClass", "PredClass", "Overall")
     # Get the percentages: 
     conf_xbg$Percentage <- 0 # Create a numeric column
     # Get the TP and FN:
     conf_xbg[c(1:2), ]$Percentage <- round(conf_xbg[c(1:2), ]$Overall/
                                              sum(conf_xbg[c(1:2), ]$Overall)*100, digits = 0) 
     
     # Get the TN and FP:
     conf_xbg[c(3:4), ]$Percentage <- round(conf_xbg[c(3:4), ]$Overall/
                                              sum(conf_xbg[c(3:4), ]$Overall)*100, digits = 0) 
     
     # Generate the plot:
     xgb_confusio_plot <- 
       ggplot(data = conf_xbg, mapping = aes(x = InitClass, y = Percentage, fill = PredClass)) +
       geom_bar(stat = "identity", colour = "black") +
       geom_text(data = conf_xbg, aes(x = InitClass, y = Percentage, label = paste(Percentage,"%",sep="")), 
                 size = 12, position = "stack", colour = c("white", "chartreuse3", "blue3", "black"), 
                 alpha = c(1, 0, 0, 1), vjust = 1.2) +
       scale_fill_manual(values = alpha(c("blue3", "chartreuse3"), 0.7), labels = c(j, "Rest")) +
       scale_y_continuous(name = "Class Predictions (%)", breaks = seq(from = 0, to = 100, by = 20),
                          expand = c(0, 1)) +
       scale_x_discrete(labels = c(j, "Rest")) + 
       labs(x = "Classes") +
       coord_cartesian(ylim = c(1, 100)) +
       theme(plot.title = element_text(size = 40, hjust = 0.5),
             axis.title = element_text(size = 40),
             axis.text.x  = element_text(size = 38, hjust = 0.5, colour = "black"),
             axis.text.y  = element_text(size = 40, hjust = 0.5, colour = "black"),
             legend.title = element_blank(),
             legend.key.size = unit(1.3, "cm"), 
             legend.text  = element_text(size = 40),
             axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
             axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
             panel.background = element_rect(fill = "white"),
             axis.ticks.length = unit(0.4, "cm"))
     
     write.table(conf_xbg, 
                 file = paste("../../output/confusion_plot/xgb_conf_", j, ".tsv", sep = ""), 
                 sep = "\t", 
                 quote = FALSE,
                 row.names = FALSE)
     
     png(file = paste("../../output/3_machine_learning/confusion_plot/xgb_conf_", j, ".png", sep = ""),
         width = 900, height = 1100)
     
     print(xgb_confusio_plot)
     
     dev.off()
  
}
  
# -----------------------> FINAL RESULTS:
  
# Here, we are going to generate the final table with the results of all the
# constructed models and create an accuracy figure for all the models.
  
# I will create another script to identify the common relevant genes for all the 
# models. 
  
# Create a result table from the final table:
result_table <- final_table
  
# Round percentages to exact numbers:
result_table[, c(4, 7:10)] <- round(result_table[, c(4, 7:10)], digits = 0)
  
# Add 2 decimals to sensitivity and specificity:
result_table[, c(5,6)] <- round(result_table[, c(5,6)], digits = 2)
  
# Write the table with the results for each dataset:
  
write.table(result_table, 
            file = paste("../../output/3_machine_learning/ml_results.tsv", sep = ""),              
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE)
  
# -----> CREATE Accuracy figure:
  
accuracy_plot <- 
  ggplot(data = result_table, aes(x = cluster, y = percentage_accuracy, fill = machine_learning)) +
    geom_bar(stat = "identity",  position = position_dodge2(preserve = 'single'), colour = "black") +
    scale_fill_manual(values=c("firebrick2","dodgerblue2",'goldenrod1')) + 
    geom_text(position = position_dodge(w = 0.9), vjust = -0.7,
              aes(label = paste(percentage_accuracy, "%", sep = "")), size = 8.5) +
    scale_y_continuous(name = "% Accuracy Score", breaks = seq(from = 0, to = 100, by = 20), 
                       expand = c(0, 1)) + # Expand to delete the space between zero and the x-axis. 
    coord_cartesian(ylim = c(1, 100)) +
    theme(plot.title = element_text(size = 40, hjust = 0.5),
          axis.title = element_text(size = 40),
          axis.title.x = element_blank(),
          axis.text.x  =  element_text(size = 40, hjust = 1, angle = 45, colour = "black"), # Put color to the labels
          axis.text.y  = element_text(size = 40, hjust = 1, colour = "black"), # Put color to the labels
          axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
          axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
          panel.background = element_rect(fill = "white"),
          legend.title = element_blank(),
          legend.position = "top",
          legend.key.size = unit(1.3, "cm"), 
          legend.text  = element_text(size = 35),
          legend.spacing.x = unit(1, "cm"),
          axis.ticks.length = unit(0.4, "cm"))
  
png(file = paste("../../output/accuracy_plot.png", sep = ""),
    width = 1200, height = 1000)
  
print(accuracy_plot)
  
dev.off()


  
