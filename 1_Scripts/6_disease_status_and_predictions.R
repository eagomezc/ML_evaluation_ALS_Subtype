#-------------- ALS Subtype Evaluation: Disease status  ----------------------#

# Since the subtypes in brain are partially found in blood, we want to explore if
# close to death samples are more likely to be classify in brain subtypes, since 
# this subtypes were found in post-mortem samples. 

# In that sense, analyzing early and late stage of the disease and seeing how well
# it works, can give us a sense of if it will work closer to death. Late stage
# should be better than early.


#-----------------------------> LIBRARY LOAD:
if (!require('betareg')) install.packages('betareg'); library('betareg')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('ggeffects')) install.packages('ggeffects'); library('ggeffects')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('scales')) install.packages('scales'); library('scales')


if (!require('tidyr')) install.packages('tidyr'); library('tidyr')




set.seed(42) # To get same results even with the random part.

# Set directory to be in the source file location (Rstudio):
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# -------------------------> DATA LOAD: 

# Grima meta data: 
grima_meta <- read.delim("GSE234297_meta_GRIMA.txt")
grima_meta <- grima_meta[, c(1,3,4,9:12)]
grima_meta$Group <- gsub("PeripheralBlood_", "", grima_meta$Group)
grima_meta <- grima_meta[grima_meta$Deceased != 0, ] # Excluded controls and alive patients
colnames(grima_meta)[1] <- "Sample"

# Clusters:
cluster <- c("Neu", "OxA", "SNs")
final_logit <- data.frame()
final_prob <- data.frame()

# ------------------------> CORRELATION ANALYSIS:
  
  # Open LDA classification:
  LDA_class <- read.delim("LDA_WG_cluster_AllGenes_in_Grima.tsv")
  
  # Add a classification or non classification column:
  LDA_class$class <- "classified"
  LDA_class[LDA_class$Prob < 0.85, ]$class <- "Non_classified"
    
    # Merge meta and classification tables: 
    long_table <- merge(LDA_class, grima_meta, by = "Sample")
    colnames(long_table)[c(9,10)] <- c("Collection_paper", "Collection_medium")
    
    # Save long table:
    write.table(long_table, 
                file = "../../output/Longitudinal_table.tsv", 
                sep = "\t", 
                quote = FALSE,
                row.names = FALSE)
    
    # Test normality for the two variables:
    prob_s <- shapiro.test(long_table$Prob)
    cp_s <- shapiro.test(long_table$Collection_point)
    
    if (prob_s$p.value < 0.05 | cp_s$p.value < 0.05) {test <- "spearman"} else {test <- "pearson"}
    
    # Run the correlation analysis:
    corr_plot <- 
      ggscatter(long_table, x = "Collection_point", y = "Prob", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = test,
                cor.coef.size = 12, cor.coef.coord = c(0.25, max(long_table$Prob) + 0.1), 
                add.params = list(size = 1.5,  color = "black", fill = "skyblue"),    
                xlab = "Collection point", ylab = "Probability of assigment") +
      geom_point(data = long_table[long_table$Cluster == "Neu", ], color = "#619CFF", 
                 size = 5, position = position_dodge(0.75), show.legend = FALSE) +
      geom_point(data = long_table[long_table$Cluster == "OxA", ], color = "#00BA38", 
                 size = 5, position = position_dodge(0.75), show.legend = FALSE) +
      geom_point(data = long_table[long_table$Cluster == "SNs", ], color = "#F8766D", 
                 size = 5, position = position_dodge(0.75), show.legend = FALSE) +
      ylim(NA, max(long_table$Prob) + 0.15) +
      theme(axis.title = element_text(size = 40),
            plot.title = element_text(size = 50, hjust = 0.5, colour = "black"),
            axis.title.x = element_text(size = 40, hjust = 0.5, colour = "black"),
            axis.text.x  = element_text(size = 35, colour = "black"), # Put color to the labels
            axis.text.y  = element_text(size = 35, hjust = 1, colour = "black"), # Put color to the labels
            axis.line = element_line(colour = 'black', linewidth = 2), # Color and thickness of axis
            axis.ticks = element_line(colour = "black", linewidth = 2), # Color and thickness of every axis sep. 
            panel.background = element_rect(fill = "white"),
            legend.position = "none",
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.key.size = unit(1.5, "cm"), 
            legend.text  = element_text(size = 50),
            axis.ticks.length = unit(0.5, "cm")) 
    
    # Save the plot:
    png(file = paste("../../output/Correlation_All_clusters_", test, ".png", sep = ""),
        width = 1000, height = 1000)
    
    print(corr_plot)
    
    dev.off()
    
    # Correlation analysis per cluster:
    
    for (k in cluster) {
      
      # Define points colors:
      if (k == "Neu") {
        point_color <- "royalblue2"} else if (k == "OxA") {
          point_color <- "green4"} else if (k == "SNs") {
            point_color <- "firebrick3"}
      
      # Run the correlation analysis:
      corr_cluster_plot <- 
        ggscatter(long_table[long_table$Cluster == k, ], x = "Collection_point", y = "Prob", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = test, 
                  cor.coef.size = 20, cor.coef.coord = c(0.1, max(long_table$Prob) + 0.1), 
                  add.params = list(size = 2,  color = "black", fill = "skyblue"),    
                  xlab = "Collection point", ylab = "Probability of assigment") +
        geom_point(data = long_table[long_table$Cluster == k, ], color = point_color, 
                   size = 8, position = position_dodge(0.75), show.legend = FALSE) +
        ylim(NA, max(long_table$Prob) + 0.15) +
        theme(axis.title = element_text(size = 60, hjust = 0.5, colour = "black"),
              axis.text.x  = element_text(size = 55, colour = "black"), # Put color to the labels
              axis.text.y  = element_text(size = 55, hjust = 1, colour = "black"), # Put color to the labels
              axis.line = element_line(colour = 'black', linewidth = 3), # Color and thickness of axis
              axis.ticks = element_line(colour = "black", linewidth = 3), # Color and thickness of every axis sep. 
              panel.background = element_rect(fill = "white"),
              legend.position = "none",
              legend.title = element_blank(),
              legend.key = element_blank(),
              legend.key.size = unit(1.5, "cm"), 
              legend.text  = element_text(size = 50),
              axis.ticks.length = unit(0.5, "cm")) 
      
      # Save the plot:
      png(file = paste("../../output/Correlation_", k, "_", test, ".png", sep = ""), width = 1000, height = 1000)
      
      print(corr_cluster_plot)
      
      dev.off()
      
    }
    
    # ------------------------> LOGISTIC REGRESSION:

    # Logistic Regression Models to test disease status effect in classification in blood:
    # The brain-trained ML model performs better in late-stage patients than early-state one.
    
    # Creates binary classification:
    long_table$class_bi <- 1
    long_table[long_table$class == "Non_classified", ]$class_bi <- 0
    
    # Convert categorical variables in factors: 
    long_table[, c(2,4,5,9:11)] <- lapply(long_table[, c(2,4,5,9:11)], as.factor)
    
    # Create logistic model (Disease status):
    # We want to build a model that predicts the classified or unclassified samples, 
    # either non_classified (0) or classified (1), based on the disease status. 
    # The response variable is the binary variable class_bi and the predictor variable is Collection point.
    # With Sex and Age at collection as covariates. We build a logit model by applying 
    # the glm() function. For the logistic regression model we specify family = 'binomial'.
    
    # Generate a loop with predictors:
    
    predictors <- c("Collection_point", "Collection_paper")
    
    for (l in predictors) {
    
      # ------> CLASSIFICATION vs NON ClASSIFICATION:
      
      # When using binary response variable and continuous explanatory variables:
      # Intercept represent log odds of classification (1), when the explanatory variables are zero.
      # The rest of Coefficients represents the effect of that specific variable over class (1)
      # considering the others as covariates. 
      
      # When using binary response and categorical explanatory variable:
      # Intercept represent log odds of class (1), when the explanatory variable is collection[early].
      # The rest of coef represent the effect of the non reference value collection[late] over class (1)
      # considering other covariates.
      
      # ------> NO Cluster:
      
      # Create formula for the glm: as class_bi ~ Collection_point + Sex + Age_at_collection:
      eq <- as.formula(paste("class_bi ~", l, "+ Sex + Age_at_collection"))
      
      # Model without considering cluster:
      model_no_cluster <- glm(eq, data = long_table, family = "binomial")
      
      # Plot the models:
      if (l == "Collection_point") {
        
        # Plot predicted probability vs Continuous variable:
        # Calculate prediction 
        pred <- ggpredict(model_no_cluster, terms = l)
        
        # Plot: 
        plot_model <-
          ggplot(pred, aes(x, predicted)) +
          geom_line(size = 1) +
          geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
          labs(y = "Predicted probability", x = "Collection Point") +
            theme(plot.title = element_text(size = 15, colour = "black", hjust = 0.5),
                  axis.title.y = element_text(size = 30, colour = "black"),
                  axis.title.x = element_text(size = 30, colour = "black"),
                  axis.text.x  =  element_text(size = 30, colour = "black"), # Put color to the labels
                  axis.text.y  = element_text(size = 30, colour = "black"), # Put color to the labels
                  axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                  axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
                  panel.background = element_rect(fill = "white"),
                  legend.title = element_text(size = 30),
                  legend.text  = element_text(size = 25),
                  legend.key = element_rect(fill = NA),
                  legend.key.size = unit(1.3, "cm"), 
                  axis.ticks.length = unit(0.4, "cm"))
      
      # Save the model:
      png(file = paste("../../output/Binary_NoCl_", l, ".png", sep = ""), width = 1000, height = 1000)
      
      print(plot_model)
      
      dev.off() } else {
        
        pred <- ggpredict(model_no_cluster, terms = l)
        
        if (nrow(pred) == 3) {pred$x <- factor(pred$x, levels = c("early", "medium", "late"))}
        
        plot_model <-
          ggplot(pred, aes(x = x, y = predicted)) +
            geom_point(size = 4) +
            geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, size = 1) +
            labs(y = "Predicted probability", x = "Collection Point") +
            theme(plot.title = element_text(size = 15, colour = "black", hjust = 0.5),
                  axis.title.y = element_text(size = 30, colour = "black"),
                  axis.title.x = element_text(size = 30, colour = "black"),
                  axis.text.x  =  element_text(size = 30, colour = "black"), # Put color to the labels
                  axis.text.y  = element_text(size = 30, colour = "black"), # Put color to the labels
                  axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                  axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
                  panel.background = element_rect(fill = "white"),
                  legend.title = element_text(size = 30),
                  legend.text  = element_text(size = 25),
                  legend.key = element_rect(fill = NA),
                  legend.key.size = unit(1.3, "cm"), 
                  axis.ticks.length = unit(0.4, "cm"))
          
        png(file = paste("../../output/Binary_NoCl_", l, ".png", sep = ""), width = 1000, height = 1000)
          
        print(plot_model)
          
        dev.off()
        
      }
      
      # Get coefficient tables:
      coef_table <- as.data.frame(summary(model_no_cluster)$coefficients)
      
      # Save all the table of coefficient for future checking:
      write.table(coef_table, 
                  file = paste("../../output/Binary_No_cluster_", l, ".tsv", sep = ""), 
                  sep = "\t", 
                  quote = FALSE,
                  col.names = NA,
                  row.names = TRUE)
      
      # -------> CLUSTER: 
      
      # When class binary and continuous explanatory adding Cluster combined effect:
      # Intercept represents Collection Point at 0 with Neu considering females. 
      # Other variables are there is with covariates.
      # And then Collection_point:ClusterOxA the interaction effect: How the effect of Status 
      # changes depending on the Subtype.
      
      # When binary response and binary explanatory, adding Cluster combined effect:
      # The interaction effect is for the Collection Point[late] since early is the Intercept. 
      
      # Create formula for the glm: class_bi ~ Collection_point*Cluster + Sex + Age_at_collection:
      eq_cl <- as.formula(paste("class_bi ~", l, "*Cluster + Sex + Age_at_collection"))
      
      # Model considering cluster:
      model_cluster <- glm(eq_cl, data = long_table, family = "binomial")
      
      # Plot the models:
      if (l == "Collection_point") {
        
        # Plot predicted probability vs Continuous variable:
        # Calculate prediction 
        pred <- ggpredict(model_cluster, terms = c(l, "Cluster"))
        
        # Plot: 
        plot_model <-
          ggplot(pred, aes(x, predicted, color = group)) +
            geom_line(size = 1) +
            geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
            labs(y = "Predicted probability", x = "Collection Point") +
            scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
            facet_wrap(~group, labeller=label_parsed, ncol = 3, nrow = 1) + 
            theme(plot.title = element_text(size = 15, colour = "black", hjust = 0.5),
                  axis.title.y = element_text(size = 30, colour = "black"),
                  axis.title.x = element_text(size = 30, colour = "black"),
                  axis.text.x  =  element_text(size = 30, colour = "black"), # Put color to the labels
                  axis.text.y  = element_text(size = 30, colour = "black"), # Put color to the labels
                  axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                  axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
                  panel.background = element_rect(fill = "white"),
                  legend.title = element_text(size = 30),
                  legend.text  = element_text(size = 25),
                  legend.key = element_rect(fill = NA),
                  legend.key.size = unit(1.3, "cm"), 
                  strip.text.x = element_text(size = 20),
                  strip.background =element_rect(fill="snow2"),
                  legend.position = "none",
                  axis.ticks.length = unit(0.4, "cm"))
        
        # Save the model:
        png(file = paste("../../output/Binary_Cl_", l, ".png", sep = ""), width = 1000, height = 1000)
        
        print(plot_model)
        
        dev.off() } else {
          
          pred <- ggpredict(model_cluster, terms = c(l, "Cluster"))
          
          if (nrow(pred) == 9) {pred$x <- factor(pred$x, levels = c("early", "medium", "late"))}
          
          plot_model <- 
            ggplot(pred, aes(x = x, y = predicted, color = group)) +
            geom_point(position = position_dodge(0.3), size = 4) +
            geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                          position = position_dodge(0.3),
                          width = 0.2, size = 1) +
              scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
              labs(y = "Predicted probability", x = "Collection Point") +
              guides(color = guide_legend(title = "Cluster")) +
              theme(plot.title = element_text(size = 15, colour = "black", hjust = 0.5),
                    axis.title.y = element_text(size = 30, colour = "black"),
                    axis.title.x = element_text(size = 30, colour = "black"),
                    axis.text.x  =  element_text(size = 30, colour = "black"), # Put color to the labels
                    axis.text.y  = element_text(size = 30, colour = "black"), # Put color to the labels
                    axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                    axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
                    panel.background = element_rect(fill = "white"),
                    legend.title = element_text(size = 30),
                    legend.text  = element_text(size = 25),
                    legend.key = element_rect(fill = NA),
                    legend.key.size = unit(1.3, "cm"), 
                    axis.ticks.length = unit(0.4, "cm"))
          
          png(file = paste("../../output/Binary_Cl_", l, ".png", sep = ""), width = 1000, height = 1000)
          
          print(plot_model)
          
          dev.off()
          
        }
      
      # Get coefficient tables:
      coef_table <- as.data.frame(summary(model_cluster)$coefficients)
      
      # Save all the table of coefficient for future checking:
      write.table(coef_table, 
                  file = paste("../../output/Binary_Cluster_", l, ".tsv", sep = ""), 
                  sep = "\t", 
                  quote = FALSE,
                  col.names = NA,
                  row.names = TRUE)
      
      # ------> PROBABILITY OF CLASSIFICATION:
    
      
      # ------> NO Cluster:
      
      # When both are probability, Intercept represents Collection_point at zero.
      # Coeff refers to that specific variable effect over the probability (pos value 
      # indicates pos correlation and reverse). 
      
      # When continuous response and categorical explanatory: Intercept represents 
      # "The predicted probability for the early status group".
      # Status[Late]: This means the "large" group has a probability 2 times higher than the "early" group.
      
      # Create formula for the glm: as Prob ~ Collection_point + Sex + Age_at_collection:
      eq <- as.formula(paste("Prob ~", l, "+ Sex + Age_at_collection"))
      
      # Model without considering cluster:
      model_no_cluster <- glm(eq, data = long_table, Gamma(link = "log"))
      
      # Plot the models:
      if (l == "Collection_point") {
        
        # Plot predicted probability vs Continuous variable:
        # Calculate prediction 
        pred <- ggpredict(model_no_cluster, terms = l)
        
        # Plot: 
       plot_model <-
        ggplot(pred, aes(x, predicted)) +
          geom_line(size = 1) +
          geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
          labs(y = "Predicted probability", x = "Collection Point") +
          theme(plot.title = element_text(size = 15, colour = "black", hjust = 0.5),
                axis.title.y = element_text(size = 30, colour = "black"),
                axis.title.x = element_text(size = 30, colour = "black"),
                axis.text.x  =  element_text(size = 30, colour = "black"), # Put color to the labels
                axis.text.y  = element_text(size = 30, colour = "black"), # Put color to the labels
                axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
                panel.background = element_rect(fill = "white"),
                legend.title = element_text(size = 30),
                legend.text  = element_text(size = 25),
                legend.key = element_rect(fill = NA),
                legend.key.size = unit(1.3, "cm"), 
                axis.ticks.length = unit(0.4, "cm"))
        
        # Save the model:
        png(file = paste("../../output/Prob_NoCl_", l, ".png", sep = ""), width = 1000, height = 1000)
        
        print(plot_model)
        
        dev.off() } else {
          
          pred <- ggpredict(model_no_cluster, terms = l)
          
          if (nrow(pred) == 3) {pred$x <- factor(pred$x, levels = c("early", "medium", "late"))}
          
          plot_model <-
            ggplot(pred, aes(x = x, y = predicted)) +
            geom_point(size = 4) +
            geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, size = 1) +
            labs(y = "Predicted probability", x = "Collection Point") +
            theme(plot.title = element_text(size = 15, colour = "black", hjust = 0.5),
                  axis.title.y = element_text(size = 30, colour = "black"),
                  axis.title.x = element_text(size = 30, colour = "black"),
                  axis.text.x  =  element_text(size = 30, colour = "black"), # Put color to the labels
                  axis.text.y  = element_text(size = 30, colour = "black"), # Put color to the labels
                  axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                  axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
                  panel.background = element_rect(fill = "white"),
                  legend.title = element_text(size = 30),
                  legend.text  = element_text(size = 25),
                  legend.key = element_rect(fill = NA),
                  legend.key.size = unit(1.3, "cm"), 
                  axis.ticks.length = unit(0.4, "cm"))
          
          png(file = paste("../../output/Prob_NoCl_", l, ".png", sep = ""), width = 1000, height = 1000)
          
          print(plot_model)
          
          dev.off()
          
        }
      
      # Get coefficient tables:
      coef_table <- as.data.frame(summary(model_no_cluster)$coefficients)
      coef_table_save <- coef_table
      coef_table_save$perc_incr_decr <- (exp(coef_table_save$Estimate) - 1)*100
      
      # Save all the table of coefficient for future checking:
      write.table(coef_table_save, 
                  file = paste("../../output/Prob_No_cluster_", l, ".tsv", sep = ""), 
                  sep = "\t", 
                  quote = FALSE,
                  col.names = NA,
                  row.names = TRUE)
      
      # -------> CLUSTER: 
      
      # In case of using Beta regression, make sure there is no absolute zeros and 1:
      long_table$b_prob <- (long_table$Prob * (nrow(long_table) - 1) + 0.5)/nrow(long_table)
      
      # The Intercept represents the expected log-mean when all categorical predictors 
      # are at their reference level.
      # It is the mean probability for a subject that is both Early and Neu.
      
      # The interaction terms (e.g., Status[med]:Subtype[beta]) are the most complex part. 
      # They answer: "Is the effect of being Medium the same for Alpha as it is for Beta?"
      # If the interaction coefficient is significant (p < 0.05), it means the categories 
      # don't just "add up" independently—they have a unique synergy.
      # To calculate the probability for a non-reference cell (e.g., Medium + Beta):
      # You must add the Intercept + Main Effect (Med) + Main Effect (Beta) + Interaction (Med:Beta).
      
      # For beta Regression:
      eq_cl <- as.formula(paste("b_prob ~", l, "*Cluster + Sex + Age_at_collection"))
      
      # Beta regression:
      model_cluster <- betareg(eq_cl, data = long_table)
      
      # Plot the models:
      if (l == "Collection_point") {
        
        # Plot predicted probability vs Continuous variable:
        # Calculate prediction 
        pred <- ggpredict(model_cluster, terms = c(l, "Cluster"))
        
        # Save prediction table:
        write.table(pred, 
                    file = paste("../../output/beta_Prob_Cl_", l, "_pred.tsv", sep = ""), 
                    sep = "\t", 
                    quote = FALSE,
                    row.names = FALSE)
        
        # Plot: 
        plot_model <-
         ggplot(pred, aes(x, predicted, color = group)) +
          geom_line(size = 2) +
          geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = "white", 
                      fill = "skyblue") +
          labs(y = "Predicted probability", x = "Collection Point") +
          scale_color_manual(values = c("royalblue2", "green4", "firebrick3")) +
          scale_y_continuous(labels = label_number(accuracy = 0.01)) +
          scale_x_continuous(labels = label_number(accuracy = 0.1)) +
          facet_wrap(~group, labeller=label_parsed, ncol = 3, nrow = 1, scales = "free") + 
          theme(axis.title.y = element_text(size = 50, colour = "black"),
                axis.title.x = element_text(size = 50, colour = "black"),
                axis.text.x  =  element_text(size = 40, colour = "black"), # Put color to the labels
                axis.text.y  = element_text(size = 40, colour = "black"), # Put color to the labels
                axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
                panel.background = element_rect(fill = "white"),
                legend.title = element_text(size = 30),
                legend.text  = element_text(size = 25),
                legend.key = element_rect(fill = NA),
                legend.key.size = unit(1.3, "cm"), 
                strip.text.x = element_text(size = 50),
                strip.background =element_rect(fill="snow2"),
                legend.position = "none",
                axis.ticks.length = unit(0.4, "cm"))
        
        # Save the model:
        png(file = paste("../../output/beta_Prob_Cl_", l, ".png", sep = ""), width = 1600, height = 800)
        
        print(plot_model)
        
        dev.off() } else {
          
          pred <- ggpredict(model_cluster, terms = c(l, "Cluster"))
          
          if (nrow(pred) == 9) {pred$x <- factor(pred$x, levels = c("early", "medium", "late"))}
          
          # Save prediction table:
          write.table(pred, 
                      file = paste("../../output/beta_Prob_Cl_", l, "_pred.tsv", sep = ""), 
                      sep = "\t", 
                      quote = FALSE,
                      row.names = FALSE)
          
          plot_model <- 
            ggplot(pred, aes(x = x, y = predicted, color = group)) +
            geom_point(position = position_dodge(0.3), size = 6) +
            geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                          position = position_dodge(0.3),
                          width = 0.2, size = 2) +
            scale_color_manual(values = c("royalblue2", "green4", "firebrick3")) +
            labs(y = "Predicted probability", x = "Collection Point") +
            guides(color = guide_legend(title = "Cluster")) +
            theme(axis.title.y = element_text(size = 50, colour = "black"),
                  axis.title.x = element_text(size = 50, colour = "black"),
                  axis.text.x  =  element_text(size = 40, colour = "black"), # Put color to the labels
                  axis.text.y  = element_text(size = 40, colour = "black"), # Put color to the labels
                  axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                  axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
                  panel.background = element_rect(fill = "white"),
                  legend.title = element_text(size = 45),
                  legend.text  = element_text(size = 35),
                  legend.key = element_rect(fill = NA),
                  legend.key.size = unit(1.3, "cm"), 
                  axis.ticks.length = unit(0.4, "cm"))
          
          png(file = paste("../../output/beta_Prob_Cl_", l, ".png", sep = ""), width = 1000, height = 800)
          
          print(plot_model)
          
          dev.off()
          
        }
      
      # Get coefficient tables:
      coef_table <- as.data.frame(summary(model_cluster)$coefficients)
      coef_table_save <- coef_table
      coef_table_save$odds_ratio <- exp(coef_table_save$mean.Estimate)
      
      # Save all the table of coefficient for future checking:
      write.table(coef_table_save, 
                  file = paste("../../output/beta_Prob_Cluster_", l, ".tsv", sep = ""), 
                  sep = "\t", 
                  quote = FALSE,
                  col.names = NA,
                  row.names = TRUE)
    
    }
  


