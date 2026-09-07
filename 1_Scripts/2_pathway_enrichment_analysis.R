#----------- ALS Subtype Evaluation: Pathway Enrichment Analysis --------------#

# After performing unsupervised clustering and classified samples in different
# ALS subtypes and extracted the relevant genes for clustering, now we explored
# this genes to see which enriched pathways are associated with each subtype. 

# -------------------------> LIBRARY LOAD: 

if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if (!require('clusterProfiler')) BiocManager::install('clusterProfiler', update = FALSE); library('clusterProfiler')
if (!require('org.Hs.eg.db')) BiocManager::install("org.Hs.eg.db", character.only = TRUE); library("org.Hs.eg.db")
if (!require('ReactomePA')) BiocManager::install("ReactomePA", character.only = TRUE); library("ReactomePA")
if (!require('stringr')) install.packages('stringr'); library('stringr')

# Set directory to be in the source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
op <- par(no.readonly = TRUE) # Image display

# -------------------------> DATA LOAD: 

# Open table with ensembl and symbol names:
ens_to_sym <- read.delim("sex_chr/Ens_to_symbol.txt", header = FALSE, row.names = 1)
colnames(ens_to_sym) <- "Symbol"

# Open table with backgroun genes:
bck_genes <- rownames(read.delim("Ens_all_Genes.tsv", row.names = 1))

# -----------------------> ENRICHMENT ANALYSIS: 
  
# The list of informative genes for each cluster was then used 
# to characterise their molecular phenotypes by performing gene enrichment 
# analysis using the clusterProfiler and ReactomePA R package. 
# Genes from the whole expression matrix were used as a custom gene background. 
# Evaluated databases:
# Gene Ontology (Biological Process (GO:BP), Molecular Function (GO:MF) and Cellular Component (GO:CC)).
# Reactome.
  
# -----> CLUSTER genes:
  
rel_genes <- read.delim("../../output/2_unsupervised_clustering/All_Genes/Cluster_features_All_Genes.tsv")
  
# Get the genes from the RELEVANT GENE TABLE: 
rel_cl_genes <- rel_genes[, c(1,3)]
colnames(rel_cl_genes)[1] <- "Symbol"
  
# Get the genes in the useful format: 
rel_inf_genes <- split(rel_cl_genes$Symbol, rel_cl_genes$Group)
  
# -----> ENRICHMENT analysis:
  
# Enrichment analysis GO terms:
go_enrich <- lapply(names(rel_inf_genes), function(cluster) {
  simplify(
  enrichGO(
    gene = rel_inf_genes[[cluster]],
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "ALL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe = bck_genes,
    minGSSize = 10,
    maxGSSize = 500))})
  
# Enrichment analysis Reactome:
# Get background genes in the right nomenclature:
bck_entriz <- bitr(bck_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  
# Run Reactome enrichment anaylysis:
react_enrich <- lapply(names(rel_inf_genes), function(cluster) {
  # Get the right nomenclature:
  entrez_id <- bitr(rel_inf_genes[[cluster]], fromType = "ENSEMBL", toType = "ENTREZID",
                    OrgDb = "org.Hs.eg.db")
  
  enrichPathway(
    gene = entrez_id$ENTREZID,
    universe = bck_entriz$ENTREZID,
    organism = "human",
    minGSSize = 10,
    maxGSSize = 500,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE)})
  
# Get plots and tables: 
  for (n in 1:length(go_enrich)) {
    
    # Get the table with enrich results: 
    go_results <- go_enrich[[n]]@result
    react_results <- react_enrich[[n]]@result
    
    if (!is.null(go_results)) { # only generate plot if there is info in the table
      
      # Go paper:
      go_paper <- go_results[, c(1:4, 7, 9, 10, 12, 13)]
      go_paper$Description <- str_to_sentence(go_paper$Description)
      go_paper$final_term <- paste(go_paper$ONTOLOGY, go_paper$Description)
      go_paper$generatio <- sapply(go_paper$GeneRatio, function(x) eval(parse(text = x)))
      go_paper$logpval <- -log10(go_paper$p.adjust)
      go_paper <- go_paper[go_paper$Count >= 3, ]
      go_paper <- go_paper[, c(1,2,10,11,5:7,12,9,8)]
      colnames(go_paper) <- c("Database", "ID", "Final Term", "Gene Ratio", "Fold Enrichment",
                              "p value", "adjust p value", "-Log(Adj P)", "Count", "Gene IDs")
      
      # Some table work:
      go_results <- go_results[ , c(1, 3, 4, 10, 13)]
      go_results$final_term <- paste(go_results$ONTOLOGY, go_results$Description)
      go_results$generatio <- sapply(go_results$GeneRatio, function(x) eval(parse(text = x)))
      go_results$logpval <- -log10(go_results$p.adjust)
      go_results <- go_results[, c(6:8,5)]
      colnames(go_results)[c(2:4)] <- c("GeneRatio", "-Log(Adj P)", "Count")
      go_results <- go_results[order(go_results$`-Log(Adj P)`, decreasing = TRUE), ]
      #go_results <- go_results[go_results$Count >= 3, ]
      
      write.table(go_results, 
                  file = paste("../../output/GO_Cluster_inf_", n, ".tsv", sep = ""), 
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE) 
      
      # If the table has a lot of results, the figure cant be read, so I reduce to 100: 
      if (nrow(go_results > 100)) {go_results <- head(go_results, 100)}
      
      # Get the enrich plot: 
      go_plot <-
        ggplot(go_results, aes(x = `-Log(Adj P)`, y = reorder(final_term, `-Log(Adj P)`))) +
        geom_point(aes(size = Count, color = GeneRatio)) +
        scale_color_gradient(low = "red", high = "blue", name = "Gene Ratio") +
        scale_size_continuous(range = c(1, 10), breaks = seq(0, max(go_results$Count), by = 20), name = "Count") +
        labs(y = NULL, x = "-log(Adj P)") +
        theme(axis.title = element_text(size = 25, hjust = 0.5),
              plot.title = element_text(size = 28, hjust = 0.5),
              axis.text.x  =  element_text(size = 22, colour = "black"),
              axis.text.y  =  element_text(colour = "black"),
              legend.title = element_text(size = 22),
              legend.text = element_text(size = 20),
              legend.position = "right",
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black")) 
      
      png(file = paste("../../output/GO_Cluster_inf_", n, ".png", sep = ""), 
          width = 1500, height = 800)
      
      print(go_plot)
      
      dev.off()
      
    }
      
      if (!is.null(react_results)) { # only generate plot if there is info in the table
        
        # Go paper:
        react_paper <- react_results[, c(1:3, 6, 8, 9, 11, 12)]
        react_paper$Database <- "React"
        react_paper$final_term <- paste(react_paper$Database, react_paper$Description)
        react_paper$generatio <- sapply(react_paper$GeneRatio, function(x) eval(parse(text = x)))
        react_paper$logpval <- -log10(react_paper$p.adjust)
        react_paper <- react_paper[react_paper$Count >= 3, ]
        react_paper <- react_paper[react_paper$p.adjust <= 0.05, ]
        react_paper <- react_paper[, c(9,1,10,11,4:6,12,8,7)]
        colnames(react_paper) <- c("Database", "ID", "Final Term", "Gene Ratio", "Fold Enrichment",
                                   "p value", "adjust p value", "-Log(Adj P)", "Count", "Gene IDs")
        
        # Some table work:
        react_results <- react_results[ , c(2, 3, 9, 12)]
        react_results <- react_results[react_results$p.adjust <= 0.08, ]
        react_results$final_term <- paste("React", react_results$Description)
        react_results$generatio <- sapply(react_results$GeneRatio, function(x) eval(parse(text = x)))
        react_results$logpval <- -log10(react_results$p.adjust)
        react_results <- react_results[, c(5:7,4)]
        colnames(react_results)[c(2:4)] <- c("GeneRatio", "-Log(Adj P)", "Count")
        react_results <- react_results[order(react_results$`-Log(Adj P)`, decreasing = TRUE), ]
        # react_results <- react_results[react_results$Count >= 3, ]
        
        write.table(react_results, 
                    file = paste("../../output/REACT_Cluster_inf_", n, ".tsv", sep = ""), 
                    sep = "\t",
                    quote = FALSE,
                    row.names = FALSE)
        
        # If the table has a lot of results, the figure cant be read, so I reduce to 100: 
        if (nrow(react_results > 100)) {react_results <- head(react_results, 100)}
        
        # Get the enrich plot:
        react_plot <-
          ggplot(react_results, aes(x = `-Log(Adj P)`, y = reorder(final_term, `-Log(Adj P)`))) +
          geom_point(aes(size = Count, color = GeneRatio)) +
          scale_color_gradient(low = "red", high = "blue", name = "Gene Ratio") +
          scale_size_continuous(range = c(1, 10), breaks = seq(0, max(react_results$Count), by = 20), name = "Count") +
          labs(y = NULL, x = "-log(Adj P)") +
          theme(axis.title = element_text(size = 25, hjust = 0.5),
                plot.title = element_text(size = 28, hjust = 0.5),
                axis.text.x  =  element_text(size = 22, colour = "black"),
                axis.text.y  =  element_text(colour = "black"),
                legend.title = element_text(size = 22),
                legend.text = element_text(size = 20),
                legend.position = "right",
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black")) 
        
        png(file = paste("../../output/REACT_Cluster_inf_", n, ".png", sep = ""), 
            width = 1500, height = 800)
        
        print(react_plot)
        
        dev.off()
    
      }
    
    # Combined table and plots:
    go_and_react <- rbind(go_paper, react_paper)
    
    # Save all table: 
    write.table(go_and_react, 
                file = paste("../../output/All_G&R_Cluster_inf_", n, ".tsv", sep = ""), 
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
    
    # Prepare table for plot:
    go_and_react <- go_and_react[order(go_and_react$`-Log(Adj P)`, decreasing = TRUE), ] # Order by pvalue
    go_and_react <- go_and_react[go_and_react$Database == "BP" | go_and_react$Database == "React", ] # Keep just BP and Reactome
    
    all_plot <-
      ggplot(go_and_react, aes(x = `-Log(Adj P)`, y = reorder(`Final Term`, `-Log(Adj P)`))) +
      geom_point(aes(size = Count, color = `Gene Ratio`)) +
      scale_color_gradient(low = "red", high = "blue", name = "Gene Ratio") +
      scale_size_continuous(range = c(1, 10), breaks = seq(0, max(go_and_react$Count), by = 10), name = "Count") +
      labs(y = NULL, x = "-log(Adj P)") +
      xlim(NA, max(go_and_react$`-Log(Adj P)`) + 0.2) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
      theme(axis.title = element_text(size = 55, hjust = 0.5),
            axis.text.x  =  element_text(size = 45, colour = "black"),
            axis.text.y  =  element_text(size = 15, colour = "black", lineheight = 0.6),
            legend.title = element_text(size = 45),
            legend.text = element_text(size = 30),
            legend.position = "right",
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black")) 
    
    png(file = paste("../../output/All_G&R_Cluster_inf_", n, ".png", sep = ""), 
        width = 1500, height = 1000)
    
    print(all_plot)
    
    dev.off()

  }


    
    
    
  
  
  
  
