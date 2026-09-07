#-------- ALS Subtype Evaluation: 3-way Differential Gene Expression ----------#

# Since we didn't know the real classification for the blood samples, there is no
# way to evaluate the LDA models. A indirect way to evaluate them is to run 
# differential gene expression analysis to identify up and down-regulated genes 
# for each subtype and then run and enrichment analysis. Enriched pathways are then
# compared with the identified pathways from the unsupervised clustering and see
# if we got a match. 

# -------------------------> LIBRARY LOAD: 

if (!require('ggplot2')) install.packages('gplots'); library('ggplot2')
if (!require('volcano3D')) install.packages('volcano3D'); library('volcano3D')
if (!require('htmlwidgets')) install.packages('htmlwidgets'); library('htmlwidgets')
if (!require('Hmisc')) install.packages('Hmisc'); library('Hmisc')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('plotly')) install.packages('plotly'); library('plotly')
if (!require('stringr')) install.packages('stringr'); library('stringr')
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if (!require('DESeq2')) BiocManager::install('DESeq2', update = FALSE); library('DESeq2')
if (!require('clusterProfiler')) BiocManager::install('clusterProfiler', update = FALSE); library('clusterProfiler')
if (!require('org.Hs.eg.db')) BiocManager::install("org.Hs.eg.db", character.only = TRUE); library("org.Hs.eg.db")
if (!require('ReactomePA')) BiocManager::install("ReactomePA", character.only = TRUE); library("ReactomePA")

# Set directory to be in the source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
op <- par(no.readonly = TRUE) # Image display

# -------------------------> DATA LOAD: 

# Open the tsv file the meta data information:

grima_meta <- read.delim("GSE234297_meta_GRIMA.txt")
grima_meta <- grima_meta[grepl("Case", grima_meta$Group), ] # Get only cases
grima_meta <- grima_meta[, c(1, 3:8)] # Keep just the relevant columns.
colnames(grima_meta)[1] <- "Sample"
grima_meta$Sample <- gsub("PeripheralBlood_", "", grima_meta$Sample)

# Create categorical variables:
grima_meta$AgeCat <- cut(grima_meta$Age_at_collection, breaks = 5, labels = c(1, 2, 3, 4, 5))
grima_meta$RINCat <- cut(grima_meta$RIN, breaks = 5, labels = c(1, 2, 3, 4, 5))

# Open the tsv file with the classification information:
grima_class <- read.delim("LDA_WG_cluster_AllGenes_in_Grima.tsv") # All these samples are cases.

# Open the tsv file with the RAW read counts: 
grima_counts <- read.delim("GSE234297_gene_raw_GRIMA.txt")
# Get only cases:
grima_counts <- grima_counts[ ,colnames(grima_counts) == "EntrezGeneID" | colnames(grima_counts)
                              %in% grima_meta$Sample]

# Open the txt file with genes in sex chromosomes:
sex_genes <- read.delim("sex_with_symbols.txt", sep = " ", header = FALSE)
sex_genes <- sex_genes$V2

# Open the table with the entriz to ensembl names:
ent_to_ens <- read.delim("GeneID_to_Ens_GRIMA.tsv")
ent_to_ens <- ent_to_ens[, c(7,8,2)]

# -----------------------> DATA MANIPULATION:
  
  # ------> MERGE meta data and class file:
  
  # Merge by sample name:
  grima_meta_cases <- merge(grima_meta, grima_class[, c(1:2)], by = "Sample")
  rownames(grima_meta_cases) <- grima_meta_cases$Sample
  grima_meta_cases[c(2,10)] <- lapply(grima_meta_cases[c(2,10)], factor)

  # ------> Change the name of the genes in the raw file:
  grima_counts_cases <- merge(grima_counts, ent_to_ens[, c(1,3)], by = "EntrezGeneID")
  grima_counts_cases <- grima_counts_cases[!duplicated(grima_counts_cases$symbol), ]
  rownames(grima_counts_cases) <- grima_counts_cases$symbol
  grima_counts_cases <- grima_counts_cases[, -c(1,98)] 
  grima_counts_cases <- round(grima_counts_cases, 0)
  
  # Remove sex cromosomes: 
  grima_counts_cases <- grima_counts_cases[rownames(grima_counts_cases) %nin% sex_genes, ] # Exclude sex genes
  
  # PCA analysis for the grima dataset:
  # The DESeq class will contain the reads, class info and "design = ~1" indicates no design
  # for analysis (meaning not based on groups, or adding covariants, etc.)
  
  dds <- DESeqDataSetFromMatrix(countData = grima_counts_cases, colData = grima_meta_cases, design = ~Sex + AgeCat + RINCat + Cluster)
  
  # Calculate the Size Factors:
  dds <- estimateSizeFactors(dds)

  # ------> FILTERING and normalization: 

  # Gene to be expressed at a reasonable level in a sample if it has 5 counts in that 
  # sample. Gene should be expressed in at least 10 samples (as indicated but Marriott script).  

  keep <- rowSums(counts(dds, normalized = TRUE) >= 5) >= 10
  filter_dds <- dds[keep, ]

  # Normalization was performed using VST.
  # Explanation of VST (Variance-stabilizing transformation):
  # Data transformation used to make the variance approximately constant across its range
  # of values. Fix a mean-dependent variance model that shrink or increase the variance
  # of the different features so the data become constant and somehow normalize without
  # loosing the variability between groups. 

  vsd <- vst(filter_dds, blind=TRUE) # This indicates that it does not consider groups for the calculations.
  
  # Get the backgroun genes for enrichment analysis:
  bck_genes <- rownames(assay(vsd))

  # -------> PCA analysis: 

  # PCA plot
  pca_data <- plotPCA(vsd, intgroup = c("Sex", "Age_at_collection", "Age_of_onset",
                                        "Disease_duration_months", "RIN", "Cluster"), 
                                      returnData = TRUE, ntop = 5000)
  percentVar <- round(100 * attr(pca_data, "percentVar")) # Calculate variance
  
  # Reorganize table:
  pca_data <- pca_data[, c(10,1,2,9,4:8)] 

  # Create empty list:
  plot_list <- list()

  for (j in 4:9) {
    
    feature <- colnames(pca_data)[j]
  
    pca <- 
      ggplot(pca_data, aes(PC1, PC2, color = !!sym(feature))) +
      geom_point(size = 3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      ggtitle(paste("PCA of VST-counts (", feature, ")", sep = "")) +
      {if (feature == "Cluster") scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D"))} +
      {if (feature == "Sex") scale_color_manual(values = c("red", "blue"))} +
      {if (feature == "Age_at_collection") scale_color_gradient(low = "lavenderblush4", high = "coral")} +
      {if (feature == "Disease_duration_months") scale_color_gradient(low = "grey8", high = "grey50")} +
      {if (feature == "Age_of_onset") scale_color_gradient(low = "plum2", high = "purple4")} +
      {if (feature == "RIN") scale_color_gradient(low = "lightblue", high = "dodgerblue4")} +
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
    plot_list[[j]] <- pca # From the previous empty list, save the plots.
  
  }
  
  # Since the for start with, we need to delete the empty 3 elements in the list:
  plot_list <- plot_list[lengths(plot_list) > 0]
  
  png(file = "../../../output/PCA_VST_grima.png", width = 1500, height = 800)
  
  # Save all the plots in one grid: 
  gridExtra::grid.arrange(grobs = plot_list, ncol = 3)
  
  dev.off()
  
  # -----------------------> VOLCANO 3D PLOT:
  
  # Three-way differential expression. RLE counts were Z-score normalized and mean Z scores
  # calculated for the three groups. This three dimensional data were reduced to two dimensional
  # polar coordinate system analogous to color space conversion.
  
  # Fold change can be used as an alternative to Z score for the radial scale. z axis shows -
  # log10 p value for likelihood ratio test comparing all three groups. 
  
  # Intial analysis run:
  dds_DE <- DESeq(filter_dds)
  
  # Likelihood ratio test on "Cluster': 
  dds_LTR <- DESeq(filter_dds, test = "LRT", reduced = ~Sex + AgeCat + RINCat, parallel = TRUE) 
  
  # Create 'volc3d' class object for plotting:
  res <- deseq_polar(dds_DE, dds_LTR, 'Cluster', padj.method = "BH")
  
  # Get most significant genes per each group:
  rel_genes <- res@df$scaled # Get a table with coordinates and gene info
  
  # Output result table:
  write.table(rel_genes, 
              file = "../../../output/Volcano3D_table.tsv", 
              sep = "\t", 
              quote = FALSE,
              col.names = NA,
              row.names = TRUE)
  
  rel_genes <- rel_genes[rel_genes$lab != "ns", ] # Dont want to label not significant genes
  rel_genes$symbol <- rownames(rel_genes)
  
  # Get the first 5 most significant genes for each color:
  rel_genes_top <- rel_genes %>%
    group_by(col) %>%
    slice_min(order_by = pvalue, n = 5) %>%
    ungroup()
  
  # Just to higlight the ones associated only with one subtype:
  rel_genes_top <- rel_genes_top[rel_genes_top$lab %in% c("S+", "N+", "O+"), ]
  
  # Then just to double check that selected genes are significant after adjusted p values:
  
  padj <- as.data.frame(res@padj) # Get adjusted p values for each of the comparisons
  padj <- padj %>%
    filter(if_any(everything(), ~ . <= 0.05))
  
  # Get final list of genes: 
  final_padj <- rownames(padj[rownames(padj) %in% rel_genes_top$symbol, ])
  
  # Plot 3d volcano plot:
  # Define color scheme:
  res@scheme <- c("grey60", "green4", "red", "firebrick3", "purple", "royalblue2", "yellow")
  
  # Since Type is not specified. Here Type 1 is used which is Z-score (scaled Fold Change)
  v3D <- 
    volcano3D(res, axis_width = 4, grid_width = 3, z_axis_title_size = 30, radial_axis_title_size = 30)
  
  # Save as Interactive HTML:
  htmlwidgets::saveWidget(as_widget(v3D), "../../../output/Volcano3D.html")
  
  # Volcano3D version from the top:
  # Type 2 indicates that is showing FOLD CHANGE.
  top_3D <- 
    radial_plotly(res, type = 2,
                  axis_width = 4, grid_width = 3, axis_title_size = 30, axis_label_size = 20)
  
  # Save as Interactive HTML:
  htmlwidgets::saveWidget(as_widget(top_3D), "../../../output/Volcano3D_top.html")
  
  # Save table with significant genes for specfic groups: 
  rel_genes_padj <- rel_genes[rel_genes$symbol %in% rownames(padj) & 
                                (rel_genes$lab == "N+" | rel_genes$lab == "O+" |
                                   rel_genes$lab == "S+"), ] # All LTR significant
  rel_genes_padj <- rel_genes_padj[, c(12, 11, 8, 9)]
  
  # Genes that were significant for every subtype: Specific pairwise:
  # A = Neu, B = Oxa, C = SNs.
  
  rel_gens_padj_neu <- rel_genes[rel_genes$symbol %in% rownames(padj[padj$AvB <= 0.05 & padj$AvC <= 0.05, ]), ]
  rel_gens_padj_oxa <- rel_genes[rel_genes$symbol %in% rownames(padj[padj$AvB <= 0.05 & padj$BvC <= 0.05, ]), ]
  rel_gens_padj_sns <- rel_genes[rel_genes$symbol %in% rownames(padj[padj$AvC <= 0.05 & padj$BvC <= 0.05, ]), ]
  
  rel_gens_padj_pair <- rbind(rel_gens_padj_neu, rel_gens_padj_oxa, rel_gens_padj_sns)
  rel_gens_padj_pair <- rel_gens_padj_pair[rel_gens_padj_pair$lab == "N+" | rel_gens_padj_pair$lab == "O+" |
                                             rel_gens_padj_pair$lab == "S+", ]
  rel_gens_padj_pair <- rel_gens_padj_pair[, c(12, 11, 9)]
  
  # Output tables:
  write.table(rel_genes_padj, 
              file = "../../../output/V3D_genes_padj_LTR.tsv", 
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  write.table(rel_gens_padj_pair, 
              file = "../../../output/V3D_genes_padj_allcomp.tsv", 
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  # -----------------------> ENRICHMENT ANALYSIS: 
  
  # The list of informative genes for each cluster was then used 
  # to characterise their molecular phenotypes by performing gene enrichment 
  # analysis using the clusterProfiler and ReactomePA R package. 
  # Genes from the whole expression matrix were used as a custom gene background. 
  # Evaluated databases:
  # Gene Ontology (Biological Process (GO:BP), Molecular Function (GO:MF) and Cellular Component (GO:CC)).
  # Reactome.
  
  # ------> LTR Genes:
  
  # Get the significant genes from the volcano 3D plot analysis: 
  rel_cl_genes <- rel_genes_padj[, c(1,2)]
  
  # Get the genes in the useful format: 
  rel_inf_genes <- split(rel_cl_genes$symbol, rel_cl_genes$lab)
  rel_inf_genes <- rel_inf_genes[lengths(rel_inf_genes) > 0] # Just get factos with values
  
  # -----> ENRICHMENT analysis:
  
  # Enrichment analysis REL GENES:
  go_enrich <- lapply(names(rel_inf_genes), function(cluster) {
    ego <- enrichGO(
      gene = rel_inf_genes[[cluster]],
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "ALL",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      universe = bck_genes,
      minGSSize = 10,
      maxGSSize = 500)
    
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      simplify(ego)
    } 
  })
  
  # Enrichment analysis REL GENES:
  # Get background genes in the right nomenclature:
  bck_entriz <- bitr(bck_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  
  # Run Reactome enrichment anaylysis:
  react_enrich <- lapply(names(rel_inf_genes), function(cluster) {
    
    entrez_id <- bitr(rel_inf_genes[[cluster]], fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = "org.Hs.eg.db")
    
    enrichPathway(
      gene          = entrez_id$ENTREZID,
      universe      = bck_entriz$ENTREZID,
      organism      = "human",
      minGSSize     = 10,
      maxGSSize     = 500,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      readable      = TRUE)})
  
  for (n in 1:length(go_enrich)) {
    
    # Get the table with enrich results: 
    if (!is.null(go_enrich[[n]])) {go_results <- go_enrich[[n]]@result } else {go_results <- NULL}
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
      go_results <- go_results[go_results$Count >= 3, ]
      
      write.table(go_results, 
                  file = paste("../../../output/GO_Cluster_", names(rel_inf_genes)[[n]], ".tsv", sep = ""), 
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
        ggtitle(paste("Cluster ", names(rel_inf_genes)[[n]], sep = "")) +
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
      
      png(file = paste("../../../output/GO_Cluster_", names(rel_inf_genes)[[n]], ".png", sep = ""), 
          width = 1500, height = 800)
      
      print(go_plot)
      
      dev.off()
      
    }
    
    if (!is.null(react_results) & nrow(react_results[react_results$p.adjust <= 0.05, ]) > 0) { # only generate plot if there is info in the table
      
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
      react_results <- react_results[react_results$p.adjust <= 0.05, ]
      react_results$final_term <- paste("React", react_results$Description)
      react_results$generatio <- sapply(react_results$GeneRatio, function(x) eval(parse(text = x)))
      react_results$logpval <- -log10(react_results$p.adjust)
      react_results <- react_results[, c(5:7,4)]
      colnames(react_results)[c(2:4)] <- c("GeneRatio", "-Log(Adj P)", "Count")
      react_results <- react_results[order(react_results$`-Log(Adj P)`, decreasing = TRUE), ]
      react_results <- react_results[react_results$Count >= 3, ]
      
      write.table(react_results, 
                  file = paste("../../../output/REACT_Cluster_", names(rel_inf_genes)[[n]], ".tsv", sep = ""), 
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
        ggtitle(paste("Cluster ", names(rel_inf_genes)[[n]], sep = "")) +
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
      
      png(file = paste("../../../output/REACT_Cluster_", names(rel_inf_genes)[[n]], ".png", sep = ""), 
          width = 1500, height = 800)
      
      print(react_plot)
      
      dev.off()
      
    }
    
    if (!is.null(go_results) & !is.null(react_results)) {
    
    # Combined table and plots:
    go_and_react <- rbind(go_paper, react_paper)
    
    # Save all table: 
    write.table(go_and_react, 
                file = paste("../../../output/All_G&R_Cluster_inf_", names(rel_inf_genes)[[n]], ".tsv", sep = ""), 
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
    
    # Prepare table for plot:
    go_and_react <- go_and_react[order(go_and_react$`-Log(Adj P)`, decreasing = TRUE), ] # Order by pvalue
    go_and_react <- go_and_react[go_and_react$Database == "BP" | go_and_react$Database == "React", ] # Keep just BP and Reactome
    go_and_react <- go_and_react[go_and_react$Database == "React", ] # Keep Reactome
    
    # Get the enrich plot:
    go_and_react <- go_and_react[go_and_react$Count >= 20, ]
    go_and_react <- go_and_react[go_and_react$ID != "R-HSA-198933", ]
    
    # Terms you want to highlight
    highlight_terms <- c("React Neuronal System",
      "React Neurotransmitter receptors and postsynaptic signal transmission",
      "React Transmission across Chemical Synapses")
    
    # Create colored labels
    go_and_react$label_col <- ifelse(
      go_and_react$`Final Term` %in% highlight_terms,
      paste0("<span style='color:firebrick3;'>", go_and_react$`Final Term`, "</span>"),
      go_and_react$`Final Term`)
    
    # Create the plot:
    all_plot <-
      ggplot(go_and_react, aes(x = `-Log(Adj P)`, y = reorder(label_col, `-Log(Adj P)`))) +
      geom_point(aes(size = Count, color = `Gene Ratio`)) +
      scale_color_gradient(low = "red", high = "blue", name = "Gene Ratio") +
      scale_size_continuous(range = c(1, 10), breaks = seq(0, max(react_results$Count), by = 20), name = "Count") +
      xlim(NA, max(go_and_react$`-Log(Adj P)`) + 0.2) +
      labs(y = NULL, x = "-log(Adj P)") +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
      # ggtitle(paste("Cluster ", n, " Inf ", names(sample_list)[i], sep = "")) +
      theme(axis.title = element_text(size = 55, hjust = 0.5),
            axis.text.x  =  element_text(size = 45, colour = "black"),
            axis.text.y  =  ggtext::element_markdown(size = 15, colour = "black", lineheight = 0.6),
            legend.title = element_text(size = 45),
            legend.text = element_text(size = 30),
            legend.position = "right",
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black")) 
    
    png(file = paste("../../../output/All_G&R_Cluster_inf_", names(rel_inf_genes)[[n]], ".png", sep = ""), 
        width = 1500, height = 800)
    
    print(all_plot)
    
    dev.off()
    
    }
    
  }
  
    
    
  
  
