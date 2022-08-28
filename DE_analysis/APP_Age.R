
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


library(statmod)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(purrr)
library(tidyverse)
library(biomaRt)
library("pheatmap")
library(dendextend)
library(RColorBrewer)
library('fgsea')
library(dplyr)
library(EnhancedVolcano)
library(ggrepel)
library(ggfortify)


ensembl_id_to_gene_symbol <- function(values) {
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
    mgi_ids <- getBM(attributes=c("ensembl_gene_id","mgi_symbol"), filter="ensembl_gene_id", value=values, mart=mouse, uniqueRows=T)
    return(data.frame(mgi_ids))
}   

run_deseq_APP_Age <- function(cnts, samples,  count_filter, condition) {
  # Convert col types to factors
  sample_cols <- colnames(samples)
  samples[sample_cols] <- lapply(samples[sample_cols], factor)  
  class(samples[sample_cols])
  #summary(samples)
  samples
  
  # reorder rownames of metadata table and colnames of counts file so that they are
  # in the same order. Otherwise DESeqDataSet will throw an error.
  samples2 <- samples[order(rownames(samples)), ]
  cnts2 <- cnts[, order(colnames(cnts))]
  
  # create the DESeq object
  dds <- DESeqDataSetFromMatrix(countData=cnts2, colData=samples2, design= ~ Age + APP) 
  keep <- rowSums(counts(dds)) >= count_filter
  dds <- dds[keep,]
  dds$APP <- relevel(dds$APP, ref = "WT")
  dds <- DESeq(dds)
  #Check which comparisions were run
  resultsNames(dds)
  message("resultsNames:")
  print(resultsNames(dds))
  # Generate Gene Expression table
  results <- results(dds)
  summary(results)
  #differential expression results sorted by adjusted p-value.
  results <- results[order(results$padj), ]
  return(list('res' = results, 'dds' = dds))
}

# Histogram
deseq_FC_hist_plot <- function(deseq_result, title) {
  #histogram
  log2fc_plot <- ggplot(deseq_result, aes(log2FoldChange))+
    geom_histogram(color="black", fill="lightblue", bins=100) +
    theme_bw() +
    xlab("log2FoldChange") + ylab("count") +
    labs(title = title) +
    theme(plot.title = element_text(size = 8))
  return(log2fc_plot)
  
}
# Scatter plot
deseq_DE_scat_plot <- function(deseq_result, title) {
  #Scatter Plot
  scat_plot <- ggplot(deseq_result, aes(x=log2FoldChange, y=pvalue)) + 
    geom_point(size=1.0) + 
    theme_bw() +
    xlab("log2FoldChange") + ylab("pvalue") +
    labs(title = title) +
    theme(plot.title = element_text(size = 8))
  
  return(scat_plot)
}  


plot_enhanced_volcano <- function(deseq_res_with_genename, padj_threshold,  FC_threshold, plot_subtitle) {
  labeled_results <- label_res(deseq_res_with_genename, padj_threshold, FC_threshold)
  labeled_results <- labeled_results[order(labeled_results$padj, decreasing = FALSE), ]
  labeled_results %>% relocate(genes, volc_plot_status, log2FoldChange, padj)
  Up_genes <- filter(labeled_results, (labeled_results$volc_plot_status == 'UP'))
  Down_genes <- filter(labeled_results, (labeled_results$volc_plot_status == 'DOWN'))
  no_of_UP <- nrow(Up_genes)
  no_of_DOWN <- nrow(Down_genes)
  no_of_NS <- nrow(labeled_results) - no_of_UP - no_of_DOWN
  toptable <- labeled_results
  
  # create custom key-value for 'Up', 'Down' and 'NS' Genes
  keyvals <- ifelse(
    toptable$log2FoldChange > FC_threshold & toptable$padj < padj_threshold, '#00BFC4',
    ifelse(toptable$log2FoldChange < FC_threshold & toptable$padj < padj_threshold, '#F8766D',
           'blue'))
  keyvals[is.na(keyvals)] <- 'blue'
  names(keyvals)[keyvals == '#00BFC4'] <- 'Up Regulated'
  names(keyvals)[keyvals == 'blue'] <- 'Not Significant'
  names(keyvals)[keyvals == '#F8766D'] <- 'Down Regulated'
  
  title <- "Enhanced Volcano Plot"
  caption <- paste0("\nTotal Number of Genes:  ", nrow(labeled_results),
                    "\nNumber of Non Significant Genes:  ", no_of_NS,
                    "\nNumber of Up Regulated Genes Significant at padj < ", padj_threshold, " and log2FC > ", FC_threshold, ":  ", no_of_UP, 
                    "\nNumber of Down Regulated Genes Significant at padj < ", padj_threshold, " and log2FC < ", FC_threshold, ":  ", no_of_DOWN)
  
  EnhancedVolcano_plot <- 
    EnhancedVolcano(toptable = toptable, 
                    lab = toptable$mgi_symbol, 
                    x = 'log2FoldChange', 
                    y = 'padj',
                    xlim = c(min(toptable$log2FoldChange, na.rm = TRUE) - 1.5, max(toptable$log2FoldChange, na.rm = TRUE) + 1.5),
                    ylim = c(0, max(-log10(toptable$padj), na.rm = TRUE) + 5),
                    xlab = bquote(~Log[2] ~ "FoldChange"),
                    ylab = bquote(~-Log[10] ~ (padj)),
                    axisLabSize = 12,
                    pCutoff = padj_threshold,
                    FCcutoff =  FC_threshold,
                    pointSize = 0.30,
                    labSize = 2.0,
                    colCustom = keyvals,
                    colAlpha = 1,
                    title = title,
                    subtitle = subtitle,
                    caption = caption,
                    titleLabSize = 11,
                    subtitleLabSize = 8,
                    captionLabSize = 8,
                    legendPosition = 'right',
                    legendLabSize = 9,
                    legendIconSize = 3.0,
                    drawConnectors = TRUE,
                    widthConnectors = 0.25,
                    colConnectors = "grey10",
                    directionConnectors = 'both',
                    max.overlaps = 15
    )
  return(EnhancedVolcano_plot)
}


label_res <- function(deseq2_res, padj_threshold, FC_threshold) {
  #Using the results generated by DESeq2, convert it to a `tibble` and add one additional column that 
  #denotes whether a gene is either: 1. Upregulated and significant at a given padj threshold, 
  #2.Downregulated and significant at a given padj threshold, or 3. Not significant at a given 
  #padj threshold. Name this column `volc_plot_status`. Make sure the values for these labels are 
  #UP, DOWN, and NS, respectively. Have your function take the DESeq2 results as well as a
  #user-defined threshold for padj and set it at padj_threshold.
  
  labeled <- deseq2_res %>%
    as_tibble(rownames='genes') %>%
    mutate(volc_plot_status = case_when(log2FoldChange > FC_threshold & padj < padj_threshold ~ 'UP', 
                                        log2FoldChange < FC_threshold & padj < padj_threshold ~ 'DOWN', 
                                        TRUE ~ 'NS'))
  return(labeled)
  
}

# Generate heatmap with both Ensemble Ids and Gene names (mgi_symbols)
generate_heatmap <- function(dds, results, title){
  dds_counts <- counts(dds,normalized=TRUE)
  dds_counts <- as.matrix(dds_counts)
  ## can now scale the rows via z-score
  z_score <- function(x) {(x- mean(x)) / sd(x)}
  ## Create list of top 50 DEGs from the deseq results 
  top50 <- row.names(results[1:50,])
  ## now subset the dds counts so that we only have the top 50 genes remaining
  top50_dds <- dds_counts[rownames(dds_counts) %in% top50, ]
  top50_dds_norm <- t(apply(top50_dds, 1, z_score))
  cols <- colData(dds)[,c("APP", "Age")]
  df <- as.data.frame(cols)
  ## Now perform hierarchical clustering to obtain gene cluseters
  hclust_gene <- hclust(dist(top50_dds_norm), method = "complete")
  # Now we can form 2 clusters using cutree() function
  gene_col <- cutree(tree = as.dendrogram(hclust_gene), k=2 )
  gene_col <- data.frame(Cluster= ifelse(test = gene_col ==1, yes = "Cluster 1", no= "Cluster 2"))
  
  # pheatmap with Gene Names/mgi_symbol
  #Add gene name to deseq results
  top50_dds_norm_genename <- data.frame(top50_dds_norm)
  # Convert rownames (which are the geneids) to a column
  top50_dds_norm_genename <- rownames_to_column(top50_dds_norm_genename, "Geneid")
  # add a new col of mgi_symbols and relocate the mgi_symbol col to first col
  # Following is the resulting col names samples
  # "mgi_symbol", "Geneid", "X102", "X103", "X106", "X107", "X117", "X122", "X124", "X125", "X127", "X129", 
  # "X130", "X132", "X134", "X139", "X141", "X147", "X149", "X152", "X153", "X155", "X157", "X160", etc.,
  top50_dds_norm_genename$mgi_symbol <- results$mgi_symbol[match(top50_dds_norm_genename$Geneid, rownames(results))]
  top50_dds_norm_genename <- top50_dds_norm_genename %>% relocate(mgi_symbol)
  # now make the mgi_symbols as rownames
  row.names(top50_dds_norm_genename) <- top50_dds_norm_genename$mgi_symbol
  # we no longer need the mgi_symbol and Geneid columns so lets get rid of them
  top50_dds_norm_genename <- subset(top50_dds_norm_genename, select = -c(1,2))
  # lastly we need to remove character X from all the column names
  colnames(top50_dds_norm_genename) <- gsub("^X", "",  colnames(top50_dds_norm_genename))
  
  # Now perform hierarchical clustering to obtain gene cluseters
  hclust_gene_2 <- hclust(dist(top50_dds_norm_genename), method = "complete")
  
  # Form 2 clusters using cutree() function
  gene_col_2 <- cutree(tree = as.dendrogram(hclust_gene_2), k=2 )
  gene_col_2 <- data.frame(Cluster= ifelse(test = gene_col_2 ==1, yes = "Cluster 1", no= "Cluster 2"))
  
  # Heatmap with Ensemble Id
  # pheat <- pheatmap(top50_dds_norm, annotation_row = gene_col, annotation_col=df, fontsize_col=3.5, fontsize=5.5, main=title)
  # Heatmap with Gene names/ mgi_symbol
  pheat_map <- pheatmap(top50_dds_norm_genename, annotation_row = gene_col_2, annotation_col=df,
                        fontsize_col=4, fontsize_row=4.8, fontsize=7, main=title)
  
  return(pheat_map)
}

# Function to run fgsea on DESeq2 results
run_gsea <- function(labeled_results, gmt) {
  # Strip the digit extensions in the gene names in labled_results and pull the 
  # updated/modified gene names to be used in getBM to retrieve information
  labeled_results_modified <- separate(labeled_results, genes, sep='\\.', into='genes', remove=TRUE)
  # Pull the modified gene names
  genes <- labeled_results_modified$genes
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mgi_symbols <- getBM(attributes=c("ensembl_gene_id","mgi_symbol"), filter="ensembl_gene_id", value=genes, mart=mouse, uniqueRows=T)
  # Add mgi symbols to modified labeled_results
  mgi_results <- left_join(labeled_results_modified, mgi_symbols, by=c('genes' = 'ensembl_gene_id'))
  
  # Create ranked list of log2FoldChange values
  #log2FC <- filter(mgi_results, !is.na(HGNC.symbol) & !is.na(log2FoldChange))
  log2FC <- filter(mgi_results, !is.na(mgi_symbol) & !is.na(log2FoldChange))
  log2FC <- log2FC %>% arrange(desc(log2FoldChange))
  log2FC_rank <- log2FC[c("mgi_symbol", "log2FoldChange")]
  log2FoldChange_rank <- deframe(log2FC_rank)
  # Get M2 Canonical Pathways gene set collection
  gmt_pathways <- gmtPathways(gmt)
  # Run gsea on ranked log2FoldChange values
  fgsea_results <- fgsea(gmt_pathways, log2FoldChange_rank)
  fgsea_results <- as_tibble(fgsea_results)
  return(fgsea_results)
}

# Function to plot top ten positive NES and top ten negative NES pathways
# in a barchart
top_pathways <- function(fgsea_results, num_paths, title){
  #fgsea_results is already sorted in the aescending order of NES
  # make a tibble with only pathways with top 10 positive and top 10 negative NES values
  # gather top 10 negative values
  top_neg <- fgsea_results[1:num_paths,]
  # gather top 10 positive values
  start_pos <- nrow(fgsea_results) - (num_paths-1)
  stop_pos  <- nrow(fgsea_results)
  top_pos <- fgsea_results[start_pos:stop_pos,]
  # filter out the rows wtih top 10 pos and neg valued pathways
  pos_pathway <- top_pos$pathway
  neg_pathway <- top_neg$pathway
  fgsea_subset <- fgsea_results %>% filter(pathway %in% c(pos_pathway, neg_pathway)) 
  # In fgsea_subset create a new col of pathway names (x_axis_names) without '-' 
  # in the name and reorder them for x-axis lables 
  fgsea_subset <- fgsea_subset %>% mutate(x_axis_names = str_replace_all(pathway, '_', ' ')) %>% 
    mutate(x_axis_names = forcats::fct_reorder(factor(x_axis_names), NES))
  
  # Create a bar plot of top 10 pos and neg NES pathways
  nes_plot <- 
    ggplot(fgsea_subset) +
    geom_bar(aes(x=x_axis_names, y=NES, fill = NES > 0), stat='identity', show.legend = FALSE) +
    scale_fill_manual(values = c('TRUE' = 'blue', 'FALSE' = 'red')) + 
    theme_minimal(base_size = 8) +
    ggtitle(label= "fgsea results for M2 CP gene sets derived from: \nWikiPathways database (m2.cp.wikipathways.v0.3.symbols.gmt)",
            subtitle = title) +
    theme(plot.title = element_text(hjust=0.5, color="blue", size=8, face="bold"),
          plot.subtitle = element_text(hjust=0.5, color="black", size=6, face="bold")) +
    ylab('Normalized Enrichment Score (NES)') +
    xlab('') 
  
  # Since some of the the pathway names are very long let us wrap the names so that we can
  # display the bar plot correctly
  wrap_pathway_name <- function(x) str_wrap(x, width = 80)
  nes_plot <- nes_plot +
    scale_x_discrete(labels = wrap_pathway_name) 
  
  # Filp the plot from vertical to horizontal layout
  nes_plot <- nes_plot + coord_flip()
  return(nes_plot)
}

# Define a function to calculate the proportion of variance explained by each PC
calculate_variance_explained <- function(pca_results) {
  variance <- pca_results$sdev^2 / sum((pca_results$sdev)^2)
  return(variance)
}

# Define a function that takes in the variance values and the PCA results to
# make a tibble with PCA names, variance explained by each PC, and the
# cumulative sum of variance explained
make_variance_tibble <- function(pca_ve, pca_results) {
  pca_cum <- cumsum(pca_results$sdev^2 / sum(pca_results$sdev^2))
  variance_explained = c(pca_ve)
  cumulative = c(pca_cum)
  principal_components = c(colnames(pca_results$rotation))
  var_tibble <- tibble(variance_explained, principal_components, cumulative)
  return(var_tibble)
}

# Define a function that creates a bar plot of the variance explained by each
# PC along with a scatter plot showing the cumulative sum of variance explained
# using ggplot2
plot_pca_variance <- function(variance_tibble) {
  variance_df <- data.frame(variance_tibble)
  # we need to order the pincipal components according to nuemerical order (PC1, PC2, PC3 ...)
  # instead of alphabetical order(PC1, PC10, PC11,...,PC2...)
  variance_tibble$principal_components <- 
    factor(variance_tibble$principal_components, levels=c(variance_tibble$principal_components))
  
  bplot <- 
    ggplot(variance_tibble, aes(x=principal_components)) +
    # Need geom_point and geom_line to draw a dot conneced by line plot
    geom_point(mapping=aes(y=cumulative, colour="Cumulative", group=1)) + 
    geom_line(mapping=aes(y=cumulative, colour="Cumulative", group=1)) + 
    # set color and legend text for colour parameter used in geom_point and geom_line
    scale_colour_manual(values="black", name="Cumulative") +
    geom_bar(mapping=aes(y=variance_explained, fill="Variance Explained", 
                         group=variance_explained), colour="#006000", stat='identity') +
    # set color and legend text for fill parameter used in geom_bar above
    scale_fill_manual(values="lightblue", name="Variance Explained") +
    theme_bw() +
    xlab("PC") + 
    ylab("% variance") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.50, hjust = 0.50, size=4)) +
    ylim(0, 1)
  
  return(bplot)
}

# Define a function to create a biplot of PC1 vs. PC2 labeled by
# APP(APP vs WT)
make_biplot <- function(metadata, pca_results, title) {
  meta_data <- read.csv(metadata)
  meta_data <- data.frame(meta_data)
  # select only Sample and APP which is the GSM sample data from meta_data
  meta_df <- meta_data %>% dplyr::select(APP, Sample)   
  # filter the row names of meta_df using row names of pca_results$x
  meta_df_sample <- filter(meta_df, meta_df[,"Sample"] %in% rownames(pca_results$x))
  
  # Now change the row names in pca_results$x to APP names from meta data
  # which we will use to color code the PC1 and PC2 data in the plotting
  row.names(pca_results$x) <- meta_df_sample[,"APP"]
  APP = c(rownames(pca_results$x))
  s <- summary(pca_results)

  bplot <- 
    pca_results$x %>% as.data.frame %>%
    ggplot(aes(x=PC1,y=PC2, col=APP)) + 
    geom_point() + 
    #theme_classic() +
    theme_bw() +
    labs(title=title) +
    theme(plot.title = element_text(size = 9)) +
    geom_hline(yintercept = 0, linetype=2, colour="grey50") +
    geom_vline(xintercept = 0, linetype=2, colour="grey50") +
    xlab(paste0("PC1 (", round(s$importance[2,1]*100, 1), "%)")) + 
    ylab(paste0("PC2 (", round(s$importance[2,2]*100, 1), "%)"))
  
  return(bplot)
}


