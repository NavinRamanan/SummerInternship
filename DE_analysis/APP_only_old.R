
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')

#BiocManager::install('EnhancedVolcano')

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
#library(ggfortify)

ensembl_id_to_gene_symbol <- function(values) {
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
    mgi_ids <- getBM(attributes=c("ensembl_gene_id","mgi_symbol"), filter="ensembl_gene_id", value=values, mart=mouse, uniqueRows=T)
    return(data.frame(mgi_ids))
}   

run_deseq_APP_only <- function(cnts, samples,  count_filter, condition) {
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
  dds <- DESeqDataSetFromMatrix(countData=cnts2, colData=samples2, design= ~ APP) 
  keep <- rowSums(counts(dds)) >= count_filter
  dds <- dds[keep,]
  dds$APP <- relevel(dds$APP, ref = "WT")
  dds <- DESeq(dds)
  #Check which comparisions were run
  resultsNames(dds)
  #view(resultsNames(dds))
  message("resultsNames:")
  print(resultsNames(dds))
  # Generate Gene Expression table
  results <- results(dds)
  #view(results)
  #results <- lfcShrink(dds, contrast = c('dex','trt','untrt'), res=results, type = 'normal')
  #results <- lfcShrink(dds, coef=2)
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
    theme(plot.title = element_text(size = 9))
  return(log2fc_plot)
  
}
# Scatter plot
#Cortex_deseq_scat_plot1 <- deseq_DE_scat_plot(Cortex_deseq_results_APP_padj_filtered, Cortex_deseq_scat_title1)
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


plot_volcano <- function(deseq2_results, padj_threshold, FC_threshold, plot_title) {

  labeled_results <- label_res(deseq2_results, padj_threshold, FC_threshold)
  labeled_results <- labeled_results[order(labeled_results$padj, decreasing = FALSE), ]
  labeled_results %>% relocate(genes, volc_plot_status, log2FoldChange, padj)
  #view(labeled_results)
  Up_genes <- filter(labeled_results, (labeled_results$volc_plot_status == 'UP'))
  Down_genes <- filter(labeled_results, (labeled_results$volc_plot_status == 'DOWN'))
  no_of_UP <- nrow(Up_genes)
  no_of_DOWN <- nrow(Down_genes)
  no_of_NS <- nrow(labeled_results) - no_of_UP - no_of_DOWN
  
  cols <- c("UP"="#00BFC4", "DOWN"="#F8766D", "NS"="blue")
  new_label <- c("Up Regulated", "Down Regulated", "Not Significant")
  caption <- paste0("\nTotal Number of DE Genes:  ", nrow(labeled_results),
                    "\nNumber of Non Significant DE Genes:  ", no_of_NS,
                    "\nNumber of Up Regulated Genes Significant at padj < ", padj_threshold, " and log2FC < ", FC_threshold, ":  ", no_of_UP, 
                    "\nNumber of Down Regulated Genes Genes Significant at padj < ", padj_threshold, " and log2FC < ", FC_threshold, ":  ", no_of_DOWN)
  #############################################################
  volcano_scat_plot <- ggplot(labeled_results) + 
    geom_point(aes(x=log2FoldChange, y=-log10(padj), color=volc_plot_status), size=0.9) + 
    scale_color_manual(name="Gene Status", values = cols, label = new_label) +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed")  + 
    theme_minimal() +
    labs(title = plot_title, subtitle = caption) +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 9, color = "blue")) +
    xlab("log2FoldChange") + 
    ylab("-log10(padj)")
  
  #Following code will show the colors used by ggplot
  #view(ggplot_build(volcano_scat_plot)$data[[1]])
  #colour         x         y PANEL group shape size fill alpha stroke
  #1   #619CFF 7.6359387 195.90690     1     3    19  0.9   NA    NA    0.5
  #2   #619CFF 6.2351925 162.77491     1     3    19  0.9   NA    NA    0.5
  #3   #619CFF 6.4682492 146.89611     1     3    19  0.9   NA    NA    0.5
  
  message("")
  message("**************************************************")
  message(paste0("Total Number of DE Genes:   ", nrow(labeled_results)))
  message(paste0("Number of UP Regulated Genes:   ", no_of_UP))
  message(paste0("Number of DOWN Regulated Genes: ", no_of_DOWN))
  message(paste0("Number of Non Significant Genes: ", no_of_NS))
  message("**************************************************")
  message("")
  
  return(volcano_scat_plot)
}

plot_volcano_all <- function(deseq2_results, padj_threshold, FC_threshold, plot_title) {
  #labeled_results <- label_res(deseq2_results, padj_threshold)
  labeled_results <- deseq2_results %>%
    as_tibble(rownames='genes') %>%
    mutate(volc_plot_status = case_when(log2FoldChange > FC_threshold & padj < padj_threshold ~ 'UP',
                                        log2FoldChange < FC_threshold & padj < padj_threshold ~ 'DOWN',
                                        TRUE ~ 'NS',))
  
  labeled_results <- labeled_results[order(labeled_results$padj, decreasing = FALSE), ]
  labeled_results %>% relocate(genes, volc_plot_status, log2FoldChange, padj)
  
  #view(labeled_results)
  Up_genes <- filter(labeled_results, (labeled_results$volc_plot_status == 'UP'))
  Down_genes <- filter(labeled_results, (labeled_results$volc_plot_status == 'DOWN'))
  no_of_UP <- nrow(Up_genes)
  no_of_DOWN <- nrow(Down_genes)
  no_of_NS <- nrow(labeled_results) - no_of_UP - no_of_DOWN
  
  cols <- c("UP"="#00BFC4", "DOWN"="#F8766D", "NS"="blue")
  new_label <- c("Up Regulated", "Down Regulated", "Not Significant")
  caption <- paste0("\nTotal Number of DE Genes:  ", nrow(labeled_results),
                    "\nNumber of Non Significant DE Genes:  ", no_of_NS,
                    "\nNumber of Up Regulated Genes Significant at padj < ", padj_threshold, " and log2FC < ", FC_threshold, ":  ", no_of_UP, 
                    "\nNumber of Down Regulated Genes Genes Significant at padj < ", padj_threshold, " and log2FC < ", FC_threshold, ":  ", no_of_DOWN)
  #############################################################
  volcano_scat_plot <- ggplot(labeled_results) + 
    geom_point(aes(x=log2FoldChange, y=-log10(padj), color=volc_plot_status), size=0.9) + 
    geom_hline(yintercept = -log10(0.1), linetype = "dashed")  + 
    theme_minimal() +
    scale_color_manual(name="Gene Status", values = cols, label = new_label) +
    labs(title = plot_title, subtitle = caption) +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 9, color = "blue")) +
    xlab("log2FoldChange") + 
    ylab("-log10(padj)")
  
  #Following code will show the colors used by ggplot
  #view(ggplot_build(volcano_scat_plot)$data[[1]])
  #colour         x         y PANEL group shape size fill alpha stroke
  #1   #619CFF 7.6359387 195.90690     1     3    19  0.9   NA    NA    0.5
  #2   #619CFF 6.2351925 162.77491     1     3    19  0.9   NA    NA    0.5
  #3   #619CFF 6.4682492 146.89611     1     3    19  0.9   NA    NA    0.5
  
  message("")
  message("**************************************************")
  message(paste0("Total Number of DE Genes:   ", nrow(labeled_results)))
  message(paste0("Number of UP Regulated Genes:   ", no_of_UP))
  message(paste0("Number of DOWN Regulated Genes: ", no_of_DOWN))
  message(paste0("Number of Non Significant Genes: ", no_of_NS))
  message("**************************************************")
  message("")
  
  return(volcano_scat_plot)
}

plot_volcano_all_genename_4 <- function(deseq2_results, padj_threshold,  FC_threshold, plot_title) {
  labeled_results <- label_res(deseq2_results, padj_threshold, FC_threshold)
  labeled_results <- labeled_results[order(labeled_results$padj, decreasing = FALSE), ]
  labeled_results %>% relocate(mgi_symbol, volc_plot_status, log2FoldChange, padj)
  
  Up_genes <- filter(labeled_results, (labeled_results$volc_plot_status == 'UP'))
  Down_genes <- filter(labeled_results, (labeled_results$volc_plot_status == 'DOWN'))
  no_of_UP <- nrow(Up_genes)
  no_of_DOWN <- nrow(Down_genes)
  no_of_NS <- nrow(labeled_results) - no_of_UP - no_of_DOWN
  
  UP_and_Down_genes <- filter(labeled_results, (labeled_results$volc_plot_status == 'UP' |
                                                  labeled_results$volc_plot_status == 'DOWN'))
  cols <- c("Up Regulated"="#00BFC4", "Down Regulated"="#F8766D", "Not Significant"="blue")
  caption <- paste0("\nTotal Number of DE Genes:  ", nrow(labeled_results),
                    "\nNumber of Non Significant DE Genes:  ", no_of_NS,
                    "\nNumber of Up Regulated Genes Significant at padj < ", padj_threshold, " and log2FC < ", FC_threshold, ":  ", no_of_UP, 
                    "\nNumber of Down Regulated Genes Genes Significant at padj < ", padj_threshold, " and log2FC < ", FC_threshold, ":  ", no_of_DOWN)
  
  #############################################################
  volcano_scat_plot <- 
    labeled_results %>%
    ggplot(mapping=aes(x=log2FoldChange, y=-log10(padj))) + 
    geom_point(aes(color=volc_plot_status), size=0.5) + 
    geom_text_repel(data = UP_and_Down_genes, aes(label = mgi_symbol, fill=mgi_symbol)) +
    #                size = 2, # font size in the text labels
    #                color = 'black', # text color on label
    #                box.padding = unit(.7, "lines"),
    #                hjust= 0.2,
    #                max.overlaps = 20,
    #                segment.color = 'grey50'
    #) +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size=3))) +
    scale_color_manual(values = c("DOWN"="#F8766D", "NS"="blue", "UP"="#00BFC4")) +
    #scale_color_manual(values = c("Up Regulated"="#00BFC4", "Down Regulated"="#F8766D", "Not Significant"="blue")) +
    labs(title = plot_title,
         subtitle = paste0("Total Number of DE Genes:  ", nrow(labeled_results),
                           "\nNumber of Non Significant DE Genes:  ", no_of_NS,
                           "\nNumber of Up Regulated Genes Significant at padj < ", padj_threshold, ":  ", no_of_UP, 
                           "\nNumber of Down Regulated Genes Genes Significant at padj < ", padj_threshold, ":  ", no_of_DOWN)) + 
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 9, color = "blue")) +
    geom_hline(yintercept = 0, linetype=2, colour="grey50") +
    geom_vline(xintercept = 0, linetype=2, colour="grey50") +
    xlab("log2FoldChange") + 
    ylab("-log10(padj)")
  
  print("plot_volcano_all_genename - 2")
  #volcano_scat_plot <- volcano_scat_plot + geom_text_repel()
  
  #Following code will show the colors used by ggplot
  #view(ggplot_build(volcano_scat_plot)$data[[1]])
  #colour         x         y PANEL group shape size fill alpha stroke
  #1   #619CFF 7.6359387 195.90690     1     3    19  0.9   NA    NA    0.5
  #2   #619CFF 6.2351925 162.77491     1     3    19  0.9   NA    NA    0.5
  #3   #619CFF 6.4682492 146.89611     1     3    19  0.9   NA    NA    0.5
  
  message("")
  message("**************************************************")
  message(paste0("Total Number of DE Genes:   ", nrow(labeled_results)))
  message(paste0("Number of UP Regulated Genes:   ", no_of_UP))
  message(paste0("Number of DOWN Regulated Genes: ", no_of_DOWN))
  message(paste0("Number of Non Significant Genes: ", no_of_NS))
  message("**************************************************")
  message("")
  
  return(volcano_scat_plot)
}

plot_volcano_all_genename <- function(deseq2_results, padj_threshold,  FC_threshold, plot_title) {
  labeled_results <- label_res(deseq2_results, padj_threshold, FC_threshold)
  labeled_results <- labeled_results[order(labeled_results$padj, decreasing = FALSE), ]
  labeled_results %>% relocate(mgi_symbol, volc_plot_status, log2FoldChange, padj)
  
  #view(labeled_results)
  Up_genes <- filter(labeled_results, (labeled_results$volc_plot_status == 'UP'))
  Down_genes <- filter(labeled_results, (labeled_results$volc_plot_status == 'DOWN'))
  no_of_UP <- nrow(Up_genes)
  no_of_DOWN <- nrow(Down_genes)
  no_of_NS <- nrow(labeled_results) - no_of_UP - no_of_DOWN
  UP_and_Down_genes <- filter(labeled_results, (labeled_results$volc_plot_status == 'UP' |
                                                labeled_results$volc_plot_status == 'DOWN'))
  #cols <- c("Up Regulated"="#00BFC4", "Down Regulated"="#F8766D", "Not Significant"="blue")
  cols <- c("UP"="#00BFC4", "DOWN"="#F8766D", "NS"="blue")
  new_label <- c("Up Regulated", "Down Regulated", "Not Significant")
  caption <- paste0("\nTotal Number of DE Genes:  ", nrow(labeled_results),
                    "\nNumber of Non Significant DE Genes:  ", no_of_NS,
                    "\nNumber of Up Regulated Genes Significant at padj < ", padj_threshold, " and log2FC < ", FC_threshold, ":  ", no_of_UP, 
                    "\nNumber of Down Regulated Genes Genes Significant at padj < ", padj_threshold, " and log2FC < ", FC_threshold, ":  ", no_of_DOWN)
  #############################################################
  volcano_scat_plot <- 
    #labeled_results %>%
    ggplot(labeled_results, mapping=aes(x=log2FoldChange, y=-log10(padj))) + 
    geom_point(labeled_results, mapping=aes(color=volc_plot_status), size=0.5) + 
    scale_color_manual(name="Gene Status", values = cols, labels=new_label) +
    #geom_label_repel(data = UP_and_Down_genes, aes(label = mgi_symbol)) + 
    geom_text_repel(data = UP_and_Down_genes, aes(label = mgi_symbol), 
                     size = 3, # font size in the text labels
                     color = 'black', # text color on label
    #                box.padding = unit(.45, "lines"),
                     hjust= 1,
                     max.overlaps = 15,
                     segment.color = 'grey50',
                     nudge_y = -3,
                     nudge_x = -3,
                     #force = 1,    #force of repulsion between overlapping text labels
                     #force_pull = 10, #force of attraction between each text label and its data point
                     #direction    = "both", #move text labels “both” (default), “x”, or “y” directions
                     #segment.size = 0.5,
                     #segment.linetype = 1
                     #max.iter = 1000000,
                     #max.time = 10,
                     #min.segment.length = 0
                     #fontface = "bold"
                     #family = "Times"
                     # Repel away from the left edge, not from the right.
                     #xlim = c(NA, Inf),
                     xlim = c(-Inf, NA),
                     # Do not repel from top or bottom edges.
                     #ylim = c(-Inf, Inf)
                     ylim = c(NA, NA)
    
                      ) +
    theme_minimal() +
    guides(color = guide_legend(override.aes = list(size=3))) + # make the size of the legend icon bigger
    labs(title = plot_title, subtitle = caption) + 
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 9, color = "blue"),
          axis.line = element_line(colour = 'black', size = 1)) +
    xlim(c(min(labeled_results$log2FoldChange, na.rm = TRUE) - 1.5, max(labeled_results$log2FoldChange, na.rm = TRUE) + 1.5)) +
    ylim(c(0, max(-log10(labeled_results$padj), na.rm = TRUE) + 5)) +
    geom_hline(yintercept = 0, linetype=2, colour="grey50") +
    geom_vline(xintercept = 0, linetype=2, colour="grey50") +
    xlab("log2FoldChange") + 
    ylab("-log10(padj)")
  
  #Following code will show the colors used by ggplot
  #view(ggplot_build(volcano_scat_plot)$data[[1]])
  #colour         x         y PANEL group shape size fill alpha stroke
  #1   #619CFF 7.6359387 195.90690     1     3    19  0.9   NA    NA    0.5
  #2   #619CFF 6.2351925 162.77491     1     3    19  0.9   NA    NA    0.5
  #3   #619CFF 6.4682492 146.89611     1     3    19  0.9   NA    NA    0.5
  
  message("")
  message("**************************************************")
  message(paste0("Total Number of DE Genes:   ", nrow(labeled_results)))
  message(paste0("Number of UP Regulated Genes:   ", no_of_UP))
  message(paste0("Number of DOWN Regulated Genes: ", no_of_DOWN))
  message(paste0("Number of Non Significant Genes: ", no_of_NS))
  message("**************************************************")
  message("")
  
  return(volcano_scat_plot)
}

plot_volcano_all_genename_2 <- function(deseq2_results, padj_threshold,  FC_threshold, plot_title) {
  #labeled_results <- label_res(deseq2_results, padj_threshold)
  #view(deseq2_results)
  labeled_results <- label_res(deseq2_results, padj_threshold, FC_threshold)
  labeled_results <- labeled_results[order(labeled_results$padj, decreasing = FALSE), ]
  labeled_results %>% relocate(mgi_symbol, volc_plot_status, log2FoldChange, padj)
  
  #view(labeled_results)
  Up_genes <- filter(labeled_results, (labeled_results$volc_plot_status == 'UP'))
  Down_genes <- filter(labeled_results, (labeled_results$volc_plot_status == 'DOWN'))
  no_of_UP <- nrow(Up_genes)
  no_of_DOWN <- nrow(Down_genes)
  no_of_NS <- nrow(labeled_results) - no_of_UP - no_of_DOWN
  UP_and_Down_genes <- filter(labeled_results, (labeled_results$volc_plot_status == 'UP' |
                                                  labeled_results$volc_plot_status == 'DOWN'))
  cols <- c("Up Regulated"="#00BFC4", "Down Regulated"="#F8766D", "Not Significant"="blue")
  caption <- paste0("\nTotal Number of DE Genes:  ", nrow(labeled_results),
                    "\nNumber of Non Significant DE Genes:  ", no_of_NS,
                    "\nNumber of Up Regulated Genes Significant at padj < ", padj_threshold, " and log2FC < ", FC_threshold, ":  ", no_of_UP, 
                    "\nNumber of Down Regulated Genes Genes Significant at padj < ", padj_threshold, " and log2FC < ", FC_threshold, ":  ", no_of_DOWN)
  #############################################################
  volcano_scat_plot <- 
    #labeled_results %>%
    ggplot(labeled_results, mapping=aes(x=log2FoldChange, y=-log10(padj), color=volc_plot_status)) + 
    geom_point(labeled_results, mapping=aes(color=volc_plot_status), size=0.5) + 
    #guides(color = guide_legend(override.aes = list(size=3))) +
    scale_color_manual(name="Gene Status", values = cols) +
    geom_text_repel(data = UP_and_Down_genes, mapping=aes(label = mgi_symbol, fill=mgi_symbol), 
                    size = 2, # font size in the text labels
                    color = 'black', # text color on label
                    box.padding = unit(.45, "lines"),
                    hjust= 0.2,
                    max.overlaps = 15) +
    #segment.color = 'grey50') +
    scale_fill_manual(values = cols) + # Modify point colour
    geom_hline(yintercept = 0, linetype=2, colour="grey50") +
    geom_vline(xintercept = 0, linetype=2, colour="grey50") +
    theme_bw() +
    labs(title = plot_title, subtitle = caption) + 
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 9, color = "blue")) +
    xlim(c(min(labeled_results$log2FoldChange, na.rm = TRUE) - 1.5, max(labeled_results$log2FoldChange, na.rm = TRUE) + 1.5)) +
    ylim(c(0, max(-log10(labeled_results$padj), na.rm = TRUE) + 5)) +
    xlab("log2FoldChange") + 
    ylab("-log10(padj)")
  
  #Following code will show the colors used by ggplot
  #view(ggplot_build(volcano_scat_plot)$data[[1]])
  #colour         x         y PANEL group shape size fill alpha stroke
  #1   #619CFF 7.6359387 195.90690     1     3    19  0.9   NA    NA    0.5
  #2   #619CFF 6.2351925 162.77491     1     3    19  0.9   NA    NA    0.5
  #3   #619CFF 6.4682492 146.89611     1     3    19  0.9   NA    NA    0.5
  
  message("")
  message("**************************************************")
  message(paste0("Total Number of DE Genes:   ", nrow(labeled_results)))
  message(paste0("Number of UP Regulated Genes:   ", no_of_UP))
  message(paste0("Number of DOWN Regulated Genes: ", no_of_DOWN))
  message(paste0("Number of Non Significant Genes: ", no_of_NS))
  message("**************************************************")
  message("")
  
  return(volcano_scat_plot)
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
  
  #view(toptable)
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
                    titleLabSize = 13,
                    subtitleLabSize = 11,
                    captionLabSize = 9,
                    legendPosition = 'right',
                    legendLabSize = 10,
                    legendIconSize = 3.0,
                    drawConnectors = TRUE,
                    widthConnectors = 0.25,
                    colConnectors = "grey10",
                    directionConnectors = 'both',
                    max.overlaps = 15
     )
  return(EnhancedVolcano_plot)
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < padj_threshold and has a positive log
#' fold change, 2. Significant at padj < padj_threshold and has a negative log fold change,
#' 3. Not significant at padj < padj_threshold. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`.
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < padj_threshold.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
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
  #cols <- colData(dds)[,c("APP", "Diet", "Age", "Sex")]
  #df <- as.data.frame(cols)
  APP <- colData(dds)[,"APP"]
  df <- as.data.frame(APP)
  ## Now perform hierarchical clustering to obtain gene cluseters
  hclust_gene <- hclust(dist(top50_dds_norm), method = "complete")
  #as.dendrogram(hclust_gene) %>% plot(horiz= T)
  # Now we can form 2 clusters using cutree() function
  gene_col <- cutree(tree = as.dendrogram(hclust_gene), k=2 )
  gene_col <- data.frame(cluster= ifelse(test = gene_col ==1, yes = "cluster 1", no= "cluster 2"))

  
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
  #as.dendrogram(hclust_gene) %>% plot(horiz= T)

  # Form 2 clusters usingh cutree() function
  gene_col_2 <- cutree(tree = as.dendrogram(hclust_gene_2), k=2 )
  gene_col_2 <- data.frame(cluster= ifelse(test = gene_col_2 ==1, yes = "cluster 1", no= "cluster 2"))
  
  # Heatmap with Gene names/ mgi_symbol
  pheat_map <- pheatmap(top50_dds_norm_genename, annotation_row = gene_col_2, annotation_col=df,
                        #cellwidth=4, cellheight=4, 
                        fontsize_col=4, fontsize=4.8, main=title)
  
  # Heatmap with Ensemble Id
  pheat <- pheatmap(top50_dds_norm, annotation_row = gene_col, annotation_col=df,
                    fontsize_col=3.5, fontsize=5.5, main=title)
  
  return(pheat)
}

#' Function to run fgsea on DESeq2 results
#' @param labeled_results (tibble): the labeled results from DESeq2
#' @param gmt (str): the path to the GMT file
#' @return tibble containing the results from running fgsea using descending
#' log2foldchange as a ranking metric
#' @examples fgsea_results <- run_gsea(labeled_results, 'm2.cp.wikipathways.v0.3.symbols.gmt')
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

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths, title){
  #view(fgsea_results)
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
  #view(fgsea_subset)
  # In fgsea_subset create a new col of pathway names (x_axis_names) without '-' 
  # in the name and reorder them for x-axis lables 
  fgsea_subset <- fgsea_subset %>% mutate(x_axis_names = str_replace_all(pathway, '_', ' ')) %>% 
    mutate(x_axis_names = forcats::fct_reorder(factor(x_axis_names), NES))
  
  fgsea_subset2 <- fgsea_subset
  #names(rownames(fgsea_subset2)) <- paste0(rownames(fgsea_subset2), "-", fgsea_subset2$pval)
  names(fgsea_subset2$pathway) <- paste0(names(fgsea_subset2$pathway), ":", fgsea_subset2$pval)
  view(fgsea_subset2)
  print(fgsea_subset2$pathway)
  view(fgsea_subset)
  
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

top_pathways_padj <- function(fgsea_results, num_paths, title){
  view(fgsea_results)
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
  #view(fgsea_subset)
  # In fgsea_subset create a new col of pathway names (x_axis_names) without '-' 
  # in the name and reorder them for x-axis lables 
  fgsea_subset <- fgsea_subset %>% mutate(x_axis_names = str_replace_all(pathway, '_', ' ')) %>% 
    mutate(x_axis_names = forcats::fct_reorder(factor(x_axis_names), pval))
  
  # Create a bar plot of top 10 pos and neg NES pathways
  padj_plot <- 
    ggplot(fgsea_subset) +
    geom_bar(aes(x=x_axis_names, y=padj, fill = pval > 0), stat='identity', show.legend = FALSE) +
    scale_fill_manual(values = c('TRUE' = 'blue', 'FALSE' = 'red')) + 
    theme_minimal(base_size = 8) +
    ggtitle(label= "fgsea results for M2 CP gene sets derived from: \nWikiPathways database (m2.cp.wikipathways.v0.3.symbols.gmt)",
            subtitle = title) +
    theme(plot.title = element_text(hjust=0.5, color="blue", size=8, face="bold"),
          plot.subtitle = element_text(hjust=0.5, color="black", size=6, face="bold")) +
    xlim(-1,1) +
    ylab('Adjusted p-value') +
    xlab('') 
  
  # Since some of the the pathway names are very long let us wrap the names so that we can
  # display the bar plot correctly
  wrap_pathway_name <- function(x) str_wrap(x, width = 80)
  padj_plot <- padj_plot +
    scale_x_discrete(labels = wrap_pathway_name) 
  
  # Filp the plot from vertical to horizontal layout
  padj_plot <- padj_plot + coord_flip()
  return(padj_plot)
}

# Generate barplot for Total DEGs for each of the test cases run
generate_barplot_for_DEGs <- function(DEGs_table, is_Cortex) {
  # make sure the test variant names are ordered correctly
  Variants <- factor(DEGs_table$Comparison, levels=DEGs_table$Comparison)
  if (is_Cortex == TRUE) {
    Cortex_bplot_title <- "Cortex:Total DE Genes After Padj_threshold < 0.05 filtering"
    Cortex_DEGs <- DEGs_table %>% ggplot() + 
      geom_bar(aes(x=Comparison, y=Cortex_DEGs_after_Padj_filtering, fill=Comparison), stat="identity") +
      #geom_bar(aes(x=Variants, y=Cortex_DEGs_after_Padj_filtering, fill=Variants), stat="identity") +
      geom_text(aes(Comparison, Cortex_DEGs_after_Padj_filtering, label = Cortex_DEGs_after_Padj_filtering), vjust = 0, size=2) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=7), 
            legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 9)) +
      xlab("Variants") + ylab("Number of DE Genes") + labs(title=Cortex_bplot_title)
    return (Cortex_DEGs)
  }
  
  if (is_Cortex == FALSE) {
    Hippocampus_bplot_title <- "Hippocampus:Total DE Genes After Padj_threshold < 0.05 filtering"
    Hippocampus_DEGs <- DEGs_table %>% ggplot() + 
      geom_bar(aes(x=Comparison, y=Hippocampus_DEGs_after_Padj_filtering, fill=Comparison), stat="identity") +
      #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
      geom_text(aes(Comparison, Hippocampus_DEGs_after_Padj_filtering, label = Hippocampus_DEGs_after_Padj_filtering), vjust = 0, size=2) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=7), 
            legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 9)) +
      xlab("Variants") + ylab("Number of DE Genes") + labs(title=Hippocampus_bplot_title)
    return (Hippocampus_DEGs)
  }
}

# Generate barplot for Up Regulated and Down Regulated DEGs
generate_barplot_for_Up_Down_DEGs <- function(DEGs_table, is_Cortex) {
  #print("generate_barplot_for_DEGs")
  # make sure the test variant names are ordered correctly
  Variants <- factor(DEGs_table$Comparison, levels=DEGs_table$Comparison)
  if (is_Cortex == TRUE) {
    Cortex_bplot_title <- "Cortex:Up and Down Regulated DE Genes After Padj_threshold < 0.05 filtering"
    Cortex_Up_Down_DEGs <-
      DEGs_table %>% group_by(Comparison) %>%
      summarise(n=n(),
                Up_Regulated = Cortex_Up_Regulated_DEGs,
                Down_Regulated = Cortex_Down_Regulated_DEGs) %>% 
      gather("DEGs", "value", - c(Comparison, n)) %>%
      ggplot(aes(x = Comparison, y = value, group = DEGs, fill = DEGs, label = value)) + geom_col() +
      geom_text(size = 3, position = position_stack(vjust = 0.5)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=7), 
            legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 9)) +
      xlab("Variants") + ylab("# of Up and Down Regulated DE Genes") + labs(title=Cortex_bplot_title)
    
    return(Cortex_Up_Down_DEGs) 
  }
  
  # This plot request must be for Hippocampus
  if (is_Cortex == FALSE) {
    Hippocampus_bplot_title <- "Hippocampus:Up and Down Regulated DE Genes After Padj_threshold < 0.05 filtering"
    Hippocampus_Up_Down_DEGs <- 
      DEGs_table %>% group_by(Comparison) %>%
      summarise(n=n(),
                Up_Regulated = Hippocampus_Up_Regulated_DEGs,
                Down_Regulated = Hippocampus_Down_Regulated_DEGs) %>% 
      gather("DEGs", "value", - c(Comparison, n)) %>%
      ggplot(aes(x = Comparison, y = value, group = DEGs, fill = DEGs, label = value)) + geom_col() +
      geom_text(size = 3, position = position_stack(vjust = 0.5)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=7), 
            legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 9)) +
      xlab("Variants") + ylab("# of Up and Down Regulated DE Genes") + labs(title=Hippocampus_bplot_title)
    
    return(Hippocampus_Up_Down_DEGs) 
  }
}

#' Define a function to calculate the proportion of variance explained by each PC
#'
#' @param pca_results (obj): the results returned by `prcomp()`
#'
#' @return A vector containing the values of the variance explained by each PC
#' @export
#'
#' @examples
calculate_variance_explained <- function(pca_results) {
  variance <- pca_results$sdev^2 / sum((pca_results$sdev)^2)
  return(variance)
}

#' Define a function that takes in the variance values and the PCA results to
#' make a tibble with PCA names, variance explained by each PC, and the
#' cumulative sum of variance explained
#' @param pca_ve (vector): the vector generated in the previous function with
#'   the variance explained values
#' @param pca_results (object): the results returned by `prcomp()`
#'
#' @return A tibble that contains the names of the PCs, the individual variance
#'   explained and the cumulative variance explained
#' @export
#'
#' @examples 
make_variance_tibble <- function(pca_ve, pca_results) {
  
  #str(summary(pca_results))
  pca_cum <- cumsum(pca_results$sdev^2 / sum(pca_results$sdev^2))
  variance_explained = c(pca_ve)
  cumulative = c(pca_cum)
  principal_components = c(colnames(pca_results$rotation))
  var_tibble <- tibble(variance_explained, principal_components, cumulative)
  return(var_tibble)
}

#' Define a function that creates a bar plot of the variance explained by each
#' PC along with a scatter plot showing the cumulative sum of variance explained
#' using ggplot2
#' @param variance_tibble (tibble): the tibble gnerated in the previous function
#' that contains each PC label, the variance explained by each PC, and the 
#' cumulative sum of variance explained
#'
#' @return A ggplot with a barchart representing individual variance
#'   explained and a scatterplot (connected with a line) that represents the
#'   cumulative sum of PCs
#' @export
#'
#' @examples

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

#' Define a function to create a biplot of PC1 vs. PC2 labeled by
#' SixSubTypesClassification
#'
#' @param metadata (str): The path to the proj_metadata.csv file
#' @param pca_results (obj): The results returned by `prcomp()`
#'
#' @return A ggplot consisting of a scatter plot of PC1 vs PC2 labeled by
#'   SixSubTypesClassification found in the metadata
#' @export
#'
#' @examples

make_biplot <- function(metadata, pca_results, title) {
  meta_data <- read.csv(metadata)
  meta_data <- data.frame(meta_data)
  #view(meta_data)
  # select only SixSubtypesClassification and geo_accession which is the GSM sample data from meta_data
  meta_df <- meta_data %>% dplyr::select(APP, Sample)   
  #meta_df <- meta_data %>% select(Sample)   
  #view(meta_df)
  #str(summary(pca_results))
  #print(rownames(pca_results$x))
  # filter the row names of meta_df using row names of pca_results$x
  meta_df_sample <- filter(meta_df, meta_df[,"Sample"] %in% rownames(pca_results$x))
  #view(meta_df_sample)
  
  # Now change the row names in pca_results$x to SixSubtypesClassification names from meta data
  # which we will use to color code the PC1 and PC2 data in the plotting
  row.names(pca_results$x) <- meta_df_sample[,"APP"]
  APP = c(rownames(pca_results$x))
  s <- summary(pca_results)
  #str(summary(pca_results))
  
  bplot <- 
    pca_results$x %>% as.data.frame %>%
    ggplot(aes(x=PC1,y=PC2, col=APP)) + 
    geom_point() + 
    #theme_classic() +
    theme_bw() +
    labs(title=title) +
    theme(plot.title = element_text(size = 10)) +
    geom_hline(yintercept = 0, linetype=2, colour="grey50") +
    geom_vline(xintercept = 0, linetype=2, colour="grey50") +
    xlab(paste0("PC1 (", round(s$importance[2,1]*100, 1), "%)")) + 
    ylab(paste0("PC2 (", round(s$importance[2,2]*100, 1), "%)"))
  
  return(bplot)
}

make_biplot2 <- function(metadata, pca_results, title) {
  meta_data <- read.csv(metadata)
  meta_data <- data.frame(meta_data)
  #view(meta_data)
  # select only APP and Sample  from meta_data
  meta_df <- meta_data %>% dplyr::select(APP, Sample)   
  #view(meta_df)
  str(summary(pca_results))
  summary(pca_results)

  # filter the row names of meta_df using row names of pca_results$x
  meta_df_sample <- filter(meta_df, meta_df[,"Sample"] %in% rownames(pca_results$x))
  #view(meta_df_sample)
  
  # Now change the row names in pca_results$x to APP names from meta data
  # which we will use to color code the PC1 and PC2 data in the plotting
  #row.names(pca_results$x) <- meta_df_sample[,"APP"]
  #SixSubtypesClassification = c(rownames(pca_results$x))
  #view(pca_results$x)
  #APP = c(rownames(pca_results$x))
  APP <- meta_df_sample[,"APP"]
  #str(summary(pca_results))
  
  #####################
  a <- pca_results$rotation
  a1 <- as.data.frame(a)
  a2 <- rownames_to_column(a1, var="genes")
  # a %>% as.data.frame %>% rownames_to_column %>% dplyr::select(rowname, PC1, PC2) %>% arrange(desc(PC1^2+PC2^2)) %>% 
  #  head(10)
  a3 <- a2 %>% dplyr::select(genes, PC1, PC2)
  a4 <- a3 %>% arrange(desc(PC1^2+PC2^2))
  #view(a4)
  #####################
  s <- summary(pca_results)
  bplot <- 
    pca_results$x %>% as.data.frame %>%
      ggplot(aes(x=PC1,y=PC2, col=APP)) + 
      geom_point() + 
      geom_text_repel(
        aes(#label =  rownames(pca_results$x), fill = factor(APP)),
            label =  rownames(pca_results$x)),
            color = 'black', # text color on label
            size = 3,
            #box.padding = unit(0.25, "lines"), # how far box is away from point + "lines" means connect a line between point and label
            max.overlaps = 20,
            segment.color = 'grey50' # color of said line
            ) +
      theme_bw() +
      labs(title=title) +
      theme(plot.title = element_text(size = 9)) +
      geom_hline(yintercept = 0, linetype=2, colour="grey50") +
      geom_vline(xintercept = 0, linetype=2, colour="grey50") +
      xlab(paste0("PC1 (", round(s$importance[2,1]*100, 1), "%)")) + 
      ylab(paste0("PC2 (", round(s$importance[2,2]*100, 1), "%)"))

    return(bplot)
}

