#Useful MSigDB (From ZMIZ1 package)
library(org.Hs.eg.db)
library(vulcan)
library(ZMIZ1)

#Load Data
res<-readRDS("007_res.RDS")
# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Convert the DESeq2 results object to a data frame
res_df <- as.data.frame(res)

# Add gene names and row names (entrez IDs)
res_df$gene_name <- mapIds(org.Hs.eg.db, keys=row.names(res_df), column="SYMBOL", keytype="ENTREZID", multiVals="first")
res_df$entrez_id <- row.names(res_df)

#Load MSigDB genes
SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP <- msigdb[["c2_cgp;_;SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP"]]

# Define the MSigDB gene list "SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP"
msigdb_genes <- c(SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP) # Replace with the actual gene symbols

# Mark the genes from the MSigDB gene list in a new column
res_df$in_msigdb <- res_df$entrez_id %in% msigdb_genes

# Add a column to identify SCN5A
res_df$SCN5A <- res_df$gene_name == "SCN5A"

# Sort the data frame by adjusted p-value and select the top 10 genes
top_genes <- res_df %>% arrange( padj ) %>% head(15)

# Create the volcano plot
volcano_plot <- ggplot() +
  geom_point(data=res_df[!res_df$in_msigdb & !res_df$SCN5A,], aes(x=log2FoldChange, y=-log10(padj), color="Others"), alpha=0.6) +
  geom_point(data=res_df[res_df$in_msigdb,], aes(x=log2FoldChange, y=-log10(padj), color="MSigDB"), alpha=0.6) +
  geom_point(data=res_df[res_df$SCN5A,], aes(x=log2FoldChange, y=-log10(padj), color="SCN5A"), alpha=0.6) +
  geom_text_repel(show.legend = FALSE, data=top_genes, aes(x=log2FoldChange, y=-log10(padj), label=gene_name, color=factor(ifelse(in_msigdb, "MSigDB", ifelse(SCN5A, "SCN5A", "OthersText")))), size=4, box.padding=0.2) +
  scale_color_manual(values=c("Others" = "gray", "MSigDB" = "blue", "SCN5A" = "red", "OthersText" = "black"), labels = c("MSigDB Gene Set", "SCN5A"), breaks = c("MSigDB", "SCN5A")) +
  labs(title="RNA-seq expression changes in response to SCN5a shRNA knockdown ",
       x="log2 fold change",
       y="-log10 adjusted p-value") +
  theme_minimal() +
  guides(color = guide_legend(title = "Key"))

# Show the plot
volcano_plot

#
pval_to_asterisk_R <- function(pval) {
  if (pval > 0.05) {
    return(" ")
  } else if (pval <= 0.05 & pval > 0.01) {
    return("*   ")
  } else if (pval <= 0.01 & pval > 0.001) {
    return("**  ")
  } else if (pval <= 0.001 & pval > 0.0001) {
    return("*** ")
  } else {
    return("****")
  }
}

pval_to_asterisk_L <- function(pval) {
  if (pval > 0.05) {
    return(" ")
  } else if (pval <= 0.05 & pval > 0.01) {
    return("   *")
  } else if (pval <= 0.01 & pval > 0.001) {
    return("  **")
  } else if (pval <= 0.001 & pval > 0.0001) {
    return(" ***")
  } else {
    return("****")
  }
}

# Add a new column with the asterisk notation for p-values in res_df
res_df <- res_df %>%
  mutate(pval_asterisk_L = sapply(padj, pval_to_asterisk_L))

res_df <- res_df %>%
  mutate(pval_asterisk_R = sapply(padj, pval_to_asterisk_R))

#Create top 20 gene plots
top_genes_up <- res_df %>% arrange(desc(log2FoldChange)) %>% head(20)
top_genes_down <- res_df %>% arrange(log2FoldChange) %>% head(20)
top20_genes <- rbind(top_genes_up, top_genes_down)

# Create a ggplot object for genes with increased fold change
increased_fc_plot <- ggplot(data = top_genes_up, aes(x = log2FoldChange, y = reorder(gene_name, log2FoldChange), fill = factor(ifelse(in_msigdb, "MSigDB", ifelse(SCN5A, "SCN5A", "Others"))))) +
  geom_col() +
  scale_fill_manual(values = c("Others" = "black", "MSigDB" = "blue", "SCN5A" = "red")) +
  geom_text(data=top_genes_up, aes(x = log2FoldChange, y = reorder(gene_name, log2FoldChange), label = pval_asterisk_L, vjust = 0.8, hjust=1.3,size=10, color="white")) +
  labs(title = "Increased expression (Top 20)",
       x = "log2 fold change",
       y = "Gene name") +
  scale_color_identity() +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(hjust = -0.2))
 

# Create a ggplot object for genes with decreased fold change
decreased_fc_plot <- ggplot(data = top_genes_down, aes(x = log2FoldChange, y = reorder(gene_name, log2FoldChange), fill = factor(ifelse(in_msigdb, "MSigDB", ifelse(SCN5A, "SCN5A", "Others"))))) +
  geom_col() +
  scale_fill_manual(values = c("Others" = "black", "MSigDB" = "blue", "SCN5A" = "red")) +
  geom_text(data=top_genes_down, aes(x = log2FoldChange, y = reorder(gene_name, log2FoldChange), label = pval_asterisk_R, vjust = 0.8, hjust=-0.5,size=10, color="white")) +
  labs(title = "Decreased expression (Top 20)",
       x = "log2 fold change",
       y = "Gene name") +
  scale_color_identity() +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(hjust = -0.2))

# Show the plots
increased_fc_plot
decreased_fc_plot


library(cowplot)
# Combine the three plots using plot_grid() and label them A, B, and C
combined_plot <- plot_grid(
  plot_grid(volcano_plot, NULL, ncol = 1, rel_heights = c(1, 0.05), labels = c("A", "")),
  plot_grid(increased_fc_plot + theme(plot.title = element_text(hjust = -2.5)), 
            decreased_fc_plot + theme(plot.title = element_text(hjust = -2.3)), ncol = 2, labels = c("B", "C")),
  ncol = 1, rel_heights = c(1, 1)
)

# Show the combined plot
combined_plot


#Plots GSEA
geneList<-res$log2FoldChange
names(geneList) <- rownames(res)

SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP <- msigdb[["c2_cgp;_;SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP"]]
obj <- gsea(sort(geneList, decreasing = TRUE), set = SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP, 
            method = "pareto")
plot_gsea(obj, bottomYtitle = "sh70c/ctrl", title = "SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP")



pdf("figures/008_FigA-C.pdf",width=8,heigh= 6)
combined_plot
dev.off()

pdf("figures/008_FigD.pdf",width=8,heigh= 4)
plot_gsea(obj, bottomYtitle = "sh70c/ctrl", title = "SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP")
dev.off()


