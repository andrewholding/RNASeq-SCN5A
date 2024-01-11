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
res_df$gene_name <- mapIds(org.Hs.eg.db, keys=row.names(res_df),
                           column="SYMBOL", keytype="ENTREZID",
                           multiVals="first")
res_df$entrez_id <- row.names(res_df)
 
library(clusterProfiler)
library(ggplot2)


plotGO<-function(log2foldChange=0, padj=1,terms="BP") {
  sigResults<-res_df[res_df$log2FoldChange < log2foldChange & res_df$padj < padj,]$entrez_id
  print(paste0("Genes: ", length(  sigResults)))
  goResults<-enrichGO(gene=sigResults, OrgDb="org.Hs.eg.db", keyType="ENTREZID", ont=terms) 
  print(paste0("Go Terms: ", nrow(as.data.frame(goResults))))
  if(terms=="BP") {title="Biological Process"} 
  if(terms=="MF") {title="Molecular Function"} 
  if(terms=="CC") {title="Cellular Compartment"} 
  return(dotplot(goResults, showCategory = 20)   
      + ggtitle(paste0(title,", Log2 fold change < ", log2foldChange,", Padj < ", padj ))
        )   
}

pdf(file="GO terms", width=8, height=11)
for (terms in c("BP","MF","CC")) {
  print(terms)
  p<-plotGO(log2foldChange=0,padj=0.05,terms=terms) #note log2 is <
  plot(p)
}
dev.off()
