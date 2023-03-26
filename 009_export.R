library(org.Hs.eg.db)

#Load Data
res<-readRDS("007_res.RDS")

res_df<-as.data.frame(res)

# Add gene names and row names (entrez IDs)
res_df$gene_name <- mapIds(org.Hs.eg.db, keys=row.names(res_df), column="SYMBOL", keytype="ENTREZID", multiVals="first")
res_df$entrez_id <- row.names(res_df)
write.csv(res_df, file="processed_data/009_DeSeq_Results.csv")


countMatrix<-as.data.frame(readRDS(file='processed_data/countMatrix.RDS')$counts)
countMatrix
countMatrix$gene_name <- mapIds(org.Hs.eg.db, keys=row.names(countMatrix), column="SYMBOL", keytype="ENTREZID", multiVals="first")
countMatrix$entrez_id <- row.names(countMatrix)
write.csv(res_df, file="processed_data/009_CountMatrix.csv")