#!/usr/bin/env Rscript



countMatrix<-readRDS(file='processed_data/countMatrix.RDS')$counts
countMatrix



coldata<-cbind(rep(1:6,2), c(rep('ctrl',6),rep('sh70c',6)))
rownames(coldata)<-colnames(countMatrix)
colnames(coldata)<-c('rep','condition')


#Check SCN5a counts
df<-data.frame(Counts=countMatrix['6331',],Condition=c(rep('Control',6),rep('shRNA SCN5a',6)),Replicate=rep(1:6,2))
library(ggpubr)           
my_comparisons<-list(c('Control','shRNA SCN5a'))
ggboxplot(df, x='Condition', y='Counts',color='Condition', add="jitter", legend.title="Raw SCN5a RNA-seq counts") + stat_compare_means(comparisons = my_comparisons)


library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = coldata,
                              design = ~ condition)

#Ensure that transcripts have at least 3 samples with a count of 10 or more.
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

#Ensure control is Reference
dds$condition <- relevel(dds$condition, ref = "ctrl")


#Run!
dds <- DESeq(dds)
res <- results(dds)
res



#Useful function to give gene symbols (From ZMIZ1 package)
library(org.Hs.eg.db)
eg2sym<-function(x){
    eg2symmap<-as.list( org.Hs.egSYMBOL[mappedkeys( org.Hs.egSYMBOL)])
  x<-as.character(x)
  out<-eg2symmap[x]
  names(out)<-x
  out2<-sapply(out,function(x){
    if(is.null(x)){
      return(NA)
    } else {
      return(x[1])
    }
  })
  out3<-unlist(out2)
  out<-setNames(out3,names(out))
  return(out)
}

#Add gene symbols to results
res$symbol<-eg2sym(rownames(res))


#Some quick plots of the data
plotMA(res, ylim=c(-5,5))
plotCounts(dds, gene='6331', intgroup="condition")


#PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd)
plotPCA(vsd,intgroup="rep")
