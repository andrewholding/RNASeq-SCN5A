#!/usr/bin/env Rscript



countMatrix<-readRDS(file='processed_data/countMatrix.RDS')$counts
countMatrix


cage=c("H228", "H228", "H210", "H178", "H210", "H240", 
       "H288", "H288", "H240", "H178", "H210", "H210")

coldata<-cbind(rep(1:6,2), c(rep('ctrl',6),rep('sh70c',6)),cage)
rownames(coldata)<-colnames(countMatrix)
colnames(coldata)<-c('rep','condition','cage')



#Check SCN5a counts
df<-data.frame(Counts=countMatrix['6331',],Condition=c(rep('Control',6),rep('shRNA SCN5a',6)),Replicate=rep(1:6,2))
library(ggpubr)           
my_comparisons<-list(c('Control','shRNA SCN5a'))
ggboxplot(df, x='Condition', y='Counts',color='Condition', add="jitter", legend.title="Raw SCN5a RNA-seq counts") + stat_compare_means(comparisons = my_comparisons)


library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = coldata,
                              design = ~ cage + condition)

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
plotPCA(vsd,intgroup="cage")+ggtitle("Uncorrected for Cage Effect")


mat <- assay(vsd)
mm <- model.matrix(~condition, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$cage, design=mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup="cage")+ggtitle("Correct for Cage Effect")


library(vulcan)
library(ZMIZ1)
geneList<-res$log2FoldChange
names(geneList) <- rownames(res)

SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP <- msigdb[["c2_cgp;_;SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP"]]
obj <- gsea(sort(geneList, decreasing = TRUE), set = SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP, 
            method = "pareto")
plot_gsea(obj, bottomYtitle = "sh70c/ctrl", title = "SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP")


SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_DN <- msigdb[["c2_cgp;_;SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_DN"]]
obj <- gsea(sort(geneList, decreasing = TRUE), set = SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_DN, 
            method = "pareto")
plot_gsea(obj, bottomYtitle = "sh70c/ctrl", title = "SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_DN")


VANTVEER_BREAST_CANCER_METASTASIS_DN <- msigdb[["c2_cgp;_;WANG_TUMOR_INVASIVENESS_UP"]]
obj <- gsea(sort(geneList, decreasing = TRUE), set = VANTVEER_BREAST_CANCER_METASTASIS_DN, 
            method = "pareto")
plot_gsea(obj, bottomYtitle = "sh70c/ctrl", title = "VANTVEER_BREAST_CANCER_METASTASIS_DN")


GOBP_ACID_SECRETION<-c(1066,1133,115111,123264,18,1813,1814,189,200931,2550,29887,320,3274,3351,3784,389015,407039,4763,4923,4987,5020,5021,5027,51738,5244,54659,5562,5733,6343,6446,65010,6529,6554,7032,7124,7200,7223,7349,7369,7442,785,846,85413,8647,887,9026,9368,9900)
GOBP_LACTATE_METABOLIC_PROCESS<-c(11315,197257,25953,3029,3091,347862,3939,3945,3948,406992,5208,57103,7157,8604,8864,89,92483)
GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT<-c(133418,159963,160728,23539,6566,9123,9194)

obj <- gsea(sort(geneList, decreasing =TRUE), set = GOBP_ACID_SECRETION, 
            method = "pareto")
plot_gsea(obj, bottomYtitle = "sh70c/ctrl", title = "GOBP_ACID_SECRETION")

obj <- gsea(sort(geneList, decreasing = TRUE), set = GOBP_LACTATE_METABOLIC_PROCESS, 
            method = "pareto")
plot_gsea(obj, bottomYtitle = "sh70c/ctrl", title = "GOBP_LACTATE_METABOLIC_PROCESS")

obj <- gsea(sort(geneList, decreasing = TRUE), set = GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT, 
            method = "pareto")
plot_gsea(obj, bottomYtitle = "sh70c/ctrl", title = "GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT")



