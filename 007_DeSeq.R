#!/usr/bin/env Rscript

countMatrix<-readRDS(file='processed_data/countMatrix.RDS')$counts
countMatrix

cage=c("H228", "H228", "H210", "H178", "H210", "H240", 
       "H288", "H288", "H240", "H178", "H210", "H210")

coldata<-cbind(rep(1:6,2), c(rep('ctrl',6),rep('sh70c',6)),cage)
rownames(coldata)<-colnames(countMatrix)
colnames(coldata)<-c('rep','condition','cage')

#### DeSeq2 Analysis
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
saveRDS(res,file="007_res.RDS")


#### Suplimentary Plot
library(ggpubr)

#Check SCN5a counts
df<-data.frame(Counts=countMatrix['6331',],
              Condition=c(rep('Control',6),
              rep('shRNA SCN5a',6)),
              Replicate=rep(1:6,2))
my_comparisons<-list(c('Control','shRNA SCN5a'))
SCN5aCounts<-ggboxplot(df, x='Condition', 
                           y='Counts',
                           ylim=c(0,1750),
                           color='Condition',
                          add="jitter",
                           legend.title="Raw counts") +
                          ylab("SCN5a Counts") +
            stat_compare_means(comparisons = my_comparisons) + 
            theme(legend.text=element_text(size=8),
                  legend.title=element_text(size=10))


#Check PCA corrected for cage.
vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)
mm <- model.matrix(~condition, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$cage, design=mm)
assay(vsd) <- mat
cagePlot<-plotPCA(vsd, intgroup="cage")+
                              ggtitle("PCA of samples coloured by cage")+
                              guides(color=guide_legend(title="Cage"))
condPlot<-plotPCA(vsd, intgroup="condition")+
                              ggtitle("PCA of sample coloured by condition")+
                              guides(color=guide_legend(title="Condition"))

#GSEA SCHUETZ Ductal Down 
library(vulcan)
library(ZMIZ1) #from Andrew Holding's github
geneList<-res$log2FoldChange
names(geneList) <- rownames(res)
SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_DN <- msigdb[["c2_cgp;_;SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_DN"]]
obj <- gsea(sort(geneList, decreasing = TRUE), set = SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_DN, 
            method = "pareto")

gseaPlotDN<-""

pdf("figures/007_suppA-C.pdf",width=8,heigh= 5)
ggarrange(SCN5aCounts, condPlot, cagePlot,  gseaPlotDN,
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)
dev.off()

pdf("figures/007_suppD.pdf",width=8,heigh= 4)
plot_gsea(obj, bottomYtitle = "sh70c/ctrl", title = "SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_DN")
dev.off()
