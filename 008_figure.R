

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

# 
# #Some quick plots of the data
# plotMA(res, ylim=c(-5,5))
# plotCounts(dds, gene='6331', intgroup="condition")
# 

res<-readRDS("007_res.RDS")

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
