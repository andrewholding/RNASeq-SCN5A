#!/usr/bin/env Rscript
library(Rsubread)


condition<-c(rep("ctrl",6),rep("sh70c",6))
rep<-c(1:6,1:6)
BAMFiles<-paste0('processed_data/',condition,'_',rep,'.hg38.bam')

print(BAMFiles)


countMatrix<-featureCounts(files=BAMFiles,annot.inbuilt="hg38", nthreads=10, isPairedEnd=TRUE)


saveRDS(file='processed_data/countMatrix.RDS',countMatrix)

