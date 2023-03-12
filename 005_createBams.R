#!/usr/bin/env Rscript
library(Rsubread)

#ctrl_1  ctrl_2  ctrl_3  ctrl_4  ctrl_5  ctrl_6  sh70c_1  sh70c_2  sh70c_3  sh70c_4  sh70c_5  sh70c_6
#sh70c_6.hg38.R1.fq.gz

for (condition in c("ctrl","sh70c"))
  {
  for (rep in c(1:6))
    {
    directory=paste0('processed_data','/',condition,'_',rep)
    filenameR1=paste0(condition,'_',rep,'.hg38.R', 1,'.fq.gz')
    filenameR2=paste0(condition,'_',rep,'.hg38.R', 2,'.fq.gz')
    sampleR1=paste0(directory,'/',filenameR1)
    sampleR2=paste0(directory,'/',filenameR2)
    BAM=paste0('processed_data/',condition,'_',rep,'.hg38.bam')
    align(index="RsubreadIndexHg",readfile1=sampleR1,readfile2=sampleR2,type="rna", sortReadsByCoordinates=TRUE, output_file=BAM)
    print(sampleR1)
    print(sampleR2)
    print(BAM)
    }
  }


