#!/bin/bash
FILE=csv
if [ -d "$FILE" ]; then
    echo "$FILE exists."
else 
    echo "$FILE does not exist. Creating"
    mkdir $FILE
fi


#Differential Expression analysis from Novagene 
cp ../20230222-tess-leslie/X204SC20090619-Z01-F001/X204SC20090619-Z01-F001/Result_X204SC20090619-Z01-F001_hg38/4.DiffExprAnalysis/DEGlist/sh70cvsctrl/sh70cvsctrl.DEG.xls csv/sh70cvsctrl.DEG.csv
cp ../20230222-tess-leslie/X204SC20090619-Z01-F001/X204SC20090619-Z01-F001/Result_X204SC20090619-Z01-F001_hg38/3.Quantification/Count/readcount_genename.xls csv/readcount_genename.csv 

FILE=txt
if [ -d "$FILE" ]; then
    echo "$FILE exists."
else
    echo "$FILE does not exist. Creating"
    mkdir $FILE
fi


#Ouput? But from where?
cp ../20230222-tess-leslie/ambiguous.genes.txt txt/ambiguous.genes_OldPipeline.txt


#Output off BBSplit on the .FA files from Novagene
cp ../disambiguate_reads/ambiguous_read_counts.tsv txt/ambiguous_read_counts_OldPipeline.tsv



FILE=scripts
if [ -d "$FILE" ]; then
    echo "$FILE exists."
else 
    echo "$FILE does not exist. Creating"
    mkdir $FILE
fi

#Get Script for BBSplit 
cp ../disambiguate_reads/run_bbsplit.sh scripts

#Get Raw Data from Novagene
FILE=rawdata
if [ -d "$FILE" ]; then
    echo "$FILE exists."
else
    echo "$FILE does not exist. Creating"
    mkdir $FILE
fi
cp  ../20230222-tess-leslie/X204SC20090619-Z01-F001/X204SC20090619-Z01-F001/rawdata/*.fq.gz rawdata



FILE=references
if [ -d "$FILE" ]; then
    echo "$FILE exists."
else
    echo "$FILE does not exist. Creating"
    mkdir $FILE
fi

 cp ../disambiguate_reads/mm10.fa references
 cp ../disambiguate_reads/hg38.fa references
