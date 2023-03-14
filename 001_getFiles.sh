#!/bin/bash

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
