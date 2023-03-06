#!/bin/bash

module load lang/Java/13.0.2

mkdir processed_data

for i in ctrl sh70c
do
    for j in 1 2 3 4 5 6
    do
        sample=${i}_$j
        mkdir processed_data/$sample
	bbmap/bbsplit.sh build=1 threads=50 \
					   ref=references/hg38.fa,references/mm10.fa \
					   in=rawdata/${sample}_1.fq.gz  \
					   in2=rawdata/${sample}_2.fq.gz  \
					   scafstats=processed_data/$sample/$sample.scafstats.txt \
					   refstats=processed_data/$sample/$sample.refstats.txt \
 					   ambiguous2=toss \
					   basename=processed_data/$sample/$sample.%.R#.fq.gz
	done
done

