#!/usr/bin/env bash
for i in ctrl sh70c
do
    for j in 1 2 3 4 5 6
    do
        sample=${i}_$j
        mkdir $sample
        bbsplit.sh build=1 threads=50 in=../X204SC20090619-Z01-F001/X204SC20090619-Z01-F001/rawdata/${sample}_1.fq.gz in2=../X204SC20090619-Z01-F001/X204SC20090619-Z01-F001/rawdata/${sample}_2.fq.gz scafstats=$sample/$sample.scafstats.txt refstats=$sample/$sample.refstats.txt ambiguous2=split basename=$sample/$sample.%.R#.fq.gz &> $sample/$sample.bbsplit.log
    done
done
