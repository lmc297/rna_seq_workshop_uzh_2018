#!/bin/bash

# sh sam2sortedbam.sh <path to directory with .sam files>
# a simple loop to convert sam to sorted bam files
# example for RNA-Seq Workshop, University of Zurich, April 2018
# Laura Carroll, lmc297@cornell.edu

# go to the directory with our *.sam files
cd $1

# loop through each sam file
for f in *.sam
do
# convert sam file to bam file, skipping alignments with mapping quality smaller than 1
samtools view -b -q 1 $f -o ${f%.sam}.bam

# sort the resulting bam file
samtools sort ${f%.sam}.bam -o ${f%.sam}.sorted.bam

done
