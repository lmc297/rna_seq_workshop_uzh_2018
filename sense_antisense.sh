#!/bin/bash
# sh sense_antisense.sh <path to directory with sorted bam files>
# a simple loop to separate sorted bam files into mapped sense and anti-sense reads and index these files
# example for RNA-Seq Workshop, University of Zurich, April 2018
# Laura Carroll, lmc297@cornell.edu

# go to the directory with our sorted bam files
cd $1
# loop through each sorted bam file
for f in *.sorted.bam
do
# separate out reads mapped to the anti-sense strand
samtools view -f 16 $f -b -o ${f%.sorted.bam}.as.sorted.bam
# separate out reads mapped to the sense strand
samtools view -F 16 $f -b -o ${f%.sorted.bam}.s.sorted.bam
# index sorted anti-sense strand bam file
samtools index ${f%.sorted.bam}.as.sorted.bam
# index sorted sense strand bam file
samtools index ${f%.sorted.bam}.s.sorted.bam
done
