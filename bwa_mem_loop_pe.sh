#!/bin/bash


# sh bwa_mem_loop_pe.sh <path to reference.fasta> <path to directory with paired-end .fastq.gz reads>
# a simple loop to create bwa indexes and run bwa mem for paired-end reads ending in _1.fastq.gz and _2.fastq.gz
# example for RNA-Seq Workshop, University of Zurich, April 2018
# Laura Carroll, lmc297@cornell.edu

# first, let's index our reference
bwa index $1

# now, let's move to our directory with .fastq.gz files
cd $2
# let's loop through all of our forward reads
for f in *_1.fastq.gz
do
# run bwa mem using both forward and reverse reads 
# save our output as a sam file with a name that matches our fastq.gz file
bwa mem $1 $f ${f%_1.fastq.gz}_2.fastq.gz > ${f%_1.fastq.gz}.sam
done
