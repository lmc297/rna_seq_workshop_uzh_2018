#!/usr/bin/python

# python gene_list2protein_fasta.py <list of significant up- or down-regulated genes> <reference protein fasta faa file>
# Given a list of significant genes and reference protein fasta file, pull out the fasta sequences of significant genes
# RNA-Seq Workshop, University of Zurich, April 2018
# Laura Carroll, lmc297@cornell.edu

import sys
from Bio import SeqIO

# Read list of genes
significant_genes = []
gene_list = open(str(sys.argv[1]), "r")
for line in gene_list:
        if "annotation" not in line and "FDR.DE" not in line:
                protein = line.split("\t")[-1]
                significant_genes.append(protein.strip())
gene_list.close()

# Parse fasta and output protein sequence if protein sequence is present in significant gene list
infile = open(str(sys.argv[2]), "r")
prefix = str(sys.argv[1]).replace("_gene_list.txt", ".faa")
outfile = open(prefix, "a")
detected = []
for record in SeqIO.parse(infile, "fasta"):
        seqid = str(record.id).strip()
        seqdes = str(record.description).strip()
        seqseq = str(record.seq).strip()
        if seqid in significant_genes:
                detected.append(seqid)
                print >> outfile, ">"+seqdes
                print >> outfile, seqseq
infile.close()
outfile.close()
undetected = 0
for gene in significant_genes:
        if gene not in detected:
                undetected += 1
if undetected > 0:
        print "Unable to find protein accession numbers for "+str(undetected)+" gene(s). Ignoring "+str(undetected)+" gene(s)."
