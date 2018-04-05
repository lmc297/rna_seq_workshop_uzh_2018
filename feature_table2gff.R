#!/usr/bin/env Rscript

# Rscript feature_table2gff.R <path to NCBI feature table>
# Convert NCBI feature table to 5-column GFF file
# RNA-Seq Workshop, University of Zurich, April 2018
# Laura Carroll, lmc297@cornell.edu

# Read command line arguments
args = commandArgs(trailingOnly=TRUE)

# Read feature table
feature.table <- read.delim(file = args[1],
                            header = T, sep = "\t")

# Save protein coding genes
gff.coding <- feature.table[which(feature.table$X..feature=="CDS"),]

# Create 5-column table
gff <- data.frame(gff.coding$locus_tag, gff.coding$genomic_accession, gff.coding$start, gff.coding$end, gff.coding$strand)

# Write table
outfile = gsub(pattern = "_feature_table.txt",
               replacement = "_final.gff", x = args[1])
write.table(x = gff, file = outfile,
            append = F, quote = F,
            sep = "\t", row.names = F, col.names = F)
