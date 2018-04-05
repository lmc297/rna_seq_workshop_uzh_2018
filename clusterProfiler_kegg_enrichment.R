# clusterProfiler Workflow for Detecting Significant Gene Clusters
# RNA-Seq Workshop, University of Zurich, April 2018
# Laura Carroll, lmc297@cornell.edu

# Step 0: Set working directory
# setwd("/path/to/directory/with/feature_table.txt)

# Step 1: Install clusterProfiler, if necessary
#source("https://bioconductor.org/biocLite.R")
#biocLite("clusterProfiler")

# Step 2: Load clusterProfiler
library(clusterProfiler)

# Step 3: Load our feature table
feature.table <- read.delim(file = "GCF_000210855.2_ASM21085v2_feature_table.txt",
                            header = T, sep = "\t")

# Step 4: Format data for clusterProfiler

# Currently, clusterProfiler uses locus tags under "old_locus_tag" to use with KEGG
# We need to convert our current locus IDs for our up and down regulated genes to the old version

# Let's get the old IDs from our feature table
old <- feature.table[grepl(pattern = "old_locus_tag", x = feature.table$attributes),]
old.loci <- gsub(pattern = ".*=",
                 replacement = "",
                 x = old$attributes)

# This data frame has our current locus tags and the old one associated with it
# We can think of it as a "dictionary"
old.df <- data.frame(old$locus_tag, old.loci)

# Now let's load our upregulated genes
upregulated <- read.delim(file = "upregulated_gene_list.txt", 
                          header = T, sep = "\t")

# Let's find the rows of our dictionary that are in our upregulated gene list
old.up <- old.df[which(old.df$old.locus_tag%in%upregulated$annotation),]

# Now let's get our upregulated genes that are in our dictionary
# Some of our upregulated genes may not have an old locus ID, sadly
upregulated.subset <- upregulated[which(upregulated$annotation%in%old.up$old.locus_tag),]

# Let's check if our subsetted dictionary is in the same order as our upregulated gene list
table(droplevels(old.up$old.locus_tag)==droplevels(upregulated.subset$annotation))

# Let's reorder our subsetted dictionary 
old.upreordered <- old.up[ order(match(old.up$old.locus_tag, upregulated.subset$annotation)), ]

# NOW let's make sure our subsetted dictionary and upregulated gene list are in the same order
table(droplevels(old.upreordered$old.locus_tag)==droplevels(upregulated.subset$annotation))

# Our lists are in the same order!
# Now we can add the old locus column to our list of upregulated genes!
upregulated.subset$old_locus <- old.upreordered$old.loci

# Step 5: Find enriched (overrepresented) KEGG categories  

up.enriched <- enrichKEGG(gene = upregulated.subset$old_locus,
                 organism     = 'sey',
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "fdr")
View(up.enriched)

#  We can view our enriched pathways 

browseKEGG(up.enriched, 'sey02040')


# Step 6: Find enriched KEGG modules

up.module <- enrichMKEGG(gene = upregulated.subset$old_locus,
                   organism = 'sey', pvalueCutoff = 0.05, 
                   pAdjustMethod = "fdr")

View(up.module)


# Step 7: Repeat Steps 4-6 for Downregulated Genes

downregulated <- read.delim(file = "downregulated_gene_list.txt", 
                          header = T, sep = "\t")

old.down <- old.df[which(old.df$old.locus_tag%in%downregulated$annotation),]

downregulated.subset <- downregulated[which(downregulated$annotation%in%old.down$old.locus_tag),]

table(droplevels(old.down$old.locus_tag)==droplevels(downregulated.subset$annotation))

old.downreordered <- old.down[ order(match(old.down$old.locus_tag, downregulated.subset$annotation)), ]

table(droplevels(old.downreordered$old.locus_tag)==droplevels(downregulated.subset$annotation))

downregulated.subset$old_locus <- old.downreordered$old.loci

down.enriched <- enrichKEGG(gene = downregulated.subset$old_locus,
                          organism     = 'sey',
                          pvalueCutoff = 0.05, 
                          pAdjustMethod = "fdr")
View(down.enriched)

# We can also view this pathway
browseKEGG(up.enriched, 'sey01130')


down.module <- enrichMKEGG(gene = downregulated.subset$old_locus,
                         organism = 'sey', pvalueCutoff = 0.05, 
                         pAdjustMethod = "fdr")

View(down.module)



