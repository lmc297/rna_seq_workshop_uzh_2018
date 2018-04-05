library(Rsubread)

setwd("/Volumes/My Passport for Mac/rnaseq2018/typhimurium/bam/")
bamlist <- Sys.glob(as.character("*s.sorted.bam"))

ann <- read.delim(file = "/Volumes/My Passport for Mac/rnaseq2018/typhimurium/GCF_000210855.2_ASM21085v2_final.gff",
                  header = F, sep = "\t", stringsAsFactors = F)
colnames(ann) <- c("GeneID", "Chr", "Start", "End", "Strand")


fc <- featureCounts(files = bamlist, annot.ext = ann, isPairedEnd = T)

fc.counts <- as.data.frame(fc$counts)

table(rownames(fc.counts)==ann$GeneID)

fc.counts$Start <- ann$Start
fc.counts$End <- ann$End
fc.counts$Strand <- ann$Strand
fc.counts$geneID <- rownames(fc.counts)

write.table(x = fc.counts, file = "final_coverage_table.txt", 
            append = F, quote = F, row.names = F, col.names = T, sep = "\t")






