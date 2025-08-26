
library(BaalChIP)
library(MASS)

sampleSheet <- 'sampleSheet.tsv'

hets <- c('essential'='gained_essential_SNPs.BaalChIPhets',
'nonEssential'='all.nonEssential_SNPs.BaalChIPhets')

##object
res <- BaalChIP(samplesheet=samplesheet, hets=hets)

##This function is a wrapper of the following functions: alleleCounts, QCfilter, mergePerGroup,filter1allele, getASB
res <- BaalChIP.run(res, cores=4)

summaryQC(res)
summaryASB(res)

counts <- BaalChIP.get(res, 'mergedCounts')

##essential
head(counts[[1]])
write.table(counts[[1]], "/path/to/your/output/gained_essential_SNPs.allelicCounts.matrix.txt", 
sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)


##nonEssential
head(counts[[2]])
write.table(counts[[2]], "/path/to/your/output/nonEssential_SNPs.allelicCounts.matrix.txt", 
sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)



##get ASB
ASB <- BaalChIP.report(res)
##essential
write.table(ASB[[1]], "/path/to/your/output/gained_essential_SNPs.mergedAllelicCounts.matrix.txt", 
sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
##nonEssential
write.table(ASB[[2]], "/path/to/your/output/nonEssential_SNPs.mergedAllelicCounts.matrix.txt", 
sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
