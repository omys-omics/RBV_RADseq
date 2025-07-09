setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/RBVs/09.5plates.phylogeography/mito.SNPs/all")

library(vcfR)
library(SNPfiltR)
library(readr)


############################################################
####### Visualize and filter SNPs by quality and missingness
############################################################


# Read in vcf
all_species <- read.vcfR("populations.snps.vcf")

#generate popmap file. Two column popmap with the same format as stacks, and the columns must be named 'id' and 'pop'
detailed.popmap <- read_delim("~/Documents/Ben's Stuff/0 KU/Dissertation/RBVs/detailed.popmap.20240229.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
mito.names <- colnames(extract.gt(all_species))
popmap <- data.frame(detailed.popmap[detailed.popmap$id %in% colnames(extract.gt(all_species)),])
popmap <- popmap[,c("id", "pop")]

# make sure popmap is in same order as ids in the vcf
identical(colnames(extract.gt(all_species)), popmap$id)
popmap <- popmap[order(match(popmap$id, colnames(extract.gt(all_species)))),]
identical(colnames(extract.gt(all_species)), popmap$id)


# get rid of samples with all missing genotypes
all_species_genotypes <- extract.gt(all_species)
n <- 0
for (x in colnames(all_species_genotypes)) {
  n <- n + 1
  if (all(is.na(all_species_genotypes[,x]))) {
    print(x)
    print(n)
  }
}


#plot depth per snp and per sample
dp <- extract.gt(all_species, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)

#plot genotype quality per snp and per sample
gq <- extract.gt(all_species, element = "GQ", as.numeric=TRUE)
heatmap.bp(gq, rlabels = FALSE)

# average depth per sample
depth <- dp
depth[is.na(depth)] <- 0
mean.depth <- colSums(depth)/nrow(depth)
hist(mean.depth, breaks = 113)
sort(mean.depth)

#visualize distributions
hard_filter(vcfR=all_species)

# hard filter to minimum depth of 5, minimum genotype quality of 30
all_species_filter <- hard_filter(vcfR=all_species, depth = 4, gq = 30)

# allele balance filter (by setting min to 0.6 and max to 0.4, all hets are removed)
all_species_filter <- filter_allele_balance(all_species_filter, min.ratio = 0.6, max.ratio = 0.4)

#remove invariant sites
all_species_filter<-min_mac(all_species_filter, min.mac = 1)

#visualize missing data by SNP and the effect of various cutoffs on the missingness of each sample
missing_by_snp(all_species_filter)

#choose a value that retains an acceptable amount of missing data in each sample, and maximizes SNPs retained while minimizing overall missing data, and filter vcf
all_species_filter<-missing_by_snp(all_species_filter, cutoff = 0.1)

#subset popmap to only include retained individuals
popmap_filter<-popmap[popmap$id %in% colnames(all_species_filter@gt),]

#see which individuals and pops were dropped
setdiff(popmap$id, popmap_filter$id)
setdiff(popmap$pop, popmap_filter$pop)

#remove invariant sites
all_species_filter<-min_mac(all_species_filter, min.mac = 1)

#run function to visualize samples above the threshold we want from the vcf
missing_by_sample(vcfR=all_species_filter, popmap = popmap_filter)

#plot depth per snp and per sample
dp <- extract.gt(all_species_filter, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)

#plot genotype quality per snp and per sample
gq <- extract.gt(all_species_filter, element = "GQ", as.numeric=TRUE)
heatmap.bp(gq, rlabels = FALSE)

# Loci kept:
all_species_filter@fix
snps <- all_species_filter@fix[,"ID"]
loci <- unique(sapply(strsplit(snps, ":"), '[', 1))
loci
length(loci)

#write out thinned vcf
vcfR::write.vcf(all_species_filter, "mito.SNPs.all.vcf.gz")
