setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/RBVs/09.5plates.phylogeography/mito.SNPs/both")

library(vcfR)
library(SNPfiltR)
library(readr)


############################################################
####### Visualize and filter SNPs by quality and missingness
############################################################


# Read in vcf
all_species <- read.vcfR("populations.snps.vcf")

#generate popmap file. Two column popmap with the same format as stacks, and the columns must be named 'id' and 'pop'
detailed.popmap <- read_delim("../all/popmap.mt_species_3.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
mito.names <- colnames(extract.gt(all_species))
popmap <- data.frame(detailed.popmap[detailed.popmap$id %in% colnames(extract.gt(all_species)),])
popmap <- popmap[,c("id", "mt_species")]
colnames(popmap) <- c("id", "pop")

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

# hard filter to minimum depth of 4, minimum genotype quality of 30
all_species_filter <- hard_filter(vcfR=all_species, depth = 4, gq = 30)

# allele balance filter (by setting min to 0.6 and max to 0.4, all hets are removed)
all_species_filter <- filter_allele_balance(all_species_filter, min.ratio = 0.6, max.ratio = 0.4)

#remove invariant sites
all_species_filter<-min_mac(all_species_filter, min.mac = 1)







#######################
####    RUTILUS    ####
#######################
rut.pop <- popmap[popmap$pop=="rutilus",]
rutilus <- all_species_filter[samples=rut.pop$id]


#visualize missing data by SNP and the effect of various cutoffs on the missingness of each sample
missing_by_snp(rutilus)

#choose a value that retains an acceptable amount of missing data in each sample, and maximizes SNPs retained while minimizing overall missing data, and filter vcf
rutilus<-missing_by_snp(rutilus, cutoff = 0.3)

#remove invariant sites
rutilus<-min_mac(rutilus, min.mac = 1)

#run function to visualize samples above the threshold we want from the vcf
missing_by_sample(vcfR=rutilus, popmap = rut.pop)

#plot depth per snp and per sample
dp <- extract.gt(rutilus, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)

#plot genotype quality per snp and per sample
gq <- extract.gt(rutilus, element = "GQ", as.numeric=TRUE)
heatmap.bp(gq, rlabels = FALSE)

# Loci kept:
rutilus@fix
rut.snps <- rutilus@fix[,"ID"]
rut.loci <- unique(sapply(strsplit(rut.snps, ":"), '[', 1))
rut.loci
length(rut.loci)


# Identify rutilus with <= 20% missing data on rutilus-specific SNPs

# make sure popmap is in same order as ids in the vcf
identical(colnames(extract.gt(rutilus)), rut.pop$id)

# calculate percent missing data for each individual
percent_missing <- colSums(is.na(rutilus@gt))/nrow(rutilus@gt)
identical(names(percent_missing[-1]), rut.pop$id)
rut.pop$percent_missing <- percent_missing[-1]

rut.pop <- rut.pop[rut.pop$percent_missing<=0.2,]




#######################
####    GAPPERI    ####
#######################
gap.pop <- popmap[popmap$pop=="gapperi",]
gapperi <- all_species_filter[samples=gap.pop$id]


#visualize missing data by SNP and the effect of various cutoffs on the missingness of each sample
missing_by_snp(gapperi)

#choose a value that retains an acceptable amount of missing data in each sample, and maximizes SNPs retained while minimizing overall missing data, and filter vcf
gapperi<-missing_by_snp(gapperi, cutoff = 0.3)

#remove invariant sites
gapperi<-min_mac(gapperi, min.mac = 1)

#run function to visualize samples above the threshold we want from the vcf
missing_by_sample(vcfR=gapperi, popmap = gap.pop)

#plot depth per snp and per sample
dp <- extract.gt(gapperi, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)

#plot genotype quality per snp and per sample
gq <- extract.gt(gapperi, element = "GQ", as.numeric=TRUE)
heatmap.bp(gq, rlabels = FALSE)

# Loci kept:
gapperi@fix
gap.snps <- gapperi@fix[,"ID"]
gap.loci <- unique(sapply(strsplit(gap.snps, ":"), '[', 1))
gap.loci
length(gap.loci)






############################
####        BOTH        ####
############################

# all unique snps in rut and gap
snps <- c(rut.snps, gap.snps)
snps <- unique(snps)

# all snps from lightly filtered vcf
all.snps <- all_species_filter@fix[,"ID"]
snp.indices <- which(all.snps %in% snps)

# make new vcf of all inds, but only at the SNPs identified for rutilus and gapperi
both <- all_species_filter[snp.indices]


# Only keep rutilus with low missing data
both <- both[samples=c(rut.pop$id, gap.pop$id, "UAM23201", "UAM34216")]
both <- min_mac(both, min.mac = 1)
both
popmap <- popmap[popmap$id %in% c(rut.pop$id, gap.pop$id, "UAM23201", "UAM34216"),]

#run function to visualize samples above the threshold we want from the vcf
missing_by_sample(vcfR=both, popmap = popmap)

#plot depth per snp and per sample
dp <- extract.gt(both, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)

#plot genotype quality per snp and per sample
gq <- extract.gt(both, element = "GQ", as.numeric=TRUE)
heatmap.bp(gq, rlabels = FALSE)

# Loci kept:
both@fix
both.snps <- both@fix[,"ID"]
both.loci <- unique(sapply(strsplit(both.snps, ":"), '[', 1))
both.loci
length(both.loci)


#write out thinned vcf
vcfR::write.vcf(both, "mito.SNPs.both.vcf.gz")
