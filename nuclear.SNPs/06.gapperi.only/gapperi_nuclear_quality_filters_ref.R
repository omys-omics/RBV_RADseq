setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/RBVs/09.5plates.phylogeography/nuclear.SNPs/06.gapperi.only")

library(vcfR)
library(SNPfiltR)
library(readr)


############################################################
####### Visualize and filter SNPs by quality and missingness
############################################################


# Read in vcf
all_species <- read.vcfR("../01.quality.filters/all.nuc.recode.vcf.gz")

# read in popmap with nuc and mito calls
detailed.popmap <- read_delim("../02.nuc.analyses/no.russ.popmap.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

# keep only gapperi
detailed.popmap <- detailed.popmap[detailed.popmap$nuc_species=="gapperi",]


# make popmap
popmap <- detailed.popmap[,c("id", "pop")]

# filter vcfR for correct gapperi individuals
gapperi <- all_species[samples=popmap$id]


# make sure popmap is in same order as ids in the vcf
identical(colnames(extract.gt(gapperi)), popmap$id)
popmap <- popmap[order(match(popmap$id, colnames(extract.gt(gapperi)))),]
identical(colnames(extract.gt(gapperi)), popmap$id)


#plot depth per snp and per sample
dp <- extract.gt(gapperi, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)

#plot genotype quality per snp and per sample
gq <- extract.gt(gapperi, element = "GQ", as.numeric=TRUE)
heatmap.bp(gq, rlabels = FALSE)

# average depth per sample
depth <- dp
depth[is.na(depth)] <- 0
mean.depth <- colSums(depth)/nrow(depth)
hist(mean.depth, breaks = 265)

#visualize distributions
hard_filter(vcfR=gapperi)

# hard filter to minimum depth of 4, minimum genotype quality of 30
gapperi_filter <- hard_filter(vcfR=gapperi, depth = 4, gq = 30)

# allele balance filter
filter_allele_balance(gapperi_filter)

#visualize and pick appropriate max depth cutoff
max_depth(gapperi_filter)
gapperi_filter <- max_depth(gapperi_filter, maxdepth = 400)

#remove invariant sites
gapperi_filter<-min_mac(gapperi_filter, min.mac = 1)


#visualize missing data by SNP and the effect of various cutoffs on the missingness of each sample
missing_by_snp(gapperi_filter)

#verify that missing data is not driving clustering patterns among samples
miss<-assess_missing_data_pca(vcfR=gapperi_filter, popmap = popmap, thresholds = 0.8, clustering = FALSE)


#choose a value that retains an acceptable amount of missing data in each sample, and maximizes SNPs retained while minimizing overall missing data, and filter vcf
gapperi_filter_cutoff<-missing_by_snp(gapperi_filter, cutoff = 0.8)
popmap_filter<-popmap[popmap$id %in% colnames(gapperi_filter@gt),]


#plot depth per snp and per sample
dp <- extract.gt(gapperi_filter_cutoff, element = "DP", as.numeric=TRUE)
depth <- dp
depth[is.na(depth)] <- 0
mean.depth <- colSums(depth)/nrow(depth)
hist(mean.depth, breaks = 250)
hist(mean.depth, breaks = 250, xlim = c(0,200))
sort(mean.depth)
summary(mean.depth)
heatmap.bp(dp, rlabels = FALSE)

#missing data for each sample
sort(colMeans(is.na(dp)))
colMeans(is.na(dp))[colMeans(is.na(dp))>0.5]


#plot genotype quality per snp and per sample
gq <- extract.gt(gapperi_filter_cutoff, element = "GQ", as.numeric=TRUE)
heatmap.bp(gq, rlabels = FALSE)

#write out thinned vcf
vcfR::write.vcf(gapperi_filter_cutoff, "gapperi.nuc.vcf.gz")


