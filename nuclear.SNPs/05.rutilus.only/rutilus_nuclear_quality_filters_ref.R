setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/RBVs/09.5plates.phylogeography/nuclear.SNPs/05.rutilus.only")

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

# keep only rutilus
detailed.popmap <- detailed.popmap[detailed.popmap$nuc_species=="rutilus",]


# make popmap
popmap <- detailed.popmap[,c("id", "pop")]

# filter vcfR for correct rutilus individuals
rutilus <- all_species[samples=popmap$id]


# make sure popmap is in same order as ids in the vcf
identical(colnames(extract.gt(rutilus)), popmap$id)
popmap <- popmap[order(match(popmap$id, colnames(extract.gt(rutilus)))),]
identical(colnames(extract.gt(rutilus)), popmap$id)


#plot depth per snp and per sample
dp <- extract.gt(rutilus, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)

#plot genotype quality per snp and per sample
gq <- extract.gt(rutilus, element = "GQ", as.numeric=TRUE)
heatmap.bp(gq, rlabels = FALSE)

# average depth per sample
depth <- dp
depth[is.na(depth)] <- 0
mean.depth <- colSums(depth)/nrow(depth)
hist(mean.depth, breaks = 265)

#visualize distributions
hard_filter(vcfR=rutilus)

# hard filter to minimum depth of 4, minimum genotype quality of 30
rutilus_filter <- hard_filter(vcfR=rutilus, depth = 4, gq = 30)

# allele balance filter
filter_allele_balance(rutilus_filter)

#visualize and pick appropriate max depth cutoff
max_depth(rutilus_filter)
rutilus_filter <- max_depth(rutilus_filter, maxdepth = 400)

#remove invariant sites
rutilus_filter<-min_mac(rutilus_filter, min.mac = 1)


#visualize missing data by SNP and the effect of various cutoffs on the missingness of each sample
missing_by_snp(rutilus_filter)

#verify that missing data is not driving clustering patterns among samples
miss<-assess_missing_data_pca(vcfR=rutilus_filter, popmap = popmap, thresholds = 0.8, clustering = FALSE)


#choose a value that retains an acceptable amount of missing data in each sample, and maximizes SNPs retained while minimizing overall missing data, and filter vcf
rutilus_filter_cutoff<-missing_by_snp(rutilus_filter, cutoff = 0.8)
popmap_filter<-popmap[popmap$id %in% colnames(rutilus_filter@gt),]

#run function to visualize missing data by individual
dev.off()
pdf(file = "missing_by_sample.pdf", width = 6, height = 10)
missing_by_sample(vcfR=rutilus_filter_cutoff)
dev.off()


#plot depth per snp and per sample
dp <- extract.gt(rutilus_filter_cutoff, element = "DP", as.numeric=TRUE)
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
gq <- extract.gt(rutilus_filter_cutoff, element = "GQ", as.numeric=TRUE)
heatmap.bp(gq, rlabels = FALSE)

#write out thinned vcf
vcfR::write.vcf(rutilus_filter_cutoff, "rutilus.nuc.vcf.gz")


