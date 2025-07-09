setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/RBVs/09.5plates.phylogeography/nuclear.SNPs/01.quality.filters")

library(vcfR)
library(SNPfiltR)
library(readr)


############################################################
####### Visualize and filter SNPs by quality and missingness
############################################################


# Read in vcf
all_species <- read.vcfR("all.nuc.recode.vcf.gz")

#generate popmap file. Two column popmap with the same format as stacks, and the columns must be named 'id' and 'pop'
detailed.popmap <- read_delim("~/Documents/Ben's Stuff/0 KU/Dissertation/RBVs/detailed.popmap.20240229.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

# popmap of what is sequenced
sequenced.popmap <- detailed.popmap[detailed.popmap$id %in% colnames(extract.gt(all_species)),]
write.table(sequenced.popmap, file = "RADseq.samples.txt", quote = F, sep = "\t")

# remove accidental Microtus
detailed.popmap <- detailed.popmap[!detailed.popmap$id %in% c("FN428", "UAM50301", "UAM68273"),]

# remove embryos
detailed.popmap <- detailed.popmap[!detailed.popmap$embryos %in% c("FN611", "FN572", "FN577", "FN2587"),]

# make popmap
popmap <- data.frame(detailed.popmap[detailed.popmap$id %in% colnames(extract.gt(all_species)),])
popmap <- popmap[,c("id", "pop")]

# remove accidental Microtus and embryos
all_species <- all_species[samples=popmap$id]


# make sure popmap is in same order as ids in the vcf
identical(colnames(extract.gt(all_species)), popmap$id)
popmap <- popmap[order(match(popmap$id, colnames(extract.gt(all_species)))),]
identical(colnames(extract.gt(all_species)), popmap$id)


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
hist(mean.depth, breaks = 265)

#visualize distributions
hard_filter(vcfR=all_species)

# hard filter to minimum depth of 4, minimum genotype quality of 30
all_species_filter <- hard_filter(vcfR=all_species, depth = 4, gq = 30)

# allele balance filter
filter_allele_balance(all_species_filter)

#visualize and pick appropriate max depth cutoff
max_depth(all_species_filter)
all_species_filter <- max_depth(all_species_filter, maxdepth = 400)

#remove invariant sites
all_species_filter<-min_mac(all_species_filter, min.mac = 1)


#run function to visualize samples above the threshold we want from the vcf
pdf(file = "missing_by_sample.pdf", width = 2, height = 30)
missing_by_sample(vcfR=all_species_filter, popmap = popmap, cutoff = 0.984)
dev.off()


#0.984 gets rid of top 14
all_species_filter <- missing_by_sample(vcfR=all_species_filter, cutoff = 0.984) 


#subset popmap to only include retained individuals
popmap_filter<-popmap[popmap$id %in% colnames(all_species_filter@gt),]

#see which individuals and pops were dropped
setdiff(popmap$id, popmap_filter$id)
setdiff(popmap$pop, popmap_filter$pop)

#remove invariant sites
all_species_filter<-min_mac(all_species_filter, min.mac = 1)

#visualize missing data by SNP and the effect of various cutoffs on the missingness of each sample
missing_by_snp(all_species_filter)

#verify that missing data is not driving clustering patterns among samples
miss<-assess_missing_data_pca(vcfR=all_species_filter, popmap = popmap_filter, thresholds = 0.8, clustering = FALSE)


#choose a value that retains an acceptable amount of missing data in each sample, and maximizes SNPs retained while minimizing overall missing data, and filter vcf
all_species_filter_cutoff<-missing_by_snp(all_species_filter, cutoff = 0.8)
popmap_filter<-popmap[popmap$id %in% colnames(all_species_filter@gt),]

#run function to visualize samples above the threshold we want from the vcf
pdf(file = "missing_by_sample_2.pdf", width = 8, height = 30)
missing_by_sample(vcfR=all_species_filter_cutoff)
dev.off()


#plot depth per snp and per sample
dp <- extract.gt(all_species_filter_cutoff, element = "DP", as.numeric=TRUE)
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
gq <- extract.gt(all_species_filter_cutoff, element = "GQ", as.numeric=TRUE)
heatmap.bp(gq, rlabels = FALSE)

#write out thinned vcf
vcfR::write.vcf(all_species_filter_cutoff, "rbv.nuc.all.vcf.gz")


