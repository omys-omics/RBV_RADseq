setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/RBVs/09.5plates.phylogeography/mito.SNPs/both")


library(vcfR)
library(adegenet)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(diveRsity)
library(RColorBrewer)
library(gtools)
library(ggplot2) 
library(vecsets)
library(readr)
library(SNPfiltR)
library(readxl)


# VCF file for all samples:

# read in vcf file
mito.SNPs.both <- read.vcfR("mito.SNPs.both.vcf.gz")
# read in popmap
detailed.popmap <- read_delim("../all/popmap.mt_species_3.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
popmap <- data.frame(detailed.popmap[detailed.popmap$id %in% colnames(extract.gt(mito.SNPs.both)),])

mito.SNPs.both
table(popmap$mt_species)


# Only keep rut inds that passed previous rut threshold

#rut.pop.low.miss <- read_delim("../rutilus/popmap.mt_rutilus_4.txt", 
#                               delim = "\t", escape_double = FALSE, 
#                               trim_ws = TRUE)

#popmap <- popmap[popmap$id %in% rut.pop.low.miss$id | !popmap$mt_species=="rutilus",]
#mito.SNPs.both <- mito.SNPs.both[samples=popmap$id]
#mito.SNPs.both <- min_mac(mito.SNPs.both, min.mac = 1)




# make sure popmap is in same order as ids in the vcf
identical(colnames(extract.gt(mito.SNPs.both)), popmap$id)
popmap <- popmap[order(match(popmap$id, colnames(extract.gt(mito.SNPs.both)))),]
identical(colnames(extract.gt(mito.SNPs.both)), popmap$id)


# calculate percent missing data for each individual
percent_missing <- colSums(is.na(mito.SNPs.both@gt))/nrow(mito.SNPs.both@gt)
sort(percent_missing)
identical(names(percent_missing[-1]), popmap$id)
popmap$percent_missing <- percent_missing[-1]



# build a matrix of mitochondrial genotypes
genotype_matrix <- extract.gt(mito.SNPs.both)


#plot depth per snp and per sample
dp <- extract.gt(mito.SNPs.both, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)

# Loci kept:
mito.SNPs.both@fix
both.snps <- mito.SNPs.both@fix[,"ID"]
both.loci <- unique(sapply(strsplit(both.snps, ":"), '[', 1))
both.loci
length(both.loci)

# Write vcf of low miss samples
vcfR::write.vcf(mito.SNPs.both, "mito.SNPs.both.vcf.gz")




# Build IQTREE




###########################
##   DISTANCE MATRIX   ####
###########################

# Make a genlight of subsamples
genlight.rbv <- vcfR2genlight(mito.SNPs.both)
# Assign populations
pop(genlight.rbv) <- popmap$pop
# Make a distance matrix for splitstree
dist.rbv <- stamppNeisD(genlight.rbv, pop = FALSE)
# Replace NAs with 0
dist.rbv[is.na(dist.rbv)] <- 0
# Replace infinite values
dist.rbv[dist.rbv==Inf] <- 1
#export for splitstree
stamppPhylip(distance.mat=dist.rbv, file="dist.both.txt")

# name inds by pop
dist.rbv.pop <- dist.rbv
row.names(dist.rbv.pop) <- paste(popmap$pop, 1:nrow(dist.rbv), sep = "_")
stamppPhylip(distance.mat=dist.rbv.pop, file="dist.both.pop.txt")






