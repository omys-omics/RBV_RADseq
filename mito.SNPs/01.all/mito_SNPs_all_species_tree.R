setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/RBVs/09.5plates.phylogeography/mito.SNPs/all")


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
mito.SNPs.all <- read.vcfR("mito.SNPs.all.vcf.gz")
# read in popmap
detailed.popmap <- read_delim("~/Documents/Ben's Stuff/0 KU/Dissertation/RBVs/detailed.popmap.20240229.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
popmap <- data.frame(detailed.popmap[detailed.popmap$id %in% colnames(extract.gt(mito.SNPs.all)),])
# make sure popmap is in same order as ids in the vcf
identical(colnames(extract.gt(mito.SNPs.all)), popmap$id)
popmap <- popmap[order(match(popmap$id, colnames(extract.gt(mito.SNPs.all)))),]
identical(colnames(extract.gt(mito.SNPs.all)), popmap$id)



###########################
##   DISTANCE MATRIX   ####
###########################

# Make a genlight of subsamples
genlight.rbv <- vcfR2genlight(mito.SNPs.all)
# Assign populations
pop(genlight.rbv) <- popmap$pop
# Make a distance matrix for splitstree
dist.rbv <- stamppNeisD(genlight.rbv, pop = FALSE)
# Replace NAs with 0
dist.rbv[is.na(dist.rbv)] <- 0
# Replace infinite values
dist.rbv[dist.rbv==Inf] <- 1
#export for splitstree
stamppPhylip(distance.mat=dist.rbv, file="dist.rbv.txt")

# name inds by pop
dist.rbv.pop <- dist.rbv
row.names(dist.rbv.pop) <- paste(popmap$pop, 1:nrow(dist.rbv), sep = "_")
stamppPhylip(distance.mat=dist.rbv.pop, file="dist.rbv.pop.txt")

plot(nj(dist.rbv), "u")



##########################################
###   Assign mitochondrial species  ######
##########################################
# First, based on whether closer to Minn or MAK in dist matrix

mito_assigner <- function(dist.matrix = NULL, popmap = NULL, round = NULL) {
  if (round == 1) {
    # Get ID of MINN and MAK individuals
    MINN <- popmap[popmap$pop == "MINN" & !is.na(popmap$pop),]$id
    MAK <- popmap[popmap$pop == "MAK" & !is.na(popmap$pop),]$id
    print(MINN)
    print(MAK)
    # make new popmap object
    pm <- popmap
    
    # initiate mito species column
    pm$mt_species <- rep(NA, nrow(pm))
    print(pm)
    # make new distance matrix and assign column names
    dist <- dist.matrix
    colnames(dist) <- rownames(dist)
    print(dist)
    
    for (rbv in pm$id) {
      print(rbv)
      if (is.na(rbv)) {
        next
      }
      if (!rbv %in% rownames(dist.matrix)) {
        next
      }
      # find distance to MINN and MAK
      gap.dist <- dist[rbv,MINN]
      rut.dist <- dist[rbv,MAK]
      print(gap.dist)
      print(rut.dist)
      # assign mito species based on whichever distance is smaller
      if (gap.dist < rut.dist) {
        pm[pm$id==rbv & !is.na(pm$id),"mt_species"] <- "gapperi"
      } else if (gap.dist > rut.dist) {
        pm[pm$id==rbv & !is.na(pm$id),"mt_species"] <- "rutilus"
      } else if (gap.dist == rut.dist) {
        pm[pm$id==rbv & !is.na(pm$id),"mt_species"] <- "SAME DISTANCE?????"
      }
    }
    return(pm)
  }

  
  if (round > 1) {
    # Get IDs of rutilus and gapperi mito identified from first round
    rut <- popmap[popmap$mt_species == "rutilus",]$id
    gap <- popmap[popmap$mt_species == "gapperi",]$id
    
    # make new popmap object
    pm <- popmap
    
    # initiate mito species column
    pm$mt_species <- rep(NA, nrow(pm))
    print(pm)
    # make new distance matrix and assign column names
    dist <- dist.matrix
    colnames(dist) <- rownames(dist)
    print(dist)
    
    for (rbv in pm$id) {
      print(rbv)
      if (is.na(rbv)) {
        next
      }
      if (!rbv %in% rownames(dist.matrix)) {
        next
      }
      # find average distance from previously identified gapperi and rutilus
      gap.dist <- mean(dist[rbv,gap])
      rut.dist <- mean(dist[rbv,rut])
      print(gap.dist)
      print(rut.dist)
      # assign mito species based on whichever distance is smaller
      if (gap.dist < rut.dist) {
        pm[pm$id==rbv & !is.na(pm$id),"mt_species"] <- "gapperi"
      } else if (gap.dist > rut.dist) {
        pm[pm$id==rbv & !is.na(pm$id),"mt_species"] <- "rutilus"
      } else if (gap.dist == rut.dist) {
        pm[pm$id==rbv & !is.na(pm$id),"mt_species"] <- "SAME DISTANCE?????"
      }
    }
    return(pm)
  }
  
}

# calculate percent missing data for each individual
percent_missing <- colSums(is.na(mito.SNPs.all@gt))/nrow(mito.SNPs.all@gt)
sort(percent_missing)

# build a matrix of mitochondrial genotypes
genotype_matrix <- extract.gt(mito.SNPs.all)

popmap.mt_species_1  <- mito_assigner(dist.matrix = dist.rbv, popmap = popmap, round = 1)
popmap.mt_species_2 <- mito_assigner(dist.matrix = dist.rbv, popmap = popmap.mt_species_1, round = 2)
popmap.mt_species_3 <- mito_assigner(dist.matrix = dist.rbv, popmap = popmap.mt_species_2, round = 3)

identical(popmap.mt_species_2$mt_species, popmap.mt_species_3$mt_species)

identical(names(percent_missing[-1]), popmap.mt_species_3$id)

popmap.mt_species_3$percent_missing <- percent_missing[-1]

popmap.mt_species_3[popmap.mt_species_3$id=="FN428","mt_species"] <- "Microtus"
popmap.mt_species_3[popmap.mt_species_3$id=="UAM50301","mt_species"] <- "Microtus"
popmap.mt_species_3[popmap.mt_species_3$id=="UAM68273","mt_species"] <- "Microtus"
popmap.mt_species_3[popmap.mt_species_3$id=="UAM23201","mt_species"] <- "Alticola"
popmap.mt_species_3[popmap.mt_species_3$id=="UAM34216","mt_species"] <- "Alticola"


# name inds by pop
dist.rbv.mito <- dist.rbv
row.names(dist.rbv.mito) <- paste(popmap.mt_species_3$mt_species, 1:nrow(dist.rbv.mito), sep = "_")
stamppPhylip(distance.mat=dist.rbv.mito, file="dist.rbv.mt_species.txt")


# remove individuals with >90% missing data
popmap.mt_species_4 <- popmap.mt_species_3[popmap.mt_species_3$percent_missing < 0.9,]

# remove embryos
popmap.mt_species_4 <- popmap.mt_species_4[!popmap.mt_species_4$embryos %in% c("FN611", "FN577", "FN572", "FN2587"),]

# remove Alticola and accidental Microtus
popmap.mt_species_4 <- popmap.mt_species_4[!popmap.mt_species_4$mt_species %in% c("Alticola", "Microtus"),]



# write popmap with mito IDs and lat long to look at in QGIS
write.table(popmap.mt_species_4, file = "popmap.mt_species_4.txt", sep = "\t", col.names = T, row.names = F, quote = F)


##### Redo distance matrix with filtered samples
mito.SNPs.filtered <- mito.SNPs.all[samples=popmap.mt_species_4$id]


# Make a genlight of subsamples
genlight.rbv <- vcfR2genlight(mito.SNPs.filtered)
# Assign populations
pop(genlight.rbv) <- popmap.mt_species_4$pop
# Make a distance matrix for splitstree
dist.rbv <- stamppNeisD(genlight.rbv, pop = FALSE)
# Replace NAs with 0
dist.rbv[is.na(dist.rbv)] <- 0
# Replace infinite values
dist.rbv[dist.rbv==Inf] <- 1
#export for splitstree
stamppPhylip(distance.mat=dist.rbv, file="dist.rbv.filtered.txt")

# name inds by pop
dist.rbv.pop <- dist.rbv
row.names(dist.rbv.pop) <- paste(popmap.mt_species_4$mt_species, 1:nrow(dist.rbv), sep = "_")
stamppPhylip(distance.mat=dist.rbv.pop, file="dist.rbv.mito.filtered.txt")


# MAKE LIST OF SAMPLES WITH RUTILUS MITOGENOMES FOR ALIGNING TO RUTILUS MITOGENOME
rutilus.mt <- popmap.mt_species_4[popmap.mt_species_4$mt_species=="rutilus","id"]
write.table(rutilus.mt, file = "rutilus.mt.txt", quote = F, sep = "\t", col.names = F, row.names = F)

# MAKE LIST OF SAMPLES WITH GAPPERI MITOGENOMES FOR ALIGNING TO GAPPERI MITOGENOME
gapperi.mt <- popmap.mt_species_4[popmap.mt_species_4$mt_species=="gapperi","id"]
write.table(gapperi.mt, file = "gapperi.mt.txt", quote = F, sep = "\t", col.names = F, row.names = F)



