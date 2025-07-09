setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/RBVs/09.5plates.phylogeography/nuclear.SNPs/02.nuc.analyses")


library(hierfstat)
library(adegenet)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(diveRsity)
library(dartR)
library(RColorBrewer)
library(gtools)
library(reshape2)
library(ggplot2)
library(cowplot)
library(dplyr)
library(apex)
library(vecsets)
library(readr)
library(gridExtra)
library(SNPfiltR)
library(readxl)
library(strataG)
library(pophelper)
library(scales)
library(triangulaR)
library(ggpubr)
library(vcfR)


##########################################
###      Define helpful functions     ####
##########################################
vcfR2splitstree <- function(vcfR = NULL, popmap = NULL, file.prefix = NULL) {
  if (!identical(colnames(extract.gt(vcfR)), popmap$id)) {
    stop("order of inds in vcfR and popmap doesn't match")
  }
  # Make a genlight of subsamples
  gl <- vcfR2genlight(vcfR)
  # Assign populations
  pop(gl) <- popmap$pop
  # Make a distance matrix for splitstree
  dist <- stamppNeisD(gl, pop = FALSE)
  #export for splitstree
  stamppPhylip(distance.mat=dist, file=paste0(file.prefix, ".txt"))

  # name inds by pop
  dist.pop <- dist
  row.names(dist.pop) <- paste(popmap$pop, 1:nrow(popmap), sep = "_")
  stamppPhylip(distance.mat=dist.pop, file=paste0(file.prefix, ".pop.txt"))
  return(dist)
}


pcaSubset <- function(vcfr = NULL, popmap = NULL, morph_type = NA, morphs = NA, species_type = NA, species = NA) {

  # Make a popmap of only requested groups
  pm <- popmap
  if (!is.na(morph_type)) {
    pm <- pm[pm[[morph_type]] %in% morphs,]
  }
  if (!is.na(species_type)) {
    pm <- pm[pm[[species_type]] %in% species,]
  }
  # Make a vcfR of only requested groups
  v <- min_mac(vcfr[samples = pm$id], min.mac = 1)
  # Make a genlight of only requested groups
  gl <- vcfR2genlight(v)
  # Check that inds line up
  if(!identical(gl@ind.names, pm$id)) {
    stop("genlight and popmap ids not the same")
  }
  # Assign populations
  pop(gl) <- pm$pop
  #perform PCA
  pca<-glPca(gl, nf = 10)
  #pull pca scores out of df
  pca.scores<-as.data.frame(pca$scores)
  pca.scores$id <- pm$id
  pca.scores$pop <- pm$pop
  pca.scores$date <- pm$date
  pca.scores$locality <- pm$locality
  pca.scores$lat <- pm$lat
  pca.scores$nuc_species <- pm$nuc_species
  pca.scores$nuc_morph <- pm$nuc_morph
  pca.scores$mt_species <- pm$mt_species
  pca.scores$mt_morph <- pm$mt_morph
  #calculate missingness by individual
  miss<-colSums(is.na(v@gt))/nrow(v@gt)
  #add missingness to df
  pca.scores$missing<-miss[-1]
  #record variance percentages explained
  var.frac<- pca$eig/sum(pca$eig)*100
  # get heterozygosity for each individual
  het <- is.het(extract.gt(v))
  for (rbv in rownames(pca.scores)) {
    pca.scores[rbv, "het"] <- sum(het[,rbv])/length(het[,rbv])
  }

  # put pca.scores and var.frac into a list to use outside of function
  pca.list <- list(pca.scores, var.frac)
  # return list
  return(pca.list)
}


#  function for plotting structure results
betterStructurePlot <- function(
  q.mat, 
  ind.col = 1,
  pop.col = 3, 
  prob.col = 4, 
  sort.probs = TRUE,
  label.pops = FALSE,
  label.inds = TRUE,
  col = NULL, 
  horiz = TRUE, 
  type = NULL,
  legend.position = c("top", "left", "right", "bottom", "none"),
  plot = TRUE
) {
  
  legend.position <- match.arg(legend.position)
  
  # convert q.mat to sorted data.table
  prob.cols <- prob.col:ncol(q.mat)
  qm <- as.data.frame(q.mat)[, c(ind.col, pop.col, prob.cols), drop = FALSE]
  qm[, 2] <- factor(
    qm[, 2], 
    levels = sort(unique(qm[, 2]), decreasing = !horiz)
  )
  sort.cols <- c(2, if(sort.probs) 3:ncol(qm) else NULL)
  i <- do.call(
    order, 
    c(as.list(qm[, sort.cols, drop = FALSE]), decreasing = TRUE)
  )
  qm <- qm[i, ] # this reverses the order of pops for plotting
  qm$x <- 1:nrow(qm)
  
  # Get population frequencies, centers and dividing points
  pop.freq <- table(qm[, 2])
  levels(qm[, 2]) <- paste(
    levels(qm[, 2]), "\n(n = ", pop.freq, ")", sep = ""
  )
  pop.cntr <- tapply(qm$x, qm[, 2], mean)
  ind.cntr <- tapply(qm$x, qm[, 1], mean)
  pop.div <- rev(tapply(qm$x, qm[, 2], min))[-1] - 0.5
  
  # Create data.frame for plotting
  df <- melt(qm[, c("x", colnames(qm)[-ncol(qm)])], id.vars = c(1:3),
             variable.name = "Group", value.name = "probability")
  colnames(df)[1:3] <- c("x", "id", "population")
  df <- df[order(-as.numeric(df$Group), df$probability), ]
  
  type <- if(is.null(type)) {
    if(nrow(df) <= 100) "bar" else "area"
  } else {
    match.arg(type, c("bar", "area"))
  }
  
  # Plot stacked bar graphs
  g <- ggplot2::ggplot(df, ggplot2::aes_string("x", "probability")) +  
    switch(
      type,
      area = ggplot2::geom_area(
        ggplot2::aes_string(fill = "Group"), 
        stat = "identity"
      ),
      bar = ggplot2::geom_bar(
        ggplot2::aes_string(fill = "Group"), 
        stat = "identity"
      )
    ) +
    ggplot2::ylab("Pr(Group Membership)") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = legend.position,
      legend.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    )
  if(label.pops) {
    g <- g + 
      ggplot2::geom_vline(xintercept = pop.div, size = 1.5) +
      ggplot2::scale_x_continuous(
        name = "", 
        breaks = pop.cntr, 
        labels = names(pop.cntr),
        expand = c(0, 0)
      )
  } else if (label.inds) {
    g <- g + 
      ggplot2::geom_vline(xintercept = pop.div, size = 1.5) +
      ggplot2::scale_x_continuous(
        name = "", 
        breaks = ind.cntr, 
        labels = names(ind.cntr),
        expand = c(0, 0)
      ) +
      ggplot2::theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  } else {
    g <- g + 
      ggplot2::xlab("") + 
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }
  if(horiz) g <- g + ggplot2::coord_flip()
  if(!is.null(col)) g <- g + ggplot2::scale_fill_manual(values = col)
  
  if(plot) print(g)
  invisible(g)
}


#########################################
######        Setup data       ##########
#########################################

# VCF file for all samples:

# read in vcf file
rbv_filtered <- read.vcfR("../01.quality.filters/rbv.nuc.all.vcf.gz")
# read in full popmap
detailed.popmap <- read_delim("~/Documents/Ben's Stuff/0 KU/Dissertation/RBVs/detailed.popmap.20240229.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

# get mito popmap
mito.popmap <- read_delim("../../mito.SNPs/all/popmap.mt_species_4.txt",
                              delim = "\t", escape_double = FALSE,
                              trim_ws = TRUE)


popmap <- data.frame(detailed.popmap[detailed.popmap$id %in% colnames(extract.gt(rbv_filtered)),])



# add mito calls to full popmap
popmap$mt_species <- rep("UNKNOWN", nrow(popmap))
for (rbv in popmap$id) {
  if (rbv %in% mito.popmap$id) {
    popmap[popmap$id==rbv, "mt_species"] <- mito.popmap[mito.popmap$id==rbv, "mt_species"]
  }
}
  
  

# make sure popmap is in same order as ids in the vcf
identical(colnames(extract.gt(rbv_filtered)), popmap$id)
popmap <- popmap[order(match(popmap$id, colnames(extract.gt(rbv_filtered)))),]
identical(colnames(extract.gt(rbv_filtered)), popmap$id)


# choose only one of the duplicated individuals, whichever has least missing data
popmap$nuc_missing_data <- (colSums(is.na(rbv_filtered@gt))/nrow(rbv_filtered@gt))[-1]
popmap <- popmap[!popmap$id %in% c("FN2564_2", "FN2586_2"),]
rbv_filtered <- min_mac(rbv_filtered[samples=popmap$id], min.mac = 1)
identical(colnames(extract.gt(rbv_filtered)), popmap$id)



######################################
##   Remove Alticola for structure  ##
######################################

no.russ.popmap <- popmap[popmap$pop != "RUSS",]
vcfR.for.structure <- min_mac(rbv_filtered[samples=no.russ.popmap$id], min.mac = 1)
# write out vcf
vcfR::write.vcf(vcfR.for.structure, "rbv.for.structure.vcf.gz")
# write popmap for structure
write.table(no.russ.popmap[,c("id", "pop")], "popmap.structure.txt", quote = F, sep = "\t", row.names = F, col.names = T)

# run structure outside of R



##################################
#     Plot structure results     #
##################################

for (file in list.files("../03.structure/01.output.1M.K2_8/", pattern = "*_f", full.names = T)) {
  # Get number to name structure run
  num_of_k <- paste(sapply(strsplit(file, split='_', fixed=TRUE), `[`, 3), sapply(strsplit(file, split='_', fixed=TRUE), `[`, 4), sep = "_")
  # read in structure results
  raw <- structureRead(file)
  # Assign pop names
  for (i in 1:nrow(raw$q.mat)) {
    current_ind <- raw$q.mat[i,"id"]
    raw$q.mat[i,"orig.pop"] <- popmap[popmap$id == current_ind,"pop"]
  }
  # Rename to structure run
  assign(paste0(num_of_k, "_structure"), raw)
}


# create color palette
color_ramp <- colorRampPalette(c("orange", "blue", "green", "red1", "yellow", "purple"))
cols <- color_ramp(20)



# Evanno method
slist <- readQ(list.files("../03.structure/01.output.1M.K2_8/", pattern = "*_f", full.names = T))
tabQ <- tabulateQ(slist)
sumQ <- summariseQ(tabQ)
evanno.data <- evannoMethodStructure(sumQ)
evanno.plot <- evannoMethodStructure(sumQ, returnplot = T, exportplot = F, returndata = F, basesize=12, linesize=0.7)
grid.arrange(evanno.plot)


# view structure results labeled by pop
betterStructurePlot(K2_rep1_structure$q.mat, horiz = F, label.pops = T, legend.position = "top", col = color_ramp(2), sort.probs = T, type = "bar")

# view structure results labeled by individual
betterStructurePlot(K2_rep1_structure$q.mat, horiz = F, legend.position = "top", col = color_ramp(2), sort.probs = T, type = "bar")


# save plots as image
for (k in ls(pattern = "K*_structure")) {
  # get number of k
  num <- sapply(strsplit(paste(sapply(strsplit(k, split='_', fixed=TRUE), `[`, 1), sep = "_"), split = "K", fixed =TRUE), `[`, 2)
  # get number of rep
  repp <- sapply(strsplit(paste(sapply(strsplit(k, split='_', fixed=TRUE), `[`, 2), sep = "_"), split = "p", fixed =TRUE), `[`, 2)
  # assign temp object for data
  kdata <- get(k)
  # plot with individuals labels
  plot <- betterStructurePlot(kdata$q.mat, horiz = F, legend.position = "top", col = color_ramp(num), sort.probs = T, type = "bar")
  ggplot2::ggsave(filename = paste0("../03.structure/all_structure_plots/", k, "_plot_ind_labels.png"), plot, width = 15, height = 5, path = getwd())
  # plot with pop labels
  plot <- betterStructurePlot(kdata$q.mat, label.pops = T, horiz = F, legend.position = "top", col = color_ramp(num), sort.probs = T, type = "bar")
  ggplot2::ggsave(filename = paste0("../03.structure/all_structure_plots/", k, "_plot_pops_labels.png"), plot, width = 15, height = 5, path = getwd())
}


#################################
##   Assign structure species  ##
#################################
species.q.mat <- K2_rep1_structure$q.mat
identical(no.russ.popmap$id, species.q.mat$id)
no.russ.popmap <- no.russ.popmap[order(match(no.russ.popmap$id, species.q.mat$id)),]
identical(no.russ.popmap$id, species.q.mat$id)

no.russ.popmap$q1 <- species.q.mat$Group.1
no.russ.popmap$q2 <- species.q.mat$Group.2

no.russ.popmap[no.russ.popmap$q1 > 0.8 & !is.na(no.russ.popmap$q1), "nuc_species"] <- "rutilus"
no.russ.popmap[no.russ.popmap$q2 > 0.8 & !is.na(no.russ.popmap$q1), "nuc_species"] <- "gapperi"
no.russ.popmap[is.na(no.russ.popmap$nuc_species),"nuc_species"] <- "hybrid"


######################################
####     Check on missing data    ####
######################################

# Assess correlation between depth and missing data
read.count <- read.table("reads.txt", col.names = c("id", "reads"))
read.count <- read.count[read.count$id %in% no.russ.popmap$id,]

identical(no.russ.popmap$id, read.count$id)
read.count <- read.count[order(match(read.count$id, no.russ.popmap$id)),]
identical(no.russ.popmap$id, read.count$id)

no.russ.popmap$read.count <- read.count$reads

summary(lm(nuc_missing_data ~ read.count, data = no.russ.popmap))

ggplot(no.russ.popmap, aes(y=nuc_missing_data, x=read.count)) +
  geom_point(cex = 2, alpha=1) +
  #geom_smooth(method = "lm", se = FALSE) +
  ylab(paste("Proportion nuclear missing data")) +
  xlab(paste("Read count")) +
  scale_x_log10(breaks=c(300000, 1000000, 3000000, 10000000, 30000000)) +
  labs(title = "")

# Assess correlation between missing data and q1
summary(lm(q1 ~ nuc_missing_data, no.russ.popmap))

q1_r <- c(summary(lm(q1 ~ nuc_missing_data, no.russ.popmap))$r.squared,
          summary(lm(q1 ~ nuc_species, no.russ.popmap))$r.squared,
          summary(lm(q1 ~ nuc_species + nuc_missing_data, no.russ.popmap))$r.squared)

q1_r <- data.frame(t(q1_r))

colnames(q1_r) <- c("missing_data", "species", "species+missing_data")
q1_r$improvement <- q1_r$`species+missing_data` - q1_r$species
q1_r <- round(q1_r, digits = 8)

write.table(q1_r, file = "variance_from_missingness.txt", quote = F, sep = "\t")

q1_miss <- ggplot(no.russ.popmap, aes(x=nuc_missing_data, y=q1, color=nuc_species)) +
  geom_point(cex = 2, alpha=1) +
  #geom_smooth(method = "lm", se = FALSE) +
  #annotate("text", x=0.7, y=1, label= "R2=0.003067") +
  #annotate("text", x=0.7, y=0.95, label= "p=0.4084") +
  xlab(paste("Proportion nuclear missing data")) +
  ylab(paste("Q1 from structure")) +
  labs(title = "")

ggplot2::ggsave(filename = "missing_q1.png", q1_miss, width = 7, height = 4, path = getwd())


ggplot(no.russ.popmap, aes(x=nuc_missing_data)) +
  geom_histogram(bins = 40) +
  xlab(paste("Proportion nuclear missing data")) +
  ylab(paste("Count")) +
  labs(title = "")

summary(no.russ.popmap$nuc_missing_data)
sum(no.russ.popmap$nuc_missing_data < 0.15) / nrow(no.russ.popmap)

#####################################
##   whitelist loci for iqtree   ####
#####################################

################
#   all inds  ##
################

# get locus names in vcfR object
whitelist.rbv <- rbv_filtered@fix[,"ID"]
whitelist.rbv <- sapply(strsplit(whitelist.rbv, split = ":"), "[", 1)
# make sure each locus name is unique
length(unique(whitelist.rbv)) == length(whitelist.rbv)
# write out table of locus names
write.table(whitelist.rbv, file = "whitelist.rbv.txt", quote = F, row.names = F, col.names = F)
# write popmap for iq tree (each ind is its own pop)
iqtree.popmap.rbv <- cbind(popmap$id, popmap$id)
write.table(iqtree.popmap.rbv, file = "iqtree.popmap.rbv.txt", sep = "\t", quote = F, row.names = F, col.names = F)

########################
#  only inds in mito  ##
#  trees + Alticola   ##
########################
rutilus.mito.popmap <- read_delim("../../mito.SNPs/rutilus/popmap.mt_rutilus_4.txt", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)

gapperi.mito.popmap <- read_delim("../../mito.SNPs/gapperi/popmap.mt_gapperi.txt", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)


mt.inds <- c(rutilus.mito.popmap$id, gapperi.mito.popmap$id, "UAM23201", "UAM34216")
# make popmap
iqtree.popmap.rbv.mito.inds <- cbind(mt.inds, mt.inds)
write.table(iqtree.popmap.rbv.mito.inds, file = "iqtree.popmap.rbv.mito.inds.txt", sep = "\t", quote = F, row.names = F, col.names = F)


# run IQ-TREE outside of R



#############################################
##   sort structure by order on iqtree   ####
#############################################

# get q matrix
q <- K2_rep1_structure$q.mat

# add species column
q$species <- rep(NA, nrow(q))
q[q$Group.1 > q$Group.2, "species"] <- "rutilus"
q[q$Group.2 > q$Group.1, "species"] <- "gapperi"

#melt into dataframe for plotting
df <- melt(q, id.vars = c(1:3,6), variable.name = "Group", value.name = "ancestry")


# get order of inds on tree
tree.order <- read.table(file = "../04.iqtree/rbv.tree.order.txt")[,1]
# without the Alticola
tree.order <- tree.order[-1:-2]

# sort ids by order in iqtree
# make sure popmap is in same order as ids in the vcf
identical(df$id, tree.order)
df <- df[order(match(df$id, tree.order)),]
identical(unique(df$id), tree.order)


# in order to add black outlines to each bar
pop.div <- 1:nrow(df) - 0.5


# plot big with labels
big.structure.plot <- ggplot(df, aes(x=ancestry, y=id, fill = Group)) +
  geom_bar(stat = "identity") +
  xlab("Ancestry") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(factor(df$id, tree.order)))) + 
  theme(
    axis.ticks.y = ggplot2::element_blank(),
    legend.position = "none",
    legend.title = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()
    ) +
  geom_hline(yintercept = pop.div, size = 0.1, color = "white") +
  scale_fill_manual(values = c("#D8B365", "#5AB4AC"))


# plot small, both species together
structure.plot <- ggplot(df, aes(x=ancestry, y=id, fill = Group)) +
  geom_bar(stat = "identity") +
  xlab("Ancestry") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(factor(df$id, tree.order)))) + 
  theme(
    axis.ticks.y = ggplot2::element_blank(),
    legend.position = "none",
    legend.title = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.text.y = element_blank()
  ) +
  geom_hline(yintercept = pop.div, size = 0.1, color = "white") +
  scale_fill_manual(values = c("#D8B365", "#5AB4AC"))

rutilus.structure <- ggplot(df[df$species=="rutilus",], aes(x=ancestry, y=id, fill = Group)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  xlab("Ancestry") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(factor(df$id, tree.order)))) + 
  theme(
    axis.ticks.y = ggplot2::element_blank(),
    legend.position = "none",
    legend.title = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.text.y = element_blank()
  ) +
  geom_hline(yintercept = pop.div, size = 0.1, color = "white") +
  scale_fill_manual(values = c("#D8B365", "#5AB4AC"))

gapperi.structure <- ggplot(df[df$species=="gapperi",], aes(x=ancestry, y=id, fill = Group)) +
  geom_bar(stat = "identity") +
  xlab("Ancestry") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(factor(df$id, tree.order)))) + 
  theme(
    axis.ticks.y = ggplot2::element_blank(),
    legend.position = "none",
    legend.title = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.text.y = element_blank()
  ) +
  geom_hline(yintercept = pop.div, size = 0.1, color = "white") +
  scale_fill_manual(values = c("#D8B365", "#5AB4AC"))

ggplot2::ggsave(filename = "big.structure.plot.pdf", big.structure.plot, width = 2, height = 24, path = getwd())
ggplot2::ggsave(filename = "structure.plot.pdf", structure.plot, width = 1, height = 8, path = getwd())
ggplot2::ggsave(filename = "rutilus.structure.plot.pdf", rutilus.structure, width = 1, height = 8, path = getwd())
ggplot2::ggsave(filename = "gapperi.structure.plot.pdf", gapperi.structure, width = 1, height = 8, path = getwd())



#################################
##   Make popmap for mapping   ##
#################################

nuc.mito.map <- no.russ.popmap[no.russ.popmap$id %in% iqtree.popmap.rbv.mito.inds,]
setdiff(iqtree.popmap.rbv.mito.inds[,1], nuc.mito.map$id)

#add UAM44610 to map, confirmed nuclear rutilus with IQTREE
UAM44610 <- detailed.popmap[detailed.popmap$id=="UAM44610" & !is.na(detailed.popmap$id),]
UAM44610$mt_species <- "rutilus"
UAM44610$nuc_missing_data <- NA
UAM44610$q1 <- NA
UAM44610$q2 <- NA
UAM44610$nuc_species <- "rutilus"
nuc.mito.map <- rbind(nuc.mito.map, UAM44610)

#add column for unique combinations
nuc.mito.map$class <- paste(nuc.mito.map$pop, nuc.mito.map$nuc_species, nuc.mito.map$mt_species, sep = "_")

write.table(nuc.mito.map, file = "nuc.mito.map.txt", sep = "\t", quote = F, col.names = T, row.names = F)


#################################################
##    structure gapperi and rutilus separate   ##
#################################################

# write full popmap with every individual in structure run
write.table(no.russ.popmap, file = "no.russ.popmap.txt", sep = "\t", quote = F, col.names = T, row.names = F)

# refilter SNPs for rutilus and gapperi for separate analyses (other R scripts)







######################
##   TRIANGLE PLOTS ##
######################

# make popmap for triangle plot
tri.popmap <- no.russ.popmap
tri.popmap$tri.pop <- tri.popmap$pop
tri.popmap[tri.popmap$id %in% c("MSB147123", "MSB147516", "MSB158085", "MSB198875"),]$tri.pop <- "rutilus.parental"
tri.popmap[tri.popmap$id %in% c("MSB157962", "MSB199211", "MSB199275", "UAM59916"),]$tri.pop <- "gapperi.parental"
tri.popmap <- tri.popmap[,c("id", "tri.pop", "nuc_species")]
colnames(tri.popmap) <- c("id", "pop", "species")
tri.popmap <- tri.popmap[!tri.popmap$pop %in% c("MAK", "MINN"),]

# make vcfR for triangle plot
rbv_tri <- rbv_filtered[samples = tri.popmap$id]
rbv_tri <- min_mac(rbv_tri, min.mac = 1)

# find AIMs and calculate hybrid index and interclass heterozygosity
rbv_diff <- alleleFreqDiff(vcfR = rbv_tri, pm = tri.popmap, p1 = "gapperi.parental", p2 = "rutilus.parental", difference = 1)
hi_het <- hybridIndex(vcfR = rbv_diff, pm = tri.popmap, p1 = "gapperi.parental", p2 = "rutilus.parental")
triangulaR::triangle.plot(hi_het)
freq_diff <- specFreqDiff(vcfR = rbv_tri, pm = tri.popmap, p1 = "gapperi.parental", p2 = "rutilus.parental")
hist(freq_diff$diff)

spec_freq_diff <- ggplot(data = freq_diff,  aes(x=diff)) +
  geom_histogram(color="#FFFFFF", alpha=1, position = 'identity', bins = 9) +
  xlab(paste("Observed Allele Frequency Difference")) +
  ylab(paste("Count")) +
  #ylim(0,8000) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position = "none") +  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

ggsave("spec_freq_diff.pdf", plot = spec_freq_diff, path = "./", width = 6, height = 4)




# write triangle vcfR and popmap
vcfR::write.vcf(rbv_tri, file = "rbv_triangle.vcf.gz")
write.table(tri.popmap, file = "rbv_tri_popmap.txt", sep = "\t", quote = F, row.names = F, col.names = T, append = F)

# add nuc species to hi_het
identical(hi_het$id, tri.popmap$id)
hi_het$species <- tri.popmap$species


tri.BC <- ggplot(hi_het[hi_het$pop == "BC",], aes(x = hybrid.index, y = heterozygosity, color = as.factor(species))) +
  geom_segment(aes(x = 0.5, xend = 1, y = 1, yend = 0), color = "black") +
  geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 1), color = "black") +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
  stat_function(fun = function(hi) 2 * hi * (1 - hi), xlim = c(0, 1), color = "black", linetype = "dashed") +
  geom_jitter(cex = 2, alpha = 1, width = 0, height = 0) +
  guides(shape = guide_legend(override.aes = list(size = 5), order = 2, label.theme = element_text(face = "italic"))) +
  xlab(paste("Hybrid Index")) +
  ylab(paste("Interclass Heterozygosity")) +
  labs(title = "British Columbia") +
  scale_color_manual("species", values = c("#5AB4AC", "#D8B365")) +
  ylim(c(-0.05, 1.05)) + xlim(c(-0.05, 1.05)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")


tri.NWT <- ggplot(hi_het[hi_het$pop == "NWT",], aes(x = hybrid.index, y = heterozygosity, color = as.factor(species))) +
  geom_segment(aes(x = 0.5, xend = 1, y = 1, yend = 0), color = "black") +
  geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 1), color = "black") +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
  stat_function(fun = function(hi) 2 * hi * (1 - hi), xlim = c(0, 1), color = "black", linetype = "dashed") +
  geom_jitter(cex = 2, alpha = 1, width = 0, height = 0) +
  guides(shape = guide_legend(override.aes = list(size = 5), order = 2, label.theme = element_text(face = "italic"))) +
  xlab(paste("Hybrid Index")) +
  ylab(paste("Interclass Heterozygosity")) +
  labs(title = "Northwest Territories") +
  scale_color_manual("species", values = c("#5AB4AC", "black", "#D8B365")) +
  ylim(c(-0.05, 1.05)) + xlim(c(-0.05, 1.05)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")


tri.SEAK <- ggplot(hi_het[hi_het$pop == "SEAK",], aes(x = hybrid.index, y = heterozygosity, color = as.factor(species))) +
  geom_segment(aes(x = 0.5, xend = 1, y = 1, yend = 0), color = "black") +
  geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 1), color = "black") +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
  stat_function(fun = function(hi) 2 * hi * (1 - hi), xlim = c(0, 1), color = "black", linetype = "dashed") +
  geom_jitter(cex = 2, alpha = 1, width = 0, height = 0) +
  guides(shape = guide_legend(override.aes = list(size = 5), order = 2, label.theme = element_text(face = "italic"))) +
  xlab(paste("Hybrid Index")) +
  ylab(paste("Interclass Heterozygosity")) +
  labs(title = "Southeast Alaska") +
  scale_color_manual("species", values = c("#5AB4AC", "#D8B365")) +
  ylim(c(-0.05, 1.05)) + xlim(c(-0.05, 1.05)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")


all.tri <- ggarrange(tri.NWT, tri.BC, tri.SEAK,
          labels = c("A", "B", "C"),
          font.label = list(size = 14),
          #label.y = -0.2,
          label.x = 0.1,
          ncol = 1, nrow = 3)


ggplot2::ggsave(filename = "triangle.plots.pdf", all.tri, width = 3, height = 8, path = getwd())




###################################
########      Dsuite BC     #######
###################################

# make popmap for dsuite
dsuite.popmap <- popmap
dsuite.popmap <- popmap[popmap$pop %in% c("BC", "RUSS"),]
dsuite.popmap$dsuite_pop <- rep(NA, nrow(dsuite.popmap))

# add groupings for BC populations
dsuite.popmap[dsuite.popmap$lat > 57.52, "dsuite_pop"] <- "CLRU"
dsuite.popmap[dsuite.popmap$lat < 56.32, "dsuite_pop"] <- "CLGA"
dsuite.popmap[dsuite.popmap$locality == "South of Kinaskan Lake", "dsuite_pop"] <- "KIN_LAKE"
dsuite.popmap[dsuite.popmap$locality == "Cassiar Highway (Hwy 37), Burrage River Crossing", "dsuite_pop"] <- "BURRAGE"
dsuite.popmap[dsuite.popmap$lat == 57.1698, "dsuite_pop"] <- "SLATE_N"
dsuite.popmap[dsuite.popmap$lat == 57.15511667, "dsuite_pop"] <- "SLATE_S"
dsuite.popmap[dsuite.popmap$locality == "Iskut River Valley", "dsuite_pop"] <- "ISKUT"
dsuite.popmap[dsuite.popmap$locality == "Cassiar Highway (Hwy 37), Thomas Creek Ccrossing", "dsuite_pop"] <- "THOMAS"
dsuite.popmap[dsuite.popmap$locality == "Cassiar Highway (Hwy 37), Thomas Creek Crossing", "dsuite_pop"] <- "THOMAS"
dsuite.popmap[dsuite.popmap$locality == "Cassiar Highway (Hwy 37), Devil Creek Crossing", "dsuite_pop"] <- "DEVIL"
dsuite.popmap[dsuite.popmap$locality == "Cassiar Highway (Hwy 37), Junction with Gamma Creek", "dsuite_pop"] <- "GAMMA"
dsuite.popmap[dsuite.popmap$locality == "Cassiar Highway (Hwy 37), Nigunsaw Summit", "dsuite_pop"] <- "NIGUNSAW"
dsuite.popmap[dsuite.popmap$locality == "Cassiar Highway (Hwy 37), 5km N of Bell II Crossing", "dsuite_pop"] <- "BELL2"
dsuite.popmap[dsuite.popmap$locality == "S side Glacier creek", "dsuite_pop"] <- "GLACIER"
dsuite.popmap[dsuite.popmap$locality == "S side Dettaki Creek", "dsuite_pop"] <- "DETTAKI"
dsuite.popmap[dsuite.popmap$locality == "Stewart-Cassiar Highway, British Columbia Highway 37, Dease Lake Highway, 36.1 km N Meziadin Junction", "dsuite_pop"] <- "Mez36N"
dsuite.popmap[dsuite.popmap$locality == "Spruce Creek", "dsuite_pop"] <- "SPRUCE"
dsuite.popmap[dsuite.popmap$pop == "RUSS", "dsuite_pop"] <- "Outgroup"

# just two columns for input to dsuite
dsuite.popmap.all <- dsuite.popmap
dsuite.popmap <- dsuite.popmap[!is.na(dsuite.popmap$dsuite_pop), c("id", "dsuite_pop")]

# write out popmap and vcf
for (i in unique(dsuite.popmap$dsuite_pop)) {
  if (i %in% c("Outgroup", "CLRU", "CLGA")) {
    next
  }
  p <- dsuite.popmap[dsuite.popmap$dsuite_pop %in% c("Outgroup", "CLRU", "CLGA", i),]
  v <- min_mac(rbv_filtered[samples=p$id], min.mac = 1)
  
  print(i)
  print(v)
  
  vcfR::write.vcf(v, file = paste0("dsuite/dsuite.", i, ".vcf.gz"))
  write.table(p, file = paste0("dsuite/dsuite.", i, ".popmap.txt"), quote = F, append = F, col.names = F, row.names = F, sep = "\t")
}



# read in each output file
dstats <- data.frame(matrix(nrow = 0, ncol = 11))
for (i in unique(dsuite.popmap$dsuite_pop)) {
  if (i %in% c("Outgroup", "CLRU", "CLGA")) {
    next
  }
  o <- read.table(file = paste0("dsuite_output/", "_", i, "_BBAA.txt"), header = T, sep = "\t")
  o <- cbind(i, o)
  dstats <- rbind(o, dstats)
}
colnames(dstats) <- c("group", "P1",	"P2",	"P3",	"Dstatistic",	"Z-score",	"p-value",	"f4-ratio",	"BBAA",	"ABBA",	"BABA")

str(dstats)

# make D negative if CLGA is P2
for (i in 1:nrow(dstats)) {
  if (dstats[i,"P2"] == "CLGA") {
    dstats[i,"Dstatistic"] <- dstats[i,"Dstatistic"] * -1
  }
}


# add latitude 
dstats$lat <- rep(NA, nrow(dstats))
for (i in 1:nrow(dstats)) {
  g <- dstats[i,"group"]
  dstats[i, "lat"] <- mean(dsuite.popmap.all[dsuite.popmap.all$dsuite_pop == g, "lat"], na.rm = T)
}

# add 95% confidence intervals
# Dsuite paper:  Z = D/SE(D)
# So SE(D) = D/Z
# get standard error, multiply by 1.96, add and subtract from D
se <- abs(dstats$Dstatistic / dstats$`Z-score`)
dstats$CI.low <- dstats$Dstatistic - (1.96*se)
dstats$CI.high <- dstats$Dstatistic + (1.96*se)

# make plot
bc_d <- ggplot(dstats, aes(x=Dstatistic, y=lat)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbar(aes(xmin = CI.low, xmax = CI.high), width = 0.02, colour = "grey") +
  geom_point(cex = 2, alpha=1) +
  scale_x_continuous(n.breaks = 8) +
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  xlab("D statistic")+
  ylab("Latitude")+
  labs(title = ) +
  theme_classic()


ggplot2::ggsave(filename = "BC_D.pdf", bc_d, width = 4, height = 8, path = getwd())





#####################################
########      Dsuite SEAK     #######
#####################################


# make popmap for dsuite
dsuite.popmap <- popmap
dsuite.popmap <- popmap[popmap$pop %in% c("SEAK", "BC", "RUSS"),]
dsuite.popmap$dsuite_pop <- rep(NA, nrow(dsuite.popmap))

# add groupings for BC and SEAK populations
dsuite.popmap[dsuite.popmap$lat > 57.52 & dsuite.popmap$pop == "BC", "dsuite_pop"] <- "CLRU_BC"
dsuite.popmap[dsuite.popmap$locality == "Indian Point, north of LeConte Bay", "dsuite_pop"] <- "CLRU_SEAK"
dsuite.popmap[dsuite.popmap$pop == "SEAK" & dsuite.popmap$mt_species == "gapperi", "dsuite_pop"] <- "CLGA"
dsuite.popmap[dsuite.popmap$pop == "RUSS", "dsuite_pop"] <- "Outgroup"
dsuite.popmap[dsuite.popmap$locality == "Mallard Slough", "dsuite_pop"] <- "Stikine"
dsuite.popmap[dsuite.popmap$locality == "Mallard Slough, Petersburg Quad, Alexander Archipelago", "dsuite_pop"] <- "Stikine"
dsuite.popmap[dsuite.popmap$locality == "Stikine River", "dsuite_pop"] <- "Stikine"
dsuite.popmap[dsuite.popmap$locality == "Stikine River, Petersburg Quad, Alexander Archipelago, Farm Island", "dsuite_pop"] <- "Stikine"
dsuite.popmap[dsuite.popmap$locality == "Stikine River, Ellis cabin, Petersburg Quad, Alexander Archipelago, Farm Island", "dsuite_pop"] <- "Stikine"
dsuite.popmap[dsuite.popmap$locality == "Wrangell District, Stikine River and Andrew's Creek junction, Mount Rynda Cabin", "dsuite_pop"] <- "Stikine"
dsuite.popmap[dsuite.popmap$locality == "Stikine River, Figure Eight Lake", "dsuite_pop"] <- "Stikine"
dsuite.popmap[dsuite.popmap$locality == "Shake's Slough, north side of the Stikine River", "dsuite_pop"] <- "Stikine"
dsuite.popmap[dsuite.popmap$locality == "Rynda's Cabin, south side of Stikine River off Andrew Creek", "dsuite_pop"] <- "Stikine"
dsuite.popmap[dsuite.popmap$locality == "Rynda's Cabin, south side of Stikine River off of Andrew's Slough", "dsuite_pop"] <- "Stikine"
dsuite.popmap[dsuite.popmap$locality == "Rynda's Cabin, South side of the Stikine River off of Andrew's Slough", "dsuite_pop"] <- "Stikine"
dsuite.popmap[dsuite.popmap$locality == "Stikine River, Petersburg Quad, Alexander Archipelago, Limb Island", "dsuite_pop"] <- "Stikine"
dsuite.popmap[dsuite.popmap$locality == "Limb Island, Stikine River, Petersburg Quad, Alexander Archipelago, Limb Island", "dsuite_pop"] <- "Stikine"


# just two columns for input to dsuite
dsuite.popmap.all <- dsuite.popmap
dsuite.popmap <- dsuite.popmap[!is.na(dsuite.popmap$dsuite_pop), c("id", "dsuite_pop")]

# write out popmap and vcf
for (i in unique(dsuite.popmap$dsuite_pop)) {
  if (i %in% c("Outgroup", "Stikine", "CLGA")) {
    next
  }
  p <- dsuite.popmap[dsuite.popmap$dsuite_pop %in% c("Outgroup", "Stikine", "CLGA", i),]
  v <- min_mac(rbv_filtered[samples=p$id], min.mac = 1)
  
  print(i)
  print(v)
  
  vcfR::write.vcf(v, file = paste0("dsuite_SEAK/dsuite.", i, ".vcf.gz"))
  write.table(p, file = paste0("dsuite_SEAK/dsuite.", i, ".popmap.txt"), quote = F, append = F, col.names = F, row.names = F, sep = "\t")
}




# read in each output file
dstats_seak <- data.frame(matrix(nrow = 0, ncol = 11))
for (i in unique(dsuite.popmap$dsuite_pop)) {
  if (i %in% c("Outgroup", "Stikine", "CLGA")) {
    next
  }
  o <- read.table(file = paste0("dsuite_SEAK_output/", "_", i, "_BBAA.txt"), header = T, sep = "\t")
  o <- cbind(i, o)
  dstats_seak <- rbind(o, dstats_seak)
}
colnames(dstats_seak) <- c("group", "P1",	"P2",	"P3",	"Dstatistic",	"Z-score",	"p-value",	"f4-ratio",	"BBAA",	"ABBA",	"BABA")

str(dstats_seak)


# add 95% confidence intervals
# Dsuite paper:  Z = D/SE(D)
# So SE(D) = D/Z
# get standard error, multiply by 1.96, add and subtract from D
se <- abs(dstats_seak$Dstatistic / dstats_seak$`Z-score`)
dstats_seak$CI.low <- dstats_seak$Dstatistic - (1.96*se)
dstats_seak$CI.high <- dstats_seak$Dstatistic + (1.96*se)

# make plot
seak_d <- ggplot(dstats_seak[2,], aes(x=Dstatistic, y=group)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbar(aes(xmin = CI.low, xmax = CI.high), width = 0.02, colour = "grey") +
  geom_point(cex = 2, alpha=1) +
  scale_x_continuous(n.breaks = 8, limits = c(-0.6, 0.6)) +
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  xlab("D statistic")+
  ylab("")+
  labs(title = ) +
  theme_classic()


ggplot2::ggsave(filename = "SEAK_D.pdf", seak_d, width = 4, height = 1, path = getwd())






######################################
##   Distance matrix only gapperi
######################################
gapperi.popmap <- popmap[popmap$q1 > 0.8 & !is.na(popmap$q1),]
gapperi.vcfR <- min_mac(rbv_filtered[samples=gapperi.popmap$id], min.mac = 1)
gapperi.vcfR
dist.gapperi <- vcfR2splitstree(vcfR = gapperi.vcfR, popmap = gapperi.popmap, file.prefix = "dist.gapperi")
plot(nj(dist.gapperi),"u")

#####################################
##   Distance matrix only rutilus
#####################################
rutilus.popmap <- popmap[popmap$q1 < 0.7 & !is.na(popmap$q1),]
rutilus.vcfR <- min_mac(rbv_filtered[samples=rutilus.popmap$id], min.mac = 1)
rutilus.vcfR
dist.rutilus <- vcfR2splitstree(vcfR = rutilus.vcfR, popmap = rutilus.popmap, file.prefix = "dist.rutilus")
plot(nj(dist.rutilus),"u")


##############################
######   Depth 
##############################
depth <- extract.gt(rbv_filtered, element = "DP", as.numeric=TRUE)
depth[is.na(depth)] <- 0
identical(colnames(depth), popmap$id)
popmap$depth <- colMeans(depth)


#################################
####### PCA
#################################


# Make a genlight of subsamples
genlight.rbv <- vcfR2genlight(rbv_filtered)
# Assign populations
pop(genlight.rbv) <- popmap$pop
#perform PCA
pca.rbv<-glPca(genlight.rbv, nf = 10)
#pull pca scores out of df
pca.scores.rbv<-as.data.frame(pca.rbv$scores)
pca.scores.rbv$id <- popmap$id
pca.scores.rbv$pop <- popmap$pop
pca.scores.rbv$date <- popmap$date
pca.scores.rbv$locality <- popmap$locality
pca.scores.rbv$nuc_morph <- popmap$morph


pca.scores.rbv$nuc_species <- popmap$nuc_species
pca.scores.rbv$nuc_morph <- popmap$nuc_morph
pca.scores.rbv$mt_species <- popmap$mt_species
pca.scores.rbv$mt_morph <- popmap$mt_morph

pca.scores.rbv$depth <- popmap$depth
#calculate missingness by individual
miss.rbv<-colSums(is.na(rbv_filtered@gt))/nrow(rbv_filtered@gt)
#add missingness to df
pca.scores.rbv$missing<-miss.rbv[-1]
#record variance percentages explained
var.frac.rbv<- pca.rbv$eig/sum(pca.rbv$eig)*100
#ggplot color by species
ggplot(pca.scores.rbv, aes(x=PC1, y=PC2, color=pop)) +
  geom_point(cex = 2, alpha=1)+
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  xlab(paste("PC1, ", round(var.frac.rbv[1], digits = 2), "% variance explained", sep = ""))+
  ylab(paste("PC2, ", round(var.frac.rbv[2], digits = 2), "% variance explained", sep = ""))+
  labs(title = ) +
  theme_classic()

#ggplot color by species, PC3
ggplot(pca.scores.rbv, aes(x=PC1, y=PC3, color=pop)) +
  geom_point(cex = 2, alpha=1)+
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  xlab(paste("PC1, ", round(var.frac.rbv[1], digits = 2), "% variance explained", sep = ""))+
  ylab(paste("PC3, ", round(var.frac.rbv[3], digits = 2), "% variance explained", sep = ""))+
  labs(title = ) +
  theme_classic()

#ggplot color by species, PC4
ggplot(pca.scores.rbv, aes(x=PC1, y=PC4, color=pop)) +
  geom_point(cex = 2, alpha=1)+
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  xlab(paste("PC1, ", round(var.frac.rbv[1], digits = 2), "% variance explained", sep = ""))+
  ylab(paste("PC4, ", round(var.frac.rbv[4], digits = 2), "% variance explained", sep = ""))+
  labs(title = ) +
  theme_classic()

#ggplot color by missingness
ggplot(pca.scores.rbv, aes(x=PC1, y=PC2, color=missing)) +
  geom_point(cex = 2, alpha=1)+
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  xlab(paste("PC1, ", round(var.frac.rbv[1], digits = 2), "% variance explained", sep = ""))+
  ylab(paste("PC2, ", round(var.frac.rbv[2], digits = 2), "% variance explained", sep = ""))+
  labs(title = ) +
  theme_classic()

#ggplot color by missingness
ggplot(pca.scores.rbv, aes(x=PC1, y=PC4, color=missing)) +
  geom_point(cex = 2, alpha=1)+
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  xlab(paste("PC1, ", round(var.frac.rbv[1], digits = 2), "% variance explained", sep = ""))+
  ylab(paste("PC3, ", round(var.frac.rbv[3], digits = 2), "% variance explained", sep = ""))+
  labs(title = ) +
  theme_classic()

#ggplot color by depth
ggplot(pca.scores.rbv, aes(x=PC1, y=PC4, color=depth)) +
  geom_point(cex = 2, alpha=1)+
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  xlab(paste("PC1, ", round(var.frac.rbv[1], digits = 2), "% variance explained", sep = ""))+
  ylab(paste("PC3, ", round(var.frac.rbv[3], digits = 2), "% variance explained", sep = ""))+
  labs(title = ) +
  theme_classic()

#ggplot missing vs depth
ggplot(pca.scores.rbv, aes(x=depth, y=missing, color=pop)) +
  geom_point(cex = 2, alpha=1)+
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  #xlab(paste("PC1, ", round(var.frac.rbv[1], digits = 2), "% variance explained", sep = ""))+
  #ylab(paste("PC3, ", round(var.frac.rbv[3], digits = 2), "% variance explained", sep = ""))+
  labs(title = ) +
  geom_segment()
theme_classic()

#ggplot missing vs depth
ggplot(pca.scores.rbv, aes(x=depth, y=missing, color=pop)) +
  geom_point(cex = 2, alpha=1)+
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  #xlab(paste("PC1, ", round(var.frac.rbv[1], digits = 2), "% variance explained", sep = ""))+
  #ylab(paste("PC3, ", round(var.frac.rbv[3], digits = 2), "% variance explained", sep = ""))+
  labs(title = ) +
  scale_x_continuous(limits = c(0,100)) +
  geom_segment(aes(x = 8.8, xend = 8.8, y = 0, yend = 0.8), color = "black") +
  theme_classic()


#ggplot color by species
ggplot(pca.scores.rbv, aes(x=PC1, y=PC2, color=nuc_species)) +
  geom_point(cex = 2, alpha=1)+
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  xlab(paste("PC1, ", round(var.frac.rbv[1], digits = 2), "% variance explained", sep = ""))+
  ylab(paste("PC2, ", round(var.frac.rbv[2], digits = 2), "% variance explained", sep = ""))+
  labs(title = ) +
  theme_classic()

#ggplot color by date
ggplot(pca.scores.rbv, aes(x=PC1, y=PC2, color=date)) +
  geom_point(cex = 2, alpha=1)+
  guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
  xlab(paste("PC1, ", round(var.frac.rbv[1], digits = 2), "% variance explained", sep = ""))+
  ylab(paste("PC2, ", round(var.frac.rbv[2], digits = 2), "% variance explained", sep = ""))+
  labs(title = ) +
  theme_classic()









##################################################
#  WRITE OUT GEOGRAPHICALLY SEPERATED VCF FILES  #
##################################################


#vcfR::write.vcf(vcfR.seak, "rbv.seak.vcf.gz")
#vcfR::write.vcf(vcfR.bc, "rbv.bc.vcf.gz")
#vcfR::write.vcf(vcfR.nwt, "rbv.nwt.vcf.gz")
#write.table(popmap.seak, file = "popmap.seak.txt", quote = F, sep = "\t", col.names = T, row.names = F)
#write.table(popmap.bc, file = "popmap.bc.txt", quote = F, sep = "\t", col.names = T, row.names = F)
#write.table(popmap.nwt, file = "popmap.nwt.txt", quote = F, sep = "\t", col.names = T, row.names = F)







