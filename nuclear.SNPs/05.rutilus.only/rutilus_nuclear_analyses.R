setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/RBVs/09.5plates.phylogeography/nuclear.SNPs/05.rutilus.only")


#library("tinytex")
#library("poppr")
library("hierfstat")
#library("mmod")
library(vcfR)
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
library(rEEMSplots)
library(rgdal)
library(rworldmap)
library(rworldxtra)
library(ggplot2)
library(dplyr)

###############################
####    Define funcitons    ###
###############################

vcf2str <- function(vcfR = NULL, popmap = NULL, filename = NULL, missing.character = -9) {
  
  # extract genotype matrix
  m <- extract.gt(vcfR)
  m <- t(m)
  
  # recode missing data
  m[is.na(m)] <- (paste0(missing.character, "/", missing.character))
  
  # get unique values in matrix
  snp.values <- names(table(m))
  if (!all(snp.values %in% c("-9/-9", "0/0", "0/1", "1/1"))) {
    stop("unexpected values in genotype matrix")
  }
  
  # put popmap in same order as rownames in genotype matrix
  pm <- popmap[order(match(popmap$id, colnames(m))),]
  if (!identical(rownames(m), pm$id)) {
    stop("the individuals in the popmap and genotype matrix cannot be put in the same order... are some missing?")
  }
  
  # add id and pop columns to matrix
  id <- paste0(pm$id, "/", pm$id)
  pop <- paste0(pm$pop, "/", pm$pop)
  int.pop <- pop
  # convert pop to integer
  int <- 2
  for (p in unique(pop)) {
    int.pop[int.pop==p] <- paste0(int, "/", int)
    int <- int + 1
  }
  
  m <- cbind(id, int.pop, m)
  
  # make new matrix in structure format
  n <- matrix(nrow = 0, ncol = ncol(m))
  colnames(n) <- colnames(m)
  
  # populate new matrix
  for(i in 1:nrow(m)) {
    allele1 <- sapply(strsplit(m[i,], "/"), '[', 1)
    n <- rbind(n, allele1)
    
    allele2 <- sapply(strsplit(m[i,], "/"), '[', 2)
    n <- rbind(n, allele2)
  }
  
  # rename columns for structure
  colnames(n) <- c(" ", " ", colnames(n)[-1:-2])
  
  # write file
  write.table(n, file = filename, append = F, quote = F, sep = "\t", col.names = T, row.names = F)
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
rutilus_filtered <- read.vcfR("rutilus.nuc.vcf.gz")

# read in popmap with nuc and mito calls
detailed.popmap <- read_delim("../02.nuc.analyses/no.russ.popmap.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

# keep only rutilus
popmap <- detailed.popmap[detailed.popmap$nuc_species=="rutilus",]


# make sure popmap is in same order as ids in the vcf
identical(colnames(extract.gt(rutilus_filtered)), popmap$id)



######################################
##    Write popmap for structure    ##
######################################

# run structure outside of R
vcf2str(vcfR = rutilus_filtered, popmap = popmap, filename = "rutilus.nuc.str")
rutilus_filtered







##################################
#     Plot structure results     #
##################################

for (file in list.files("01.output.1M.K2_8/", pattern = "*_f", full.names = T)) {
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
slist <- readQ(list.files("01.output.1M.K2_8/", pattern = "*_f", full.names = T))
tabQ <- tabulateQ(slist)
sumQ <- summariseQ(tabQ)
evanno.data <- evannoMethodStructure(sumQ)
evanno.plot <- evannoMethodStructure(sumQ, returnplot = T, exportplot = F, returndata = F, basesize=12, linesize=0.7)
grid.arrange(evanno.plot)


# view structure results labeled by pop
betterStructurePlot(K3_rep1_structure$q.mat, horiz = F, label.pops = T, legend.position = "top", col = color_ramp(3), sort.probs = T, type = "bar")

# view structure results labeled by individual
betterStructurePlot(K6_rep9_structure$q.mat, horiz = F, legend.position = "top", col = color_ramp(6), sort.probs = T, type = "bar")


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
  ggplot2::ggsave(filename = paste0("all_structure_plots/", k, "_plot_ind_labels.png"), plot, width = 15, height = 5, path = getwd())
  # plot with pop labels
  plot <- betterStructurePlot(kdata$q.mat, label.pops = T, horiz = F, legend.position = "top", col = color_ramp(num), sort.probs = T, type = "bar")
  ggplot2::ggsave(filename = paste0("all_structure_plots/", k, "_plot_pops_labels.png"), plot, width = 15, height = 5, path = getwd())
}





#############################################
##   sort structure by order on iqtree   ####
#############################################

# get q matrix and melt into dataframe for plotting
q <- K3_rep1_structure$q.mat
df <- melt(q, id.vars = c(1:3), variable.name = "Group", value.name = "ancestry")

# get order of inds on tree
tree.order <- read.table(file = "../04.iqtree/rbv.tree.order.txt")[,1]
# order of only rutilus
tree.order <- tree.order[-1:-2]
tree.order <- tree.order[tree.order %in% df$id]

# sort ids by order in iqtree
# make sure popmap is in same order as ids in the vcf
identical(df$id, tree.order)
df <- df[order(match(df$id, tree.order)),]
identical(unique(df$id), tree.order)


# in order to add black outlines to each bar
pop.div <- 1:nrow(df) - 0.5

# plot
structure <- ggplot(df, aes(x=ancestry, y=id, fill = Group)) +
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
  scale_fill_manual(values = c("#FC9272", "#A50F15", "#EF3B2C"))

ggplot2::ggsave(filename = "rutilus.alone.structure.plot.pdf", structure, width = 2, height = 8, path = getwd())





#################################
##   Make popmap for mapping   ##
#################################



rutilus.q.mat <- K3_rep1_structure$q.mat
identical(popmap$id, rutilus.q.mat$id)

popmap$rutilus.q1 <- rutilus.q.mat$Group.1
popmap$rutilus.q2 <- rutilus.q.mat$Group.2
popmap$rutilus.q3 <- rutilus.q.mat$Group.3

write.table(popmap, file = "rutilus.nuc.map.txt", sep = "\t", quote = F, col.names = T, row.names = F)



###################################
##   Missing data and q values   ##
###################################

# make sure popmap is in same order as ids in the vcf
identical(colnames(extract.gt(rutilus_filtered)), popmap$id)
# calculate missing data for rutilus-specific dataset
popmap$rut_nuc_missing_data <- (colSums(is.na(rutilus_filtered@gt))/nrow(rutilus_filtered@gt))[-1]

# Q1
summary(lm(rutilus.q1 ~ rut_nuc_missing_data, popmap, subset = (pop == "SEAK")))
summary(lm(rutilus.q1 ~ rut_nuc_missing_data, popmap, subset = (pop == "BC")))
summary(lm(rutilus.q1 ~ rut_nuc_missing_data, popmap, subset = (pop == "NWT")))

q1_r <- c(summary(lm(rutilus.q1 ~ rut_nuc_missing_data, popmap))$r.squared,
          summary(lm(rutilus.q1 ~ pop, popmap))$r.squared,
          summary(lm(rutilus.q1 ~ rut_nuc_missing_data + pop, popmap))$r.squared)

q1_miss <- ggplot(popmap, aes(x=rut_nuc_missing_data, y=rutilus.q1, color=pop)) +
  geom_point(cex = 3, alpha = 1, shape=1) +
  #geom_smooth(method = "lm", se = FALSE) +
  #annotate("text", x=0.5, y=1, label= "R2=0.012") +
  #annotate("text", x=0.5, y=0.95, label= "p=0.2933") +
  xlab(paste("Proportion nuclear missing data")) +
  ylab(paste("Q1 from structure (C. rutilus only)")) +
  labs(title = "")

ggplot2::ggsave(filename = "rut_missing_q1.png", q1_miss, width = 7, height = 4, path = getwd())

# Q2
summary(lm(rutilus.q2 ~ rut_nuc_missing_data, popmap, subset = (pop == "SEAK")))
summary(lm(rutilus.q2 ~ rut_nuc_missing_data, popmap, subset = (pop == "BC")))
summary(lm(rutilus.q2 ~ rut_nuc_missing_data, popmap, subset = (pop == "NWT")))

q2_r <- c(summary(lm(rutilus.q2 ~ rut_nuc_missing_data, popmap))$r.squared,
          summary(lm(rutilus.q2 ~ pop, popmap))$r.squared,
          summary(lm(rutilus.q2 ~ rut_nuc_missing_data + pop, popmap))$r.squared)


q2_miss <- ggplot(popmap, aes(x=rut_nuc_missing_data, y=rutilus.q2, color=pop)) +
  geom_point(cex = 3, alpha = 1, shape=1) +
  #geom_smooth(method = "lm", se = FALSE) +
  #annotate("text", x=0.48, y=1, label= "R2=6.899e-05") +
  #annotate("text", x=0.48, y=0.95, label= "p=0.9367") +
  xlab(paste("Proportion nuclear missing data")) +
  ylab(paste("Q2 from structure (C. rutilus only)")) +
  labs(title = "")

ggplot2::ggsave(filename = "rut_missing_q2.png", q2_miss, width = 7, height = 4, path = getwd())

#Q3
summary(lm(rutilus.q3 ~ rut_nuc_missing_data, popmap, subset = (pop == "SEAK")))
summary(lm(rutilus.q3 ~ rut_nuc_missing_data, popmap, subset = (pop == "BC")))
summary(lm(rutilus.q3 ~ rut_nuc_missing_data, popmap, subset = (pop == "NWT")))

q3_r <- c(summary(lm(rutilus.q3 ~ rut_nuc_missing_data, popmap))$r.squared,
          summary(lm(rutilus.q3 ~ pop, popmap))$r.squared,
          summary(lm(rutilus.q3 ~ rut_nuc_missing_data + pop, popmap))$r.squared)

q3_miss <- ggplot(popmap, aes(x=rut_nuc_missing_data, y=rutilus.q3, color=pop)) +
  geom_point(cex = 3, alpha = 1, shape=1) +
  #geom_smooth(method = "lm", se = FALSE) +
  #annotate("text", x=0.5, y=1, label= "R2=0.0094") +
  #annotate("text", x=0.5, y=0.95, label= "p=0.3533") +
  xlab(paste("Proportion nuclear missing data")) +
  ylab(paste("Q3 from structure (C. rutilus only)")) +
  labs(title = "")

ggplot2::ggsave(filename = "rut_missing_q3.png", q3_miss, width = 7, height = 4, path = getwd())

rsquared <- data.frame(rbind(q1_r, q2_r, q3_r))
colnames(rsquared) <- c("missing_data", "geography", "geography+missing_data")
rsquared$improvement <- rsquared$`geography+missing_data` - rsquared$geography
rsquared <- round(rsquared, digits = 6)

write.table(rsquared, file = "variance_from_missingness.txt", quote = F, sep = "\t")

##############################################
###    SET UP: EEMS only SEAK BC rutilus   ###
##############################################

# make sure popmap is in same order as ids in the vcf
identical(colnames(extract.gt(rutilus_filtered)), popmap$id)

# just SEAK and BC
eems.popmap <- popmap[popmap$pop %in% c("SEAK", "BC"),]

# get rid of sample with missing coordinates
eems.popmap <- eems.popmap[!is.na(eems.popmap$lat),]

# make new vcfR
eems.vcfR <- min_mac(rutilus_filtered[samples=eems.popmap$id], min.mac = 1)
eems.vcfR

# make sure popmap is in same order as ids in the vcf
identical(colnames(extract.gt(eems.vcfR)), eems.popmap$id)

# Write out coordinates for each sample, tab separated
write.table(eems.popmap[,c("lat", "long")], "datapath.rutilus.coord", col.names = F, row.names = F, quote = F, sep='\t')


###### make genetic input file


#define function to convert a genotype matrix into a pairwise difference matrix
bed2diffs_v2 <- function(genotypes) {
  nIndiv <- nrow(genotypes)
  nSites <- ncol(genotypes)
  missing <- is.na(genotypes)
  ## Impute NAs with the column means (= twice the allele frequencies)
  geno_means <- colMeans(genotypes, na.rm = TRUE)
  # nIndiv rows of genotype means
  geno_means <- matrix(geno_means, nrow = nIndiv, ncol = nSites, byrow = TRUE)
  ## Set the means which correspond to observed genotypes to 0
  geno_means[missing == FALSE] <- 0
  ## Set the missing genotypes to 0 (used to be NA)
  genotypes[missing == TRUE] <- 0
  genotypes <- genotypes + geno_means
  similarities <- genotypes %*% t(genotypes) / nSites
  self_similarities <- diag(similarities)
  vector1s <- rep(1, nIndiv)
  diffs <- self_similarities %*% t(vector1s) + vector1s %*% t(self_similarities) - 2 * similarities
  diffs
}

gl<-vcfR2genlight(eems.vcfR) 
geno.mat<-as.matrix(gl)
diff.mat<-bed2diffs_v2(genotypes = geno.mat)
#make sure order matches
identical(rownames(geno.mat), eems.popmap$id)

write.table(diff.mat, "datapath.rutilus.diffs", col.names = F, row.names = F, quote = F, sep='\t')

#INFO:
eems.vcfR



########################
###   rutilus plots  ###
########################

mcmcpath = "rep100.rutilus/"
plotpath = "./rep100.rutilus/rep100_rutilus"

projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"
#projection_mercator <- "+init=epsg:3338"

coords <- eems.popmap[,c("long", "lat")]
#colors <- c("red", "green", "blue", "purple", "orange")
#labels <- LETTERS[1:5]
coords_merc <- sp::spTransform(
  SpatialPoints(coords, CRS(projection_none)),
  CRS(projection_mercator)
)
# `coords_merc` is a SpatialPoints structure
# but we only need the coordinates themselves
coords_merc <- coords_merc@coords
# Coordinates in latitude/longitude
coords
# The same coordinates projected
coords_merc



eems.plots(
  mcmcpath = mcmcpath,
  plotpath = plotpath,
  longlat = FALSE,
  projection.in = projection_none,
  projection.out = projection_mercator,
  add.map = TRUE,
  col.map = "black",
  lwd.map = 1,
  out.png = FALSE,
  m.plot.xy = {
    points(coords_merc, col = "black", pch = 16)
  },
  q.plot.xy = {
    points(coords_merc, col = "black", pch = 16)
  }
  
)


























