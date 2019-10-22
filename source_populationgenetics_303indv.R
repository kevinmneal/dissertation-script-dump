#### source file to load packages, files, objects, and custom functions
#### run this before running other analysis scripts!

##### required packages
library(adegenet)
library(SNPRelate)
library(vcfR)
library(poppr)
library(PopGenome)
library(cluster)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ape)
library(vegan)
library(adespatial)
library(mmod)
library(RColorBrewer)
library(igraph)
library(gplots)
library(viridis)
library(tess3r) # installed from github; make sure it's up-to-date
library(LEA) # installed from github; make sure it's up-to-date
library(rworldmap)
library(maps)
library(gstudio)
library(hierfstat)
library(pegas)
library(radiator)
library(unPC)
library(rhierbaps)
library(pophelper)
library(popgraph)
library(gstudio)
library(gridExtra)
library(latticeExtra)
library(pcadapt) # installed from github; make sure it's up-to-date
library(DataCombine)
library(StAMPP)
library(dartR)
library(gdistance)
library(ResistanceGA) # installed from github; make sure it's up-to-date





parentdir.303indiv <- "G:/My Drive/Illumina Sequencing Data/20181004_ALLsocal_303indiv_highqual"
setwd(parentdir.303indiv)
vcfpaths.303indiv <- list.files(parentdir.303indiv, pattern=".vcf$", recursive=T)

write.csv(popcoords.303, "popcoords.303.csv")
popcoords.303 <- read.csv("G:/My Drive/Illumina Sequencing Data/20181004_ALLsocal_303indiv_highqual/allsocal_303indiv_popcoords_multilevelpops_betternames_20181017.csv",
                          header=T,
                          stringsAsFactors = F)

#colnames(vcf.50pctmsng@gt)[-1] == popcoords.303$Ind

vcf.10pctmsng <- read.vcfR(vcfpaths.303indiv[8]) # 1378 SNPs
genind.10pctmsng <- vcfR2genind(vcf.10pctmsng)
indNames(genind.10pctmsng) <- popcoords.303$BetterName
pop(genind.10pctmsng) <- popcoords.303$Pop
genind.10pctmsng@other$xy <- popcoords.303[,c("Longitude", "Latitude")]
genind.10pctmsng.seppops <- seppop(genind.10pctmsng)
genpop.10pctmsng <- genind2genpop(genind.10pctmsng)
gl.10pctmsng <- vcfR2genlight(vcf.10pctmsng)
ploidy(gl.10pctmsng) <- 2
pop(gl.10pctmsng) <- popcoords.303$Pop
gl.10pctmsng@other$xy <- popcoords.303[,c("Longitude", "Latitude")]



vcf.20pctmsng <- read.vcfR(vcfpaths.303indiv[2])  # 2135 SNPs
genind.20pctmsng <- vcfR2genind(vcf.20pctmsng)
indNames(genind.20pctmsng) <- popcoords.303$BetterName
pop(genind.20pctmsng) <- popcoords.303$Pop
genind.20pctmsng@other$xy <- popcoords.303[,c("Longitude", "Latitude")]
genind.20pctmsng.seppops <- seppop(genind.20pctmsng)
genpop.20pctmsng <- genind2genpop(genind.20pctmsng)
gl.20pctmsng <- vcfR2genlight(vcf.20pctmsng)
ploidy(gl.20pctmsng) <- 2
pop(gl.20pctmsng) <- popcoords.303$Pop
gl.20pctmsng@other$xy <- popcoords.303[,c("Longitude", "Latitude")]


vcf.50pctmsng <- read.vcfR(vcfpaths.303indiv[3]) # 3910 SNPs
genind.50pctmsng <- vcfR2genind(vcf.50pctmsng)
indNames(genind.50pctmsng) <- popcoords.303$BetterName
pop(genind.50pctmsng) <- popcoords.303$Pop
genind.50pctmsng@other$xy <- popcoords.303[,c("Longitude", "Latitude")]
genind.50pctmsng.seppops <- seppop(genind.50pctmsng)
genpop.50pctmsng <- genind2genpop(genind.50pctmsng)
gl.50pctmsng <- vcfR2genlight(vcf.50pctmsng)
ploidy(gl.50pctmsng) <- 2
pop(gl.50pctmsng) <- popcoords.303$Pop
gl.50pctmsng@other$xy <- popcoords.303[,c("Longitude", "Latitude")]

# 75 pct missing data
vcf.75pctmsng <- read.vcfR(vcfpaths.303indiv[8]) # 6056 SNPs; excludes singletons and private doubletons
# 8894 if we only exclude singletons:
vcf.75pctmsng <- read.vcfR("G:/My Drive/Illumina Sequencing Data/20181004_ALLsocal_303indiv_highqual/nopopfile_75pctmsng/vcftools_allsocal303indiv_75pctmsng_nosingl_20181126.recode.vcf")
genind.75pctmsng <- vcfR2genind(vcf.75pctmsng)
indNames(genind.75pctmsng) <- popcoords.303$IndPop2
pop(genind.75pctmsng) <- popcoords.303$PopLevel3geogenclusters
genind.75pctmsng@other$xy <- popcoords.303[,c("Longitude3", "Latitude3")]
genind.75pctmsng.seppops <- seppop(genind.75pctmsng)
genpop.75pctmsng <- genind2genpop(genind.75pctmsng)
gl.75pctmsng <- vcfR2genlight(vcf.75pctmsng)
ploidy(gl.75pctmsng) <- 2
pop(gl.75pctmsng) <- popcoords.303$PopLevel3geogenclusters
gl.75pctmsng@other$xy <- popcoords.303[,c("Longitude3", "Latitude3")]


vcf.SNPinallpops.10pctmsng <- read.vcfR(vcfpaths.303indiv[4]) # 761 SNPs
genind.SNPinallpops.10pctmsng <- vcfR2genind(vcf.SNPinallpops.10pctmsng)
indNames(genind.SNPinallpops.10pctmsng) <- popcoords.303$BetterName
pop(genind.SNPinallpops.10pctmsng) <- popcoords.303$Pop
genind.SNPinallpops.10pctmsng@other$xy <- popcoords.303[,c("Longitude", "Latitude")]
genind.SNPinallpops.10pctmsng.seppops <- seppop(genind.SNPinallpops.10pctmsng)
genpop.SNPinallpops.10pctmsng <- genind2genpop(genind.SNPinallpops.10pctmsng)
gl.SNPinallpops.10pctmsng <- vcfR2genlight(vcf.SNPinallpops.10pctmsng)
ploidy(gl.SNPinallpops.10pctmsng) <- 2
pop(gl.SNPinallpops.10pctmsng) <- popcoords.303$Pop
gl.SNPinallpops.10pctmsng@other$xy <- popcoords.303[,c("Longitude", "Latitude")]


vcf.SNPinallpops <- read.vcfR(vcfpaths.303indiv[5]) # 928 SNPs
genind.SNPinallpops <- vcfR2genind(vcf.SNPinallpops)
indNames(genind.SNPinallpops) <- popcoords.303$BetterName
pop(genind.SNPinallpops) <- popcoords.303$Pop
genind.SNPinallpops@other$xy <- popcoords.303[,c("Longitude", "Latitude")]
genind.SNPinallpops.seppops <- seppop(genind.SNPinallpops)
genpop.SNPinallpops <- genind2genpop(genind.SNPinallpops)
gl.SNPinallpops <- vcfR2genlight(vcf.SNPinallpops)
ploidy(gl.SNPinallpops) <- 2
pop(gl.SNPinallpops) <- popcoords.303$Pop
gl.SNPinallpops@other$xy <- popcoords.303[,c("Longitude", "Latitude")]


vcf.nomsng <- read.vcfR(vcfpaths.303indiv[6]) # 131 SNPs
genind.nomsng <- vcfR2genind(vcf.nomsng)
indNames(genind.nomsng) <- popcoords.303$BetterName
pop(genind.nomsng) <- popcoords.303$Pop
genind.nomsng@other$xy <- popcoords.303[,c("Longitude", "Latitude")]
genind.nomsng.seppops <- seppop(genind.nomsng)
genpop.nomsng <- genind2genpop(genind.nomsng)
gl.nomsng <- vcfR2genlight(vcf.nomsng)
ploidy(gl.nomsng) <- 2
pop(gl.nomsng) <- popcoords.303$Pop
gl.nomsng@other$xy <- popcoords.303[,c("Longitude", "Latitude")]

popcoords.303 <- read.csv("G:/My Drive/Illumina Sequencing Data/20181004_ALLsocal_303indiv_highqual/allsocal_303indiv_popcoords_multiplelevels_betternames.csv",
                          header=T,
                          stringsAsFactors = F)

#### 303 individuals, 50 percent and 20 percent missing data

popcoords <- popcoords.303
vcf.50pctmsng <- read.vcfR("G:/My Drive/Illumina Sequencing Data/20181004_ALLsocal_303indiv_highqual/nopopfile_50pctmsng/vcftools_allsocal_303indv_highqual_50pctmsng_oneSNP_nosingletons_20181004.recode.vcf") # 3498 SNPs
genind.50pctmsng <- vcfR2genind(vcf.50pctmsng) # 3498 SNPs
indNames(genind.50pctmsng) <- popcoords$IndPop2
pop(genind.50pctmsng) <- popcoords$PopLevel2_prox
genind.50pctmsng@other$xy <- popcoords[,c("Longitude2", "Latitude2")]
genind.50pctmsng.seppops <- seppop(genind.50pctmsng)
genpop.50pctmsng <- genind2genpop(genind.50pctmsng)
gl.50pctmsng <- vcfR2genlight(vcf.50pctmsng)
ploidy(gl.50pctmsng) <- 2
pop(gl.50pctmsng) <- popcoords$PopLevel2_prox
gl.50pctmsng@other$xy <- popcoords[,c("Longitude2", "Latitude2")]



vcf.20pctmsng <- read.vcfR("G:/My Drive/Illumina Sequencing Data/20181004_ALLsocal_303indiv_highqual/nopopfile_20pctmsng/vcftools_allsocal_303indv_highqual_20pctmsng_oneSNP_nosingletons_20181004.recode.vcf")
genind.20pctmsng <- vcfR2genind(vcf.20pctmsng) # 1987 SNPs
indNames(genind.20pctmsng) <- popcoords$IndPop2
pop(genind.20pctmsng) <- popcoords$PopLevel2_prox
genind.20pctmsng@other$xy <- popcoords[,c("Longitude2", "Latitude2")]
genind.20pctmsng.seppops <- seppop(genind.20pctmsng)
genpop.20pctmsng <- genind2genpop(genind.20pctmsng)
gl.20pctmsng <- vcfR2genlight(vcf.20pctmsng)
ploidy(gl.20pctmsng) <- 2
pop(gl.20pctmsng) <- popcoords$PopLevel2_prox
gl.20pctmsng@other$xy <- popcoords[,c("Longitude2", "Latitude2")]


palettelist <- list(
  "shiny"=c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E"),
  "strong"=c("#11A4C8","#63C2C5","#1D4F9F","#0C516D","#2A2771","#396D35","#80C342","#725DA8","#B62025","#ED2224","#ED1943","#ED3995","#7E277C","#F7EC16","#F8941E","#8C2A1C","#808080"),
  "oceanfive"=c("#00A0B0", "#6A4A3C", "#CC333F", "#EB6841", "#EDC951"),
  "keeled"=c("#48B098", "#91CB62", "#FFEE3B", "#FB9013", "#FF3C28"),
  "vintage"=c("#400F13", "#027368", "#A3BF3F", "#F2B950", "#D93A2B"),
  "muted"=c("#46BDDD","#82DDCE","#F5F06A","#F5CC6A","#F57E6A"),
  "teal"=c("#CFF09E","#A8DBA8","#79BD9A","#3B8686","#0B486B"),
  "merry"=c("#5BC0EB","#FDE74C","#9BC53D","#E55934","#FA7921"),
  "funky"=c("#A6CEE3", "#3F8EAA", "#79C360", "#E52829", "#FDB762","#ED8F47","#9471B4"),
  "retro"=c("#01948E","#A9C4E2","#E23560","#01A7B3","#FDA963","#323665","#EC687D"),
  "cb_paired"=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),
  "cb_set3"=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"),
  "morris"=c("#4D94CC","#34648A","#8B658A","#9ACD32","#CC95CC","#9ACD32","#8B3A39","#CD6601","#CC5C5B","#8A4500"),
  "wong"=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#006699","#D55E00","#CC79A7"),
  "krzywinski"=c("#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE"))
# probably stick to shiny
combinedpalettelist <- unique(unlist(c(palettelist)))
palettesub <- sample(combinedpalettelist, size=50)


plot(c(1:50), col=palettesub, cex=2, pch=16)

distToDataframe <- function(inDist) {
  if (class(inDist) != "dist") stop("wrong input type")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  data.frame(
    row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
    col = rep(B[-length(B)], (length(B)-1):1),
    value = as.vector(inDist))
}

# r2.high2 <- FindReplace(dat.sort2,
#                         Var="INDV2",
#                         replaceData=popcoords.225,
#                         from="Ind",
#                         to="Pop")




