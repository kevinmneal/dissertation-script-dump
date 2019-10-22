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





parentdir.303indiv <- "G:/My Drive/Illumina Sequencing Data/20181128_allsocal303indiv_finaldatasets/"
setwd(parentdir.303indiv)
vcfpaths.303indiv <- list.files(parentdir.303indiv, pattern=".vcf$", recursive=T)

write.csv(popcoords.303, "popcoords.303.csv")
popcoords.303 <- read.csv("G:/My Drive/Illumina Sequencing Data/20181128_allsocal303indiv_finaldatasets/allsocal_303indiv_popcoords_multilevelpops_betternames_20181017.csv",
                          header=T,
                          stringsAsFactors = F)

popcoords.alloc <- read.csv("G:/My Drive/Illumina Sequencing Data/20181128_allsocal303indiv_finaldatasets/allOC_286indiv_popcoords_multilevelpops_betternames_20181211.csv",
                            header=T,
                            stringsAsFactors = F)
popcoords.artif <- popcoords.303[popcoords.303$PopLevel7K3=="ARTIF",]
popcoords.coastal <- popcoords.303[popcoords.303$PopLevel7K3=="COAST",]
popcoords.inland <- popcoords.303[popcoords.303$PopLevel7K3=="INLAND",]

#colnames(vcf.75pctmsng.artif@gt)[-1] == popcoords.artif$Indfile


# 75 pct missing data - allOC (so both native and artif. pops in Orange County only)
vcf.75pctmsng.allOC <- read.vcfR("75pctmsng_nosingletons/vcftools_allOC286indiv_75pctmsng_oneSNPnosngl_20181127.recode.vcf") # 6056 SNPs; excludes singletons and private doubletons
# 5961 if we only exclude singletons
genind.75pctmsng.allOC <- vcfR2genind(vcf.75pctmsng.allOC)
indNames(genind.75pctmsng.allOC) <- popcoords.alloc$IndPop2
pop(genind.75pctmsng.allOC) <- popcoords.alloc$PopLevel3geogenclusters
genind.75pctmsng.allOC@other$xy <- popcoords.alloc[,c("LongitudePopLevel3", "LatitudePopLevel3")]
genind.75pctmsng.allOC.seppops <- seppop(genind.75pctmsng.allOC)
genpop.75pctmsng.allOC <- genind2genpop(genind.75pctmsng.allOC)
gl.75pctmsng.allOC <- vcfR2genlight(vcf.75pctmsng.allOC)
ploidy(gl.75pctmsng.allOC) <- 2
indNames(gl.75pctmsng.allOC) <- popcoords.alloc$IndPop2
pop(gl.75pctmsng.allOC) <- popcoords.alloc$PopLevel7K3
gl.75pctmsng.allOC@other$xy <- popcoords.alloc[,c("LongitudePopLevel3", "LatitudePopLevel3")]

# 75 pct missing data + HWE filtering - allOC (so both native and artif. pops in Orange County only)
vcf.75pctmsng.allOChwe <- read.vcfR("75pctmsng_nosingletons/vcftools_allOC_286indiv_pophwept001_75pctmsng_oneSNPnosngl_20181212.recode.vcf") # 6056 SNPs; excludes singletons and private doubletons
# 5339 SNPs if we only exclude singletons
genind.75pctmsng.allOChwe <- vcfR2genind(vcf.75pctmsng.allOChwe)
indNames(genind.75pctmsng.allOChwe) <- popcoords.alloc$IndPop2
pop(genind.75pctmsng.allOChwe) <- popcoords.alloc$PopLevel3geogenclusters
genind.75pctmsng.allOChwe@other$xy <- popcoords.alloc[,c("LongitudePopLevel3", "LatitudePopLevel3")]
genind.75pctmsng.allOChwe.seppops <- seppop(genind.75pctmsng.allOChwe)
genpop.75pctmsng.allOChwe <- genind2genpop(genind.75pctmsng.allOChwe)
gl.75pctmsng.allOChwe <- vcfR2genlight(vcf.75pctmsng.allOChwe)
ploidy(gl.75pctmsng.allOChwe) <- 2
indNames(gl.75pctmsng.allOChwe) <- popcoords.alloc$IndPop2
pop(gl.75pctmsng.allOChwe) <- popcoords.alloc$PopLevel7K3
gl.75pctmsng.allOChwe@other$xy <- popcoords.alloc[,c("LongitudePopLevel3", "LatitudePopLevel3")]


# 75 pct missing data + HWE filtering - allsocal
vcf.75pctmsng.allsocalhwe <- read.vcfR("75pctmsng_nosingletons/vcftools_allsocal_303indiv_pophwept001_75pctmsng_oneSNPnosngl_20181212.recode.vcf") # 6056 SNPs; excludes singletons and private doubletons
# 5946 SNPs if we only exclude singletons
genind.75pctmsng.allsocalhwe <- vcfR2genind(vcf.75pctmsng.allsocalhwe)
indNames(genind.75pctmsng.allsocalhwe) <- popcoords.303$IndPop2
pop(genind.75pctmsng.allsocalhwe) <- popcoords.303$PopLevel3geogenclusters
genind.75pctmsng.allsocalhwe@other$xy <- popcoords.303[,c("LongitudePopLevel3", "LatitudePopLevel3")]
genind.75pctmsng.allsocalhwe.seppops <- seppop(genind.75pctmsng.allsocalhwe)
genpop.75pctmsng.allsocalhwe <- genind2genpop(genind.75pctmsng.allsocalhwe)
gl.75pctmsng.allsocalhwe <- vcfR2genlight(vcf.75pctmsng.allsocalhwe)
ploidy(gl.75pctmsng.allsocalhwe) <- 2
indNames(gl.75pctmsng.allsocalhwe) <- popcoords.303$IndPop2
pop(gl.75pctmsng.allsocalhwe) <- popcoords.303$PopLevel7K4
gl.75pctmsng.allsocalhwe@other$xy <- popcoords.303[,c("LongitudePopLevel3", "LatitudePopLevel3")]




vcf.75pctmsng.artif <- read.vcfR("75pctmsng_nosingletons/vcftools_artifcluster_75pctmsng_oneSNPnosngl_20181210.recode.vcf") # 6056 SNPs; excludes singletons and private doubletons
# 3389 if we only exclude singletons
genind.75pctmsng.artif <- vcfR2genind(vcf.75pctmsng.artif)
indNames(genind.75pctmsng.artif) <- popcoords.artif$IndPop2
pop(genind.75pctmsng.artif) <- popcoords.artif$PopLevel1
genind.75pctmsng.artif@other$xy <- popcoords.artif[,c("Longitude2", "Latitude2")]
genind.75pctmsng.artif.seppops <- seppop(genind.75pctmsng.artif)
genpop.75pctmsng.artif <- genind2genpop(genind.75pctmsng.artif)
gl.75pctmsng.artif <- vcfR2genlight(vcf.75pctmsng.artif)
ploidy(gl.75pctmsng.artif) <- 2
indNames(gl.75pctmsng.artif) <- popcoords.artif$IndPop2
pop(gl.75pctmsng.artif) <- popcoords.artif$PopLevel1
gl.75pctmsng.artif@other$xy <- popcoords.artif[,c("Longitude2", "Latitude2")]

vcf.75pctmsng.coastal <- read.vcfR("75pctmsng_nosingletons/vcftools_coastalcluster_75pctmsng_oneSNPnosngl_20181210.recode.vcf") # 6056 SNPs; excludes singletons and private doubletons
# 3073 if we only exclude singletons
genind.75pctmsng.coastal <- vcfR2genind(vcf.75pctmsng.coastal)
indNames(genind.75pctmsng.coastal) <- popcoords.coastal$IndPop2
pop(genind.75pctmsng.coastal) <- popcoords.coastal$PopLevel1
genind.75pctmsng.coastal@other$xy <- popcoords.coastal[,c("Longitude2", "Latitude2")]
genind.75pctmsng.coastal.seppops <- seppop(genind.75pctmsng.coastal)
genpop.75pctmsng.coastal <- genind2genpop(genind.75pctmsng.coastal)
gl.75pctmsng.coastal <- vcfR2genlight(vcf.75pctmsng.coastal)
ploidy(gl.75pctmsng.coastal) <- 2
indNames(gl.75pctmsng.coastal) <- popcoords.coastal$IndPop2
pop(gl.75pctmsng.coastal) <- popcoords.coastal$PopLevel1
gl.75pctmsng.coastal@other$xy <- popcoords.coastal[,c("Longitude2", "Latitude2")]

vcf.75pctmsng.inland <- read.vcfR("75pctmsng_nosingletons/vcftools_inlandcluster_75pctmsng_oneSNPnosngl_20181210.recode.vcf") # 6056 SNPs; excludes singletons and private doubletons
# 5508 if we only exclude singletons
genind.75pctmsng.inland <- vcfR2genind(vcf.75pctmsng.inland)
indNames(genind.75pctmsng.inland) <- popcoords.inland$IndPop2
pop(genind.75pctmsng.inland) <- popcoords.inland$PopLevel1
genind.75pctmsng.inland@other$xy <- popcoords.inland[,c("Longitude2", "Latitude2")]
genind.75pctmsng.inland.seppops <- seppop(genind.75pctmsng.inland)
genpop.75pctmsng.inland <- genind2genpop(genind.75pctmsng.inland)
gl.75pctmsng.inland <- vcfR2genlight(vcf.75pctmsng.inland)
ploidy(gl.75pctmsng.inland) <- 2
indNames(gl.75pctmsng.inland) <- popcoords.inland$IndPop2
pop(gl.75pctmsng.inland) <- popcoords.inland$PopLevel1
gl.75pctmsng.inland@other$xy <- popcoords.inland[,c("Longitude2", "Latitude2")]
gl.75pctmsng.inland.seppops <- seppop()




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

`%notin%` <- Negate(`%in%`)


