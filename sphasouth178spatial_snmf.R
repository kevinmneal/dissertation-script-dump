# sNMF - sphasouth178spatial
# 20190322
library(vcfR)
library(pals)
library(adegenet)
library(LEA)
library(tess3r)
library(hierfstat)
library(poppr)
library(radiator)
library(assigner) # use to run global and pairwise WCFst with confidence intervals

setwd("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/snmf")
geno.sphasouth178spatialpath <- vcf2geno(input.file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial.radiator.50pctmsng.oneSNPmac3.vcf")
#geno.sphasouth178spatial.shortpath <- "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/sphasouth178spatial__oneSNPmac2_20190302.recode.geno"
#sometimes if the filename is too long it won't be able to read the geno, so have to manually shorten it
geno.sphasouth178spatial <- read.geno(input.file=geno.sphasouth178spatialpath)
geno.sphasouth178spatial[geno.sphasouth178spatial==9] <- NA

popcoords.sphasouth178spatial <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/sphasouth178spatial_popfile_coords.csv",
                                          header=T,
                                          stringsAsFactors = F)

colnames(popcoords.sphasouth178spatial)
# snmf also hates long folder paths... 

# try running with alpha=1 and alpha=10 and see which has the lower CE values
# do 3 runs: full SNP set; putatitvely neutral SNP set; and putative outlier SNP set
project.sphasouth178spatial.50pctmsng <- NULL
project.sphasouth178spatial.50pctmsng <- snmf("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial.radiator.50pctmsng.oneSNPmac3.geno",
                                      K = 1:20,
                                      entropy = TRUE,
                                      percentage = 0.05,
                                      tolerance=0.0000001,
                                      repetitions = 10,
                                      iterations = 500,
                                      project = "new")

project.sphasouth178spatial.50pctmsng <- load.snmfProject("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial.radiator.50pctmsng.oneSNPmac3.snmfProject")
plot(project.sphasouth178spatial.50pctmsng)

snmf.sphasouth178spatial.50pctmsng.ce <- NULL
for (i in 1:20){
  snmf.sphasouth178spatial.50pctmsng.ce <- cbind(snmf.sphasouth178spatial.50pctmsng.ce, cross.entropy(project.sphasouth178spatial.50pctmsng, K=i))
}
# looks like K=14 and K=21 are best according to sNMF; maybe K=8 as well

boxplot(snmf.sphasouth178spatial.50pctmsng.ce, xlab="Cross-entropy", main="sphanorth124.spatial.50pctmsng.snmf")

#plot(project.sphasouth178spatial.50pctmsng)

# par(mfrow = c(9,1),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1) # change to (4,0,1,1) if using a unique axis for each chart

cairo_pdf(file="sNMF_K2-20_sortbyQ_sphasouth178spatial_radiator_50pctmsng.pdf",
    width=20,
    height=50)

par(mfrow = c(19,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(8,0,1,1) + 0.1) # change to (4,0,1,1) if using a unique axis for each chart

for (i in 2:20){
  LEA::barchart(project.sphasouth178spatial.50pctmsng,
                K=i,
                run=which.min(cross.entropy(project.sphasouth178spatial.50pctmsng, K=i)),
                border=NA,
                space=0,
                col=cols25(19),
                main=paste0("sNMF sorted by Q, K=", i),
                sort.by.Q=TRUE) -> sortbyq10
  axis(1, at = 1:length(sortbyq10$order),
       labels = popcoords.sphasouth178spatial$poplevel2[sortbyq10$order], las=2,
       cex.axis = 1)
}
dev.off() # only use if using pdf() before generating the plots

cairo_pdf(file="sNMF_K2-20_sphasouth178spatial_radiator_50pctmsng.pdf",
    width=20,
    height=50)

par(mfrow = c(19,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(8,0,1,1) + 0.1) # change to (4,0,1,1) if using a unique axis for each chart

for (i in 2:20){
  LEA::barchart(project.sphasouth178spatial.50pctmsng,
                K=i,
                run=which.min(cross.entropy(project.sphasouth178spatial.50pctmsng, K=i)),
                border=NA,
                space=0,
                col=cols25(19),
                main=paste0("sNMF, K=", i),
                sort.by.Q=FALSE) -> sortbyq10
  axis(1, at = 1:length(sortbyq10$order),
       labels = popcoords.sphasouth178spatial$poplevel2, las=2,
       cex.axis = 1)
}
dev.off() # only use if using pdf() before generating the plots





cairo_pdf(file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sNMF_bestKs_sortbyQ_sphasouth178spatial_radiator_50pctmsng.pdf",
          width=20,
          height=20)

par(mfrow = c(6,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(8,0,1,1) + 0.1) # change to (4,0,1,1) if using a unique axis for each chart

for (i in c(2,3,4,5,6,10)){
  LEA::barchart(project.sphasouth178spatial.50pctmsng,
                K=i,
                run=which.min(cross.entropy(project.sphasouth178spatial.50pctmsng, K=i)),
                border=NA,
                space=0,
                col=cols25(19),
                main=paste0("sNMF sorted by Q, K=", i),
                sort.by.Q=TRUE) -> sortbyq10
  axis(1, at = 1:length(sortbyq10$order),
       labels = popcoords.sphasouth178spatial$poplevel2[sortbyq10$order], las=2,
       cex.axis = 1)
}
dev.off() # only use if using pdf() before generating the plots

cairo_pdf(file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sNMF_bestKs_2-10_sphasouth178spatial_radiator_50pctmsng.pdf",
          width=18,
          height=15)
#CairoSVG(file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sNMF_bestKs_2-9_sphasouth178spatial_radiator_50pctmsng.svg",
#          width=14,
#          height=21)

par(mfrow = c(6,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(6,0,1,1) + 0.1) # change to (4,0,1,1) if using a unique axis for each chart

for (i in c(2,3,4,5,6,10)){
  LEA::barchart(project.sphasouth178spatial.50pctmsng,
                K=i,
                run=which.min(cross.entropy(project.sphasouth178spatial.50pctmsng, K=i)),
                border=NA,
                space=0,
                col=cols25(19),
                main=paste0("sNMF, K=", i),
                sort.by.Q=FALSE) -> sortbyq10
  axis(1, at = 1:length(sortbyq10$order),
       labels = popcoords.sphasouth178spatial$poplevel2, las=2,
       cex.axis = 1)
}
dev.off() # only use if using pdf() before generating the plots
















#### try with cv-10


project.sphasouth178spatial.50pctmsng.pt1 <- NULL
project.sphasouth178spatial.50pctmsng.pt1 <- snmf("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial.radiator.50pctmsng.oneSNPmac3.geno",
                                              K = 1:20,
                                              entropy = TRUE,
                                              percentage = 0.1,
                                              tolerance=0.000001,
                                              repetitions = 10,
                                              iterations = 200,
                                              project = "new")

#project.sphasouth178spatial.50pctmsng.pt1 <- load.snmfProject("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial.radiator.50pctmsng.oneSNPmac3.snmfProject")
plot(project.sphasouth178spatial.50pctmsng.pt1)

snmf.sphasouth178spatial.50pctmsng.ce <- NULL
for (i in 1:20){
  snmf.sphasouth178spatial.50pctmsng.ce <- cbind(snmf.sphasouth178spatial.50pctmsng.ce, which.min(cross.entropy(project.sphasouth178spatial.50pctmsng.pt1, K=i)))
}
# looks like K=14 and K=21 are best according to sNMF; maybe K=8 as well

boxplot(snmf.sphasouth178spatial.50pctmsng.ce, xlab="Cross-entropy", main="sphanorth124.spatial.50pctmsng.snmf")

#plot(project.sphasouth178spatial.50pctmsng.pt1)

# par(mfrow = c(9,1),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1) # change to (4,0,1,1) if using a unique axis for each chart

Cairo(file="sNMF_K2-3-7_sphasouth178spatial_radiator_50pctmsng_20190607.png",
          width=20,
          height=10,
         units="in",
         res=320,
      type="png",
      bg="white")

CairoPDF(file="sNMF_K2-3-7_sphasouth178spatial_radiator_50pctmsng_20190607.png",
      width=20,
      height=10)
par(mfrow = c(3,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(8,0,1,1) + 0.1) # change to (4,0,1,1) if using a unique axis for each chart

for (i in c(2,3,7)){
  LEA::barchart(project.sphasouth178spatial.50pctmsng.pt1,
                K=i,
                run=which.min(cross.entropy(project.sphasouth178spatial.50pctmsng.pt1, K=i)),
                border=NA,
                space=0,
                col=cols25(19),
                main=paste0("SPHA-SOUTH sNMF clusters, K=", i),
                sort.by.Q=FALSE) -> sortbyq10
  axis(1, at = 1:length(sortbyq10$order),
       labels = popcoords.sphasouth178spatial$poplevel2[sortbyq10$order], las=2,
       cex.axis = 1)
}
dev.off() # only use if using pdf() before generating the plots

LEA::barchart(project.sphasouth178spatial.50pctmsng.pt1,
              K=2,
              run=which.min(cross.entropy(project.sphasouth178spatial.50pctmsng.pt1, K=2)),
              border=NA,
              space=0,
              col=cols25(19),
              main=paste0("sNMF sorted by Q, K=", i),
              sort.by.Q=TRUE) -> sortbyq10
axis(1, at = 1:length(sortbyq10$order),
     labels = popcoords.sphasouth178spatial$poplevel2[sortbyq10$order], las=2,
     cex.axis = 1)





### snmf has fst outlier tests...

p.snmf = snmf.pvalues(project.sphasouth178spatial.50pctmsng.pt1, entropy = TRUE, ploidy = 2, K = 13)
p.snmf$GIF

par(mfrow = c(2,1))
hist(p.snmf$pvalues, col = "orange")

plot(-log10(p.snmf$pvalues), pch = 19, col = "blue", cex = .7)




snmf.sphasouth178spatial.50pctmsng.pt1.ce <- NULL
for (i in 1:20){
  snmf.sphasouth178spatial.50pctmsng.pt1.ce <- rbind(snmf.sphasouth178spatial.50pctmsng.pt1.ce, cross.entropy(project.sphasouth178spatial.50pctmsng.pt1, K=i, run=which.min(cross.entropy(project.sphasouth178spatial.50pctmsng.pt1, K=i))))
}

plot(snmf.sphasouth178spatial.50pctmsng.pt1.ce)

Q(project.sphasouth178spatial.50pctmsng.pt1, K=2, run=which.min(cross.entropy(project.sphasouth178spatial.50pctmsng.pt1, K=2)))

### can interpolate any q values using tess3r. Try it with snmf results:

# tess3
ggtess3Q(qmatrix(tess3.obj, K=5, rep="best"),
         coord=coordinates.sphasouth178spatial,
         raster.filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/cmiRGA_SPHAsouth_noimperviousmask_maxentmask10.asc",
         col.palette=CreatePalette(cols25(10), 20))

# snmf
ggtess3Q(Q(project.sphasouth178spatial.50pctmsng.pt1, K=3, run=which.min(cross.entropy(project.sphasouth178spatial.50pctmsng.pt1, K=3))),
         coord=coordinates.sphasouth178spatial,
         raster.filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/cmiRGA_SPHAsouth_noimperviousmask_maxentmask10.asc",
         col.palette=CreatePalette(cols25(3)[c(1,2,3)], 10))
plot(Q(project.sphasouth178spatial.50pctmsng.pt1, K=3, run=which.min(cross.entropy(project.sphasouth178spatial.50pctmsng.pt1, K=3))),
         coord=coordinates.sphasouth178spatial,
         raster.filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/cmiRGA_SPHAsouth_noimperviousmask_maxentmask10.asc",
         col.palette=CreatePalette(cols25(3)[c(1,2,3)], 10))




### plot cross entropy for snmf, tess3, and BIC from snapclust.choosek.bic

par(mfrow=c(2,2))
boxplot(snmf.sphasouth178spatial.50pctmsng.ce, ylab="Cross-entropy", main="A) sphasouth178spatial.50pctmsng.snmf")
plot(choosek.bic, main="B) sphasouth178spatial.choosek.bic") # sphasouth178spatial: k=3 (lowest)
abline(v=2, col="darkred")
plot(choosek.kic, main="C) sphasouth178spatial.choosek.kic") # sphasouth178spatial: k=3 (elbow) or k=5 (lowest)
abline(v=5, col="darkred")
plot(choosek.aic, main="D) sphasouth178spatial.choosek.aic") # sphasouth178spatial: k=5 (elbow) or 8 (lowest)
abline(v=10, col="darkred")
dev.off()


#### plot north and south choice of K figure
tess3.obj.north <- readRDS("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphanorth124spatial/sphanorth124spatial.choosek.rds")
choosek.bic.north <- readRDS("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphanorth124spatial/sphanorth124spatial.choosekbic.rds")
snmfcrossval.north <- readRDS("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphanorth124spatial/sphanorth124spatial.snmfcrossval.rds")

CairoPNG(file = "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphanorth124spatial/sphanorthsouth.onefig.choiceofk.png",
         units="in", height=7, width=10, dpi=320)
par(mfrow=c(2,3))
plot(snmfcrossval.north, ylab="Cross-entropy", xlab="K", main="A) SPHA-NORTH sNMF cross-entropy")
plot(tess3.obj.north, ylab="Cross-entropy", main="B) SPHA-NORTH TESS3 cross-entropy")
plot(choosek.bic.north, main="C) SPHA-NORTH BIC", xlab="K", ylab="BIC")
abline(v=3, col="red")

plot(snmf.sphasouth178spatial.50pctmsng.pt1.ce, ylab="Cross-entropy", xlab="K", main="D) SPHA-SOUTH sNMF cross-entropy")
plot(tess3.obj, ylab="Cross-entropy", main="E) SPHA-SOUTH TESS3 cross-entropy")
plot(choosek.bic, main="F) SPHA-SOUTH BIC", xlab="K", ylab="BIC")
abline(v=2, col="red")

dev.off()

