### Spea genetic in R - outlier loci with pcadapt
### .R file creation date: 2018-08-28
### Last updated: 2018-11-30

### vignettes at https://bcm-uga.github.io/pcadapt/articles/index.html
### paper at https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.12592

# since the artificial pops may have undergone a bottleneck, maybe identifying outlier loci doesn't make sense,
# so maybe only run this on native populations?
# and then run resistanceGA using both neutral and outlier loci, and also faststructure


##### required packages
library(pcadapt)
library(qvalue)
library(adegenet)
library(SNPRelate)
library(vcfR)
library(poppr)
library(PopGenome)
library(dplyr)
library(ape)
library(RColorBrewer)
library(igraph)
library(gplots)
library(viridis)
library(tess3r)
library(rworldmap)
library(maps)
library(gstudio)
library(hierfstat)
library(pegas)
library(radiator)
library(pals)

#####

vcf2pcadapt("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial.radiator.50pctmsng.oneSNPmac3.vcf",
              out="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial.radiator.50pctmsng.oneSNPmac3.pcadapt")

#vcfpcadapt <- read.pcadapt("vcftools_allsocal_303indv_highqual_20pctmsng_SNPsinallbignativepops_oneSNP_nosingletons_20181009.recode.vcf", type="vcf")
vcfpcadapt.sphasouth178spatial.50pctmsng <- read.pcadapt("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial.radiator.50pctmsng.oneSNPmac3.pcadapt", type="pcadapt")

popcoords.sphasouth178spatial <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/sphasouth178spatial_popfile_coords.csv",
                                          header=T,
                                          stringsAsFactors = F)

# K should be chosen as n (number of pops/clusters) - 1
# so can use your snmf/tess3/whatever optimal K and subtract 1 (because need n-1 PCs to discriminate among n clusters)
# use K = 3 (i.e. 4 pops) or 14, based on scree plot
plot(pcadapt(input=vcfpcadapt.sphasouth178spatial.50pctmsng, K=20, min.maf=0.0001), option="screeplot") # think about choice of K and min.maf here
pcadapt.sphasouth178spatial.50pctmsng <- pcadapt(input=vcfpcadapt.sphasouth178spatial.50pctmsng, K=3, min.maf=0.0001) # think about choice of K and min.maf here
x.50pctmsng <- pcadapt.sphasouth178spatial.50pctmsng
plot(x.50pctmsng, option="screeplot") # 2 or 5 at min.maf=0.02 or 0.05; 2 or 4 at 0.0001 or 0.01; choice should be optimal minus 1 using Cattell's Rule
plot(x.50pctmsng, option="scores", pop=popcoords.sphasouth178spatial$poplevel2, col=glasbey(n=length(unique(popcoords.sphasouth178spatial$poplevel2))))
plot(x.50pctmsng, option="scores", pop=popcoords.sphasouth178spatial$poplevel2, i=1,j=3, col=glasbey(n=length(unique(popcoords.sphasouth178spatial$poplevel2))))
plot(x.50pctmsng, option="scores", pop=popcoords.sphasouth178spatial$poplevel2, i=5,j=6, col=glasbey(n=length(unique(popcoords.sphasouth178spatial$poplevel2))))
plot(x.50pctmsng, option="scores", pop=popcoords.sphasouth178spatial$poplevel2, i=7,j=8, col=glasbey(n=length(unique(popcoords.sphasouth178spatial$poplevel2))))
plot(x.50pctmsng, option="scores", pop=popcoords.sphasouth178spatial$poplevel2, i=9,j=10, col=glasbey(n=length(unique(popcoords.sphasouth178spatial$poplevel2))))
plot(x.50pctmsng, option="scores", pop=popcoords.sphasouth178spatial$poplevel2)
plot(x.50pctmsng, option="manhattan")
plot(x.50pctmsng, option="qqplot")
hist(x.50pctmsng$pvalues, xlab = "p-values", main = NULL, breaks = 50, col="orange")
plot(x.50pctmsng, option = "stat.distribution")
plot(x.50pctmsng$loadings[,1])

# benjamini-yekutieli - use this one. Has been cited in other papers using pcadapt
padj.by.50pctmsng <- p.adjust(x.50pctmsng$pvalues,method="BY")
alpha <- 0.05 # use 0.01 if being conservative?
outliers.by.pt01 <- which(padj.by.50pctmsng < 0.01)
outliers.by.pt05 <- which(padj.by.50pctmsng < 0.05)
outliers.by.pt10 <- which(padj.by.50pctmsng < 0.10)
neutral.by.pt10 <- which(padj.by.50pctmsng >= 0.10) # bear in mind these are neutral SNPs above maf=0.05

length(outliers.by.pt10)
length(neutral.by.pt10)
snp_pc.50pctmsng <- get.pc(x.50pctmsng, outliers.by.pt01)

# to be conservative, treat < 0.01 as outliers and > 0.1 as neutral

vcfpcadapt.sphasouth178spatial.50pctmsng


vcfR.sphasouth178spatial.50pctmsng <- read.vcfR("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial.radiator.50pctmsng.oneSNPmac3.vcf")
vcfR.sphasouth178spatial.50pctmsng.chroms <- getCHROM(vcfR.sphasouth178spatial.50pctmsng)
length(vcfR.sphasouth178spatial.50pctmsng.chroms)
pcadapt.sphasouth178spatial.50pctmsng.outlierspt01 <- vcfR.sphasouth178spatial.50pctmsng.chroms[outliers.by.pt01]
pcadapt.sphasouth178spatial.50pctmsng.neutralpt01 <- vcfR.sphasouth178spatial.50pctmsng.chroms[!(vcfR.sphasouth178spatial.50pctmsng.chroms %in% pcadapt.sphasouth178spatial.50pctmsng.outlierspt01)]
pcadapt.sphasouth178spatial.50pctmsng.outlierspt05 <- vcfR.sphasouth178spatial.50pctmsng.chroms[outliers.by.pt05]
pcadapt.sphasouth178spatial.50pctmsng.neutralpt05 <- vcfR.sphasouth178spatial.50pctmsng.chroms[!(vcfR.sphasouth178spatial.50pctmsng.chroms %in% pcadapt.sphasouth178spatial.50pctmsng.outlierspt05)]
pcadapt.sphasouth178spatial.50pctmsng.outlierspt10 <- vcfR.sphasouth178spatial.50pctmsng.chroms[outliers.by.pt10]
pcadapt.sphasouth178spatial.50pctmsng.neutralpt10 <- vcfR.sphasouth178spatial.50pctmsng.chroms[!(vcfR.sphasouth178spatial.50pctmsng.chroms %in% pcadapt.sphasouth178spatial.50pctmsng.outlierspt10)]
pcadapt.sphasouth178spatial.50pctmsng.neutralpt10 <- vcfR.sphasouth178spatial.50pctmsng.chroms[neutral.by.pt10]


pcadapt(read.pcadapt(pcadapt.sphasouth178spatial.50pctmsng.neutralpt10))


length(pcadapt.sphasouth178spatial.50pctmsng.neutralpt10)
# use these to whitelist an outlier and neutral set with radiator, and run sNMF on both sets

# tess3 outliers

tess3outliers <- pvalue(tess3.obj, K=3)
snmfoutliers <- snmf.pvalues(project.sphasouth178spatial.50pctmsng.pt1, 
                             genomic.control=TRUE,
                             entropy=TRUE,
                             K=3,
                             ploidy=2)
hist(tess3outliers)

# benjamini-yekutieli - use this one. Has been cited in other papers using pcadapt
padj.tess3.by <- p.adjust(tess3outliers,method="BY")
alpha <- 0.05 # use 0.01 if being conservative?
outliers.tess3.by.pt01 <- which(padj.tess3.by < 0.01)
outliers.tess3.by.pt05 <- which(padj.tess3.by < 0.05)
outliers.tess3.by.pt10 <- which(padj.tess3.by < 0.10)
length(outliers.tess3.by.pt10)
plot(padj.tess3.by)

intersect(outliers.tess3.by.pt10, outliers.by.pt01)

padj.snmf.by <- p.adjust(snmfoutliers$pvalues,method="BY")
alpha <- 0.05 # use 0.01 if being conservative?
outliers.snmf.by.pt01 <- which(padj.snmf.by < 0.01)
outliers.snmf.by.pt05 <- which(padj.snmf.by < 0.05)
outliers.snmf.by.pt10 <- which(padj.snmf.by < 0.10)
neutral.snmf.by.pt10 <- which(padj.snmf.by >= 0.10)
length(outliers.snmf.by.pt10)
plot(padj.snmf.by)

intersect(outliers.snmf.by.pt10, outliers.by.pt10)


plot(tess3.obj)
plot(project.sphasouth178spatial.50pctmsng.pt1)




#gl.sphasouth178spatial.radiator.50pctmsng <- readRDS("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85_radiator_final/sphasouth178spatial_radiator_mac3/19_radiator_genomic_converter_20190322@2235/radiator_genlight_20190322@2236.RData")
`%notin%` <- Negate(`%in%`)

gl.sphasouth178spatial.radiator.50pctmsng <- vcfR2genlight(vcfR.sphasouth178spatial.50pctmsng)
indNames(gl.sphasouth178spatial.radiator.50pctmsng)
colnames(sphasouth178spatial.vcf@gt)[2:179] == indNames(gl.sphasouth178spatial.radiator.50pctmsng) # ALWAYS run this check
pop(gl.sphasouth178spatial.radiator.50pctmsng) <- popcoords.sphasouth178spatial$poplevel2
ploidy(gl.sphasouth178spatial.radiator.50pctmsng) <- 2

#gl.sphasouth178spatial.radiator.50pctmsng.neutral <- gl.sphasouth178spatial.radiator.50pctmsng[,(gl.sphasouth178spatial.radiator.50pctmsng@loc.names %notin% gl.sphasouth178spatial.radiator.50pctmsng@loc.names[outliers.by.pt05])]
#gl.sphasouth178spatial.radiator.50pctmsng.outliers <- gl.sphasouth178spatial.radiator.50pctmsng[,gl.sphasouth178spatial.radiator.50pctmsng@loc.names[outliers.by.pt01]]
gl.sphasouth178spatial.radiator.50pctmsng.neutralpt10 <- gl.sphasouth178spatial.radiator.50pctmsng[,gl.sphasouth178spatial.radiator.50pctmsng@loc.names[neutral.by.pt10]]
gl.sphasouth178spatial.radiator.50pctmsng.outlierspt01 <- gl.sphasouth178spatial.radiator.50pctmsng[,gl.sphasouth178spatial.radiator.50pctmsng@loc.names[outliers.by.pt01]]
gl.sphasouth178spatial.radiator.50pctmsng.outlierspt10 <- gl.sphasouth178spatial.radiator.50pctmsng[,gl.sphasouth178spatial.radiator.50pctmsng@loc.names[outliers.by.pt10]]

gl.sphasouth178spatial.radiator.50pctmsng.snmfoutlierspt10 <- gl.sphasouth178spatial.radiator.50pctmsng[,gl.sphasouth178spatial.radiator.50pctmsng@loc.names[outliers.snmf.by.pt10]]
gl.sphasouth178spatial.radiator.50pctmsng.snmfneutralpt10 <- gl.sphasouth178spatial.radiator.50pctmsng[,gl.sphasouth178spatial.radiator.50pctmsng@loc.names[neutral.snmf.by.pt10]]



pca.sphasouth178spatial.pcadapt.k3.outlierspt10 <- glPca(gl.sphasouth178spatial.radiator.50pctmsng.outlierspt10, nf=8)
col <- glasbey(length(unique(pop(gl.sphasouth178spatial.radiator.50pctmsng.outliers))))
s.class(pca.sphasouth178spatial.pcadapt.k3.outlierspt10$scores, pop(gl.sphasouth178spatial.radiator.50pctmsng.outliers), xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)

pca.sphasouth178spatial.pcadapt.k3.neutralpt10 <- glPca(gl.sphasouth178spatial.radiator.50pctmsng.neutralpt10, nf=8)
col <- glasbey(length(unique(pop(gl.sphasouth178spatial.radiator.50pctmsng.neutral))))
s.class(pca.sphasouth178spatial.pcadapt.k3.neutralpt10$scores, pop(gl.sphasouth178spatial.radiator.50pctmsng.neutralpt10), xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=0.5, grid=FALSE)

pca.sphasouth178spatial.pcadapt.k3.snmfoutlierspt10 <- glPca(gl.sphasouth178spatial.radiator.50pctmsng.snmfoutlierspt10, nf=8)
col <- glasbey(length(unique(pop(gl.sphasouth178spatial.radiator.50pctmsng.outliers))))
s.class(pca.sphasouth178spatial.pcadapt.k3.snmfoutlierspt10$scores, pop(gl.sphasouth178spatial.radiator.50pctmsng.outliers), xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)

pca.sphasouth178spatial.pcadapt.k3.snmfneutralpt10 <- glPca(gl.sphasouth178spatial.radiator.50pctmsng.snmfneutralpt10, nf=8)
col <- glasbey(length(unique(pop(gl.sphasouth178spatial.radiator.50pctmsng.neutral))))
s.class(pca.sphasouth178spatial.pcadapt.k3.snmfneutralpt10$scores, pop(gl.sphasouth178spatial.radiator.50pctmsng.neutralpt10), xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=0.5, grid=FALSE)

