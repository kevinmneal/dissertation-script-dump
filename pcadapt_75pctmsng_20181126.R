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

#####

#vcfpcadapt <- vcf2pcadapt("vcftools_allsocal_303indv_highqual_20pctmsng_SNPsinallbignativepops_oneSNP_nosingletons_20181009.recode.vcf", type="vcf")
#vcfpcadapt <- read.pcadapt("vcftools_allsocal_303indv_highqual_20pctmsng_SNPsinallbignativepops_oneSNP_nosingletons_20181009.recode.vcf", type="vcf")
vcfpcadapt <- read.pcadapt("75pctmsng_nosingletons/vcftools_allsocal303indiv_75pctmsng_oneSNPnosngl_20181127.recode.vcf", type="vcf")

read.pcadapt

# K should be chosen as n (number of pops/clusters) - 1
# use K = 3 (i.e. 4 pops) or 14, based on scree plot
pcadapt.k3 <- pcadapt(input=vcfpcadapt, K=3, min.maf=0.0001) # think about choice of K and min.maf here
x <- pcadapt.k3
plot(pcadapt.k3, option="screeplot")
plot(pcadapt.k3, option="scores", pop=popcoords.303$PopLevel3geogenclusters)
plot(pcadapt.k3, option="manhattan")
plot(pcadapt.k3, option="qqplot")
hist(pcadapt.k3$pvalues, xlab = "p-values", main = NULL, breaks = 50, col="orange")
plot(pcadapt.k3, option = "stat.distribution")

pcadapt.k14 <- pcadapt(input=vcfpcadapt, K=14, min.maf=0.0001) # think about choice of K and min.maf here
x <- pcadapt.k14
plot(pcadapt.k14, option="screeplot")
plot(pcadapt.k14, option="scores", pop=popcoords.303$PopLevel3geogenclusters)
plot(pcadapt.k14, option="manhattan")
plot(pcadapt.k14, option="qqplot")
hist(pcadapt.k14$pvalues, xlab = "p-values", main = NULL, breaks = 50, col="orange")
plot(pcadapt.k14, option = "stat.distribution")


# ### choosing outlier cutoff
# # q-values with no correction
# library(qvalue)
# qval <- qvalue(x$pvalues)$qvalues
# alpha <- 0.05
# outliers.qv <- which(qval < alpha)
# length(outliers.qv)
# 
# # benjamini-hochberg
# padj.bh <- p.adjust(x$pvalues,method="BH")
# alpha <- 0.05
# outliers.bh <- which(padj.bh < alpha)
# length(outliers.bh)
# 
# # hochberg
# padj.hoch <- p.adjust(x$pvalues,method="hochberg")
# alpha <- 0.05
# outliers.hoch <- which(padj.hoch < alpha)
# length(outliers.hoch)

# benjamini-yekutieli - use this one. Has been cited in other papers using pcadapt
padj.by.k3 <- p.adjust(pcadapt.k3$pvalues,method="BY")
alpha <- 0.05
outliers.by.k3 <- which(padj.by.k3 < alpha)
length(outliers.by.k3)

padj.by.k14 <- p.adjust(pcadapt.k14$pvalues,method="BY")
alpha <- 0.05
outliers.by.k14 <- which(padj.by.k14 < alpha)
length(outliers.by.k14)

# 
# # holm - dominates the bonferroni method and is valid under arbitrary assumptions
# padj.holm <- p.adjust(x$pvalues,method="holm")
# alpha <- 0.05
# outliers.holm <- which(padj.holm < alpha)
# length(outliers.holm)
# 
# # hommel - more powerful than hochberg's but difference is usually small and computationally slower
# padj.hommel <- p.adjust(x$pvalues,method="hommel")
# alpha <- 0.05
# outliers.hommel <- which(padj.hommel < alpha)
# length(outliers.hommel)
# 
# outliers.qv == outliers.bh
# 
# # bonferroni correction - likely overly conservative, according to qvalue package info
# padj.bf <- p.adjust(pcadapt.run1$pvalues,method="bonferroni")
# alpha <- 0.05
# outliers.bf <- which(padj.bf < alpha)
# length(outliers.bf)
# 
# outliers.intersect <- intersect(outliers.holm, outliers.by)
# length(outliers.intersect)

gl.75pctmsng.allsocal.k3.neutral <- gl.75pctmsng.allsocal[,(gl.75pctmsng.allsocal@loc.names %notin% gl.75pctmsng.allsocal@loc.names[outliers.by.k3])]
gl.75pctmsng.allsocal.k3.outliers <- gl.75pctmsng.allsocal[,gl.75pctmsng.allsocal@loc.names[outliers.by.k3]]

gl.75pctmsng.allsocal.k14.neutral <- gl.75pctmsng.allsocal[,(gl.75pctmsng.allsocal@loc.names %notin% gl.75pctmsng.allsocal@loc.names[outliers.by.k14])]
gl.75pctmsng.allsocal.k14.outliers <- gl.75pctmsng.allsocal[,gl.75pctmsng.allsocal@loc.names[outliers.by.k14]]

testpca.k3.outliers.centered <- glPca(gl.75pctmsng.allsocal.k3.outliers, nf=7, scale=F, center=T) # 7 axes
col <- viridis(length(unique(pop(gl.75pctmsng.allsocal))))
s.class(testpca.k3.outliers.centered$scores, pop(gl.75pctmsng.allsocal.outliers), xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)
testpca.k14.outliers.centered <- glPca(gl.75pctmsng.allsocal.k14.outliers, nf=7, scale=F, center=T)
#col <- viridis(length(unique(pop(gl.75pctmsng.allsocal))))
s.class(testpca.k14.outliers.centered$scores, pop(gl.75pctmsng.allsocal.outliers), xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)

testpca.k3.neutral.centered <- glPca(gl.75pctmsng.allsocal.k3.neutral, nf=7, scale=F, center=T)
s.class(testpca.k3.neutral.centered$scores, pop(gl.75pctmsng.allsocal.k3.neutral), xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)
testpca.k14.neutral.centered <- glPca(gl.75pctmsng.allsocal.k14.neutral, nf=7, scale=F, center=T)
s.class(testpca.k14.neutral.centered$scores, pop(gl.75pctmsng.allsocal.k14.neutral), xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)
# k=14 looks a little more reasonable than K=3
#scatter(testpca.neutral, ratio=0.2)
#col <- funky(length(unique(pop(gl.75pctmsng.allsocal))))
col <- viridis(length(unique(pop(gl.75pctmsng.allsocal))))

#col <- funky(length(unique(pop(gl.75pctmsng.allsocal))))
#col <- viridis(50)

# switch out poplevel1/site labels for labeling the pca plot
s.class(testpca.neutral$scores, pop(gl.75pctmsng.allsocal.neutral), xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)


##### outlier analysis on native socal pops only
vcfpcadapt.native <- read.pcadapt("75pctmsng_nosingletons/vcftools_nativesocal198indiv_75pctmsng_oneSNPnosngl_20181127.recode.vcf", type="vcf")
gl.75pctmsng.nativesocal <- popsub(gl.75pctmsng.allsocal, 
                                   blacklist=c("IM01", "IM02", "IM03", "IM06", "IM07", "IM09", "IM12", "IM14", "SHOE"),
                                   drop=TRUE)

pcadapt.nativesocal <- pcadapt(input=vcfpcadapt.native, K=20, min.maf=0.0001) # think about choice of K and min.maf here
x <- pcadapt.nativesocal
plot(pcadapt.nativesocal, option="screeplot")
plot(pcadapt.nativesocal, option="scores", pop=pop(gl.75pctmsng.nativesocal))
plot(pcadapt.nativesocal, option="manhattan")
plot(pcadapt.nativesocal, option="qqplot")
hist(pcadapt.nativesocal$pvalues, xlab = "p-values", main = NULL, breaks = 50, col="orange")
plot(pcadapt.nativesocal, option = "stat.distribution")

# benjamini-yekutieli - use this one. Has been cited in other papers using pcadapt
padj.nativesocal.by <- p.adjust(pcadapt.nativesocal$pvalues,method="BY")
alpha <- 0.05
outliers.nativesocal.by <- which(padj.nativesocal.by < alpha)
length(outliers.nativesocal.by)

gl.75pctmsng.nativesocal.neutral <- gl.75pctmsng.nativesocal[,(gl.75pctmsng.nativesocal@loc.names %notin% gl.75pctmsng.nativesocal@loc.names[outliers.nativesocal.by])]
gl.75pctmsng.nativesocal.outliers <- gl.75pctmsng.nativesocal[,gl.75pctmsng.nativesocal@loc.names[outliers.nativesocal.by]]

testpca.nativesocal.outliers <- glPca(gl.75pctmsng.nativesocal.outliers)
col2 <- viridis(length(unique(pop(gl.75pctmsng.nativesocal))))
s.class(testpca.nativesocal.outliers$scores, pop(gl.75pctmsng.nativesocal.neutral), xax=1,yax=2, col=transp(col2,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)

testpca.nativesocal.neutral <- glPca(gl.75pctmsng.nativesocal.neutral)
s.class(testpca.nativesocal.neutral$scores, pop(gl.75pctmsng.nativesocal.neutral), xax=1,yax=2, col=transp(col2,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)



### outlier analysis on inland pops only
vcfpcadapt.inland <- read.pcadapt("G:/My Drive/Illumina Sequencing Data/20181128_allsocal303indiv_finaldatasets/75pctmsng_nosingletons/vcftools_inlandcluster_75pctmsng_oneSNPnosngl_20181210.recode.vcf", type="vcf")

pcadapt.inland <- pcadapt(input=vcfpcadapt.inland, K=6, min.maf=0.0001) # think about choice of K and min.maf here
x <- pcadapt.inland
plot(pcadapt.inland, option="screeplot")
plot(pcadapt.inland, option="scores", pop=pop(gl.75pctmsng.inland))
plot(pcadapt.inland, option="manhattan")
plot(pcadapt.inland, option="qqplot")
hist(pcadapt.inland$pvalues, xlab = "p-values", main = NULL, breaks = 50, col="orange")
plot(pcadapt.inland, option = "stat.distribution")

# benjamini-yekutieli - use this one. Has been cited in other papers using pcadapt
padj.inland.by <- p.adjust(pcadapt.inland$pvalues,method="BY")
alpha <- 0.05
outliers.inland.by <- which(padj.inland.by < alpha)
length(outliers.inland.by)

gl.75pctmsng.inland.neutral <- gl.75pctmsng.inland[,(gl.75pctmsng.inland@loc.names %notin% gl.75pctmsng.inland@loc.names[outliers.inland.by])]
gl.75pctmsng.inland.outliers <- gl.75pctmsng.inland[,gl.75pctmsng.inland@loc.names[outliers.inland.by]]

testpca.inland.outliers <- glPca(gl.75pctmsng.inland.outliers)
colinland <- viridis(length(unique(pop(gl.75pctmsng.inland))))
s.class(testpca.inland.outliers$scores, pop(gl.75pctmsng.inland.neutral), xax=1,yax=2, col=transp(colinland,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)

testpca.inland.neutral <- glPca(gl.75pctmsng.inland.neutral)
s.class(testpca.inland.neutral$scores, pop(gl.75pctmsng.inland.neutral), xax=1,yax=2, col=transp(col2,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)




### try lfmm too, in package LEA, and I think there are tess/snmf functions to do something without environmental variables