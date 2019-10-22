### pca and dapc - adegenet

library(adegenet)
library(poppr)
library(PopGenome)
library(adegenet)
library(hierfstat)
library(StAMPP)
library(vcfR)
library(cluster)
library(Matrix)

popcoords.303 <- read.csv("G:/My Drive/Illumina Sequencing Data/20181004_ALLsocal_303indiv_highqual/allsocal_303indiv_popcoords_multiplelevels_betternames.csv",
                          header=T,
                          stringsAsFactors = F)
popcoords <- popcoords.303

setwd("G:/My Drive/Illumina Sequencing Data/20181004_ALLsocal_303indiv_highqual/popfile_allbignativepops/")

vcf.20pctmsng <- read.vcfR("vcftools_allsocal_303indv_highqual_20pctmsng_SNPsinallbignativepops_oneSNP_nosingletons_20181009.recode.vcf")
genind.20pctmsng <- vcfR2genind(vcf.20pctmsng) # __ SNPs
indNames(genind.20pctmsng) <- popcoords$IndPop2
pop(genind.20pctmsng) <- popcoords$PopLevel2_prox
genind.20pctmsng@other$xy <- popcoords[,c("Longitude2", "Latitude2")]
genind.20pctmsng.seppops <- seppop(genind.20pctmsng)
genpop.20pctmsng <- genind2genpop(genind.20pctmsng)
gl.20pctmsng <- vcfR2genlight(vcf.20pctmsng)
ploidy(gl.20pctmsng) <- 2
pop(gl.20pctmsng) <- popcoords$PopLevel2_prox
gl.20pctmsng@other$xy <- popcoords[,c("Longitude2", "Latitude2")]
indNames(gl.20pctmsng) <- popcoords$IndPop2

spha.genind <- genind.20pctmsng

xvaldapc.20pctmsng <- xvalDapc(genind.75pctmsng.allsocal, grp=popcoords$PopLevel3geogenclusters)

spha.pca <- find.clusters(genind.75pctmsng.allsocal, max.n.clust=20)
pca.dapc <- dapc(genind.75pctmsng.allsocal, spha.pca$grp)
compoplot(pca.dapc, show.lab=TRUE, lab=pop(genind.75pctmsng.allsocal), legend=FALSE, main="DAPC K=3, 303 individuals")
scatter(pca.dapc)

par(mfrow=c(1,2))
spea.dudipca <- dudi.pca(df = tab(genind.75pctmsng.allsocal, NA.method = "mean"), scannf = TRUE, scale=F)

col <- funky(20)
s.class(spea.dudipca$li, pop(genind.75pctmsng.allsocal),xax=1,yax=2, col=transp(col,.6), axesell=TRUE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)

pca.20pctmsng <- find.clusters(gl.178indiv)
dapc.20pctmsng <- dapc(gl.178indiv, pca.20pctmsng$grp)
compoplot(dapc.20pctmsng, show.lab=TRUE, lab=pop(gl.178indiv), legend=FALSE, main="DAPC, 303 individuals")
scatter(dapc.20pctmsng, label=pop(gl.181indiv), scree.da=FALSE)
temp <- optim.a.score(dapc.20pctmsng) # optimal looks to be 13 PCs

xval <- xvalDapc(tab(genind.181indiv, NA.method="mean"), pop(genind.181indiv), n.pca.max = 200, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
scatter(xval$DAPC, label=pop(genind.20pctmsng)) # well this doesn't seem useful...


# snapclust k selection might still be useful...
spha.choosek.aic <- snapclust.choose.k(30, genind.75pctmsng.allOChwe, IC=AIC)
spha.choosek.aicc <- snapclust.choose.k(30, genind.75pctmsng.allOChwe, IC=AICc)
spha.choosek.bic <- snapclust.choose.k(30, genind.75pctmsng.allOChwe, IC=BIC)
spha.choosek.kic <- snapclust.choose.k(30, genind.75pctmsng.allOChwe, IC=KIC)

plot(spha.choosek.aic) # allsocal: best K=3, 10, or 14; allOC: 3, 8, 10, 11, 16
plot(spha.choosek.aicc)
plot(spha.choosek.bic) # allsocal: K=3; also 3 for allOC
plot(spha.choosek.kic) # allsocal: K=3 or 10; 3 or 8 if looking at allOC

snapk2 <- snapclust(genind.75pctmsng.allOChwe, k=2)
snapk3 <- snapclust(genind.75pctmsng.allOChwe, k=3)
snapk4 <- snapclust(genind.75pctmsng.allOChwe, k=4)
snapk8 <- snapclust(genind.75pctmsng.allOChwe, k=8)
snapk10 <- snapclust(genind.75pctmsng.allOChwe, k=10)


par(mfrow=c(3,1))

compoplot(snapk2$proba, 
          #border="black",
          main="Spea hammondii socal, K=2 (snapclust)",
          lab=pop(genind.75pctmsng.allOChwe),
          show.lab = TRUE,
          legend=FALSE)
compoplot(snapk3$proba, 
          #border="black",
          lab=pop(genind.75pctmsng.allOChwe),
          main="Spea hammondii socal, K=3 (snapclust)",
          show.lab = TRUE,
          legend=FALSE)
compoplot(snapk4$proba, 
          #border="black",
          main="Spea hammondii socal, K=10 (snapclust)",
          show.lab=TRUE)


testpca.75pctmsng.allsocal <- glPca(gl.75pctmsng.allsocal)
#scatter(testpca, ratio=0.2)
#col <- funky(length(unique(pop(gl.75pctmsng))))
col <- viridis(length(unique(pop(gl.75pctmsng))))
s.class(testpca.75pctmsng.allsocal$scores, pop(gl.75pctmsng.allsocal),xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)
# try clusplot or just ggplot

# SIAP
testpca.SIAP.allsocal <- glPca(gl.SIAP.allsocal)
#col <- viridis(length(unique(pop(gl.75pctmsng.allsocal))))
s.class(testpca.SIAP.allsocal$scores, pop(gl.75pctmsng.allsocal),xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)

# 20pctmsng - no SIAP
testpca.20pctmsng.allsocal <- glPca(gl.20pctmsng.allsocal)
#col <- viridis(length(unique(pop(gl.75pctmsng))))
s.class(testpca.20pctmsng.allsocal$scores, pop(gl.75pctmsng.allsocal),xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)

# 5pctmsng - no SIAP
testpca.5pctmsng.allsocal <- glPca(gl.5pctmsng.allsocal)
#col <- viridis(length(unique(pop(gl.75pctmsng))))
s.class(testpca.5pctmsng.allsocal$scores, pop(gl.75pctmsng.allsocal),xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)

# nomsng - no SIAP
testpca.nomsng.allsocal <- glPca(gl.nomsng.allsocal)
#col <- viridis(length(unique(pop(gl.75pctmsng))))
s.class(testpca.nomsng.allsocal$scores, pop(gl.75pctmsng.allsocal),xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)


# test glPca params
testpca.75pctmsng.allsocal <- glPca(gl.75pctmsng.allsocal, scale=F, center=F, nf=5)
s.class(testpca.75pctmsng.allsocal$scores, pop(gl.75pctmsng.allsocal),xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)
# interestingly, when scale and center=F, plotting axes 2 and 3 looks the way axes 1 and 2 look when center=T
# also looks like, by not centering or scaling, pc1 actually separates LAVENTURA pops from the rest;
# pc2 separates coastal, inland OC, and artificial pools
# BUT it isn't really a pca if you don't center it...
testpca.75pctmsng.nativesocal <- glPca(gl.75pctmsng.nativesocal, scale=F, center=F, nf=5)
col2 <- viridis(length(unique(pop(gl.75pctmsng.nativesocal))))
s.class(testpca.75pctmsng.nativesocal$scores, pop(gl.75pctmsng.nativesocal),xax=1,yax=2, col=transp(col2,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)


# artificial pools alone
testpca.75pctmsng.artif.centered <- glPca(gl.75pctmsng.artif, scale=F, center=T, nf=5)
colartif <- viridis(length(unique(pop(gl.75pctmsng.artif))))
s.class(testpca.75pctmsng.artif.centered$scores, pop(gl.75pctmsng.artif),xax=1,yax=2, col=transp(colartif,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)

# coastal pops alone
# note: you should center your data... centering is a part of how PCA is done, otherwise it's a different sort of analysis?
testpca.75pctmsng.coastal <- glPca(gl.75pctmsng.coastal, scale=F, center=F, nf=5)
colcoastal <- viridis(length(unique(pop(gl.75pctmsng.coastal))))
s.class(testpca.75pctmsng.coastal$scores, pop(gl.75pctmsng.coastal),xax=1,yax=2, col=transp(colcoastal,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)
testpca.75pctmsng.coastal$eig[[1]]/sum(testpca.75pctmsng.coastal$eig)*100 
testpca.75pctmsng.coastal.centered <- glPca(gl.75pctmsng.coastal, scale=F, center=T, nf=5)
#colcoastal <- viridis(length(unique(pop(gl.75pctmsng.coastal))))
s.class(testpca.75pctmsng.coastal.centered$scores, pop(gl.75pctmsng.coastal),xax=1,yax=2, col=transp(colcoastal,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)

# inland pops alone
testpca.75pctmsng.inland.centered <- glPca(gl.75pctmsng.inland, scale=F, center=T, nf=5)
colinland <- viridis(length(unique(pop(gl.75pctmsng.inland))))
s.class(testpca.75pctmsng.inland.centered$scores, pop(gl.75pctmsng.inland),xax=1,yax=2, col=transp(colinland,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)
title("pca.75pctmsng.inland, PC1 PC2")
s.class(testpca.75pctmsng.inland.centered$scores, pop(gl.75pctmsng.inland),xax=3,yax=4, col=transp(colinland,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)
title("pca.75pctmsng.inland, PC3 PC4")

gl.75pctmsng.inlandoc <- popsub(gl.75pctmsng.inland, blacklist=c("MISSION", "BOUQ", "SOBOB1", "SOBOB2", "HCC1", "HCC2", "SHUL", "LEM", "CAST"))
pca.75pctmsng.inlandoc <- glPca(gl.75pctmsng.inlandoc, scale=F, center=T, nf=8)
colinlandoc <- viridis(length(unique(pop(gl.75pctmsng.inlandoc))))
par(mfrow=c(2,2))
s.class(pca.75pctmsng.inlandoc$scores, pop(gl.75pctmsng.inlandoc),xax=1,yax=2, col=transp(colinlandoc,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=0.7, grid=FALSE, chullSize=1)
title("pca.75pctmsng.inlandOC, PC1 PC2")
s.class(pca.75pctmsng.inlandoc$scores, pop(gl.75pctmsng.inlandoc),xax=3,yax=4, col=transp(colinlandoc,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=0.7, grid=FALSE)
title("pca.75pctmsng.inlandOC, PC3 PC4")
s.class(pca.75pctmsng.inlandoc$scores, pop(gl.75pctmsng.inlandoc),xax=5,yax=6, col=transp(colinlandoc,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=0.7, grid=FALSE)
title("pca.75pctmsng.inlandOC, PC5 PC6")
s.class(pca.75pctmsng.inlandoc$scores, pop(gl.75pctmsng.inlandoc),xax=7,yax=8, col=transp(colinlandoc,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=0.7, grid=FALSE)
title("pca.75pctmsng.inlandOC, PC7 PC8")
#pairs(pca.75pctmsng.inlandoc$scores[, c(1:3)], pch=19)
s.chull(pca.75pctmsng.inlandoc$scores, pop(gl.75pctmsng.inlandoc),xax=1,yax=2, col=transp(colinlandoc,.6),
        cpoint=3, clabel=0.5, grid=FALSE, optchull=1)
title("pca.75pctmsng.inlandOC, PC1 PC2")
s.chull(pca.75pctmsng.inlandoc$scores, pop(gl.75pctmsng.inlandoc),xax=3,yax=4, col=transp(colinlandoc,.6),
        cpoint=3, clabel=0.5, grid=FALSE, optchull=1)
title("pca.75pctmsng.inlandOC, PC3 PC4")
s.chull(pca.75pctmsng.inlandoc$scores, pop(gl.75pctmsng.inlandoc),xax=5,yax=6, col=transp(colinlandoc,.6),
        cpoint=3, clabel=0.5, grid=FALSE, optchull=1)
title("pca.75pctmsng.inlandOC, PC5 PC6")
s.chull(pca.75pctmsng.inlandoc$scores, pop(gl.75pctmsng.inlandoc),xax=7,yax=8, col=transp(colinlandoc,.6),
        cpoint=3, clabel=0.5, grid=FALSE, optchull=1)
title("pca.75pctmsng.inlandOC, PC7 PC8")

# percent variation explained by each PC
plot(pca.75pctmsng.inlandoc$eig/sum(pca.75pctmsng.inlandoc$eig)*100)

gl.75pctmsng.inlandoc2 <- popsub(gl.75pctmsng.inland, blacklist=c("STARR", "MISSION", "BOUQ", "SOBOB1", "SOBOB2", "HCC1", "HCC2", "SHUL", "LEM", "CAST"))
pca.75pctmsng.inlandoc2 <- glPca(gl.75pctmsng.inlandoc2, scale=F, center=T, nf=6)
colinlandoc2 <- viridis(length(unique(pop(gl.75pctmsng.inlandoc2))))
s.class(pca.75pctmsng.inlandoc2$scores, pop(gl.75pctmsng.inlandoc2),xax=1,yax=2, col=transp(colinlandoc2,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)
#pairs(pca.75pctmsng.inlandoc2$scores[, c(1:3)], pch=19)
gl.75pctmsng.allochwe <- gl.75pctmsng.allOChwe

### did some HWE filtering - trying these out now
pca.75pctmsng.allochwe2 <- glPca(gl.75pctmsng.allOChwe, scale=F, center=T)
#pca.75pctmsng.allochwek3 <- pca.75pctmsng.allochwe 
#pca.75pctmsng.allochwe.sitesaspops <- pca.75pctmsng.allochwe
pca.75pctmsng.allochwe.p001 <- pca.75pctmsng.allochwe
colalloc <- viridis(length(unique(pop(gl.75pctmsng.allOChwe))))
s.class(pca.75pctmsng.allochwe$scores, pop(gl.75pctmsng.allOChwe),xax=1,yax=2, col=transp(colalloc,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)
s.chull(pca.75pctmsng.allochwe$scores, pop(gl.75pctmsng.allochwe),xax=1,yax=2, col=transp(colalloc,.6),
        cpoint=3, clabel=1, grid=FALSE, optchull=1)
title("pca.75pctmsng.allochwe, PC1 PC2")
s.chull(pca.75pctmsng.allochwe$scores, pop(gl.75pctmsng.allochwe),xax=3,yax=4, col=transp(colalloc,.6),
        cpoint=3, clabel=1, grid=FALSE, optchull=1)
title("pca.75pctmsng.allochwe, PC3 PC4")
s.chull(pca.75pctmsng.allochwe$scores, pop(gl.75pctmsng.allochwe),xax=5,yax=6, col=transp(colalloc,.6),
        cpoint=3, clabel=1, grid=FALSE, optchull=1)
title("pca.75pctmsng.allochwe, PC5 PC6")
s.chull(pca.75pctmsng.allochwe$scores, pop(gl.75pctmsng.allochwe),xax=7,yax=8, col=transp(colalloc,.6),
        cpoint=3, clabel=1, grid=FALSE, optchull=1)
title("pca.75pctmsng.allochwe, PC7 PC8")

pca.75pctmsng.allsocalhwe <- glPca(gl.75pctmsng.allsocalhwe, scale=F, center=T, nf=8)
#pca.75pctmsng.allsocalhwe.k3 <- pca.75pctmsng.allsocalhwe
#pca.75pctmsng.allsocalhwe.sitesaspops <- pca.75pctmsng.allsocalhwe
pca.75pctmsng.allsocalhwe.p001 <- pca.75pctmsng.allsocalhwe
colallsocal <- viridis(length(unique(pop(gl.75pctmsng.allsocalhwe))))
s.class(pca.75pctmsng.allsocalhwe$scores, pop(gl.75pctmsng.allsocalhwe),xax=1,yax=2, col=transp(colallsocal,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=1, grid=FALSE)

s.chull(pca.75pctmsng.allsocalhwe$scores, pop(gl.75pctmsng.allsocalhwe),xax=1,yax=2, col=transp(colallsocal,.6),
        cpoint=3, clabel=1, grid=FALSE, optchull=1)
title("pca.75pctmsng.allsocalhwe, PC1 PC2")
