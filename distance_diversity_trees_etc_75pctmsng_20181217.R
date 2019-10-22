### genetic distance, distance tree, and diversity

library(adegenet)
library(poppr)
#library(PopGenome)
library(adegenet)
library(hierfstat)
library(StAMPP)
library(vcfR)
library(cluster)
library(Matrix)

popcoords.303 <- read.csv("G:/My Drive/Illumina Sequencing Data/20181004_ALLsocal_303indiv_highqual/allsocal_303indiv_popcoords_multilevelpops_betternames_20181017.csv",
                          header=T,
                          stringsAsFactors = F)
popcoords <- popcoords.303

setwd("G:/My Drive/Illumina Sequencing Data/20181128_allsocal303indiv_finaldatasets/75pctmsng_nosingletons")

vcf.75pctmsng <- read.vcfR("G:/My Drive/Illumina Sequencing Data/20181128_allsocal303indiv_finaldatasets/75pctmsng_nosingletons/vcftools_allsocal303indiv_75pctmsng_oneSNPnosngl_20181127.recode.vcf")
genind.75pctmsng <- vcfR2genind(vcf.75pctmsng) # 1672 SNPs
indNames(genind.75pctmsng) <- popcoords$IndPop2
pop(genind.75pctmsng) <- popcoords$PopLevel3geogenclusters
# genind.75pctmsng@other$xy <- popcoords[,c("LongitudePopLevel3", "LatitudePopLevel3")]
# genind.75pctmsng.seppops <- seppop(genind.75pctmsng)
# genpop.75pctmsng <- genind2genpop(genind.75pctmsng)
# gl.75pctmsng <- vcfR2genlight(vcf.75pctmsng)
# ploidy(gl.75pctmsng) <- 2
# indNames(gl.75pctmsng) <- popcoords$IndPop2
# pop(gl.75pctmsng) <- popcoords
# gl.75pctmsng@other$xy <- popcoords[,c("LongitudePopLevel3", "LatitudePopLevel3")]
# 


privatealleles.75pctmsng.allOChwe.k3 <- poppr::private_alleles(genind.75pctmsng.allOChwe.k3,
                                                   alleles ~ .,
                                                   report="table",
                                                   level="population",
                                                   count.alleles=FALSE,
                                                   drop=TRUE)
privatealleles.counts.75pctmsng.allOChwe.k3 <- rowSums(privatealleles.75pctmsng.allOChwe.k3)
# out of 11020 alleles (5510 loci):
# INLAND 1117 private alleles; ARTIF, 319; COAST, 158

propshared.75pctmsng.allOChwe.k3 <- pairwise.propShared(genind.75pctmsng.allOChwe.k3)
propshared.75pctmsng.allsocalhwe.k4 <- pairwise.propShared(genind.75pctmsng.allsocalhwe.k4)
propshared.75pctmsng <- pairwise.propShared(genind.75pctmsng.bigpops)

#### Genetic diversity
# Might just use DnaSP 6 - nevermind, I don't trust it... numbers don't make sense

hierfstat.75pctmsng <- genind2hierfstat(genind.75pctmsng)
basicstats.75pctmsng <- basic.stats(hierfstat.75pctmsng)
ho.pops.75pctmsng <- colMeans(basicstats.75pctmsng$Ho, na.rm=T)
#he_adegenet.pops.75pctmsng <- Hs(genind.75pctmsng) 
hs.pops.75pctmsng <- colMeans(basicstats.75pctmsng$Hs, na.rm=T)
fis.pops.75pctmsng <- colMeans(basicstats.75pctmsng$Fis, na.rm=T) # they're all negative... does this make sense?
#inbreeding.pops.75pctmsng <- inbreeding(genind.75pctmsng, pop=pop(genind.75pctmsng), res.type="estimate") # adegenet's way of calculating it; doesn't give pop estimates?
ar.75pctmsng <- allelic.richness(hierfstat.75pctmsng)
ar.pops.75pctmsng <- colMeans(ar.75pctmsng$Ar, na.rm=T)

hohe.pops.75pctmsng <- cbind(ho.pops.75pctmsng, he_adegenet.pops.75pctmsng)
arhohsfis.pops.75pctmsng <- cbind(ar.pops.75pctmsng, ho.pops.75pctmsng, hs.pops.75pctmsng, fis.pops.75pctmsng) # this is "Nei's gene diversity", which is arguably superior to standard Hexp - accounts for sampling error?
# in any case, both He and Hs give similar results. Ho in almost all cases is > He or Hs...

colnames(arhohsfis.pops.75pctmsng) <- c("Ar", "Ho", "He", "Fis")
#summarise(popcoords$PopLevel2_prox)
popcounts <- count(popcoords, PopLevel3geogenclusters)

arhohsfissamplesize.pops.75pctmsng <- cbind(popcounts, arhohsfis.pops.75pctmsng[order(row.names(arhohsfis.pops.75pctmsng)),])
View(arhohsfissamplesize.pops.75pctmsng)
write.csv(arhohsfissamplesize.pops.75pctmsng, file="arhohsfis_75pctmsng_allsocal303.csv")

#ppfisboot.75pctmsng <- boot.ppfis(hierfstat.75pctmsng)
#ppfis.75pctmsng <- boot.ppfis(hierfstat.75pctmsng, nboot=1) # gives same general numbers as basic.stats, so don't bother

popprdiv <- poppr(genind.75pctmsng) # turns out, these numbers basically just reflect sample size. Useless!
# don't bother going further with this.

# global fst
globalfst.75pctmsng <- fstat(genind.75pctmsng.allsocal) # global Fst for all 303 samples is 0.2535; Fit is 0.1614

#### Pairwise Fst (Weir-Cockerman 1984) and Nei's Da (1983), and maybe Jost's D (2008)
# remove pops with 1 or otherwise very few individuals - gets weird results. use poppr popsub

## make subset of populations with >3 individuals

fst.75pctmsng.k3 <- fstat(genind.75pctmsng.allOChwe.k3)
pwfst.75pctmsng.k3 <- pairwise.WCfst(hierfstat.75pctmsng.allOC.k3)

### intrapop genetic diversity (by taking mean or median genetic distances among individuals within a pop)
# genet.dist reads the pop() info; would have to make each individual its own population
# make a function... note, WCfst doesn't work
# easiest way is probably to run on the full 303 individuals and use dplyr or something to find means by population


genind.75pctmsng.allOChwe.k3 <- genind.75pctmsng.allOChwe
pop(genind.75pctmsng.allOChwe.k3) <- popcoords.alloc$PopLevel7K3
hierfstat.75pctmsng.allOC.k3 <- genind2hierfstat(genind.75pctmsng.allOChwe.k3)
genind.75pctmsng.allOChwe.seppops <- seppop(genind.75pctmsng.allOChwe.k3)
genind.INLAND <- genind.75pctmsng.allOChwe.seppops$INLAND
genind.COAST <- genind.75pctmsng.allOChwe.seppops$COAST
genind.ARTIF <- genind.75pctmsng.allOChwe.seppops$ARTIF

pop(genind.INLAND) <- indNames(genind.INLAND)
hierfstat.INLAND <- genind2hierfstat(genind.INLAND)
par(mfrow=c(2,2))
distance.reynold.INLAND <- genet.dist(hierfstat.INLAND, method="Fst") # equivalent (mostly) to Reynold's distance
mean(distance.reynold.INLAND) 
hist(distance.reynold.INLAND, breaks=100)
distancemat.reynold.INLAND <- as.matrix(distance.reynold.INLAND)
rownames(distancemat.reynold.INLAND) <- as.character(unique(pop(genind.INLAND))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.reynold.INLAND) <- as.character(unique(pop(genind.INLAND)))
#diag(distancemat.reynold.INLAND) <- NA
heatmap.2(distancemat.reynold.INLAND, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Reynold's distance (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")
distance.Da.INLAND <- genet.dist(hierfstat.INLAND, method="Da") # equivalent (mostly) to Reynold's distance
mean(distance.Da.INLAND) 
hist(distance.Da.INLAND, breaks=100)
distancemat.reynold.COAST <- as.matrix(distance.reynold.COAST)
rownames(distancemat.reynold.COAST) <- as.character(unique(pop(genind.COAST))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.reynold.COAST) <- as.character(unique(pop(genind.COAST)))
#diag(distancemat.reynold.COAST) <- NA
heatmap.2(distancemat.Da.INLAND, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")

pop(genind.COAST) <- indNames(genind.COAST)
hierfstat.COAST <- genind2hierfstat(genind.COAST)
par(mfrow=c(2,2))
distance.reynold.COAST <- genet.dist(hierfstat.COAST, method="Fst") # equivalent (mostly) to Reynold's distance
mean(distance.reynold.COAST) 
hist(distance.reynold.COAST, breaks=100)
distancemat.reynold.COAST <- as.matrix(distance.reynold.COAST)
rownames(distancemat.reynold.COAST) <- as.character(unique(pop(genind.COAST))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.reynold.COAST) <- as.character(unique(pop(genind.COAST)))
#diag(distancemat.reynold.COAST) <- NA
heatmap.2(distancemat.reynold.COAST, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")
distance.Da.COAST <- genet.dist(hierfstat.COAST, method="Da") # equivalent (mostly) to Reynold's distance
mean(distance.Da.COAST) 
hist(distance.Da.COAST, breaks=100)
distancemat.Da.COAST <- as.matrix(distance.Da.COAST)
rownames(distancemat.Da.COAST) <- as.character(unique(pop(genind.COAST))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.COAST) <- as.character(unique(pop(genind.COAST)))
#diag(distancemat.Da.COAST) <- NA
heatmap.2(distancemat.Da.COAST, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")

pop(genind.ARTIF) <- indNames(genind.ARTIF)
hierfstat.ARTIF <- genind2hierfstat(genind.ARTIF)
par(mfrow=c(2,2))
distance.reynold.ARTIF <- genet.dist(hierfstat.ARTIF, method="Fst") # equivalent (mostly) to Reynold's distance
mean(distance.reynold.ARTIF) 
hist(distance.reynold.ARTIF, breaks=100)
distance.Da.ARTIF <- genet.dist(hierfstat.ARTIF, method="Da") # equivalent (mostly) to Reynold's distance
mean(lower(distance.Da.ARTIF))
hist(distance.Da.ARTIF, breaks=100)
distancemat.reynold.ARTIF <- as.matrix(distance.reynold.ARTIF)
rownames(distancemat.reynold.ARTIF) <- as.character(unique(pop(genind.ARTIF))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.reynold.ARTIF) <- as.character(unique(pop(genind.ARTIF)))
#diag(distancemat.reynold.ARTIF) <- NA
heatmap.2(distancemat.reynold.ARTIF, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")
distancemat.Da.ARTIF <- as.matrix(distance.Da.ARTIF)
mean(lower(distancemat.Da.ARTIF))
rownames(distancemat.Da.ARTIF) <- as.character(unique(pop(genind.ARTIF))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.ARTIF) <- as.character(unique(pop(genind.ARTIF)))
#diag(distancemat.Da.ARTIF) <- NA
heatmap.2(distancemat.Da.ARTIF, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")


distance.Da.75pctmsng.allOC.k3.allochwe <- genet.dist(hierfstat.75pctmsng.allOC.k3, method="Da") # doesn't retain pop names, ugh!
distancemat.Da.75pctmsng.allOC.k3.allochwe <- as.matrix(distance.Da.75pctmsng.allOC.k3.allochwe)
rownames(distancemat.Da.75pctmsng.allOC.k3.allochwe) <- as.character(unique(pop(genind.75pctmsng.allOChwe.k3))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.75pctmsng.allOC.k3.allochwe) <- as.character(unique(pop(genind.75pctmsng.allOChwe.k3)))
#diag(distancemat.Da.75pctmsng.allOC.k3.allochwe) <- NA

distancemat.Da.75pctmsng.allOC.k3.allochwe[1,1] <- mean(lower(distancemat.Da.INLAND))
distancemat.Da.75pctmsng.allOC.k3.allochwe[2,2] <- mean(lower(distancemat.Da.ARTIF))
distancemat.Da.75pctmsng.allOC.k3.allochwe[3,3] <- mean(lower(distancemat.Da.COAST))


heatmap.2(distancemat.Da.75pctmsng.allOC.k3.allochwe, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")


mean(genet.dist(hierfstat.MORO, method="Da"))
hist(genet.dist(hierfstat.MORO, method="Da"), breaks=100)
mean(genet.dist(hierfstat.MORO, method="WC84"))
hist(genet.dist(hierfstat.MORO, method="WC84"), breaks=100)
mean(genet.dist(hierfstat.MORO, method="Nei87"))
hist(genet.dist(hierfstat.MORO, method="Nei87"), breaks=100)


###### genetic distances
pdf(file="geneticdistance_heatmaps_allpops_75pctmsng.pdf",
    width=15,
    height=15)

par(mfrow=c(1,1))


## Nei's Da (1983) - attributes distance to both genetic drift and mutation
distance.Da.75pctmsng.bigpops <- genet.dist(hierfstat.75pctmsng.bigpops, method="Da") # doesn't retain pop names, ugh!
distancemat.Da.75pctmsng.bigpops <- as.matrix(distance.Da.75pctmsng.bigpops)
rownames(distancemat.Da.75pctmsng.bigpops) <- as.character(unique(pop(genind.75pctmsng.bigpops))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.75pctmsng.bigpops) <- as.character(unique(pop(genind.75pctmsng.bigpops)))
#diag(distancemat.Da.75pctmsng.bigpops) <- NA

heatmap.2(distancemat.Da.75pctmsng.bigpops, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")

## Nei's Da (1983) - only on native OC populations (input for ResistanceGA)
genind.75pctmsng.nativeOC <- popsub(genind.75pctmsng.allOC, blacklist=c("FREM", "SADDLE", "IM0103", "IM679", "IM12", "IM02", "IM14", "SHOE", "TORO"), drop=TRUE)
genpop.75pctmsng.bigpops <- genind2genpop(genind.75pctmsng.bigpops)


hierfstat.75pctmsng.nativeOC <- genind2hierfstat(genind.75pctmsng.nativeOC)
distance.Da.75pctmsng.nativeOC <- genet.dist(hierfstat.75pctmsng.nativeOC, method="Da") # doesn't retain pop names, ugh!
distancemat.Da.75pctmsng.nativeOC <- as.matrix(distance.Da.75pctmsng.nativeOC)
rownames(distancemat.Da.75pctmsng.nativeOC) <- as.character(unique(pop(genind.75pctmsng.nativeOC))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.75pctmsng.nativeOC) <- as.character(unique(pop(genind.75pctmsng.nativeOC)))
#diag(distancemat.Da.75pctmsng.nativeOC) <- NA

heatmap.2(distancemat.Da.75pctmsng.nativeOC, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")

## Nei's Da on all pops
distance.Da.75pctmsng <- genet.dist(hierfstat.75pctmsng, method="Da") # doesn't retain pop names, ugh!
distancemat.Da.75pctmsng <- as.matrix(distance.Da.75pctmsng)
rownames(distancemat.Da.75pctmsng) <- as.character(unique(pop(genind.75pctmsng.allsocal))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.75pctmsng) <- as.character(unique(pop(genind.75pctmsng.allsocal)))

heatmap.2(distancemat.Da.75pctmsng, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), all pops, hclust=ward.D2, up to 75% missing data")

# Reynold's distance 1983 / approximated by genet.dist(method="Fst")(see Takezaki and Nei 1996 eq 3) - may not approximate it as closely as they say...
# or, use adegenet dist.genpop(method=3) to use the precise Reynold's distance equation - results between the two are 0.999 correlated, see below
# plotting them shows they're almost perfectly correlated, but dist.genpop returns higher values (e.g. same pair is 0.3 in hierfstat but 0.5 in dist.genpop)
# gonna use the real Reynold's for consistency's sake

# distance.reynold.75pctmsng.bigpops <- genet.dist(hierfstat.75pctmsng.bigpops, method="Fst")
# distancemat.reynold.75pctmsng.bigpops <- as.matrix(distance.reynold.75pctmsng.bigpops)
# rownames(distancemat.reynold.75pctmsng.bigpops) <- as.character(unique(pop(genind.75pctmsng.bigpops))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
# colnames(distancemat.reynold.75pctmsng.bigpops) <- as.character(unique(pop(genind.75pctmsng.bigpops)))

# heatmap.2(distancemat.reynold.75pctmsng.bigpops, 
#           col=magma(256),
#           trace="none", symm=T,
#           #hclustfun=diana)
#           hclustfun=function(x)hclust(x, method="ward.D2"),
#           main="'Reynold's distance', (1983, via hierfstat::genet.dist), hclust=ward.D2, up to 20% missing data")


## Nei's standard Ds (1972) - attributes distance to both genetic drift and mutation
distance.Ds.75pctmsng.bigpops <- genet.dist(hierfstat.75pctmsng.bigpops, method="Ds") # doesn't retain pop names, ugh!
distancemat.Ds.75pctmsng.bigpops <- as.matrix(distance.Ds.75pctmsng.bigpops)
rownames(distancemat.Ds.75pctmsng.bigpops) <- as.character(unique(pop(genind.75pctmsng.bigpops))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Ds.75pctmsng.bigpops) <- as.character(unique(pop(genind.75pctmsng.bigpops)))
diag(distancemat.Ds.75pctmsng.bigpops) <- NA

heatmap.2(distancemat.Ds.75pctmsng.bigpops, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Ds (1972, via hierfstat), hclust=ward.D2, up to 20% missing data")



distance.reynolds.75pctmsng.bigpops <- dist.genpop(genpop.75pctmsng.bigpops, method=3) # use poppr::reynolds.dist for interindividual distances
distancemat.reynolds.75pctmsng.bigpops <- as.matrix(distance.reynolds.75pctmsng.bigpops)
diag(distancemat.reynolds.75pctmsng.bigpops) <- 0 # can fill in with the intrapop distances at some point
# oddly, changing the diag from 0 to NA changes the topology of the hclust tree (for the worse, at least for Reynold's dist)
# these heatmaps are just for visualization though, so they don't really matter

heatmap.2(distancemat.reynolds.75pctmsng.bigpops, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=agnes,
          hclustfun=function(x)hclust(x, method="ward.D2"), # ward.D, ward.D2, single, complete, average, mcquitty, median, centroid; defaults to "complete"
          main="Reynold's distance (1983, via adegenet::dist.genpop), hclust=ward.D2, up to 75% missing data")

#distancegenpop.reynold.75pctmsng.bigpops <- dist.genpop(genpop.75pctmsng.bigpops, method=3)
# heatmap.2(as.matrix(distancegenpop.reynold.75pctmsng.bigpops), 
#           col=magma(256),
#           trace="none", symm=T,
#           #hclustfun=diana)
#           hclustfun=function(x)hclust(x, method="ward.D2"),
#           main="Reynold's distance (1983, via adegenet), hclust=ward.D2, up to 20% missing data")
# 
# cor(distance.reynold.75pctmsng.bigpops, distancegenpop.reynold.75pctmsng.bigpops, method="spearman") # cor = 0.99349

## Weir and Cockerman Fst (1984)
distance.wcfst.75pctmsng.bigpops <- pairwise.WCfst(hierfstat.75pctmsng.bigpops) # standard output is a matrix
distancemat.wcfst.75pctmsng.bigpops <- distance.wcfst.75pctmsng.bigpops
distance.wcfst.75pctmsng.bigpops <- as.dist(distancemat.wcfst.75pctmsng.bigpops)
diag(distancemat.wcfst.75pctmsng.bigpops) <- 0

heatmap.2(distancemat.wcfst.75pctmsng.bigpops, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=agnes,
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Weir-Cockerman Fst (1984, via hierfstat), hclust=ward.D2, up to 20% missing data")


## Jost's D (2008)
distance.jostD.75pctmsng.bigpops <- pairwise_D(genind.75pctmsng.bigpops, linearized=FALSE)
heatmap.2(as.matrix(distance.jostD.75pctmsng.bigpops), 
          col=viridis(256),
          trace="none", symm=T,
          #hclustfun=agnes,
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Jost's D (2008, via mmod), hclust=ward.D2, up to 75% missing data")
# can linearize by doing D/(1-D)
# distance.jostDlin.75pctmsng.bigpops <- pairwise_D(genind.75pctmsng.bigpops, linearized=TRUE)
# heatmap.2(as.matrix(distance.jostDlin.75pctmsng.bigpops), 
#           col=magma(256),
#           trace="none", symm=T,
#           #hclustfun=agnes,
#           hclustfun=function(x)hclust(x, method="ward.D2"),
#           main="Jost's D (2008, via mmod), 20% missing data")






df.distance.Da.75pctmsng.bigpops <- distToDataframe(distance.Da.75pctmsng.bigpops)
df.distance.reynold.75pctmsng.bigpops <- distToDataframe(distance.reynold.75pctmsng.bigpops)
df.distance.wcfst.75pctmsng.bigpops <- distToDataframe(distance.wcfst.75pctmsng.bigpops)
df.distance.jostD.75pctmsng.bigpops <- distToDataframe(distance.jostD.75pctmsng.bigpops)

df.distancemeasures <- cbind(df.distance.Da.75pctmsng.bigpops$value,
      df.distance.reynold.75pctmsng.bigpops$value,
      df.distance.wcfst.75pctmsng.bigpops$value,
      df.distance.jostD.75pctmsng.bigpops$value)
colnames(df.distancemeasures) <- c("Da", "Reynold", "WCFst", "JostD")

corrplot(cor(df.distancemeasures), method="number", number.digits=4)
cor(df.distancemeasures)

mantel.rtest(distance.Da.75pctmsng.bigpops, distance.reynold.75pctmsng.bigpops)
mantel.rtest(distance.reynold.75pctmsng.bigpops, distancegenpop.reynold.75pctmsng.bigpops)

# dist for pop coordinates
#distance.eucgeographic.bigpops <- dist(unique(genind.75pctmsng.bigpops$other$xy))
distancemat.geographic.bigpops <- spDists(as.matrix(unique(genind.75pctmsng.bigpops$other$xy)), longlat=TRUE) # uses WGS84 ellipsoid to get distances in km
distance.geographic.bigpops <- as.dist(distancemat.geographic.bigpops)

pdf(file="geneticdistance_vs_geographic_bigpops_20pctmsng.pdf",
    width=10,
    height=15)

# to linearize: dist/(1-dist)
par(mfrow=c(3,2))
plot(distance.geographic.bigpops, distance.Da.75pctmsng.bigpops,
     main="Nei's Da (1983) distance vs geographic distance", xlim=c(0,30))
abline(lm(distance.Da.75pctmsng.bigpops ~ distance.geographic.bigpops))
#plot(distance.geographic.bigpops, distance.Da.75pctmsng.bigpops/(1-distance.Da.75pctmsng.bigpops), xlim=c(0,35)) # linearized
plot(distance.geographic.bigpops, distance.reynold.75pctmsng.bigpops,
     main="Reynold's distance (1983) vs geographic distance")
abline(lm(distance.reynold.75pctmsng.bigpops ~ distance.geographic.bigpops))
#plot(distance.geographic.bigpops, distance.reynold.75pctmsng.bigpops/(1-distance.reynold.75pctmsng.bigpops), xlim=c(0,35)) # linearized
plot(distance.geographic.bigpops, distance.wcfst.75pctmsng.bigpops,
     main="Weir-Cockerman Fst (1984) vs geographic distance")
abline(lm(distance.wcfst.75pctmsng.bigpops ~ distance.geographic.bigpops))
plot(distance.geographic.bigpops, distance.wcfst.75pctmsng.bigpops/(1-distance.wcfst.75pctmsng.bigpops),
     main="Linearized Weir-Cockerman Fst (1984) vs geographic distance") # linearized
abline(lm(distance.wcfst.75pctmsng.bigpops/(1-distance.wcfst.75pctmsng.bigpops) ~ distance.geographic.bigpops))
plot(distance.geographic.bigpops, distance.jostD.75pctmsng.bigpops,
     main="Jost's D (2008) vs geographic distance")
abline(lm(distance.jostD.75pctmsng.bigpops ~ distance.geographic.bigpops))
plot(distance.geographic.bigpops, distance.jostD.75pctmsng.bigpops/(1-distance.jostD.75pctmsng.bigpops),
     main="Linearized Jost's D (2008) vs geographic distance")# linearized
abline(lm(distance.jostD.75pctmsng.bigpops/(1-distance.jostD.75pctmsng.bigpops) ~ distance.geographic.bigpops))










gl.75pctmsng.poplabel <- gl.75pctmsng
indNames(gl.75pctmsng.poplabel) <- popcoords$PopLevel2prox
bitwisedist.75pctmsng <- bitwise.dist(gl.75pctmsng.poplabel)
msn.75pctmsng <- poppr.msn(gl.75pctmsng.poplabel, distmat=bitwisedist.75pctmsng)
plot_poppr_msn(gl.75pctmsng.poplabel, msn.75pctmsng, palette=palettesub) #, cutoff=0.065, beforecut=TRUE) # refine this figure, or just make it in SplitsTree
# if a node connects to other pops before its own at a certain cutoff, could be indicative of gene flow, I think

########################### make a tree
## bitwise distance is equivalent to provesti's distance?

bitwise.euclidean <- function(x) bitwise.dist(x, percent = TRUE, mat = FALSE, missing_match = FALSE,
                                              scale_missing = FALSE, euclidean = TRUE, differences_only = FALSE,
                                              threads = 0L)

tree.75pctmsng <- aboot(gl.75pctmsng.allsocal, 
                        tree="fastme.bal",
                        distance="bitwise.dist", 
                        sample=100,
                        cutoff=50,
                        showtree=FALSE,
                        missing="ignore" #"mean" or "ignore"
)

tree.75pctmsng.euc <- aboot(gl.75pctmsng.allsocal, 
                            tree="fastme.bal",
                            distance="bitwise.euclidean", 
                            sample=100,
                            cutoff=50,
                            showtree=FALSE,
                            missing="ignore" #"mean" or "ignore"
                        )

tree.75pctmsng.meanmiss <- aboot(gl.75pctmsng.allsocal, 
                        tree="fastme.bal",
                        distance="bitwise.euclidean", 
                        sample=100,
                        cutoff=50,
                        showtree=FALSE,
                        missing="mean" #"mean" or "ignore"
)

# use if using fastme or whatever, or anything I guess
#tree.50pct.og.maf2 <- root(tree.50pct.og.maf2, outgroup=indNames(spha.genind.seppops$TEJON), resolve.root=TRUE)

# for tree=, can enter: "upgma", function(x)(upgma, method="ward"), 
# "fastme.bal" (takes awhile), "NJ", "bionj", etc

cols <- palettesub #palettelist$shiny #brewer.pal(n=nPop(gl.75pctmsng), name="Spectral")
cols <- colorRampPalette(cols)(nPop(gl.75pctmsng.allsocal))
plot.phylo(tree.75pctmsng, use.edge.length = FALSE, cex=0.5, font=2, adj=0, tip.color=cols[pop(gl.75pctmsng.allsocal)],
           type="unrooted")
nodelabels(tree.75pctmsng$node.label, adj=c(1.3, -0.5), frame="n", cex=0.5, font=3, xpd=T)
plot.phylo(tree.75pctmsng, use.edge.length = FALSE, cex=0.5, font=2, adj=0, tip.color=cols[pop(gl.75pctmsng.allsocal)],
           type="phylogram")
nodelabels(tree.75pctmsng$node.label, adj=c(1.3, -0.5), frame="n", cex=0.5, font=3, xpd=T)

legend('topleft', 
       legend=popNames(gl.75pctmsng),
       fill=cols,
       border=F,
       bty="n",
       cex=0.3)
axis(side=1)
title(xlab="Genetic distance (proportion of loci that are different/bitwise.dist), fastme.bal tree")

