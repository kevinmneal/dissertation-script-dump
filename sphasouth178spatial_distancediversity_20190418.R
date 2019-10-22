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
library(dplyr)
library(PopGenReport)


vcf.OCreduced148.radiator.50pctmsng <- vcfR2genlight("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85_radiator_final/OCreduced148_radiator_mac3/19_radiator_genomic_converter_20190322@2235/OCreduced148.radiator.75pctmsng.oneSNPmac3.vcf")


vcf.75pctmsng.OCreduced148 <- read.vcfR("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/OCreduced148_85clust_75pctmsng_oneSNPmac2_20190205.recode.vcf")
vcf.75pctmsng.OCall208 <- read.vcfR("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/OCall208_85clust_75pctmsng_oneSNPmac2_20190205.recode.vcf")

#genind.75pctmsng.OCreduced148 <- vcfR2genind(vcf.75pctmsng.OCreduced148)

popcoords.OCreduced148 <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/OCreduced148_85clust_popcoords_20190205.csv")
popcoords.OCall208 <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/OCall208_85clust_popcoords_20190205.csv")
saveRDS(popcoords.OCreduced148, file="G:/My Drive/WORKING FOLDER Jan 2015/STUFF THAT ISNT NEWCAL/Spea/RAD_data/allsocal_303indv/popcoords.OCreduced148.rds")

genind.75pctmsng.OCall208 <- vcfR2genind(vcf.75pctmsng.OCall208)
pop(genind.75pctmsng.OCall208) <- popcoords.OCall208$Site
genind.75pctmsng.OCall208.k3 <- genind.75pctmsng.OCall208
pop(genind.75pctmsng.OCall208.k3) <- popcoords.OCall208$K3hier

gl.75pctmsng.OCall208 <- vcfR2genlight(vcf.75pctmsng.OCall208)
ploidy(gl.75pctmsng.OCall208) <- 2
pop(gl.75pctmsng.OCall208) <- popcoords.OCall208$SiteIMcondensed
gl.75pctmsng.OCall208.k3 <- gl.75pctmsng.OCall208
ploidy(gl.75pctmsng.OCall208.k3) <- 2
pop(gl.75pctmsng.OCall208.k3) <- popcoords.OCall208$K3hier

genind.75pctmsng.OCreduced148 <- vcfR2genind(vcf.75pctmsng.OCreduced148)
pop(genind.75pctmsng.OCreduced148) <- popcoords.OCreduced148$SiteIMcondensed
genind.75pctmsng.OCreduced148.k3 <- genind.75pctmsng.OCreduced148
pop(genind.75pctmsng.OCreduced148.k3) <- popcoords.OCreduced148$K3hier

gl.75pctmsng.OCreduced148 <- vcfR2genlight(vcf.75pctmsng.OCreduced148)
ploidy(gl.75pctmsng.OCreduced148) <- 2
pop(gl.75pctmsng.OCreduced148) <- popcoords.OCreduced148$SiteIMcondensed
gl.75pctmsng.OCreduced148.k3 <- gl.75pctmsng.OCreduced148
ploidy(gl.75pctmsng.OCreduced148.k3) <- 2
pop(gl.75pctmsng.OCreduced148.k3) <- popcoords.OCreduced148$K3hier



privatealleles.75pctmsng.OCreduced148.k3 <- poppr::private_alleles(genind.75pctmsng.OCreduced148.k3,
                                                   alleles ~ .,
                                                   report="table",
                                                   level="population",
                                                   count.alleles=FALSE,
                                                   drop=TRUE)
privatealleles.counts.75pctmsng.OCreduced148.k3 <- rowSums(privatealleles.75pctmsng.OCreduced148.k3, na.rm=TRUE)

propshared.75pctmsng.OCall208 <- pairwise.propShared(genind.75pctmsng.OCall208)
propshared.75pctmsng.OCreduced148 <- pairwise.propShared(genind.75pctmsng.OCreduced148)
propshared.75pctmsng.OCreduced148.k3 <- pairwise.propShared(genind.75pctmsng.OCreduced148.k3)

privatealleles.75pctmsng.OCall208 <- poppr::private_alleles(genind.75pctmsng.OCall208,
                                                                   alleles ~ .,
                                                                   report="table",
                                                                   level="population",
                                                                   count.alleles=TRUE,
                                                                   drop=TRUE)
privatealleles.counts.75pctmsng.OCall208 <- rowSums(privatealleles.75pctmsng.OCall208)
# Get percentages.
sweep(privatealleles.counts.75pctmsng.OCall208, 2, nAll(genind.75pctmsng.OCall208)[colnames(privatealleles.counts.75pctmsng.OCall208)], FUN = "/")

unique(pop(genind.75pctmsng.OCall208))
genind.75pctmsng.OCall208.bigpops <- popsub(genind.75pctmsng.OCall208, blacklist=c("SADDLE", "TORO"))

propshared.75pctmsng.OCall208.bigpops <- pairwise.propShared(genind.75pctmsng.OCall208.bigpops)
propshared.75pctmsng.OCreduced148 <- pairwise.propShared(genind.75pctmsng.OCreduced148)
propshared.75pctmsng.OCall208.k3 <- pairwise.propShared(genind.75pctmsng.OCall208.k3)


#### Genetic diversity
# Might just use DnaSP 6 - nevermind, I don't trust it... numbers don't make sense

hierfstat.75pctmsng.OCreduced148 <- genind2hierfstat(genind.75pctmsng.OCreduced148)
basicstats.75pctmsng.OCreduced148 <- basic.stats(hierfstat.75pctmsng.OCreduced148)
ho.pops.75pctmsng.OCreduced148 <- colMeans(basicstats.75pctmsng.OCreduced148$Ho, na.rm=T)
#he_adegenet.pops.75pctmsng.OCreduced148 <- Hs(genind.75pctmsng.OCreduced148) 
hs.pops.75pctmsng.OCreduced148 <- colMeans(basicstats.75pctmsng.OCreduced148$Hs, na.rm=T)
fis.pops.75pctmsng.OCreduced148 <- colMeans(basicstats.75pctmsng.OCreduced148$Fis, na.rm=T) # they're all negative... does this make sense?
inbreeding.pops.75pctmsng.OCreduced148 <- inbreeding(genind.75pctmsng.OCreduced148, pop=pop(genind.75pctmsng.OCreduced148), res.type="estimate") # adegenet's way of calculating it; doesn't give pop estimates?
ar.75pctmsng.OCreduced148 <- allelic.richness(hierfstat.75pctmsng.OCreduced148)
ar.pops.75pctmsng.OCreduced148 <- colMeans(ar.75pctmsng.OCreduced148$Ar, na.rm=T)

#hohe.pops.75pctmsng.OCreduced148 <- cbind(ho.pops.75pctmsng.OCreduced148, he_adegenet.pops.75pctmsng.OCreduced148)
arhohsfis.pops.75pctmsng.OCreduced148 <- cbind(ar.pops.75pctmsng.OCreduced148, ho.pops.75pctmsng.OCreduced148, hs.pops.75pctmsng.OCreduced148, fis.pops.75pctmsng.OCreduced148) # this is "Nei's gene diversity", which is arguably superior to standard Hexp - accounts for sampling error?
# in any case, both He and Hs give similar results. Ho in almost all cases is > He or Hs...

colnames(arhohsfis.pops.75pctmsng.OCreduced148) <- c("Ar", "Ho", "He", "Fis")
#summarise(popcoords$PopLevel2_prox)
popcounts <- count(popcoords.OCreduced148, SiteIMcondensed)

arhohsfissamplesize.pops.75pctmsng.OCreduced148 <- cbind(popcounts, arhohsfis.pops.75pctmsng.OCreduced148[order(row.names(arhohsfis.pops.75pctmsng.OCreduced148)),])
View(arhohsfissamplesize.pops.75pctmsng.OCreduced148)
write.csv(arhohsfissamplesize.pops.75pctmsng.OCreduced148, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/arhohsfis_75pctmsng.OCreduced148.csv")


### for the 3 clusters
hierfstat.75pctmsng.OCreduced148.k3 <- genind2hierfstat(genind.75pctmsng.OCreduced148.k3)
basicstats.75pctmsng.OCreduced148.k3 <- basic.stats(hierfstat.75pctmsng.OCreduced148.k3)
ho.pops.75pctmsng.OCreduced148.k3 <- colMeans(basicstats.75pctmsng.OCreduced148.k3$Ho, na.rm=T)
#he_adegenet.pops.75pctmsng.OCreduced148.k3 <- Hs(genind.75pctmsng.OCreduced148.k3) 
hs.pops.75pctmsng.OCreduced148.k3 <- colMeans(basicstats.75pctmsng.OCreduced148.k3$Hs, na.rm=T)
fis.pops.75pctmsng.OCreduced148.k3 <- colMeans(basicstats.75pctmsng.OCreduced148.k3$Fis, na.rm=T) # they're all negative... does this make sense?
inbreeding.pops.75pctmsng.OCreduced148.k3 <- inbreeding(genind.75pctmsng.OCreduced148.k3, pop=pop(genind.75pctmsng.OCreduced148.k3), res.type="estimate") # adegenet's way of calculating it; doesn't give pop estimates?
ar.75pctmsng.OCreduced148.k3 <- allelic.richness(hierfstat.75pctmsng.OCreduced148.k3)
ar.pops.75pctmsng.OCreduced148.k3 <- colMeans(ar.75pctmsng.OCreduced148.k3$Ar, na.rm=T)

#hohe.pops.75pctmsng.OCreduced148.k3 <- cbind(ho.pops.75pctmsng.OCreduced148.k3, he_adegenet.pops.75pctmsng.OCreduced148.k3)
arhohsfis.pops.75pctmsng.OCreduced148.k3 <- cbind(ar.pops.75pctmsng.OCreduced148.k3, ho.pops.75pctmsng.OCreduced148.k3, hs.pops.75pctmsng.OCreduced148.k3, fis.pops.75pctmsng.OCreduced148.k3) # this is "Nei's gene diversity", which is arguably superior to standard Hexp - accounts for sampling error?
# in any case, both He and Hs give similar results. Ho in almost all cases is > He or Hs...

colnames(arhohsfis.pops.75pctmsng.OCreduced148.k3) <- c("Ar", "Ho", "He", "Fis")
#summarise(popcoords$PopLevel2_prox)
popcounts <- count(popcoords.OCreduced148, K3hier)

arhohsfissamplesize.pops.75pctmsng.OCreduced148.k3 <- cbind(popcounts, arhohsfis.pops.75pctmsng.OCreduced148.k3[order(row.names(arhohsfis.pops.75pctmsng.OCreduced148.k3)),])
View(arhohsfissamplesize.pops.75pctmsng.OCreduced148.k3)
write.csv(arhohsfissamplesize.pops.75pctmsng.OCreduced148.k3, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/arhohsfis_75pctmsng.OCreduced148.k3.csv")


hierfstat.75pctmsng.OCall208 <- genind2hierfstat(genind.75pctmsng.OCall208)
basicstats.75pctmsng.OCall208 <- basic.stats(hierfstat.75pctmsng.OCall208)
ho.pops.75pctmsng.OCall208 <- colMeans(basicstats.75pctmsng.OCall208$Ho, na.rm=T)
#he_adegenet.pops.75pctmsng.OCall208 <- Hs(genind.75pctmsng.OCall208) 
hs.pops.75pctmsng.OCall208 <- colMeans(basicstats.75pctmsng.OCall208$Hs, na.rm=T)
fis.pops.75pctmsng.OCall208 <- colMeans(basicstats.75pctmsng.OCall208$Fis, na.rm=T) # they're all negative... does this make sense?
#inbreeding.pops.75pctmsng.OCall208 <- inbreeding(genind.75pctmsng.OCall208, pop=pop(genind.75pctmsng.OCall208), res.type="estimate") # adegenet's way of calculating it; doesn't give pop estimates?
ar.75pctmsng.OCall208 <- allelic.richness(hierfstat.75pctmsng.OCall208)
ar.pops.75pctmsng.OCall208 <- colMeans(ar.75pctmsng.OCall208$Ar, na.rm=T)

#hohe.pops.75pctmsng.OCall208 <- cbind(ho.pops.75pctmsng.OCall208, he_adegenet.pops.75pctmsng.OCall208)
arhohsfis.pops.75pctmsng.OCall208 <- cbind(ar.pops.75pctmsng.OCall208, ho.pops.75pctmsng.OCall208, hs.pops.75pctmsng.OCall208, fis.pops.75pctmsng.OCall208) # this is "Nei's gene diversity", which is arguably superior to standard Hexp - accounts for sampling error?
# in any case, both He and Hs give similar results. Ho in almost all cases is > He or Hs...

colnames(arhohsfis.pops.75pctmsng.OCall208) <- c("Ar", "Ho", "He", "Fis")
#summarise(popcoords$PopLevel2_prox)
popcounts <- count(popcoords.OCall208, Site)

arhohsfissamplesize.pops.75pctmsng.OCall208 <- cbind(popcounts, arhohsfis.pops.75pctmsng.OCall208[order(row.names(arhohsfis.pops.75pctmsng.OCall208)),])
View(arhohsfissamplesize.pops.75pctmsng.OCall208)
write.csv(arhohsfissamplesize.pops.75pctmsng.OCall208, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/arhohsfis_75pctmsng.OCall208_IMsites.csv")

# try fis_summary in radiator, or use ppfis in hierfstat
# dangit, fis_summary gives locus-level specs so gotta do colMeans again...
OCall208.75pctmsng.fissummary <- fis_summary(gl.75pctmsng.OCall208) # took 806 seconds on laptop
OCall208.75pctmsng.fissummary$fis.summary


boot.ppfis
boot.ppfst

OCall208.gendiv <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/arhohsfis_75pctmsng.OCall208_IMsites_outliersremoved.csv",
                            header=T,
                            stringsAsFactors = T)
View(OCall208.gendiv)
pairs(OCall208.gendiv)



#par(mfrow=c(1,1))
OCall208.gendiv <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/stacks_populations/OCall208.75pctmsng.fullsumstats_evensmallersummary_forcorrelationsinR_20190311.tsv.csv",
                            header=T,
                            stringsAsFactors = T)
View(OCall208.gendiv)
pairs(OCall208.gendiv)

svg("G:/My Drive/Manuscripts/Spea_socal_OC_thesischapter/figures_20190211/OCall208_gendiv_correlations_20190218.svg",
    width=14,
    height=21,
    pointsize=12,
    antialias = "default",
    onefile=TRUE)
cairo_pdf("G:/My Drive/Manuscripts/Spea_socal_OC_thesischapter/figures_20190211/OCall208_gendiv_correlations_20190311.pdf",
          width=14,
          height=21,
          pointsize=12,
          antialias="default",
          fallback_resolution = 600,
          onefile=TRUE)
par(mfrow=c(3,2),
    mar=c(1,1,1,1))
corm.spearman <- cor.mtest(OCall208.gendiv[,-c(1,2,3,4)], 
                      alternative="two.sided", # consider using "greater" to find significant *positive* associations
                      method="spearman", 
                      exact=FALSE)
corm.kendall <- cor.mtest(OCall208.gendiv[,-c(1,2,3,4)], 
                           alternative="two.sided", # consider using "greater" to find significant *positive* associations
                           method="kendall", 
                           exact=FALSE)
OCall208.cor.spearman <- cor(OCall208.gendiv[,-c(1,2,3,4)], use="pairwise.complete.obs", method="spearman")
OCall208.cor.spearman[corm.spearman$p > 0.05] <- 0
corrplot(OCall208.cor.spearman, 
         title="Genetic diversity correlations, all sites",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
         #p.mat=corm.spearman$p, insig="pch")
OCall208.cor.kendall <- cor(OCall208.gendiv[,-c(1,2,3,4)], use="pairwise.complete.obs", method="kendall")
OCall208.cor.kendall[corm.kendall$p > 0.05] <- 0
corrplot(OCall208.cor.kendall, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, all sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9)
         #p.mat=corm.kendall$p, insig="pch")

corm.spearman.ARTIF <- cor.mtest(OCall208.gendiv[OCall208.gendiv$Artifnatural=="ARTIF",-c(1,2,3,4)], 
                           alternative="two.sided", # consider using "greater" to find significant *positive* associations
                           method="spearman", 
                           exact=FALSE)
corm.kendall.ARTIF <- cor.mtest(OCall208.gendiv[OCall208.gendiv$Artifnatural=="ARTIF",-c(1,2,3,4)], 
                          alternative="two.sided", # consider using "greater" to find significant *positive* associations
                          method="kendall", 
                          exact=FALSE)
OCall208.cor.spearman.ARTIF <- cor(OCall208.gendiv[OCall208.gendiv$Artifnatural=="ARTIF",-c(1,2,3,4)], use="pairwise.complete.obs", method="spearman")
OCall208.cor.spearman.ARTIF[corm.spearman.ARTIF$p > 0.05] <- 0
corrplot(OCall208.cor.spearman.ARTIF, 
         title="Genetic diversity correlations, artificial sites",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
OCall208.cor.kendall.ARTIF <- cor(OCall208.gendiv[OCall208.gendiv$Artifnatural=="ARTIF",-c(1,2,3,4)], use="pairwise.complete.obs", method="kendall")
OCall208.cor.kendall.ARTIF[corm.kendall.ARTIF$p > 0.05] <- 0
corrplot(OCall208.cor.kendall.ARTIF, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, all sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9)
#p.mat=corm.kendall$p, insig="pch")

corm.spearman.NATURAL <- cor.mtest(OCall208.gendiv[OCall208.gendiv$Artifnatural=="NATURAL",-c(1,2,3,4)], 
                                 alternative="two.sided", # consider using "greater" to find significant *positive* associations
                                 method="spearman", 
                                 exact=FALSE)
corm.kendall.NATURAL <- cor.mtest(OCall208.gendiv[OCall208.gendiv$Artifnatural=="NATURAL",-c(1,2,3,4)], 
                                alternative="two.sided", # consider using "greater" to find significant *positive* associations
                                method="kendall", 
                                exact=FALSE)
OCall208.cor.spearman.NATURAL <- cor(OCall208.gendiv[OCall208.gendiv$Artifnatural=="NATURAL",-c(1,2,3,4)], use="pairwise.complete.obs", method="spearman")
OCall208.cor.spearman.NATURAL[corm.spearman.NATURAL$p > 0.05] <- 0
corrplot(OCall208.cor.spearman.NATURAL, 
         title="Genetic diversity correlations, all natural sites",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
OCall208.cor.kendall.NATURAL <- cor(OCall208.gendiv[OCall208.gendiv$Artifnatural=="NATURAL",-c(1,2,3,4)], use="pairwise.complete.obs", method="kendall")
OCall208.cor.kendall.NATURAL[corm.kendall.NATURAL$p > 0.05] <- 0
corrplot(OCall208.cor.kendall.NATURAL, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, all sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9)
#p.mat=corm.kendall$p, insig="pch")

corm.spearman.COAST <- cor.mtest(OCall208.gendiv[OCall208.gendiv$ClusterK4=="Coast",-c(1,2,3,4)], 
                                 alternative="two.sided", # consider using "greater" to find significant *positive* associations
                                 method="spearman", 
                                 exact=FALSE)
corm.kendall.COAST <- cor.mtest(OCall208.gendiv[OCall208.gendiv$ClusterK4=="Coast",-c(1,2,3,4)], 
                                alternative="two.sided", # consider using "greater" to find significant *positive* associations
                                method="kendall", 
                                exact=FALSE)
OCall208.cor.spearman.COAST <- cor(OCall208.gendiv[OCall208.gendiv$ClusterK4=="Coast",-c(1,2,3,4)], use="pairwise.complete.obs", method="spearman")
OCall208.cor.spearman.COAST[corm.spearman.COAST$p > 0.05] <- 0
corrplot(OCall208.cor.spearman.COAST, 
         title="Genetic diversity correlations, Coast sites",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
OCall208.cor.kendall.COAST <- cor(OCall208.gendiv[OCall208.gendiv$ClusterK4=="Coast",-c(1,2,3,4)], use="pairwise.complete.obs", method="kendall")
OCall208.cor.kendall.COAST[corm.kendall$p > 0.05] <- 0
corrplot(OCall208.cor.kendall.COAST, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, Inland sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9)

corm.spearman.INLAND <- cor.mtest(OCall208.gendiv[OCall208.gendiv$ClusterK4=="Inlandnatural",-c(1,2,3,4)], 
                                 alternative="two.sided", # consider using "greater" to find significant *positive* associations
                                 method="spearman", 
                                 exact=FALSE)
corm.kendall.INLAND <- cor.mtest(OCall208.gendiv[OCall208.gendiv$ClusterK4=="Inlandnatural",-c(1,2,3,4)], 
                                alternative="two.sided", # consider using "greater" to find significant *positive* associations
                                method="kendall", 
                                exact=FALSE)
OCall208.cor.spearman.INLAND <- cor(OCall208.gendiv[OCall208.gendiv$ClusterK4=="Inlandnatural",-c(1,2,3,4)], use="pairwise.complete.obs", method="spearman")
OCall208.cor.spearman.INLAND[corm.spearman.INLAND$p > 0.05] <- 0
corrplot(OCall208.cor.spearman.INLAND, 
         title="Genetic diversity correlations, natural Inland sites",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
OCall208.cor.kendall.INLAND <- cor(OCall208.gendiv[OCall208.gendiv$ClusterK4=="Inlandnatural",-c(1,2,3,4)], use="pairwise.complete.obs", method="kendall")
OCall208.cor.kendall.INLAND[corm.kendall$p > 0.05] <- 0
corrplot(OCall208.cor.kendall.INLAND, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, Inland sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9)


dev.off()


OCreduced148.choosek.aic <- snapclust.choose.k(30, genind.75pctmsng.OCreduced148, IC=AIC)
OCreduced148.choosek.aicc <- snapclust.choose.k(30, genind.75pctmsng.OCreduced148, IC=AICc)
OCreduced148.choosek.bic <- snapclust.choose.k(20, genind.75pctmsng.OCreduced148, IC=BIC)
OCreduced148.choosek.kic <- snapclust.choose.k(20, genind.75pctmsng.OCreduced148, IC=KIC)

plot(OCreduced148.choosek.aic) # allsocal: best K=3, 10, or 14; allOC: 3, 8, 10, 11, 16
plot(OCreduced148.choosek.aicc)
plot(OCreduced148.choosek.bic) # allsocal: K=3; also 3 for allOC
plot(OCreduced148.choosek.kic)

kmeans.75pctmsng.OCreduced148 <- find.clusters(gl.75pctmsng.OCreduced148)
dapc.75pctmsng.OCreduced148 <- dapc(gl.75pctmsng.OCreduced148, kmeans.75pctmsng.OCreduced148$grp$grp)
compoplot(dapc.75pctmsng.OCreduced148, show.lab=TRUE, lab=pop(gl.75pctmsng.OCreduced148), legend=FALSE, main="DAPC, 303 individuals")
scatter(dapc.75pctmsng.OCreduced148, scree.da=FALSE)
temp <- optim.a.score(dapc.75pctmsng.OCreduced148) # optimal looks to be 13 PCs




ppfisboot.75pctmsng.OCall208 <- boot.ppfis(hierfstat.75pctmsng.OCall208, nboot=1000) # only gives confidence intervals...
# run with no bootstraps to get a point estimate. What does it mean if it's outside the CI???
ppfisboot.75pctmsng.OCall208.boot100 <- boot.ppfis(hierfstat.75pctmsng.OCall208, nboot=1)

ppfisboot.75pctmsng.OCreduced148.k3.1boot <- boot.ppfis(hierfstat.75pctmsng.OCreduced148.k3, nboot=1)
ppfisboot.75pctmsng.OCreduced148.k3.1000boot <- boot.ppfis(hierfstat.75pctmsng.OCreduced148.k3, nboot=1000)

popprdiv.75pctmsng.OCreduced148.k3 <- poppr(genind.75pctmsng.OCreduced148.k3) # turns out, these numbers basically just reflect sample size. Useless!
View(popprdiv.75pctmsng.OCreduced148.k3)
# don't bother going further with this.

# adegenet popgen summary
adegenetsummary.OCall208.75pctmsng <- summary(genind.75pctmsng.OCall208)
adegenetbartlett.OCall208.75pctmsng <- bartlett.test(list(adegenetsummary.OCall208.75pctmsng$Hexp, adegenetsummary.OCall208.75pctmsng$Hobs))

# global fst
globalfst.75pctmsng.OCreduced148 <- fstat(genind.75pctmsng.OCreduced148) # global Fst for all 208 samples is 0.2752; Fit is 0.2561
globalfst.75pctmsng.OCreduced148.k3 <- fstat(genind.75pctmsng.OCreduced148.k3) # global Fst for all 208 samples is 0.2752; Fit is 0.2561
# wc estimates of global fst/fis
wc.75pctmsng.OCreduced148 <- wc(hierfstat.75pctmsng.OCreduced148) # FST: 0.2467; FIS: 0.0121. Same results for fstat function...
wc.75pctmsng.OCreduced148.k3 <- wc(hierfstat.75pctmsng.OCreduced148.k3)


vcf.75pctmsng.OCartif79 <- read.vcfR("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/OCartif79_85clust_75pctmsng_oneSNPmac2_20190205.recode.vcf")
genind.75pctmsng.OCartif79 <- vcfR2genind(vcf.75pctmsng.OCartif79)
popcoords.OCartif79 <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/OCartif79_85clust_popcoords_20190205.csv")
pop(genind.75pctmsng.OCartif79) <- popcoords.OCartif79$Site

#genind.75pctmsng.OCartif <- popsub(genind.75pctmsng.OCall208, sublist=c("IM01", "IM02", "IM03", "IM06", "IM07", "IM09", "IM12", "IM14", "SHOE"))
globalfst.75pctmsng.OCartif79 <- fstat(genind.75pctmsng.OCartif79) # global Fst for all 208 samples is 0.2752; Fit is 0.2561


#### Pairwise Fst (Weir-Cockerman 1984) and Nei's Da (1983), and maybe Jost's D (2008)
# remove pops with 1 or otherwise very few individuals - gets weird results. use poppr popsub

## make subset of populations with >3 individuals

#fst.75pctmsng.k3 <- fstat(genind.75pctmsng.allOChwe.k3)
#pwfst.75pctmsng.k3 <- pairwise.WCfst(hierfstat.75pctmsng.allOC.k3)

### intrapop genetic diversity (by taking mean or median genetic distances among individuals within a pop)
# genet.dist reads the pop() info; would have to make each individual its own population
# make a function... note, WCfst doesn't work
# easiest way is probably to run on the full 303 individuals and use dplyr or something to find means by population


genind.75pctmsng.OCreduced148.k3.seppops <- seppop(genind.75pctmsng.OCreduced148.k3)
genind.INLAND <- genind.75pctmsng.OCreduced148.k3.seppops$INLAND
genind.COAST <- genind.75pctmsng.OCreduced148.k3.seppops$COAST
genind.SANTENA <- genind.75pctmsng.OCreduced148.k3.seppops$SANTENA

gl.75pctmsng.OCreduced148.k3.seppops <- seppop(gl.75pctmsng.OCreduced148.k3)
gl.INLAND <- gl.75pctmsng.OCreduced148.k3.seppops$INLAND
gl.COAST <- gl.75pctmsng.OCreduced148.k3.seppops$COAST
gl.SANTENA <- gl.75pctmsng.OCreduced148.k3.seppops$SANTENA

pop(genind.INLAND) <- indNames(genind.INLAND)
hierfstat.INLAND <- genind2hierfstat(genind.INLAND)
#par(mfrow=c(2,2))
distance.Da.INLAND <- genet.dist(hierfstat.INLAND, method="Da") 
mean(distance.Da.INLAND) 
hist(distance.Da.INLAND, breaks=100)
distancemat.Da.INLAND <- as.matrix(distance.Da.INLAND)
rownames(distancemat.Da.INLAND) <- as.character(unique(pop(genind.INLAND))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.INLAND) <- as.character(unique(pop(genind.INLAND)))
#diag(distancemat.Da.INLAND) <- NA
heatmap.2(distancemat.Da.INLAND, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")

pop(gl.INLAND) <- indNames(gl.INLAND)
bitwise.ia.INLAND <- bitwise.ia(gl.INLAND)
bitwise.dist.euc.INLAND <- bitwise.dist(gl.INLAND, euclidean=TRUE)
pop(genind.COAST) <- indNames(genind.COAST)
hierfstat.COAST <- genind2hierfstat(genind.COAST)
distance.Da.COAST <- genet.dist(hierfstat.COAST, method="Da") # equivalent (mostly) to Reynold's distance
mean(distance.Da.COAST) 
hist(distance.Da.COAST, breaks=100)
distancemat.Da.COAST <- as.matrix(distance.Da.COAST)
rownames(distancemat.Da.COAST) <- as.character(unique(pop(genind.COAST))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.COAST) <- as.character(unique(pop(genind.COAST)))
#diag(distancemat.Da.COAST) <- NA
bionj.COAST <- nj(distancemat.Da.COAST)
plot(bionj.COAST)
bitwise.ia.COAST <- bitwise.ia()

pop(gl.COAST) <- indNames(gl.COAST)
bitwise.ia.COAST <- bitwise.ia(gl.COAST)
bitwise.dist.euc.COAST <- bitwise.dist(gl.COAST, euclidean=TRUE)
heatmap.2(distancemat.Da.COAST, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=function(x)agnes(x, diss=TRUE, metric="euclidean", method="ward"),
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")

pop(genind.SANTENA) <- indNames(genind.SANTENA)
hierfstat.SANTENA <- genind2hierfstat(genind.SANTENA)
distance.Da.SANTENA <- genet.dist(hierfstat.SANTENA, method="Da") # equivalent (mostly) to Reynold's distance
mean(distance.Da.SANTENA)
#mean(lower(distance.Da.SANTENA))
hist(distance.Da.SANTENA, breaks=100)
distancemat.Da.SANTENA <- as.matrix(distance.Da.SANTENA)
#mean(lower(distancemat.Da.SANTENA))
rownames(distancemat.Da.SANTENA) <- as.character(unique(pop(genind.SANTENA))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.SANTENA) <- as.character(unique(pop(genind.SANTENA)))
#diag(distancemat.Da.SANTENA) <- NA
heatmap.2(distancemat.Da.SANTENA, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")

genind.75pctmsng.OCreduced148.indiv <- genind.75pctmsng.OCreduced148
pop(genind.75pctmsng.OCreduced148.indiv) <- indNames(genind.75pctmsng.OCreduced148.indiv)
hierfstat.75pctmsng.OCreduced148.indiv <- genind2hierfstat(genind.75pctmsng.OCreduced148.indiv)
distance.Da.75pctmsng.OCreduced148.indiv <- genet.dist(hierfstat.75pctmsng.OCreduced148.indiv, method="Da") # doesn't retain pop names, ugh!
distancemat.Da.75pctmsng.OCreduced148.indiv <- as.matrix(distance.Da.75pctmsng.OCreduced148.indiv)
rownames(distancemat.Da.75pctmsng.OCreduced148.indiv) <- as.character(pop(genind.75pctmsng.OCreduced148.indiv)) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.75pctmsng.OCreduced148.indiv) <- as.character(pop(genind.75pctmsng.OCreduced148.indiv))
#diag(distancemat.Da.75pctmsng.OCreduced148.k3) <- NA

distancemat.Da.75pctmsng.OCreduced148.k3["INLAND","INLAND"] <- mean(lower(distancemat.Da.INLAND))
distancemat.Da.75pctmsng.OCreduced148.k3["COAST","COAST"] <- mean(lower(distancemat.Da.COAST))
distancemat.Da.75pctmsng.OCreduced148.k3["SANTENA","SANTENA"] <- mean(lower(distancemat.Da.SANTENA))


heatmap.2(distancemat.Da.75pctmsng.OCreduced148.indiv, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")

stamppPhylip(distancemat.Da.75pctmsng.OCreduced148.indiv, 
             file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/OCreduced148_NeiDa_75pctmsng_indivdist_popnames.phy.dst")
# splitstree doesn't like duplicate names for individuals, so can't name them as their populations...

distance.X2.75pctmsng.OCreduced148.indiv <- genet.dist(hierfstat.75pctmsng.OCreduced148.indiv, method="X2") # doesn't retain pop names, ugh!
distancemat.X2.75pctmsng.OCreduced148.indiv <- as.matrix(distance.X2.75pctmsng.OCreduced148.indiv)
rownames(distancemat.X2.75pctmsng.OCreduced148.indiv) <- as.character(unique(pop(genind.75pctmsng.OCreduced148.indiv))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.X2.75pctmsng.OCreduced148.indiv) <- as.character(unique(pop(genind.75pctmsng.OCreduced148.indiv)))
#diag(distancemat.X2.75pctmsng.OCreduced148.k3) <- NA
heatmap.2(distancemat.X2.75pctmsng.OCreduced148.indiv, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")
stamppPhylip(distancemat.X2.75pctmsng.OCreduced148.indiv, 
             file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/OCreduced148_X2_75pctmsng_indivdist.phy.dst")
stamppPhylip(as.matrix(distance.jostD.75pctmsng.OCreduced148), 
             file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/OCreduced148_JostD_75pctmsng_dist.phy.dst")

genind.75pctmsng.OCall208.indiv <- genind.75pctmsng.OCall208
pop(genind.75pctmsng.OCall208.indiv) <- indNames(genind.75pctmsng.OCall208.indiv)
hierfstat.75pctmsng.OCall208.indiv <- genind2hierfstat(genind.75pctmsng.OCall208.indiv)
distance.Da.75pctmsng.OCall208.indiv <- genet.dist(hierfstat.75pctmsng.OCall208.indiv, method="Da") # doesn't retain pop names, ugh!
distancemat.Da.75pctmsng.OCall208.indiv <- as.matrix(distance.Da.75pctmsng.OCall208.indiv)
rownames(distancemat.Da.75pctmsng.OCall208.indiv) <- as.character(unique(pop(genind.75pctmsng.OCall208.indiv))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.75pctmsng.OCall208.indiv) <- as.character(unique(pop(genind.75pctmsng.OCall208.indiv)))
#diag(distancemat.Da.75pctmsng.OCall208.k3) <- NA

#distancemat.Da.75pctmsng.OCall208.k3[1,1] <- mean(lower(distancemat.Da.INLAND))
#distancemat.Da.75pctmsng.OCall208.k3[2,2] <- mean(lower(distancemat.Da.SANTENA))
#distancemat.Da.75pctmsng.OCall208.k3[3,3] <- mean(lower(distancemat.Da.COAST))


heatmap.2(distancemat.Da.75pctmsng.OCall208.indiv, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")


hierfstat.75pctmsng.OCreduced148.k3 <- genind2hierfstat(genind.75pctmsng.OCreduced148.k3)
distance.Da.75pctmsng.OCreduced148.k3 <- genet.dist(hierfstat.75pctmsng.OCreduced148.k3, method="Da") # doesn't retain pop names, ugh!
distancemat.Da.75pctmsng.OCreduced148.k3 <- as.matrix(distance.Da.75pctmsng.OCreduced148.k3)
rownames(distancemat.Da.75pctmsng.OCreduced148.k3) <- as.character(unique(pop(genind.75pctmsng.OCreduced148.k3))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.75pctmsng.OCreduced148.k3) <- as.character(unique(pop(genind.75pctmsng.OCreduced148.k3)))
#diag(distancemat.Da.75pctmsng.OCreduced148.k3) <- NA

#distancemat.Da.75pctmsng.OCreduced148.k3[1,1] <- mean(lower(distancemat.Da.INLAND))
#distancemat.Da.75pctmsng.OCreduced148.k3[2,2] <- mean(lower(distancemat.Da.SANTENA))
#distancemat.Da.75pctmsng.OCreduced148.k3[3,3] <- mean(lower(distancemat.Da.COAST))


heatmap.2(distancemat.Da.75pctmsng.OCreduced148.k3, 
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


####### amova
#bitwise.dist.euc.75pctmsng.OCreduced148.indiv <- bitwise.dist(gl.75pctmsng.OCreduced148, euclidean=TRUE)
strata(gl.75pctmsng.OCreduced148) <- data.frame(popcoords.OCreduced148[,c("SiteIMcondensed", "K3hier")])
strata(gl.75pctmsng.OCreduced148)
strata(genind.75pctmsng.OCreduced148) <- data.frame(popcoords.OCreduced148[,c("SiteIMcondensed", "K3hier")])
strata(genind.75pctmsng.OCreduced148, formula=~SiteIMcondensed/K3hier, combine=FALSE)

amova.75pctmsng.OCreduced148 <- poppr.amova(genind.75pctmsng.OCreduced148, hier=~SiteIMcondensed/K3hier, 
                                            #dist=bitwise.dist.euc.75pctmsng.OCreduced148.indiv,
                                            within=FALSE,
                                            method="ade4",
                                            missing="mean")



#######





## Nei's Da (1983) - attributes distance to both genetic drift and mutation
distance.Da.75pctmsng.OCreduced148 <- genet.dist(hierfstat.75pctmsng.OCreduced148, method="Da") # doesn't retain pop names, ugh!
distancemat.Da.75pctmsng.OCreduced148 <- as.matrix(distance.Da.75pctmsng.OCreduced148)
rownames(distancemat.Da.75pctmsng.OCreduced148) <- as.character(unique(pop(genind.75pctmsng.OCreduced148))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.75pctmsng.OCreduced148) <- as.character(unique(pop(genind.75pctmsng.OCreduced148)))
#diag(distancemat.Da.75pctmsng.OCreduced148) <- NA

heatmap.2(distancemat.Da.75pctmsng.OCreduced148, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")
saveRDS(distancemat.Da.75pctmsng.OCreduced148, file="G:/My Drive/WORKING FOLDER Jan 2015/STUFF THAT ISNT NEWCAL/Spea/RAD_data/allsocal_303indv/distancemat.Da.75pctmsng.OCreduced148.rds")
distancemat.Da.75pctmsng.OCreduced148


distance.Da.75pctmsng.OCreduced148.k3 <- genet.dist(hierfstat.75pctmsng.OCreduced148.k3, method="Da") # doesn't retain pop names, ugh!
distancemat.Da.75pctmsng.OCreduced148.k3 <- as.matrix(distance.Da.75pctmsng.OCreduced148.k3)
rownames(distancemat.Da.75pctmsng.OCreduced148.k3) <- as.character(unique(pop(genind.75pctmsng.OCreduced148.k3))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.75pctmsng.OCreduced148.k3) <- as.character(unique(pop(genind.75pctmsng.OCreduced148.k3)))
#diag(distancemat.Da.75pctmsng.OCreduced148.k3) <- NA

heatmap.2(distancemat.Da.75pctmsng.OCreduced148.k3, 
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


genpop.75pctmsng.OCreduced148 <- genind2genpop(genind.75pctmsng.OCreduced148)
distance.reynolds.75pctmsng.OCreduced148 <- dist.genpop(genpop.75pctmsng.OCreduced148, method=3) # use poppr::reynolds.dist for interindividual distances
distancemat.reynolds.75pctmsng.OCreduced148 <- as.matrix(distance.reynolds.75pctmsng.OCreduced148)
diag(distancemat.reynolds.75pctmsng.OCreduced148) <- 0 # can fill in with the intrapop distances at some point
# oddly, changing the diag from 0 to NA changes the topology of the hclust tree (for the worse, at least for Reynold's dist)
# these heatmaps are just for visualization though, so they don't really matter

heatmap.2(distancemat.reynolds.75pctmsng.OCreduced148, 
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
distance.wcfst.75pctmsng.OCreduced148.k3 <- pairwise.WCfst(hierfstat.75pctmsng.OCreduced148.k3) # standard output is a matrix
distancemat.wcfst.75pctmsng.OCreduced148.k3 <- distance.wcfst.75pctmsng.OCreduced148.k3
distance.wcfst.75pctmsng.OCreduced148.k3 <- as.dist(distancemat.wcfst.75pctmsng.OCreduced148.k3)
diag(distancemat.wcfst.75pctmsng.OCreduced148.k3) <- 0

heatmap.2(distancemat.wcfst.75pctmsng.OCreduced148.k3, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=agnes,
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Weir-Cockerman Fst (1984, via hierfstat), hclust=ward.D2, up to 75% missing data")

distance.wcfst.75pctmsng.OCreduced148 <- pairwise.WCfst(hierfstat.75pctmsng.OCreduced148) # standard output is a matrix
distancemat.wcfst.75pctmsng.OCreduced148 <- as.matrix(distance.wcfst.75pctmsng.OCreduced148)
distance.wcfst.75pctmsng.OCreduced148 <- as.dist(distancemat.wcfst.75pctmsng.OCreduced148)
diag(distancemat.wcfst.75pctmsng.OCreduced148) <- 0

heatmap.2(distancemat.wcfst.75pctmsng.OCreduced148, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=agnes,
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Weir-Cockerman Fst (1984, via hierfstat), hclust=ward.D2, up to 75% missing data")

heatmap.2(distancemat.wcfst.75pctmsng.OCall208, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=agnes,
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Weir-Cockerman Fst (1984, via hierfstat), hclust=ward.D2, up to 75% missing data")

write.table(distancemat.wcfst.75pctmsng.OCall208, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/GeneticDistance_tables/distancemat.wcfst.75pctmsng.OCall208.txt")
write.table(distancemat.wcfst.75pctmsng.OCreduced148, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/GeneticDistance_tables/distancemat.wcfst.75pctmsng.OCreduced148.txt")
stamppPhylip(as.matrix(distance.jostD.75pctmsng.OCall208), file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/GeneticDistance_tables/distancemat.jostD.75pctmsng.OCall208.phy.dst")

distancemat.wcfst.75pctmsng.OCreduced148


distancemat.wcfst.boot.75pctmsng.OCall208 <- boot.ppfst(genind.75pctmsng.OCall208)


distancemat.wcfst.75pctmsng.OCreduced148.COAST <- distancemat.wcfst.75pctmsng.OCreduced148[c("BBEND", "LAGUN", "MORO", "CCSP1", "CCSP2", "CCSP3"),c("BBEND", "LAGUN", "MORO", "CCSP1", "CCSP2", "CCSP3")]
range(lower(distancemat.wcfst.75pctmsng.OCreduced148.COAST))
mean(lower(distancemat.wcfst.75pctmsng.OCreduced148.COAST))
median(lower(distancemat.wcfst.75pctmsng.OCreduced148.COAST))

#distancemat.wcfst.75pctmsng.OCreduced148.INLAND <- distancemat.wcfst.75pctmsng.OCreduced148[c("FREM","IMSHOE",  "THOM", "LOMA", "LIME", "STARR", "GREAT"),c("FREM", "IMSHOE",  "THOM", "LOMA", "LIME", "STARR", "GREAT")]
distancemat.wcfst.75pctmsng.OCreduced148.INLAND
range(lower(distancemat.wcfst.75pctmsng.OCreduced148.INLAND))
mean(lower(distancemat.wcfst.75pctmsng.OCreduced148.INLAND))
median(lower(distancemat.wcfst.75pctmsng.OCreduced148.INLAND))

distancemat.wcfst.75pctmsng.OCreduced148.INLANDbigpops <- distancemat.wcfst.75pctmsng.OCreduced148[c("FREM","IMSHOE",  "THOM", "LOMA", "LIME", "STARR", "GREAT"),c("FREM", "IMSHOE",  "THOM", "LOMA", "LIME", "STARR", "GREAT")]
distancemat.wcfst.75pctmsng.OCreduced148.INLANDbigpops
range(lower(distancemat.wcfst.75pctmsng.OCreduced148.INLANDbigpops))
mean(lower(distancemat.wcfst.75pctmsng.OCreduced148.INLANDbigpops))
median(lower(distancemat.wcfst.75pctmsng.OCreduced148.INLANDbigpops))


distancemat.wcfst.75pctmsng.OCreduced148.SANTENA <- distancemat.wcfst.75pctmsng.OCreduced148[c("SANCAN", "TENA"),c("SANCAN", "TENA")]
distancemat.wcfst.75pctmsng.OCreduced148.SANTENA
range(lower(distancemat.wcfst.75pctmsng.OCreduced148.SANTENA))
mean(lower(distancemat.wcfst.75pctmsng.OCreduced148.SANTENA))
median(lower(distancemat.wcfst.75pctmsng.OCreduced148.SANTENA))

distancemat.Da.75pctmsng.OCreduced148.k3
range(lower(distancemat.Da.75pctmsng.OCreduced148.k3))
mean(lower(distancemat.Da.75pctmsng.OCreduced148.k3))
median(lower(distancemat.Da.75pctmsng.OCreduced148.k3))

distancemat.Da.75pctmsng.OCreduced148.COAST <- distancemat.Da.75pctmsng.OCreduced148[c("BBEND", "LAGUN", "MORO", "CCSP1", "CCSP2", "CCSP3"),c("BBEND", "LAGUN", "MORO", "CCSP1", "CCSP2", "CCSP3")]
distancemat.Da.75pctmsng.OCreduced148.COAST
range(lower(distancemat.Da.75pctmsng.OCreduced148.COAST))
mean(lower(distancemat.Da.75pctmsng.OCreduced148.COAST))
median(lower(distancemat.Da.75pctmsng.OCreduced148.COAST))

distancemat.Da.75pctmsng.OCreduced148.INLAND <- distancemat.Da.75pctmsng.OCreduced148[c("FREM", "SADDLE", "IMSHOE", "THOM", "LOMA", "LIME", "STARR", "GREAT", "TORO"),c("FREM", "SADDLE", "IMSHOE",  "THOM", "LOMA", "LIME", "STARR", "GREAT", "TORO")]
distancemat.Da.75pctmsng.OCreduced148.INLAND
range(lower(distancemat.Da.75pctmsng.OCreduced148.INLAND))
mean(lower(distancemat.Da.75pctmsng.OCreduced148.INLAND))
median(lower(distancemat.Da.75pctmsng.OCreduced148.INLAND))

distancemat.Da.75pctmsng.OCreduced148.INLANDbigpops <- distancemat.Da.75pctmsng.OCreduced148[c("FREM","IMSHOE", "THOM", "LOMA", "LIME", "STARR", "GREAT"),c("FREM", "IMSHOE",  "THOM", "LOMA", "LIME", "STARR", "GREAT")]
distancemat.Da.75pctmsng.OCreduced148.INLANDbigpops
range(lower(distancemat.Da.75pctmsng.OCreduced148.INLANDbigpops))
mean(lower(distancemat.Da.75pctmsng.OCreduced148.INLANDbigpops))
median(lower(distancemat.Da.75pctmsng.OCreduced148.INLANDbigpops))


distancemat.Da.75pctmsng.OCreduced148.SANTENA <- distancemat.Da.75pctmsng.OCreduced148[c("SANCAN", "TENA"),c("SANCAN", "TENA")]
distancemat.Da.75pctmsng.OCreduced148.SANTENA
range(lower(distancemat.Da.75pctmsng.OCreduced148.SANTENA))
mean(lower(distancemat.Da.75pctmsng.OCreduced148.SANTENA))
median(lower(distancemat.Da.75pctmsng.OCreduced148.SANTENA))



## Jost's D (2008)
library(mmod)
distance.jostD.75pctmsng.OCreduced148.k3 <- pairwise_D(genind.75pctmsng.OCreduced148.k3, linearized=FALSE)
heatmap.2(as.matrix(distance.jostD.75pctmsng.OCreduced148.k3), 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=agnes,
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Jost's D (2008, via mmod), hclust=ward.D2, up to 75% missing data")
stamppPhylip(as.matrix(distance.jostD.75pctmsng.OCreduced148.k3), 
             file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/OCreduced148_JostD_75pctmsng.k3_dist.phy.dst")

distance.jostD.75pctmsng.OCreduced148 <- pairwise_D(genind.75pctmsng.OCreduced148, linearized=FALSE)
heatmap.2(as.matrix(distance.jostD.75pctmsng.OCreduced148), 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana,
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Jost's D (2008, via mmod), hclust=ward.D2, up to 75% missing data")

distance.jostD.75pctmsng.OCall208 <- pairwise_D(genind.75pctmsng.OCall208, linearized=FALSE)
heatmap.2(as.matrix(distance.jostD.75pctmsng.OCall208), 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana,
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Jost's D (2008, via mmod), hclust=ward.D2, up to 75% missing data")

stamppPhylip(as.matrix(distance.jostD.75pctmsng.OCreduced148), 
             file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/OCreduced148_JostD_75pctmsng_dist.phy.dst")

OCreduced148.amova <- stamppAmova(distance.Da.75pctmsng.OCreduced148.indiv, gl.75pctmsng.OCreduced148)
propshared.75pctmsng.OCall208.bigpops
propshared.75pctmsng.OCall208.bigpops2 <- as.matrix(propshared.75pctmsng.OCall208.bigpops)
diag(propshared.75pctmsng.OCall208.bigpops2) <- NA

heatmap.2(propshared.75pctmsng.OCall208.bigpops2, 
          col=inferno(256, direction=-1),
          trace="none", symm=T,
          #hclustfun=diana,
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="proportion of shared loci, up to 75% missing data")
stamppPhylip(as.matrix(propshared.75pctmsng.OCall208.bigpops), 
             file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/OCall208_propsharedloci_75pctmsng_dist.phy.dst")


### diveRsity - might be faster (this doesn't f'in work.... UGH)
# 
# ip <- rgp(infile)
# diveRsity.diffcalc.75pctmsng.OCall208 <- diffCalc(infile="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/OCall208.75pctmsng_genepop.gen",
#                                                    outfile="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/diveRsity_results/",
#                                                    pairwise=TRUE,
#                                                    fst=TRUE,
#                                                    bs_pairwise = TRUE,
#                                                    boots=100,
#                                                    para=TRUE)
# 


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


pca.75pctmsng.OCreduced148 <- glPca(gl.75pctmsng.OCreduced148)
col.75pctmsng.OCreduced148 <- viridis(length(unique(pop(gl.75pctmsng.OCreduced148))))
col.75pctmsng.OCreduced148.k3 <- viridis(length(unique(pop(gl.75pctmsng.OCreduced148.k3))))
col.75pctmsng.OCreduced148.k3.sitecols <- c("#440154FF","#440154FF","#440154FF","#440154FF","#297B8EFF","#297B8EFF","#1F978BFF","#440154FF","#297B8EFF","#297B8EFF","#297B8EFF","#440154FF","#297B8EFF","#89D548FF","#297B8EFF","#89D548FF","#297B8EFF","#297B8EFF")

# [1] "#440154FF" "#21908CFF" "#FDE725FF"
# > col.75pctmsng.OCreduced148
# [1] "#440154FF" "#481769FF" "#472A7AFF" "#433D84FF" "#3D4E8AFF" "#355E8DFF" "#2E6D8EFF" "#297B8EFF" "#23898EFF" "#1F978BFF"
# [11] "#21A585FF" "#2EB37CFF" "#46C06FFF" "#65CB5EFF" "#89D548FF" "#B0DD2FFF" "#D8E219FF" "#FDE725FF"
# > unique(pop(gl.75pctmsng.OCreduced148)))
# Error: unexpected ')' in "unique(pop(gl.75pctmsng.OCreduced148)))"
# > unique(pop(gl.75pctmsng.OCreduced148))
# [1] "#440154FF"   "#440154FF" "#440154FF" "#21908CFF"  "#21908CFF"  "#21908CFF"   "#440154FF"   "#440154FF"   "#440154FF"   "#440154FF"   "#440154FF"  "#21908CFF"  "#21908CFF"  "#21908CFF"  "#440154FF"  "#440154FF"   "#FDE725FF"   "#FDE725FF"
# Levels: "#440154FF" "#440154FF" "#440154FF" "#440154FF" "#21908CFF" "#21908CFF" "#21A585FF" "#440154FF" "#21908CFF" "#21908CFF" "#21908CFF" "#440154FF" "#21908CFF" "#FDE725FF" "#21908CFF" "#FDE725FF" "#21908CFF" "#21908CFF"
# > unique(pop(gl.75pctmsng.OCreduced148.k3))
# [1] INLAND=  COAST=   SANTENA=
# Levels: COAST="#440154FF" INLAND="#297B8EFF" SANTENA="#89D548FF"
# use a slightly different shade for IMSHOE within Inland: "#21A585FF"


scatter(pca.75pctmsng.OCreduced148)
loadingplot(pca.75pctmsng.OCreduced148)
s.class(pca.75pctmsng.OCreduced148$scores, pop(gl.75pctmsng.OCreduced148),xax=1,yax=2, col=transp(col.75pctmsng.OCreduced148,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=0.5, grid=FALSE)#,

### for vector graphics, use svg or pdf; can convert to eps or whatever in inkscape

tiff("G:/My Drive/Manuscripts/Spea_socal_OC_thesischapter/figures_20190211/pca.75pctmsng.OCreduced148_pc1-8_20190211.tif",
    width=10,
    height=10,
    pointsize=12,
    units="in",
    res=600,
    antialias="default",
    type="cairo")
cairo_pdf("G:/My Drive/Manuscripts/Spea_socal_OC_thesischapter/figures_20190211/pca.75pctmsng.OCreduced148_pc1-8_20190211.pdf",
    width=10,
    height=10,
    pointsize=12,
    antialias="default",
    fallback_resolution = 600,
    onefile=TRUE)

options(bitmapType="cairo")

# cairo_ps("G:/My Drive/Manuscripts/Spea_socal_OC_thesischapter/figures_20190211/pca.75pctmsng.OCreduced148_pc1-8_20190211.eps",
#            width=10,
#            height=10,
#            pointsize=12,
#            antialias = "default",
#            fallback_resolution = 600,
#          onefile=TRUE)

svg("G:/My Drive/Manuscripts/Spea_socal_OC_thesischapter/figures_20190211/pca.75pctmsng.OCreduced148_pc1-8_20190211.svg",
         width=10,
         height=10,
         pointsize=12,
         antialias = "default",
    onefile=FALSE)

par(mfrow=c(2,2))
parmardefault <- par("mar")
par(mar=c(1,1,1,1))
par("mar")
fix(s.chull) # comment out the par args
#s.chull3(pca.75pctmsng.OCreduced148$scores, pop(gl.75pctmsng.OCreduced148),xax=1,yax=2, col=transp(col.75pctmsng.OCreduced148.k3.sitecols,.6),
#        cpoint=3, clabel=0.5, grid=FALSE, optchull=1, xlab="PC1", ylab="PC2", main="PC1, PC2")
par(mfrow=c(2,2))
s.chull(pca.75pctmsng.OCreduced148$scores, pop(gl.75pctmsng.OCreduced148),xax=1,yax=2, col=transp(col.75pctmsng.OCreduced148.k3.sitecols,.6),
        cpoint=3, clabel=0.5, grid=FALSE, optchull=1)
title("PCA, OCreduced148, up to 75% missing data per locus", line=1, xlab="PC1 (6.42% of variance)", ylab="PC2 (4.36% of variance")
s.chull(pca.75pctmsng.OCreduced148$scores, pop(gl.75pctmsng.OCreduced148),xax=3,yax=4, col=transp(col.75pctmsng.OCreduced148.k3.sitecols,.6),
        cpoint=3, clabel=0.5, grid=FALSE, optchull=1)
title("PCA, OCreduced148, up to 75% missing data per locus", line=1, xlab="PC3 (3.65% of variance)", ylab="PC4 (3.36% of variance")
s.chull(pca.75pctmsng.OCreduced148$scores, pop(gl.75pctmsng.OCreduced148),xax=5,yax=6, col=transp(col.75pctmsng.OCreduced148.k3.sitecols,.6),
        cpoint=3, clabel=0.5, grid=FALSE, optchull=1)
title("PCA, OCreduced148, up to 75% missing data per locus", line=1, xlab="PC5 (2.88% of variance)", ylab="PC6 (2.51% of variance")
s.chull(pca.75pctmsng.OCreduced148$scores, pop(gl.75pctmsng.OCreduced148),xax=7,yax=8, col=transp(col.75pctmsng.OCreduced148.k3.sitecols,.6),
        cpoint=3, clabel=0.5, grid=FALSE, optchull=1)
title("PCA, OCreduced148, up to 75% missing data per locus", line=1, xlab="PC7 (2.19% of variance)", ylab="PC8 (1.82% of variance")

cairo_pdf("G:/My Drive/Manuscripts/Spea_socal_OC_thesischapter/figures_20190211/pca.75pctmsng.OCreduced148_pc1-8_20190308.pdf",
    width=10,
    height=10,
    pointsize=12,
    antialias = "default",
    onefile=TRUE)
par(mar=c(1,1,1,1))
par(mfrow=c(2,2))
s.chull(pca.75pctmsng.OCreduced148$scores, pop(gl.75pctmsng.OCreduced148),xax=1,yax=2, col=transp(col.75pctmsng.OCreduced148.k3.sitecols,.6),
        cpoint=3, clabel=0.5, grid=FALSE, optchull=1)
title(line=0, xlab="PC1 (6.42% of variance)", ylab="PC2 (4.36% of variance")
s.chull(pca.75pctmsng.OCreduced148$scores, pop(gl.75pctmsng.OCreduced148),xax=3,yax=4, col=transp(col.75pctmsng.OCreduced148.k3.sitecols,.6),
        cpoint=3, clabel=0.5, grid=FALSE, optchull=1)
title(line=0, xlab="PC3 (3.65% of variance)", ylab="PC4 (3.36% of variance")
s.chull(pca.75pctmsng.OCreduced148$scores, pop(gl.75pctmsng.OCreduced148),xax=5,yax=6, col=transp(col.75pctmsng.OCreduced148.k3.sitecols,.6),
        cpoint=3, clabel=0.5, grid=FALSE, optchull=1)
title(line=0, xlab="PC5 (2.88% of variance)", ylab="PC6 (2.51% of variance")
s.chull(pca.75pctmsng.OCreduced148$scores, pop(gl.75pctmsng.OCreduced148),xax=7,yax=8, col=transp(col.75pctmsng.OCreduced148.k3.sitecols,.6),
        cpoint=3, clabel=0.5, grid=FALSE, optchull=1)
title(line=0, xlab="PC7 (2.19% of variance)", ylab="PC8 (1.82% of variance")


dev.off()

# % contribution - need to retain all PCs?
plot(pca.75pctmsng.OCreduced148$eig/sum(pca.75pctmsng.OCreduced148$eig)*100)
sum(pca.75pctmsng.OCreduced148$eig/sum(pca.75pctmsng.OCreduced148$eig)*100) # should be 100



s.chull(pca.75pctmsng.OCreduced148$scores, pop(gl.75pctmsng.OCreduced148),xax=7,yax=8, col=transp(col.75pctmsng.OCreduced148,.6),
        cpoint=3, clabel=0.5, grid=FALSE, optchull=1)

par(mfrow=c(2,2))
s.chull(pca.75pctmsng.OCreduced148$scores, pop(gl.75pctmsng.OCreduced148.k3),xax=1,yax=2, col=transp(col.75pctmsng.OCreduced148,.6),
        cpoint=3, clabel=0.2, grid=FALSE, optchull=1)
title("pca.75pctmsng.OCreduced148, PC1 PC2")
s.chull(pca.75pctmsng.OCreduced148$scores, pop(gl.75pctmsng.OCreduced148.k3),xax=3,yax=4, col=transp(col.75pctmsng.OCreduced148.k3,.6),
        cpoint=3, clabel=0.2, grid=FALSE, optchull=1)
title("pca.75pctmsng.OCreduced148, PC3 PC4")
s.chull(pca.75pctmsng.OCreduced148$scores, pop(gl.75pctmsng.OCreduced148.k3),xax=5,yax=6, col=transp(col.75pctmsng.OCreduced148.k3,.6),
        cpoint=3, clabel=0.5, grid=FALSE, optchull=1)
title("pca.75pctmsng.OCreduced148, PC5 PC6")
s.chull(pca.75pctmsng.OCreduced148$scores, pop(gl.75pctmsng.OCreduced148.k3),xax=7,yax=8, col=transp(col.75pctmsng.OCreduced148.k3,.6),
        cpoint=3, clabel=0.5, grid=FALSE, optchull=1)
title("pca.75pctmsng.OCreduced148, PC7 PC8")
# come up with a way to split out individual sites, but coloring them the same within the 3 FastStructure clusters

findclust1 <- find.clusters(gl.75pctmsng.OCreduced148.k3)
dapc1 <- dapc(gl.75pctmsng.OCreduced148, 
              as.factor(pop(gl.75pctmsng.OCreduced148.k3)),
              #pop=findclust1$grp, 
              glPca=pca.75pctmsng.OCreduced148,
              n.pca=10,
              n.da=2)
dapc1.optima <- optim.a.score(dapc1) # 8
dapc1.ascore <- a.score(dapc1)
scatter(dapc1)
compoplot(dapc1)




pop=as.factor(pop(gl.75pctmsng.OCreduced148.k3))


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

library(phangorn)
netcols <- palettesub #palettelist$shiny #brewer.pal(n=nPop(gl.75pctmsng), name="Spectral")
netcols <- colorRampPalette(netcols)(nPop(gl.75pctmsng.OCall208.k3))
netcols <- colorRampPalette(netcols)(nPop(gl.75pctmsng.OCall208.k3))
netcols <- col.75pctmsng.OCreduced148.k3

OCall208.75pctmsng.jostD.neighbornet <- neighborNet(distance.jostD.75pctmsng.OCall208) #input is a dist object, not a matrix, I think
plot(OCall208.75pctmsng.jostD.neighbornet) # makes an interactive 3D plot!
plot(OCall208.75pctmsng.jostD.neighbornet, type="2D",
     tip.color=netcols[unique(cbind(popcoords.OCall208$Site, popcoords.OCall208$K3hier))[,2]])
OCall208.75pctmsng.jostD.splitsnet <- splitsNetwork(distance.jostD.75pctmsng.OCall208) #input is a dist object, not a matrix, I think
plot(OCall208.75pctmsng.jostD.neighbornet, type="2D",
     tip.color=netcols[unique(pop(gl.75pctmsng.OCall208))])

# look up ?as.network to see how to manipulate the network plot
