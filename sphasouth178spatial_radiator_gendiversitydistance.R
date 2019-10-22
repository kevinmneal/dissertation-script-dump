

library(vcfR)
library(adegenet)
library(LEA)
library(tess3r)
library(hierfstat)
library(poppr)
library(radiator)
library(assigner) # use to run global and pairwise WCFst with confidence intervals

popcoords.sphasouth178spatial <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/sphasouth178spatial_popfile_coords.csv",
                                          header=T,
                                          stringsAsFactors = F)

sphasouth178spatial.vcf <- read.vcfR("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial.radiator.50pctmsng.oneSNPmac3.vcf")
sphasouth178spatial.allsnps.vcf <- read.vcfR("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial.radiator.50pctmsng.")

gl.sphasouth178spatial.radiator.50pctmsng <- vcfR2genlight(sphasouth178spatial.vcf)
genind.sphasouth178spatial.radiator.50pctmsng <- vcfR2genind(sphasouth178spatial.vcf)
pop(genind.sphasouth178spatial.radiator.50pctmsng) <- popcoords.sphasouth178spatial$poplevel2
ploidy(genind.sphasouth178spatial.radiator.50pctmsng) <- 2

tidy.sphasouth178spatial.radiator.50pctmsng <- tidy_genomic_data(data="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial.radiator.50pctmsng.oneSNPmac3.vcf",
                                                        strata="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial_popcoords_radiator_strata.txt",
                                                        parallel.core=1,
                                                        filter.common.markers=FALSE)


genomic_converter(data=tidy.sphasouth178spatial.radiator.50pctmsng,
                  strata="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial_popcoords_radiator_strata.txt",
                  parallel.core=1,
                  filter.common.markers=FALSE,
                  output="faststructure",
                  filename="sphasouth178spatial.radiator.50pctmsng")



#tidy.test <- read_rad(data="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial.radiator.50pctmsng.oneSNPmac3.rad")
#strata.test <- read_strata(strata="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial_popcoords_radiator_strata.txt",
#                           pop.id=TRUE)

indNames(gl.sphasouth178spatial.radiator.50pctmsng)
colnames(sphasouth178spatial.vcf@gt)[2:179] == indNames(gl.sphasouth178spatial.radiator.50pctmsng) # ALWAYS run this check
pop(gl.sphasouth178spatial.radiator.50pctmsng) <- popcoords.sphasouth178spatial$poplevel2
ploidy(gl.sphasouth178spatial.radiator.50pctmsng) <- 2


pca.sphasouth178spatial.radiator.50pctmsng <- glPca(gl.sphasouth178spatial.radiator.50pctmsng)
col.sphasouth178spatial.radiator.50pctmsng <- kovesi.rainbow(length(unique(pop(gl.sphasouth178spatial.radiator.50pctmsng))))
#col.sphasouth178spatial.radiator.50pctmsng <- viridis(length(unique(pop(gl.sphasouth178spatial.radiator.50pctmsng))))
#col.sphasouth178spatial.radiator.50pctmsng.sitecols <- c("#440154FF","#440154FF","#440154FF","#440154FF","#297B8EFF","#297B8EFF","#1F978BFF","#440154FF","#297B8EFF","#297B8EFF","#297B8EFF","#440154FF","#297B8EFF","#89D548FF","#297B8EFF","#89D548FF","#297B8EFF","#297B8EFF")
# coltests
#par(mfrow=c(1,1))
#pal.heatmap(parula(length(unique(pop(gl.sphasouth178spatial.radiator.50pctmsng)))))
#pal.scatter(parula(length(unique(pop(gl.sphasouth178spatial.radiator.50pctmsng)))))
#pal.safe(cols25(length(unique(pop(gl.sphasouth178spatial.radiator.50pctmsng)))))

scatter(pca.sphasouth178spatial.radiator.50pctmsng)
loadingplot(pca.sphasouth178spatial.radiator.50pctmsng)
s.class(pca.sphasouth178spatial.radiator.50pctmsng$scores, pop(gl.sphasouth178spatial.radiator.50pctmsng),xax=1,yax=2, col=transp(col.sphasouth178spatial.radiator.50pctmsng,.6), axesell=FALSE,
        cstar=1, cpoint=3, clabel=0.5, grid=FALSE)#,

cairo_pdf("G:/My Drive/Manuscripts/Spea_rangewide_thesischapter/figures_20190418/pca.sphasouth178spatial.radiator.50pctmsng_pc1-8_20190418.pdf",
          width=12,
          height=12,
          pointsize=12,
          antialias = "default",
          onefile=TRUE)
par(mar=c(1,1,1,1))
par(mfrow=c(2,2))
s.chull(pca.sphasouth178spatial.radiator.50pctmsng$scores, pop(gl.sphasouth178spatial.radiator.50pctmsng),xax=1,yax=2, col=transp(col.sphasouth178spatial.radiator.50pctmsng,.6),
        cpoint=2, clabel=0.5, grid=FALSE, optchull=1)
title(line=0, xlab="PC1 (10.89% of variance)", ylab="PC2 (6.30% of variance")
s.chull(pca.sphasouth178spatial.radiator.50pctmsng$scores, pop(gl.sphasouth178spatial.radiator.50pctmsng),xax=3,yax=4, col=transp(col.sphasouth178spatial.radiator.50pctmsng,.6),
        cpoint=2, clabel=0.5, grid=FALSE, optchull=1)
title(line=0, xlab="PC3 (4.83% of variance)", ylab="PC4 (3.95% of variance")
s.chull(pca.sphasouth178spatial.radiator.50pctmsng$scores, pop(gl.sphasouth178spatial.radiator.50pctmsng),xax=5,yax=6, col=transp(col.sphasouth178spatial.radiator.50pctmsng,.6),
        cpoint=2, clabel=0.5, grid=FALSE, optchull=1)
title(line=0, xlab="PC5 (2.74% of variance)", ylab="PC6 (2.59% of variance")
s.chull(pca.sphasouth178spatial.radiator.50pctmsng$scores, pop(gl.sphasouth178spatial.radiator.50pctmsng),xax=7,yax=8, col=transp(col.sphasouth178spatial.radiator.50pctmsng,.6),
        cpoint=2, clabel=0.5, grid=FALSE, optchull=1)
title(line=0, xlab="PC7 (2.39% of variance)", ylab="PC8 (2.14% of variance")
dev.off()

### DAPC - seems completely pointless! Optimal K is 3-5, doesn't add anything that sNMF, TESS3, and PCA don't already show
# keep all PCs for find.clusters(); keep optimal number of PCs for dapc()
par(mar=c(5,4,1,1))
findclust1 <- find.clusters(gl.sphasouth178spatial.radiator.50pctmsng,
                            n.pca=200)
findclust2 <- find.clusters(gl.sphasouth178spatial.radiator.50pctmsng,
                            n.pca=200,
                            method="ward")
findclust3 <- find.clusters(gl.sphasouth178spatial.radiator.50pctmsng,
                            n.pca=200,
                            method="kmeans")
dapc1 <- dapc(gl.sphasouth178spatial.radiator.50pctmsng, 
              #pop=as.factor(pop(gl.sphasouth178spatial.radiator.50pctmsng)),
              pop=findclust2$grp, 
              glPca=pca.sphasouth178spatial.radiator.50pctmsng)#,
              #n.pca=10,
              #n.da=2)
dapc2 <- dapc(gl.sphasouth178spatial.radiator.50pctmsng, 
              pop=findclust2$grp)
dapc3 <- dapc(gl.sphasouth178spatial.radiator.50pctmsng, 
              pop=findclust3$grp)

dapc2.optima <- optim.a.score(dapc2) # 11
dapc2 <- dapc(gl.sphasouth178spatial.radiator.50pctmsng, 
              pop=findclust1$grp,
              pca.select="percVar",
              perc.pca=50)
              #n.pca=dapc2.optima$best)

#scatter(dapc1)
#compoplot(dapc1)
par(mfrow=c(1,1))
scatter(dapc3)
compoplot(dapc3)
dapc3$grp


choosek.bic <- snapclust.choose.k(genind.sphasouth178spatial.radiator.50pctmsng, IC=BIC, max=20)
plot(choosek.bic)
choosek.kic <- snapclust.choose.k(genind.sphasouth178spatial.radiator.50pctmsng, IC=KIC, max=20)
plot(choosek.kic)
choosek.aic <- snapclust.choose.k(genind.sphasouth178spatial.radiator.50pctmsng, IC=AIC, max=20)
plot(choosek.aic)

svg(file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sNMF_crossentropy_sphasouth178spatial_radiator_50pctmsng.svg",
          width=10,
          height=10)
#cairo_pdf(file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sNMF_crossentropy_sphasouth178spatial_radiator_50pctmsng.pdf",
#    width=10,
#    height=10)
#CairoPNG(filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sNMF_crossentropy_sphasouth178spatial_radiator_50pctmsng.pdf",
#          width=10,
#          height=10)
par(mfrow=c(2,2))
boxplot(snmf.sphasouth178spatial.50pctmsng.ce, ylab="Cross-entropy", main="A) sphasouth178spatial.50pctmsng.snmf")
plot(choosek.bic, main="B) sphasouth178spatial.choosek.bic") # sphasouth178spatial: k=3 (lowest)
abline(v=2, col="darkred")
plot(choosek.kic, main="C) sphasouth178spatial.choosek.kic") # sphasouth178spatial: k=3 (elbow) or k=5 (lowest)
abline(v=5, col="darkred")
plot(choosek.aic, main="D) sphasouth178spatial.choosek.aic") # sphasouth178spatial: k=5 (elbow) or 8 (lowest)
abline(v=10, col="darkred")
dev.off()


#### genetic distance

# propshared - popgenreport
propshared.sphasouth178spatial.radiator.50pctmsng <- pairwise.propShared(genind.sphasouth178spatial.radiator.50pctmsng)
# requires 2 or more individuals per pop. Use clustering results to combine?

# basic stats

#### Genetic diversity
# summary_rad() in radiator calculates strata-level Ho, He, and FIS - use this!
summary.sphasouth178spatial.radiator.50pctmsng <- summary_rad(tidy.sphasouth178spatial.radiator.50pctmsng)
write.csv(summary.sphasouth178spatial.radiator.50pctmsng$summary.stats.pop, file="summary_rad.sphasouth178spatial.radiator.50pctmsng.csv")

# try in adegenet: Hs(), HS.test(), gengraph(); 
hs.adegenet.sphasouth178spatial.radiator.50pctmsng <- Hs(genind.sphasouth178spatial.radiator.50pctmsng)
genpop.sphasouth178spatial.radiator.50pctmsng <- genind2genpop(genind.sphasouth178spatial.radiator.50pctmsng)
gengraph.sphasouth178spatial.radiator.50pctmsng <- gengraph(genpop.sphasouth178spatial.radiator.50pctmsng, method=2)



basicstats.sphasouth178spatial.radiator.50pctmsng <- basic.stats(hierfstat.sphasouth178spatial.radiator.50pctmsng)
ho.pops.sphasouth178spatial.radiator.50pctmsng <- colMeans(basicstats.sphasouth178spatial.radiator.50pctmsng$Ho, na.rm=T)
#he_adegenet.pops.sphasouth178spatial.radiator.50pctmsng <- Hs(genind.sphasouth178spatial.radiator.50pctmsng) 
hs.pops.sphasouth178spatial.radiator.50pctmsng <- colMeans(basicstats.sphasouth178spatial.radiator.50pctmsng$Hs, na.rm=T)
fis.pops.sphasouth178spatial.radiator.50pctmsng <- colMeans(basicstats.sphasouth178spatial.radiator.50pctmsng$Fis, na.rm=T) # they're all negative... does this make sense?
#inbreeding.pops.sphasouth178spatial.radiator.50pctmsng <- inbreeding(genind.sphasouth178spatial.radiator.50pctmsng, pop=pop(genind.sphasouth178spatial.radiator.50pctmsng), res.type="estimate") # adegenet's way of calculating it; doesn't give pop estimates?
ar.sphasouth178spatial.radiator.50pctmsng <- allelic.richness(hierfstat.sphasouth178spatial.radiator.50pctmsng)
ar.pops.sphasouth178spatial.radiator.50pctmsng <- colMeans(ar.sphasouth178spatial.radiator.50pctmsng$Ar, na.rm=T)
arhohsfis.pops.sphasouth178spatial.radiator.50pctmsng <- cbind(ar.pops.sphasouth178spatial.radiator.50pctmsng, ho.pops.sphasouth178spatial.radiator.50pctmsng, hs.pops.sphasouth178spatial.radiator.50pctmsng, fis.pops.sphasouth178spatial.radiator.50pctmsng) # this is "Nei's gene diversity", which is arguably superior to standard Hexp - accounts for sampling error?
colnames(arhohsfis.pops.sphasouth178spatial.radiator.50pctmsng) <- c("Ar", "Ho", "He", "Fis")
popcounts <- count(popcoords.sphasouth178spatial, poplevel2)
arhohsfissamplesize.pops.sphasouth178spatial.radiator.50pctmsng <- cbind(popcounts, arhohsfis.pops.sphasouth178spatial.radiator.50pctmsng[order(row.names(arhohsfis.pops.sphasouth178spatial.radiator.50pctmsng)),])
View(arhohsfissamplesize.pops.sphasouth178spatial.radiator.50pctmsng)
write.csv(arhohsfissamplesize.pops.sphasouth178spatial.radiator.50pctmsng, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/arhohsfis_sphasouth178spatial.radiator.50pctmsng.csv")



# also try ibdg_fh() in radiator for the FH measure of IBDg (measures inbreeding) - "realized prop. of genome that is identical by descent relative to current population under random mating"
# IBDg is tailored for RADseq apparently
ibdg.sphasouth178spatial.radiator.50pctmsng2 <- ibdg_fh(tidy.sphasouth178spatial.radiator.50pctmsng)
ibdg.sphasouth178spatial.radiator.50pctmsng2$fh.box.plot
ibdg.sphasouth178spatial.radiator.50pctmsng$fh.stats


# use pi() in radiator to get nucleotide diversity
# requires the GT_BIN column be named GT_BIN!
# should probably use a dataset that includes invariant sites and all the SNPs before calculating pi...
pi.sphasouth178spatial.radiator.50pctmsng <- pi(tidy.sphasouth178spatial.radiator.50pctmsng, read.length=60)

betas.sphasouthspatial.radiator.50pctmsng <- betas_estimator(tidy.sphasouth178spatial.radiator.50pctmsng)
privatehaps.sphasouthspatial.radiator.50pctmsng <- private_haplotypes(tidy.sphasouth178spatial.radiator.50pctmsng)
privatealleles.sphasouthspatial.radiator.50pctmsng <- radiator::private_alleles(tidy.sphasouth178spatial.radiator.50pctmsng)


View(betas.sphasouthspatial.radiator.50pctmsng$betaiovl)
View(betas.sphasouthspatial.radiator.50pctmsng$Hw)
View(betas.sphasouthspatial.radiator.50pctmsng$Hb)
# need tidy_data_frame or GDS for radiator functions. Can convert vcf to tidy. use tidy_vcf()

hierfstat.sphasouth178spatial.radiator.50pctmsng <- genind2hierfstat(genind.sphasouth178spatial.radiator.50pctmsng)
basicstats.sphasouth178spatial.radiator.50pctmsng <- basic.stats(hierfstat.sphasouth178spatial.radiator.50pctmsng)
ho.pops.sphasouth178spatial.radiator.50pctmsng <- colMeans(basicstats.sphasouth178spatial.radiator.50pctmsng$Ho, na.rm=T)
#he_adegenet.pops.sphasouth178spatial.radiator.50pctmsng <- Hs(genind.sphasouth178spatial.radiator.50pctmsng) 
hs.pops.sphasouth178spatial.radiator.50pctmsng <- colMeans(basicstats.sphasouth178spatial.radiator.50pctmsng$Hs, na.rm=T)
fis.pops.sphasouth178spatial.radiator.50pctmsng <- colMeans(basicstats.sphasouth178spatial.radiator.50pctmsng$Fis, na.rm=T) # they're all negative... does this make sense?
inbreeding.pops.sphasouth178spatial.radiator.50pctmsng <- inbreeding(genind.sphasouth178spatial.radiator.50pctmsng, pop=pop(genind.sphasouth178spatial.radiator.50pctmsng), res.type="estimate") # adegenet's way of calculating it; doesn't give pop estimates?
ar.sphasouth178spatial.radiator.50pctmsng <- allelic.richness(hierfstat.sphasouth178spatial.radiator.50pctmsng)
ar.pops.sphasouth178spatial.radiator.50pctmsng <- colMeans(ar.sphasouth178spatial.radiator.50pctmsng$Ar, na.rm=T)

#hohe.pops.sphasouth178spatial.radiator.50pctmsng <- cbind(ho.pops.sphasouth178spatial.radiator.50pctmsng, he_adegenet.pops.sphasouth178spatial.radiator.50pctmsng)
arhohsfis.pops.sphasouth178spatial.radiator.50pctmsng <- cbind(ar.pops.sphasouth178spatial.radiator.50pctmsng, ho.pops.sphasouth178spatial.radiator.50pctmsng, hs.pops.sphasouth178spatial.radiator.50pctmsng, fis.pops.sphasouth178spatial.radiator.50pctmsng) # this is "Nei's gene diversity", which is arguably superior to standard Hexp - accounts for sampling error?
# in any case, both He and Hs give similar results. Ho in almost all cases is > He or Hs...

colnames(arhohsfis.pops.sphasouth178spatial.radiator.50pctmsng) <- c("Ar", "Ho", "He", "Fis")
#summarise(popcoords$PopLevel2_prox)
popcounts <- count(popcoords.sphasouth178spatial, SiteIMcondensed)

arhohsfissamplesize.pops.sphasouth178spatial.radiator.50pctmsng <- cbind(popcounts, arhohsfis.pops.sphasouth178spatial.radiator.50pctmsng[order(row.names(arhohsfis.pops.sphasouth178spatial.radiator.50pctmsng)),])
View(arhohsfissamplesize.pops.sphasouth178spatial.radiator.50pctmsng)
write.csv(arhohsfissamplesize.pops.sphasouth178spatial.radiator.50pctmsng, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/arhohsfis_sphasouth178spatial.radiator.50pctmsng.csv")


# also try ibdg_fh() in radiator for the FH measure of IBDg (measures inbreeding) - "realized prop. of genome that is identical by descent relative to current population under random mating"
# IBDg is tailored for RADseq apparently
ibdg.sphasouth178spatial.radiator.50pctmsng <- ibdg_fh(tidy.sphasouth178spatial.radiator.50pctmsng)
ibdg.sphasouth178spatial.radiator.50pctmsng$fh.box.plot
ibdg.sphasouth178spatial.radiator.50pctmsng$fh


# use pi() in radiator to get nucleotide diversity
# requires the GT_BIN column be named GT_BIN!
# should probably use a dataset that includes invariant sites and all the SNPs before calculating pi...
pi.sphasouth178spatial.radiator.50pctmsng <- pi(tidy.sphasouth178spatial.radiator.50pctmsng, read.length=60)





# global FST, ppfis, ppfst
boot.ppfis
boot.ppfst

hierfstat.sphasouth178spatial.radiator.50pctmsng <- genind2hierfstat(genind.sphasouth178spatial.radiator.50pctmsng)
wc.sphasouth178spatial.radiator.50pctmsng <- wc(hierfstat.sphasouth178spatial.radiator.50pctmsng)
# hierfstat global wcfst: 0.358; global WCFIS: 0.0617

## WCfst (use assigner package - has fast implementation)
tidy.sphasouth178spatial.radiator.50pctmsng2 <- tidy.sphasouth178spatial.radiator.50pctmsng
colnames(tidy.sphasouth178spatial.radiator.50pctmsng2)[11] <- "GT" # need to change GT_BIN column name to "GT"
wcfst.sphasouth178spatial.radiator.50pctmsng <- fst_WC84(data=tidy.sphasouth178spatial.radiator.50pctmsng,
                                                         ci=TRUE,
                                                         digits=5,
                                                         verbose=TRUE,
                                                         heatmap=TRUE,
                                                         filename="wcfst.sphasouth178spatial.radiator.50pctmsng")
# Fst (overall): 0.27494 [0.26976 - 0.28169]
# Fis (overall): -0.02881
# mean individual pi: 0.001352764
# median individual pi: 0.001300557
# overall pi: 0.000907633 (NOT the same as average pop or indiv pi)
# overall FH: 0.000009387084 (actually the average population FH)


# okay, I would just avoid using assigner... these values don't make sense and differ from hierfstat
pwwcfst.sphasouth178spatial.radiator.50pctmsng <- fst_WC84(data=tidy.sphasouth178spatial.radiator.50pctmsng,
                                                         ci=TRUE,
                                                         digits=5,
                                                         verbose=TRUE,
                                                         pairwise = TRUE,
                                                         heatmap=TRUE,
                                                         filename="pwwcfst.sphasouth178spatial.radiator.50pctmsng")



pwwcfstmatrix.sphasouth178spatial.radiator.50pctmsng <- `class<-`(pwwcfst.sphasouth178spatial.radiator.50pctmsng$pairwise.fst.full.matrix, 'numeric')
pwwcfstmatrix.sphasouth178spatial.radiator.50pctmsng[pwwcfstmatrix.sphasouth178spatial.radiator.50pctmsng==0] <- NA
heatmap.2(pwwcfstmatrix.sphasouth178spatial.radiator.50pctmsng,
         symm=T,
         trace="none",
         #hclustfun=diana,
         hclustfun=function(x)hclust(x, method="ward.D2"),
         col=inferno(255))
stamppPhylip(pwwcfstmatrix.sphasouth178spatial.radiator.50pctmsng, 
             file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial_wcfst_50pctmsng_dist.phy.dst")

smean.cl.boot(lower(pwwcfstmatrix.sphasouth178spatial.radiator.50pctmsng), na.rm=T)
smedian.hilow(lower(pwwcfstmatrix.sphasouth178spatial.radiator.50pctmsng), na.rm=T)

distancemat.Da.sphasouth178spatial.radiator.50pctmsng[distancemat.Da.sphasouth178spatial.radiator.50pctmsng==0] <- NA
smean.cl.boot(lower(distancemat.Da.sphasouth178spatial.radiator.50pctmsng), na.rm=T)
smedian.hilow(lower(distancemat.Da.sphasouth178spatial.radiator.50pctmsng), na.rm=T)


colMeans(pwwcfstmatrix.sphasouth178spatial.radiator.50pctmsng, na.rm=T)

pwwcfstmatrix.sphasouth178spatial.radiator.50pctmsng
distancemat.wcfst.sphasouth178spatial.radiator.50pctmsng2 <- pwwcfstmatrix.sphasouth178spatial.radiator.50pctmsng
diag(distancemat.wcfst.sphasouth178spatial.radiator.50pctmsng2) <- NA
meanwcfst.sphasouth178spatial.radiator.50pctmsng <- colMeans(distancemat.wcfst.sphasouth178spatial.radiator.50pctmsng2, na.rm=T)
sort(meanwcfst.sphasouth178spatial.radiator.50pctmsng)
View(meanwcfst.sphasouth178spatial.radiator.50pctmsng)

wcfst.vals <- NULL
wcfst.vals$wcfst.min <- apply(distancemat.wcfst.sphasouth178spatial.radiator.50pctmsng2, 2, function(x)(min(x, na.rm=TRUE)))   
wcfst.vals$wcfst.max <- apply(distancemat.wcfst.sphasouth178spatial.radiator.50pctmsng2, 2, function(x)(max(x, na.rm=TRUE)))
wcfst.vals$wcfst.mean <- apply(distancemat.wcfst.sphasouth178spatial.radiator.50pctmsng2, 2, function(x)(mean(x, na.rm=TRUE)))
wcfst.vals$wcfst.sd <- apply(distancemat.wcfst.sphasouth178spatial.radiator.50pctmsng2, 2, function(x)(sd(x, na.rm=TRUE)))
wcfst.vals$wcfst.median <- apply(distancemat.wcfst.sphasouth178spatial.radiator.50pctmsng2, 2, function(x)(median(x, na.rm=TRUE)))
wcfst.vals$wcfst.mad <- apply(distancemat.wcfst.sphasouth178spatial.radiator.50pctmsng2, 2, function(x)(mad(x, na.rm=TRUE)))

write.csv(wcfst.vals, file="wcfstvals.sphasouth178spatial.radiator.50pctmsng.csv")

# Nei's Da
## Nei's Da (1983) - attributes distance to both genetic drift and mutation
distance.Da.sphasouth178spatial.radiator.50pctmsng <- genet.dist(hierfstat.sphasouth178spatial.radiator.50pctmsng, method="Da") # doesn't retain pop names, ugh!
distancemat.Da.sphasouth178spatial.radiator.50pctmsng <- as.matrix(distance.Da.sphasouth178spatial.radiator.50pctmsng)
rownames(distancemat.Da.sphasouth178spatial.radiator.50pctmsng) <- as.character(unique(pop(genind.sphasouth178spatial.radiator.50pctmsng))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.sphasouth178spatial.radiator.50pctmsng) <- as.character(unique(pop(genind.sphasouth178spatial.radiator.50pctmsng)))
#diag(distancemat.Da.sphasouth178spatial.radiator.50pctmsng) <- NA

heatmap.2(distancemat.Da.sphasouth178spatial.radiator.50pctmsng, 
          col=inferno(255),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="sphasouth178spatial.radiator.50pctmsng Nei's Da (1983)")
saveRDS(distancemat.Da.sphasouth178spatial.radiator.50pctmsng, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/distancemat.Da.sphasouth178spatial.radiator.50pctmsng.rds")
distancemat.Da.sphasouth178spatial.radiator.50pctmsng # cailliez() confirms Da distance is Euclidean
stamppPhylip(distancemat.Da.sphasouth178spatial.radiator.50pctmsng, 
             file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial_Da_50pctmsng_dist.phy.dst")

distancemat.Da.sphasouth178spatial.radiator.50pctmsng2 <- distancemat.Da.sphasouth178spatial.radiator.50pctmsng
diag(distancemat.Da.sphasouth178spatial.radiator.50pctmsng2) <- NA
meanDa.sphasouth178spatial.radiator.50pctmsng <- colMeans(distancemat.Da.sphasouth178spatial.radiator.50pctmsng2, na.rm=T)
sort(meanDa.sphasouth178spatial.radiator.50pctmsng)
View(meanDa.sphasouth178spatial.radiator.50pctmsng)

Da.vals <- NULL
Da.vals$Da.min <- apply(distancemat.Da.sphasouth178spatial.radiator.50pctmsng2, 2, function(x)(min(x, na.rm=TRUE)))   
Da.vals$Da.max <- apply(distancemat.Da.sphasouth178spatial.radiator.50pctmsng2, 2, function(x)(max(x, na.rm=TRUE)))
Da.vals$Da.mean <- apply(distancemat.Da.sphasouth178spatial.radiator.50pctmsng2, 2, function(x)(mean(x, na.rm=TRUE)))
Da.vals$Da.sd <- apply(distancemat.Da.sphasouth178spatial.radiator.50pctmsng2, 2, function(x)(sd(x, na.rm=TRUE)))
Da.vals$Da.median <- apply(distancemat.Da.sphasouth178spatial.radiator.50pctmsng2, 2, function(x)(median(x, na.rm=TRUE)))
Da.vals$Da.mad <- apply(distancemat.Da.sphasouth178spatial.radiator.50pctmsng2, 2, function(x)(mad(x, na.rm=TRUE)))

write.csv(data.frame(Da.vals), file="sphasouth178spatial.meanDastats.radiator.50pctmsng.csv")

## chord distance
distance.Dch.sphasouth178spatial.radiator.50pctmsng <- genet.dist(hierfstat.sphasouth178spatial.radiator.50pctmsng, method="Dch") # doesn't retain pop names, ugh!
distancemat.Dch.sphasouth178spatial.radiator.50pctmsng <- as.matrix(distance.Dch.sphasouth178spatial.radiator.50pctmsng)
rownames(distancemat.Dch.sphasouth178spatial.radiator.50pctmsng) <- as.character(unique(pop(genind.sphasouth178spatial.radiator.50pctmsng))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Dch.sphasouth178spatial.radiator.50pctmsng) <- as.character(unique(pop(genind.sphasouth178spatial.radiator.50pctmsng)))
#diag(distancemat.Dch.sphasouth178spatial.radiator.50pctmsng) <- NA

heatmap.2(distancemat.Dch.sphasouth178spatial.radiator.50pctmsng, 
          col=inferno(255),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="sphasouth178spatial.radiator.50pctmsng chord distance")
saveRDS(distancemat.Dch.sphasouth178spatial.radiator.50pctmsng, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/distancemat.Dch.sphasouth178spatial.radiator.50pctmsng.rds")
cailliez(dist(distancemat.Da.sphasouth178spatial.radiator.50pctmsng)) # I don't believe Da 1983 is euclidean but cailliez() says it is
stamppPhylip(distancemat.Dch.sphasouth178spatial.radiator.50pctmsng, 
             file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial_chord_50pctmsng_dist.phy.dst")



## Weir and Cockerman Fst (1984)
distance.wcfst.sphasouth178spatial.radiator.50pctmsng <- pairwise.WCfst(hierfstat.sphasouth178spatial.radiator.50pctmsng) # standard output is a matrix
distancemat.wcfst.sphasouth178spatial.radiator.50pctmsng <- distance.wcfst.sphasouth178spatial.radiator.50pctmsng
distance.wcfst.sphasouth178spatial.radiator.50pctmsng <- as.dist(distancemat.wcfst.sphasouth178spatial.radiator.50pctmsng)
diag(distancemat.wcfst.sphasouth178spatial.radiator.50pctmsng) <- 0
saveRDS(distancemat.wcfst.sphasouth178spatial.radiator.50pctmsng)

heatmap.2(distancemat.wcfst.sphasouth178spatial.radiator.50pctmsng, 
          col=inferno(255),
          trace="none", symm=T,
          #hclustfun=agnes,
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="sphasouth178spatial.radiator.50pctmsng WCFst")


###### Jost's D (2008)
# need at least 2 individuals in each population for this to run; combine based on sNMF or tess3r clustering
library(mmod)
#distance.jostD.sphasouth153spatial.nohybrids.radiator.50pctmsng.old <- distance.jostD.sphasouth153spatial.nohybrids.radiator.50pctmsng
distance.jostD.sphasouth153spatial.nohybrids.radiator.50pctmsng <- pairwise_D(genind.sphasouth153spatial.nohybrids.radiator.50pctmsng, linearized=FALSE)
heatmap.2(as.matrix(distance.jostD.sphasouth153spatial.nohybrids.radiator.50pctmsng), 
          col=inferno(255),
          trace="none", symm=T,
          #hclustfun=agnes,
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="sphasouth153spatial.nohybrids.radiator.50pctmsng Jost's D")
stamppPhylip(as.matrix(distance.jostD.sphasouth153spatial.nohybrids.radiator.50pctmsng), 
             file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth153spatial.nohybrids_jostD_50pctmsng_dist.phy.dst")


### do it for K=15 clusters as pops, with individuals <50% maxQ removed
pwwcfst.sphasouth153spatial.nohybrids.radiator.50pctmsng <- fst_WC84(data=tidy.sphasouth153spatial.nohybrids.radiator.50pctmsng,
                                                                     ci=TRUE,
                                                                     digits=5,
                                                                     verbose=TRUE,
                                                                     pairwise = TRUE,
                                                                     filename="pwwcfst.sphasouth153spatial.nohybrids.radiator.50pctmsng")

pwwcfstmatrix.sphasouth153spatial.nohybrids.radiator.50pctmsng <- `class<-`(pwwcfst.sphasouth153spatial.nohybrids.radiator.50pctmsng$pairwise.fst.full.matrix, 'numeric')
heatmap.2(pwwcfstmatrix.sphasouth153spatial.nohybrids.radiator.50pctmsng,
          symm=T,
          trace="none",
          #hclustfun=diana,
          hclustfun=function(x)hclust(x, method="ward.D2"),
          col=inferno(255))
stamppPhylip(pwwcfstmatrix.sphasouth153spatial.nohybrids.radiator.50pctmsng, 
             file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth153spatial.nohybrids_wcfst_50pctmsng_dist.phy.dst")







# can export to phylip distance matrix for SplitsTree viewing, or use phangorn:
library(phangorn)
netcols <- palettesub #palettelist$shiny #brewer.pal(n=nPop(gl.50pctmsng), name="Spectral")
netcols <- colorRampPalette(netcols)(nPop(gl.sphasouth178spatial.radiator.50pctmsng))
netcols <- colorRampPalette(netcols)(nPop(gl.sphasouth178spatial.radiator.50pctmsng))
netcols <- col.sphasouth178spatial.radiator.50pctmsng

sphasouth178spatial.50pctmsng.jostD.neighbornet <- neighborNet(distance.jostD.sphasouth178spatial.radiator.50pctmsng) #input is a dist object, not a matrix, I think
plot(sphasouth178spatial.50pctmsng.jostD.neighbornet) # makes an interactive 3D plot!
plot(sphasouth178spatial.50pctmsng.jostD.neighbornet, type="2D",
     tip.color=netcols[unique(cbind(popcoords.sphasouth178spatial$Site, popcoords.sphasouth178spatial$K3hier))[,2]])
sphasouth178spatial.50pctmsng.jostD.splitsnet <- splitsNetwork(distance.jostD.sphasouth178spatial.radiator.50pctmsng) #input is a dist object, not a matrix, I think
plot(sphasouth178spatial.50pctmsng.jostD.neighbornet, type="2D",
     tip.color=netcols[unique(pop(gl.sphasouth178spatial.radiator.50pctmsng))])

# look up ?as.network to see how to manipulate the network plot




