distancemat.wcfst.75pctmsng.OCall208
distancemat.wcfst.75pctmsng.OCall208.2 <- distancemat.wcfst.75pctmsng.OCall208
diag(distancemat.wcfst.75pctmsng.OCall208.2) <- NA
distancemat.wcfst.75pctmsng.bigpops

library(gplots)

library(pals)

CairoEPS(filename="G:/My Drive/Manuscripts/Spea_socal_OC_thesischapter/figures_20190211/wcfst.OCall208.75pctmsng_higherres2.tif",
         height=10, width=10, res=320, units="in", bg="white")
heatmap.2(distancemat.wcfst.75pctmsng.OCall208, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          #hclustfun=function(x)hclust(x, method="complete"),
          #hclustfun=function(x)hclust(x, method="single"),
          #hclustfun=function(x)hclust(x, method="average"),
          #hclustfun=function(x)hclust(x, method="mcquitty"),
          main="Weir-Cockerham Pairwise Fst, OC_all_208")
dev.off()







OCall2088.gendiv.old <- OCall208.gendiv
OCall208.gendiv <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/stacks_populations/OCall208.75pctmsng.fullsumstats_veryfewstatsatall_20190610.csv",
                            header=T,
                            stringsAsFactors = T)
View(OCall208.gendiv)
pairs(OCall208.gendiv)


CairoPNG(filename="G:/My Drive/Manuscripts/Spea_socal_OC_thesischapter/figures_20190211/sphaOCclust85_popgencorr_all_higherres.png",
         height=10,
         width=7,
         res=320,
         units="in",
         bg="white")
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
OCall208.cor.kendall.COAST[corm.kendall.COAST$p > 0.05] <- 0
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
OCall208.cor.kendall.INLAND[corm.kendall.INLAND$p > 0.05] <- 0
corrplot(OCall208.cor.kendall.INLAND, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, Inland sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9)


dev.off()








par(mfrow=c(3,2))

#### with adjusted p-values. not sure if this is the right thing to do though but might as well use it

corm.spearman <- cor.mtest(OCall208.gendiv[,-c(1,2,3,4)], 
                           alternative="two.sided", # consider using "greater" to find significant *positive* associations
                           method="spearman", 
                           exact=FALSE)
corm.kendall <- cor.mtest(OCall208.gendiv[,-c(1,2,3,4)], 
                          alternative="two.sided", # consider using "greater" to find significant *positive* associations
                          method="kendall", 
                          exact=FALSE)
OCall208.cor.spearman <- cor(OCall208.gendiv[,-c(1,2,3,4)], use="pairwise.complete.obs", method="spearman")
OCall208.cor.spearman[p.adjust(corm.spearman$p) > 0.1] <- 0
#OCall208.cor.spearman[corm.spearman$p > 0.05] <- 0
corrplot(OCall208.cor.spearman, 
         title="Genetic diversity correlations, all sites",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
OCall208.cor.kendall <- cor(OCall208.gendiv[,-c(1,2,3,4)], use="pairwise.complete.obs", method="kendall")
OCall208.cor.kendall[p.adjust(corm.kendall$p) > 0.1] <- 0
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
OCall208.cor.spearman.ARTIF[p.adjust(corm.spearman.ARTIF$p) > 0.1] <- 0
corrplot(OCall208.cor.spearman.ARTIF, 
         title="Genetic diversity correlations, artificial sites",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
OCall208.cor.kendall.ARTIF <- cor(OCall208.gendiv[OCall208.gendiv$Artifnatural=="ARTIF",-c(1,2,3,4)], use="pairwise.complete.obs", method="kendall")
OCall208.cor.kendall.ARTIF[p.adjust(corm.kendall.ARTIF$p) > 0.1] <- 0
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
OCall208.cor.spearman.NATURAL[p.adjust(corm.spearman.NATURAL$p) > 0.1] <- 0
corrplot(OCall208.cor.spearman.NATURAL, 
         title="Genetic diversity correlations, all natural sites",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
OCall208.cor.kendall.NATURAL <- cor(OCall208.gendiv[OCall208.gendiv$Artifnatural=="NATURAL",-c(1,2,3,4)], use="pairwise.complete.obs", method="kendall")
OCall208.cor.kendall.NATURAL[p.adjust(corm.kendall.NATURAL$p) > 0.1] <- 0
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
OCall208.cor.spearman.COAST[p.adjust(corm.spearman.COAST$p) > 0.1] <- 0
corrplot(OCall208.cor.spearman.COAST, 
         title="Genetic diversity correlations, Coast sites",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
OCall208.cor.kendall.COAST <- cor(OCall208.gendiv[OCall208.gendiv$ClusterK4=="Coast",-c(1,2,3,4)], use="pairwise.complete.obs", method="kendall")
OCall208.cor.kendall.COAST[p.adjust(corm.kendall.COAST$p) > 0.1] <- 0
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
OCall208.cor.spearman.INLAND[p.adjust(corm.spearman.INLAND$p) > 0.1] <- 0
#OCall208.cor.spearman.INLAND[corm.spearman.INLAND$p > 0.05] <- 0
corrplot(OCall208.cor.spearman.INLAND, 
         title="Genetic diversity correlations, natural Inland sites",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
OCall208.cor.kendall.INLAND <- cor(OCall208.gendiv[OCall208.gendiv$ClusterK4=="Inlandnatural",-c(1,2,3,4)], use="pairwise.complete.obs", method="kendall")
OCall208.cor.kendall.INLAND[p.adjust(corm.kendall.INLAND$p) > 0.1] <- 0
#OCall208.cor.kendall.INLAND[corm.kendall.INLAND$p > 0.05] <- 0
corrplot(OCall208.cor.kendall.INLAND, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, Inland sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9)




