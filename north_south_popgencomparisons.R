## north and south popgen comparisons

library(ResistanceGA)
library(reshape2)
library(ggplot2)
library(ggpubr)

# plot Da vs distance on same graph


sphasouth178spatial.geogdist
sphanorth124spatial.geogdist <- as.matrix(read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphanorth124spatial/north_ssoptim_pcarasters_aggfac2/Results/Distance_commuteDistance_distMat.csv",
                                         header=F,
                                         stringsAsFactors = F))
diag(sphanorth124spatial.geogdist) <- NA
sphanorth124spatial.Dadist <- readRDS("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/distancemat.Da.75pctmsng.sphanorth124spatial.rds")
diag(sphanorth124spatial.Dadist) <- NA

sphasouth178spatial.geogdist <- as.matrix(read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/south_ssoptim_pcarasters_hires/Distance_commuteDistance_distMat.csv",
                                                   header=F,
                                                   stringsAsFactors = F))
diag(sphasouth178spatial.geogdist) <- NA

sphasouth178spatial.Dadist <- readRDS("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/distancemat.Da.75pctmsng.sphasouth178spatial.rds")
diag(sphasouth178spatial.Dadist) <- NA


range(c(sphasouth178spatial.geogdist, sphanorth124spatial.geogdist), na.rm=T)
range(c(sphasouth178spatial.Dadist, sphanorth124spatial.Dadist), na.rm=T)

plot(sphasouth178spatial.geogdist, sphasouth178spatial.Dadist,
     col="red",
     pch=19,
     xlim=range(c(sphasouth178spatial.geogdist, sphanorth124spatial.geogdist), na.rm=T),
     ylim=range(c(sphasouth178spatial.Dadist, sphanorth124spatial.Dadist), na.rm=T),
     ylab="Genetic distance (Nei's Da)",
     xlab="Geographic distance")

points(sphanorth124spatial.geogdist, sphanorth124spatial.Dadist, col="blue", pch=19)

cor(lower(sphanorth124spatial.Dadist), lower(sphanorth124spatial.geogdist))
#mantel.randtest(dist(sphanorth124spatial.Dadist), dist(sphanorth124spatial.geogdist))
mantel(sphanorth124spatial.geogdist, sphanorth124spatial.Dadist, method="pearson")

cor(lower(sphasouth178spatial.Dadist), lower(sphasouth178spatial.geogdist))
#mantel.randtest(dist(sphasouth178spatial.Dadist), dist(sphasouth178spatial.geogdist))
mantel(sphasouth178spatial.geogdist, sphasouth178spatial.Dadist, method="pearson")


# run this to plot north against south in different subgraphs
require(ggplot2)

northsouth.popgenstats.allpops <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphanorth_south_popgencomparisons_sharedcolumns_allpops.csv",
                                   header=T,
                                   stringsAsFactors = F)

df.m.allpops <- melt(northsouth.popgenstats.allpops)

p <- ggplot(data = df.m.allpops, aes(x=variable, y=value)) 
p <- p + geom_boxplot(aes(fill = clade))
# if you want color for points replace group with colour=Label
p <- p + geom_point(aes(y=value, group=clade), position = position_dodge(width=0.75))
p <- p + facet_wrap( ~ variable, scales="free")
p <- p + xlab("x-axis") + ylab("y-axis") + ggtitle("North-South SPHA comparisons, all populations")
p <- p + guides(fill=guide_legend(title="Clade"))
p 

northsouth.popgenstats.bigpops <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphanorth_south_popgencomparisons_sharedcolumns_bigpops.csv",
                                   header=T,
                                   stringsAsFactors = F)

df.m.bigpops <- melt(northsouth.popgenstats.bigpops)

northsouth.popgenstats.bigpops.fewervars <- northsouth.popgenstats.bigpops[,c("n", "Ar", "Ho", "He", "Fis", "")]

q <- ggplot(data = df.m.bigpops, aes(x=variable, y=value)) 
q <- q + geom_boxplot(aes(fill = clade))
# if you want color for points replace group with colour=Label
q <- q + geom_point(aes(y=value, group=clade), position = position_dodge(width=0.75))
q <- q + facet_wrap( ~ variable, scales="free")
q <- q + xlab("x-axis") + ylab("y-axis") + ggtitle("North-South SPHA comparisons, populations with n>=3")
q <- q + guides(fill=guide_legend(title="Clade"))
q 



northsouth.popgenstats.fewervars <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphanorth_south_popgencomparisons_sharedcolumns_fewervars_allorn4sites3.csv",
                                           header=T,
                                           stringsAsFactors = F)

northsouth.popgenstats.fewervars <- northsouth.popgenstats.fewervars[, !names(northsouth.popgenstats.fewervars) %in% c("n")]
df.m.bigpops2 <- melt(northsouth.popgenstats.fewervars)
#hdata$Mono<- factor(hdata$Mon, levels = c(12, 1:8))
df.m.bigpops2$clade <- factor(df.m.bigpops2$clade, levels=c("SPHANorth.allsites", "SPHASouth.allsites", "SPHANorth.n4sites", "SPHASouth.n4sites"))
#northsouth.popgenstats.bigpops.fewervars <- northsouth.popgenstats.bigpops[, c("n", "Ar", "Ho", "He", "Fis", "")]


# test colors
plot(c(1,2,3,4), col=parula(100)[c(15,35,75,85)], cex=10, pch=19)

qq <- ggplot(data = df.m.bigpops2, aes(x=allor4, y=value, fill=clade)) 
qq <- qq + geom_boxplot(aes(x=allor4, fill=clade)) + scale_fill_manual(values=parula(100)[c(15,35,85,75)])
# if you want color for points replace group with colour=Label
qq <- qq + geom_point(aes(x=allor4, y=value, group=clade), position = position_dodge(width=0.75))
qq <- qq + facet_wrap( ~ variable, scales="free")# + theme(axis.text.x = element_blank())
qq <- qq + ggtitle("North-South SPHA comparisons")
qq <- qq + guides(fill=guide_legend(title="Clade"))
qq 


### for sites with 3 or more samples
northsouth.popgenstats.fewervars.n3.1 <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphanorth_south_popgencomparisons_sharedcolumns_fewervars_n3sites.csv",
                                             header=T,
                                             stringsAsFactors = F)

northsouth.popgenstats.fewervars.n3 <- northsouth.popgenstats.fewervars.n3.1[, !names(northsouth.popgenstats.fewervars.n3.1) %in% c("n")]
df.m.n3 <- melt(northsouth.popgenstats.fewervars.n3)
#hdata$Mono<- factor(hdata$Mon, levels = c(12, 1:8))
#df.m.n3$clade <- factor(df.m.bigpops2$clade, levels=c("SPHANorth.allsites", "SPHASouth.allsites", "SPHANorth.n4sites", "SPHASouth.n4sites"))
#northsouth.popgenstats.bigpops.fewervars <- northsouth.popgenstats.bigpops[, c("n", "Ar", "Ho", "He", "Fis", "")]

View(northsouth.popgenstats.fewervars.n3)
northsouth.popgenstats.fewervars.n3.2 <- northsouth.popgenstats.fewervars.n3.1[,c("clade", "pop", "n", "Ho", "He", "Fis", 
                                                                                  "IBDg.FH", "Neb", "Da.median", 
                                                                                  "wcfst.median", "maxentmedian.2km", "elev.median.2km",
                                                                                  "imperv.max.2km")]
View(northsouth.popgenstats.fewervars.n3.2)

df.m.n3 <- melt(northsouth.popgenstats.fewervars.n3.2)
View(df.m.n3)

# test colors
plot(c(1,2,3,4), col=parula(100)[c(15,35,75,85)], cex=10, pch=19)

uu <- ggplot(data = df.m.n3, aes(x=variable, y=value, fill=clade)) 
uu <- uu + geom_boxplot(aes(fill=clade)) + scale_fill_manual(values=parula(100)[c(15,85)])
# if you want color for points replace group with colour=Label
uu <- uu + geom_point(aes(y=value, group=clade), position = position_dodge(width=0.75))
uu <- uu + facet_wrap( ~ variable, scales="free") + 
                                          theme(axis.title.x=element_blank(),
                                            axis.text.x=element_blank(),
                                            axis.ticks.x=element_blank(),
                                            axis.title.y=element_blank())
uu <- uu + ggtitle("North-South SPHA comparisons, site-level")
uu <- uu + guides(fill=guide_legend(title="Clade"))
#uu + stat_compare_means(method="wilcox.test", label="p.format")
uuv <- uu + stat_compare_means(method="wilcox.test", label="p.format")
uuv
#uu + stat_compare_means(method="wilcox.test", aes(label=..p.adj..))

ggsave(plot=uuv, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/northsouth_popgencomparisons.northsouth_nopi.eps", device="eps",
       dpi=320, width=10, height=7, units="in")
ggsave(plot=uuv, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/northsouth_popgencomparisons.northsouth_nopi.png", device=png(),
       dpi=320, width=10, height=7, units="in")

compare_means(value ~ clade, group.by="variable", data=uu$data, p.adjust.method="BY")
# I'm honestly not sure if these need p-value adjustments... I want to say these are independent comparisons so they don't,
# while looking at the correlations among values within north and south are not independent and so do need p.adjust?

#rr + stat_compare_means(method="t.test")

# do a wilcox means test using ggpubr? or manually do MADAM tests with bootstrapping? too much work for the latter...

install.packages("ggpubr")
library(ggpubr)

northsouth.popgenstats.bigpops 
northsouth.popgenstats.fewervars.n3
northsouth.popgenstats.fewervars.n3.1
ns.n3.latlon <- left_join(northsouth.popgenstats.fewervars.n3.1, northsouth.popgenstats.bigpops, by="pop")
####### testing using a BY correction on p-values to control for multiple testing
corm.spearman <- cor.mtest(northsouth.popgenstats.fewervars.n3.1[,-c(1,2)], 
                           alternative="two.sided", # consider using "greater" to find significant *positive* associations
                           method="spearman", 
                           exact=FALSE)
corm.spearman.by <- p.adjust(corm.spearman$p, method="BY")
# benjamini-yekutieli - use this one. Has been cited in other papers using pcadapt

corm.kendall <- cor.mtest(northsouth.popgenstats.fewervars.n3.1[,-c(1,2)], 
                          alternative="two.sided", # consider using "greater" to find significant *positive* associations
                          method="kendall", 
                          exact=FALSE)
sphasouth178spatial.cor.spearman <- cor(northsouth.popgenstats.fewervars.n3.1[,-c(1,2)], use="pairwise.complete.obs", method="spearman")
sphasouth178spatial.cor.spearman[p.adjust(corm.spearman$p, method="BY") > 0.1] <- 0
corrplot(sphasouth178spatial.cor.spearman, 
         title="Genetic diversity correlations, sites with n>=3 samples, SPHA North and South",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
sphasouth178spatial.cor.kendall <- cor(northsouth.popgenstats.fewervars.n3.1[,-c(1,2)], use="pairwise.complete.obs", method="kendall")
sphasouth178spatial.cor.kendall[p.adjust(corm.kendall$p, method="BY") > 0.1] <- 0
corrplot(sphasouth178spatial.cor.kendall, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, all sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9)
#p.mat=corm.kendall$p, insig="pch")




### make final figure


CairoPNG(filename = "gendivcorrelations.northsouthspha.png",
    height=7,
    width = 10,
    res=320,
    units="in")
par(mfrow=c(1,2))
####### testing using a BY correction on p-values to control for multiple testing
corm.spearman <- cor.mtest(northsouth.popgenstats.fewervars.n3.1[northsouth.popgenstats.fewervars.n3.1$clade=="SPHANorth",-c(1,2)], 
                           alternative="two.sided", # consider using "greater" to find significant *positive* associations
                           method="spearman", 
                           exact=FALSE)
corm.spearman.by <- p.adjust(corm.spearman$p, method="BY")
# benjamini-yekutieli - use this one. Has been cited in other papers using pcadapt

corm.kendall <- cor.mtest(northsouth.popgenstats.fewervars.n3.1[northsouth.popgenstats.fewervars.n3.1$clade=="SPHANorth",-c(1,2)], 
                          alternative="two.sided", # consider using "greater" to find significant *positive* associations
                          method="kendall", 
                          exact=FALSE)
sphasouth178spatial.cor.spearman <- cor(northsouth.popgenstats.fewervars.n3.1[northsouth.popgenstats.fewervars.n3.1$clade=="SPHANorth",-c(1,2)], use="pairwise.complete.obs", method="spearman")
sphasouth178spatial.cor.spearman[p.adjust(corm.spearman$p, method="BY") > 0.1] <- 0
corrplot(sphasouth178spatial.cor.spearman, 
         title="Genetic diversity correlations, sites with n>=3 samples, SPHA North only",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
sphasouth178spatial.cor.kendall <- cor(northsouth.popgenstats.fewervars.n3.1[northsouth.popgenstats.fewervars.n3.1$clade=="SPHANorth",-c(1,2)], use="pairwise.complete.obs", method="kendall")
sphasouth178spatial.cor.kendall[p.adjust(corm.kendall$p, method="BY") > 0.1] <- 0
corrplot(sphasouth178spatial.cor.kendall, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, all sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9)
#p.mat=corm.kendall$p, insig="pch")




corm.spearman <- cor.mtest(northsouth.popgenstats.fewervars.n3.1[northsouth.popgenstats.fewervars.n3.1$clade=="SPHASouth",-c(1,2)], 
                           alternative="two.sided", # consider using "greater" to find significant *positive* associations
                           method="spearman", 
                           exact=FALSE)
corm.spearman.by <- p.adjust(corm.spearman$p, method="BY")
# benjamini-yekutieli - use this one. Has been cited in other papers using pcadapt

corm.kendall <- cor.mtest(northsouth.popgenstats.fewervars.n3.1[northsouth.popgenstats.fewervars.n3.1$clade=="SPHASouth",-c(1,2)], 
                          alternative="two.sided", # consider using "greater" to find significant *positive* associations
                          method="kendall", 
                          exact=FALSE)
sphasouth178spatial.cor.spearman <- cor(northsouth.popgenstats.fewervars.n3.1[northsouth.popgenstats.fewervars.n3.1$clade=="SPHASouth",-c(1,2)], use="pairwise.complete.obs", method="spearman")
sphasouth178spatial.cor.spearman[p.adjust(corm.spearman$p, method="BY") > 0.1] <- 0
corrplot(sphasouth178spatial.cor.spearman, 
         title="Genetic diversity correlations, sites with n>=3 samples, SPHA South only",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
sphasouth178spatial.cor.kendall <- cor(northsouth.popgenstats.fewervars.n3.1[northsouth.popgenstats.fewervars.n3.1$clade=="SPHANorth",-c(1,2)], use="pairwise.complete.obs", method="kendall")
sphasouth178spatial.cor.kendall[p.adjust(corm.kendall$p, method="BY") > 0.1] <- 0
corrplot(sphasouth178spatial.cor.kendall, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, all sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9)
#p.mat=corm.kendall$p, insig="pch")
dev.off()


#### consider removing the pairwise distances from this chart and making independent charts for those
#### that just show the overall distribution of pairwise values per species?
#### can also do the same for pairwise resistance distances before and after impervious/maxent filtering? that'd be interesting

boxplot(lower(sphasouth178spatial.Dadist))
boxplot(lower(sphasouth178spatial.geogdist))

sphasouth178spatial.cmidist <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/south_ssoptim_pcarasters_hires/Results/climaticMoistureIndexSPHAsouth30s_commuteDistance_distMat.csv",
                                        header=F,
                                        stringsAsFactors = F)

sphanorth124spatial.maxentdist <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphanorth124spatial/north_ssoptim_pcarasters_aggfac2/Results/maxentcloglogx100SPHANorth_commuteDistance_distMat.csv",
                                           header=F,
                                           stringsAsFactors = F)

range(c(sphasouth178spatial.cmidist, sphanorth124spatial.maxentdist), na.rm=T)
boxplot(lower(sphasouth178spatial.cmidist))

plot(lower(sphasouth178spatial.cmidist), lower(sphasouth178spatial.Dadist),
     col="red",
     pch=19,
     xlim=c(0,1500),#range(c(sphasouth178spatial.cmidist, sphanorth124spatial.maxentdist), na.rm=T),
     ylim=range(c(sphasouth178spatial.Dadist, sphanorth124spatial.Dadist), na.rm=T),
     ylab="Genetic distance (Nei's Da)",
     xlab="Resistance distance (best model per clade)")

points(lower(sphanorth124spatial.maxentdist), lower(sphanorth124spatial.Dadist), col="blue", pch=19)
points(lower(sphanorth124spatial.geogdist), lower(sphanorth124spatial.Dadist), col="lightblue", pch=19)
points(lower(sphasouth178spatial.geogdist), lower(sphasouth178spatial.Dadist), col="pink", pch=19)

boxplot(lower(sphanorth124spatial.maxentdist), lower(sphasouth178spatial.cmidist))

cor(lower(sphanorth124spatial.Dadist), lower(sphanorth124spatial.geogdist))
#mantel.randtest(dist(sphanorth124spatial.Dadist), dist(sphanorth124spatial.geogdist))
mantel(sphanorth124spatial.geogdist, sphanorth124spatial.Dadist, method="pearson")

cor(lower(sphasouth178spatial.Dadist), lower(sphasouth178spatial.geogdist))
#mantel.randtest(dist(sphasouth178spatial.Dadist), dist(sphasouth178spatial.geogdist))
mantel(sphasouth178spatial.geogdist, sphasouth178spatial.Dadist, method="pearson")






northsouth.popgenstats.n3pops.fewervars2 <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphanorth_south_popgencomparisons_sharedcolumns_bigpops_lessvars.csv",
                                                     header=T,
                                                     stringsAsFactors = F)

CairoPNG(filename = "gendivcorrelations.northsouthspha.n3pops3.png",
         height=15,
         width = 20,
         res=320,
         units="in")
par(mfrow=c(1,1))
####### testing using a BY correction on p-values to control for multiple testing
corm.spearman <- cor.mtest(northsouth.popgenstats.n3pops.fewervars2[northsouth.popgenstats.n3pops.fewervars2$clade=="SPHANorth",-c(1,2)], 
                           alternative="two.sided", # consider using "greater" to find significant *positive* associations
                           method="spearman", 
                           exact=FALSE)
corm.spearman.by <- p.adjust(corm.spearman$p, method="BY")
# benjamini-yekutieli - use this one. Has been cited in other papers using pcadapt

corm.kendall <- cor.mtest(northsouth.popgenstats.n3pops.fewervars2[northsouth.popgenstats.n3pops.fewervars2$clade=="SPHANorth",-c(1,2)], 
                          alternative="two.sided", # consider using "greater" to find significant *positive* associations
                          method="kendall", 
                          exact=FALSE)
sphasouth178spatial.cor.spearman <- cor(northsouth.popgenstats.n3pops.fewervars2[northsouth.popgenstats.n3pops.fewervars2$clade=="SPHANorth",-c(1,2)], use="pairwise.complete.obs", method="spearman")
sphasouth178spatial.cor.spearman[p.adjust(corm.spearman$p, method="BY") > 0.1] <- 0
corrplot(sphasouth178spatial.cor.spearman, 
         title="Genetic diversity correlations, sites with n>=3 samples, SPHA North only",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
sphasouth178spatial.cor.kendall <- cor(northsouth.popgenstats.n3pops.fewervars2[northsouth.popgenstats.n3pops.fewervars2$clade=="SPHANorth",-c(1,2)], use="pairwise.complete.obs", method="kendall")
sphasouth178spatial.cor.kendall[p.adjust(corm.kendall$p, method="BY") > 0.1] <- 0
corrplot(sphasouth178spatial.cor.kendall, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, all sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9)
#p.mat=corm.kendall$p, insig="pch")




corm.spearman <- cor.mtest(northsouth.popgenstats.n3pops.fewervars2[northsouth.popgenstats.n3pops.fewervars2$clade=="SPHASouth",-c(1,2)], 
                           alternative="two.sided", # consider using "greater" to find significant *positive* associations
                           method="spearman", 
                           exact=FALSE)
corm.spearman.by <- p.adjust(corm.spearman$p, method="BY")
# benjamini-yekutieli - use this one. Has been cited in other papers using pcadapt

corm.kendall <- cor.mtest(northsouth.popgenstats.n3pops.fewervars2[northsouth.popgenstats.n3pops.fewervars2$clade=="SPHASouth",-c(1,2)], 
                          alternative="two.sided", # consider using "greater" to find significant *positive* associations
                          method="kendall", 
                          exact=FALSE)
sphasouth178spatial.cor.spearman <- cor(northsouth.popgenstats.n3pops.fewervars2[northsouth.popgenstats.n3pops.fewervars2$clade=="SPHASouth",-c(1,2)], use="pairwise.complete.obs", method="spearman")
sphasouth178spatial.cor.spearman[p.adjust(corm.spearman$p, method="BY") > 0.1] <- 0
corrplot(sphasouth178spatial.cor.spearman, 
         title="Genetic diversity correlations, sites with n>=3 samples, SPHA South only",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
sphasouth178spatial.cor.kendall <- cor(northsouth.popgenstats.n3pops.fewervars2[northsouth.popgenstats.n3pops.fewervars2$clade=="SPHANorth",-c(1,2)], use="pairwise.complete.obs", method="kendall")
sphasouth178spatial.cor.kendall[p.adjust(corm.kendall$p, method="BY") > 0.1] <- 0
corrplot(sphasouth178spatial.cor.kendall, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, all sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9)
#p.mat=corm.kendall$p, insig="pch")
dev.off()





### just using spearman; north and south on one plot - can put correlations in the diagonal!!!
CairoPNG(filename = "gendivcorrelations.northsouthspha.n3pops.spearman.png",
         height=9,
         width = 12,
         res=320,
         units="in")
par(mfrow=c(1,1))
####### testing using a BY correction on p-values to control for multiple testing
corm.spearman <- cor.mtest(northsouth.popgenstats.n3pops.fewervars2[northsouth.popgenstats.n3pops.fewervars2$clade=="SPHANorth",-c(1,2)], 
                           alternative="two.sided", # consider using "greater" to find significant *positive* associations
                           method="spearman", 
                           exact=FALSE)
corm.spearman.by <- p.adjust(corm.spearman$p, method="BY")
# benjamini-yekutieli - use this one. Has been cited in other papers using pcadapt

corm.kendall <- cor.mtest(northsouth.popgenstats.n3pops.fewervars2[northsouth.popgenstats.n3pops.fewervars2$clade=="SPHASouth",-c(1,2)], 
                          alternative="two.sided", # consider using "greater" to find significant *positive* associations
                          method="kendall", 
                          exact=FALSE)
sphasouth178spatial.cor.spearman <- cor(northsouth.popgenstats.n3pops.fewervars2[northsouth.popgenstats.n3pops.fewervars2$clade=="SPHASouth",-c(1,2)], use="pairwise.complete.obs", method="spearman")
sphasouth178spatial.cor.spearman[p.adjust(corm.spearman$p, method="BY") > 0.1] <- 0
corrplot(sphasouth178spatial.cor.spearman, 
         title="Genetic diversity correlations, sites with n>=3 samples",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9,
         diag=TRUE)
#p.mat=corm.spearman$p, insig="pch")
sphasouth178spatial.cor.kendall <- cor(northsouth.popgenstats.n3pops.fewervars2[northsouth.popgenstats.n3pops.fewervars2$clade=="SPHASouth",-c(1,2)], use="pairwise.complete.obs", method="kendall")
sphasouth178spatial.cor.kendall[p.adjust(corm.kendall$p, method="BY") > 0.1] <- 0
corrplot(sphasouth178spatial.cor.kendall, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, all sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9,
         diag=FALSE)
#p.mat=corm.kendall$p, insig="pch")


dev.off()



