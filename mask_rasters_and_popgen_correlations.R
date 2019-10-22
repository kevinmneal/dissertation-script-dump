library(raster)

# mask a resistance surface using impervious surface layer and maxent suitability, then run through circuitscape, 
# to get a realistic connectivity surface for management

maskbyimpervmaxent <- function(rastertomask,
                               impervmask, maxentmask, impervthreshold=50, maxentthreshold=0.1, filename)
{
  rastertomask <- rastertomask
  impervmask <- impervmask
  maxentmask <- maxentmask
  impervmask[impervmask>impervthreshold] <- NA
  rastermask1 <- mask(rastertomask, impervmask)
  maxentmask[maxentmask<maxentthreshold] <- NA
  rastermask2 <- mask(rastermask1, maxentmask)
  plot(rastermask2, col=parula(255), main="raster masked by impervious and maxent values")
  writeRaster(rastermask2, filename = filename)
}

topowetsouth <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/south_ssoptim_pcarasters_hires/Results/topoWetSPHAsouth30s.asc")
maxentresissouth <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/south_ssoptim_pcarasters_hires/Results/maxentcloglogSPHAsouth30s.asc")

maskbyimpervmaxent(rastertomask=maxentresissouth,
                   impervmask=impervsouth,
                   impervthreshold = 50,
                   maxentmask = maxentsouth,
                   maxentthreshold = 0.1,
                   filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/maxentRGA_SPHAsouth_imperviousmask50_maxentmask10.asc")

maskbyimpervmaxent(rastertomask=cminorth,
                   impervmask=impervnorth,
                   impervthreshold = 100,
                   maxentmask = maxentnorth,
                   maxentthreshold = 0.1,
                   filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/cmiRGA_SPHAnorth_noimperviousmask_maxentmask10.asc")
maskbyimpervmaxent(rastertomask=maxentnorthresis,
                   impervmask=impervnorth,
                   impervthreshold = 20,
                   maxentmask = maxentnorth,
                   maxentthreshold = 0.1,
                   filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/maxentRGA_SPHAnorth_imperviousmask20_maxentmask10.asc")


impervnorth <- raster("G:/My Drive/CHELSA_1.2/envlayers_north_aggfac2/impervious_SPHANorth_aggfac2.asc")
plot(impervnorth)

impervnorthmask50 <- impervnorth
impervnorthmask20 <- impervnorth
impervnorthmask50[impervnorthmask50>50] <- NA
plot(impervnorthmask50)
impervnorthmask20[impervnorthmask20>20] <- NA
plot(impervnorthmask20)


cminorth <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphanorth124spatial/north_ssoptim_pcarasters_aggfac2/Results/climaticMoistureIndexSPHANorth.asc")
cminorthmasked <- mask(cminorth, impervnorthmask)
plot(cminorth)
plot(cminorthmasked)

cminorthmasked50 <- mask(cminorth, impervnorthmask50)
cminorthmasked20 <- mask(cminorth, impervnorthmask20)
plot(cminorthmasked20)
plot(cminorthmasked50)

writeRaster(cminorthmasked50, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/north_ssoptim_pcarasters_hires/climaticMoistureIndexSPHAnorth30s_imperviousmask50.asc")
writeRaster(cminorthmasked20, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/north_ssoptim_pcarasters_hires/climaticMoistureIndexSPHAnorth30s_imperviousmask20.asc")


# threshold for maxent suitability: use 0.1? Or take mean "maximum test sensitivity plus specificity" threshold across all 10 reps?
meanmaxSS

maxentnorth <- raster("G:/My Drive/CHELSA_1.2/envlayers_north_aggfac2/maxent_SPHANorth_aggfac2.asc")
maxentnorth10thresh <- maxentnorth
maxentnorth10thresh[maxentnorth10thresh<0.1] <- NA
maxentnorth20thresh <- maxentnorth
maxentnorth20thresh[maxentnorth20thresh<0.2] <- NA
maxentnorth14thresh <- maxentnorth
maxentnorth14thresh[maxentnorth14thresh<0.14] <- NA
maxentnorth30thresh <- maxentnorth
maxentnorth30thresh[maxentnorth30thresh<0.3] <- NA

plot(stack(maxentnorth10thresh, maxentnorth14thresh, maxentnorth20thresh, maxentnorth30thresh), col=parula(255))

# extract values from raster at the locality points to find a minimum to use for the threshold that doesn't omit any (or omits some; make a histogram)
maxentnorthpoints <- read.table("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/sphanorth124spatial_allmaxentcoords_CAonly_circuitscape_nodes.txt",
                                header=F)

maxentnorthresis <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphanorth124spatial/north_ssoptim_pcarasters_aggfac2/Results/maxentcloglogx100SPHANorth.asc")

resisvals <- extract(maxentnorth, maxentnorthpoints[,2:3])
sort(resisvals) # use higher of lowest minus 5%, and 0.1
resisvals.samplesites <- extract(maxentnorth, unique(coordinates.sphanorth124spatial))
sort(resisvals.samplesites) # use lowest, minus 5%
norththres <- min(resisvals.samplesites, na.rm=T) - 0.05

maxentnorththres <- maxentnorth
maxentnorththres[maxentnorththres<0.1] <- NA
plot(maxentnorththres, col=parula(255))


#cminorthmasked50.maxent <- mask(cminorthmasked50, maxentnorththres)
#cminorthmasked20.maxent <- mask(cminorthmasked20, maxentnorththres)
cminorthmasked50.maxent <- mask(cminorthmasked50, maxentnorththres)
cminorthmasked20.maxent <- mask(cminorthmasked20, maxentnorththres)

maxentresisnorthmasked50 <- mask(maxentnorthresis, impervnorthmask50)
maxentresisnorthmasked20 <- mask(maxentnorthresis, impervnorthmask20)
maxentresisnorthmasked50.maxent <- mask(maxentresisnorthmasked50, maxentnorththres)
maxentresisnorthmasked20.maxent <- mask(maxentresisnorthmasked50, maxentnorththres)

plot(cminorthmasked20.maxent)
plot(cminorthmasked50.maxent)
plot(maxentresisnorthmasked50.maxent)
plot(maxentresisnorthmasked20.maxent)

points(maxentnorthpoints[,2:3])
plot(stack(cminorthmasked50.maxent, cminorthmasked50.maxent2, cminorthmasked20.maxent, cminorthmasked20.maxent2), col=parula(255))

writeRaster(cminorthmasked50.maxent, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/cmiRGA_SPHAnorth_imperviousmask50_maxentmask10.asc")
writeRaster(cminorthmasked20.maxent, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/cmiRGA_SPHAnorth_imperviousmask20_maxentmask10.asc")
writeRaster(maxentresisnorthmasked50.maxent, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/maxentRGA_SPHAnorth_imperviousmask50_maxentmask10.asc")
writeRaster(maxentresisnorthmasked20.maxent, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/maxentRGA_SPHAnorth_imperviousmask20_maxentmask10.asc")


cminorthcircuitscape.cum <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/north_ssoptim_pcarasters_aggfac2/circuitscape_outputs/climaticmoisturerga_maxentpoints_sphanorth124spatial_cum_curmap.asc")
cminorthcircuitscape.max <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/north_ssoptim_pcarasters_aggfac2/circuitscape_outputs/climaticmoisturerga_maxentpoints_sphanorth124spatial_curmap_max.asc")

maxentnorthresiscircuitscape.cum <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/north_ssoptim_pcarasters_aggfac2/circuitscape_outputs/maxentrga_maxentpoints_sphanorth124spatial_cum_curmap.asc")
maxentnorthresiscircuitscape.max <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/north_ssoptim_pcarasters_aggfac2/circuitscape_outputs/maxentrga_maxentpoints_sphanorth124spatial_curmap_max.asc")


# resistance and maxent values at populations
# put pops in alphabetical order
popcoords.sphanorth124spatial.unique.alphabetical <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/sphanorth124spatial_uniquepopcoords_forrga.csv",
                                                              header=T,
                                                              stringsAsFactors = F)
popcoords.sphanorth124spatial.unique.alphabetical <- popcoords.sphanorth124spatial.unique.alphabetical[order(popcoords.sphanorth124spatial.unique.alphabetical$poplevel2),]
sphanorth124spatial.popgenlandscapestats <- cbind(popcoords.sphanorth124spatial.unique.alphabetical, arhohsfissamplesize.pops.sphanorth124spatial.radiator.50pctmsng)
# sphanorth124spatial.popgenlandscapestats$maxentmedian.aggfac3 <- extract(aggregate(maxentnorth, fun=median, fact=2), popcoords.sphanorth124spatial.unique.alphabetical[,2:3])
# sphanorth124spatial.popgenlandscapestats$cmiresismedian.aggfac3 <- extract(aggregate(cminorth, fun=median, fact=2), popcoords.sphanorth124spatial.unique.alphabetical[,2:3])
# sphanorth124spatial.popgenlandscapestats$maxentresismedian.aggfac3 <- extract(aggregate(maxentnorthresis, fun=median, fact=2), popcoords.sphanorth124spatial.unique.alphabetical[,2:3])
# sphanorth124spatial.popgenlandscapestats$cmicirccum.median.aggfac3 <- extract(aggregate(cminorthcircuitscape.cum, fun=median, fact=2), popcoords.sphanorth124spatial.unique.alphabetical[,2:3])
# sphanorth124spatial.popgenlandscapestats$cmicircmax.median.aggfac3 <- extract(aggregate(cminorthcircuitscape.max, fun=median, fact=2), popcoords.sphanorth124spatial.unique.alphabetical[,2:3])
# sphanorth124spatial.popgenlandscapestats$maxentresiscirccum.median.aggfac3 <- extract(aggregate(maxentnorthresiscircuitscape.cum, fun=median, fact=2), popcoords.sphanorth124spatial.unique.alphabetical[,2:3])
# sphanorth124spatial.popgenlandscapestats$maxentresiscircmax.median.aggfac3 <- extract(aggregate(maxentnorthresiscircuitscape.max, fun=median, fact=2), popcoords.sphanorth124spatial.unique.alphabetical[,2:3])
# sphanorth124spatial.popgenlandscapestats$imperv.max.aggfac2 <- extract(aggregate(impervnorth, fun=max, fact=2), popcoords.sphanorth124spatial.unique.alphabetical[,2:3])
# sphanorth124spatial.popgenlandscapestats$imperv.max.aggfac3 <- extract(aggregate(impervnorth, fun=max, fact=3), popcoords.sphanorth124spatial.unique.alphabetical[,2:3])
# sphanorth124spatial.popgenlandscapestats$imperv.max.aggfac4 <- extract(aggregate(impervnorth, fun=max, fact=4), popcoords.sphanorth124spatial.unique.alphabetical[,2:3])
# sphanorth124spatial.popgenlandscapestats$imperv.mean.aggfac4 <- extract(aggregate(impervnorth, fun=mean, fact=4), popcoords.sphanorth124spatial.unique.alphabetical[,2:3])
# sphanorth124spatial.popgenlandscapestats$pi.allsites.60bp <- pi.sphanorth124spatial.radiator.50pctmsng$pi.populations$PI_NEI[order(as.character(pi.sphanorth124spatial.radiator.50pctmsng$pi.populations$POP_ID))]
# 

rawelevnorth <- raster("G:/My Drive/CHELSA_1.2/envlayers_north_croppedmasked/elev_SPHANorth.asc")


sphanorth124spatial.popgenlandscapestats$maxentmedian.2km <- extract(x=maxentnorth, y=popcoords.sphanorth124spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
sphanorth124spatial.popgenlandscapestats$cmiresismedian.2km <- extract(cminorth, y=popcoords.sphanorth124spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
#sphanorth124spatial.popgenlandscapestats$topowetresismedian.2km <- extract(topowetnorth, y=popcoords.sphanorth124spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
sphanorth124spatial.popgenlandscapestats$maxentresismedian.2km <- extract(maxentnorthresis, y=popcoords.sphanorth124spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
sphanorth124spatial.popgenlandscapestats$cmicirccum.median.2km <- extract(cminorthcircuitscape.cum, y=popcoords.sphanorth124spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
sphanorth124spatial.popgenlandscapestats$cmicircmax.median.2km <- extract(cminorthcircuitscape.max, y=popcoords.sphanorth124spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
sphanorth124spatial.popgenlandscapestats$rawelev.median.2km <- extract(rawelevnorth, y=popcoords.sphanorth124spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
sphanorth124spatial.popgenlandscapestats$maxentresiscirccum.median.2km <- extract(maxentnorthresiscircuitscape.cum, y=popcoords.sphanorth124spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
sphanorth124spatial.popgenlandscapestats$maxentresiscircmax.median.2km <- extract(maxentnorthresiscircuitscape.max, y=popcoords.sphanorth124spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
sphanorth124spatial.popgenlandscapestats$imperv.mean.2km <- extract(impervnorth, y=popcoords.sphanorth124spatial.unique.alphabetical[,2:3], buffer=2000, fun=mean)
sphanorth124spatial.popgenlandscapestats$imperv.max.2km <- extract(impervnorth, y=popcoords.sphanorth124spatial.unique.alphabetical[,2:3], buffer=2000, fun=max)
sphanorth124spatial.popgenlandscapestats$imperv.mean.5km <- extract(impervnorth, y=popcoords.sphanorth124spatial.unique.alphabetical[,2:3], buffer=5000, fun=mean)
sphanorth124spatial.popgenlandscapestats$imperv.max.5km <- extract(impervnorth, y=popcoords.sphanorth124spatial.unique.alphabetical[,2:3], buffer=5000, fun=max)




fh.nooverall <- ibdg.sphanorth124spatial.radiator.50pctmsng$fh.stats
fh.nooverall <- fh.nooverall[-28,]
#sphanorth124spatial.popgenlandscapestats$FH <- ibdg.sphanorth124spatial.radiator.50pctmsng$fh.stats$FH[order(as.character(ibdg.sphanorth124spatial.radiator.50pctmsng$fh.stats$POP_ID))]
sphanorth124spatial.popgenlandscapestats$FH <- fh.nooverall$FH[order(as.character(fh.nooverall$POP_ID))]

write.csv(sphanorth124spatial.popgenlandscapestats, file="sphanorth124spatial.radiator.50pctmsng.popgenlandscapestats.csv")

pi.sphanorth124spatial.individuals <- as.data.frame(pi.sphanorth124spatial.radiator.50pctmsng$pi.individuals)
pi.sphanorth124spatial.individuals$POP_ID <- revalue(pi.sphanorth124spatial.individuals$POP_ID, c("TEJONSSIDE"="TEJON", "TEJONNSIDE"="TEJON"))
pi.mean.indiv <- pi.sphanorth124spatial.individuals %>%
  group_by(POP_ID) %>%
  summarise(meanpi=mean(PI), sdpi=sd(PI), medianpi=median(PI), madpi=mad(PI))

pi.mean.indiv <- pi.mean.indiv[order(as.character(pi.mean.indiv$POP_ID)),]
View(pi.mean.indiv)

sphanorth124spatial.popgenlandscapestats2 <- cbind(sphanorth124spatial.popgenlandscapestats, pi.mean.indiv)
write.csv(sphanorth124spatial.popgenlandscapestats2, file="sphanorth124spatial.radiator.50pctmsng.popgenlandscapestats3.csv")

# manually created a table in excel with all popgen stats, Neb, maxent, and resistance values



sphanorth124spatial.bigpops.allpopgenstats <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphanorth124spatial/sphanorth124spatial.radiator.50pctmsng.popgenlandscapestats_allpops.csv",
                                                       header=T,
                                                       stringsAsFactors = F,
                                                       row.names = 1)

corm.spearman <- cor.mtest(sphanorth124spatial.bigpops.allpopgenstats, 
                           alternative="two.sided", # consider using "greater" to find significant *positive* associations
                           method="spearman", 
                           exact=FALSE)
corm.kendall <- cor.mtest(sphanorth124spatial.bigpops.allpopgenstats, 
                          alternative="two.sided", # consider using "greater" to find significant *positive* associations
                          method="kendall", 
                          exact=FALSE)
sphanorth124spatial.cor.spearman <- cor(sphanorth124spatial.bigpops.allpopgenstats, use="pairwise.complete.obs", method="spearman")
sphanorth124spatial.cor.spearman[corm.spearman$p > 0.05] <- 0
corrplot(sphanorth124spatial.cor.spearman, 
         title="Genetic diversity correlations, all sites",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
sphanorth124spatial.cor.kendall <- cor(sphanorth124spatial.bigpops.allpopgenstats, use="pairwise.complete.obs", method="kendall")
sphanorth124spatial.cor.kendall[corm.kendall$p > 0.05] <- 0
corrplot(sphanorth124spatial.cor.kendall, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, all sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9)
#p.mat=corm.kendall$p, insig="pch")




####### testing using a BY correction on p-values to control for multiple testing
corm.spearman <- cor.mtest(sphanorth124spatial.bigpops.allpopgenstats, 
                           alternative="two.sided", # consider using "greater" to find significant *positive* associations
                           method="spearman", 
                           exact=FALSE)
corm.spearman.by <- p.adjust(corm.spearman$p, method="BY")
# benjamini-yekutieli - use this one. Has been cited in other papers using pcadapt

corm.kendall <- cor.mtest(sphanorth124spatial.bigpops.allpopgenstats, 
                          alternative="two.sided", # consider using "greater" to find significant *positive* associations
                          method="kendall", 
                          exact=FALSE)
sphanorth124spatial.cor.spearman <- cor(sphanorth124spatial.bigpops.allpopgenstats, use="pairwise.complete.obs", method="spearman")
sphanorth124spatial.cor.spearman[p.adjust(corm.spearman$p, method="BY") > 0.1] <- 0
corrplot(sphanorth124spatial.cor.spearman, 
         title="Genetic diversity correlations, sites n>3, sphanorth124spatial",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
sphanorth124spatial.cor.kendall <- cor(sphanorth124spatial.bigpops.allpopgenstats, use="pairwise.complete.obs", method="kendall")
sphanorth124spatial.cor.kendall[p.adjust(corm.kendall$p, method="BY") > 0.1] <- 0
corrplot(sphanorth124spatial.cor.kendall, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, all sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9)
#p.mat=corm.kendall$p, insig="pch")



boxplot(sphanorth124spatial.bigpops.allpopgenstats$maxentmedian.aggfac3)


# run this to plot north against north in different subgraphs
require(ggplot2)
p <- ggplot(data = df.m, aes(x=variable, y=value)) 
p <- p + geom_boxplot(aes(fill = Label))
# if you want color for points replace group with colour=Label
p <- p + geom_point(aes(y=value, group=Label), position = position_dodge(width=0.75))
p <- p + facet_wrap( ~ variable, scales="free")
p <- p + xlab("x-axis") + ylab("y-axis") + ggtitle("Title")
p <- p + guides(fill=guide_legend(title="Legend_Title"))
p 