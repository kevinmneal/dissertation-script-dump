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

maskbyimpervmaxent(rastertomask=topowetsouth,
                   impervmask=impervsouth,
                   impervthreshold = 100,
                   maxentmask = maxentsouth,
                   maxentthreshold = 0.1,
                   filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/topowetRGA_SPHAsouth_noimperviousmask_maxentmask10.asc")


impervsouth <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/south_ssoptim_pcarasters_hires/RAW_impervious_SPHAsouth.asc")
plot(impervsouth)

impervsouthmask50 <- impervsouth
impervsouthmask20 <- impervsouth
impervsouthmask50[impervsouthmask50>50] <- NA
plot(impervsouthmask50)
impervsouthmask20[impervsouthmask20>20] <- NA
plot(impervsouthmask20)


cmisouth <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/south_ssoptim_pcarasters_hires/climaticMoistureIndexSPHAsouth30s.asc")
cmisouthmasked <- mask(cmisouth, impervsouthmask)
plot(cmisouth)
plot(cmisouthmasked)

cmisouthmasked50 <- mask(cmisouth, impervsouthmask50)
cmisouthmasked20 <- mask(cmisouth, impervsouthmask20)
plot(cmisouthmasked20)
plot(cmisouthmasked50)

writeRaster(cmisouthmasked50, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/south_ssoptim_pcarasters_hires/climaticMoistureIndexSPHAsouth30s_imperviousmask50.asc")
writeRaster(cmisouthmasked20, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/south_ssoptim_pcarasters_hires/climaticMoistureIndexSPHAsouth30s_imperviousmask20.asc")


# threshold for maxent suitability: use 0.1? Or take mean "maximum test sensitivity plus specificity" threshold across all 10 reps?
meanmaxSS

maxentsouth <- raster("G:/My Drive/CHELSA_1.2/envlayers_south_croppedmasked/maxent_SPHAsouth.asc")
maxentsouth10thresh <- maxentsouth
maxentsouth10thresh[maxentsouth10thresh<0.1] <- NA
maxentsouth20thresh <- maxentsouth
maxentsouth20thresh[maxentsouth20thresh<0.2] <- NA
maxentsouth14thresh <- maxentsouth
maxentsouth14thresh[maxentsouth14thresh<0.14] <- NA
maxentsouth30thresh <- maxentsouth
maxentsouth30thresh[maxentsouth30thresh<0.3] <- NA

plot(stack(maxentsouth10thresh, maxentsouth14thresh, maxentsouth20thresh, maxentsouth30thresh), col=parula(255))

# extract values from raster at the locality points to find a minimum to use for the threshold that doesn't omit any (or omits some; make a histogram)
maxentsouthpoints <- read.table("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/sphasouth178spatial_allmaxentcoords_CAonly_circuitscape_nodes.txt",
                                header=F)

resisvals <- extract(maxentsouth, maxentsouthpoints[,2:3])
sort(resisvals) # use higher of lowest minus 5%, and 0.1
resisvals.samplesites <- extract(maxentsouth, unique(coordinates.sphasouth178spatial))
sort(resisvals.samplesites) # use lowest, minus 5%
souththres <- min(resisvals.samplesites, na.rm=T) - 0.05

maxentsouththres <- maxentsouth
maxentsouththres[maxentsouththres<souththres] <- NA
plot(maxentsouththres, col=parula(255))


cmisouthmasked50.maxent <- mask(cmisouthmasked50, maxentsouththres)
cmisouthmasked20.maxent <- mask(cmisouthmasked20, maxentsouththres)
cmisouthmasked50.maxent2 <- mask(cmisouthmasked50, maxentsouth10thresh)
cmisouthmasked20.maxent2 <- mask(cmisouthmasked20, maxentsouth10thresh)
plot(cmisouthmasked20.maxent)
plot(cmisouthmasked50.maxent)
points(maxentsouthpoints[,2:3])
plot(stack(cmisouthmasked50.maxent, cmisouthmasked50.maxent2, cmisouthmasked20.maxent, cmisouthmasked20.maxent2), col=parula(255))

writeRaster(cmisouthmasked50.maxent, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/south_ssoptim_pcarasters_hires/climaticMoistureIndexSPHAsouth30s_imperviousmask50_maxentmask29.asc")
writeRaster(cmisouthmasked20.maxent, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/south_ssoptim_pcarasters_hires/climaticMoistureIndexSPHAsouth30s_imperviousmask20_maxentmask29.asc")
writeRaster(cmisouthmasked50.maxent2, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/south_ssoptim_pcarasters_hires/climaticMoistureIndexSPHAsouth30s_imperviousmask50_maxentmask10.asc")
writeRaster(cmisouthmasked20.maxent2, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/south_ssoptim_pcarasters_hires/climaticMoistureIndexSPHAsouth30s_imperviousmask20_maxentmask10.asc")


cmisouthcircuitscape.cum <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/south_ssoptim_pcarasters_hires/circuitscape_outputs/climaticmoisturerga_maxentpoints_sphasouth178spatial_cum_curmap.asc")
cmisouthcircuitscape.max <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/south_ssoptim_pcarasters_hires/circuitscape_outputs/climaticmoisturerga_maxentpoints_sphasouth178spatial_curmap_max.asc")

# resistance and maxent values at populations
# put pops in alphabetical order
popcoords.sphasouth178spatial.unique.alphabetical <- unique(popcoords.sphasouth178spatial[,c(8,10,11)])
popcoords.sphasouth178spatial.unique.alphabetical <- popcoords.sphasouth178spatial.unique.alphabetical[order(popcoords.sphasouth178spatial.unique.alphabetical$poplevel2),]
sphasouth178spatial.popgenstats.maxent <- cbind(extract(aggregate(maxentsouth, fun=median, fact=2), popcoords.sphasouth178spatial.unique.alphabetical[,2:3]), arhohsfissamplesize.pops.sphasouth178spatial.radiator.50pctmsng)
sphasouth178spatial.popgenstats.maxent.resis <- cbind(extract(aggregate(cmisouth, fun=median, fact=2), popcoords.sphasouth178spatial.unique.alphabetical[,2:3]), sphasouth178spatial.popgenstats.maxent)
sphasouth178spatial.popgenstats.maxent.resis.circ1 <- cbind(extract(aggregate(cmisouthcircuitscape.cum, fun=median, fact=2), popcoords.sphasouth178spatial.unique.alphabetical[,2:3]), sphasouth178spatial.popgenstats.maxent.resis)
sphasouth178spatial.popgenstats.maxent.resis.circ2 <- cbind(extract(aggregate(cmisouthcircuitscape.max, fun=median, fact=2), popcoords.sphasouth178spatial.unique.alphabetical[,2:3]), sphasouth178spatial.popgenstats.maxent.resis.circ1)

# sphasouth178spatial.popgenlandscapestats <- cbind(popcoords.sphasouth178spatial.unique.alphabetical, arhohsfissamplesize.pops.sphasouth178spatial.radiator.50pctmsng)
# sphasouth178spatial.popgenlandscapestats$maxentmedian.aggfac3 <- extract(aggregate(maxentsouth, fun=median, fact=2), popcoords.sphasouth178spatial.unique.alphabetical[,2:3])
# sphasouth178spatial.popgenlandscapestats$cmiresismedian.aggfac3 <- extract(aggregate(cmisouth, fun=median, fact=2), popcoords.sphasouth178spatial.unique.alphabetical[,2:3])
# sphasouth178spatial.popgenlandscapestats$maxentresismedian.aggfac3 <- extract(aggregate(maxentsouthresis, fun=median, fact=2), popcoords.sphasouth178spatial.unique.alphabetical[,2:3])
# sphasouth178spatial.popgenlandscapestats$cmicirccum.median.aggfac3 <- extract(aggregate(cmisouthcircuitscape.cum, fun=median, fact=2), popcoords.sphasouth178spatial.unique.alphabetical[,2:3])
# sphasouth178spatial.popgenlandscapestats$cmicircmax.median.aggfac3 <- extract(aggregate(cmisouthcircuitscape.max, fun=median, fact=2), popcoords.sphasouth178spatial.unique.alphabetical[,2:3])
# sphasouth178spatial.popgenlandscapestats$maxentresiscirccum.median.aggfac3 <- extract(aggregate(maxentsouthresiscircuitscape.cum, fun=median, fact=2), popcoords.sphasouth178spatial.unique.alphabetical[,2:3])
# sphasouth178spatial.popgenlandscapestats$maxentresiscircmax.median.aggfac3 <- extract(aggregate(maxentsouthresiscircuitscape.max, fun=median, fact=2), popcoords.sphasouth178spatial.unique.alphabetical[,2:3])
# sphasouth178spatial.popgenlandscapestats$imperv.max.aggfac2 <- extract(aggregate(impervsouth, fun=max, fact=2), popcoords.sphasouth178spatial.unique.alphabetical[,2:3])
# sphasouth178spatial.popgenlandscapestats$imperv.max.aggfac3 <- extract(aggregate(impervsouth, fun=max, fact=3), popcoords.sphasouth178spatial.unique.alphabetical[,2:3])
# sphasouth178spatial.popgenlandscapestats$pi.allsites.60bp <- pi.sphasouth178spatial.radiator.50pctmsng$pi.populations$PI_NEI[order(as.character(pi.sphasouth178spatial.radiator.50pctmsng$pi.populations$POP_ID))]

maxentsouthresis <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/south_ssoptim_pcarasters_hires/maxentcloglogx100cv10SPHAsouth30s.asc")
topowetsouth <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/south_ssoptim_pcarasters_hires/topoWetSPHAsouth30s.asc")
rawelevsouth <- raster("G:/My Drive/CHELSA_1.2/envlayers_south_croppedmasked/elev_SPHAsouth.asc")

sphasouth178spatial.popgenlandscapestats <- cbind(popcoords.sphasouth178spatial.unique.alphabetical, arhohsfissamplesize.pops.sphasouth178spatial.radiator.50pctmsng)
sphasouth178spatial.popgenlandscapestats$maxentmedian.2km <- extract(x=maxentsouth, y=popcoords.sphasouth178spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
sphasouth178spatial.popgenlandscapestats$cmiresismedian.2km <- extract(cmisouth, y=popcoords.sphasouth178spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
sphasouth178spatial.popgenlandscapestats$topowetresismedian.2km <- extract(topowetsouth, y=popcoords.sphasouth178spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
sphasouth178spatial.popgenlandscapestats$maxentresismedian.2km <- extract(maxentsouthresis, y=popcoords.sphasouth178spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
sphasouth178spatial.popgenlandscapestats$cmicirccum.median.2km <- extract(cmisouthcircuitscape.cum, y=popcoords.sphasouth178spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
sphasouth178spatial.popgenlandscapestats$cmicircmax.median.2km <- extract(cmisouthcircuitscape.max, y=popcoords.sphasouth178spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
sphasouth178spatial.popgenlandscapestats$rawelev.median.2km <- extract(rawelevsouth, y=popcoords.sphasouth178spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
#sphasouth178spatial.popgenlandscapestats$maxentresiscirccum.median.2km <- extract(maxentsouthresiscircuitscape.cum, y=popcoords.sphasouth178spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
#sphasouth178spatial.popgenlandscapestats$maxentresiscircmax.median.2km <- extract(maxentsouthresiscircuitscape.max, y=popcoords.sphasouth178spatial.unique.alphabetical[,2:3], buffer=2000, fun=median)
sphasouth178spatial.popgenlandscapestats$imperv.mean.2km <- extract(impervsouth, y=popcoords.sphasouth178spatial.unique.alphabetical[,2:3], buffer=2000, fun=mean)
sphasouth178spatial.popgenlandscapestats$imperv.max.2km <- extract(impervsouth, y=popcoords.sphasouth178spatial.unique.alphabetical[,2:3], buffer=2000, fun=max)
sphasouth178spatial.popgenlandscapestats$imperv.mean.5km <- extract(impervsouth, y=popcoords.sphasouth178spatial.unique.alphabetical[,2:3], buffer=5000, fun=mean)
sphasouth178spatial.popgenlandscapestats$imperv.max.5km <- extract(impervsouth, y=popcoords.sphasouth178spatial.unique.alphabetical[,2:3], buffer=5000, fun=max)

pi.nooverall <- pi.sphasouth178spatial.radiator.50pctmsng$pi.populations
pi.nooverall <- pi.nooverall[-32,]
sphasouth178spatial.popgenlandscapestats$pi.pop <- pi.nooverall$PI_NEI[order(as.character(pi.nooverall$POP_ID))]


fh.nooverall <- ibdg.sphasouth178spatial.radiator.50pctmsng$fh.stats
fh.nooverall <- fh.nooverall[-32,]
#sphasouth178spatial.popgenlandscapestats$FH <- ibdg.sphasouth178spatial.radiator.50pctmsng$fh.stats$FH[order(as.character(ibdg.sphasouth178spatial.radiator.50pctmsng$fh.stats$POP_ID))]
sphasouth178spatial.popgenlandscapestats$FH <- fh.nooverall$FH[order(as.character(fh.nooverall$POP_ID))]

#write.csv(sphasouth178spatial.popgenlandscapestats, file="sphasouth178spatial.radiator.50pctmsng.popgenlandscapestats.csv")

pi.sphasouth178spatial.individuals <- as.data.frame(pi.sphasouth178spatial.radiator.50pctmsng$pi.individuals)
pi.mean.indiv <- pi.sphasouth178spatial.individuals %>%
  group_by(POP_ID) %>%
  summarise(meanpi=mean(PI), sdpi=sd(PI), medianpi=median(PI), madpi=mad(PI))

pi.mean.indiv <- pi.mean.indiv[order(as.character(pi.mean.indiv$POP_ID)),]
View(pi.mean.indiv)

sphasouth178spatial.popgenlandscapestats2 <- cbind(sphasouth178spatial.popgenlandscapestats, pi.mean.indiv)
write.csv(sphasouth178spatial.popgenlandscapestats2, file="sphasouth178spatial.radiator.50pctmsng.popgenlandscapestats.csv")


#write.csv(sphasouth178spatial.popgenstats.maxent.resis.circ2, file="sphasouth178spatial.popgenstats.maxent.cmiresis.circ_aggfac2.csv")


# manually created a table in excel with all popgen stats, Neb, maxent, and resistance values

sphasouth178spatial.bigpops.allpopgenstats <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial.radiator.50pctmsng.popgenstats_landscapevalues_bigpops.csv",
                                                       header=T,
                                                       stringsAsFactors = F)

corm.spearman <- cor.mtest(sphasouth178spatial.bigpops.allpopgenstats[,-c(1,12,13,14)], 
                           alternative="two.sided", # consider using "greater" to find significant *positive* associations
                           method="spearman", 
                           exact=FALSE)
corm.kendall <- cor.mtest(sphasouth178spatial.bigpops.allpopgenstats[,-c(1,12,13,14)], 
                          alternative="two.sided", # consider using "greater" to find significant *positive* associations
                          method="kendall", 
                          exact=FALSE)
sphasouth178spatial.cor.spearman <- cor(sphasouth178spatial.bigpops.allpopgenstats[,-c(1,12,13,14)], use="pairwise.complete.obs", method="spearman")
sphasouth178spatial.cor.spearman[corm.spearman$p > 0.05] <- 0
corrplot(sphasouth178spatial.cor.spearman, 
         title="Genetic diversity correlations, all sites",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
sphasouth178spatial.cor.kendall <- cor(sphasouth178spatial.bigpops.allpopgenstats[,-c(1,12,13,14)], use="pairwise.complete.obs", method="kendall")
sphasouth178spatial.cor.kendall[corm.kendall$p > 0.05] <- 0
corrplot(sphasouth178spatial.cor.kendall, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, all sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9)
#p.mat=corm.kendall$p, insig="pch")




####### testing using a BY correction on p-values to control for multiple testing
corm.spearman <- cor.mtest(sphasouth178spatial.bigpops.allpopgenstats[,-c(1,12,13,14)], 
                           alternative="two.sided", # consider using "greater" to find significant *positive* associations
                           method="spearman", 
                           exact=FALSE)
corm.spearman.by <- p.adjust(corm.spearman$p, method="BY")
# benjamini-yekutieli - use this one. Has been cited in other papers using pcadapt

corm.kendall <- cor.mtest(sphasouth178spatial.bigpops.allpopgenstats[,-c(1,12,13,14)], 
                          alternative="two.sided", # consider using "greater" to find significant *positive* associations
                          method="kendall", 
                          exact=FALSE)
sphasouth178spatial.cor.spearman <- cor(sphasouth178spatial.bigpops.allpopgenstats[,-c(1,12,13,14)], use="pairwise.complete.obs", method="spearman")
sphasouth178spatial.cor.spearman[p.adjust(corm.spearman$p, method="BY") > 0.1] <- 0
corrplot(sphasouth178spatial.cor.spearman, 
         title="Genetic diversity correlations, all sites",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt",
         number.cex=0.9)
#p.mat=corm.spearman$p, insig="pch")
sphasouth178spatial.cor.kendall <- cor(sphasouth178spatial.bigpops.allpopgenstats[,-c(1,12,13,14)], use="pairwise.complete.obs", method="kendall")
sphasouth178spatial.cor.kendall[p.adjust(corm.kendall$p, method="BY") > 0.1] <- 0
corrplot(sphasouth178spatial.cor.kendall, 
         method="number", 
         mar=c(0,0,2,0),
         #title="Genetic diversity kendall correlations, all sites",
         add=TRUE, type="upper",
         tl.pos="n",
         number.cex=0.9)
#p.mat=corm.kendall$p, insig="pch")



boxplot(sphasouth178spatial.bigpops.allpopgenstats$maxentmedian.aggfac3)


# run this to plot north against south in different subgraphs
require(ggplot2)
p <- ggplot(data = df.m, aes(x=variable, y=value)) 
p <- p + geom_boxplot(aes(fill = Label))
# if you want color for points replace group with colour=Label
p <- p + geom_point(aes(y=value, group=Label), position = position_dodge(width=0.75))
p <- p + facet_wrap( ~ variable, scales="free")
p <- p + xlab("x-axis") + ylab("y-axis") + ggtitle("Title")
p <- p + guides(fill=guide_legend(title="Legend_Title"))
p 