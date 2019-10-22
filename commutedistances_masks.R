### making commuteDistance on masked rasters to see differences in commute distances between 
### standard resistnace layers and those masked by impervious surfaces


# test

#southresistancelayers <- stack(list.files("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/", pattern=".*SPHAsouth_.*asc$", full.names = TRUE))
southresistancelayers <- stack(list.files("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/", pattern="cmiRGA.*SPHAsouth_.*asc$", full.names = TRUE))
southresistancelayers2 <- southresistancelayers
#southresistancelayers.10000NA <- southresistancelayers2
#southresistancelayers2[is.na(southresistancelayers)] <- 10000
#southresistancelayers$cmiRGA_SPHAsouth_imperviousmask20_maxentmask10 <- mask(southresistancelayers$cmiRGA_SPHAsouth_imperviousmask20_maxentmask10, cmisouth)
southresistancelayers2 <- setMinMax(southresistancelayers)
maxValue(southresistancelayers2)

for (i in 1:length(names(southresistancelayers2))){
  southresistancelayers2[[i]][is.na(southresistancelayers2[[i]])] <- 10000 #<- maxValue(southresistancelayers2[[i]]*10)
  southresistancelayers2[[i]] <- mask(southresistancelayers2[[i]], cmisouth)
}
plot(southresistancelayers2, col=parula(100))

northresistancelayers <- stack(list.files("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/", pattern="maxentRGA.*SPHAnorth_.*asc$", full.names = TRUE))
northresistancelayers2 <- northresistancelayers
#northresistancelayers.10000NA <- northresistancelayers2
northresistancelayers2 <- setMinMax(northresistancelayers)
maxValue(northresistancelayers2)

maxentresisnorth <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/maxentcloglogx100SPHANorth.asc")
for (i in 1:length(names(northresistancelayers2))){
  northresistancelayers2[[i]][is.na(northresistancelayers2[[i]])] <- 10000 #<- maxValue(northresistancelayers2[[i]]*10)
  northresistancelayers2[[i]] <- mask(northresistancelayers2[[i]], maxentresisnorth)
}
plot(northresistancelayers2, col=parula(100))

# use longlat=TRUE, directions=8, transitionFunction=function(x) 1 / mean(x)

unique(popcoords.sphasouth178spatial[,c(8,10,11)])

sphasouth178.rgacoords <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/sphasouth178spatial_uniquepopcoords_forrga.csv",
         header=T,
         stringsAsFactors = F)
spdf.sphasouth178spatial.pops <- SpatialPoints(sphasouth178.rgacoords[,c(2,3)])

hist(southresistancelayers$cmiRGA_SPHAsouth_imperviousmask20_maxentmask10)
southresistancelayers$cmiRGA_SPHAsouth_imperviousmask20_maxentmask10[is.na(southresistancelayers$cmiRGA_SPHAsouth_imperviousmask20_maxentmask10)] <- 10000
southresistancelayers$cmiRGA_SPHAsouth_imperviousmask20_maxentmask10 <- mask(southresistancelayers$cmiRGA_SPHAsouth_imperviousmask20_maxentmask10, cmisouth)
testtrans <- transition(southresistancelayers$cmiRGA_SPHAsouth_imperviousmask20_maxentmask10,
                        transitionFunction = function(x) 1/mean(x),
                        directions=8)
testtransgeo <- geoCorrection(testtrans, type="r")
testcommute4 <- commuteDistance(testtransgeo,
                               coords=spdf.sphasouth178spatial.pops)
testcommute3
View(as.matrix(testcommute4))
hist(testcommute3, breaks=100)
hist(testcommute4, breaks=100)

sphanorth124.rgacoords <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/sphanorth124spatial_uniquepopcoords_forrga.csv",
                                   header=T,
                                   stringsAsFactors = F)
spdf.sphanorth124spatial.pops <- SpatialPoints(sphanorth124.rgacoords[,c(2,3)])


southcommutedists.nax10 <- list()
northcommutedists.nax10 <- list()

plot(northresistancelayers2, col=parula(20))
plot(southresistancelayers2, col=parula(20))

for (i in 1:length(names(southresistancelayers2))){
  transsouth <- transition(southresistancelayers2[[i]],
             transitionFunction = function(x) 1/mean(x),
             directions=8)
  transgeo <- geoCorrection(transsouth, type="r")
  commdistout <- commuteDistance(transgeo, coords=spdf.sphasouth178spatial.pops)
  southcommutedists.nax10[[i]] <- as.matrix(commdistout)
  rownames(southcommutedists[[1]]) <- sphasouth178.rgacoords$SPHASpoplvl2
  colnames(southcommutedists[[1]]) <- sphasouth178.rgacoords$SPHASpoplvl2
  names(southcommutedists)[i] <- paste0(names(southresistancelayers2[[i]]),"_","commuteDistance")
}


for (i in 1:length(names(northresistancelayers2))){
  transnorth <- transition(northresistancelayers2[[i]],
                           transitionFunction = function(x) 1/mean(x),
                           directions=8)
  transgeo <- geoCorrection(transnorth, type="r")
  commdistout <- commuteDistance(transgeo, coords=spdf.sphanorth124spatial.pops)
  northcommutedists.nax10[[i]] <- as.matrix(commdistout)
  rownames(northcommutedists[[1]]) <- sphanorth124.rgacoords$poplevel2
  colnames(northcommutedists[[1]]) <- sphanorth124.rgacoords$poplevel2
  names(northcommutedists)[i] <- paste0(names(northresistancelayers2[[i]]),"_","commuteDistance")
}

#save old settings
op <- par(no.readonly = TRUE)
#change settings
par(mar=c(8, 4, 2, 2) + 0.1)

### get named pairs
# keep only the lower triangle of the matrix and then melt it; only need to do this for one matrix and make this a column in the df
m <- southcommutedists$cmiRGA_SPHAsouth_imperviousmask20_maxentmask10_commuteDistance
m[!lower.tri(m)] <- NA
View(m)
southpairmelt <- cbind(reshape2::melt(m, varnames = c('pop1', 'pop2'), na.rm = TRUE), lower(southcommutedists$cmiRGA_SPHAsouth_imperviousmask20_maxentmask10_commuteDistance))
View(southpairmelt)
southpairmelt.concat <- paste0(southpairmelt$pop1,".",southpairmelt$pop2)
View(southpairmelt.concat)


# test that order is preserved as it is with using lower()
testpairmelt$value==testpairmelt$`lower(southcommutedists$cmiRGA_SPHAsouth_imperviousmask20_maxentmask10_commuteDistance)`


south.commutedists.df <- NULL
south.commutedists.df$none <- lower(southcommutedists$cmiRGA_SPHAsouth_noimperviousmask_maxentmask10_commuteDistance)
south.commutedists.df$imp50 <- lower(southcommutedists$cmiRGA_SPHAsouth_imperviousmask50_maxentmask10_commuteDistance)
south.commutedists.df$imp20 <- lower(southcommutedists$cmiRGA_SPHAsouth_imperviousmask20_maxentmask10_commuteDistance)
#south.commutedists.df$cmiRGA.maxent10 <- lower(southcommutedists$cmiRGA_SPHAsouth_noimperviousmask_maxentmask10_commuteDistance)
#south.commutedists.df$cmiRGA.maxent10.imp50 <- lower(southcommutedists$cmiRGA_SPHAsouth_imperviousmask50_maxentmask10_commuteDistance)
#south.commutedists.df$cmiRGA.maxent10.imp20 <- lower(southcommutedists$cmiRGA_SPHAsouth_imperviousmask20_maxentmask10_commuteDistance)
#south.commutedists.df$maxentRGA.maxent10 <- lower(southcommutedists$maxentRGA_SPHAsouth_noimperviousmask_maxentmask10_commuteDistance)
#south.commutedists.df$maxentRGA.maxent10.imp50 <- lower(southcommutedists$maxentRGA_SPHAsouth_imperviousmask50_maxentmask10_commuteDistance)
#south.commutedists.df$maxentRGA.maxent10.imp20 <- lower(southcommutedists$maxentRGA_SPHAsouth_imperviousmask20_maxentmask10_commuteDistance)
#south.commutedists.df$topowetRGA.maxent10 <- lower(southcommutedists$topowetRGA_SPHAsouth_noimperviousmask_maxentmask10_commuteDistance)
#south.commutedists.df$topowetRGA.maxent10.imp50 <- lower(southcommutedists$topowetRGA_SPHAsouth_imperviousmask50_maxentmask10_commuteDistance)
#south.commutedists.df$topowetRGA.maxent10.imp20 <- lower(southcommutedists$topowetRGA_SPHAsouth_imperviousmask20_maxentmask10_commuteDistance)
#south.commutedists.df$populationpair <- factor(c(1:length(south.commutedists.df$cmiRGA.maxent10)))
south.commutedists.df$populationpair <- factor(southpairmelt.concat)
south.commutedists.df$clade <- c(rep("SPHASouth", times=length(south.commutedists.df$imp50)))
#View(south.commutedists.df$clade)
south.commutedists.df2 <- data.frame(south.commutedists.df)
View(south.commutedists.df2)

south.commutedists.melt <- melt(south.commutedists.df2, value.name="resis.dist", variable.name="mask", id.vars=c("clade", "populationpair"))
View(south.commutedists.melt)
#south.commutedists.melt <- cbind(rep("SPHASouth", n=length(south.commutedists.melt$resis.dist)), south.commutedists.melt)
#colnames(south.commutedists.melt) <- c("clade", "resis.dist", "layer")

write.csv(south.commutedists.df2, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/south.commuteDists.cmiRGA.impervmaxentmasks.poppairs.df.csv")


### get named pairs
# keep only the lower triangle of the matrix and then melt it; only need to do this for one matrix and make this a column in the df
m <- northcommutedists$maxentRGA_SPHAnorth_imperviousmask20_maxentmask10_commuteDistance
m[!lower.tri(m)] <- NA
View(m)
northpairmelt <- cbind(reshape2::melt(m, varnames = c('pop1', 'pop2'), na.rm = TRUE), lower(northcommutedists$maxentRGA_SPHAnorth_imperviousmask20_maxentmask10_commuteDistance))
View(northpairmelt)
northpairmelt.concat <- paste0(northpairmelt$pop1,".",northpairmelt$pop2)
View(northpairmelt.concat)


# test that order is preserved as it is with using lower()
northpairmelt$value==northpairmelt$`lower(northcommutedists$maxentRGA_SPHAnorth_imperviousmask20_maxentmask10_commuteDistance)`


north.commutedists.df <- NULL
north.commutedists.df$none <- lower(northcommutedists$maxentRGA_SPHAnorth_noimperviousmask_maxentmask10_commuteDistance)
north.commutedists.df$imp50 <- lower(northcommutedists$maxentRGA_SPHAnorth_imperviousmask50_maxentmask10_commuteDistance)
north.commutedists.df$imp20 <- lower(northcommutedists$maxentRGA_SPHAnorth_imperviousmask20_maxentmask10_commuteDistance)
#north.commutedists.df$maxentRGA.maxent10 <- lower(northcommutedists$maxentRGA_SPHAnorth_noimperviousmask_maxentmask10_commuteDistance)
#north.commutedists.df$maxentRGA.maxent10.imp50 <- lower(northcommutedists$maxentRGA_SPHAnorth_imperviousmask50_maxentmask10_commuteDistance)
#north.commutedists.df$maxentRGA.maxent10.imp20 <- lower(northcommutedists$maxentRGA_SPHAnorth_imperviousmask20_maxentmask10_commuteDistance)
#north.commutedists.df$maxentRGA.maxent10 <- lower(northcommutedists$maxentRGA_SPHAnorth_noimperviousmask_maxentmask10_commuteDistance)
#north.commutedists.df$maxentRGA.maxent10.imp50 <- lower(northcommutedists$maxentRGA_SPHAnorth_imperviousmask50_maxentmask10_commuteDistance)
#north.commutedists.df$maxentRGA.maxent10.imp20 <- lower(northcommutedists$maxentRGA_SPHAnorth_imperviousmask20_maxentmask10_commuteDistance)
#north.commutedists.df$topowetRGA.maxent10 <- lower(northcommutedists$topowetRGA_SPHAnorth_noimperviousmask_maxentmask10_commuteDistance)
#north.commutedists.df$topowetRGA.maxent10.imp50 <- lower(northcommutedists$topowetRGA_SPHAnorth_imperviousmask50_maxentmask10_commuteDistance)
#north.commutedists.df$topowetRGA.maxent10.imp20 <- lower(northcommutedists$topowetRGA_SPHAnorth_imperviousmask20_maxentmask10_commuteDistance)
#north.commutedists.df$populationpair <- factor(c(1:length(north.commutedists.df$maxentRGA.maxent10)))
north.commutedists.df$populationpair <- factor(northpairmelt.concat)
north.commutedists.df$clade <- c(rep("SPHANorth", times=length(north.commutedists.df$imp50)))
#View(north.commutedists.df$clade)
north.commutedists.df2 <- data.frame(north.commutedists.df)
View(north.commutedists.df2)

north.commutedists.melt <- melt(north.commutedists.df2, value.name="resis.dist", variable.name="mask", id.vars=c("clade", "populationpair"))
View(north.commutedists.melt)
#north.commutedists.melt <- cbind(rep("SPHASouth", n=length(north.commutedists.melt$resis.dist)), north.commutedists.melt)
#colnames(north.commutedists.melt) <- c("clade", "resis.dist", "layer")

write.csv(north.commutedists.df2, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/north.commuteDists.maxentRGA.impervmaxentmasks.poppairs.df.csv")






#north.commutedists.melt <- cbind(rep("SPHANorth", n=length(north.commutedists.melt$resis.dist)), north.commutedists.melt)
#colnames(north.commutedists.melt) <- c("clade", "resis.dist", "layer")

northsouth.commutedists.melt <- dplyr::bind_rows(north.commutedists.melt, south.commutedists.melt)
write.csv(northsouth.commutedists.melt, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/northsouth.commuteDists.impervmaxentmasks.poppairs.bestlayer.melt.fixed.csv")
View(northsouth.commutedists.melt)

#hist(northsouth.commutedists.melt$resis.dist, breaks=1000, xlim=c(0,10000000), ylim=c(0,100))
#quantile(northsouth.commutedists.melt$resis.dist, probs=c(0, 0.025, 0.05, 0.95, 0.975, 1))

#northsouth.commutedists.melt$resis.dist[northsouth.commutedists.melt$resis.dist > 10000000] <- NA
#write.csv(northsouth.commutedists.melt, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/northsouth.commuteDists.impervmaxentmasks_NAs.melt.csv")

#View(northsouth.commutedists.melt)

#ns.comm.m <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/northsouth.commuteDists.impervmaxentmasks.melt.csv")
#ns.comm.m$maskedby <- factor(ns.comm.m$maskedby, levels=c("maxent10", "maxent10imp50", "maxent10imp20"))

ns.comm.m.bestlayer <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/northsouth.commuteDists.impervmaxentmasks.poppairs.bestlayerperclade.melt.csv",
                                header=T,
                                stringsAsFactors = F)
ns.comm.m.bestlayer$mask <- factor(ns.comm.m.bestlayer$mask, levels=c("none", "imperv>50%", "imperv>20%"))


ppp <- ggplot(data = ns.comm.m.bestlayer, aes(x=layer, y=resis.dist)) 
ppp <- ppp + geom_boxplot(aes(fill=mask)) + scale_fill_manual(values=parula(100)[c(15,85,50,16,17,80,82,85)])
#ppp <- ppp + geom_violin(aes(fill=maskedby)) + scale_fill_manual(values=parula(100)[c(15,85,16,17,80,82,85)])

# if you want color for points replace group with colour=Label
#ppp <- ppp + geom_point(aes(y=resis.dist, group=clade), position = position_dodge(width=0.75))
ppp <- ppp + facet_wrap( ~ clade, scales="free") + 
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
ppp <- ppp + ggtitle("North-South SPHA resistance distances")
ppp <- ppp + guides(fill=guide_legend(title="mask"))
ppp

?ggpaired
ggpairedplotted <- ggpaired(data=ns.comm.m.bestlayer[ns.comm.m.bestlayer$clade=="SPHASouth",],
                            x="mask",
                            y="resis.dist",
                            id="populationpair",
                            fill="mask",
                            repel=TRUE) + facet_wrap( ~ clade, scales="free") # + stat_compare_means(method="wilcox.test")
ggpairedplotted
compare_means(resis.dist ~ mask, data=ggpairedplotted$data, p.adjust.method="BY", paired=TRUE, method.args=list(conf.int=TRUE), method="wilcox.test")


ggpairedplotted2 <- ggpaired(data=ns.comm.m.bestlayer,
                             x="mask",
                             y="resis.dist",
                             id="populationpair",
                             fill="mask",
                             facet.by="clade",
                             repel=TRUE,
                             line.color = "black",
                             #point.size=0.6,
                             #outlier.size=0.6,
                             line.size=0.05,
                             yscale="log10",
                             ylab="Resistance distance in meters (log scale)",
                             xlab="mask",
                             palette=parula(4),
                             ggtheme=theme_gray()) + facet_wrap( ~ clade, scales="free") #+ stat_compare_means(method="wilcox.test")
ggpairedplotted2
testcompare <- compare_means(resis.dist ~ mask, group.by="clade", data=ggpairedplotted2$data, p.adjust.method="BY", paired=TRUE, method="wilcox.test")
testcompare
ns.comm.m.bestlayer.south <- subset(ns.comm.m.bestlayer, clade=="SPHASouth", drop=TRUE)
pairedwilcox.resisdist.south <- pairwise.wilcox.test(x=ns.comm.m.bestlayer.south$resis.dist, g=ns.comm.m.bestlayer.south$mask, p.adjust.method="BY", paired=TRUE, exact=FALSE)
pairedwilcox.resisdist.south

ns.comm.m.bestlayer.north <- subset(ns.comm.m.bestlayer, clade=="SPHANorth", drop=TRUE)
pairedwilcox.resisdist.north <- pairwise.wilcox.test(x=ns.comm.m.bestlayer.north$resis.dist, g=ns.comm.m.bestlayer.north$mask, p.adjust.method="BY", paired=TRUE, exact=FALSE)
# if differences are zero/there are ties, must set exact=FALSE
pairedwilcox.resisdist.north


testcompare
ggsave(plot=ggpairedplotted2, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/ggpaired_resistance.impervmask.northsouth.eps", device="eps",
       dpi=320, width=10, height=5, units="in")

# when set to paired, doesnt show value for north imperv20 vs 50 - I think because th pvalue is 1

#ppp + stat_compare_means(label="p.format", resis.dist ~ maskedby, group.by="clade", data=ppp$data, p.adjust.method="BY", paired=FALSE, method.args=list(conf.int=TRUE), method="wilcox.test")
#uu + stat_compare_means(method="wilcox.test", aes(label=..p.adj..))

compare_means(resis.dist ~ mask, group.by="clade", data=ppp$data, p.adjust.method="BY", paired=TRUE, method.args=list(conf.int=TRUE), method="wilcox.test")
# I'm honestly not sure if these need p-value adjustments... I want to say these are independent comparisons so they don't,
# while looking at the correlations among values within north and south are not independent and so do need p.adjust?
# for intrapop with the "masking" as different "dosages" for the same treatment, should probbaly use a paired test? But need to make sure the 
# values are actually paired in the input, in that case... as long as I didn't sort the values, they should be paired, I think
# or should i do a kruskal.test if I have multiple comparisons in each clade (e.g. imp20, imp50, imp100, etc.? Maybe ask what Brad thinks...
# nah, stick with paired wilcox test (it's not really multiple comparisons since they're all simulated dosages/treatments?)
# looks like it is smart enough to know they're paired even without labels of the pairs? 

ggplot_build(ppp)$data # use to look up medians plotted on the map - "middle" column






































pp <- ggplot(data = ns.comm.m, aes(x=clade, y=resis.dist)) 
pp <- pp + geom_boxplot(aes(fill=maskedby)) + scale_fill_manual(values=parula(100)[c(15,85,50,16,17,80,82,85)])
#pp <- pp + geom_violin(aes(fill=maskedby)) + scale_fill_manual(values=parula(100)[c(15,85,16,17,80,82,85)])

# if you want color for points replace group with colour=Label
#pp <- pp + geom_point(aes(y=resis.dist, group=clade), position = position_dodge(width=0.75))
pp <- pp + facet_wrap( ~ layer, scales="free") + 
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
pp <- pp + ggtitle("North-South SPHA resistance distances")
pp <- pp + guides(fill=guide_legend(title="maskedby"))
pp


tt <- ggplot(data = ns.comm.m, aes(x=maskedby, y=resis.dist)) 
tt <- tt + geom_boxplot(aes(fill=clade)) + scale_fill_manual(values=parula(100)[c(15,85,16,17,80,82,85)])
#tt <- tt + geom_violin(aes(fill=maskedby)) + scale_fill_manual(values=parula(100)[c(15,85,16,17,80,82,85)])

# if you want color for points replace group with colour=Label
#tt <- tt + geom_point(aes(y=resis.dist, group=clade), position = position_dodge(width=0.75))
tt <- tt + facet_wrap( ~ layer, scales="free") + 
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
tt <- tt + ggtitle("North-South SPHA resistance distances")
tt <- tt + guides(fill=guide_legend(title="Clade"))
tt


#rr + stat_compare_means(method="wilcox.test", label="p.format")


#### only include the best resistnace surface for each clade: cmiRGA for south, maxentRGA for north. Make it simpler.

northsouth.commutedists.melt
colnames(northsouth.commutedists.melt)

ns.comm.m.bestlayer <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/northsouth.commuteDists.impervmaxentmasks.poppairs.bestlayerperclade.melt.csv",
         header=T,
         stringsAsFactors = F)
ns.comm.m.bestlayer$mask <- factor(ns.comm.m.bestlayer$mask, levels=c("none", "imperv>50%", "imperv>20%"))


ppp <- ggplot(data = ns.comm.m.bestlayer, aes(x=layer, y=resis.dist)) 
ppp <- ppp + geom_boxplot(aes(fill=mask)) + scale_fill_manual(values=parula(100)[c(15,85,50,16,17,80,82,85)])
#ppp <- ppp + geom_violin(aes(fill=maskedby)) + scale_fill_manual(values=parula(100)[c(15,85,16,17,80,82,85)])

# if you want color for points replace group with colour=Label
#ppp <- ppp + geom_point(aes(y=resis.dist, group=clade), position = position_dodge(width=0.75))
ppp <- ppp + facet_wrap( ~ clade, scales="free") + 
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
ppp <- ppp + ggtitle("North-South SPHA resistance distances")
ppp <- ppp + guides(fill=guide_legend(title="mask"))
ppp

?ggpaired
ggpairedplotted <- ggpaired(data=ns.comm.m.bestlayer[ns.comm.m.bestlayer$clade=="SPHASouth",],
         x="mask",
         y="resis.dist",
         id="populationpair",
         fill="mask",
         repel=TRUE) + facet_wrap( ~ clade, scales="free") # + stat_compare_means(method="wilcox.test")
ggpairedplotted
compare_means(resis.dist ~ mask, data=ggpairedplotted$data, p.adjust.method="BY", paired=TRUE, method.args=list(conf.int=TRUE), method="wilcox.test")


ggpairedplotted2 <- ggpaired(data=ns.comm.m.bestlayer,
         x="mask",
         y="resis.dist",
         id="populationpair",
         fill="mask",
         facet.by="clade",
         repel=TRUE,
         line.color = "black",
         #point.size=0.6,
         #outlier.size=0.6,
         line.size=0.05,
         yscale="log10",
         ylab="Resistance distance in meters (log scale)",
         xlab="mask",
         palette=parula(4),
         ggtheme=theme_gray()) + facet_wrap( ~ clade, scales="free") #+ stat_compare_means(method="wilcox.test")
ggpairedplotted2
testcompare <- compare_means(resis.dist ~ mask, group.by="clade", data=ggpairedplotted2$data, p.adjust.method="BY", paired=TRUE, method="wilcox.test")
testcompare
ns.comm.m.bestlayer.south <- subset(ns.comm.m.bestlayer, clade=="SPHASouth", drop=TRUE)
pairedwilcox.resisdist.south <- pairwise.wilcox.test(x=ns.comm.m.bestlayer.south$resis.dist, g=ns.comm.m.bestlayer.south$mask, p.adjust.method="BY", paired=TRUE, exact=FALSE)
pairedwilcox.resisdist.south

ns.comm.m.bestlayer.north <- subset(ns.comm.m.bestlayer, clade=="SPHANorth", drop=TRUE)
pairedwilcox.resisdist.north <- pairwise.wilcox.test(x=ns.comm.m.bestlayer.north$resis.dist, g=ns.comm.m.bestlayer.north$mask, p.adjust.method="BY", paired=TRUE, exact=FALSE)
# if differences are zero/there are ties, must set exact=FALSE
pairedwilcox.resisdist.north


testcompare
ggsave(plot=ggpairedplotted2, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/ggpaired_resistance.impervmask.northsouth.eps", device="eps",
       dpi=320, width=10, height=5, units="in")

# when set to paired, doesnt show value for north imperv20 vs 50 - I think because th pvalue is 1

#ppp + stat_compare_means(label="p.format", resis.dist ~ maskedby, group.by="clade", data=ppp$data, p.adjust.method="BY", paired=FALSE, method.args=list(conf.int=TRUE), method="wilcox.test")
#uu + stat_compare_means(method="wilcox.test", aes(label=..p.adj..))

compare_means(resis.dist ~ mask, group.by="clade", data=ppp$data, p.adjust.method="BY", paired=TRUE, method.args=list(conf.int=TRUE), method="wilcox.test")
# I'm honestly not sure if these need p-value adjustments... I want to say these are independent comparisons so they don't,
# while looking at the correlations among values within north and south are not independent and so do need p.adjust?
# for intrapop with the "masking" as different "dosages" for the same treatment, should probbaly use a paired test? But need to make sure the 
# values are actually paired in the input, in that case... as long as I didn't sort the values, they should be paired, I think
# or should i do a kruskal.test if I have multiple comparisons in each clade (e.g. imp20, imp50, imp100, etc.? Maybe ask what Brad thinks...
# nah, stick with paired wilcox test (it's not really multiple comparisons since they're all simulated dosages/treatments?)
# looks like it is smart enough to know they're paired even without labels of the pairs? 

ggplot_build(ppp)$data # use to look up medians plotted on the map - "middle" column

(13132027.1-405134.1)/405134.1
(467312.7-459934.4)/459934.4
### remove pairwise comparisons that are disconnected in *any* raster (the way circuitscape would treat them)
north.commutedists.df.intact <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/north.commuteDists.impervmaxentmasks.intactpopsonly.df.csv",
                                         header=T,
                                         stringsAsFactors = F)
south.commutedists.df.intact <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/south.commuteDists.impervmaxentmasks.intactpopsonly.df.csv",
                                         header=T,
                                         stringsAsFactors = F)
View(south.commutedists.df.intact)
south.commutedists.melt.intact <- melt(south.commutedists.df.intact, value.name="resis.dist")
View(south.commutedists.melt.intact)
south.commutedists.melt.intact <- cbind(rep("SPHASouth", n=length(south.commutedists.melt.intact$resis.dist)), south.commutedists.melt.intact)
colnames(south.commutedists.melt.intact) <- c("clade", "layer", "resis.dist")

north.commutedists.melt.intact <- melt(north.commutedists.df.intact, value.name="resis.dist")
View(north.commutedists.melt.intact)
north.commutedists.melt.intact <- cbind(rep("SPHANorth", n=length(north.commutedists.melt.intact$resis.dist)), north.commutedists.melt.intact)
colnames(north.commutedists.melt.intact) <- c("clade", "layer", "resis.dist")
northsouth.commutedists.melt.intact <- dplyr::bind_rows(north.commutedists.melt.intact, south.commutedists.melt.intact)
View(northsouth.commutedists.melt.intact)
write.csv(northsouth.commutedists.melt.intact, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/northsouth.commuteDists.impervmaxentmasks.melt.intact.csv")

ns.comm.m.intact <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/northsouth.commuteDists.impervmaxentmasks.melt.intact.csv",
                             header=T,
                             stringsAsFactors = F)
ns.comm.m.intact$maskedby <- factor(ns.comm.m.intact$maskedby, levels=c("maxent10", "maxent10imp50", "maxent10imp20"))


ss <- ggplot(data = ns.comm.m.intact, aes(x=clade, y=resis.dist)) 
ss <- ss + geom_boxplot(aes(fill=maskedby)) + scale_fill_manual(values=parula(100)[c(15,85,50,17,80,82,95)])
#ss <- ss + geom_violin(aes(fill=maskedby)) + scale_fill_manual(values=parula(100)[c(15,85,16,17,80,82,85)])

# if you want color for points replace group with colour=Label
#ss <- ss + geom_point(aes(y=resis.dist, group=clade), position = position_dodge(width=0.75))
ss <- ss + facet_wrap( ~ layer, scales="free") + 
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
ss <- ss + ggtitle("North-South SPHA resistance distances")
ss <- ss + guides(fill=guide_legend(title="maskedby"))
ss





























##### different parameterization of impervious layers: maximum resistance value within the layer * 10

southcommutedists.nax10 <- list()
northcommutedists.nax10 <- list()

plot(northresistancelayers2, col=parula(20))
plot(southresistancelayers2, col=parula(20))

for (i in 1:length(names(southresistancelayers2))){
  transsouth <- transition(southresistancelayers2[[i]],
                           transitionFunction = function(x) 1/mean(x),
                           directions=8)
  transgeo <- geoCorrection(transsouth, type="r")
  commdistout <- commuteDistance(transgeo, coords=spdf.sphasouth178spatial.pops)
  southcommutedists.nax10[[i]] <- as.matrix(commdistout)
  rownames(southcommutedists.nax10[[1]]) <- sphasouth178.rgacoords$SPHASpoplvl2
  colnames(southcommutedists.nax10[[1]]) <- sphasouth178.rgacoords$SPHASpoplvl2
  names(southcommutedists.nax10)[i] <- paste0(names(southresistancelayers2[[i]]),"_","commuteDistance")
}


for (i in 1:length(names(northresistancelayers2))){
  transnorth <- transition(northresistancelayers2[[i]],
                           transitionFunction = function(x) 1/mean(x),
                           directions=8)
  transgeo <- geoCorrection(transnorth, type="r")
  commdistout <- commuteDistance(transgeo, coords=spdf.sphanorth124spatial.pops)
  northcommutedists.nax10[[i]] <- as.matrix(commdistout)
  rownames(northcommutedists.nax10[[1]]) <- sphanorth124.rgacoords$poplevel2
  colnames(northcommutedists.nax10[[1]]) <- sphanorth124.rgacoords$poplevel2
  names(northcommutedists.nax10)[i] <- paste0(names(northresistancelayers2[[i]]),"_","commuteDistance")
}

#save old settings
op <- par(no.readonly = TRUE)
#change settings
par(mar=c(8, 4, 2, 2) + 0.1)


south.commutedists.nax10.df <- NULL
south.commutedists.nax10.df$cmiRGA.maxent10 <- lower(southcommutedists.nax10$cmiRGA_SPHAsouth_noimperviousmask_maxentmask10_commuteDistance)
south.commutedists.nax10.df$cmiRGA.maxent10.imp50 <- lower(southcommutedists.nax10$cmiRGA_SPHAsouth_imperviousmask50_maxentmask10_commuteDistance)
south.commutedists.nax10.df$cmiRGA.maxent10.imp20 <- lower(southcommutedists.nax10$cmiRGA_SPHAsouth_imperviousmask20_maxentmask10_commuteDistance)
south.commutedists.nax10.df$maxentRGA.maxent10 <- lower(southcommutedists.nax10$maxentRGA_SPHAsouth_noimperviousmask_maxentmask10_commuteDistance)
south.commutedists.nax10.df$maxentRGA.maxent10.imp50 <- lower(southcommutedists.nax10$maxentRGA_SPHAsouth_imperviousmask50_maxentmask10_commuteDistance)
south.commutedists.nax10.df$maxentRGA.maxent10.imp20 <- lower(southcommutedists.nax10$maxentRGA_SPHAsouth_imperviousmask20_maxentmask10_commuteDistance)
south.commutedists.nax10.df$topowetRGA.maxent10 <- lower(southcommutedists.nax10$topowetRGA_SPHAsouth_noimperviousmask_maxentmask10_commuteDistance)
south.commutedists.nax10.df$topowetRGA.maxent10.imp50 <- lower(southcommutedists.nax10$topowetRGA_SPHAsouth_imperviousmask50_maxentmask10_commuteDistance)
south.commutedists.nax10.df$topowetRGA.maxent10.imp20 <- lower(southcommutedists.nax10$topowetRGA_SPHAsouth_imperviousmask20_maxentmask10_commuteDistance)

south.commutedists.nax10.melt <- melt(south.commutedists.nax10.df, value.name="resis.dist")
View(south.commutedists.nax10.melt)
south.commutedists.nax10.melt <- cbind(rep("SPHASouth", n=length(south.commutedists.nax10.melt$resis.dist)), south.commutedists.nax10.melt)
colnames(south.commutedists.nax10.melt) <- c("clade", "resis.dist", "layer")

write.csv(north.commutedists.nax10.df, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/north.commuteDists.nax10.impervmaxentmasks.df.csv")
write.csv(south.commutedists.nax10.df, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/south.commuteDists.nax10.impervmaxentmasks.df.csv")



north.commutedists.nax10.df <- NULL
north.commutedists.nax10.df$cmiRGA.maxent10 <- lower(northcommutedists.nax10$cmiRGA_SPHAnorth_noimperviousmask_maxentmask10_commuteDistance)
north.commutedists.nax10.df$cmiRGA.maxent10.imp50 <- lower(northcommutedists.nax10$cmiRGA_SPHAnorth_imperviousmask50_maxentmask10_commuteDistance)
north.commutedists.nax10.df$cmiRGA.maxent10.imp20 <- lower(northcommutedists.nax10$cmiRGA_SPHAnorth_imperviousmask20_maxentmask10_commuteDistance)
north.commutedists.nax10.df$maxentRGA.maxent10 <- lower(northcommutedists.nax10$maxentRGA_SPHAnorth_noimperviousmask_maxentmask10_commuteDistance)
north.commutedists.nax10.df$maxentRGA.maxent10.imp50 <- lower(northcommutedists.nax10$maxentRGA_SPHAnorth_imperviousmask50_maxentmask10_commuteDistance)
north.commutedists.nax10.df$maxentRGA.maxent10.imp20 <- lower(northcommutedists.nax10$maxentRGA_SPHAnorth_imperviousmask20_maxentmask10_commuteDistance)

north.commutedists.nax10.melt <- melt(north.commutedists.nax10.df, value.name="resis.dist")
View(north.commutedists.nax10.melt)
north.commutedists.nax10.melt <- cbind(rep("SPHANorth", n=length(north.commutedists.nax10.melt$resis.dist)), north.commutedists.nax10.melt)
colnames(north.commutedists.nax10.melt) <- c("clade", "resis.dist", "layer")

northsouth.commutedists.nax10.melt <- dplyr::bind_rows(north.commutedists.nax10.melt, south.commutedists.nax10.melt)

write.csv(northsouth.commutedists.nax10.melt, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/northsouth.commuteDists.nax10.impervmaxentmasks.melt.csv")

hist(northsouth.commutedists.nax10.melt$resis.dist, breaks=1000, xlim=c(0,10000000), ylim=c(0,100))
quantile(northsouth.commutedists.nax10.melt$resis.dist, probs=c(0, 0.025, 0.05, 0.95, 0.975, 1))

#northsouth.commutedists.nax10.melt$resis.dist[northsouth.commutedists.nax10.melt$resis.dist > 10000000] <- NA
write.csv(northsouth.commutedists.nax10.melt, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/northsouth.commuteDists.nax10.impervmaxentmasks_NAs.melt.csv")

View(northsouth.commutedists.nax10.melt)

ns.comm.m.nax10 <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/northsouth.commuteDists.nax10.impervmaxentmasks.melt.csv")
ns.comm.m.nax10$maskedby <- factor(ns.comm.m.nax10$maskedby, levels=c("maxent10", "maxent10imp50", "maxent10imp20"))


rr <- ggplot(data = ns.comm.m.nax10, aes(x=maskedby, y=resis.dist)) 
rr <- rr + geom_boxplot(aes(fill=clade)) + scale_fill_manual(values=parula(100)[c(15,85,16,17,80,82,85)])
#rr <- rr + geom_violin(aes(fill=maskedby)) + scale_fill_manual(values=parula(100)[c(15,85,16,17,80,82,85)])

# if you want color for points replace group with colour=Label
#rr <- rr + geom_point(aes(y=resis.dist, group=clade), position = position_dodge(width=0.75))
rr <- rr + facet_wrap( ~ layer, scales="free") + 
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
rr <- rr + ggtitle("North-South SPHA resistance distances")
rr <- rr + guides(fill=guide_legend(title="maskedby"))
rr
#rr + stat_compare_means(method="wilcox.test", label="p.format")

rrr <- ggplot(data = ns.comm.m.nax10, aes(x=clade, y=resis.dist)) 
rrr <- rrr + geom_boxplot(aes(fill=maskedby)) + scale_fill_manual(values=parula(100)[c(15,85,50,16,17,80,82,85)])
#rrr <- rrr + geom_violin(aes(fill=maskedby)) + scale_fill_manual(values=parula(100)[c(15,85,16,17,80,82,85)])

# if you want color for points replace group with colour=Label
#rrr <- rrr + geom_point(aes(y=resis.dist, group=clade), position = position_dodge(width=0.75))
rrr <- rrr + facet_wrap( ~ layer, scales="free") + 
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
rrr <- rrr + ggtitle("North-South SPHA resistance distances")
rrr <- rrr + guides(fill=guide_legend(title="maskedby"))
rrr
#rrr + stat_compare_means(method="wilcox.test", label="p.format")


### remove pairwise comparisons that are disconnected in *any* raster (the way circuitscape would treat them)
north.commutedists.nax10.df.intact <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/north.commuteDists.nax10.impervmaxentmasks.intactpopsonly.df.csv",
                                         header=T,
                                         stringsAsFactors = F)
south.commutedists.nax10.df.intact <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/south.commuteDists.nax10.impervmaxentmasks.intactpopsonly.df.csv",
                                         header=T,
                                         stringsAsFactors = F)
View(south.commutedists.nax10.df.intact)
south.commutedists.nax10.melt.intact <- melt(south.commutedists.nax10.df.intact, value.name="resis.dist")
View(south.commutedists.nax10.melt.intact)
south.commutedists.nax10.melt.intact <- cbind(rep("SPHASouth", n=length(south.commutedists.nax10.melt.intact$resis.dist)), south.commutedists.nax10.melt.intact)
colnames(south.commutedists.nax10.melt.intact) <- c("clade", "layer", "resis.dist")

north.commutedists.nax10.melt.intact <- melt(north.commutedists.nax10.df.intact, value.name="resis.dist")
View(north.commutedists.nax10.melt.intact)
north.commutedists.nax10.melt.intact <- cbind(rep("SPHANorth", n=length(north.commutedists.nax10.melt.intact$resis.dist)), north.commutedists.nax10.melt.intact)
colnames(north.commutedists.nax10.melt.intact) <- c("clade", "layer", "resis.dist")
northsouth.commutedists.nax10.melt.intact <- dplyr::bind_rows(north.commutedists.nax10.melt.intact, south.commutedists.nax10.melt.intact)
View(northsouth.commutedists.nax10.melt.intact)
write.csv(northsouth.commutedists.nax10.melt.intact, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/northsouth.commuteDists.nax10.impervmaxentmasks.melt.intact.csv")

ns.comm.m.intact <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/northsouth.commuteDists.nax10.impervmaxentmasks.melt.intact.csv",
                             header=T,
                             stringsAsFactors = F)
ns.comm.m.intact$maskedby <- factor(ns.comm.m.intact$maskedby, levels=c("maxent10", "maxent10imp50", "maxent10imp20"))


ss <- ggplot(data = ns.comm.m.intact, aes(x=clade, y=resis.dist)) 
ss <- ss + geom_boxplot(aes(fill=maskedby)) + scale_fill_manual(values=parula(100)[c(15,85,50,17,80,82,95)])
#ss <- ss + geom_violin(aes(fill=maskedby)) + scale_fill_manual(values=parula(100)[c(15,85,16,17,80,82,85)])

# if you want color for points replace group with colour=Label
#ss <- ss + geom_point(aes(y=resis.dist, group=clade), position = position_dodge(width=0.75))
ss <- ss + facet_wrap( ~ layer, scales="free") + 
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
ss <- ss + ggtitle("North-South SPHA resistance distances")
ss <- ss + guides(fill=guide_legend(title="maskedby"))
ss

sss <- ggplot(data = ns.comm.m.intact, aes(x=clade, y=resis.dist)) 
sss <- sss + geom_boxplot(aes(fill=maskedby)) + scale_fill_manual(values=parula(100)[c(15,85,50,17,80,82,95)])
#sss <- sss + geom_violin(aes(fill=maskedby)) + scale_fill_manual(values=parula(100)[c(15,85,16,17,80,82,85)])

# if you want color for points replace group with colour=Label
#sss <- sss + geom_point(aes(y=resis.dist, group=clade), position = position_dodge(width=0.75))
sss <- sss + facet_wrap( ~ layer, scales="free") + 
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
sss <- sss + ggtitle("North-South SPHA resistance distances")
sss <- sss + guides(fill=guide_legend(title="maskedby"))
sss
