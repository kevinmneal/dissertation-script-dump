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
south.commutedists.df$impervious50 <- lower(southcommutedists$cmiRGA_SPHAsouth_imperviousmask50_maxentmask10_commuteDistance)
south.commutedists.df$impervious20 <- lower(southcommutedists$cmiRGA_SPHAsouth_imperviousmask20_maxentmask10_commuteDistance)
#south.commutedists.df$cmiRGA.maxent10 <- lower(southcommutedists$cmiRGA_SPHAsouth_noimperviousmask_maxentmask10_commuteDistance)
#south.commutedists.df$cmiRGA.maxent10.impervious50 <- lower(southcommutedists$cmiRGA_SPHAsouth_imperviousmask50_maxentmask10_commuteDistance)
#south.commutedists.df$cmiRGA.maxent10.impervious20 <- lower(southcommutedists$cmiRGA_SPHAsouth_imperviousmask20_maxentmask10_commuteDistance)
#south.commutedists.df$maxentRGA.maxent10 <- lower(southcommutedists$maxentRGA_SPHAsouth_noimperviousmask_maxentmask10_commuteDistance)
#south.commutedists.df$maxentRGA.maxent10.impervious50 <- lower(southcommutedists$maxentRGA_SPHAsouth_imperviousmask50_maxentmask10_commuteDistance)
#south.commutedists.df$maxentRGA.maxent10.impervious20 <- lower(southcommutedists$maxentRGA_SPHAsouth_imperviousmask20_maxentmask10_commuteDistance)
#south.commutedists.df$topowetRGA.maxent10 <- lower(southcommutedists$topowetRGA_SPHAsouth_noimperviousmask_maxentmask10_commuteDistance)
#south.commutedists.df$topowetRGA.maxent10.impervious50 <- lower(southcommutedists$topowetRGA_SPHAsouth_imperviousmask50_maxentmask10_commuteDistance)
#south.commutedists.df$topowetRGA.maxent10.impervious20 <- lower(southcommutedists$topowetRGA_SPHAsouth_imperviousmask20_maxentmask10_commuteDistance)
#south.commutedists.df$populationpair <- factor(c(1:length(south.commutedists.df$cmiRGA.maxent10)))
south.commutedists.df$populationpair <- factor(southpairmelt.concat)
south.commutedists.df$clade <- c(rep("SPHASouth", times=length(south.commutedists.df$impervious50)))
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
north.commutedists.df$impervious50 <- lower(northcommutedists$maxentRGA_SPHAnorth_imperviousmask50_maxentmask10_commuteDistance)
north.commutedists.df$impervious20 <- lower(northcommutedists$maxentRGA_SPHAnorth_imperviousmask20_maxentmask10_commuteDistance)
#north.commutedists.df$maxentRGA.maxent10 <- lower(northcommutedists$maxentRGA_SPHAnorth_noimperviousmask_maxentmask10_commuteDistance)
#north.commutedists.df$maxentRGA.maxent10.impervious50 <- lower(northcommutedists$maxentRGA_SPHAnorth_imperviousmask50_maxentmask10_commuteDistance)
#north.commutedists.df$maxentRGA.maxent10.impervious20 <- lower(northcommutedists$maxentRGA_SPHAnorth_imperviousmask20_maxentmask10_commuteDistance)
#north.commutedists.df$maxentRGA.maxent10 <- lower(northcommutedists$maxentRGA_SPHAnorth_noimperviousmask_maxentmask10_commuteDistance)
#north.commutedists.df$maxentRGA.maxent10.impervious50 <- lower(northcommutedists$maxentRGA_SPHAnorth_imperviousmask50_maxentmask10_commuteDistance)
#north.commutedists.df$maxentRGA.maxent10.impervious20 <- lower(northcommutedists$maxentRGA_SPHAnorth_imperviousmask20_maxentmask10_commuteDistance)
#north.commutedists.df$topowetRGA.maxent10 <- lower(northcommutedists$topowetRGA_SPHAnorth_noimperviousmask_maxentmask10_commuteDistance)
#north.commutedists.df$topowetRGA.maxent10.impervious50 <- lower(northcommutedists$topowetRGA_SPHAnorth_imperviousmask50_maxentmask10_commuteDistance)
#north.commutedists.df$topowetRGA.maxent10.impervious20 <- lower(northcommutedists$topowetRGA_SPHAnorth_imperviousmask20_maxentmask10_commuteDistance)
#north.commutedists.df$populationpair <- factor(c(1:length(north.commutedists.df$maxentRGA.maxent10)))
north.commutedists.df$populationpair <- factor(northpairmelt.concat)
north.commutedists.df$clade <- c(rep("SPHANorth", times=length(north.commutedists.df$impervious50)))
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
#ns.comm.m$maskedby <- factor(ns.comm.m$maskedby, levels=c("maxent10", "maxent10impervious50", "maxent10impervious20"))





# ns.comm.m.bestlayer <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/northsouth.commuteDists.impervmaxentmasks.poppairs.bestlayer.melt.fixed.csv",
#                                 header=T,
#                                 stringsAsFactors = F)
ns.comm.m.bestlayer <- northsouth.commutedists.melt
#View(northsouth.commutedists.melt)
ns.comm.m.bestlayer$mask <- factor(ns.comm.m.bestlayer$mask, levels=c("none", "impervious50", "impervious20"))


ppp <- ggplot(data = ns.comm.m.bestlayer, aes(x="mask", y=resis.dist)) 
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
ggsave(plot=ggpairedplotted2, filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/ggpaired_resistance.impervmask.northsouth.svg", device=svg(),
       dpi=320, width=10, height=5, units="in")

# when set to paired, doesnt show value for north imperv20 vs 50 - I think because th pvalue is 1

#ppp + stat_compare_means(label="p.format", resis.dist ~ maskedby, group.by="clade", data=ppp$data, p.adjust.method="BY", paired=FALSE, method.args=list(conf.int=TRUE), method="wilcox.test")
#uu + stat_compare_means(method="wilcox.test", aes(label=..p.adj..))

compare_means(resis.dist ~ mask, group.by="clade", data=ppp$data, p.adjust.method="BY", paired=TRUE, method.args=list(conf.int=TRUE), method="wilcox.test")
# I'm honestly not sure if these need p-value adjustments... I want to say these are independent comparisons so they don't,
# while looking at the correlations among values within north and south are not independent and so do need p.adjust?
# for intrapop with the "masking" as different "dosages" for the same treatment, should probbaly use a paired test? But need to make sure the 
# values are actually paired in the input, in that case... as long as I didn't sort the values, they should be paired, I think
# or should i do a kruskal.test if I have multiple comparisons in each clade (e.g. impervious20, impervious50, imp100, etc.? Maybe ask what Brad thinks...
# nah, stick with paired wilcox test (it's not really multiple comparisons since they're all simulated dosages/treatments?)
# looks like it is smart enough to know they're paired even without labels of the pairs? 

ggplot_build(ppp)$data # use to look up medians plotted on the map - "middle" column

# north, difference in medians from none to impervious50:
(1140935.2 - 1146806.0)/1146806.0*100 # -0.5% decrease - attribute to the idiosyncricies of random-walk-based commute distance

(13132027.1 - 405134.1)/405134.1*100 # 3141% increase






### plot masked rasters

# example from github:
# breakpoints <- c(94,100,120,140,160,180,195)
# colors <- c("red","white","white","white","white","blue")
# plot(volcanoR,breaks=breakpoints,col=colors)


# get minmax of unmasked resistance layer and set breaks in the plot
#maxValue(northresistancelayers2$maxentRGA_SPHAnorth_noimperviousmask_maxentmask10[maxentnorth>=0.1]) # 2496.366 is highest
maxentnorth.mask10.2 <- northresistancelayers2$maxentRGA_SPHAnorth_noimperviousmask_maxentmask10
maxentnorth.mask10.2[maxentnorth.mask10.2>9999] <- NA
#plot(maxentnorth.mask10.2)
maxValue(setMinMax(maxentnorth.mask10.2)) # 1201.015

plot(northresistancelayers2, 
     breaks=c(seq(from=0, to=1200, by=100),10000), #c(0,250,500,750,1000,1250,1500,1750,2000,2250,2500,10000),
     col=c(parula(length(seq(from=0, to=1200, by=100))-1), "black"), 
     nc=3, nr=1)

south.at <- seq(0,1200,100)


parulaTheme <- function (region = c(pals::parula(10), "black"), ...) 
{
  rasterTheme(region = region, ...)
}
png(filename="northresistance.cmi.imperv.png", height=5, width=10, units="in", res=320)

levelplot(northresistancelayers2[[c(3,2,1)]], par.settings=parulaTheme, 
          at=c(seq(0,1300,100),10000),
          colorkey=list(at=c(seq(0,1300,100)), labels=list(labels=c(seq(0,1200,100),"10000"), at=c(seq(0,1300,100)))),
          layout=c(3,1),
          names.attr=c("no impervious mask", "impervious>50% masked", "impervious>20% masked")) +
  layer(sp.points(spdf.sphanorth124spatial.pops, pch=5, cex=0.5, col="red"))
dev.off()

cmisouth.mask10.2 <- southresistancelayers2$cmiRGA_SPHAsouth_noimperviousmask_maxentmask10
cmisouth.mask10.2[cmisouth.mask10.2>9999] <- NA
#plot(cmisouth.mask10.2)
maxValue(setMinMax(cmisouth.mask10.2)) # 249.6816

plot(southresistancelayers2, 
     breaks=c(seq(from=0, to=250, by=25),10000), #c(0,250,500,750,1000,1250,1500,1750,2000,2250,2500,10000),
     col=c(parula(length(seq(from=0, to=250, by=25))-1), "black"), 
     nc=3, nr=1)

png(filename="southresistance.cmi.imperv.png", height=3, width=10, units="in", res=320)

levelplot(southresistancelayers2[[c(3,2,1)]], par.settings=parulaTheme, 
          at=c(seq(0,275,25),10000),
          colorkey=list(at=c(seq(0,275,25)), labels=list(labels=c(seq(0,250,25),"10000"), at=c(seq(0,275,25)))),
          layout=c(3,1),
          names.attr=c("no impervious mask", "impervious>50% masked", "impervious>20% masked")) +
  layer(sp.points(spdf.sphasouth178spatial.pops, pch=5, cex=0.5, col="red"))
dev.off()
