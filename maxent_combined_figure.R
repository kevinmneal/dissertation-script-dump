### sum and normalize maxent models of north and south


maxentstack <- stack(c("G:/My Drive/CHELSA_1.2/SPHA_N_20190308_all49cropped_redo3_nolandcover_cv10/SPHA_N_avg.asc",
                        "G:/My Drive/CHELSA_1.2/SPHA_S_20190308_all49cropped_redo3_nolandcover_cv10/SPHA_S_avg.asc"))
maxentstack2 <- maxentstack

maxentstack2$SPHA_N_avg <- setMinMax(maxentstack$SPHA_N_avg)
maxentstack2$SPHA_S_avg <- setMinMax(maxentstack$SPHA_S_avg)

#maxentstack2$SPHA_N_avg[maxentstack2$SPHA_N_avg<0.1] <- 0
#maxentstack2$SPHA_N_avg[maxentstack2$SPHA_N_avg<0.1] <- 0

maxentsum2 <- sum(maxentstack)
maxentsum2 <- maxentsum2/maxValue(maxentsum2)
plot(maxentsum2, col=parula(100))
maxentsum3 <- maxentsum2
maxentsum3[maxentsum3<0.1] <- 0
plot(maxentsum3, col=c("black", parula(100)))


sphanorth.maxentcoords <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/Maxent_20190302/SPHA_N_maxent_coords_20190302.csv")
sphasouth.maxentcoords <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/Maxent_20190302/SPHA_S_maxent_coords_20190302.csv")

sphamaxentcoords <- rbind(sphanorth.maxentcoords, sphasouth.maxentcoords)[,2:3]

sphamaxentcoords.sub <- cbind(sphamaxentcoords, raster::extract(maxentsum2, SpatialPoints(sphamaxentcoords)))
sphamaxentcoords.sub <- sphamaxentcoords.sub[complete.cases(sphamaxentcoords.sub),]
View(sphamaxentcoords.sub)
#sphamaxentcoords.sub[-c(9,138),] # these points stand out as seeming unreasonable... remove?
sphamaxentcoords.sub <- SpatialPoints(sphamaxentcoords.sub)
popcoords.sphanorth124spatial.unique.alphabetical <- readRDS("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphanorth124spatial/popcoords.sphanorth124spatial.unique.alphabetical.rds")


CairoPNG(filename = "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphanorthsouth.maxent.withpoints2.png",
    height=9, width=8, units="in", res=320)
plot(maxentsum2, col=c("black",viridis(100)), main="Maxent habitat suitability")
points(sphamaxentcoords.sub, col="black", bg="wheat", pch=21, cex=0.8)
#points(sphanorth124.rgacoords[,2:3], col="white", pch=21, cex=0.5)
#points(sphasouth178.rgacoords[,2:3], col="white", pch=21, cex=0.5)
points(popcoords.sphasouth178spatial.unique.alphabetical[,2:3], col="black", bg="firebrick2", pch=21, cex=0.6)
points(popcoords.sphanorth124spatial.unique.alphabetical[,2:3], col="black", bg="dodgerblue", pch=21, cex=0.6)
#points(popcoords.sphasouth178spatial.unique.alphabetical[,2:3], col="firebrick2", pch=20, cex=0.5)
#points(popcoords.sphanorth124spatial.unique.alphabetical[,2:3], col="firebrick2", pch=20, cex=0.5)
#points(sphamaxentcoords.sub, col="red", pch=20, cex=0.6)

dev.off()


ggsave(device="png", filename="sphanorthsouth.maxent.withpoints.png",
       height=9, width=7, dpi=320)


sphamaxentcoords.sub



library(phangorn)
library(ape)
library(phytools)


?read.phylo
