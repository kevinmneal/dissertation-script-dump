
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
library(gplots)


# subsampled the sphasouth dataset to get things a bit more spatially even for running in ResistanceGA;
# pooled sancan and tena into one, removed ccsp2 and 3, pooled SOBOB1 and 2
# most pops are singleton pops so this might not be the best analysis...

vcf.75pctmsng.sphanorth124spatial <- read.vcfR("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/sphanorth124spatial_75pctmsng_oneSNPmac2_20190302.recode.vcf")
vcf.50pctmsng.sphanorth124spatial <- read.vcfR("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/sphanorth124spatial_50pctmsng_oneSNPmac2_20190302.recode.vcf")


popcoords.sphanorth124spatial <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/sphanorth124spatial_popcoords.csv")
#saveRDS(popcoords.OCreduced148, file="G:/My Drive/WORKING FOLDER Jan 2015/STUFF THAT ISNT NEWCAL/Spea/RAD_data/allsocal_303indv/popcoords.OCreduced148.rds")

genind.75pctmsng.sphanorth124spatial <- vcfR2genind(vcf.75pctmsng.sphanorth124spatial)
pop(genind.75pctmsng.sphanorth124spatial) <- popcoords.sphanorth124spatial$poplevel2
#genind.75pctmsng.sphanorth124spatial.k3 <- genind.75pctmsng.sphanorth124spatial
#pop(genind.75pctmsng.sphanorth124spatial.k3) <- popcoords.sphanorth124spatial$Site

gl.75pctmsng.sphanorth124spatial <- vcfR2genlight(vcf.75pctmsng.sphanorth124spatial)
ploidy(gl.75pctmsng.sphanorth124spatial) <- 2
pop(gl.75pctmsng.sphanorth124spatial) <- popcoords.sphanorth124spatial$poplevel2
#gl.75pctmsng.sphanorth124spatial.k3 <- gl.75pctmsng.sphanorth124spatial
#ploidy(gl.75pctmsng.sphanorth124spatial.k3) <- 2
#pop(gl.75pctmsng.sphanorth124spatial.k3) <- popcoords.sphanorth124spatial$K3hier

hierfstat.75pctmsng.sphanorth124spatial <- genind2hierfstat(genind.75pctmsng.sphanorth124spatial)
distance.Da.75pctmsng.sphanorth124spatial <- genet.dist(hierfstat.75pctmsng.sphanorth124spatial, method="Da") # doesn't retain pop names, ugh!
distancemat.Da.75pctmsng.sphanorth124spatial <- as.matrix(distance.Da.75pctmsng.sphanorth124spatial)
rownames(distancemat.Da.75pctmsng.sphanorth124spatial) <- as.character(unique(pop(genind.75pctmsng.sphanorth124spatial))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.75pctmsng.sphanorth124spatial) <- as.character(unique(pop(genind.75pctmsng.sphanorth124spatial)))
#diag(distancemat.Da.75pctmsng.sphanorth124spatial) <- NA

heatmap.2(distancemat.Da.75pctmsng.sphanorth124spatial, 
          col=linearl(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 75% missing data")

saveRDS(distancemat.Da.75pctmsng.sphanorth124spatial, file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/distancemat.Da.75pctmsng.sphanorth124spatial.rds")
stamppPhylip(distancemat.Da.75pctmsng.sphanorth124spatial, file="distancemat.Da.75pctmsng.sphanorth124spatial.phy.dst")



distance.jostD.75pctmsng.sphanorth124spatial <- pairwise_D(genind.75pctmsng.sphanorth124spatial, linearized=FALSE)
heatmap.2(as.matrix(distance.jostD.75pctmsng.sphanorth124spatial), 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=diana,
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Jost's D (2008, via mmod), hclust=ward.D2, up to 75% missing data")
stamppPhylip(as.matrix(distance.jostD.75pctmsng.sphanorth124spatial), 
             file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/sphanorth124spatial_jostD_75pctmsng_dist.phy.dst")

distance.wcfst.75pctmsng.sphanorth124spatial <- pairwise.WCfst(hierfstat.75pctmsng.sphanorth124spatial) # standard output is a matrix
distancemat.wcfst.75pctmsng.sphanorth124spatial <- as.matrix(distance.wcfst.75pctmsng.sphanorth124spatial)
distance.wcfst.75pctmsng.sphanorth124spatial <- as.dist(distancemat.wcfst.75pctmsng.sphanorth124spatial)
diag(distancemat.wcfst.75pctmsng.sphanorth124spatial) <- 0

heatmap.2(distancemat.wcfst.75pctmsng.sphanorth124spatial, 
          col=magma(256),
          trace="none", symm=T,
          #hclustfun=agnes,
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="sphanorth124spatial WCFst, up to 75% missing data")
stamppPhylip(distancemat.wcfst.75pctmsng.sphanorth124spatial, 
             file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/sphanorth124spatial_wcfst_75pctmsng_dist.phy.dst")





genind.50pctmsng.sphanorth124spatial <- vcfR2genind(vcf.50pctmsng.sphanorth124spatial)
pop(genind.50pctmsng.sphanorth124spatial) <- popcoords.sphanorth124spatial$poplevel2
#genind.50pctmsng.sphanorth124spatial.k3 <- genind.50pctmsng.sphanorth124spatial
#pop(genind.50pctmsng.sphanorth124spatial.k3) <- popcoords.sphanorth124spatial$Site

gl.50pctmsng.sphanorth124spatial <- vcfR2genlight(vcf.50pctmsng.sphanorth124spatial)
ploidy(gl.50pctmsng.sphanorth124spatial) <- 2
pop(gl.50pctmsng.sphanorth124spatial) <- popcoords.sphanorth124spatial$poplevel2
#gl.50pctmsng.sphanorth124spatial.k3 <- gl.50pctmsng.sphanorth124spatial
#ploidy(gl.50pctmsng.sphanorth124spatial.k3) <- 2
#pop(gl.50pctmsng.sphanorth124spatial.k3) <- popcoords.sphanorth124spatial$K3hier

hierfstat.50pctmsng.sphanorth124spatial <- genind2hierfstat(genind.50pctmsng.sphanorth124spatial)

distance.Da.50pctmsng.sphanorth124spatial <- genet.dist(hierfstat.50pctmsng.sphanorth124spatial, method="Da") # doesn't retain pop names, ugh!
distancemat.Da.50pctmsng.sphanorth124spatial <- as.matrix(distance.Da.50pctmsng.sphanorth124spatial)
rownames(distancemat.Da.50pctmsng.sphanorth124spatial) <- as.character(unique(pop(genind.50pctmsng.sphanorth124spatial))) #rownames(distancemat.neifst.50pctmsng) # pull from a matrix that does preserve pop names. genet.dist(method="Nei87") keeps pop names, so if you have one of those handy,...
colnames(distancemat.Da.50pctmsng.sphanorth124spatial) <- as.character(unique(pop(genind.50pctmsng.sphanorth124spatial)))
#diag(distancemat.Da.50pctmsng.sphanorth124spatial) <- NA

heatmap.2(distancemat.Da.50pctmsng.sphanorth124spatial, 
          col=rainbow(256),
          trace="none", symm=T,
          #hclustfun=diana)
          hclustfun=function(x)hclust(x, method="ward.D2"),
          main="Nei's Da (1983, via hierfstat), hclust=ward.D2, up to 50% missing data")

rga.popcoords <- read.csv("spha.genpop.15coords_fordistance.csv")
rga.points <- SpatialPoints(rga.popcoords[,2:3])
rga.neimatrix <- as.matrix(read.csv("spha.genpop.15pops_neidistancematrix.csv", header=T, row.names=1))
colnames(rga.neimatrix) <- rownames(rga.neimatrix)
rga.neimatrix[rga.neimatrix==0] <- NA
rga.neivector <- as.vector(rga.neimatrix)
rga.neivector <- rga.neivector[complete.cases(rga.neivector)]


corrplot.resis.spearman.north <- readRDS("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth184spatial/corrplot.resis.spearman.north.rds")
corrplot.resis.spearman.south <- readRDS("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth184spatial/corrplot.resis.spearman.south.rds")
mutinfo.north <- readRDS("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth184spatial/mutinfo.north.rds")
mutinfo.south <- readRDS("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth184spatial/mutinfo.south.rds")
mutinfo.south2 <- mutinfo.south
mutinfo.south2[mutinfo.south2<0] <- 0
heatmap.2(mutinfo.south2, symm=T, main="rga south input rasters, mutual information", col=viridis(256, direction=-1),
          mar=c(5,5), trace="none")
heatmap(mutinfo.north, symm=T, col=viridis(256, direction=-1))
heatmap.2(mutinfo.south2, symm=T, main="rga south input rasters, mutual information", col=viridis(256, direction=-1),
          mar=c(5,5), trace="none")
heatmap.2(abs(corrplot.resis.spearman.north), symm=T, main="rga north input rasters, abs spearman", col=viridis(256, direction=-1),
          mar=c(5,5), trace="none")
heatmap.2(abs(corrplot.resis.spearman.south), symm=T, main="rga south input rasters, abs spearman", col=viridis(256, direction=-1),
          mar=c(5,5), trace="none")


OCall208.cor.spearman.INLAND[corm.spearman.INLAND$p > 0.05] <- 0
spearman.south <- corrplot.resis.spearman.south
spearman.south[abs(spearman.south)<0.7] <- NA

spearman.north <- corrplot.resis.spearman.north
spearman.north[abs(spearman.north)<0.7] <- NA

corrplot(abs(spearman.south), 
         title="Spearman south abs spearman",
         mar=c(0,0,2,0),
         method="number",
         type="lower",
         tl.pos="lt")

spearman.south.removing <- spearman.south
abs(spearman.south)["bio03",] # nothing
abs(spearman.south)["bio11",] # bio01, 6, 8, continentality, elev, growingdeg, maxtempcoldest, monthcount, petcoldest, thermicity
abs(spearman.south)["bio10",] # bio01, 5, 9, gdd, petwettest
abs(spearman.south)["bio19",] # bio12, bio13, bio16
abs(spearman.south)["bio15",] # bio18
abs(spearman.south)["toporuggedness",] # slope, topowet
abs(spearman.south)["bio18",] # bio02, 14, 15, 17
abs(spearman.south)["slopedeg",] # toporuggedness, topowet
abs(spearman.south)["topoWet",] # toporuggedness
abs(spearman.south)["landcover",] # impervious weakly
abs(spearman.south)["bpsprehumanvegmodel",] # 
abs(spearman.south)["depthtobedrockrhorizon",] # 
abs(spearman.south)["climaticMoistureIndex",] # embergerQ
abs(spearman.south)["embergerQ",] # CMI
abs(spearman.south)["aridityIndexThornthwaite",] # nothing
abs(spearman.south)["bulkdensity5cm",] # nothing
abs(spearman.south)["claycontent5cm",] # nothing
abs(spearman.south)["siltcontent5cm",] # bio15 and bio18 weakly
abs(spearman.south)["thermicityIndex",] # bio01, 06, 08, 11, elev, gdd, maxtempcoldest, monthcount, petcoldest
abs(spearman.south)["impervious",] # landcover weakly
abs(spearman.south)["taxousda",] # nothing
abs(spearman.south)["canopy",] # nothing
abs(spearman.south)["continentality",] # nothing


maxtest <- raster("G:/My Drive/CHELSA_1.2/maxent3.4.1_SPHA_N_20190303_withenviremvars/SPHA_N_avg.asc")
maxtestagg2 <- aggregate(maxtest, fact=2, fun=modal)
maxtestagg3 <- aggregate(maxtest, fact=3, fun=modal)
maxtest
maxtestagg
par(mfrow=c(1,2))
plot(maxtest, col=viridis(256))
plot(maxtestagg3, col=viridis(256))
par(mfrow=c(1,1))

topotestnorth <- raster("G:/My Drive/CHELSA_1.2/topoWet_SPHANorth.asc")
topotestnorthextent <- topotestnorth
plot(topotestnorth, col=viridis(256))
topotestnorth.2 <- aggregate(topotestnorthextent, fact=2, fun=modal, na.rm=TRUE)
plot(topotestnorth.2)
topotestnorth.2[topotestnorth.2==-9999] <- NA
plot(topotestnorth.2)
plot(topotestnorth)
topotestnorth.2.useasmask2 <- aggregate(topotestnorthextent, fact=2, fun=modal, na.rm=TRUE, extend=FALSE)
plot(topotestnorth.2.useasmask2)
topotestnorth.2.useasmask2[topotestnorth.2.useasmask2==-9999] <- NA
plot(topotestnorth.3.useasmask)
topotestnorth.3again <- aggregate(topotestnorth, fact=3, fun=median, na.rm=TRUE)
plot(topotestnorth.3)
topotestnorth.3[topotestnorth.3<0] <- NA
plot(topotestnorth.3)

#topotestnorthextent <- crop(topotestnorth, topotestnorth.3)

topotestmasking3 <- mask(topotestnorth.3again, topotestnorth.3)
topotestnorthextent[is.na(topotestnorthextent)] <- -9999
plot(topotestnorthextent)
# some areas with NA values get lost with aggregation, would be good to keep those as they're barriers though, try na.rm=FALSE
# use median for continuous vars, mode for categorical
par(mfrow=c(1,3))
plot(topotestnorth.3.useasmask, col=viridis(256))
writeRaster(topotestnorth.3.useasmask, file="G:/My Drive/CHELSA_1.2/aggregated_SPHANorth_useasmask.asc", format="ascii")
plot(topotestnorth.3, col=viridis(256))
plot(topotestnorth.3again, col=viridis(256))
par(mfrow=c(1,1))
topotestnorth
topotestnorth.2
topotestnorth.3
# standard res: 1428207 pixels; agg fac=2: 357650; fac=3: ; fac=4: 89568

bedrock <- raster("G:/My Drive/CHELSA_1.2/depthtobedrockrhorizon_SPHANorth.asc")
bio18 <- raster("G:/My Drive/CHELSA_1.2/bio18_SPHANorth.asc")
elev <- raster("G:/My Drive/CHELSA_1.2/elev_SPHANorth.asc")

plot(bedrock, col=viridis(256))
plot(bio18, col=viridis(256))
plot(elev, col=viridis(256))

plot(raster("G:/My Drive/CHELSA_1.2/CHELSA_CA/CHELSA_bio01_ca.asc"))

topotestsouth <- raster("G:/My Drive/CHELSA_1.2/topoWet_SPHAsouth.asc")
topotestsouthextent <- topotestsouth
plot(topotestsouth, col=viridis(256))
topotestsouth.2.useasmask2 <- aggregate(topotestsouthextent, fact=2, fun=modal, na.rm=TRUE, extend=FALSE)
plot(topotestsouth.2)
topotestsouth.2.useasmask2[topotestsouth.2.useasmask2==-9999] <- NA
topotestsouth.2[topotestsouth.2<0] <- NA
plot(topotestsouth.2.useasmask)
plot(topotestsouth)
topotestsouth.3.useasmask <- aggregate(topotestsouthextent, fact=3, fun=modal, na.rm=TRUE)
plot(topotestsouth.3.useasmask)
topotestsouth.3.useasmask[topotestsouth.3.useasmask==-9999] <- NA
plot(topotestsouth.3.useasmask)
topotestsouth.3again <- aggregate(topotestsouth, fact=3, fun=median, na.rm=TRUE)
plot(topotestsouth.3)
topotestsouth.3[topotestsouth.3<0] <- NA
plot(topotestsouth.3)
writeRaster(topotestnorth.2.useasmask, file="G:/My Drive/CHELSA_1.2/aggregatedfac2_SPHANorth_useasmask.asc", format="ascii")

#topotestsouthextent <- crop(topotestsouth, topotestsouth.3)

topotestmasking3 <- mask(topotestsouth.3again, topotestsouth.3)
topotestsouthextent[is.na(topotestsouthextent)] <- -9999
plot(topotestsouthextent)
# some areas with NA values get lost with aggregation, would be good to keep those as they're barriers though, try na.rm=FALSE
# use median for continuous vars, mode for categorical
par(mfrow=c(1,3))
plot(topotestsouth.3.useasmask, col=viridis(256))
writeRaster(topotestsouth.2.useasmask, file="G:/My Drive/CHELSA_1.2/aggregatedfac2_SPHAsouth_useasmask.asc", format="ascii")
plot(topotestsouth.2.useasmask, col=viridis(256))
plot(topotestsouth.3again, col=viridis(256))
par(mfrow=c(1,1))
topotestsouth
topotestsouth.2
topotestsouth.3

teststack <- stack(bio18, bedrock, elev)
testagg <- mask(crop(maskobject, teststack), teststack)

testresamp <- resample(ras, maskobject, method="bilinear")
testresampmask <- mask(testresamp, maskobject)
arid <- raster("G:/My Drive/CHELSA_1.2/aridityIndexThornthwaite_SPHANorth.asc")
plot(arid)

ras <- raster("G:/My Drive/CHELSA_1.2/bio18_SPHANorth.asc")
maskobject <- raster("G:/My Drive/CHELSA_1.2/aggregatedfac2_SPHANorth_useasmask.asc")
ras.aggfac2 <- aggregate(arid, fact=2, fun=median, na.rm=TRUE)
ras.aggfac2 <- resample(arid, maskobject, method="bilinear")
ras.aggfac2 <- mask(ras.aggfac2, maskobject)


ga.test <- readRDS("G:/My Drive/CHELSA_1.2/GA.inputs.rds")
surfacetype <- ga.test$surface.type
surfacename <- ga.test$layer.names
surfacetype[ga.test$layer.names[ga.test$layer.names=="bio03SPHAsouth"]]
ga.test$layer.names[46]
surfacetype[which(surfacename=="landcoverSPHAsouth")]
surfacename[22]
surfacetype[22] <- "cat"

which(surfacename == "annualPETSPHAsouth")
maxentras <- raster("G:/My Drive/CHELSA_1.2/envlayers_south_aggfac2/maxent_SPHAsouth_aggfac2.asc")
maxentras
plot(maxentras)
landcovras <- raster("G:/My Drive/CHELSA_1.2/envlayers_south_aggfac2/landcover_SPHAsouth_aggfac2.asc")
landcovras
plot(landcovras)

maxentSraw <- raster("G:/My Drive/CHELSA_1.2/maxent3.4.1_SPHA_S_rawoutput/SPHA_S.asc")
maxentSrawresample <- resample(maxentSraw, landcovras, method="bilinear")
maxentSrawresample
plot(maxentSrawresample)
maxentSrawresample <- mask(maxentSrawresample, landcovras)
maxentSrawresample <- maxentSrawresample*1000
writeRaster(maxentNrawresample, file="G:/My Drive/CHELSA_1.2/envlayers_north_aggfac2/maxentraw_SPHANorth_aggfac2.asc")

landcovrasN <- raster("G:/My Drive/CHELSA_1.2/envlayers_north_aggfac2/landcover_SPHANorth_aggfac2.asc")

maxentNraw <- raster("G:/My Drive/CHELSA_1.2/maxent3.4.1_SPHA_N_rawoutput/SPHA_N.asc")
maxentNrawresample <- resample(maxentNraw, landcovrasN, method="bilinear")
maxentNrawresample
plot(maxentNrawresample)
maxentNrawresample <- mask(maxentNrawresample, landcovrasN)
maxentNrawresample <- maxentNrawresample*100
maxentNrawresample[maxentNrawresample<1] <- 1

pcarasters <- readRDS("G:/My Drive/CHELSA_1.2/rasters.rga.south.pca.standardized.rds")
plot(pcarasters$map$PC1)
plot(pcarasters$map$PC2)
plot(pcarasters$map$PC3)
plot(pcarasters$map$PC4)
plot(pcarasters$map$PC5)
plot(pcarasters$map$PC6)

pcarasters$model$loadings
write.csv(pcarasters$model$loadings, file="G:/My Drive/CHELSA_1.2/rasterPCAloadings.csv")

writeRaster(pcarasters$map$PC1, file="G:/My Drive/CHELSA_1.2/envlayers_south_aggfac2/PC1ofrastersSPHAsouth_aggfac2.asc", format="ascii")
writeRaster(pcarasters$map$PC2, file="G:/My Drive/CHELSA_1.2/envlayers_south_aggfac2/PC2ofrastersSPHAsouth_aggfac2.asc", format="ascii")
writeRaster(pcarasters$map$PC3, file="G:/My Drive/CHELSA_1.2/envlayers_south_aggfac2/PC3ofrastersSPHAsouth_aggfac2.asc", format="ascii")
writeRaster(pcarasters$map$PC4, file="G:/My Drive/CHELSA_1.2/envlayers_south_aggfac2/PC4ofrastersSPHAsouth_aggfac2.asc", format="ascii")
writeRaster(pcarasters$map$PC5, file="G:/My Drive/CHELSA_1.2/envlayers_south_aggfac2/PC5ofrastersSPHAsouth_aggfac2.asc", format="ascii")


pcarasters2 <- readRDS("G:/My Drive/CHELSA_1.2/rasters.rga.south.pca.rds")
plot(pcarasters2$map$PC1)
plot(pcarasters2$map$PC2)
plot(pcarasters2$map$PC3)
plot(pcarasters2$map$PC4)
plot(pcarasters2$map$PC5)

par(mfrow=c(2,2))

monthcountras <- raster("G:/My Drive/CHELSA_1.2/CHELSA_CA_croppedmasked/monthCountByTemp10_ca_30s.asc")
plot(monthcountras, col=viridis(256))
plot(monthcountras, col=inferno(256))
plot(monthcountras, col=magma(256))
plot(monthcountras, col=plasma(256))
plot(monthcountras, col=cividis(256))

# looking at output of rasterPCA on 30arcsec rasters
pcahiresrasters <- readRDS("G:/My Drive/CHELSA_1.2/rasters.rga.south.hires.pca.standardized.rds")
pcahiresrasters2 <- readRDS("G:/My Drive/CHELSA_1.2/rasters.rga.south.hires.pca.rds")
par(mfrow=c(1,2))
pcahiresstack <- stack(pcahiresrasters$map)
plot(pcahiresrasters$map$PC1)
plot(pcahiresrasters$map$PC2)
plot(pcahiresrasters$map$PC3)
plot(pcahiresrasters$map$PC4)
plot(pcahiresrasters$map$PC5)
plot(pcahiresrasters$map$PC6)

plot(pcahiresstack[[1:12]], col=viridis(256))
plot(pcahiresstack[[6]], col=inferno(256))

pcahiresrasters$model$loadings
write.csv(pcahiresrasters$model$loadings, file="G:/My Drive/CHELSA_1.2/rasterPCAhiresloadings.csv")

writeRaster(pcahiresrasters$map$PC1, file="G:/My Drive/CHELSA_1.2/envlayers_south_croppedmasked/PC1ofrastersSPHAsouth30s.asc", format="ascii")
writeRaster(pcahiresrasters$map$PC2, file="G:/My Drive/CHELSA_1.2/envlayers_south_croppedmasked/PC2ofrastersSPHAsouth30s.asc", format="ascii")
writeRaster(pcahiresrasters$map$PC3, file="G:/My Drive/CHELSA_1.2/envlayers_south_croppedmasked/PC3ofrastersSPHAsouth30s.asc", format="ascii")
writeRaster(pcahiresrasters$map$PC4, file="G:/My Drive/CHELSA_1.2/envlayers_south_croppedmasked/PC4ofrastersSPHAsouth30s.asc", format="ascii")
writeRaster(pcahiresrasters$map$PC5, file="G:/My Drive/CHELSA_1.2/envlayers_south_croppedmasked/PC5ofrastersSPHAsouth30s.asc", format="ascii")


# maxent cloglog * 100
maxentrasS <- raster("G:/My Drive/CHELSA_1.2/envlayers_south_croppedmasked/maxent_SPHAsouth.asc")
maxentrasS100 <- maxentrasS*100
plot(stack(maxentrasS, maxentrasS100))
par(mfrow=c(1,2))
plot(maxentrasS100, col=inferno(256))
plot(maxentrasS100, col=magma(256))
writeRaster(maxentrasS100, file="G:/My Drive/CHELSA_1.2/envlayers_south_croppedmasked/maxentcloglogx100SPHAsouth30s.asc", format="ascii")


# top layers ranked in low-res RGA run for south
toplayerssouth <- stack("G:/My Drive/CHELSA_1.2/envlayers_south_croppedmasked/topoWet_SPHAsouth.asc",
                        "G:/My Drive/CHELSA_1.2/envlayers_south_croppedmasked/toporuggedness_SPHAsouth.asc",
                        "G:/My Drive/CHELSA_1.2/envlayers_south_croppedmasked/slopedeg_SPHAsouth.asc",
                        "G:/My Drive/CHELSA_1.2/envlayers_south_croppedmasked/continentality_SPHAsouth.asc",
                        "G:/My Drive/CHELSA_1.2/envlayers_south_croppedmasked/climaticMoistureIndex_SPHAsouth.asc",
                        "G:/My Drive/CHELSA_1.2/envlayers_south_croppedmasked/landcover_SPHAsouth.asc")
plot(toplayerssouth[[2:3]], col=inferno(256))


rgatable <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth184spatial/All_Results_Table_commuteDistance_50reps_fac2layers.csv",
                     row.names = 1)
rownames(rgatable)
rgatable <- rgatable[-c(49:54),]

par(mfrow=c(1,2))
rgatable.kendall <- corrplot(cor(rgatable, method="kendall"), method="number")
rgatable.spearman <- corrplot(cor(rgatable, method="spearman"), method="number")
par(mfrow=c(1,1))
View(rgatable)


# hires, 100 reps, RGA-optimized CMI
cmi.hires.rgaoptimized <- raster("G:/My Drive/CHELSA_1.2/climaticMoistureIndexSPHAsouth30s.asc")
cmi.rgaoptimized <- raster("G:/My Drive/CHELSA_1.2/climaticMoistureIndexSPHAsouth.asc")
par(mfrow=c(1,2))
plot(cmi.hires.rgaoptimized, col=viridis(256))
plot(toplayerssouth$climaticMoistureIndex_SPHAsouth, col=viridis(256))


# looking at output of rasterPCA on 30arcsec rasters
pcaaggfac2rasters <- readRDS("G:/My Drive/CHELSA_1.2/rasters.rga.north.aggfac2.pca.standardized.rds")
pcaaggfac2rasters2 <- readRDS("G:/My Drive/CHELSA_1.2/rasters.rga.north.aggfac2.pca.rds")
par(mfrow=c(1,1))
pcaaggfac2stack <- stack(pcaaggfac2rasters$map)
plot(pcaaggfac2rasters$map$PC1)
plot(pcaaggfac2rasters$map$PC2)
plot(pcaaggfac2rasters$map$PC3)
plot(pcaaggfac2rasters$map$PC4)
plot(pcaaggfac2rasters$map$PC5)
plot(pcaaggfac2rasters$map$PC6)

plot(pcaaggfac2stack[[1:9]], col=inferno(256))
plot(pcaaggfac2stack[[6]], col=inferno(256))

pcaaggfac2rasters$model$loadings
write.csv(pcaaggfac2rasters$model$loadings, file="G:/My Drive/CHELSA_1.2/rasterPCAnorth_aggfac2loadings.csv")

writeRaster(pcaaggfac2rasters$map$PC1, file="G:/My Drive/CHELSA_1.2/envlayers_north_aggfac2/PC1ofrastersSPHANorth_aggfac2.asc", format="ascii")
writeRaster(pcaaggfac2rasters$map$PC2, file="G:/My Drive/CHELSA_1.2/envlayers_north_aggfac2/PC2ofrastersSPHANorth_aggfac2.asc", format="ascii")
writeRaster(pcaaggfac2rasters$map$PC3, file="G:/My Drive/CHELSA_1.2/envlayers_north_aggfac2/PC3ofrastersSPHANorth_aggfac2.asc", format="ascii")
writeRaster(pcaaggfac2rasters$map$PC4, file="G:/My Drive/CHELSA_1.2/envlayers_north_aggfac2/PC4ofrastersSPHANorth_aggfac2.asc", format="ascii")
writeRaster(pcaaggfac2rasters$map$PC5, file="G:/My Drive/CHELSA_1.2/envlayers_north_aggfac2/PC5ofrastersSPHANorth_aggfac2.asc", format="ascii")
writeRaster(pcaaggfac2rasters$map$PC6, file="G:/My Drive/CHELSA_1.2/envlayers_north_aggfac2/PC6ofrastersSPHANorth_aggfac2.asc", format="ascii")


# maxent cloglog * 100
maxentrasN <- raster("G:/My Drive/CHELSA_1.2/envlayers_north_aggfac2/maxent_SPHANorth_aggfac2.asc")
maxentrasN100 <- maxentrasN*100
plot(stack(maxentrasS, maxentrasS100))
par(mfrow=c(1,2))
plot(maxentrasN100, col=inferno(256))
plot(maxentrasN100, col=magma(256))

writeRaster(maxentrasN100, file="G:/My Drive/CHELSA_1.2/envlayers_north_aggfac2/maxentcloglogx100SPHANorth_aggfac2.asc", format="ascii")


# top layers ranked in low-res RGA run for north
toplayersnorth <- stack("G:/My Drive/CHELSA_1.2/envlayers_north_croppedmasked/aridityIndexThornthwaite_SPHANorth.asc",
                        "G:/My Drive/CHELSA_1.2/envlayers_north_croppedmasked/climaticMoistureIndex_SPHANorth.asc",
                        "G:/My Drive/CHELSA_1.2/envlayers_north_croppedmasked/continentality_SPHANorth.asc",
                        "G:/My Drive/CHELSA_1.2/envlayers_north_croppedmasked/landcover_SPHANorth.asc",
                        "G:/My Drive/CHELSA_1.2/envlayers_north_croppedmasked/toporuggedness_SPHANorth.asc",
                        "G:/My Drive/CHELSA_1.2/envlayers_north_croppedmasked/topoWet_SPHANorth.asc")
plot(toplayersnorth[[2]], col=inferno(256))

CMI.north.rgaoptimized <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/north_ssoptim_pcarasters_aggfac2/climaticMoistureIndexSPHANorth.asc")
plot(CMI.north.rgaoptimized, col=viridis(256))
CMI.north.circuitscape <- raster("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/north_ssoptim_pcarasters_aggfac2/circuitscape_outputs/climaticMoistureIndexRGA_north_cum_curmap.asc")
CMI.north.circuitscape <- mask(CMI.north.circuitscape, CMI.north.rgaoptimized)
plot(CMI.north.circuitscape, breaks=c(1,5,10,12,15,20,25,32,33,34,35), col=inferno(10))
levelplot(CMI.north.circuitscape, col.regions=spectral(255))



maxentsouth2 <- raster("G:/My Drive/CHELSA_1.2/SPHA_S_20190308_all49cropped_redo/SPHA_S.asc")
maxentsouth2 <- mask(crop(maxentsouth2, maxentrasS), maxentrasS)
maxentsouth2x100 <- maxentsouth2*100
plot(maxentsouth2x100, col=viridis(255))
writeRaster(maxentsouth2x100, file="G:/My Drive/CHELSA_1.2/SPHA_S_20190308_all49cropped_redo/maxentcloglogx100SPHAsouth30s.asc")


maxentsouth3 <- raster("G:/My Drive/CHELSA_1.2/SPHA_S_20190308_all49cropped_redo3_nolandcover_cv10/SPHA_S_avg.asc")
maxentsouth3
maxentsouth3 <- mask(crop(maxentsouth3, maxentrasS), maxentrasS)
maxentsouth3x100 <- maxentsouth3*100
par(mfrow=c(1,3))
plot(maxentsouth3x100, col=viridis(255))
writeRaster(maxentsouth3x100, file="G:/My Drive/CHELSA_1.2/SPHA_S_20190308_all49cropped_redo3_nolandcover_cv10/maxentcloglogx100cv10SPHAsouth30s.asc")

maxentsouth4 <- raster("G:/My Drive/CHELSA_1.2/SPHA_S_20190308_all49cropped_redo3_cv10/SPHA_S_avg.asc")
maxentsouth4
maxentsouth4 <- mask(crop(maxentsouth4, maxentrasS), maxentrasS)
maxentsouth4x100 <- maxentsouth4*100
par(mfrow=c(1,2))
plot(maxentsouth4x100, col=viridis(255))
writeRaster(maxentsouth4x100, file="G:/My Drive/CHELSA_1.2/SPHA_S_20190308_all49cropped_redo/maxentcloglogx100cv10SPHAsouth30s.asc")


maxentusedinRGAsouth <- raster("G:/My Drive/CHELSA_1.2/SPHA_S_20190308_all49cropped_redo/maxent_SPHAsouth.asc")
plot(maxentusedinRGAsouth, col=viridis(255))
maxentusedinRGAsouthx100 <- maxentusedinRGAsouth*100
writeRaster(maxentusedinRGAsouthx100, file="G:/My Drive/CHELSA_1.2/SPHA_S_20190308_all49cropped_redo/maxentcloglogx100bioclimandtopographySPHAsouth30s.asc")


maxentNnew <- raster("G:/My Drive/CHELSA_1.2/SPHA_N_20190308_all49cropped_redo3_nolandcover_cv10/SPHA_N_avg.asc")
maxentNnew <- resample(maxentNnew, CMI.north.rgaoptimized, method="bilinear")
maxentNnewmasked <- mask(maxentNnew, CMI.north.rgaoptimized)
maxentNnewx100 <- maxentNnewmasked*100
par(mfrow=c(1,2))
plot(maxentNnew, col=viridis(255))
plot(maxentNnewmasked, col=viridis(255))
writeRaster(maxentNnewmasked, file="G:/My Drive/CHELSA_1.2/SPHA_N_20190308_all49cropped_redo3_nolandcover_cv10/maxentcloglogx100cv10SPHANorth_aggfac2.asc")

maxentNold <- raster("G:/My Drive/CHELSA_1.2/SPHA_N_20190308_all49cropped_redo3_nolandcover_cv10/maxentcloglogx100SPHANorth_aggfac2.asc") # ran this one in rga already
plot(maxentNold, col=viridis(255))

# sphasouth178spatial - RGA-optimized layers
southrgaopt <- stack(list.files(path="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth184spatial/south_ssoptim_pcarasters_hires", pattern="*south30s.asc$", full.names=TRUE))
names(southrgaopt)

plot(southrgaopt, col=viridis(255))

# use impervious surfaces to mask your RGA-optimized rasters, then re-run them in Circuitscape
# The RGA-optimized surfaces presumably are based on historical gene flow; masking impervious
# then reveals how patterns of gene flow are affected in human-impacted landscapes
# can re-run with different %impervious cutoffs
# could then additionally mask out areas based on maxent suitability, being conservative (like cutoff at 0.1 or even lower)
# (this is only really relevant for chapter 2 in Orange County, but would still be a fun exercise)
southurbanmasks <- stack(list.files(path="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth184spatial/south_ssoptim_pcarasters_hires", pattern="^RAW.*asc$", full.names=TRUE))
names(southurbanmasks)
# landcover: 22=low intensity development (20-49% impervious); 23=medium (50-79% impervious), 24=high intensity (80-100% impervious)
southlandcover <- southurbanmasks[[2]]
southlandcover[southlandcover==c(23,24)] <- NA
southimpervious <- southurbanmasks[[1]]
southimpervious[southimpervious>50] <- NA
plot(southlandcover, col=viridis(255))
plot(southimpervious, col=viridis(255))
plot(southurbanmasks[[2]][southurbanmasks[[2]]==24])
plot(southurbanmasks, col=viridis(255))


