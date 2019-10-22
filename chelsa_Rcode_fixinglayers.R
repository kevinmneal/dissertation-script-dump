chelsa.ca <- stack(list.files(pattern=".asc$"))

plot(chelsa.ca[[1]])

#california.extent <- extent(envirem.mask.california) # xmin xmax ymin ymax ; extent for Spea hammondii

chelsa.ca.crop <- mask(crop(chelsa.ca, envirem.mask.california, snap="near"), envirem.mask.california)


chelsa.ccsm <- stack(list.files(path="/home/kevin/Rcoding/CHELSAlayers/", pattern=".tif$", full.names=TRUE))
chelsa.ca.crop <- mask(crop(chelsa.ccsm, asmask, snap="near"), asmask)

setwd("/home/kevin/Rcoding/CHELSAlayers/CHELSAlayers_croptoCA")
for (i in 1:length(names(chelsa.ca.crop))){
  writeRaster(chelsa.ca.crop[[i]], filename=paste(names(chelsa.ca.crop[[i]]), "_ca.asc", sep=""), format="ascii")
  message("done creating raster ", i)
}


bio1lgmtest2 <- (bio1lgmtest/10-273)*10 # need to convert temp layers from kelvin to celsius
# need to convert precip layers from mm to cm? Divide LGM by 10
# bio2: no change needed
# bio3: no change
# bio4: no change
# bio5: temperature convert
# bio6: temp convert
# bio7: no change
# bio8-11: temp convert
# bio12-19: precip convert

chelsa.ca.crop <- mask(crop(chelsa.ccsm, asmask, snap="near"), asmask)
chelsa.ca.crop <- crop(mask(chelsa.ccsm, asmask), asmask, snap="near")
extent(chelsa.ca.crop) <- alignExtent(chelsa.ca.crop, asmask, snap="near")


for (i in 1:length(names(chelsa.ca.crop))){
	if (i %in% c(1,5,6,8,9,10,11)){
		fixedraster <- (chelsa.ca.crop[[i]]/10-273)*10
		NAvalue(fixedraster) <- -9999
		fixedraster <- round(fixedraster, 0)
		writeRaster(fixedraster, filename=paste(names(chelsa.ca.crop[[i]]), "_fixed_ca.asc", sep=""), format="ascii")
		message("Temp transformed; done creating raster ", i)
	} else if (i %in% c(12:19)) {
		fixedraster <- chelsa.ca.crop[[i]]/10
		NAvalue(fixedraster) <- -9999
		fixedraster <- round(fixedraster, 0)
		writeRaster(fixedraster, filename=paste(names(chelsa.ca.crop[[i]]), "_fixed_ca.asc", sep=""), format="ascii")
		message("Precip transformed; done creating raster ", i)
	} else { 
		NAvalue(chelsa.ca.crop[[i]]) <- -9999
		writeRaster(chelsa.ca.crop[[i]], filename=paste(names(chelsa.ca.crop[[i]]), "_fixed_ca.asc", sep=""), format="ascii")
		message("No transformation; done creating raster ", i)
}
}

names(CArasters.uncropped) <- gsub("_ca|_|.asc|30s|CHELSA|taxousda|gmted2010bln|bilinear|nn|masked", "", names(CArasters.uncropped))) # replace underscores cuz they might fuck up functions
names(CArasters.uncropped) <- gsub("_ca_30s|_ca_ca|_|.asc|30s|CHELSA|soilgrids250m|nlcd2011|gmted2010bln|bilinear|_nn|masked", "", names(CArasters.uncropped))
names(CArasters.uncropped[[48]]) <- "toporuggedness"
northrasters.cropped <- mask(crop(CArasters.uncropped, northcrop, snap="near"), northcrop)
southrasters.cropped <- mask(crop(CArasters.uncropped, southcrop, snap="near"), southcrop)

setwd("/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_north_croppedmasked")
for (i in 1:length(names(northrasters.cropped))){
  writeRaster(northrasters.cropped[[i]], filename=paste(names(northrasters.cropped[[i]]), "_SPHANorth.asc", sep=""), format="ascii")
  message("done creating raster ", i, ", ", paste(names(northrasters.cropped[[i]])))
}

setwd("/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_south_croppedmasked")
for (i in 1:length(names(southrasters.cropped))){
  writeRaster(southrasters.cropped[[i]], filename=paste(names(southrasters.cropped[[i]]), "_SPHAsouth.asc", sep=""), format="ascii")
  message("done creating raster ", i, ", ", paste(names(southrasters.cropped[[i]])))
}

#/home/kevin/Rcoding/maxentruns/SPHA_S_avg.asc

northcrop <- raster("/home/kevin/Rcoding/CHELSAlayersCAuncropped/rastersforcropping/sphanorth_depthtobedrock_maskedbetter2_croppedcorrectly.asc")
southcrop <- raster("/home/kevin/Rcoding/CHELSAlayersCAuncropped/rastersforcropping/sphasouth_depthtobedrock_maskedbetter2_croppedcorrectly.asc")

maxent.cloglog.south <- raster("/home/kevin/Rcoding/maxentruns/SPHA_S_avg.asc")
maxent.cloglog.south.cropped <- mask(crop(maxent.cloglog.south, southcrop, snap="near"), southcrop)
writeRaster(maxent.cloglog.south.cropped, filename="/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_south_croppedmasked/maxentSPHAsouth.cloglog", format="ascii")

maxent.cloglog.north <- raster("/home/kevin/Rcoding/maxentruns/SPHA_N_avg.asc")
maxent.cloglog.north.cropped <- mask(crop(maxent.cloglog.north, northcrop, snap="near"), northcrop)
writeRaster(maxent.cloglog.north.cropped, filename="/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_north_croppedmasked/maxentSPHAnorth.cloglog", format="ascii")

corrplot.resis.spearman.north <- corrplot(cor(values(northrasters.cropped), method="spearman", use="pairwise.complete.obs"), main="spearman norcal correlations")
#corrplot.resis.kendall.north <- corrplot(cor(values(northrasters.cropped), method="kendall", use="pairwise.complete.obs"), main="kendall norcal correlations")
#raster.heat <- heatmap(as.matrix(northrasters.cropped))
mutinfo.test.north <- mutinformation(discretize(as.matrix(northrasters.cropped)))
corrplot.resis.spearman.south <- corrplot(cor(values(southrasters.cropped), method="spearman", use="pairwise.complete.obs"), main="spearman socal correlations")
corrplot.resis.kendall.south <- corrplot(cor(values(southrasters.cropped), method="kendall", use="pairwise.complete.obs"), main="kendall socal correlations")
#raster.heat <- heatmap(as.matrix(southrasters.cropped))
mutinfo.test.south <- mutinformation(discretize(as.matrix(southrasters.cropped)))


heatmap(abs(corrplot.resis.spearman), symm=T, main="rga input rasters, abs. value of spearman correlations", col=viridis(256),
        mar=c(15,15))
heatmap.2(mutinfo.test, symm=T, main="rga input rasters, mutual information", col=viridis(256),
          mar=c(15,15), trace="none")
# trying heatmap.2 in gplots package
heatmap.2(corrplot.resis.spearman, symm=T, main="rga input rasters, spearman correlations", col=viridis, trace="none")


#java maxent.jar density.MaxEnt nowarnings noprefixes -E "" -E SPHA_S responsecurves jackknife "outputdirectory=G:\My Drive\CHELSA_1.2\maxent3.4.1_SPHA_S_20190303_withenviremvars" "projectionlayers=G:\My Drive\CHELSA_1.2\CHELSA_CCSM_CA" "samplesfile=G:\My Drive\Illumina Sequencing Data\20181212_rangewide\Maxent_20190302\SPHA_S_maxent_coords_20190302.csv" "environmentallayers=G:\My Drive\CHELSA_1.2\CHELSA_CA_croppedmasked" replicates=5 threads=2


setwd("/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_north_croppedmasked")
setwd("/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_south_croppedmasked")
northrasters.cropped <- stack(list.files(pattern=".asc$"))
setwd("/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_south_croppedmasked")
southrasters.cropped <- stack(list.files(pattern=".asc$"))
names(northrasters.cropped) <- gsub("_SPHANorth", "", names(northrasters.cropped))
names(southrasters.cropped) <- gsub("_SPHASouth", "", names(southrasters.cropped))

# aggregate and crop ; aggregate with "modal" if categorical and "median" if continuous; always ignore na
topotestsouthextent[is.na(topotestsouthextent)] <- -9999
topotestsouth.2.useasmask <- aggregate(topotestsouthextent, fact=2, fun=modal, na.rm=TRUE)

depthtobedrock.aggfac2 <- aggregate(depthtobedrock, fact=2, fun=median, na.rm=TRUE)
depthtobedrock.aggfac2 <- mask(depthtobedrock.aggfac2, topotestsouth.2.useasmask)



resampmaskforRGA <- function(rasterobject, maskobject, outdir){
for (i in 1:length(names(rasterobject))){
  uniquevals <- length(unique(rasterobject[[i]]))
  rasname <- gsub("_SPHANorth|_SPHAsouth|_SPHASouth|_SPHAnorth", "", names(rasterobject[[i]]))
  if (rasname == "taxousda" || rasname == "landcover" || rasname == "bpsprehumanvegmodel" || rasname == "monthCountByTemp"){
    message(paste(names(rasterobject[[i]])), " is categorical; using ngb method")
	ras.aggfac2 <- resample(rasterobject[[i]], maskobject, method="ngb")
	ras.aggfac2 <- mask(ras.aggfac2, maskobject)
  } else { 
    message(paste(names(rasterobject[[i]])), " is continuous; using bilinear method")
	ras.aggfac2 <- resample(rasterobject[[i]], maskobject, method="bilinear")
	ras.aggfac2 <- mask(ras.aggfac2, maskobject)
  }
  writeRaster(ras.aggfac2, filename=paste(outdir,names(rasterobject[[i]]), "_aggfac2.asc", sep=""), format="ascii")
  message("done creating raster ", i, ", ", paste(names(rasterobject[[i]])))
}
}

#/home/kevin/Rcoding/CHELSAlayersCAuncropped/rastersforcropping/aggregatedfac2_SPHANorth_useasmask.asc

northmask.aggfac2 <- raster("/home/kevin/Rcoding/CHELSAlayersCAuncropped/rastersforcropping/aggregatedfac2_SPHANorth_useasmask.asc")
southmask.aggfac2 <- raster("/home/kevin/Rcoding/CHELSAlayersCAuncropped/rastersforcropping/aggregatedfac2_SPHAsouth_useasmask.asc")
northrasters2.cropped <- stack(list.files(path="/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_north_croppedmasked", pattern=".asc$", full.name=TRUE))
southrasters2.cropped <- stack(list.files(path="/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_south_croppedmasked", pattern=".asc$", full.name=TRUE))

northmask.outdir <- "/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_north_aggfac2/"
southmask.outdir <- "/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_south_aggfac2/"

resampmaskforRGA(rasterobject=northrasters2.cropped, maskobject=northmask.aggfac2, outdir=northmask.outdir)
resampmaskforRGA(rasterobject=southrasters2.cropped, maskobject=southmask.aggfac2, outdir=southmask.outdir)

for (i in 1:length(names(rasters.rga.south))){
	message(paste(length(unique(rasters.rga.south[[i]])))," ", names(rasters.rga.south[[i]]))
	}


library(ResistanceGA)
options(bitmapType='cairo')
rasters.rga.south <- stack(list.files(path="/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_south_aggfac2", pattern="aggfac2.asc$", full.name=TRUE))
names(rasters.rga.south) <- gsub("_|aggfac2|", "", names(rasters.rga.south))
rasters.rga.south.continuous <- # dropLayer(rasters.rga.south, c())

# make hires PCAs
rasters.rga.south.hires <- stack(list.files(path="/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_south_croppedmasked", pattern="SPHAsouth.asc$", full.name=TRUE))
rasters.rga.south.hires <- dropLayer(rasters.rga.south.hires, c(15,16,22,31,34,35,38,46)) # drop categorical layers and bio13, 14, monthCountByTemp, degdays0
names(rasters.rga.south.hires) <- gsub("_|aggfac2|", "", names(rasters.rga.south.hires))
mask.hires <- raster("/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_south_croppedmasked/maxent_SPHAsouth.asc")
rasters.rga.south.hires.masked <- mask(rasters.rga.south.hires, mask.hires)
rasters.rga.south.hires.pca.standardized <- rasterPCA(rasters.rga.south.hires.masked, nComp=20, spca=TRUE, maskCheck=FALSE)
rasters.rga.south.hires.pca <- rasterPCA(rasters.rga.south.hires.masked, nComp=20, spca=FALSE, maskCheck=FALSE)
saveRDS(rasters.rga.south.hires.pca.standardized, file="/home/kevin/Rcoding/CHELSAlayersCAuncropped/rasters.rga.south.hires.pca.standardized.rds")
saveRDS(rasters.rga.south.hires.pca, file="/home/kevin/Rcoding/CHELSAlayersCAuncropped/rasters.rga.south.hires.pca.rds")



# drop problematic layers
rasters.rga.south.culled <- dropLayer(rasters.rga.south, c(22,35,36,49))
names(rasters.rga.south.culled)
rasters.rga.south <- dropLayer(rasters.rga.south, c(16,31,39)) # drop bio14, degdays0, and monthCountByTemp

rga.popcoords.south <- read.csv("/home/kevin/Rcoding/ResistanceGAruns/sphasouth178spatial_uniquepopcoords_forrga_fixed.csv")
rga.points.south <- SpatialPoints(rga.popcoords.south[,2:3])
rga.Da.matrix.75pctmsng.south <- readRDS("/home/kevin/Rcoding/ResistanceGAruns/distancemat.Da.75pctmsng.sphasouth178spatial.rds")
results.dir.south <- "/home/kevin/Rcoding/ResistanceGAruns/south_realrun1_ssoptim_allrasters/"

source("/home/kevin/Rcoding/ResistanceGAruns/GA_prep_higherminforcategorical.R") # I changed the threshold for what determines if a layer is categorical to 35
source("/home/kevin/Rcoding/ResistanceGAruns/RGAinternal_helper_functions.R") # I changed the threshold for what determines if a layer is categorical to 35
library(GA)

gdist.inputs.south <- gdist.prep(n.Pops=length(rga.points.south),
                           response=lower(rga.Da.matrix.75pctmsng.south),
                           samples=rga.points.south,
                           longlat=TRUE,
                           method="commuteDistance")

GA.inputs.test.south <- GA.prep(ASCII.dir=rasters.rga.south,
                       Results.dir=results.dir.south,
                       parallel=12, 
                       maxiter=2)
GA.inputs.test.south.onesurface <- GA.prep(ASCII.dir=rasters.rga.south$PC1ofrastersSPHAsouth,
                       Results.dir=results.dir.south,
                       parallel=12, 
                       maxiter=2)
 
GA.inputs.south <- GA.prep(ASCII.dir=rasters.rga.south,
                       Results.dir=results.dir.south,
                       parallel=12, 
                       maxiter=50)

 
GA.inputs.test.south$layer.names
GA.inputs.test.south$surface.type # change any categorical layers to categorical if necessary
# GA.inputs.test.south$surface.type[22] <- "cat"

sink(paste(results.dir.south,"SS_optim_out.txt", sep=""))
SS.results.test.south.culled <- SS_optim(gdist.inputs=gdist.inputs.south,
                                   GA.inputs=GA.inputs.south,
                                   dist_mod=TRUE,
                                   null_mod=TRUE)
sink()

sink(paste(results.dir.south,"GA.inputs.test.southbpsprehumanvegmodelonly_newGAfunction.out.txt", sep=""))
GA.inputs.test.south
sink()

# run ResistanceGA on first 5 or so PCA rasters, + landcover and maxent 
rasterspca.rga.south.hires <- stack(list.files(path="/home/kevin/Rcoding/CHELSAlayersCAuncropped/PCArasters_hires_south", pattern=".asc$", full.name=TRUE))

results.dir.south.rasterspca <- "/home/kevin/Rcoding/ResistanceGAruns/south_ssoptim_pcarasters_hires/"

for (i in 1:length(names(rasterspca.rga.south.hires))){
	message(paste(length(unique(rasterspca.rga.south.hires[[i]])))," ", names(rasterspca.rga.south.hires[[i]]))
	}


GA.inputs.south.rasterspca <- GA.prep(ASCII.dir=rasterspca.rga.south.hires,
                       Results.dir=results.dir.south.rasterspca,
                       parallel=12, 
                       maxiter=100)

sink(paste(results.dir.south.rasterspca,"pcarasters_SS_optim_out.txt", sep=""))
SS.results.south.rasterspca <- SS_optim(gdist.inputs=gdist.inputs.south,
                                   GA.inputs=GA.inputs.south.rasterspca,
                                   dist_mod=TRUE,
                                   null_mod=TRUE)
sink()



# ran this and got an error "Computation failed in `stat_bin()`:`binwidth` must be positive" - related to ggplot? 

GA | iter = 1 | Mean = -32189.321 | Best =   1770.434
GA | iter = 2 | Mean = -14130.767 | Best =   1770.434
GA | iter = 3 | Mean = -5084.106 | Best =  1770.434
GA | iter = 4 | Mean = -9599.574 | Best =  1770.434
GA | iter = 5 | Mean = -20898.394 | Best =   1770.434
GA | iter = 1 | Mean = -23133.22 | Best =   1769.87
GA | iter = 2 | Mean = -32169.15 | Best =   1769.87
GA | iter = 3 | Mean = -23126.19 | Best =   1769.87
GA | iter = 4 | Mean = -18597.932 | Best =   1769.871
GA | iter = 5 | Mean = -9554.681 | Best =  1769.871
Error in :
  task 1 failed - "arguments imply differing number of rows: 35, 3"
In addition: Warning messages:
1: Computation failed in `stat_bin()`:
`binwidth` must be positive
2: Computation failed in `stat_bin()`:
`binwidth` must be positive
3: Computation failed in `stat_bin()`:
`binwidth` must be positive
4: Computation failed in `stat_bin()`:
`binwidth` must be positive

Error in  :
  task 1 failed - "arguments imply differing number of rows: 435, 406"
> traceback()
5: stop(simpleError(msg, call = expr))
4: e$fun(obj, substitute(ex), parent.frame(), e$data)
3: foreach(i. = seq_len(popSize), .combine = "c") %DO% {
       if (is.na(Fitness[i.]))
           do.call(fitness, c(list(Pop[i., ]), callArgs))
       else Fitness[i.]
   }
2: ga(type = "real-valued", fitness = Resistance.Opt_single, Resistance = r,
       population = GA.inputs$population, selection = GA.inputs$selection,
       pcrossover = GA.inputs$pcrossover, pmutation = GA.inputs$pmutation,
       crossover = GA.inputs$crossover, Min.Max = GA.inputs$Min.Max,
       GA.inputs = GA.inputs, gdist.inputs = gdist.inputs, lower = GA.inputs$min.list[[i]],
       upper = GA.inputs$max.list[[i]], parallel = GA.inputs$parallel,
       popSize = GA.inputs$pop.mult * length(GA.inputs$max.list[[i]]),
       maxiter = GA.inputs$maxiter, run = GA.inputs$run, keepBest = GA.inputs$keepBest,
       elitism = GA.inputs$percent.elite, mutation = GA.inputs$mutation,
       seed = GA.inputs$seed, iter = i, quiet = GA.inputs$quiet)
1: SS_optim(gdist.inputs = gdist.inputs.south, GA.inputs = GA.inputs.test.south.continuous,
       dist_mod = TRUE, null_mod = TRUE)

Issue seems to be with the maxentcloglog layer... and toporuggedness, and bpsprehumanvegmodel
bpsprehumanvegmodel only fails if I assign its surface.type to "cat"

SOLUTION: SOME RASTERS HAVE NA VALUES ALONG THE EDGE AND SOME POINTS ARE NOT ON THE RASTER!!!

rasters.rga.north <- stack(list.files(path="/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_north_aggfac2", pattern="aggfac2.asc$", full.name=TRUE))
names(rasters.rga.north) <- gsub("_|aggfac2|", "", names(rasters.rga.north))
rga.popcoords.north <- read.csv("/home/kevin/Rcoding/ResistanceGAruns/sphanorth124spatial_uniquepopcoords_forrga.csv")
rga.points.north <- SpatialPoints(rga.popcoords.north[,2:3])
rga.Da.matrix.75pctmsng.north <- readRDS("/home/kevin/Rcoding/ResistanceGAruns/distancemat.Da.75pctmsng.sphanorth124spatial.rds")
results.dir.north <- "/home/kevin/Rcoding/ResistanceGAruns/north_testrun_ssoptim_allrasters/"

gdist.inputs.north <- gdist.prep(n.Pops=length(rga.points.north),
                           response=lower(rga.Da.matrix.75pctmsng.north),
                           samples=rga.points.north,
                           longlat=TRUE,
                           method="commuteDistance")

GA.inputs.test.north <- GA.prep(ASCII.dir=rasters.rga.north,
                       Results.dir=results.dir.north,
                       parallel=12, 
                       maxiter=5)
  
SS.results.test.north <- SS_optim(gdist.inputs=gdist.inputs.north,
                                   GA.inputs=GA.inputs.test.north,
                                   dist_mod=TRUE,
                                   null_mod=TRUE)  


## NORTH!!!
rasters.rga.north.aggfac2 <- stack(list.files(path="/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_north_aggfac2", pattern="SPHANorth_aggfac2.asc$", full.name=TRUE))
rasters.rga.north.aggfac2 <- dropLayer(rasters.rga.north.aggfac2, c(15,16,22,31,34,35,36,39,47)) # drop categorical layers and bio13, 14, monthCountByTemp, degdays0, and maxent
names(rasters.rga.north.aggfac2) <- gsub("_|aggfac2|", "", names(rasters.rga.north.aggfac2))
mask.north.aggfac2 <- raster("/home/kevin/Rcoding/CHELSAlayersCAuncropped/envlayers_north_aggfac2/maxentcloglog_SPHANorth_aggfac2.asc")
rasters.rga.north.aggfac2.masked <- mask(rasters.rga.north.aggfac2, mask.north.aggfac2)
rasters.rga.north.aggfac2.pca.standardized <- rasterPCA(rasters.rga.north.aggfac2.masked, nComp=20, spca=TRUE, maskCheck=FALSE)
rasters.rga.north.aggfac2.pca <- rasterPCA(rasters.rga.north.aggfac2.masked, nComp=20, spca=FALSE, maskCheck=FALSE)
saveRDS(rasters.rga.north.aggfac2.pca.standardized, file="/home/kevin/Rcoding/CHELSAlayersCAuncropped/rasters.rga.north.aggfac2.pca.standardized.rds")
saveRDS(rasters.rga.north.aggfac2.pca, file="/home/kevin/Rcoding/CHELSAlayersCAuncropped/rasters.rga.north.aggfac2.pca.rds")


# run ResistanceGA on first 5 or so PCA rasters, + landcover and maxent 
rasterspca.rga.north.aggfac2 <- stack(list.files(path="/home/kevin/Rcoding/CHELSAlayersCAuncropped/PCArasters_aggfac2_north", pattern=".asc$", full.name=TRUE))
names(rasterspca.rga.north.aggfac2) <- gsub("_|aggfac2|", "", names(rasterspca.rga.north.aggfac2))
rasterspca.rga.north.aggfac2 <- mask(rasterspca.rga.north.aggfac2, mask.north.aggfac2) # seems doing this turns it into a rasterbrick, but needs to be a rasterstack
rasterspca.rga.north.aggfac2 <- stack(rasterspca.rga.north.aggfac2)
results.dir.north.rasterspca <- "/home/kevin/Rcoding/ResistanceGAruns/north_ssoptim_pcarasters_aggfac2/"

for (i in 1:length(names(rasterspca.rga.north.aggfac2))){
	message(paste(length(unique(rasterspca.rga.north.aggfac2[[i]])))," ", names(rasterspca.rga.north.aggfac2[[i]]))
	}

	
rga.popcoords.north <- read.csv("/home/kevin/Rcoding/ResistanceGAruns/sphanorth124spatial_uniquepopcoords_forrga.csv")
rga.points.north <- SpatialPoints(rga.popcoords.north[,2:3])
rga.Da.matrix.75pctmsng.north <- readRDS("/home/kevin/Rcoding/ResistanceGAruns/distancemat.Da.75pctmsng.sphanorth124spatial.rds")

gdist.inputs.north <- gdist.prep(n.Pops=length(rga.points.north),
                           response=lower(rga.Da.matrix.75pctmsng.north),
                           samples=rga.points.north,
                           longlat=TRUE,
                           method="commuteDistance")
	

GA.inputs.north.rasterspca <- GA.prep(ASCII.dir=rasterspca.rga.north.aggfac2,
                       Results.dir=results.dir.north.rasterspca,
                       parallel=12, 
                       maxiter=100)

sink(paste(results.dir.north.rasterspca,"pcarasters_SS_optim_out.txt", sep=""))
SS.results.north.rasterspca <- SS_optim(gdist.inputs=gdist.inputs.north,
                                   GA.inputs=GA.inputs.north.rasterspca,
                                   dist_mod=TRUE,
                                   null_mod=TRUE)
sink()



####
### do a run just for maxent layer in south 

#maxentsouth <- stack(list.files(path="/home/kevin/Rcoding/CHELSAlayersCAuncropped/south_croppedmasked_maxentonly", pattern=".asc$", full.name=TRUE))
maxentsouth <- raster("/home/kevin/Rcoding/CHELSAlayersCAuncropped/south_croppedmasked_maxentonly/maxentcloglogx100cv10SPHAsouth30s.asc")
results.dir.south.maxent <- "/home/kevin/Rcoding/ResistanceGAruns/south_ssoptim_maxentonly_better_inversemono_200iters/"

GA.inputs.south.maxent <- GA.prep(ASCII.dir=maxentsouth,
                       Results.dir=results.dir.south.maxent,
                       parallel=12, 
                       maxiter=200,
					   run=50,
					   select.trans=list(7)) # inverse mono and reverse mono only. Probably reverse (5). But north was inverse (7)...

sink(paste(results.dir.south.maxent,"maxent_SS_optim_out.txt", sep=""))
SS.results.south.maxent <- SS_optim(gdist.inputs=gdist.inputs.south,
                                   GA.inputs=GA.inputs.south.maxent,
                                   dist_mod=TRUE,
                                   null_mod=TRUE)
sink()

#north_ssoptim_maxentonly_better
#maxentnorth <- stack(list.files(path="/home/kevin/Rcoding/CHELSAlayersCAuncropped/north_croppedmasked_maxentonly", pattern=".asc$", full.name=TRUE))
maxentnorth <- raster("/home/kevin/Rcoding/CHELSAlayersCAuncropped/north_croppedmasked_maxentonly/maxentcloglogx100cv10SPHANorth_aggfac2.asc")
results.dir.north.maxent <- "/home/kevin/Rcoding/ResistanceGAruns/north_ssoptim_maxentcv10_200iters/"

GA.inputs.north.maxent <- GA.prep(ASCII.dir=maxentnorth,
                       Results.dir=results.dir.north.maxent,
                       parallel=12, 
                       maxiter=200,
					   run=50,
					   select.trans=list(7)) # inverse mono and reverse mono only. Probably reverse (5). But north was inverse (7)... =list(c(5,7)). If running two different rasters, must put list(7,7) or list(c(5,7),c(5,7))

sink(paste(results.dir.north.maxent,"maxent_SS_optim_out.txt", sep=""))
SS.results.north.maxent <- SS_optim(gdist.inputs=gdist.inputs.north,
                                   GA.inputs=GA.inputs.north.maxent,
                                   dist_mod=TRUE,
                                   null_mod=TRUE)
sink()
