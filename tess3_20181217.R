### TESS3

library(tess3r)


#### try tess3r, also uses geno
#spha.tess <- tess3(spha.geno.obj, coord=as.matrix(popcoords[,2:3]), K=1:10, ploidy=2, openMP.core.num=4)
# missing data must be encoded as NA. Seems they appear as "9" in here

### load 75pctmsng dataset
setwd("G:/My Drive/Illumina Sequencing Data/20181128_allsocal303indiv_finaldatasets/75pctmsng_nosingletons")

geno.75pctmsng.OCreduced148.path <- vcf2geno(input.file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/OCreduced148_85clust_75pctmsng_oneSNPmac2_20190205.recode.vcf")
geno.75pctmsng.OCreduced148.shortpath <- "OCreduced148.geno"

geno.75pctmsng.OCreduced148 <- read.geno(input.file=geno.75pctmsng.OCreduced148.shortpath)
geno.75pctmsng.OCreduced148[geno.75pctmsng.OCreduced148==9] <- NA


# load longitude and latitude
coordinates <- as.matrix(popcoords.OCreduced148[,7:8])

#### run tess3
tess.75pctmsng.OCreduced148 <- tess3(geno.75pctmsng.OCreduced148, coord=as.matrix(popcoords.OCreduced148[,7:8]), K=1:30, ploidy=2, openMP.core.num=4)


tess3.obj <- tess.75pctmsng.OCreduced148

par(mfrow=c(1,1))
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")


#my.colors <- c("tomato", "darkblue", "olivedrab", "turquoise", "green","purple", "darkred", "gray", "orange", "lightblue", "pink")
my.palette <- CreatePalette(palettesub, 8)
barplot(qmatrix(tess3.obj, K=4), border = NA, space = 0, 
        main = "Ancestry matrix", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette) -> bp
axis(1, at = 1:nrow(qmatrix(tess3.obj, K=3)), labels = popcoords.OCreduced148$Site[bp$order], las = 3, cex.axis = .4) 

#coordinates <- as.matrix(popcoords[,2:3])

# retrieve tess3 Q matrix for K = 5 clusters 
q.matrix <- qmatrix(tess3.obj, K = 3)
# STRUCTURE-like barplot for the Q-matrix 
barplot(q.matrix, border = NA, space = 0, 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "Ancestry matrix") -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 

#my.colors <- c("tomato", "darkblue", "olivedrab", "turquoise", "green","purple", "darkred", "darkgreen", "orange", "blue", "pink")
#my.palette <- CreatePalette(my.colors, 10)

par(mar=c(0,1,2,1))
par(mfrow=c(11,1))
barplot(qmatrix(tess3.obj, K=2), border = NA, space = 0, 
        main = "Ancestry matrix K=2", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette,
        sort.by.Q = FALSE) -> bp
#axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 
axis(1, at = 1:nrow(q.matrix), labels = popcoords.210$Pop, las = 3, cex.axis = .4) 


barplot(qmatrix(tess3.obj, K=3), border = NA, space = 0, 
        main = "Ancestry matrix K=3", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette,
        sort.by.Q = FALSE) -> bp
#axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 
axis(1, at = 1:nrow(q.matrix), labels = popcoords.210$Pop, las = 3, cex.axis = .4) 

barplot(qmatrix(tess3.obj, K=4), border = NA, space = 0, 
        main = "Ancestry matrix K=4", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette,
        sort.by.Q = FALSE) -> bp
#axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 
axis(1, at = 1:nrow(q.matrix), labels = popcoords.210$Pop, las = 3, cex.axis = .4) 

barplot(q.matrix5, border = NA, space = 0, 
        main = "Ancestry matrix K=5", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette,
        sort.by.Q = FALSE) -> bp
#axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 
axis(1, at = 1:nrow(q.matrix), labels = popcoords$Pop, las = 3, cex.axis = .4) 

barplot(q.matrix6, border = NA, space = 0, 
        main = "Ancestry matrix K=6", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette,
        sort.by.Q = FALSE) -> bp
#axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 
#axis(1, at = 1:nrow(q.matrix), labels = popcoords$Pop, las = 3, cex.axis = .4) 

barplot(q.matrix7, border = NA, space = 0, 
        main = "Ancestry matrix K=7", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette,
        sort.by.Q = FALSE) -> bp
#axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 
#axis(1, at = 1:nrow(q.matrix), labels = popcoords$Pop, las = 3, cex.axis = .4) 

barplot(q.matrix8, border = NA, space = 0, 
        main = "Ancestry matrix K=8", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette,
        sort.by.Q = FALSE) -> bp
#axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 
#axis(1, at = 1:nrow(q.matrix), labels = popcoords$Pop, las = 3, cex.axis = .4) 

barplot(q.matrix9, border = NA, space = 0, 
        main = "Ancestry matrix K=9", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette,
        sort.by.Q = FALSE) -> bp
#axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 
#axis(1, at = 1:nrow(q.matrix), labels = popcoords$Pop, las = 3, cex.axis = 1) 

barplot(q.matrix10, border = NA, space = 0, 
        main = "Ancestry matrix K=10", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette,
        sort.by.Q = FALSE) -> bp
#axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 

barplot(q.matrix11, border = NA, space = 0, 
        main = "Ancestry matrix K=11", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette,
        sort.by.Q = FALSE) -> bp
#axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 
#axis(1, at = 1:nrow(q.matrix), labels = popcoords$Pop, las = 3, cex.axis = 1) 

barplot(q.matrix12, border = NA, space = 0, 
        main = "Ancestry matrix K=12", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette,
        sort.by.Q = FALSE) -> bp
#axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 
axis(1, at = 1:nrow(q.matrix), labels = popcoords$Pop, las = 3, cex.axis = 1) 


coordinates <- as.matrix(cbind(popcoords$Longitude, popcoords$Latitude))

#par(mfrow=c(1,1))
par(mfrow=c(4,3))

q.matrix2 <- qmatrix(tess3.obj, K=2)
plot(q.matrix2, coordinates, method = "map.max", interpol = FieldsKrigModel(10), 
     raster.filename="G:/My Drive/Illumina Sequencing Data/20181128_allsocal303indiv_finaldatasets/75pctmsng_nosingletons/depthtobedrockRhorizon.asc",
     main = "Ancestry coefficients, K=2",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = my.palette)

q.matrix3 <- qmatrix(tess3.obj, K=3)
plot(q.matrix3, coordinates, method = "map.max", interpol = FieldsKrigModel(10), 
     raster.filename="G:/My Drive/Illumina Sequencing Data/20181128_allsocal303indiv_finaldatasets/75pctmsng_nosingletons/depthtobedrockRhorizon.asc",
     main = "Ancestry coefficients, K=3",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = my.palette)

q.matrix4 <- qmatrix(tess3.obj, K=4)
plot(q.matrix4, coordinates, method = "map.max", interpol = FieldsKrigModel(10), 
     raster.filename="G:/My Drive/Illumina Sequencing Data/20181128_allsocal303indiv_finaldatasets/75pctmsng_nosingletons/depthtobedrockRhorizon.asc",
     main = "Ancestry coefficients, K=4",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = my.palette)

q.matrix5 <- qmatrix(tess3.obj, K=5)
plot(q.matrix5, coordinates, method = "map.max", interpol = FieldsKrigModel(10), 
     raster.filename="G:/My Drive/Illumina Sequencing Data/20181128_allsocal303indiv_finaldatasets/75pctmsng_nosingletons/depthtobedrockRhorizon.asc",
     main = "Ancestry coefficients, K=5",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = my.palette)

q.matrix6 <- qmatrix(tess3.obj, K=6)
plot(q.matrix6, coordinates, method = "map.max", interpol = FieldsKrigModel(10), 
     raster.filename="G:/My Drive/Illumina Sequencing Data/20181128_allsocal303indiv_finaldatasets/75pctmsng_nosingletons/depthtobedrockRhorizon.asc",
     main = "Ancestry coefficients, K=6",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = my.palette)

q.matrix7 <- qmatrix(tess3.obj, K=7)
plot(q.matrix7, coordinates, method = "map.max", interpol = FieldsKrigModel(10), 
     raster.filename="elev_socal.asc",
     main = "Ancestry coefficients, K=7",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = my.palette)

q.matrix8 <- qmatrix(tess3.obj, K=8)
plot(q.matrix8, coordinates, method = "map.max", interpol = FieldsKrigModel(10), 
     raster.filename="depthtobedrockRhorizon.asc",
     main = "Ancestry coefficients, K=8",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = my.palette)

q.matrix9 <- qmatrix(tess3.obj, K=9)
plot(q.matrix9, coordinates, method = "map.max", interpol = FieldsKrigModel(10), 
     raster.filename="elev_socal.asc",
     main = "Ancestry coefficients, K=9",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = my.palette)
q.matrix10 <- qmatrix(tess3.obj, K=10)
plot(q.matrix10, coordinates, method = "map.max", interpol = FieldsKrigModel(10), 
     raster.filename="elev_socal.asc",
     main = "Ancestry coefficients, K=10",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = my.palette)

q.matrix11 <- qmatrix(tess3.obj, K=11)
plot(q.matrix11, coordinates, method = "map.max", interpol = FieldsKrigModel(10), 
     raster.filename="elev_socal.asc",
     main = "Ancestry coefficients, K=11",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = my.palette)

q.matrix12 <- qmatrix(tess3.obj, K=12)
plot(q.matrix11, coordinates, method = "map.max", interpol = FieldsKrigModel(10), 
     raster.filename="elev_socal.asc",
     main = "Ancestry coefficients, K=12",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = my.palette)

q.matrix13 <- qmatrix(tess3.obj, K=13)
plot(q.matrix11, coordinates, method = "map.max", interpol = FieldsKrigModel(10), 
     raster.filename="elev_socal.asc",
     main = "Ancestry coefficients, K=13",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = my.palette)




####### tess3, 75pctmsng, on 198 individuals (native socal--IM and SHOE removed)
setwd("G:/My Drive/Illumina Sequencing Data/20181004_ALLsocal_303indiv_highqual/nopopfile_75pctmsng")
geno.nativesocal.75pctmsng.path <- vcf2geno(input.file="vcftools_nativesocal_198indv_highqual_75pctmsng_oneSNP_nosingletons_20181009.recode.vcf")
geno.nativesocal.75pctmsng.shortpath <- "nativesocal_198indv_75pctmsng.recode.geno"

geno.nativesocal.75pctmsng <- read.geno(input.file=geno.nativesocal.75pctmsng.shortpath)
geno.nativesocal.75pctmsng[geno.nativesocal.75pctmsng==9] <- NA

popcoords.198 <- read.csv("nativesocal_198indiv_popcoords_multiplelevels_betternames.csv", 
                          header=T,
                          stringsAsFactors = F)
coordinates.198 <- as.matrix(popcoords.198[,3:4])

  
tess.nativesocal.75pctmsng <- tess3(geno.nativesocal.75pctmsng, coord=as.matrix(popcoords.198[,3,4]), K=1:30, ploidy=2, openMP.core.num=4)

tess3.obj <- tess.nativesocal.75pctmsng

par(mfrow=c(1,1))
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

par(mfrow=c(1,1))
my.palette <- CreatePalette(palettesub, 5)
barplot(qmatrix(tess3.obj, K=9), border = NA, space = 0, 
        main = "Ancestry matrix", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette) -> bp
axis(1, at = 1:nrow(qmatrix(tess3.obj, K=9)), labels = popcoords.198$PopLevel2_prox[bp$order], las = 3, cex.axis = .4) 

par(mfrow = c(8,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(4,0,1,1) + 0.1) # the 4 puts some space on the bottom to allow room for all the axes
par(mfrow=c(3,3))
for (i in 2:10){
  plot(qmatrix(tess3.obj, K=i), coordinates.198, method = "map.max", interpol = FieldsKrigModel(10), 
       raster.filename="elev_socal.asc",
       main = paste0("Ancestry coefficients of natural populations, up to 20% missing data, K=", i),
       xlab = "Longitude", ylab = "Latitude", 
       resolution = c(300,300), cex = .4, 
       col.palette = my.palette) -> bp
  axis(1, at = 1:nrow(qmatrix(tess3.obj, K=i)), labels = popcoords.198$PopLevel2_prox[bp$order], las = 3, cex.axis = .4) 
  
}

plot(qmatrix(tess3.obj, K=9), coordinates, method = "map.max", interpol = FieldsKrigModel(10), 
     raster.filename="elev_socal.asc",
     main = "Ancestry coefficients, K=9",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = my.palette)

str(qmatrix(tess3.obj, K=5))
str(as.matrix(unname(Q(project.75pctmsng, K=5, run=which.min(cross.entropy(project.75pctmsng, K=5))))))

#### 20 percent missing on OC only

####### tess3, 75pctmsng, on 198 individuals (native socal--IM and SHOE removed)
setwd("G:/My Drive/Illumina Sequencing Data/20181004_ALLsocal_303indiv_highqual/popfile_allbignativepops")
geno.nativesocal.75pctmsng.path <- vcf2geno(input.file="vcftools_nativesocal_198indv_highqual_75pctmsng_SNPsinallbignativepops_oneSNP_nosingletons_20181017.recode.vcf")
geno.nativesocal.75pctmsng.shortpath <- "nativesocal_198indv_75pctmsng.recode.geno"

geno.nativesocal.75pctmsng <- read.geno(input.file=geno.nativesocal.75pctmsng.shortpath)
geno.nativesocal.75pctmsng[geno.nativesocal.75pctmsng==9] <- NA

popcoords.198 <- read.csv("nativesocal_198indiv_popcoords_multiplelevels_betternames.csv", 
                          header=T,
                          stringsAsFactors = F)
coordinates.198 <- as.matrix(popcoords.198[,3:4])


tess.nativesocal.75pctmsng <- tess3(geno.nativesocal.75pctmsng, coord=as.matrix(popcoords.198[,3,4]), K=1:30, ploidy=2, openMP.core.num=4)

tess3.obj <- tess.nativesocal.75pctmsng

par(mfrow=c(1,1))
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

par(mfrow=c(1,1))
my.palette <- CreatePalette(palettesub, 5)
barplot(qmatrix(tess3.obj, K=9), border = NA, space = 0, 
        main = "Ancestry matrix", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette) -> bp
axis(1, at = 1:nrow(qmatrix(tess3.obj, K=9)), labels = popcoords.198$PopLevel2_prox[bp$order], las = 3, cex.axis = .4) 

par(mfrow = c(8,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(4,0,1,1) + 0.1) # the 4 puts some space on the bottom to allow room for all the axes
par(mfrow=c(3,3))
for (i in 2:10){
  plot(qmatrix(tess3.obj, K=i), coordinates.198, method = "map.max", interpol = FieldsKrigModel(10), 
       raster.filename="elev_socal.asc",
       main = paste0("Ancestry coefficients of natural populations, up to 20% missing data, K=", i),
       xlab = "Longitude", ylab = "Latitude", 
       resolution = c(300,300), cex = .4, 
       col.palette = my.palette) -> bp
  axis(1, at = 1:nrow(qmatrix(tess3.obj, K=i)), labels = popcoords.198$PopLevel2_prox[bp$order], las = 3, cex.axis = .4) 
  
}

####### tess3, 75pctmsng, on 181 individuals (native OC--IM and SHOE removed)
setwd("G:/My Drive/Illumina Sequencing Data/20181004_ALLsocal_303indiv_highqual/popfile_allbignativepops")
geno.nativeOC.75pctmsng.path <- vcf2geno(input.file="vcftools_nativeOC_181indv_highqual_75pctmsng_SNPsinallbignativepops_oneSNP_nosingletons_20181017.recode.vcf")
geno.nativeOC.75pctmsng.shortpath <- "nativeOC_181indv_75pctmsng_SNPsinallbignativepops.recode.geno"

geno.nativeOC.75pctmsng <- read.geno(input.file=geno.nativeOC.75pctmsng.shortpath)
geno.nativeOC.75pctmsng[geno.nativeOC.75pctmsng==9] <- NA

popcoords.nativeOC.181 <- read.csv("nativeOC_181indiv_popcoords_multiplelevels_betternames.csv", 
                          header=T,
                          stringsAsFactors = F)
coordinates.nativeOC.181 <- as.matrix(popcoords.nativeOC.181[,3:4])


tess.nativeOC.75pctmsng <- tess3(geno.nativeOC.75pctmsng, coord=as.matrix(popcoords.nativeOC.181[,3,4]), K=1:30, ploidy=2, openMP.core.num=4)

tess3.obj <- tess.nativeOC.75pctmsng

par(mfrow=c(1,1))
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

par(mfrow=c(1,1))
my.palette <- CreatePalette(palettesub, 20) # use when > 19 colors needed
shiny.palette <- CreatePalette(palettelist$shiny, 50)
barplot(qmatrix(tess3.obj, K=9), border = NA, space = 0, 
        main = "Ancestry matrix", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette) -> bp
axis(1, at = 1:nrow(qmatrix(tess3.obj, K=9)), labels = popcoords.nativeOC.181$PopLevel2_prox[bp$order], las = 3, cex.axis = .4) 

par(mfrow = c(8,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(4,0,1,1) + 0.1) # the 4 puts some space on the bottom to allow room for all the axes
par(mfrow=c(3,3))

uniquepopnames.nativeOC.181 <- unique(popcoords.nativeOC.181[,c("Longitude",
                                                                "Latitude",
                                                                "PopLevel2_prox")])
for (i in 2:10){
  plot(qmatrix(tess3.obj, K=i), coordinates.nativeOC.181, method = "map.max", interpol = FieldsKrigModel(10), 
       raster.filename="elev_oc.asc",
       main = paste0("Ancestry coefficients, K=", i),
       xlab = "Longitude", ylab = "Latitude", 
       resolution = c(300,300), cex = .4, 
       col.palette = shiny.palette) -> bp
  text(x = uniquepopnames.nativeOC.181$Longitude, 
       y = uniquepopnames.nativeOC.181$Latitude-0.003, 
       labels = uniquepopnames.nativeOC.181$PopLevel2_prox,
       cex=0.7)
}
#plot(raster("elev_oc.asc"), col=viridis(256), main="Elevation", xlab="Longitude", ylab="Latitude")

# save 3x3 as a pdf, 20x20 inches
# save 4x2 as pdf, 25x15 inches
# tess3maps_K2-9_75pctmsng_SNPsinallbignativepops.pdf


par(mfrow=c(1,1))
plot(raster("elev_oc.asc"), col=viridis(256), main="Elevation", xlab="Longitude", ylab="Latitude",
     xlim=c())

par(mfrow=c(3,2))
plot(qmatrix(tess3.obj, K=4), coordinates.nativeOC.181, method = "map.max", interpol = FieldsKrigModel(10), 
     raster.filename="elev_oc.asc",
     main = paste0("Ancestry coefficients (map.max), K=", 4),
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = shiny.palette) -> bp

par(mfrow=c(9,1))
shiny.palette <- CreatePalette(palettelist$shiny, 5) # change to 5, otherwise colors will be too dim
for (i in 2:10){
  barplot(qmatrix(tess3.obj, K=i), border = NA, space = 0, 
          #main = paste0("Ancestry coefficients, K=", i), 
          xlab = "Individuals", ylab = paste0("TESS3 ancestry, K=", i), 
          col.palette = shiny.palette,
          sort.by.Q=FALSE) -> bp
  axis(1, at = 1:nrow(qmatrix(tess3.obj, K=i)), labels = popcoords.nativeOC.181$PopLevel2_prox[bp$order], las = 3, cex.axis = .4) 
  
}

barplot(qmatrix(tess3.obj, K=9), border = NA, space = 0, 
        main = "Ancestry matrix", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette) -> bp
axis(1, at = 1:nrow(qmatrix(tess3.obj, K=9)), labels = popcoords.198$PopLevel2_prox[bp$order], las = 3, cex.axis = .4) 

