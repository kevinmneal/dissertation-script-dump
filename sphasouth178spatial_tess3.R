library(tess3r)
library(pals)
library(viridis)

####### tess3 on sphasouth178spatial
#setwd("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial_radiator/tess3")
geno.sphasouth178spatialpath <- vcf2geno(input.file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/radiator_final/sphasouth178spatial.radiator.50pctmsng.oneSNPmac3.vcf")
#geno.sphasouth178spatial.shortpath <- "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/sphasouth178spatial__oneSNPmac2_20190302.recode.geno"
#sometimes if the filename is too long it won't be able to read the geno, so have to manually shorten it
geno.sphasouth178spatial <- read.geno(input.file=geno.sphasouth178spatialpath)
geno.sphasouth178spatial[geno.sphasouth178spatial==9] <- NA

popcoords.sphasouth178spatial <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial_popcoords.csv",
                                          header=T,
                                          stringsAsFactors = F)

colnames(popcoords.sphasouth178spatial)
coordinates.sphasouth178spatial <- as.matrix(popcoords.sphasouth178spatial[,c("LongitudeRaw","LatitudeRaw")])

Sys.time()
#tess.sphasouth178spatial.old <- tess.sphasouth178spatial
tess.sphasouth178spatial <- tess3(geno.sphasouth178spatial, coord=coordinates.sphasouth178spatial, K=1:10, rep=10, ploidy=2, tolerance=1e-06, mask=0.1, openMP.core.num=4, verbose=TRUE)
Sys.time()
tess3.obj <- tess.sphasouth178spatial

par(mfrow=c(1,1))
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")
par(mfrow=c(1,1))

par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(4,0,1,1) + 0.1) # the 4 puts some space on the bottom to allow room for all the axes
par(mfrow=c(3,2))

uniquepopnames.sphasouth178spatialbettermapping <- unique(popcoords.sphasouth178spatial[,c("LongitudeRaw",
                                                                              "LatitudeRaw",
                                                                              "poplevel2")])
row.names(uniquepopnames.sphasouth178spatialbettermapping) <- NULL
#uniquepopnames.sphasouth178spatialbettermapping <- uniquepopnames.sphasouth178spatial[-c(8,10,20),] # remove coastal popnames since they all overlap
#row.names(uniquepopnames.sphasouth178spatialbettermapping) <- NULL
#uniquepopnames.sphasouth178spatialbettermapping
#uniquepopnames.sphasouth178spatialbettermapping[,3] <- as.character(uniquepopnames.sphasouth178spatialbettermapping[,3])
#uniquepopnames.sphasouth178spatialbettermapping[8,3] <- "Coastsites" # change one of the coastal cluster names to "Coast"

png(file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/TESS3maps_sphasouth178spatial_K3_map.all_cols25cols.png",
          width=10,
          height=10,
          units="in",
    res=320,
    type="cairo-png")
cairo_pdf(file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/TESS3maps_sphasouth178spatial_K3_map.all_cols25cols.pdf",
          width=20,
          height=20,
          onefile=TRUE)
par(mfrow=c(2,2))
for (i in 3){
  plot(qmatrix(tess3.obj, K=i, rep="best"), coordinates.sphasouth178spatial, method = "map.all", interpol = FieldsKrigModel(10), 
       raster.filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/cmiRGA_SPHAsouth_noimperviousmask_maxentmask10.asc",
       main = paste0("SPHA-SOUTH TESS3 ancestry coefficients, K=", i),
       xlab = "Longitude", ylab = "Latitude", 
       resolution = c(300,300), cex = 0.4, 
       col.palette = CreatePalette(cols25(3)[c(3,2,1)], 20)) -> bp
  plot(qmatrix(tess3.obj, K=i, rep="best"), coordinates.sphasouth178spatial, method = "map.max", interpol = FieldsKrigModel(10), 
       raster.filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/cmiRGA_SPHAsouth_noimperviousmask_maxentmask10.asc",
       main = paste0("SPHA-SOUTH TESS3 ancestry coefficients (max), K=", i),
       xlab = "Longitude", ylab = "Latitude", 
       resolution = c(300,300), cex = 0.4, 
       col.palette = CreatePalette(cols25(3)[c(3,2,1)], 20)) -> bp
  #text(x = uniquepopnames.sphasouth178spatialbettermapping$LongitudeRaw+0.012, 
  #     y = uniquepopnames.sphasouth178spatialbettermapping$LatitudeRaw-0.01, 
  #     labels = uniquepopnames.sphasouth178spatialbettermapping$poplevel2,
  #     cex=0.9)
}
dev.off()


# try ggtess3Q
ggtess3Q(qmatrix(tess3.obj, K=4, rep="best"),
         coord=coordinates.sphasouth178spatial,
         raster.filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/north_south_circuitscapedistances/cmiRGA_SPHAsouth_noimperviousmask_maxentmask10.asc",
         col.palette=CreatePalette(cols25(10), 20))




cairo_pdf(file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/tess3/TESS3maps_bestKs_sphasouth178spatial_K2-10_cols25cols2.pdf",
    width=20,
    height=30,
    onefile=TRUE)
par(mfrow=c(3,2))
for (i in c(2,3,4,5,6,10)){
  plot(qmatrix(tess3.obj, K=i, rep="best"), coordinates.sphasouth178spatial, method = "map.max", interpol = FieldsKrigModel(10), 
       raster.filename="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/topoWetSPHAsouth30s.asc",
       main = paste0("Ancestry coefficients, K=", i),
       xlab = "Longitude", ylab = "Latitude", 
       resolution = c(300,300), cex = 0.4, 
       col.palette = CreatePalette(cols25(10), 20)) -> bp
  text(x = uniquepopnames.sphasouth178spatialbettermapping$LongitudeRaw+0.012, 
       y = uniquepopnames.sphasouth178spatialbettermapping$LatitudeRaw-0.01, 
       labels = uniquepopnames.sphasouth178spatialbettermapping$poplevel2,
       cex=0.9)
}
dev.off()





pal.bands(parula, viridis, inferno, magma, cubicl, tol, kovesi.rainbow, ocean.phase,
          cols25, tableau20, polychrome, watlington, n=10, show.names=FALSE)


pdf(file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/tess3/TESS3barplots_sphasouth178spatial_K2-20_sortbyQcols25cols.pdf",
          width=20,
          height=70,
          onefile=TRUE)
par(mfrow=c(19,1))
for (i in 2:20){
  barplot(qmatrix(tess3.obj, K=i), border=NA, space = 0, 
          #main = paste0("TESS3 ancestry coefficients, K=", i), 
          #xlab = "Individuals", 
          ylab = paste0("TESS3 ancestry (sorted by Q), K=", i), 
          col.palette = CreatePalette(cols25(10), 5),
          sort.by.Q=TRUE) -> bp
  axis(1, at = 1:nrow(qmatrix(tess3.obj, K=i)), labels = popcoords.sphasouth178spatial$poplevel2[bp$order], las = 3, cex.axis = 1) 
  
}
dev.off()


cairo_pdf(file="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/TESS3barplots_sphasouth178spatial_K2-10_cols25cols.pdf",
          width=20,
          height=70,
          onefile=TRUE)
par(mfrow=c(19,1))
for (i in 2:20){
  barplot(qmatrix(tess3.obj, K=i), border=NA, space = 0, 
          #main = paste0("TESS3 ancestry coefficients, K=", i), 
          #xlab = "Individuals", 
          ylab = paste0("TESS3 ancestry, K=", i), 
          col.palette = CreatePalette(cols25(10), 5),
          sort.by.Q=FALSE) -> bp
  axis(1, at = 1:nrow(qmatrix(tess3.obj, K=i)), labels = popcoords.sphasouth178spatial$poplevel2[bp$order], las = 3, cex.axis = 1) 
  
}
dev.off()

# outlier

tess3outliers <- pvalue(tess3.obj, K=5)
hist(tess3outliers)

# benjamini-yekutieli - use this one. Has been cited in other papers using pcadapt
padj.tess3.by <- p.adjust(tess3outliers,method="BY")
alpha <- 0.05 # use 0.01 if being conservative?
outliers.tess3.by.pt01 <- which(padj.tess3.by < 0.01)
outliers.tess3.by.pt05 <- which(padj.tess3.by < 0.05)
outliers.tess3.by.pt1 <- which(padj.tess3.by < 0.10)
length(outliers.tess3.by.pt05)
plot(padj.tess3.by)

intersect(outliers.tess3.by.pt1, outliers.by.pt05)
length(intersect(outliers.tess3.by.pt1, outliers.by.pt01))

