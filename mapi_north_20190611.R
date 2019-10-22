### testing out MAPI - mapping average pairwise distances
### meant to compete with e.g. EEMS, SpaceMix, un-PC

library(mapi)
setwd("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/MAPI_sphanorth/")

mapi.popcoords <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphanorth124spatial/sphanorth124spatial_uniquepopcoords_forrga.csv")
Da.mat <- distancemat.Da.sphanorth124spatial.radiator.50pctmsng
diag(Da.mat) <- NA
Da.dist <- melt(Da.mat, varnames = c("ind1", "ind2"))
Da.dist <- Da.dist[complete.cases(Da.dist),]
#write.csv(Da.dist, file="neiDa_nativeOC_20181105.csv")
#Da.dist.interpop <- read.csv("neiDa_nativeOC_20181105_coastalpopscombined.csv", header=T, 
#                                       stringsAsFactors = F)

FST.mat <- distancemat.wcfst.75pctmsng.sphanorth124spatial
diag(FST.mat) <- NA
FST.dist <- melt(FST.mat, varnames = c("ind1", "ind2"))
FST.dist <- FST.dist[complete.cases(FST.dist),]

JostD.mat <- distancemat.JostD.75pctmsng.sphanorth124spatial
diag(JostD.mat) <- NA
JostD.dist <- melt(JostD.mat, varnames = c("ind1", "ind2"))
JostD.dist <- JostD.dist[complete.cases(JostD.dist),]

coords.points <- SpatialPoints(mapi.popcoords[,2:3], proj4string = CRS('+proj=longlat +datum=WGS84'))
utm11 <- CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") # ESPG/CRS is 26911
coords.utm <- spTransform(coords.points, utm11)
mapi.utm <- data.frame(coords.utm)
mapi.utm <- cbind(mapi.popcoords$poplevel2, mapi.utm)
colnames(mapi.utm) <- c("ind", "x", "y") # must use these exact labels

mapi.results <- MAPI_RunAuto(samples=mapi.utm, metric=Da.dist, crs=26911, beta=0.25, ecc=0.975, nbPermuts = 1000)
my.tails <- MAPI_Tails(mapi.results, alpha=0.05, minQ=5)
tails <- my.tails
#plot(my.tails)
uTail <- sf::as_Spatial(my.tails[my.tails$tail == "upper", 
                              ])

lTail <- sf::as_Spatial(my.tails[my.tails$tail == "lower", 
                              ])

MAPI_Plot(mapi.results, pal=magma(100), tails=tails, lower=TRUE, upper=TRUE)
plot(mapi.results)

mapi.results.reynolds <- MAPI_RunAuto(samples=mapi.utm, metric=, crs=26911, nbPermuts = 1000)
my.tails.reynolds <- MAPI_Tails(mapi.results.reynolds, alpha=0.05, minQ=5)
MAPI_Plot(mapi.results.reynolds, pal=magma(256)) #, tails=my.tails.reynolds, pal=viridis(256))


st_write(mapi.results, dsn=".", layer="MAPI.Da.75pctmsng.sphanorth124spatial_20190307",
         driver="ESRI Shapefile", update=TRUE, delete_layer=TRUE)




## south

setwd("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth184spatial/MAPI_sphasouth/")

mapi.popcoords.south <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth184spatial/sphasouth178spatial_uniquepopcoords_forrga.csv")
Da.mat.south <- dist distancemat.Da.75pctmsng.sphasouth178spatial
diag(Da.mat.south) <- NA
Da.dist.south <- melt(Da.mat.south, varnames = c("ind1", "ind2"))
Da.dist.south <- Da.dist.south[complete.cases(Da.dist.south),]
#write.csv(Da.dist, file="neiDa_nativeOC_20181105.csv")
#Da.dist.interpop <- read.csv("neiDa_nativeOC_20181105_coastalpopscombined.csv", header=T, 
#                                       stringsAsFactors = F)

FST.mat.south <- distancemat.wcfst.75pctmsng.sphasouth
diag(FST.mat.south) <- NA
FST.dist.south <- melt(FST.mat, varnames = c("ind1", "ind2"))
FST.dist <- FST.dist[complete.cases(FST.dist),]

JostD.mat <- distancemat.JostD.75pctmsng.sphanorth124spatial
diag(JostD.mat) <- NA
JostD.dist <- melt(JostD.mat, varnames = c("ind1", "ind2"))
JostD.dist <- JostD.dist[complete.cases(JostD.dist),]

coords.points.south <- SpatialPoints(mapi.popcoords.south[,2:3], proj4string = CRS('+proj=longlat +datum=WGS84'))
utm11.south <- CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") # ESPG/CRS is 26911
coords.utm.south <- spTransform(coords.points.south, utm11.south)
mapi.utm.south <- data.frame(coords.utm.south)
mapi.utm.south <- cbind(mapi.popcoords.south$SPHASpoplvl2, mapi.utm.south)
colnames(mapi.utm.south) <- c("ind", "x", "y") # must use these exact labels

mapi.results.south <- MAPI_RunAuto(samples=mapi.utm.south, metric=Da.dist.south, crs=26911, beta=0.25) #nbPermuts = 1000)
my.tails.south <- MAPI_Tails(mapi.results.south, alpha=0.05, minQ=5)
MAPI_Plot(mapi.results.south, pal=inferno(256))#, tails=my.tails, upper=TRUE, lower=TRUE)
plot(mapi.results.south)

st_write(mapi.results.south$avg_value, dsn=".", layer="MAPI.Da.75pctmsng.sphasouth178spatial_20190307",
         driver="ESRI Shapefile", update=TRUE, delete_layer=TRUE)

