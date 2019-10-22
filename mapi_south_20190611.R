library(mapi)

mapi.popcoords.south <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/sphasouth178spatial_uniquepopcoords_forrga.csv")
mapi.popcoords.south <- mapi.popcoords.south
Da.mat.south <- distancemat.Da.sphasouth178spatial.radiator.50pctmsng
diag(Da.mat.south) <- NA
Da.dist.south <- melt(Da.mat.south, varnames = c("ind1", "ind2"))
Da.dist.south <- Da.dist.south[complete.cases(Da.dist.south),]
#write.csv(Da.dist, file="neiDa_nativeOC_20181105.csv")
#Da.dist.interpop <- read.csv("neiDa_nativeOC_20181105_coastalpopscombined.csv", header=T, 
#                                       stringsAsFactors = F)

coords.points.south <- SpatialPoints(mapi.popcoords.south[,2:3], proj4string = CRS('+proj=longlat +datum=WGS84'))
utm11.south <- CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") # ESPG/CRS is 26911
coords.utm.south <- spTransform(coords.points.south, utm11.south)
mapi.utm.south <- data.frame(coords.utm.south)
mapi.utm.south <- cbind(mapi.popcoords.south$SPHASpoplvl2, mapi.utm.south)
colnames(mapi.utm.south) <- c("ind", "x", "y") # must use these exact labels


mapi.results <- MAPI_RunAuto(samples=mapi.utm.south, metric=Da.dist.south, crs=26911, beta=0.25, ecc=0.975, nbPermuts = 1000)
my.tails <- MAPI_Tails(mapi.results, alpha=0.05)
tails <- my.tails
#plot(my.tails)
uTail <- sf::as_Spatial(my.tails[my.tails$tail == "upper", 
                                 ])

lTail <- sf::as_Spatial(my.tails[my.tails$tail == "lower", 
                                 ])

MAPI_Plot(mapi.results, pal=magma(100), tails=tails, lower=TRUE, upper=TRUE)
plot(mapi.results)
