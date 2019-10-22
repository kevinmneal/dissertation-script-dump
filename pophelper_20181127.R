#install.packages(c("Cairo","devtools","ggplot2","gridExtra","gtable","tidyr"),dependencies=T)

#devtools::install_github('royfrancis/pophelper')
library(pophelper)
library(gridExtra)

setwd("G:/My Drive/Illumina Sequencing Data/20181128_allsocal303indiv_finaldatasets")

popcoords.303 <- read.csv("allsocal_303indiv_popcoords_multilevelpops_betternames_20181017.csv", header=T,
                            stringsAsFactors = F)
popcoords <- popcoords.303

# use FASTSTRUCTURE files (do not use this command to read local files)
#ffiles <- list.files(".", pattern="*225indiv_out.\\d*.meanQ", full.names=T) #match all 2-digit numbers
setwd("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/faststructure/k1-26_10reps")
ffiles <- list.files(".", pattern=".meanQ", full.names=T)
ffiles <- ffiles[ffiles=="./OCreduced148_85clust_75pctmsng_oneSNPmac2_20190205_r3_out.2.meanQ"]
flist <- readQ(files=ffiles)

flist.tabq <- tabulateQ(flist)
head(flist.tabq)

flist.sum <- summariseQ(flist.tabq)
flist.sum
# Levels: COAST="#440154FF" INLAND="#297B8EFF" SANTENA="#89D548FF"

poplabels <- data.frame(popcoords.OCreduced148$Site, stringsAsFactors = F)
poplabels[] <- lapply(poplabels, as.character)
colnames(poplabels) <- "Site"
#poplabels <- popcoords.OCreduced148[,3, drop=FALSE]
plotQ148 <- plotQ(flist, #imgoutput = "join", 
      grplab=poplabels,
      grplabpos=0.6,
      returnplot=T,
      exportplot=F,
      basesize=12,
      grplabsize=3,
      linesize=0.8,
      pointsize=4,
      ordergrp=TRUE,
      #sortind="all",
      #sharedindlab=FALSE,
      #sortind="Cluster1",
      indlabangle = 90,
      grplabangle=45,
      selgrp="Site",
      clustercol=c("#297B8EFF", "#440154FF", "#89D548FF"),
      splab=c("Orange County (1), K=2"),
      subset=c("BBEND", "CCSP1", "CCSP2", "CCSP3", "LAGUN", "MORO", "TENA", "SANCAN", 
              "FREM", "GREAT", "IM01", "IM02", "IM03", "IM06", "IM07", "IM09", "IM12", "IM14", "SHOE",
            "LAUR", "LIME", "LOMA", "SADDLE", "STARR", "THOM", "TORO"))
grid.arrange(plotQ148$plot[[1]])

### coast+santen, K=2, highest likelihood rep
"G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/faststructure/OCcoast_75pctmsng/OCcoast_85clust_75pctmsng_oneSNPmac2_20190205_r3_out.2.meanQ"
setwd("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/faststructure/OCcoast_75pctmsng")
ffilescoast <- list.files(".", pattern=".meanQ", full.names=T)
ffilescoast <- ffilescoast[ffilescoast=="./OCcoast_85clust_75pctmsng_oneSNPmac2_20190205_r3_out.2.meanQ"]
flistcoast <- readQ(files=ffilescoast)

flistcoast.tabq <- tabulateQ(flist)
head(flistcoast.tabq)

#flist.sum <- summariseQ(flist.tabq)
#flist.sum
# Levels: COAST="#440154FF" INLAND="#297B8EFF" SANTENA="#89D548FF"

popcoords.coast2 <- popcoords.OCreduced148[popcoords.OCreduced148$K2=="COAST",]
poplabelscoast <- data.frame(popcoords.coast2$Site, stringsAsFactors = F)
colnames(poplabelscoast) <- "Site"
poplabelscoast[] <- lapply(poplabelscoast, as.character)
#poplabels <- popcoords.OCreduced148[,3, drop=FALSE]
plotQcoast <- plotQ(flistcoast, #imgoutput = "join", 
             grplab=poplabelscoast,
             grplabpos=0.6,
             returnplot=T,
             exportplot=F,
             basesize=12,
             grplabsize=3,
             linesize=0.8,
             pointsize=4,
             ordergrp=TRUE,
             sharedindlab=FALSE,
             indlabangle = 90,
             grplabangle=45,
             clustercol=c("#440154FF","#89D548FF"),
             splab=c("Coast+Santena (2A), K=2"),
             subset=c("BBEND", "CCSP1", "CCSP2", "CCSP3", "LAGUN", "MORO", "TENA", "SANCAN"))
grid.arrange(plotQ148$plot[[1]], plotQcoast$plot[[1]])


### coast+santen, K=2, highest likelihood rep
"G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/faststructure/OCcoast_75pctmsng/OCinland_85clust_75pctmsng_oneSNPmac2_20190205_r3_out.2.meanQ"
setwd("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/faststructure/OCinland_75pctmsng")
ffilesinland <- list.files(".", pattern=".meanQ", full.names=T)
ffilesinland <- ffilesinland[ffilesinland=="./OCinland_85clust_75pctmsng_oneSNPmac2_20190205_r3_out.1.meanQ"]
flistinland <- readQ(files=ffilesinland)

flistinland.tabq <- tabulateQ(flist)
head(flistinland.tabq)

#flist.sum <- summariseQ(flist.tabq)
#flist.sum
# Levels: inland="#440154FF" INLAND="#297B8EFF" SANTENA="#89D548FF"

popcoords.inland2 <- popcoords.OCreduced148[popcoords.OCreduced148$K2=="INLAND",]
poplabelsinland <- data.frame(popcoords.inland2$Site, stringsAsFactors = F)
colnames(poplabelsinland) <- "Site"
poplabelsinland[] <- lapply(poplabelsinland, as.character)
#poplabels <- popcoords.OCreduced148[,3, drop=FALSE]
plotQinland <- plotQ(flistinland, #imgoutput = "join", 
                    grplab=poplabelsinland,
                    grplabpos=0.6,
                    returnplot=T,
                    exportplot=F,
                    basesize=12,
                    grplabsize=3,
                    linesize=0.8,
                    pointsize=4,
                    ordergrp=TRUE,
                    subset=c("FREM", "GREAT", "IM01", "IM02", "IM03", "IM06", "IM07", "IM09", "IM12", "IM14", "SHOE",
                             "LAUR", "LIME", "LOMA", "SADDLE", "STARR", "THOM", "TORO"),
                    sharedindlab=FALSE,
                    indlabangle = 90,
                    grplabangle=45,
                    clustercol=c("#297B8EFF","#89D548FF"),
                    splab=c("Inland (2B), K=1"))
grid.arrange(plotQ148$plot[[1]], plotQinland$plot[[1]], plotQcoast$plot[[1]])


setwd("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/faststructure/OCcoastsub48_75pctmsng_mac2_faststructureresults")
ffilescoastsub48 <- list.files(".", pattern=".meanQ", full.names=T)
ffilescoastsub48 <- ffilescoastsub48[ffilescoastsub48=="./OCcoastsub48_85clust_75pctmsng_oneSNPmac2_20190205_r1_out.1.meanQ"]
flistcoastsub48 <- readQ(files=ffilescoastsub48)

flistcoastsub48.tabq <- tabulateQ(flistcoastsub48)
head(flistcoastsub48.tabq)

popcoords.coastsub48 <- popcoords.OCreduced148[popcoords.OCreduced148$K3hier=="COAST",]
poplabelscoastsub48 <- data.frame(popcoords.coastsub48$Site, stringsAsFactors = F)
colnames(poplabelscoastsub48) <- "Site"
poplabelscoastsub48[] <- lapply(poplabelscoastsub48, as.character)
plotQcoastsub48 <- plotQ(flistcoastsub48, #imgoutput = "join", 
                     grplab=poplabelscoastsub48,
                     grplabpos=0.6,
                     subsetgrp=c("BBEND", "CCSP1", "CCSP2", "CCSP3", "LAGUN", "MORO"),
                     returnplot=T,
                     exportplot=F,
                     basesize=12,
                     grplabsize=3,
                     linesize=0.8,
                     pointsize=4,
                     ordergrp=TRUE,
                     # subsetgrp=c("BBEND", "CCSP1", "CCSP2", "CCSP3", "MORO", "LAGUN", "SANCAN",
                     #             "TENA", "TORO", "LOMA", "THOM", "LAUR", "SADDLE", "FREM", "STARR", "GREAT", "LIME",
                     #             "IM01", "IM03", "IM12", "IM02", "IM14", "IM06", "IM07", "IM09", "SHOE"),
                     #sortind="all",
                     sharedindlab=FALSE,
                     indlabangle = 90,
                     grplabangle=45,
                     clustercol=c("#440154FF"),
                     splab=c("Coast (3A), K=1"))
grid.arrange(plotQ148$plot[[1]], plotQinland$plot[[1]], plotQcoast$plot[[1]], plotQcoastsub48$plot[[1]])


setwd("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/faststructure/OCcoast_75pctmsng")
ffilessantena <- list.files(".", pattern=".meanQ", full.names=T)
ffilessantena <- ffilessantena[ffilessantena=="./OCcoast_85clust_75pctmsng_oneSNPmac2_20190205_r1_out.1.meanQ"]
flistsantena <- readQ(files=ffilessantena)

flistsantena.tabq <- tabulateQ(flistsantena)
head(flistsantena.tabq)

#flist.sum <- summariseQ(flist.tabq)
#flist.sum
# Levels: santena="#440154FF" INLAND="#297B8EFF" SANTENA="#89D548FF"

popcoords.santena <- popcoords.OCreduced148[popcoords.OCreduced148$K3hier=="SANTENA",]
poplabelssantena <- data.frame(popcoords.santena$Site, stringsAsFactors = F)
colnames(poplabelssantena) <- "Site"
poplabelssantena[] <- lapply(poplabelssantena, as.character)
#poplabels <- popcoords.OCreduced148[,3, drop=FALSE]
plotQsantena <- plotQ(flistsantena, #imgoutput = "join", 
                    grplab=poplabelscoast,
                    grplabpos=0.6,
                    subsetgrp=c("TENA", "SANCAN"),
                    returnplot=T,
                    exportplot=F,
                    basesize=12,
                    grplabsize=3,
                    linesize=0.8,
                    pointsize=4,
                    ordergrp=TRUE,
                    # subsetgrp=c("BBEND", "CCSP1", "CCSP2", "CCSP3", "MORO", "LAGUN", "SANCAN",
                    #             "TENA", "TORO", "LOMA", "THOM", "LAUR", "SADDLE", "FREM", "STARR", "GREAT", "LIME",
                    #             "IM01", "IM03", "IM12", "IM02", "IM14", "IM06", "IM07", "IM09", "SHOE"),
                    #sortind="all",
                    sharedindlab=FALSE,
                    indlabangle = 90,
                    grplabangle=45,
                    clustercol=c("#89D548FF"),
                    splab=c("Santena (3B), K=1"))
grid.arrange(plotQsantena$plot[[1]])

svg("G:/My Drive/Manuscripts/Spea_socal_OC_thesischapter/figures_20190211/OCreduced148_75pctmsng_FastStruc_hier_ordered_20190211.svg",
    width=16,
    height=9,
    pointsize=12,
    antialias = "default",
    onefile=FALSE)
grid.arrange(plotQ148$plot[[1]], plotQinland$plot[[1]], plotQcoast$plot[[1]], plotQcoastsub48$plot[[1]], plotQsantena$plot[[1]],
             layout_matrix=matrix(c(NA,1,NA,3,NA,2,4,5,NA),ncol=3,byrow=T))
dev.off()

cairo_pdf("G:/My Drive/Manuscripts/Spea_socal_OC_thesischapter/figures_20190211/OCreduced148_75pctmsng_FastStruc_hier_ordered_20190211.pdf",
          width=16,
          height=9,
          pointsize=12,
          antialias="default",
          fallback_resolution = 600,
          onefile=TRUE)
grid.arrange(plotQ148$plot[[1]], plotQinland$plot[[1]], plotQcoast$plot[[1]], plotQcoastsub48$plot[[1]], plotQsantena$plot[[1]],
             layout_matrix=matrix(c(NA,1,NA,3,NA,2,4,5,NA),ncol=3,byrow=T))
dev.off()
