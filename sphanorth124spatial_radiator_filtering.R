

setwd("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial")


library(vcfR)
library(radiator)
library(pals)
library(gplots)
library(strataG)
library(adegenet)
library(pcadapt)
library(devtools)
library(tess3r)
library(LEA)
#devtools::install_github("thierrygosselin/radiator", ref="0bfccf6c4daf5e71c2bbeac2424597d4e30a19b2")
devtools::install_github("thierrygosselin/radiator", ref="0bfccf6c4daf5e71c2bbeac2424597d4e30a19b2")
reload(pkgload::inst("radiator"))
#devtools::install_url("https://github.com/thierrygosselin/radiator/archive/6f109b59bf9d00252c9d156ea76fc40b3a9dfc95.zip") # working version before genomic_converter broke
#reload(pkgload::inst("radiator"))
# reinstall the above version of radiator if the latest one fails to work

# note: pcadapt output doesn't work; nor does hzar
# the other conversions I have work, but bayescan conversion is pretty slow

require(radiator)


popcoords.sphasouth178spatial <- read.csv("sphasouth178spatial_popcoords.csv",
                                          header=T,
                                          stringsAsFactors = F)


# sphasouth178spatial.folderpath <- "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/"
# sphasouth178spatialfilepathPaste <- function(filename="", 
#                                        path="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial/"){
#   paste0(path,filename)
# }

sphasouth178spatial.testdir <- paste0(getwd(),"/radiator_testing")
sphasouth178spatial.radiator.dir <- paste0(getwd(),"/radiator_final")
sphanorth158.inputvcf <- "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphanorth158.90pctmsng.minDP6.recode.vcf"

setwd(sphasouth178spatial.testdir)

sphanorth124.strata.allpops <- "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial_popcoords_radiator_strata.txt"
sphanorth124.strata.onepop <- "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial_popcoords_radiator_strata_onepop.txt"

### read_vcf
sphanorth124.tidy <- radiator::read_vcf(data = "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphanorth158.90pctmsng.minDP6.recode.vcf", 
                                    strata = "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial_popcoords_radiator_strata_onepop.txt",
                                    verbose = TRUE,
                                    parallel.core = 1,
                                    fig.upsetr=TRUE,
                                    filter.common.markers=FALSE,
                                    filter.monomorphic=TRUE)
                                    #path.folder="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial_radiator/")


#### filter_rad
### mac=3, shortLD=mac, missing=0.5, minmeanDP=1, maxmeanDP=200, SNPs per locus=6


# do not require SNPs to be common among all strata
sphanorth124.filtered <- radiator::filter_rad(data = sphanorth124.tidy, 
                                          strata = "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial_popcoords_radiator_strata_onepop.txt",
                                          #strata = sphanorth124.strata.onepop,
                                          verbose = TRUE,
                                          filename="sphasouth178spatial.radiator.50pctmsng.oneSNPmac3",
                                          parallel.core = 1,
                                          fig.upsetr=TRUE,
                                          filter.common.markers=FALSE,
                                          output=c("vcf", "fineradstructure", "tidy", "plink", "faststructure"))
#path.folder="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial_radiator/OCall207_radiator_mac3")
# output=c("vcf", "genepop", "fineradstructure", "tidy", "plink", "ldna", "stockr", "genind", "genlight", "faststructure", "hierfstat", "bayescan", "betadiv", "related"))


# convert to relevant formats, using sites as populations
OCreduced148.75pctmsng.sitesaspops.converted <- genomic_converter(data = "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/radiator_testing/filter_rad_20190405@1657/13_filtered/sphasouth178spatial.radiator.50pctmsng.oneSNPmac3.rad", 
                                              strata = paste0(sphaOC.radiator.dir,"/OCreduced148_85clust_75pctmsng_oneSNPmac2_20190205_radiator_strata.txt"),
                                              verbose = TRUE,
                                              filename="OCreduced148.radiator.75pctmsng.oneSNPmac3",
                                              parallel.core = 1,
                                              fig.upsetr=TRUE,
                                              filter.common.markers=FALSE,
                                              output=c("vcf", "genepop", "fineradstructure", "tidy", "plink", "ldna", "stockr", "genind", "genlight", "structure", "hierfstat", "bayescan", "betadiv", "related"),
                                              path.folder="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial_radiator_final/OCreduced148_radiator_mac3")





# make an NeEstimator-specific genepop with up to 90% missing data, then run in vcftools for interchrom-ld,
# then filter at different LD levels: 0.1, 0.5, 0.9, noLDfilter
# OCmost203.tidy.forNeEstimator <- radiator::read_vcf(data = "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial_radiator_final/sphasouth267_maxalleles2minmaxDP.recode.vcf", 
#                                                     strata = paste0(sphaOC.radiator.dir,"/OCmost203_radiator_strata.txt"),
#                                                     verbose = TRUE,
#                                                     parallel.core = 1,
#                                                     fig.upsetr=TRUE,
#                                                     filter.common.markers=FALSE,
#                                                     filter.monomorphic=TRUE)
# OCmost203.filtered.forNeEstimator <- radiator::filter_rad(data = OCmost203.tidy.forNeEstimator, 
#                                                           strata = paste0(sphaOC.radiator.dir,"/OCmost203_radiator_strata.txt"),
#                                                           verbose = TRUE,
#                                                           filename="OCmost203.radiator.90pctmsng.oneSNPmac3.forNeEstimator",
#                                                           parallel.core = 1,
#                                                           fig.upsetr=TRUE,
#                                                           filter.common.markers=FALSE,
#                                                           output=c("vcf", "tidy", "genepop", "strataG"),
#                                                           path.folder="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial_radiator_final/radiator_testing")
# 
# genomic_converter(data = "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial_radiator_final/radiator_testing/filter_rad_20190322@1141/13_filtered/OCmost203.radiator.75pctmsng.oneSNPmac3.rad", 
#                   strata = paste0(sphaOC.radiator.dir,"/OCmost203_radiator_strata.txt"),
#                   verbose = TRUE,
#                   filename="OCmost203.radiator.75pctmsng.oneSNPmac3.sitesaspops",
#                   parallel.core = 1,
#                   fig.upsetr=TRUE,
#                   filter.common.markers=FALSE,
#                   output=c("vcf", "genepop", "tidy", "plink",  "ldna", "stockr", "genind", "genlight", "structure", "hierfstat", "bayescan", "betadiv", "related"),
#                   path.folder="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphasouth178spatial_radiator_final/radiator_testing")
# 

