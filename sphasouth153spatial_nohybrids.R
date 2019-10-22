### some conversions using dartR

# treemix, demerelate

# remove hybrid individuals/combine pops into reasonable clusters based on sNMF to use for treemix (and maybe for genetic distance stuff (especially when it requires pops have more than 2 individuals), and phylogenetic analysis too?)

#### filter_rad
### mac=3, shortLD=mac, missing=0.5, minmeanDP=1, maxmeanDP=200, SNPs per locus=6
# do not require SNPs to be common among all strata

# filter_rad/genomic_converter still messes up the vcf output in 1.1.0
sphasouth153spatial.filtered <- radiator::filter_rad(data = "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaOCclust85/sphasouth267.vcf", 
                                              strata = "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth153spatial_nohybrids_radiator_strata_onepop.txt",
                                              #strata = sphanorth124.strata.onepop,
                                              verbose = TRUE,
                                              filename="sphasouth153spatial.nohybrids.radiator.50pctmsng.oneSNPmac3",
                                              parallel.core = 1,
                                              fig.upsetr=TRUE,
                                              filter.common.markers=FALSE,
                                              output=c("vcf", "tidy", "plink", "genind", "genlight", "genepop"))

# reload using appropriate strata
tidy.sphasouth153spatial.nohybrids.radiator.50pctmsng <- tidy_genomic_data(data="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth178spatial.radiator.50pctmsng.oneSNPmac3.vcf",
                                                                 strata="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth153spatial_nohybrids_radiator_strata.txt",
                                                                 parallel.core=1,
                                                                 filter.common.markers=FALSE)

genomic_converter(data="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/325_radiator_tidy_genomic_20190509@0113/03_tidy_vcf/radiator_20190509@0115.rad",
                  strata="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth153spatial_nohybrids_radiator_strata.txt",
                  parallel.core=1,
                  filter.common.markers=FALSE,
                  output=c("vcf", "tidy", "plink", "genind", "genlight", "genepop"))

popcoords.sphasouth153spatial.nohybrids <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/sphasouth153spatial_nohybrids_popfile_sNMFancestryK15_clusters_coords.csv",
                                          header=T,
                                          stringsAsFactors = F)

sphasouth153spatial.nohybrids.vcf <- read.vcfR("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth178spatial/330_radiator_genomic_converter_20190509@1017/radiator_data_20190509@1017.vcf")



gl.sphasouth153spatial.nohybrids.radiator.50pctmsng <- vcfR2genlight(sphasouth153spatial.nohybrids.vcf) #readRDS("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth153spatial.nohybrids/330_radiator_genomic_converter_20190509@1017/radiator_genlight_20190509@1018.RData") #<- vcfR2genlight(sphasouth153spatial.nohybrids.vcf)
genind.sphasouth153spatial.nohybrids.radiator.50pctmsng <- vcfR2genind(sphasouth153spatial.nohybrids.vcf) #<- readRDS("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/sphasouth153spatial.nohybrids/330_radiator_genomic_converter_20190509@1017/radiator_genind_20190509@1018.RData") 
pop(genind.sphasouth153spatial.nohybrids.radiator.50pctmsng) <- popcoords.sphasouth153spatial.nohybrids$poplevel2
ploidy(genind.sphasouth153spatial.nohybrids.radiator.50pctmsng) <- 2
indNames(gl.sphasouth153spatial.nohybrids.radiator.50pctmsng)
colnames(sphasouth153spatial.nohybrids.vcf@gt)[2:179] == indNames(gl.sphasouth153spatial.nohybrids.radiator.50pctmsng) # ALWAYS run this check
pop(gl.sphasouth153spatial.nohybrids.radiator.50pctmsng) <- popcoords.sphasouth153spatial.nohybrids$poplevel2
ploidy(gl.sphasouth153spatial.nohybrids.radiator.50pctmsng) <- 2


#pop(genind.sphasouth153spatial.nohybrids.radiator.50pctmsng) <- popcoords.sphasouth153spatial.nohybrids$poplevel2
treemix.sphasouth153spatial.nohybrids <- gl2treemix(gl.sphasouth153spatial.nohybrids.radiator.50pctmsng, 
                                          outfile="sphasouth153spatial.nohybrids.radiator.50pctmsng.treemix",
                                          outpath="G:\\My Drive\\Illumina Sequencing Data\\20181212_rangewide\\gitprojects\\sphasouth178spatial")
# try combining pops and removing hybrid clusters?
# use results from K=10 for sNMF, run 4 has lowest cross entropy; exclude individuals with below 50% ancestry in any
# one pop