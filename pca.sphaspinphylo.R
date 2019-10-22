### Spea genetic in R - outlier loci with pcadapt
### .R file creation date: 2018-08-28
### Last updated: 2018-11-30

### vignettes at https://bcm-uga.github.io/pcadapt/articles/index.html
### paper at https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.12592

# since the artificial pops may have undergone a bottleneck, maybe identifying outlier loci doesn't make sense,
# so maybe only run this on native populations?
# and then run resistanceGA using both neutral and outlier loci, and also faststructure


##### required packages
library(pcadapt)
library(qvalue)
library(adegenet)
library(SNPRelate)
library(vcfR)
library(poppr)
library(PopGenome)
library(dplyr)
library(ape)
library(RColorBrewer)
library(igraph)
library(gplots)
library(viridis)
library(tess3r)
library(rworldmap)
library(maps)
library(gstudio)
library(hierfstat)
library(pegas)
library(radiator)
library(pals)

#####

vcf2pcadapt("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/speaphylospatial/sphaspinphylo188spatial.50pctmsng.mac3.recode.vcf",
            out="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/speaphylospatial/sphaspinphylo188spatial.50pctmsng.mac3.pcadapt")

#vcfpcadapt <- read.pcadapt("vcftools_allsocal_303indv_highqual_20pctmsng_SNPsinallbignativepops_oneSNP_nosingletons_20181009.recode.vcf", type="vcf")
vcfpcadapt.sphaspinphylo188spatial.50pctmsng <- read.pcadapt("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/speaphylospatial/sphaspinphylo188spatial.50pctmsng.mac3.pcadapt", type="pcadapt")

popcoords.sphaspinphylo188spatial <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/sphaspinphylo188spatial/sphaspinphylo188spatial_popfile_coords.csv",
                                          header=T,
                                          stringsAsFactors = F)

poplabels.sphaspinphylo188spatial <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/speaphylospatial/sphaspinphylo188spatial.poplabels.csv",
                                              header=T,
                                              stringsAsFactors = F)

poplabels.spharangewide178spatial <- read.csv("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/speaphylospatial/spharangewide178spatial.poplabels.csv",
                                              header=T,
                                              stringsAsFactors = F)


# K should be chosen as n (number of pops/clusters) - 1
# so can use your snmf/tess3/whatever optimal K and subtract 1 (because need n-1 PCs to discriminate among n clusters)
# use K = 3 (i.e. 4 pops) or 14, based on scree plot
plot(pcadapt(input=vcfpcadapt.sphaspinphylo188spatial.50pctmsng, K=20, min.maf=0.0001), option="screeplot") # think about choice of K and min.maf here
pcadapt.sphaspinphylo188spatial.50pctmsng <- pcadapt(input=vcfpcadapt.sphaspinphylo188spatial.50pctmsng, K=2, min.maf=0.0001) # think about choice of K and min.maf here
x.sphaspinphylo188spatial.50pctmsng <- pcadapt.sphaspinphylo188spatial.50pctmsng
plot(x.sphaspinphylo188spatial.50pctmsng, option="screeplot") # 2 or 5 at min.maf=0.02 or 0.05; 2 or 4 at 0.0001 or 0.01; choice should be optimal minus 1 using Cattell's Rule
pcadaptplot.sphaspinphylo188 <- plot(x.sphaspinphylo188spatial.50pctmsng, option="scores", pop=poplabels.sphaspinphylo188spatial$Species, col=glasbey(n=length(unique(poplabels.sphaspinphylo188spatial$Species))))
plot(x.sphaspinphylo188spatial.50pctmsng, option="scores", pop=poplabels.sphaspinphylo188spatial$County, col=glasbey(n=length(unique(poplabels.sphaspinphylo188spatial$County))))
plot(x.sphaspinphylo188spatial.50pctmsng, option="scores", pop=poplabels.sphaspinphylo188spatial$pop, col=glasbey(n=length(unique(poplabels.sphaspinphylo188spatial$pop))))
plot(x.sphaspinphylo188spatial.50pctmsng, option="scores", i=2, j=3, pop=poplabels.sphaspinphylo188spatial$Species, col=glasbey(n=length(unique(poplabels.sphaspinphylo188spatial$Species))))

# if you want a density plot when using K=1:
ggplot(cbind(x.sphaspinphylo188spatial.50pctmsng$scores, poplabels.sphaspinphylo188spatial), aes(x=x.sphaspinphylo188spatial.50pctmsng$scores)) + 
  geom_density(aes(group=poplabels.sphaspinphylo188spatial$Species, 
                   col=poplabels.sphaspinphylo188spatial$Species, 
                   fill=poplabels.sphaspinphylo188spatial$Species))

pcadaptplot.sphaspinphylo188 <- pcadaptplot.sphaspinphylo188 + 
  xlab("PC1 (67.7% variance explained)") + 
  ylab("PC2 (34.3% variance explained)") +
  #ggtitle("PCA of S. hammondii (North and South) and S. intermontana")
  ggtitle(expression("PCA of"~italic("S. hammondii")~"(North and South) and"~italic("S. intermontana")))

# to get prop. of explained variance, can look at screeplot y axis, or square the $singular.values
# I don't know why the singular.values are the sqrt of the explained variance, but there it is
(x.sphaspinphylo188spatial.50pctmsng$singular.values)^2
# with maf=0.01: pc1: 71.3%; pc2: 36.7%; pc3: 
# with mac=3 (no additional maf filter): pc1: 67.7%, pc2: 34.3%, pc3: 6.4%

vcf2pcadapt("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/speaphylospatial/spharangewide178spatial.50pctmsng.mac3.recode.vcf",
            out="G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/speaphylospatial/spharangewide178spatial.50pctmsng.mac3.pcadapt")

#vcfpcadapt <- read.pcadapt("vcftools_allsocal_303indv_highqual_20pctmsng_SNPsinallbignativepops_oneSNP_nosingletons_20181009.recode.vcf", type="vcf")
vcfpcadapt.spharangewide178spatial.50pctmsng <- read.pcadapt("G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/speaphylospatial/spharangewide178spatial.50pctmsng.mac3.pcadapt", type="pcadapt")

plot(pcadapt(input=vcfpcadapt.spharangewide178spatial.50pctmsng, K=20, min.maf=0.0001), option="screeplot") # think about choice of K and min.maf here
pcadapt.spharangewide178spatial.50pctmsng <- pcadapt(input=vcfpcadapt.spharangewide178spatial.50pctmsng, K=2, min.maf=0.0001) # think about choice of K and min.maf here
x.spharangewide178spatial.50pctmsng <- pcadapt.spharangewide178spatial.50pctmsng
plot(x.spharangewide178spatial.50pctmsng, option="screeplot") # 2 or 5 at min.maf=0.02 or 0.05; 2 or 4 at 0.0001 or 0.01; choice should be optimal minus 1 using Cattell's Rule
pcadaptplot.spharangewide178 <- plot(x.spharangewide178spatial.50pctmsng, option="scores", pop=poplabels.spharangewide178spatial$Species, col=glasbey(2))
pcadaptplot.spharangewide178 <- pcadaptplot.spharangewide178 + 
  xlab("PC1 (46.1% variance explained)") + 
  ylab("PC2 (8.5% variance explained)") +
  #ggtitle("PCA of S. hammondii (North and South)") +
  ggtitle(expression("PCA of"~italic("S. hammondii")~"(North and South)"))
plot(x.spharangewide178spatial.50pctmsng, option="scores", pop=poplabels.spharangewide178spatial$County, col=glasbey(n=length(unique(poplabels.spharangewide178spatial$County))))
plot(x.spharangewide178spatial.50pctmsng, option="scores", pop=poplabels.spharangewide178spatial$pop, col=glasbey(n=length(unique(poplabels.spharangewide178spatial$pop))))

(x.spharangewide178spatial.50pctmsng$singular.values)^2
#maf 0.01: pc1 explains 49.6% of the variance; pc2, 8.9%; pc3, 4.9%
# with mac=3 (no additional maf filter): 46.1%, 8.5%, 4.7%
library(gridExtra)

CairoPNG
grid.arrange(pcadaptplot.sphaspinphylo188, pcadaptplot.spharangewide178, ncol=1)
CairoPNG(filename = "G:/My Drive/Illumina Sequencing Data/20181212_rangewide/gitprojects/speaphylospatial/pcadaptplot.sphaspin.smaller.png",
         height=5, width=7, units="in", res=320)
grid.arrange(pcadaptplot.sphaspinphylo188, pcadaptplot.spharangewide178, ncol=1)
dev.off()



plot(x.50pctmsng, option="scores", pop=poplabels.sphaspinphylo188spatial$Species, i=1,j=3, col=glasbey(n=length(unique(poplabels.sphaspinphylo188spatial$Species))))
plot(x.50pctmsng, option="scores", pop=poplabels.sphaspinphylo188spatial$Species, i=5,j=6, col=glasbey(n=length(unique(poplabels.sphaspinphylo188spatial$Species))))
plot(x.50pctmsng, option="scores", pop=poplabels.sphaspinphylo188spatial$Species, i=7,j=8, col=glasbey(n=length(unique(poplabels.sphaspinphylo188spatial$Species))))
plot(x.50pctmsng, option="scores", pop=poplabels.sphaspinphylo188spatial$Species, i=9,j=10, col=glasbey(n=length(unique(poplabels.sphaspinphylo188spatial$Species))))
plot(x.50pctmsng, option="scores", pop=poplabels.sphaspinphylo188spatial$Species)
plot(x.50pctmsng, option="manhattan")
plot(x.50pctmsng, option="qqplot")
hist(x.50pctmsng$pvalues, xlab = "p-values", main = NULL, breaks = 50, col="orange")
plot(x.50pctmsng, option = "stat.distribution")
plot(x.50pctmsng$loadings[,1])

# benjamini-yekutieli - use this one. Has been cited in other papers using pcadapt
padj.by.50pctmsng <- p.adjust(x.50pctmsng$pvalues,method="BY")
alpha <- 0.05 # use 0.01 if being conservative?
outliers.by.pt01 <- which(padj.by.50pctmsng < 0.01)
outliers.by.pt05 <- which(padj.by.50pctmsng < 0.05)
length(outliers.by.pt01)
snp_pc.50pctmsng <- get.pc(x.50pctmsng, outliers.by.pt01)

