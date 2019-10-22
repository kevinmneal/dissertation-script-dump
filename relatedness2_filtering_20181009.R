# > 0.2 seems best, just stick with what's been done in the literature
# now do it for allsocal dataset

# > 0.2 seems best, just stick with what's been done in the literature
# now do it for allsocal dataset

#rel2.345 <- read.table("vcftools_allsocal_nooutgroup_50pctmsng_nosingletons_20180910.relatedness2",
#                       header=T,
#                       stringsAsFactors = F)
relatedness2 <- read.table("G:/My Drive/Illumina Sequencing Data/20181004_ALLsocal_303indiv_highqual/nopopfile_75pctmsng/vcftools_allsocal_303indv_highqual_75pctmsng_oneSNP_nosingletons_20181004.relatedness2",
                           header=T)
rel2.303 <- relatedness2
# replace ind names with pop, then remove all comparisons that are INTERpop - only want to filter INTRApop

rel2.303.pop <- FindReplace(rel2.303,
                            Var="INDV2",
                            replaceData=popcoords.303,
                            from="Ind",
                            to="Pop")
rel2.303.pop <- FindReplace(rel2.303.pop,
                            Var="INDV1",
                            replaceData=popcoords.303,
                            from="Ind",
                            to="Pop")

rel2.303.bindpop <- cbind(rel2.303, rel2.303.pop$INDV1, rel2.303.pop$INDV2)
write.csv(rel2.303.bindpop, "vcftools_allsocal_303indv_highqual_75pctmsng_oneSNP_nosingletons_popnames_20181004.relatedness2.csv")

e<-which(!(rel2.303.bindpop$`rel2.303.pop$INDV1`==rel2.303.bindpop$`rel2.303.pop$INDV2`))
rel2.303.intrapop <- rel2.303[e, ]


rel.mat.303 <- as.matrix(xtabs(rel2.303$RELATEDNESS_PHI ~ rel2.303$INDV2 + rel2.303$INDV1))

rel.mat.303[!lower.tri(rel.mat.303)] <- 0
heatmap(rel.mat.303, symm=T)
toorelated303 <- apply(rel.mat.303,2,function(x) any(x > 0.2))
popcoords.303.toorelated <- cbind(popcoords.303, toorelated303)
write.csv(popcoords.303.toorelated, file="allsocal_303indv_toorelated.popcoords.csv")

data.new <- popcoords.nooutgroup$Ind[!apply(rel.mat.303,2,function(x) any(x > 0.2))] # selects those with no relatedness2 pairs greater than 0.2
data.new2 <- subset(popcoords.nooutgroup, Ind %in% data.new) # the final subset of unrelated individuals
head(data.new)
heatmap()



# rel2.303.pop <- FindReplace(rel2.303,
#                             Var="INDV2",
#                             replaceData=popcoords.nooutgroup,
#                             from="Ind",
#                             to="Pop")
# rel2.303.pop <- FindReplace(rel2.303.pop,
#                             Var="INDV1",
#                             replaceData=popcoords.nooutgroup,
#                             from="Ind",
#                             to="Pop")
# rel2.303.pop.sub <- rel2.303.pop[rel2.303.pop$RELATEDNESS_PHI > 0.15,]
plot(ecdf(relatedness2$RELATEDNESS_PHI), xlim=c(0, 0.5), ylim=c(0.9,1))
plot(ecdf(relatedness2$RELATEDNESS_PHI), xlim=c(0.25, 0.3), ylim=c(0.9958,0.9968))
# cutoff of 0.26 seems reasonable; 0.3 to be conservative...

#rel2.303.sub <- rel2.303[rel2.303$RELATEDNESS_PHI > 0.15,]
#rel2.303.pop.sub2 <- cbind(rel2.303.sub, rel2.303.pop.sub$INDV1, rel2.303.pop.sub$INDV2)
#write.csv(rel2.303.pop.sub, file="allsocal_50pctmsng_highr2pops.relatedness2")
#write.csv(rel2.303.pop.sub2, file="allsocal_50pctmsng_highr2indvpops.relatedness2")

# removing all above 0.2 removes some of the IM populations entirely...

# BEST cutoff might be to use the largest value that occurs BETWEEN two sites,
# because these CAN'T be full sibs (unless someone moved tadpoles or something!)
# highest interpop r2 for Irvine Mesa pops: 0.281931
# lowest-highest interpop r2 among 060709: 0.261188 (the highest 0607 > highest 0709 > highest 0609; use that 0609 value as cutoff)
# highest interpop r2 for non-IM060709 pops: 0.240741(only 1 pair above 0.24, IM01-03)
# highest interpop r2 for non-IM pops: 0.080616

toorelated303pt0806 <- apply(rel.mat.303,2,function(x) any(x > 0.08061))
popcoords.303.toorelated <- cbind(popcoords.303, toorelated303pt0806)
#write.csv(popcoords.303.toorelated, file="allsocal_toorelated_303inpopcoords.csv")

toorelated303pt2 <- apply(rel.mat.303,2,function(x) any(x > 0.2))
popcoords.303.toorelated <- cbind(popcoords.303.toorelated, toorelated303pt2)
#write.csv(popcoords.303.toorelated, file="allsocal_toorelated_303inpopcoords.csv")

toorelated303pt225 <- apply(rel.mat.303,2,function(x) any(x > 0.225))
popcoords.303.toorelated <- cbind(popcoords.303.toorelated, toorelated303pt225)
#write.csv(popcoords.303.toorelated, file="allsocal_toorelated_303inpopcoords.csv")

toorelated303pt2408 <- apply(rel.mat.303,2,function(x) any(x > 0.2408))
popcoords.303.toorelated <- cbind(popcoords.303.toorelated, toorelated303pt2408)
#write.csv(popcoords.303.toorelated, file="allsocal_toorelated_303inpopcoords.csv")

toorelated303pt25 <- apply(rel.mat.303,2,function(x) any(x > 0.25))
popcoords.303.toorelated <- cbind(popcoords.303.toorelated, toorelated303pt25)
#write.csv(popcoords.303.toorelated, file="allsocal_toorelated_303inpopcoords.csv")

toorelated303pt2612 <- apply(rel.mat.303,2,function(x) any(x > 0.2612))
popcoords.303.toorelated <- cbind(popcoords.303.toorelated, toorelated303pt2612)
#write.csv(popcoords.303.toorelated, file="allsocal_toorelated_303inpopcoords.csv")

toorelated303pt275 <- apply(rel.mat.303,2,function(x) any(x > 0.275))
popcoords.303.toorelated <- cbind(popcoords.303.toorelated, toorelated303pt275)
#write.csv(popcoords.303.toorelated, file="allsocal_toorelated_303inpopcoords.csv")

toorelated303pt282 <- apply(rel.mat.303,2,function(x) any(x > 0.282))
popcoords.303.toorelated <- cbind(popcoords.303.toorelated, toorelated303pt282)
#write.csv(popcoords.303.toorelated, file="allsocal_toorelated_303inpopcoords.csv")

toorelated303pt3 <- apply(rel.mat.303,2,function(x) any(x > 0.3))
popcoords.303.toorelated <- cbind(popcoords.303.toorelated, toorelated303pt3)
write.csv(popcoords.303.toorelated, file="allsocal_toorelated_303indv.popcoords.csv")


#

