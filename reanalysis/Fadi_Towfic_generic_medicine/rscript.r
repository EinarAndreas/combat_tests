



#setwd("C:\\Users\\vegard2013\\Dropbox\\div\\jobb\\combat_test\\papers\\Fadi_Towfic_generic_medicine")
setwd("C:\\Users\\vegard2013\\Dropbox\\div\\jobb\\combat_test\\papers\\Fadi_Towfic_generic_medicine")
#    setwd("/Users/vegardnygaard/Dropbox/div/jobb/combat_test/papers/Fadi_Towfic_generic_medicine")

#  source("http://bioconductor.org/biocLite.R")
#      biocLite("preprocessCore")
#      biocLite("GEOquery")
source("../../helper_functions.r")

library(Biobase)
library(GEOquery)
library(limma)
library(sva)
library(RColorBrewer)
library(preprocessCore)
library(car)

########################################################################
##   READ SAMPLEANNOTATION		
########################################################################
sampleannotation = read.table("data/sampleannotation.csv", sep="\t", header=TRUE,  stringsAsFactors=FALSE)
sampleannotation$code = make.names(sampleannotation$code)
sampleannotation$chip = as.character(sampleannotation$chip)
dimnames(sampleannotation)[[1]] = sampleannotation$code

#table(sampleannotation$group)
#table(sampleannotation$chip, sampleannotation$group)
sampleannotation = sampleannotation[!is.na(sampleannotation$group),] # take out 3 samples that are not assign to a group. Failed QC?
#sampleannotation = sampleannotation[!is.na(sampleannotation$covariate1),]
table(sampleannotation$chip, sampleannotation$covariate1) # the batch-covariate design.


########################################################################
##   READ Raw data file	and probe annotation file. Filter and match
########################################################################

rawdata = read.table("data/GSE40566_non_normalized.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

rawannotation = read.table("data/GPL6887_MouseWG-6_V2_0_R3_11278593_A.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, skip=8, comment.char="", quote="", fill=TRUE)

# annotation file has two separate tables, one for experimental and one at the end for controls.
# splitting and fixing the tables:
a =  rawannotation$Species=="Mus musculus"
experimentalannot = rawannotation[a,]
experimentalannot$Array_Address_Id = as.numeric(experimentalannot$Array_Address_Id)
controlannot = rawannotation[!a,]
dimnames(controlannot)[[2]] = rawannotation[rawannotation[,2]=="Array_Address_Id",]
controlannot$Array_Address_Id = suppressWarnings(as.numeric(controlannot$Array_Address_Id))
controlannot = controlannot[!is.na(controlannot$Array_Address_Id),]
controlannot=controlannot[,1:6]

probeannotation = merge(experimentalannot, controlannot, all=TRUE )
dim(probeannotation)
probeannotation = probeannotation[!duplicated(probeannotation$Array_Address_Id),]
probeannotation = probeannotation[probeannotation$Array_Address_Id %in% rawdata$ID_REF, ] # 
probeannotation$Symbol=tolower(probeannotation$Symbol)
dimnames(probeannotation)[[1]] = probeannotation$Probe_Id
dim(probeannotation)

#sort and filter probe and data similar.
datamatrix_raw = as.matrix(rawdata[,-1])
datamatrix_raw = datamatrix_raw[match( probeannotation$Array_Address_Id , rawdata$ID_REF), ]
dimnames(datamatrix_raw)[[1]] = probeannotation$Probe_Id
dim(datamatrix_raw)
dim(probeannotation)

#and match data to samples.
table(sampleannotation$code %in% dimnames(datamatrix_raw)[[2]])
table(dimnames(datamatrix_raw)[[2]] %in% sampleannotation$code)
datamatrix_raw = datamatrix_raw[, match(sampleannotation$code , dimnames(datamatrix_raw)[[2]])]


########################################################################
##   Reproduce the different version of the data as in the article
########################################################################


set.seed(100)
datamatrices = list()
datamatrices[["real_raw"]] = datamatrix_raw#[1:1000,]
datamatrices[["real_qnorm"]] = normalize.quantiles(datamatrices[["real_raw"]]) # used in paper, but seems to be the same as the limma version. Except it loses dimnames
dimnames(datamatrices[["real_qnorm"]]) = dimnames(datamatrices[["real_raw"]])

# take out samples that were talked about in the paper, assuming the others were left out from the start
#talkedabout = c("DP", "M", "N", "RS") # Drug product (also called GA), Medium, Generic (Natco), Reference
#datamatrices[["subset_raw"]] = datamatrices[["real_raw"]][, sampleannotation$covariate %in% talkedabout ]

mod = model.matrix(~as.factor(sampleannotation$covariate1))
mod0 = model.matrix(~1,data=sampleannotation)

date()
datamatrices[["real_combat_covariates"]]= as.matrix(ComBat(dat=datamatrices[["real_qnorm"]], batch=sampleannotation$chip, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE))
date()

##############################################################################
##########   Read supl tables to test that their results are reproduced
##############################################################################

table_s2 = read.table("data/table_s2.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
table_s2 = as.matrix(table_s2)
table_s2[,3:6] = as.numeric(table_s2[,3:6])
# not able to paste the pdf whitout a lot of gibberish clutter the data. Some probes are lost!
a =rowSums(is.na(table_s2[,3:6])) > 0
table(a)
#table_s2[a,]
table_s2=data.frame(table_s2[!a,], stringsAsFactors=FALSE)
for(n in 3:6) # data.frame makes the columns characters again.
{
  table_s2[,n]=as.numeric(table_s2[,n])
}

table_s3 = read.table("data/table_s3.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE) # one Probe ("1770446" or "ILM_1770446") was not found in the data and is removed from this textfile

table_s4 = read.table("data/table_s4.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
# not able to copy-paste the pdf whitout a lot of gibberish clutter the data. Some probes are lost!
a = is.na(as.numeric(table_s4$Var_generic_vs_Var_GA))
table(a)
table_s4=table_s4[!a,]
table_s4$Var_generic_vs_Var_GA = as.numeric(table_s4$Var_generic_vs_Var_GA)

table_s5 = read.table("data/table_s5.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
tmp = as.matrix(table_s5)
tmp[,3:14] = as.numeric(tmp[,3:14])
# not able to paste the pdf whitout a lot of gibberish clutter the data. Some probes are lost!
a =rowSums(is.na(tmp[,3:14])) > 0
table(a)
table_s5=data.frame(table_s5[!a,], stringsAsFactors=FALSE)
table_s5$Fold_Change = as.numeric(table_s5$Fold_Change) # refused to be numeric
#table_s5[!a,][1:10,]
	



#is.list(datamatrices[["series_matrix"]])
save.image("data/lastimage_real.rdata")
load("data/lastimage_real.rdata")


##########################################################################################
############### reproduce some of the results.   
##########################################################################################
outputdir = paste("output_", Sys.Date(), sep="")
if(!file.exists(outputdir))dir.create(outputdir)
plotdir = paste(outputdir, "/plots", sep="")
if(!file.exists(plotdir))dir.create(plotdir)
# Table S2
# Genes utilized for the tolerance method illustrated in Figure 1B.
# This standard of comparison was constructed by first identifying the top 1000 probes by absolute fold change of reference standard compared to the medium (Table S2). ... needed to have an average reference standard expression of 6.00 or higher and ones downregulated by reference standard needed to have an average medium expression of 6.00 or higher.


thisdataname = "real_combat_covariates"  # tested also for the alternative covariates settings
thiscovariatename = "covariate1"  # this produce best match to Table_S2
thisdata = datamatrices[[thisdataname]] # Using geometric mean gave 851 out of 909 of the Probes in table_S2

#thisdata = 2^datamatrices[[thisdataname]] # test for better match with Table_S2
#islog2=FALSE
#thisdataname = "real_combat_covariates_alt2"
#thiscovariatename = "covariate_alt"
meanM = rowMeans(thisdata[, sampleannotation[, thiscovariatename]=="M"]) # M is medium"
meanRS = rowMeans(thisdata[, sampleannotation[, thiscovariatename]=="RS"]) # RS is reference standard
meanDP = rowMeans(thisdata[, sampleannotation[, thiscovariatename]=="DP"]) # DP is drug product also called GA
meanN = rowMeans(thisdata[, sampleannotation[, thiscovariatename] =="N"]) # N is Natco the generic,  (or G for alternative covariates)

foldchange_RS_vs_M = meanRS - meanM
foldchange_RS_vs_M = foldchange_RS_vs_M[meanM>6 | meanRS>6]
top1000_RS_vs_M = names(foldchange_RS_vs_M[order(abs(foldchange_RS_vs_M), decreasing=TRUE)])[1:1000]
table(table_s2$ID %in% top1000_RS_vs_M) #  851 out of 909. For some reason it is not a 100% match.

#notin =!(table_s2$ID %in% top1000_RS_vs_M)
#notinnames = table_s2$ID[notin]
#cbind(table_s2[notin,], meanRS[notinnames], meanM[notinnames], FC=meanRS[notinnames]-meanM[notinnames])

# scatter plot correlation of the 4 columns in table_S2 vs. the reproduced numbers. String correlation.
pdf(file = paste(plotdir, "/reproduced_table_S2_", thisdataname,".pdf", sep=""), width=10, height=10, title="Table S2 AVG values vs reproduced") 
par(mfrow=c(2, 2))
plot(table_s2$AVG_Medium, meanM[table_s2$ID], pch=".")
plot(table_s2$AVG_Reference_Standard, meanRS[table_s2$ID], pch=".")
plot(table_s2$AVG_GA, meanDP[table_s2$ID], pch=".")
plot(table_s2$AVG_generic, meanN[table_s2$ID], pch=".")
dev.off() 

#Table S3.
#The highly variable probes that were significant by F-test in either GA or generic (see methods section) and are depicted in Figure 1A.   
# MEthods: In order to identify probes with variability induced specifically by activation (as opposed to experimental noise), we sought to identify probes that were significantly more variable when activated with either GA or generic than medium. Using an F-test, we compared GA against Medium for each probe and compared generic against Medium. We then took the set of probes where either treatment comes up to be more variable than medium (union, passes in at least at least one). For those set of probes only, we compared the variability of GA across 34 samples representing 30 batches, to the variability of generic across 11 samples representing 5 batches, utilizing an F-test to measure significance of the differences between the probes.
   
# do a f-test between  DP and medium, and N and medium. And see what is "significant" (significance cut-off not noted in the text)?
pvalsN=apply( datamatrices[["real_combat_covariates"]], 1, 
              function(x) { var.test(x[sampleannotation$covariate1=="N"], x[sampleannotation$covariate1=="M"],  
              alternative="greater")$p.value } )
pvalsDP=apply( datamatrices[["real_combat_covariates"]], 1, 
              function(x) { var.test(x[sampleannotation$covariate1=="DP"], x[sampleannotation$covariate1=="M"],  
              alternative="greater")$p.value } ) 
padjustedN = p.adjust(pvalsN,method="BH")
padjustedDP = p.adjust(pvalsDP,method="BH")
print(table(padjustedN < 0.05)) 
print(table(padjustedDP < 0.05))
# nothing significant acording to 5%FDR that was used in the Differentially expressed genes tests.

# look at the pvalues for the probes found in table_s3
pdf(file = paste(plotdir, "/reproduced_table_S3_pvals_distribution.pdf", sep=""), width=10, height=10) 
thiscolors = c("black", "brown", "red", "blue")
plot(density(pvalsN[table_s3$Probe[table_s3$Type=="generic"]]), col=thiscolors[3], main="P-values for F-test, N or DP vs. Medium")
lines(density(pvalsN), col=thiscolors[1])
lines(density(pvalsDP), col=thiscolors[2])
lines(density(pvalsDP[table_s3$Probe[table_s3$Type=="GA"]]), col=thiscolors[4])
legend("topright", legend=c("All probes, N vs. M", "All probes, DP vs. M", "Table S3 probes, N vs. M", "Table S3 probes, DP vs. M"), text.col=thiscolors)
dev.off()
# The PRobes from table_s3 are at least enriched for low p-values.
# But table_s3 is not fully reproduced. More detailed explanation is needed


# Table s4
# "Ranked list of probes by ratio of the variance in generic-activated samples to the variance in GA-activated samples."
# "Thus, we calculated for each probe the ratio of the variance in generic to the variance in GA"

# do a f-test between  N and DP. Correlate with table s4.
# f-test statistic is the same as ratio of variance?
fstatistics=apply( datamatrices[["real_combat_covariates"]], 1, 
              function(x) { var.test(x[sampleannotation$covariate1=="N"], x[sampleannotation$covariate1=="DP"],  
                                     alternative="greater")$statistic } )
fstatistics = fstatistics[order(fstatistics, decreasing=TRUE)]
table(names(fstatistics)[1:100] %in% table_s4$Probe) # not a complete reproduction.
# look at the f-statistic for the probes found in table_s3 vs all f-statistics
pdf(file = paste(plotdir, "/reproduced_table_S4_fstatistic_distribution.pdf", sep=""), width=10, height=10) 
par(mfrow=c(2, 1))
plot(table_s4[,3], fstatistics[table_s4$Probe], main="F-statitics from table S4 vs. reproduced")
thiscolors = c("black", "red")
plot(density(fstatistics), col=thiscolors[1], main="Reproduced F-statistics for N vs. DP")
lines(density(fstatistics[table_s4$Probe]), col=thiscolors[2])
legend("topright", legend=c("All probes", "Table S4 probes"), text.col=thiscolors)
dev.off()



# Table s5
#"Comparison of expression in GA to expression in GA for each probe, including fold change, ANOVA, LIMMA with background subtraction, comparative marker selection by signal-to-noise ratio, comparative marker selection by t-test, and the Wilcoxon non-parametric method."
#"Finally, genes that have significantly higher expression in samples activated by generic than in samples activated by GA by 4 different parametric methods (Table S5)"
#"Probes with Fold Change <  1.10 notdisplayed"

# Table s5 seems to be a comparison of expression in  GA (drug product) vs. generic (N) , listing the probes with greatest fold change (>abs(1.1)??) and the p-value calculated by several different tests. There seems to be no description of how the fold-change was calculated and after an inspection they look weird.

meanDP = rowMeans(thisdata[, sampleannotation$covariate1=="DP"]) # DP is drug product also called GA
meanN = rowMeans(thisdata[, sampleannotation$covariate1 =="N"]) # N is Natco the generic,  (or G for alternative covariates)

# fold change between N and DP, calculated as geometric mean.
x = table_s5$Fold_Change
y = (meanDP - meanN)[table_s5$Probe]
lim =c( min(   c(min(x), min(y))  ), max( c(max(x), max(y)) ))
plot(x, y , ylim=lim, xlim=lim, pch=".", xlab="Table_S5", ylab="reproduced", main="FC comparison") # Not a perfect match. And the table_s5 numbers seems offset.


# table S2 has the AVG for both DP and N listed
fc_table_s2 = table_s2$AVG_GA - table_s2$AVG_generic
names(fc_table_s2) = table_s2$ID
x = table_s5$Fold_Change
y = fc_table_s2[table_s5$Probe]
lim =c( min(   c(min(x, na.rm=T), min(y, na.rm=T))  ), max( c(max(x, na.rm=T), max(y, na.rm=T)) ))
plot(x, y , ylim=lim, xlim=lim, pch="." ,  xlab="Table_S5", ylab="Table_S2", main="FC comparison") # perfect match.

# But this is odd, we were able to reproduce table s2 and table s2 matches table s5, but table s5 fold changes does not match our reproduced values as good as the match for table s2.

# The reason that the line is not as straight for the reproduce FC's for table s5 as for the mean expressions in table S2 is because of different scale in the plot and the fact that the small difference will be higher for a ratio (combining to values with a individual deviation from the tables.)
x = (meanDP - meanN)[table_s2$ID]
y = fc_table_s2
plot(x, y , ylim=lim, xlim=lim, pch="." ,  xlab="reproduced", ylab="Table_S2", main="FC comparison") # perfect match.


# Still there is a offset in the numbers.
commonprobes = table_s5$Probe[table_s5$Probe %in% table_s2$ID]
fc_table_s5 = table_s5$Fold_Change[match(commonprobes ,table_s5$Probe)]
cbind(fc_table_s2[commonprobes], fc_table_s5, fc_table_s5-fc_table_s2[commonprobes], fc_table_s5/fc_table_s2[commonprobes])
# seems that some unacounted transformation has occured on the fold changes.


# reproduce the Nom P-value columns
limma_ret= getdifftab_limma(edata=datamatrices[["real_combat_covariates"]], condition=sampleannotation[,"covariate1"], contrast="DP-N" )
limma_p = limma_ret$p
names(limma_p) = limma_ret$gene
DP = sampleannotation$covariate1=="DP"
N = sampleannotation$covariate1=="N"
ttest_p = apply(datamatrices[["real_combat_covariates"]], 1, function(x) { t.test(x[DP], x[N])$p.value } )
anova_p = apply(datamatrices[["real_combat_covariates"]], 1, function(x){ summary(aov(x[DP|N]~sampleannotation$covariate1[DP|N]))[[1]][1,5]}  )  # aov er kanskje feil for dette, Er ikke live mange prøver i guppene.
wilcoxon_p = apply(datamatrices[["real_combat_covariates"]], 1, function(x) { wilcox.test(x[sampleannotation$covariate1=="DP"], x[sampleannotation$covariate1=="N"])$p.value } )

pdf(file = paste(plotdir, "/reproduced_table_S5_.pdf", sep=""), width=10, height=10) 
par(mfrow=c(2, 2))
#plot(table_s5$Fold_Change, (meanDP - meanN)[table_s5$Probe] ,  ylim=c(-1.5,1.5), xlim=c(-1.5,1.5), pch=".", xlab="Table_S5", ylab="reproduced", main="FC comparison")
plot(log10(table_s5$t_test_p), log10(ttest_p[table_s5$Probe]), pch=".", main="t-test p")
plot(log10(table_s5$limma_p), log10(limma_p[table_s5$Probe]), pch=".", main="limma p")
plot(log10(table_s5$anova_p), log10(anova_p[table_s5$Probe]), pch=".", main="anova p")
plot(log10(table_s5$wilcoxon_p), log10(wilcoxon_p[table_s5$Probe]), pch=".", main="wilcoxon p")
dev.off() 
# The p-values does not make a straight line, but they are low. But maybe all are low?

pdf(file = paste(plotdir, "/reproduced_table_S5_pvalue_distribution.pdf", sep=""), width=10, height=10)
par(mfrow=c(2, 2))
thiscolors=c("black", "red")
hist(ttest_p, border=thiscolors[1], main="Reproduced t-test p-value distribution.", breaks=100)
hist(ttest_p[table_s5$Probe], border=thiscolors[2], add=T, breaks=100)
legend("topright", legend=c("all probes", "S5 probes"), text.col=thiscolors)
hist(limma_p, border=thiscolors[1], main="Reproduced limma p-value distribution.", breaks=100)
hist(limma_p[table_s5$Probe], border=thiscolors[2], add=T, breaks=100)
legend("topright", legend=c("all probes", "S5 probes"), text.col=thiscolors)
hist(anova_p, border=thiscolors[1], main="Reproduced anova p-value distribution.", breaks=100)
hist(anova_p[table_s5$Probe], border=thiscolors[2], add=T, breaks=100)
legend("topright", legend=c("all probes", "S5 probes"), text.col=thiscolors)
hist(wilcoxon_p, border=thiscolors[1], main="Reproduced wilcoxon p-value distribution.", breaks=100)
hist(wilcoxon_p[table_s5$Probe], border=thiscolors[2], add=T, breaks=100)
legend("topright", legend=c("all probes", "S5 probes"), text.col=thiscolors)
dev.off()






##########################################################################################
############### The tables without the use of ComBat.
##########################################################################################

# table S5, The differentially expressed genes between DP and N can (and probably should) be calculated taking batch into the account while doing the test. Starting with the quantile normalized data and include the batch desgign in the model i.e Using BEASTs (Batch Effect Aware Statistical Test).



# reproduce the Nom P-value columns
limma_ret= getdifftab_limma(edata=datamatrices[["real_qnorm"]], condition=sampleannotation[,"covariate1"], contrast="DP-N", block=sampleannotation$chip)
limma_p = limma_ret$p
names(limma_p) = limma_ret$gene
DP = sampleannotation$covariate1=="DP"
N = sampleannotation$covariate1=="N"
#ttest_p = apply(datamatrices[["real_qnorm"]], 1, function(x) { t.test(x[DP], x[N])$p.value } )
cov1_fac = factor(sampleannotation$covariate1[DP|N])
chip_fac = factor(sampleannotation$chip[DP|N])
anova_p1 = apply(datamatrices[["real_qnorm"]], 1, function(x){ summary(aov(x[DP|N]~chip_fac*cov1_fac))[[1]][2,5]}  ) 
anova_p2 = apply(datamatrices[["real_qnorm"]], 1, function(x){ summary(aov(x[DP|N]~chip_fac+cov1_fac))[[1]][2,5]}  ) 
co = list(cov1_fac=contr.sum, chip_fac=contr.sum)
anova_p3 = apply(datamatrices[["real_qnorm"]], 1, function(x){ Anova(lm(x[DP|N] ~  cov1_fac + chip_fac,  contrasts=co), type=3)[[4]][2]}  ) 
#wilcoxon_p = apply(datamatrices[["real_combat_covariates"]], 1, function(x) { wilcox.test(x[sampleannotation$covariate1=="DP"], x[sampleannotation$covariate1=="N"])$p.value } )

pdf(file = paste(plotdir, "/BEAST_table_S5_pvalue_distribution.pdf", sep=""), width=10, height=10)
par(mfrow=c(2, 2))
thiscolors=c("black", "red")
#hist(ttest_p, border=thiscolors[1], main="Reproduced t-test p-value distribution.", breaks=100)
#hist(ttest_p[table_s5$Probe], border=thiscolors[2], add=T, breaks=100)
#legend("topright", legend=c("all probes", "S5 probes"), text.col=thiscolors)
hist(limma_p, border=thiscolors[1], main="Reproduced limma p-value distribution.", breaks=100)
hist(limma_p[table_s5$Probe], border=thiscolors[2], add=T, breaks=100)
legend("topright", legend=c("all probes", "S5 probes"), text.col=thiscolors)
hist(anova_p, border=thiscolors[1], main="Reproduced anova p-value distribution.", breaks=100)
hist(anova_p[table_s5$Probe], border=thiscolors[2], add=T, breaks=100)
legend("topright", legend=c("all probes", "S5 probes"), text.col=thiscolors)
#hist(wilcoxon_p, border=thiscolors[1], main="Reproduced wilcoxon p-value distribution.", breaks=100)
#hist(wilcoxon_p[table_s5$Probe], border=thiscolors[2], add=T, breaks=100)
#legend("topright", legend=c("all probes", "S5 probes"), text.col=thiscolors)
dev.off()


anova_p_old = anova_p
x  = datamatrices[["real_qnorm"]][100,]
cov1_fac = factor(sampleannotation$covariate1[DP|N])
chip_fac = factor(sampleannotation$chip[DP|N])
ret = summary(aov(x[DP|N]~cov1_fac + chip_fac ))
ret = summary(aov(x[DP|N]~ chip_fac *cov1_fac ))
replications(  x[DP|N]~cov1_fac + chip_fac, data = NULL, na.action)

summary(aov(x[DP|N]~sampleannotation$covariate1[DP|N]))

x  = datamatrices[["real_qnorm"]][100,]
Anova(lm(x[DP|N] ~  cov1_fac + chip_fac,  contrasts=list(cov1_fac=contr.sum, chip_fac=contr.sum)), type=3)


ret = Anova(lm(x[DP|N] ~   chip_fac + cov1_fac,  contrasts=list(cov1_fac=contr.sum, chip_fac=contr.sum)), type=3)

Anova(lm(x[DP|N] ~  cov1_fac + chip_fac,  contrasts=co), type=3)[[4]][3]

########################################################################
##   Create different version of the data as in the article for sanity checks and alternative analysis
########################################################################


datamatrices[["random_raw"]] = matrix(rnorm(length(datamatrices[["real_raw"]])), nrow=nrow(datamatrices[["real_raw"]]), ncol=ncol(datamatrices[["real_raw"]]) , dimnames=dimnames(datamatrices[["real_raw"]]))
#datamatrices[["real_qnorm2"]] = normalizeBetweenArrays(datamatrices[["real_raw"]])

mod = model.matrix(~as.factor(sampleannotation$covariate1))
mod0 = model.matrix(~1,data=sampleannotation)

date()
datamatrices[["real_combat_unsupervised"]]= as.matrix(ComBat(dat=datamatrices[["real_qnorm"]], batch=sampleannotation$chip, mod=NULL, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE))
date()

date()
datamatrices[["random_combat_covariates"]]= as.matrix(ComBat(dat=datamatrices[["random_raw"]], batch=sampleannotation$chip, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE))
datamatrices[["random_combat_unsupervised"]]= as.matrix(ComBat(dat=datamatrices[["random_raw"]], batch=sampleannotation$chip, mod=NULL, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE))
date()

###### make shuffled dataset, where N and DP data are shuffled within a batch
nshuffleddatasets = 3
for(i in 1:nshuffleddatasets)
{
  x = shufflesamplesinbatch( sampleannotation$code, sampleannotation$chip, sampleannotation$covariate1, c("N", "DP") )
  print(  sum(sampleannotation[,"covariate1"] != sampleannotation[x,"covariate1"]  ))  # see how many got change. Has 11 N's)	
  thisqnorm = datamatrices[["real_qnorm"]][, x]
  thiscombat = as.matrix(ComBat(dat=thisqnorm, batch=sampleannotation$chip, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE))
  datamatrices[[ paste("permute", i, "_qnorm", sep="")]] = thisqnorm
  datamatrices[[ paste("permute", i, "_combat_covariates", sep="")]] = thiscombat
}


##### make balanced datasets, i.e 7 N vs 7 DP on 7 chips
a = sampleannotation$covariate1 %in% c("N", "DP")
nbalanceddatasets = 3
for(i in 1:nbalanceddatasets)
{
  x = drawbatchbalanceddsamples(sampleannotation$code[a], sampleannotation$chip[a], sampleannotation$covariate1[a])
  print(table(sampleannotation$chip, sampleannotation$covariate1))
  thisqnorm = normalizeBetweenArrays( datamatrices[["real_qnorm"]][,x])
  mod_balanced = model.matrix(~as.factor(sampleannotation[x,"covariate1"]))
  mod0_balanced = model.matrix(~1,data=sampleannotation[x,])
  thiscombat = as.matrix(ComBat(dat=thisqnorm, batch=sampleannotation[x,"chip"], mod=mod_balanced, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)) # but will not be normaldistributed?
  datamatrices[[ paste("balanced", i, "_qnorm", sep="")]] = thisqnorm
  datamatrices[[ paste("balanced", i, "_combat_covariates", sep="")]] = thiscombat
}

#datamatrices[["balanced_combat"]]= as.matrix(ComBat(dat=datamatrices[["balanced_qnorm"]], batch=sampleannotation[x,"chip"], mod=mod_balanced, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE))
#datamatrices[["balanced_unsupervised"]]= as.matrix(ComBat(dat=datamatrices[["balanced_qnorm"]], batch=sampleannotation[x,"chip"], mod=NULL, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE))
#datamatrices[["balanced_meanadjusted"]] = pamr.batchadjust(   list(x=datamatrices[["balanced_qnorm"]],batchlabels=factor(sampleannotation[datanames,"chip"]))   )$x


#### add the data presented on geo as series_matrix, presumably normalized with partek.
#### Be awere that only 34010 probes are present here, and the probeorder is different than the other data sets
# geodataset, partek normalized
thisgeoseries=getGEO("GSE40566")
thiseset = thisgeoseries[[1]]
thismatrix = exprs(thiseset)
#table(dimnames(thismatrix)[[2]] %in% sampleannotation$geoaccession)
datamatrices[["series_matrix"]] = thismatrix[, sampleannotation$geoaccession]
dimnames(datamatrices[["series_matrix"]])[[2]] = sampleannotation$code
sum(duplicated(dimnames(thismatrix)[[1]]))
table(dimnames(thismatrix)[[1]] %in% probeannotation$Probe_Id)
datamatrices[["balanced_series_matrix"]] = datamatrices[["series_matrix"]][, dimnames(datamatrices[["balanced_qnorm"]])[[2]]]



table_s5




table(table_s4$Probe %in% table_s3$Probe)



table("ILMN_1237927" %in% names(pvalsDP))

a = table_s3$Probe[table_s3$Type=="GA"]


padjustedN = p.adjust(pvalsN,method="BH")
   

   
print(table(padjusted < 0.05))   
   print(table(fstats > 1))  
   #print(table( < 0.05))  
hist(Mpval,breaks=100)
   
  
x = datamatrices[["real_combat_covariates"]][1000,]
#x[sampleannotation$covariate1=="DP"]  
   var.test(x[sampleannotation$covariate1=="DP"]  , x[sampleannotation$covariate1=="M"] , alternative="greater" )$statistic

ret= getdifftab_limma(datamatrices[["real_combat_covariates"]], sampleannotation[,"covariate1"], "DP-N" )
table(ret$padjusted<0.05)
#ret = ret[order(ret$p),]
hist(ret$p,breaks=100)
hist(ret$fc,breaks=100)
found = as.character(ret[abs(ret$fc)>0.15, "gene"])
length(found)

table(found %in% table_s5$Probe)
table(table_s5$Probe %in% ret$gene[order(abs(ret$fc), decreasing =TRUE)][1:1000])

hist(ret[table_s5$Probe,"fc"],breaks=100)

plot(density(ret[,"fc"]), col="black")
lines(density(ret[table_s5$Probe,"fc"]), col="red")


## se om resultatet er reproduset
a = which(probeannotation$Symbol=="foxp3")
a = which(probeannotation$Symbol=="ifng")
ret[a,]
table(dimnames(probeannotation)[[1]]==ret$gene)

	#combat_dat = datamatrix
	#pcombat = f.pvalue(combat_dat, mod, mod0)
	#qcombat = p.adjust(pcombat,method="BH")	#
	#table(qcombat<0.05)

#thisnames = dimnames(datamatrices[["balanced_combat"]])[[2]][c(1:2,8:9)]
thisnames = dimnames(datamatrices[["balanced_combat"]])[[2]]
ret= getdifftab_limma(datamatrices[["balanced_combat"]][,thisnames], sampleannotation[thisnames,"covariate1"], "DP-N" )
table(ret$padjusted<0.05)
hist(ret$p,breaks=100)
	#sampleannotation[thisnames,]

ret = ret[order(ret$p),]
	


found = vector()
for(i in 1:nshuffleddatasets)
{
	ret= getdifftab_limma(datamatrices[[ paste("permute",i,"_combat_covariates", sep="")]], sampleannotation[,"covariate1"], "DP-N" )
	print(table(ret$padjusted<0.05))
	#ret = ret[order(ret$p),]
	found = c(found, ret[ret$padjusted<0.05, "gene"])
	#hist(ret$p,breaks=100)
}

table(table(found))




found = vector()
for(i in 1:nbalanceddatasets)
{
	sa = sampleannotation[dimnames(datamatrices[[paste("balanced", i, "_qnorm", sep="")]])[[2]],]
	M = datamatrices[[paste("balanced", i, "_qnorm", sep="")]][,  sa$covariate1=="DP"]-datamatrices[[paste("balanced", i, "_qnorm", sep="")]][,  sa$covariate1=="N"]
	#t.test(M)
	Mpval=apply(M, 1, function(x) { t.test(x)$p.value } )
	padjusted = p.adjust(Mpval,method="BH")
	print(table(padjusted < 0.05))
	#hist(Mpval,breaks=100)

}

table(table(found))





ret= getdifftab_limma(datamatrices[["real_qnorm"]], sampleannotation[,"covariate1"], "DP-N" , block=sampleannotation[,"chip"])
table(ret$padjusted<0.05)
hist(ret$p,breaks=100)
ret = ret[order(ret$p),]
ret[1:10,]

ret["ILMN_2635132",]

	#thisnames = dimnames(datamatrices[["balanced_combat"]])[[2]][c(1:2,8:9)]
thisnames = dimnames(datamatrices[["balanced_unsupervised"]])[[2]]
ret= getdifftab_limma(datamatrices[["balanced_unsupervised"]][,thisnames], sampleannotation[thisnames,"covariate1"], "DP-N" )
table(ret$padjusted<0.05)
hist(ret$p,breaks=100)
	#sampleannotation[thisnames,]

	
thisnames = dimnames(datamatrices[["balanced_raw"]])[[2]]
ret= getdifftab_limma(datamatrices[["balanced_raw"]][,thisnames], sampleannotation[thisnames,"covariate1"], "DP-N" )
table(ret$padjusted<0.05)
hist(ret$p,breaks=100)
ret = ret[order(ret$p),]

	
thisnames = dimnames(datamatrices[["balanced_series_matrix"]])[[2]]
ret= getdifftab_limma(datamatrices[["balanced_series_matrix"]][,thisnames], sampleannotation[thisnames,"covariate1"], "DP-N" )
table(ret$padjusted<0.05)
hist(ret$p,breaks=100)


thisnames = dimnames(datamatrices[["balanced_meanadjusted"]])[[2]]
tmp = sampleannotation[thisnames,"covariate1"][c(8,9,10,11,12,13,7,1,2,3,4,5,6,14)]
ret= getdifftab_limma(datamatrices[["balanced_meanadjusted"]][,thisnames], tmp, "DP-N" )
table(ret$padjusted<0.05)
hist(ret$p,breaks=100)


thisnames = dimnames(datamatrices[["balanced_unsupervised"]])[[2]]
tmp = sampleannotation[thisnames,"covariate1"][c(8,9,10,11,12,13,7,1,2,3,4,5,6,14)]
tmp = sampleannotation[thisnames,"covariate1"][c(1,2,3,11,12,6,7,   8,9,10,4,5,13,14)]
ret= getdifftab_limma(datamatrices[["balanced_unsupervised"]][,thisnames], tmp, "DP-N" )
table(ret$padjusted<0.05)
hist(ret$p,breaks=100)

a = c(1,2,3,4,5,6,7)
b = c(8,9,10,11,12,13,14)
#Mpval=apply(datamatrices[["balanced_unsupervised"]][,thisnames], 1, function(x) { t.test(x[a], x[b])$p.value } )
Mpval=apply(datamatrices[["balanced_unsupervised"]][,thisnames], 1, function(x) { wilcox.test(x[a], x[b])$p.value } )
hist(Mpval,breaks=100)

i=1
a = datamatrices[["balanced_series_matrix"]][,thisnames[i] ]
b = datamatrices[["balanced_meanadjusted"]][names(a),thisnames[i] ]
plot(a,b, pch=".")


thisnames = dimnames(datamatrices[["balanced_series_matrix"]])[[2]]
thisnames = thisnames[order(sampleannotation[thisnames, "chip"])]
ret = ret[order(ret$p),]
de= as.character(ret$gene[1:1000])
i=14
a = datamatrices[["balanced_series_matrix"]][,thisnames[i] ]
b = datamatrices[["balanced_combat"]][names(a),thisnames[i] ]
de = de[de %in% names(a)]
#plot(a[de],b[de], pch=".")
plot(density( a[de]-b[de] ))

#de= as.character(ret$gene[20000:20009])
de = de[de %in% names(a)]
datamatrices[["balanced_combat"]][de[1:3],thisnames ]
datamatrices[["balanced_series_matrix"]][de[1:3],thisnames ]
datamatrices[["balanced_meanadjusted"]][de[1:3],thisnames ]

datamatrices[["balanced_series_matrix"]][de[1:3],thisnames ]
datamatrices[["balanced_raw"]][de[1:3],thisnames ]

datamatrices[["balanced_series_matrix"]][1:5,1:5]
datamatrices[["balanced_raw"]][1:5,1:5]


data1 = datamatrices[["balanced_qnorm"]]


de= as.character(ret$gene[1:1000])

sa = sampleannotation[dimnames(datamatrices[["balanced1_qnorm"]])[[2]],]
M = datamatrices[["balanced1_qnorm"]][,  sa$covariate1=="DP"]-datamatrices[["balanced1_qnorm"]][,  sa$covariate1=="N"]
#t.test(M)
Mpval=apply(M, 1, function(x) { t.test(x)$p.value } )
padjusted = p.adjust(Mpval,method="BH")
table(padjusted < 0.05)
hist(Mpval,breaks=100)


#t.test(datamatrices[["balanced1_qnorm"]][de[1],  sa$covariate1=="DP"],   datamatrices[["balanced1_qnorm"]][de[1],  sa$covariate1=="N"] , paired=TRUE)






ret = getdifftab_sva(datamatrices[["balanced_unsupervised"]][,thisnames], sampleannotation[thisnames,], mod_balanced )




ret= getdifftab_limma(datamatrices[["combat_unsupervised"]], sampleannotation[,"covariate1"], "DP-N" )
table(ret$padjusted<0.05)
hist(ret$p,breaks=100)





ret = ret[order(ret$p),]


boxplot(datamatrices[["balanced_raw"]])
	
	
ret3 = getdifftab_limma(datamatrices[["series_matrix"]], sampleannotation$group, "Verified-Generic" )
table(ret3$padjusted<0.05)

ret= getdifftab_limma(datamatrices[["series_matrix"]], sampleannotation$covariate1, "DP-N" )
table(ret$padjusted<0.05)
hist(ret$p,breaks=100)


ret= getdifftab_limma(datamatrices[["combat_covariates"]], sampleannotation$covariate1, "DP-N" )
table(ret$padjusted<0.05)
hist(ret$p,breaks=100)

ret = getdifftab_limma(datamatrices[["combat_random"]], sampleannotation$covariate1, "DP-N" )
table(ret$padjusted<0.05)
hist(ret$p,breaks=100)
	
ret2 = getdifftab_limma(datamatrices[["combat_unsupervised"]], sampleannotation$group, "Verified-Generic" )
table(ret2$padjusted<0.05)


ret = getdifftab_limma(datamatrices[["combat_covariates"]], sampleannotation$group, "Verified-Generic" )
table(ret$padjusted<0.05)

table(ret$gene==probeannotation$Probe_Id)



probeannotation$Symbol=tolower(probeannotation$Symbol)
a = which(probeannotation$Symbol=="foxp3")
#a = which(probeannotation$Symbol=="cxcl9")
#probeannotation[a,]
ret[which(probeannotation$Symbol=="ifng"),]

outputdir = paste("output_", Sys.Date(), sep="")
if(!file.exists(outputdir))dir.create(outputdir)

colorpal = c( "green", "blue", "orange", "cyan", "black", "red", brewer.pal(9,"Set1"), brewer.pal(8,"Set2"), brewer.pal(8,"Set3"))
sum(duplicated(colorpal ))
colorpal=c(colorpal,colorpal,colorpal,colorpal,colorpal)

colorassignments = list()
for(i in 1:dim(sampleannotation)[2])
{
	colorassignments[[names(sampleannotation)[i]]] = colorpal[as.factor(unique(sampleannotation[,i]))]
	names(colorassignments[[names(sampleannotation)[i]]]) = as.character(unique(sampleannotation[,i]))
}



boxplotqc = function(datamatrices, plotdir, colorassignments, sa=NULL, plotsize=1200 )
{	
	for(thisdatatype in names(datamatrices))
	{
		print(paste("plotting ", thisdatatype))
		if(is.null(sa))
		{
			ncovariates=0
			thissa=NULL
		}
		if(!is.null(sa))
		{
			thissa = sa[dimnames(datamatrices[[thisdatatype]])[[2]],]
			ncovariates=dim(thissa)[2]
		}
		for(i in 1:ncovariates)
		{
			a = order(thissa[,i])
			print(paste("plotting ", thisdatatype, names(thissa)[i]))
			thisname=paste(plotdir, "/boxplot_", thisdatatype, "_", names(thissa)[i], ".png", sep="")
			png(file = thisname, pointsize = 12, width = plotsize, height = plotsize, units="px")
			if(is.null(sa)){thiscolor = "black"}else{thiscolor =colorassignments[[names(thissa)[i]]][thissa[a,i]]}
			boxplot(datamatrices[[thisdatatype]][,a], col=thiscolor , main=paste(names(thissa)[i], " colored", sep=""))
			if(length(unique(thissa[a,i]))<8)
			{
				annotlabels=paste(unique(thissa[a,i]),table(thissa[a,i])[unique(thissa[a,i])])
				legend("topleft" , legend=annotlabels, text.col=colorassignments[[names(thissa)[i]]][unique(thissa[a,i])], cex=3, bty="n")
			}
			dev.off()
		}
	}
}


#plotdir=outputdir
#sa=sampleannotation[,c("chip", "group", "slot")]
#plotsize=1200

boxplotqc(datamatrices, outputdir, colorassignments,sampleannotation[,c("chip", "group", "slot")])
plotmanypca(datamatrices, outputdir, sampleannotation[,c("chip", "group", "slot")], colorassignments )


boxplotqc(datamatrices["balanced_series_matrix"], outputdir, colorassignments,sampleannotation[,c("chip", "group", "slot")])
plotmanypca(datamatrices["balanced_series_matrix"], outputdir, sampleannotation[,c("chip", "group", "slot")], colorassignments )

boxplotqc(datamatrices[8:11], outputdir, colorassignments,sampleannotation[,c("chip", "group", "slot")])
plotmanypca(datamatrices[8:11], outputdir, sampleannotation[,c("chip", "group", "slot")], colorassignments )


#plotdirname=outputdir
#sa = sampleannotation[,c("chip", "group", "slot")]
#plotsize=1200

plotmanypca = function(datamatrices, plotdirname, sa, colorassignments, plotsize=1200)
{
	
	for(thisdatatype in names(datamatrices))
	{
		#### be sure to sort sampleannotation as data
		thissa = sa[dimnames(datamatrices[[thisdatatype]])[[2]],]
		# plot also the pvca plot.
		# plotpvca(matrices[[thisname]], sa, paste(plotdirname, "/pvca_", transcript_type,"_", thisname, ".png", sep=""))

		date()
		thisprcomp=prcomp( (t(na.omit(datamatrices[[thisdatatype]]))))
		date()
		
		for(i in 1:dim(thissa)[2])
		{			
			thisname=paste(thisdatatype, "_", names(thissa)[i], sep="")
			plotonepca( name=thisname, thisannot=thissa[,i], thisprcomp=thisprcomp, plotdirname=plotdirname, thiscolorassignments=colorassignments[[names(thissa)[i]]])
		}
	}
}

#name=thisname
#thisannot=thissa[,2]
#thiscolorassignments=colorassignments[[names(thissa)[2]]]

plotonepca = function(name, thisannot, thisprcomp, plotdirname, thiscolorassignments, plotsize=1200)
{


	
	
	#a=match(levels(thisannot), unique(thisannot))
	#a=a[!is.na(a)]
	#annotlabels=annotlabels[a]
	#annotlabelscolors=annotlabelscolors[a]
		
	png(file = paste(plotdirname , "/pca_",name, ".png", sep=""),pointsize = 12, width = plotsize, height = plotsize, units="px")
	
	plot(thisprcomp$x,  col=thiscolorassignments[thisannot], cex=2 , lwd=5, main=name, cex.main=4)
	if(length(unique(thisannot))<8)
	{
		annotlabels=paste(unique(thisannot),table(thisannot)[unique(thisannot)])
		annotlabelscolors=thiscolorassignments[unique(thisannot)]
		legend("bottomleft" , legend=annotlabels, text.col=annotlabelscolors, cex=3, bty="n")
	}
	dev.off()

}











#dim(pheno)
#names(pheno)
#t(pheno[1,])
#as.character(pheno$title)
#pheno[,c( "title", "characteristics_ch1", "description", "description.1")]
#table(as.character(pheno$title) %in% sampleannotation$Code)




save.image("data/lastimage.rdata")
load("data/lastimage.rdata")
t(thisgpl@dataTable@table[2,])

load("../../geodatasets/breast.rdata")








dim(probannot)

table(rawdata$ID_REF %in% experimentalannot$Array_Address_Id)
table(rawdata$ID_REF %in% controlannot$Array_Address_Id)
table(experimentalannot$Array_Address_Id %in% controlannot$Array_Address_Id)
a =experimentalannot$Array_Address_Id %in% controlannot$Array_Address_Id
experimentalannot[a,]
a=controlannot$Array_Address_Id %in% experimentalannot$Array_Address_Id
controlannot[a,]


sum(duplicated(experimentalannot$Array_Address_Id))
sum(duplicated(controlannot$Array_Address_Id))

x = table(controlannot$Array_Address_Id)
controlannot[controlannot$Array_Address_Id %in% names(x[x>1]),]
x[x>1]



sum(duplicated(rawdata$ID_REF))
dimnames(rawdata)[[1]][1:100]
dim(rawdata)
dim(experimentalannot)
dim(controlannot)


is.na(names(controlannot))


table(rawannotation$Species)

library(lumi)
x.lumi <- lumiR.batch("data/GSE40566_non_normalized.txt", sampleInfoFile="data/sampleannotation.csv")






set.seed(100)
a = bigsampleannotation$series == bigsampleannotation$series[1]
sampleannotation$dummy_geo_accession=sample (bigsampleannotation$geo_accession[a], nrow(sampleannotation) )


datamatrix = bigmatrix[,sampleannotation$dummy_geo_accession]
datamatrix = normalizeBetweenArrays(datamatrix)
#datamatrix[,]=rnorm( length(datamatrix))

	mod = model.matrix(~as.factor(sampleannotation$group))
	mod0 = model.matrix(~1,data=sampleannotation)
	combat_dat = as.matrix(ComBat(dat=datamatrix, batch=sampleannotation$chip, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE))
	#combat_dat = datamatrix
	#pcombat = f.pvalue(combat_dat, mod, mod0)
	#qcombat = p.adjust(pcombat,method="BH")	#
	#table(qcombat<0.05)
	


ret = getdifftab_limma(combat_dat, sampleannotation$group, "Verified-Generic" )

table(ret$padjusted<0.05)



ret[ret$fdr<0.05,]

	
	
#getdifftab_limma = function(edata, sa, threshold=0.05)
getdifftab_limma = function(edata, condition, contrast, threshold=0.05)
{
	require(limma)
	fac = as.factor(condition)	
	design = model.matrix(~0 + fac)
	colnames(design)=levels(fac)
	#contrast = paste(unique(sa$treatment)[1], "-", unique(sa$treatment)[2], sep="")
	cont.matrix = makeContrasts ( contrasts=contrast, levels=design)	
	fit <- lmFit(edata, design)
	fit2 = contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2)
	topTable(fit2)
	ret = data.frame(gene=names(fit2$Amean), p=fit2$p.value[,1], fdr=p.adjust(fit2$p.value[,1] , method="fdr"), fc=fit2$coefficients[,1])
	return(ret)
}

ret = getdifftab_limma(combat_dat, sampleannotation$group, "Verified-Generic" )

table(ret$fdr<0.05)
ret[ret$fdr<0.05,]



ret = getdifftab_limma(datamatrix, sampleannotation$group, "Verified-Generic" )

table(ret$fdr<0.05)


#thisprcomp=prcomp( (t(na.omit(combat_dat))))
thisprcomp=prcomp( (t(na.omit(datamatrix))))
plotsize = 1600
name = "pcaplot_unnorm"
thisannot = as.factor(sampleannotation$group)
colorPal = c("green", "blue", "yellow", "cyan", "black", "red")
thiscolors = colorPal[as.numeric(thisannot)]
annotlabels=paste(unique(thisannot),table(thisannot)[unique(thisannot)])
annotlabelscolors=unique((thiscolors))
png(file = paste("pcaplot.png", sep=""),pointsize = 12, width = plotsize, height = plotsize, units="px")
plot(thisprcomp$x,  col=thiscolors, cex=2 , lwd=5, main=name, cex.main=4, bg=thiscolors)
legend("bottomleft" , legend=annotlabels, text.col=annotlabelscolors, cex=3, bty="n")
dev.off()


plotonepca = function(name="PCAplot", thisannot, thisplatform, thisprcomp, plotdirname, thiscolors=NA)
{
	colorPal = c(brewer.pal(3,"Set1"), brewer.pal(8,"Set2"))
	plotsize = 1600
	goodpch = c(19,17,15,23,25)
	#thisannot =as.character(thisannot)
	#thisannot = factor(thisannot, levels=c("normal", "benign", "DCIS", "tumor", "mixed"))
	#annotcount = table(thisannot)
	#annotcount = annotcount[unique(thisannot)]
	
	if(!is.factor(thisannot))
		thisannot=factor(thisannot)
	
	#print(thisannot)
	if(is.na(thiscolors))
		thiscolors = colorPal[as.numeric(thisannot)]

	annotlabels=paste(unique(thisannot),table(thisannot)[unique(thisannot)])
	annotlabelscolors=unique((thiscolors))
	a=match(levels(thisannot), unique(thisannot))
	a=a[!is.na(a)]
	annotlabels=annotlabels[a]
	annotlabelscolors=annotlabelscolors[a]
		
	thispch = factor(thisplatform)
	levels(thispch) =  goodpch
	png(file = paste(plotdirname , "/pca_",name, ".png", sep=""),pointsize = 12, width = plotsize, height = plotsize, units="px")
	
	#plot(thisprcomp$x,  col=colorPal[as.factor(thisannot)], cex=2 , lwd=5, main=name, pch=as.numeric(as.character(thispch)) ,cex.main=4)
	#legend("bottomleft" , legend=paste(unique(thisannot),table(thisannot)[unique(thisannot)]), text.col=colorPal[unique(as.factor(thisannot))], cex=3, bty="n")
	plot(thisprcomp$x,  col=thiscolors, cex=2 , lwd=5, main=name, pch=as.numeric(as.character(thispch)) ,cex.main=4, bg=thiscolors)
	legend("bottomleft" , legend=annotlabels, text.col=annotlabelscolors, cex=3, bty="n")
	legend("bottomright" , pch=unique(as.numeric(as.character(thispch))), legend=paste("",unique(thisplatform)), text.col="black", cex=3, bty="n")	
	dev.off()

}


thisgeoseries=getGEO("GSE24759")
thiseset = thisgeoseries[[1]]
pheno = pData(thiseset)

t(pheno[1,])

#GSM609634


table(pheno[ , c("characteristics_ch1", "characteristics_ch1.1")])

colSums(table(pheno[ , c("characteristics_ch1", "characteristics_ch1.1")]))

rowSums(table(pheno[ , c("characteristics_ch1", "characteristics_ch1.1")]))

sampleannotation$description = pheno[,"description.1"][match(sampleannotation$Code,as.character(pheno$title) )  ]



pheno = pData(thiseset)

sampleannotation=pheno[, c( "title", "characteristics_ch1", "characteristics_ch1.1")]

names(sampleannotation)=c("geoaccession", "cell", "batch")

sampleannotation$cell = as.character(sampleannotation$cell)
sampleannotation$cell = gsub("\\+", "pluss", sampleannotation$cell)
sampleannotation$cell = gsub("\\-", "minus", sampleannotation$cell)
sampleannotation$cell = gsub(" ", "", sampleannotation$cell)
sampleannotation$cell = gsub("celltype:", "", sampleannotation$cell)
sampleannotation$cell = make.names(sampleannotation$cell)
table(sampleannotation[,2:3])


datamat= exprs(thiseset)
datamat[,] = rnorm(length(datamat))


#datamat= normalizeBetweenArrays(datamat)

mod = model.matrix(~as.factor(sampleannotation$cell))
mod0 = model.matrix(~1,data=sampleannotation)

date()
combatdata= as.matrix(ComBat(dat=datamat, batch=sampleannotation$batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE))
date()



cella="Megakaryocyte.erythroidprogenitor"
cellb="CD4plussCentralMemory"
beforecombat = getdifftab_limma(datamat, sampleannotation$cell, paste(cella, "-", cellb, sep="") )
table(beforecombat$padjusted<0.1)


aftercombat = getdifftab_limma(combatdata, sampleannotation$cell, paste(cella, "-", cellb, sep="") )
table(aftercombat$padjusted<0.1)


combatdata= as.matrix(ComBat(dat=combatdata, batch=sampleannotation$batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE))
aftercombat = getdifftab_limma(combatdata, sampleannotation$cell, paste(cella, "-", cellb, sep="") )
table(aftercombat$padjusted<0.5)

	



"cell type: Basophils"
"cell type: CD4+ Central Memory"

description.1
t([1,])












t(pheno[2,])


sampleannotation = read.table("data/table_S1.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#sampleannotation$inmatrix[sampleannotation$code %in%as.character(pheno$title)]
sampleannotation$group = pheno[,"characteristics_ch1"][match(sampleannotation$code,as.character(pheno$title) )  ]
sampleannotation$geoaccession = pheno[,"geo_accession"][match(sampleannotation$code,as.character(pheno$title) )  ]
sampleannotation$group = gsub("group: ", "",sampleannotation$group)
sampleannotation$group = gsub("Deliberately ", "",sampleannotation$group)
sampleannotation$group = gsub(" GA", "",sampleannotation$group)
sampleannotation$group = gsub("Unverified ", "",sampleannotation$group)
sampleannotation$group = gsub("Verified ", "",sampleannotation$group)

write.table(sampleannotation, "data/sampleannotation.csv", sep="\t", row.names=FALSE, quote=FALSE)


x = sampleannotation[,c(2,6)]
y = pheno[, c("title", "geo_accession")]

x[order(x[,1]),][1:10,]
y[order(y[,1]),][1:10,]




