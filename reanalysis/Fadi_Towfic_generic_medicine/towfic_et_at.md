Reproduction and alternative analysis of "Comparing the Biological Impact of Glatiramer Acetate with the Biological Impact of a Generic"
========================================================




2014-02-24 08:34:02

### Overview
This report aims to show that the use of the statistical tool ComBat from [Johnson et al.](http://biostatistics.oxfordjournals.org/content/8/1/118.abstract) led to false results in the analysis performed in [Towfic et al.'s ](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0083757) "Comparing the Biological Impact of Glatiramer Acetate with the Biological Impact of a Generic".

This document has five main parts
- Getting and formatting the data
- Reproduce some of the results to show that we are working on the same data and analysis workflow
- Remove the use of ComBat and perform an similar analysis with alternative established tools
- Estimate the error introduced by ComBat and the consequences for the conclusion of the study
- Perform a few more sanity checks to substantiate that the difference in results for the two above analysis is mainly due to error introduced by ComBat


### Getting the data



```r
includelibs = c("Biobase", "GEOquery", "sva", "limma", "preprocessCore")
lapply(includelibs, require, character.only = T)
```

The important sample information is described in Table S1 from Towfic et al. the its usage in ComBat is described briefly in the methods;
>  Each microarrays chip designation was supplied a batch label; there were 18 batches in all. The labels for the treatments (i.e. drug product, reference standard) were added as covariates.<cite> Towfic et al.

Table S1 does have a "Chip"- coloumn, unfortunatly there is no dedicated "treatment"-coloumn.
Communication with the corresponding author yielded this explanation
> ..and the only thing we used as a covariate was the unique treatment names themselves (e.g. ³Medium² or ³RS").

Based on these descriptions and the annotation for the GEO deposit, we compiled a more complete sample annotation file.


```r
sampleannotation = read.table("data/sampleannotation.csv", sep = "\t", header = TRUE, 
    stringsAsFactors = FALSE)
sampleannotation$code = make.names(sampleannotation$code)
sampleannotation$chip = as.character(sampleannotation$chip)
dimnames(sampleannotation)[[1]] = sampleannotation$code
head(sampleannotation)
```

```
##       number  code covariate       chip slot geogroup geoaccession
## RS.16      1 RS.16        RS 4634633002    A Verified    GSM996691
## DP.34      2 DP.34        DP 4634633002    B Verified    GSM996731
## DP.19      3 DP.19        DP 4634633002    C Verified    GSM996716
## DP.32      4 DP.32        DP 4634633002    D Verified    GSM996729
## DP.04      5 DP.04        DP 4634633002    E Verified    GSM996701
## MBP        6   MBP       MBP 4634633002    F    other    GSM996658
```

The "covariate"-coloumn is made based on the "code"-coloumn, but might not match  to what were actually used. The last two coloumns are from the GEO deposit and reveal that 3 of the samples do not have data and are taken out.


```r
# take out 3 samples that are not assign to a geoaccession. Failed QC?
sampleannotation = sampleannotation[!is.na(sampleannotation$geoaccession), ]
```

The batch/covariate design shows many batches and covariate groups.

```r
# the batch/covariate design
table(sampleannotation[, c("chip", "covariate")])
```

```
##             covariate
## chip         A C DM DP E GA M MANITOL MBP N RS S TV
##   4634633002 0 0  0  4 0  0 0       0   1 0  1 0  0
##   4634633043 0 0  1  0 0  1 1       0   0 2  1 0  0
##   4634633053 0 0  1  1 0  0 0       0   0 1  1 0  0
##   4637105025 0 0  0  5 0  0 0       1   0 0  0 0  0
##   4682416019 0 0  0  4 0  0 0       0   0 0  2 0  0
##   4682416043 0 0  0  1 0  2 0       0   0 1  2 0  0
##   4682416087 0 0  1  0 0  0 1       0   0 0  1 0  3
##   4682416090 0 0  0  4 0  0 1       0   0 1  0 0  0
##   4763128059 0 0  1  2 0  0 1       0   0 1  1 0  0
##   4763128060 0 0  1  2 0  0 1       0   0 1  1 0  0
##   4763646004 0 0  0  1 0  0 1       0   0 0  4 0  0
##   4763646012 0 1  0  1 0  0 1       0   0 2  1 0  0
##   4763646026 0 1  0  1 0  0 1       0   0 2  1 0  0
##   4763646030 0 0  0  2 0  0 2       0   0 0  2 0  0
##   5216898004 1 0  0  2 0  0 1       0   0 0  1 0  0
##   5216898008 0 0  0  2 0  0 1       0   0 0  1 1  1
##   5216898012 2 0  0  1 1  0 1       0   0 0  1 0  0
##   5216898024 0 0  0  1 1  0 1       0   0 0  1 1  1
```

A look at the most interesting groups, "DP"(drug product) and "N"(generic) reveals a lack of balance.

```r
table(sampleannotation[sampleannotation$covariate %in% c("DP", "N"), c("chip", 
    "covariate")])
```

```
##             covariate
## chip         DP N
##   4634633002  4 0
##   4634633043  0 2
##   4634633053  1 1
##   4637105025  5 0
##   4682416019  4 0
##   4682416043  1 1
##   4682416090  4 1
##   4763128059  2 1
##   4763128060  2 1
##   4763646004  1 0
##   4763646012  1 2
##   4763646026  1 2
##   4763646030  2 0
##   5216898004  2 0
##   5216898008  2 0
##   5216898012  1 0
##   5216898024  1 0
```


> The microarray data have been deposited in the Gene Expression Omnibus, under accession number [GSE40566](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40566)

The processed data from the above GEO deposit was deposited for another publication which did not processes the data as described in Tawfic et al. (personal communication). The raw data linked to as [GSE40566_RAW.tar](http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40566&format=file) and [GSE40566_non_normalized.txt.gz](http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40566&format=file&file=GSE40566%5Fnon%5Fnormalized%2Etxt%2Egz) were thus used, although it is not confirmed that this is the version of the data actually used in Towfic et al.

```r
# the output from Illumina software.
rawdata = read.table("data/GSE40566_non_normalized.txt", sep = "\t", header = TRUE, 
    stringsAsFactors = FALSE)

# the probe annotation file found inside GSE40566_RAW.tar
rawannotation = read.table("data/GPL6887_MouseWG-6_V2_0_R3_11278593_A.txt", 
    sep = "\t", header = TRUE, stringsAsFactors = FALSE, skip = 8, comment.char = "", 
    quote = "", fill = TRUE)
```


Then some formatting is needed to match the probe information to the measurments and the measurments to the sample information.

```r
# annotation file has two separate tables, one for experimental and one at
# the end for controls.  splitting and fixing the tables:
tmp = rawannotation$Species == "Mus musculus"
experimentalannot = rawannotation[tmp, ]
experimentalannot$Array_Address_Id = as.numeric(experimentalannot$Array_Address_Id)
controlannot = rawannotation[!tmp, ]
dimnames(controlannot)[[2]] = rawannotation[rawannotation[, 2] == "Array_Address_Id", 
    ]
controlannot$Array_Address_Id = suppressWarnings(as.numeric(controlannot$Array_Address_Id))
controlannot = controlannot[!is.na(controlannot$Array_Address_Id), ]
controlannot = controlannot[, 1:6]
probeannotation = merge(experimentalannot, controlannot, all = TRUE)
# dim(probeannotation)
rm(tmp, experimentalannot, controlannot)

probeannotation = probeannotation[!duplicated(probeannotation$Array_Address_Id), 
    ]
probeannotation = probeannotation[probeannotation$Array_Address_Id %in% rawdata$ID_REF, 
    ]  # 
probeannotation$Symbol = tolower(probeannotation$Symbol)
dimnames(probeannotation)[[1]] = probeannotation$Probe_Id
# dim(probeannotation)

# sort and filter probe and data similar.
datamatrix_raw = as.matrix(rawdata[, -1])
datamatrix_raw = datamatrix_raw[match(probeannotation$Array_Address_Id, rawdata$ID_REF), 
    ]
dimnames(datamatrix_raw)[[1]] = probeannotation$Probe_Id
# dim(datamatrix_raw) dim(probeannotation)

# and match data to samples. table(sampleannotation$code %in%
# dimnames(datamatrix_raw)[[2]])# check table(dimnames(datamatrix_raw)[[2]]
# %in% sampleannotation$code)# check
datamatrix_raw = datamatrix_raw[, match(sampleannotation$code, dimnames(datamatrix_raw)[[2]])]
```


Several tables with results are presented in the Supporting Information of in Tawfic et al. We are aiming to reproduce these. Unfortunatly these were only presented in a pdf-format not easely parsed by a computer. We had to resort to a ad-hoc method of cut-and-paste from pdf into text files which were somewhat polluted by pdf-formatting code. For some tables about 10% of the rows will be lost, but the rest will suffice.


```r

table_s2 = read.table("data/table_s2.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE, 
    fill = TRUE)
table_s2 = as.matrix(table_s2)
table_s2[, 3:6] = as.numeric(table_s2[, 3:6])
```

```
## Warning: NAs introduced by coercion
```

```r
# not able to paste the pdf whitout a lot of gibberish clutter the data.
# Some probes are lost!
a = rowSums(is.na(table_s2[, 3:6])) > 0
print(paste("Lost rows table s2: ", sum(a)))
```

```
## [1] "Lost rows table s2:  91"
```

```r
table_s2 = data.frame(table_s2[!a, ], stringsAsFactors = FALSE)
# data.frame makes the columns characters again.
for (n in 3:6) {
    table_s2[, n] = as.numeric(table_s2[, n])
}


table_s3 = read.table("data/table_s3.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE, 
    fill = TRUE)
# one Probe, ('1770446' or 'ILM_1770446') was not found in the data and was
# thus removed from this textfile

table_s4 = read.table("data/table_s4.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE, 
    fill = TRUE)
# not able to copy-paste the pdf whitout a lot of gibberish clutter the
# data. Some probes are lost!
a = is.na(as.numeric(table_s4$Var_generic_vs_Var_GA))
```

```
## Warning: NAs introduced by coercion
```

```r
print(paste("Lost rows table s4: ", sum(a)))
```

```
## [1] "Lost rows table s4:  1"
```

```r
table_s4 = table_s4[!a, ]
table_s4$Var_generic_vs_Var_GA = as.numeric(table_s4$Var_generic_vs_Var_GA)

table_s5 = read.table("data/table_s5.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE, 
    fill = TRUE)
tmp = as.matrix(table_s5)
tmp[, 3:14] = as.numeric(tmp[, 3:14])
```

```
## Warning: NAs introduced by coercion
```

```r
# not able to paste the pdf whitout a lot of gibberish clutter the data.
# Some probes are lost!
a = rowSums(is.na(tmp[, 3:14])) > 0
print(paste("Lost rows table s5: ", sum(a)))
```

```
## [1] "Lost rows table s5:  15"
```

```r
table_s5 = data.frame(table_s5[!a, ], stringsAsFactors = FALSE)
table_s5$Fold_Change = as.numeric(table_s5$Fold_Change)
```



### Reproduce the original results


```r
set.seed(100)
datamatrices = list()
datamatrices[["real_raw"]] = datamatrix_raw[1:5000, ]

# used in paper, but seems to do the same as the limma version. Except it
# loses dimnames
datamatrices[["real_qnorm"]] = normalize.quantiles(datamatrices[["real_raw"]])
dimnames(datamatrices[["real_qnorm"]]) = dimnames(datamatrices[["real_raw"]])

mod = model.matrix(~as.factor(sampleannotation$covariate))
mod0 = model.matrix(~1, data = sampleannotation)

datamatrices[["real_combat_covariates"]] = as.matrix(ComBat(dat = datamatrices[["real_qnorm"]], 
    batch = sampleannotation$chip, mod = mod, numCovs = NULL, par.prior = TRUE, 
    prior.plots = FALSE))
```

```
## Found 18 batches
## Found 12  categorical covariate(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
```

The nameing convention seeme to differ slightly between the Table S1 and the text. This is the covariate labels with its corresponging name in the text.
- **DP** refferd to as GA (but not GA as in table S1)
- **N** refferd to as "generic"
- **M** reffered to as "medium"
- **RS** reffered to as "reference standard""

First we will try to reproduce Table S2.

> Genes utilized for the tolerance method illustrated in Figure 1B. (figure text)

> This standard of comparison was constructed by first identifying the top 1000 probes by absolute fold change of reference standard compared to the medium (Table S2). The list includes both upregulated and downregulated probes compared to medium. Probes were filtered such that ones upregulated by reference standard needed to have an average reference standard expression of 6.00 or higher and ones downregulated by reference standard needed to have an average medium expression of 6.00 or higher. 

This is how the table looks

```r
head(table_s2)
```

```
##             ID         Gene AVG_Medium AVG_Reference_Standard AVG_GA
## 1 ILMN_2685712         IFNG      7.152                 11.255 11.163
## 2 ILMN_1247309          IL3      5.823                  9.406  9.384
## 3 ILMN_2791459         IFNG      6.636                  9.829  9.711
## 4 ILMN_1215862        CXCL9      7.342                  9.817  9.757
## 5 ILMN_2595732 LOC100046232      9.249                 11.509 11.547
## 6 ILMN_2931334          IL4      6.172                  8.367  8.425
##   AVG_generic
## 1      10.668
## 2       9.422
## 3       9.304
## 4       9.525
## 5      11.767
## 6       8.496
```


**Table S2**

```r
thiscovariatename = "covariate"
thisdata = datamatrices[["real_combat_covariates"]]

# the mean uses log2 numbers and will be a geometric mean.
meanM = rowMeans(thisdata[, sampleannotation[, "covariate"] == "M"])
meanRS = rowMeans(thisdata[, sampleannotation[, "covariate"] == "RS"])
meanDP = rowMeans(thisdata[, sampleannotation[, "covariate"] == "DP"])
meanN = rowMeans(thisdata[, sampleannotation[, "covariate"] == "N"])

foldchange_RS_vs_M = meanRS - meanM
foldchange_RS_vs_M = foldchange_RS_vs_M[meanM > 6 | meanRS > 6]
top1000_RS_vs_M = names(foldchange_RS_vs_M[order(abs(foldchange_RS_vs_M), decreasing = TRUE)])[1:1000]
table(table_s2$ID %in% top1000_RS_vs_M)  #  851 out of 909. For some reason it is not a 100% match.
```

```
## 
## FALSE  TRUE 
##   824    85
```

```r

# notin =!(table_s2$ID %in% top1000_RS_vs_M) notinnames = table_s2$ID[notin]
# cbind(table_s2[notin,], meanRS[notinnames], meanM[notinnames],
# FC=meanRS[notinnames]-meanM[notinnames])
```


```r
# scatter plot correlation of the 4 columns in table_S2 vs. the reproduced numbers. String correlation.
par(mfrow=c(2, 2))
plot(table_s2$AVG_Medium, meanM[table_s2$ID], pch=".")
plot(table_s2$AVG_Reference_Standard, meanRS[table_s2$ID], pch=".")
plot(table_s2$AVG_GA, meanDP[table_s2$ID], pch=".")
plot(table_s2$AVG_generic, meanN[table_s2$ID], pch=".")
```

![plot of chunk reproduced_tables2](figure/reproduced_tables2.png) 



**Table S3**


```r
# Table S3. The highly variable probes that were significant by F-test in
# either GA or generic (see methods section) and are depicted in Figure 1A.
# MEthods: In order to identify probes with variability induced specifically
# by activation (as opposed to experimental noise), we sought to identify
# probes that were significantly more variable when activated with either GA
# or generic than medium. Using an F-test, we compared GA against Medium for
# each probe and compared generic against Medium. We then took the set of
# probes where either treatment comes up to be more variable than medium
# (union, passes in at least at least one). For those set of probes only, we
# compared the variability of GA across 34 samples representing 30 batches,
# to the variability of generic across 11 samples representing 5 batches,
# utilizing an F-test to measure significance of the differences between the
# probes.

# do a f-test between DP and medium, and N and medium. And see what is
# 'significant' (significance cut-off not noted in the text)?
pvalsN = apply(datamatrices[["real_combat_covariates"]], 1, function(x) {
    var.test(x[sampleannotation$covariate == "N"], x[sampleannotation$covariate == 
        "M"], alternative = "greater")$p.value
})
pvalsDP = apply(datamatrices[["real_combat_covariates"]], 1, function(x) {
    var.test(x[sampleannotation$covariate == "DP"], x[sampleannotation$covariate == 
        "M"], alternative = "greater")$p.value
})
padjustedN = p.adjust(pvalsN, method = "BH")
padjustedDP = p.adjust(pvalsDP, method = "BH")
print(table(padjustedN < 0.05))
```

```
## 
## FALSE 
##  5000
```

```r
print(table(padjustedDP < 0.05))
```

```
## 
## FALSE 
##  5000
```

```r
# nothing significant acording to 5%FDR that was used in the Differentially
# expressed genes tests.
```



```r
# look at the pvalues for the probes found in table_s3
thiscolors = c("black", "brown", "red", "blue")
plot(density(pvalsN[table_s3$Probe[table_s3$Type=="generic"]]), col=thiscolors[3], main="P-values for F-test, N or DP vs. Medium")
```

```
## Error: 'x' contains missing values
```

```r
lines(density(pvalsN), col=thiscolors[1])
```

```
## Error: plot.new has not been called yet
```

```r
lines(density(pvalsDP), col=thiscolors[2])
```

```
## Error: plot.new has not been called yet
```

```r
lines(density(pvalsDP[table_s3$Probe[table_s3$Type=="GA"]]), col=thiscolors[4])
```

```
## Error: 'x' contains missing values
```

```r
legend("topright", legend=c("All probes, N vs. M", "All probes, DP vs. M", "Table S3 probes, N vs. M", "Table S3 probes, DP vs. M"), text.col=thiscolors)
```

```
## Error: plot.new has not been called yet
```

```r
# The PRobes from table_s3 are at least enriched for low p-values.
# But table_s3 is not fully reproduced. More detailed explanation is needed
```




**Table S4**


```r
# Table s4 'Ranked list of probes by ratio of the variance in
# generic-activated samples to the variance in GA-activated samples.'
# 'Thus, we calculated for each probe the ratio of the variance in generic
# to the variance in GA'

# do a f-test between N and DP. Correlate with table s4.  f-test statistic
# is the same as ratio of variance?
fstatistics = apply(datamatrices[["real_combat_covariates"]], 1, function(x) {
    var.test(x[sampleannotation$covariate == "N"], x[sampleannotation$covariate == 
        "DP"], alternative = "greater")$statistic
})
fstatistics = fstatistics[order(fstatistics, decreasing = TRUE)]
table(names(fstatistics)[1:100] %in% table_s4$Probe)  # not a complete reproduction.
```

```
## 
## FALSE  TRUE 
##    95     5
```

```r
# look at the f-statistic for the probes found in table_s3 vs all
# f-statistics
```



```r
par(mfrow=c(2, 1))
plot(table_s4[,3], fstatistics[table_s4$Probe], main="F-statitics from table S4 vs. reproduced")
thiscolors = c("black", "red")
plot(density(fstatistics), col=thiscolors[1], main="Reproduced F-statistics for N vs. DP")
lines(density(fstatistics[table_s4$Probe]), col=thiscolors[2])
```

```
## Error: 'x' contains missing values
```

```r
legend("topright", legend=c("All probes", "Table S4 probes"), text.col=thiscolors)
```

![plot of chunk reproduced_table_s4](figure/reproduced_table_s4.png) 


**Table S5**




```r
# Table s5
#"Comparison of expression in GA to expression in GA for each probe, including fold change, ANOVA, LIMMA with background subtraction, comparative marker selection by signal-to-noise ratio, comparative marker selection by t-test, and the Wilcoxon non-parametric method."
#"Finally, genes that have significantly higher expression in samples activated by generic than in samples activated by GA by 4 different parametric methods (Table S5)"
#"Probes with Fold Change <  1.10 notdisplayed"

# Table s5 seems to be a comparison of expression in  GA (drug product) vs. generic (N) , listing the probes with greatest fold change (>abs(1.1)??) and the p-value calculated by several different tests. There seems to be no description of how the fold-change was calculated and after an inspection they look weird.

meanDP = rowMeans(thisdata[, sampleannotation$covariate=="DP"]) # DP is drug product also called GA
meanN = rowMeans(thisdata[, sampleannotation$covariate =="N"]) # N is Natco the generic,  (or G for alternative covariates)

# fold change between N and DP, calculated as geometric mean.
x = table_s5$Fold_Change
y = (meanDP - meanN)[table_s5$Probe]
lim =c( min(   c(min(x), min(y))  ), max( c(max(x), max(y)) ))
plot(x, y , ylim=lim, xlim=lim, pch=".", xlab="Table_S5", ylab="reproduced", main="FC comparison") # Not a perfect match. And the table_s5 numbers seems offset.
```

```
## Error: need finite 'xlim' values
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 


```r
# table S2 has the AVG for both DP and N listed
fc_table_s2 = table_s2$AVG_GA - table_s2$AVG_generic
names(fc_table_s2) = table_s2$ID
x = table_s5$Fold_Change
y = fc_table_s2[table_s5$Probe]
lim =c( min(   c(min(x, na.rm=T), min(y, na.rm=T))  ), max( c(max(x, na.rm=T), max(y, na.rm=T)) ))
plot(x, y , ylim=lim, xlim=lim, pch="." ,  xlab="Table_S5", ylab="Table_S2", main="FC comparison") # perfect match.
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.png) 

# But this is odd, we were able to reproduce table s2 and table s2 matches table s5, but table s5 fold changes does not match our reproduced values as good as the match for table s2.

# The reason that the line is not as straight for the reproduce FC's for table s5 as for the mean expressions in table S2 is because of different scale in the plot and the fact that the small difference will be higher for a ratio (combining to values with a individual deviation from the tables.)

```r
x = (meanDP - meanN)[table_s2$ID]
y = fc_table_s2
plot(x, y , ylim=lim, xlim=lim, pch="." ,  xlab="reproduced", ylab="Table_S2", main="FC comparison") # perfect match.
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17.png) 



```r
# Still there is a offset in the numbers.
commonprobes = table_s5$Probe[table_s5$Probe %in% table_s2$ID]
fc_table_s5 = table_s5$Fold_Change[match(commonprobes, table_s5$Probe)]
head(cbind(fc_table_s2[commonprobes], fc_table_s5, fc_table_s5 - fc_table_s2[commonprobes], 
    fc_table_s5/fc_table_s2[commonprobes]))
```

```
##                    fc_table_s5             
## ILMN_2685712 0.495       1.410 0.9150 2.848
## ILMN_2791459 0.407       1.330 0.9230 3.268
## ILMN_1220788 0.361       1.284 0.9233 3.558
## ILMN_2856926 0.351       1.276 0.9246 3.634
## ILMN_1218037 0.332       1.258 0.9263 3.790
## ILMN_2707941 0.308       1.238 0.9302 4.020
```

```r
# seems that some unacounted transformation has occured on the fold changes.
```


```r
source("../../commonscripts/helper_functions.r")
getdifftab_limma
```

```
## function (edata, condition, contrast, block = NULL, threshold = 0.05) 
## {
##     require(limma)
##     fac = as.factor(condition)
##     design = model.matrix(~0 + fac)
##     if (!is.null(block)) {
##         block = as.factor(block)
##         design <- model.matrix(~block + fac)
##     }
##     colnames(design) = make.names(colnames(design))
##     fit <- lmFit(edata, design)
##     contrast = paste("fac", contrast, sep = "")
##     contrast = sub("-", "-fac", contrast)
##     cont.matrix = makeContrasts(contrasts = contrast, levels = design)
##     fit2 = contrasts.fit(fit, cont.matrix)
##     fit2 <- eBayes(fit2)
##     ret = data.frame(gene = names(fit2$Amean), p = fit2$p.value[, 
##         1], padjusted = p.adjust(fit2$p.value[, 1], method = "fdr"), 
##         fc = fit2$coefficients[, 1], stringsAsFactors = FALSE)
##     dimnames(ret)[[1]] = ret$gene
##     return(ret)
## }
```


```r
# reproduce the Nom P-value columns
limma_ret = getdifftab_limma(edata = datamatrices[["real_combat_covariates"]], 
    condition = sampleannotation[, "covariate"], contrast = "DP-N")
limma_p = limma_ret$p
names(limma_p) = limma_ret$gene
DP = sampleannotation$covariate == "DP"
N = sampleannotation$covariate == "N"
ttest_p = apply(datamatrices[["real_combat_covariates"]], 1, function(x) {
    t.test(x[DP], x[N])$p.value
})
anova_p = apply(datamatrices[["real_combat_covariates"]], 1, function(x) {
    summary(aov(x[DP | N] ~ sampleannotation$covariate[DP | N]))[[1]][1, 5]
})  # aov er kanskje feil for dette, Er ikke live mange pr?ver i guppene.
wilcoxon_p = apply(datamatrices[["real_combat_covariates"]], 1, function(x) {
    wilcox.test(x[sampleannotation$covariate == "DP"], x[sampleannotation$covariate == 
        "N"])$p.value
})
```

```
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
## Warning: cannot compute exact p-value with ties
```


```r
par(mfrow=c(2, 2))
#plot(table_s5$Fold_Change, (meanDP - meanN)[table_s5$Probe] ,  ylim=c(-1.5,1.5), xlim=c(-1.5,1.5), pch=".", xlab="Table_S5", ylab="reproduced", main="FC comparison")
plot(log10(table_s5$t_test_p), log10(ttest_p[table_s5$Probe]), pch=".", main="t-test p")
plot(log10(table_s5$limma_p), log10(limma_p[table_s5$Probe]), pch=".", main="limma p")
plot(log10(table_s5$anova_p), log10(anova_p[table_s5$Probe]), pch=".", main="anova p")
plot(log10(table_s5$wilcoxon_p), log10(wilcoxon_p[table_s5$Probe]), pch=".", main="wilcoxon p")
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21.png) 

# The p-values does not make a straight line, but they are low. But maybe all are low?


```r
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
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22.png) 




