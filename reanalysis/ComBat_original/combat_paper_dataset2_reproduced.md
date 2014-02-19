Reproduction and alternative analysis of data set 2 from Johnson et al
========================================================

Description how to obtain the data were found [here](https://groups.google.com/forum/#!msg/combat-user-forum/S9vBcXw8RGk/QIoD6oRBM9IJ) which had links to [dataExample2.txt](http://www.bu.edu/jlab/wp-assets/ComBat/data/dataExample2.txt)
and [sampleInfoExample2.txt](http://www.bu.edu/jlab/wp-assets/ComBat/data/sampleInfoExample2.txt)

The data consist of 35 samples of which 5 are annotated as "WT" and are not
refered to in the text or the plots. These are taken out before ComBat adjustment.

Read data and sample annotation
-------------------------------------------------------


```r
library(pheatmap)
library(RColorBrewer)
library(sva)
```

```
## Loading required package: corpcor
## Loading required package: mgcv
## Loading required package: nlme
## This is mgcv 1.7-28. For overview type 'help("mgcv-package")'.
```

```r
library(qvalue)
```


```r

orgdatamatrix = as.matrix(read.table("data/dataExample2.txt", sep = "\t", header = TRUE))
datamatrix = orgdatamatrix[1:3000, ]  # for dev.
datamatrix = orgdatamatrix  # for dev.
sampleannotation = read.table("data/sampleInfoExample2.txt", sep = "\t", header = TRUE, 
    stringsAsFactors = FALSE)
rownames(sampleannotation) = sampleannotation$ArrayName
sampleannotation$Batch = factor(as.character(sampleannotation$Batch))
sampleannotation$Cell = factor(as.character(sampleannotation$Cell))
# must be discrete for the pheatmap
# table(sampleannotation$ArrayName==dimnames(datamatrix)[[2]]) # ordercheck
datamatrix = datamatrix[, sampleannotation$Type != "WT"]
sampleannotation = sampleannotation[sampleannotation$Type != "WT", ]
sampleannotation$Type = factor(sampleannotation$Type)
print(dim(datamatrix))
```

```
## [1] 54675    30
```

```r
# dev/debug
useparprior = TRUE
```


Reproduce results
-------------------------------------------------------
First we try to reproduce the heatmap in Figure A.1, which uses 2698 genes with large variation.

```r
variatonmeasure = apply(datamatrix, 1, FUN= function(x){var(x)})
#variatonmeasure = apply(datamatrix, 1, FUN= function(x){sd(x)/mean(x)}) #CV
#clustermatrix = datamatrix[order(variatonmeasure, decreasing=TRUE),][1:50,]#[1:2698,]
clustermatrix = datamatrix[order(variatonmeasure, decreasing=TRUE),][1:2698,]
Batchcol = brewer.pal(8,"Set2")[1:3]
names(Batchcol) = levels(sampleannotation$Batch)
Typecol =brewer.pal(8,"Set2")[4:5]
names(Typecol) = levels(sampleannotation$Type)
ann_colors = list(Batch = Batchcol, Type = Typecol)
```




```r
pheatmap(clustermatrix, scale = "row", cluster_rows = T, cluster_cols = T, color = colorRampPalette(c("red", 
    "black", "green"))(n = 299), annotation = sampleannotation[, c("Batch", 
    "Type")], annotation_colors = ann_colors, cellheight = NA, cellwidth = 13, 
    fontsize = 12, drop_levels = TRUE, show_rownames = F, main = paste("Reproduced Fig A1.", 
        sep = ""), border_color = NA, treeheight_row = 0)
```

![plot of chunk reproduced_figA1](figure/reproduced_figA1.png) 


The Heatmap is not exactly as Fig. A1 in the paper. This could be due to different ways of calculating variation for a probe and many other ways the clustring could have been performed. But the purpose of the figure was to show that a batch effect is present in the data and espeically for batch 3. This is also achived in the reproduced figure.

Next 4 datasets were made,
- EB2: batches 1 and 2 adjusted with nonparametric use of ComBat
- EB3: All thre batches adjusted with nonparametric use of ComBat
- Batch12: batches 1 and 2 unadjusted
- Batch3: Only batch 3 (no adjustments)

```r
mat = datamatrix[,sampleannotation$Batch %in% c("1","2")]
EB2 = as.matrix( sva::ComBat(
            dat=mat, 
            batch=sampleannotation[colnames(mat),"Batch"], 
            mod=model.matrix( ~as.factor(sampleannotation[colnames(mat),"Type"])  ), 
            numCovs=NULL, 
            par.prior=useparprior, 
            prior.plots=FALSE))

mat = datamatrix[, ]
EB3 = as.matrix( sva::ComBat(
            dat=mat, 
            batch=sampleannotation[colnames(mat),"Batch"], 
            mod=model.matrix(~as.factor(sampleannotation[colnames(mat),"Type"])), 
            numCovs=NULL, 
            par.prior=useparprior, 
            prior.plots=FALSE))
Batch12 = datamatrix[,sampleannotation$Batch %in% c("1","2")]
Batch3 = datamatrix[,sampleannotation$Batch=="3"]
```


Next we try to re-create Figure A.2

> A heatmap diagram of 770 genes from data set 2 after applying the EB batch adjustments.<cite> Johnson et al

There is no description in how the 770 genes were selected so we use the same approach as in the reprodution of Figure A.1


```r
clustermatrix = EB3
variatonmeasure = apply(clustermatrix, 1, FUN = function(x) {
    var(x)
})
# variatonmeasure = apply(datamatrix, 1, FUN= function(x){sd(x)/mean(x)})
# #CV
clustermatrix = clustermatrix[order(variatonmeasure, decreasing = TRUE), ][1:770, 
    ]
Batchcol = brewer.pal(8, "Set2")[1:3]
names(Batchcol) = levels(sampleannotation$Batch)
Typecol = brewer.pal(8, "Set2")[4:6]
names(Typecol) = levels(sampleannotation$Type)
Cellcol = brewer.pal(8, "Set1")[1:5]
names(Cellcol) = levels(sampleannotation$Cell)
ann_colors = list(Batch = Batchcol, Type = Typecol, Cell = Cellcol)
pheatmap(clustermatrix, scale = "row", cluster_rows = T, cluster_cols = T, color = colorRampPalette(c("red", 
    "black", "green"))(n = 299), annotation = sampleannotation[, c("Batch", 
    "Type", "Cell")], cellheight = NA, cellwidth = 13, fontsize = 12, drop_levels = TRUE, 
    show_rownames = F, main = paste("Reproduced Fig A2.", sep = ""), border_color = NA, 
    treeheight_row = 0)
```

![plot of chunk reproduced_figA2](figure/reproduced_figA2.png) 




> Differential expression was assessed using Welchs t-test to determine the differential expression of RNAi versus control samples. EB2 produced at list of 86 significant genes at a false discovery (q-value) threshold of 0.05 (Storey and Tibshirani, 2003).<cite> Johnson et al


T-test for EB2 and count probes q<0.05.

```r
EB2_pvals = apply(EB2, 1 , 
                  FUN=function(x){t.test(
                                    x[sampleannotation[colnames(EB2), "Type"]=="C"],
                                    x[sampleannotation[colnames(EB2), "Type"]=="R"]
                                    )$p.value})
print(table(qvalue(EB2_pvals)$qvalue<0.05))
```

```
## 
## FALSE  TRUE 
## 54660    15
```




>The third batch alone produced a list of 37 significant genes using the same threshold. Crossing the significant gene lists, we observed 13 genes common in both lists (Fishers exact p-value < 1e-15). 


```r
Batch3_pvals = apply(Batch3, 1 , 
                  FUN=function(x){t.test(
                                    x[sampleannotation[colnames(Batch3), "Type"]=="C"],
                                    x[sampleannotation[colnames(Batch3), "Type"]=="R"]
                                    )$p.value})
print(table(qvalue(Batch3_pvals)$qvalue<0.05))
```

```
## 
## FALSE  TRUE 
## 54657    18
```




>Without any adjustment, combining these two batches produced a list of only 9 genes a q-value cutoff of 0.05


```r
Batch12_pvals = apply(Batch12, 1 , 
                  FUN=function(x){t.test(
                                    x[sampleannotation[colnames(Batch12), "Type"]=="C"],
                                    x[sampleannotation[colnames(Batch12), "Type"]=="R"]
                                    )$p.value})
print(table(qvalue(Batch12_pvals)$qvalue<0.05))
```

```
## 
## FALSE  TRUE 
## 54671     4
```



>Welchs t-test was also applied to EB3 to find differential expressed genes; yielding 1599 genes significant at a q-value cutoff of 0.05


```r
EB3_pvals = apply(EB3, 1 , 
                  FUN=function(x){t.test(
                                    x[sampleannotation[colnames(EB3), "Type"]=="C"],
                                    x[sampleannotation[colnames(EB3), "Type"]=="R"]
                                    )$p.value})
print(table(qvalue(EB3_pvals)$qvalue<0.05))
```

```
## 
## FALSE  TRUE 
## 53672  1003
```


>Reducing the q-value threshold to 0.01 yielded 488 significant genes

```r
print(table(qvalue(EB3_pvals)$qvalue < 0.01))
```

```
## 
## FALSE  TRUE 
## 54330   345
```


>..decreasing the threshold further to 0.001 yielded 161 significant genes.

```r
print(table(qvalue(EB3_pvals)$qvalue < 0.001))
```

```
## 
## FALSE  TRUE 
## 54572   103
```

The reproduced number of significant probes are not the same as the reported ones, but they are not completly off.


```r
require(limma)
```

```
## Loading required package: limma
```

```r
Type = as.factor(sampleannotation$Type)
design = model.matrix(~0 + Type)
colnames(design) = levels(Type)
fit <- lmFit(EB3, design)
cont.matrix = makeContrasts(contrasts = "R-C", levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
EB3_limma_pvalues = fit2$p.value[, 1]
print(table(qvalue(EB3_limma_pvalues)$qvalue < 0.05))
```

```
## 
## FALSE  TRUE 
## 53125  1550
```


Result when including batch in the statistical test
--------------


```r
Type = as.factor(sampleannotation$Type)
Block = as.factor(sampleannotation$Batch)
design <- model.matrix(~0 + Type + Block)
fit <- lmFit(datamatrix, design)
cont.matrix = makeContrasts(contrasts = "TypeR-TypeC", levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
datamatrix_limma_pvalues = fit2$p.value[, 1]
print(table(qvalue(datamatrix_limma_pvalues)$qvalue < 0.05))
```

```
## 
## FALSE  TRUE 
## 53905   770
```



```r
hist(EB3_limma_pvalues, border = "blue", main = "P-values, ComBat vs Limma", 
    breaks = 100, xlab = "p-value")
hist(datamatrix_limma_pvalues, border = "red", add = T, breaks = 100)
legend("topright", legend = c("ComBat adjusted", "Limma adjusted"), text.col = c("blue", 
    "red"))
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.svg) 

This plot shows that the p-values calculated from the ComBat adjusted data is skewed towards the low end compared to when batch is considered inside the statistical test in limma. However, the difference is not large and from this plot alone it could be argued that ComBat is just better.

Additional sanity checks
---------


```r
set.seed(100)
randommatrix = datamatrix
randommatrix[,] =rnorm(length(datamatrix), mean=0, sd=1)
EBrand = as.matrix( sva::ComBat(
            dat=randommatrix, 
            batch=sampleannotation[colnames(randommatrix),"Batch"], 
            mod=model.matrix(~as.factor(sampleannotation[colnames(randommatrix),"Type"])), 
            numCovs=NULL, 
            par.prior=useparprior, 
            prior.plots=FALSE))
```

```
## Found 3 batches
## Found 1  categorical covariate(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
```

```r

#source("../../helper_functions.r")
#randommatrixbatch = addbatcheffect(randommatrix, sampleannotation$Batch, thismean=0, thissd=0)# introduce batch effect
```


```r

Type = as.factor(sampleannotation$Type)
design = model.matrix(~0 + Type)
cont.matrix = makeContrasts ( contrasts="TypeR-TypeC", levels=design)  

fit <- lmFit(EBrand, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
EBrand_limma_pvalues = fit2$p.value[,1]

fit <- lmFit(randommatrix, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
randommatrix_limma_pvalues = fit2$p.value[,1]



Block = as.factor(sampleannotation$Batch)
design <- model.matrix(~0+Type+Block)

fit <- lmFit(randommatrix, design)
cont.matrix = makeContrasts ( contrasts="TypeR-TypeC", levels=design) 
fit2 = contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
randommatrix_limma_batch_pvalues = fit2$p.value[,1]

```



```r
hist(EBrand_limma_pvalues, border = "blue", main = "P-values, Random numbers", 
    breaks = 100, xlab = "p-value")
hist(randommatrix_limma_pvalues, border = "black", add = T, breaks = 100)
hist(randommatrix_limma_batch_pvalues, border = "red", add = T, breaks = 100)
legend("topright", legend = c("ComBat adjusted", "Limma adjusted", "No adjustment"), 
    text.col = c("blue", "black", "red"))
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.svg) 

The p-values for ComBat adjusted random numbers are enriched slighlty for low p-values. Indicating that the same enrichment seen for the real data is also false.

Another sanity check is to subset the real data into fictive batches. The EB2 data set contains the same number of samples with the same "Type" balance as the third batch alone. 
>The third batch was used for comparison against the EB2 analysis results because it was an identical experiment to EB2 other than the fact that it was conducted in a single batch.
Thus it is interesting to inspect the p-value from these two indentical experiments, EB2 (batch 1 and 2 adjusted with ComBat) and batch 3 (no ComBat adjustment). These p-values are computed above.

```r
hist(EB2_pvals, border = "blue", main = "P-values, real data", breaks = 100, 
    xlab = "p-value")
hist(Batch3_pvals, border = "red", add = T, breaks = 100)
legend("topright", legend = c("ComBat adjusted batch 1 and 2", "Batch 3"), text.col = c("blue", 
    "red"))
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17.svg) 

The 2-batch experiment is able to retrive more significant genes than running all in one batch. A further test of this is to split batch 3 into two imaginary batches in the same design as for batch 1 and 2, and then look at the p-value distribution compared with the real batch 3.


```r
Batch45 = Batch3
Batch45_annot = sampleannotation[dimnames(Batch45)[[2]],]
Batch45_annot$Batch="4"
Batch45_annot$Batch[Batch45_annot$Type=="C"][1:4] = "5"
Batch45_annot$Batch[Batch45_annot$Type=="R"][1:3] = "5"

table(sampleannotation[sampleannotation$Batch %in% c("1", "2"), c("Batch", "Type")])
```

```
##      Type
## Batch C R
##     1 2 6
##     2 4 3
##     3 0 0
```

```r
table(Batch45_annot[, c("Batch", "Type")])
```

```
##      Type
## Batch C R
##     4 2 6
##     5 4 3
```

```r

EB45 = as.matrix( sva::ComBat(
            dat=Batch45, 
            batch=Batch45_annot[colnames(Batch45),"Batch"], 
            mod=model.matrix( ~as.factor(Batch45_annot[colnames(Batch45),"Type"])  ), 
            numCovs=NULL, 
            par.prior=useparprior, 
            prior.plots=FALSE))
```

```
## Found 2 batches
## Found 1  categorical covariate(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
```

```r

EB45_pvals = apply(EB45, 1 , 
                  FUN=function(x){t.test(
                                    x[Batch45_annot[colnames(EB45), "Type"]=="C"],
                                    x[Batch45_annot[colnames(EB45), "Type"]=="R"]
                                    )$p.value})

hist(EB45_pvals, border="blue",  main="P-values, real data", breaks=100, xlab="p-value")
hist(Batch3_pvals, border="red", add=T, breaks=100)
legend("topright", legend=c("ComBat adjusted pseudo batches of batch", "Batch 3 not adjusted"), text.col=c("blue", "red"))
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18.svg) 



Done!


