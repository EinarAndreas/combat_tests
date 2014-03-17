
========================================================

Show that adjusting for batch without trying to retain group differeces also can induce group effects.



```r
library(sva)
```

```
## Loading required package: corpcor
## Loading required package: mgcv
## Loading required package: nlme
## This is mgcv 1.7-28. For overview type 'help("mgcv-package")'.
```

```r
source("../../commonscripts/helper_functions.r")
```


Creating a unbalances experiment, 3 conditions, 2 batches with some genes true different in condition 3 compared to 1 and 2. But 1 and 2 have no difference. No batch effect added but adjusted for anyways.

```r
n = 10
ngenes = 5000
ndiffgenes = round(ngenes/10)
sampleannotation = createsampleannotation(list(c(n, n, 0), c(n, 0, n)))
# sampleannotation = createsampleannotation( list(c(n,1), c(1,1), c(1,n)))
randomdata = createrandomdata(ngenes, sampleannotation, mean = 0, sd = 1)
conditiondata = addconditioneffect(df = randomdata, labels = sampleannotation$treatment, 
    ndiffgenes = ndiffgenes, thismean = 10, betweensd = 1, insidesd = 1, affectedconditions = "treatment3")
print(table(sampleannotation[, c("batch", "treatment")]))
```

```
##      treatment
## batch treatment1 treatment2 treatment3
##     1         10         10          0
##     2         10          0         10
```


Then run ComBat without covariates.

```r
batcnnormdata = ComBat(dat = conditiondata, batch = sampleannotation$batch, 
    mod = NULL, par.prior = TRUE, prior.plots = FALSE)
```

```
## Found 2 batches
## Found 0  categorical covariate(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
```

```r

```


And plot p-values

```r
ret=getdifftab_limma(batcnnormdata, sampleannotation$treatment, "treatment1-treatment2")
```

```
## Loading required package: limma
```

```r
#table(ret$padjusted<0.05)
hist(ret$p, breaks=100, border="red")

ret=getdifftab_limma(conditiondata, sampleannotation$treatment, "treatment1-treatment2", block=sampleannotation$batch)
#table(ret$padjusted<0.05)
hist(ret$p, breaks=100, add=T, border="black")

a = sampleannotation$batch==1
ret=getdifftab_limma(conditiondata[,a], sampleannotation$treatment[a], "treatment1-treatment2")
#table(ret$padjusted<0.05)
hist(ret$p, breaks=100, add=T, border="blue")
legend("topright", legend=c("ComBat adjusted", "Limma adjusted", "Only batch 1, no adjustment"), text.col=c("red", "black", "blue"))
```

![plot of chunk pvaluespermuted](figure/pvaluespermuted.svg) 

```r

```

Treatment1 and Treatment2 differs. Some of the effect in Treatment3 has been transfered overt to Treatment1 since they both are in batch2.
This problem is not limited to ComBat.






