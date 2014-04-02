
# source("figuretowfic.r")
starttime = Sys.time()
debug = FALSE
downloaddata=FALSE
set.seed(100)


includelibs = c("Biobase", "GEOquery", "sva", "limma",
                "preprocessCore")
lapply(includelibs, require, character.only=T)


sampleannotation = read.table("data/sampleannotation.csv", sep="\t", header=TRUE,  stringsAsFactors=FALSE)
sampleannotation$code = make.names(sampleannotation$code)
sampleannotation$chip = as.character(sampleannotation$chip)
dimnames(sampleannotation)[[1]] = sampleannotation$code
# take out 3 samples that are not assign to a geoaccession. Failed QC?
sampleannotation = sampleannotation[!is.na(sampleannotation$geoaccession),] 


source("../../commonscripts/helper_functions.r")


ret = loadtowfic(downloaddata)
datamatrices = list()
datamatrices[["real_raw"]] = ret[["data"]]
probeannotation = ret[["probes"]]

if(debug)
  datamatrices[["real_raw"]] = datamatrices[["real_raw"]][1:1000,]


datamatrices[["real_qnorm"]] = normalize.quantiles(datamatrices[["real_raw"]]) 
dimnames(datamatrices[["real_qnorm"]]) = dimnames(datamatrices[["real_raw"]])



# combat adjusted
combatmod = model.matrix(~as.factor(sampleannotation$covariate))
datamatrices[["real_combat_covariates"]]= as.matrix(ComBat(dat=datamatrices[["real_qnorm"]],
                                                           batch=sampleannotation$chip,
                                                           mod=combatmod,
                                                           numCovs=NULL,
                                                           par.prior=TRUE,
                                                           prior.plots=FALSE))

group = factor(sampleannotation$covariate)
design = model.matrix(~0 + group)
fit = lmFit(datamatrices[["real_combat_covariates"]], design)
cont.matrix = makeContrasts ( contrasts="groupDP-groupN", levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
limma_ret_alt = eBayes(fit2)
limma_p_alt = limma_ret_alt$p.value
print(paste("ComBat adjusted real data, significant probes: ",  sum(p.adjust(limma_p_alt, method="fdr")<0.05)))


#Limma adjusted
group = factor(sampleannotation$covariate)
block = as.factor(sampleannotation$chip)
design = model.matrix(~0+group+block)
fit = lmFit(datamatrices[["real_qnorm"]], design)
cont.matrix = makeContrasts ( contrasts="groupDP-groupN", levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
limma_ret_woc = eBayes(fit2)
limma_p_woc = limma_ret_woc$p.value[,1]
print(paste("Limma adjusted real data, significant probes: ",  sum(p.adjust(limma_p_woc, method="fdr")<0.05)))



# radnom, ComBat adjusted
set.seed(100)
datamatrices[["random_raw"]] = matrix(rnorm(length(datamatrices[["real_raw"]]), mean=0, sd=1), 
                                      nrow=nrow(datamatrices[["real_raw"]]), 
                                      ncol=ncol(datamatrices[["real_raw"]]), 
                                      dimnames=dimnames(datamatrices[["real_raw"]]))
datamatrices[["random_combat_covariates"]]= as.matrix(ComBat(dat=datamatrices[["random_raw"]],
                                                             batch=sampleannotation$chip, 
                                                             mod=combatmod, 
                                                             numCovs=NULL, 
                                                             par.prior=TRUE, 
                                                             prior.plots=FALSE))
group = factor(sampleannotation$covariate)
design = model.matrix(~0 + group)
cont.matrix = makeContrasts ( contrasts="groupDP-groupN", levels=design)
fit = lmFit(datamatrices[["random_combat_covariates"]], design)  
fit2 = contrasts.fit(fit, cont.matrix)
limma_p_rand_combat = eBayes(fit2)$p.value[,1]
print(paste("ComBat adjusted random data, significant probes: ",  sum(p.adjust(limma_p_rand_combat, method="fdr")<0.05)))






pdf(file="figure/towficpvalues.pdf")
thiscolors = c("red", "blue", "black")
xrange=1:20
xrange=1:25
# reproduced ComBat + limma
plot((xrange)/100, hist(limma_p_alt, breaks=100, plot=FALSE)$counts[xrange], 
     main="P-values",
     xlab="p-value", 
     ylab="frequency",
     type="l",
     lwd=2,
     col=thiscolors[1])

# batch handled in limma
lines((xrange)/100, hist(limma_p_woc, breaks=100, plot=FALSE)$counts[xrange], 
      col=thiscolors[2], lwd=2)

# random ComBat + limma
lines((xrange)/100, hist(limma_p_rand_combat, breaks=100, plot=F)$counts[xrange],
      col=thiscolors[3], lwd=2)

legend("topright", legend=c("Real data, ComBat adjusted",
                            "Real data, batch handled by Limma",
                            "Random data, ComBat adjusted"),
       text.col=thiscolors)

dev.off()





print(paste( "Figure generated for Towfic et al data set. Time spent ", 
             as.integer(round(difftime(Sys.time(),starttime, units="mins"))  ), 
             " minutes", sep="") )

