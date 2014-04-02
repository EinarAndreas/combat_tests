

starttime = Sys.time()

includelibs = c("sva", "qvalue", "limma")
tmp = lapply(includelibs, require, character.only=T)
print(tmp)
if(any(!unlist(tmp)))
{
  stop( paste("Not able to find all packages. Please install ",
              paste(includelibs[!unlist(tmp)], collapse=", ") )
  )
}
rm(tmp)
debug=FALSE
useparprior=TRUE


datamatrix = as.matrix(read.table("data/dataExample2.txt", sep="\t", header=TRUE))
if(debug)datamatrix = datamatrix[1:3000,]
sampleannotation = read.table("data/sampleInfoExample2.txt", 
                              sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(sampleannotation)=sampleannotation$ArrayName
sampleannotation$Batch=factor(as.character(sampleannotation$Batch)) 
datamatrix=datamatrix[,sampleannotation$Type!="WT"]
sampleannotation=sampleannotation[sampleannotation$Type!="WT",]
sampleannotation$Type=factor(sampleannotation$Type)

log2data = datamatrix
# flooring to 1
log2data[log2data<1]=1
# take out data with to much low/missing values.
negativeprobesfilter =( rowSums(log2data>1) >= (0.9*ncol(log2data)) )
log2data = log2data[negativeprobesfilter,]
# quantilenormalize
log2data=normalizeBetweenArrays(log2(log2data))

log2dataComBat = as.matrix( sva::ComBat(
  dat=log2data, 
  batch=sampleannotation[colnames(log2data),"Batch"], 
  mod=model.matrix(~as.factor(sampleannotation[colnames(log2data),"Type"])), 
  numCovs=NULL, 
  par.prior=useparprior, 
  prior.plots=FALSE))

# DE pvalues ComBat adjusted real data
Type = as.factor(sampleannotation$Type)
design = model.matrix(~0 + Type)
fit = lmFit(log2dataComBat, design)
cont.matrix = makeContrasts ( contrasts="TypeR-TypeC", levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
log2ComBat_res = eBayes(fit2)
rm(Type, design, fit, cont.matrix, fit2)
print(paste("ComBat adjusted real data, significant probes: ",  sum(qvalue(log2ComBat_res$p.value[,1])$qvalue<0.05)))

# DE pvalues not ComBat adjusted, 
# batch blocked in limma
Type = as.factor(sampleannotation$Type)
Block = as.factor(sampleannotation$Batch)
design = model.matrix(~0+Type+Block)
fit = lmFit(log2data, design)
cont.matrix = makeContrasts ( contrasts="TypeR-TypeC", levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
log2limmablock_res = eBayes(fit2)
rm(Type, Block, design, fit, cont.matrix, fit2)
print(paste("Limma batch adjusted real data, significant probes: ",  sum(qvalue(log2limmablock_res$p.value[,1])$qvalue<0.05)))



set.seed(100)
randomdata = log2data
randomdata[,] =rnorm(length(randomdata), mean=0, sd=1)
randomdataComBat = as.matrix( sva::ComBat(
  dat=randomdata, 
  batch=sampleannotation[colnames(randomdata),"Batch"], 
  mod=model.matrix(~as.factor(sampleannotation[colnames(randomdata),"Type"])), 
  numCovs=NULL, 
  par.prior=useparprior, 
  prior.plots=FALSE))

Type = as.factor(sampleannotation$Type)
design = model.matrix(~0 + Type)
cont.matrix = makeContrasts ( contrasts="TypeR-TypeC", levels=design)  

fit = lmFit(randomdataComBat, design)
fit2 = contrasts.fit(fit, cont.matrix)
randomdataComBat_pvalues = eBayes(fit2)$p.value[,1]
print(paste("ComBat adjusted random data, significant probes: ",  sum(qvalue(randomdataComBat_pvalues)$qvalue<0.05)))

pdf(file="figure/dataset2pvalues.pdf")

# reproduced ComBat + limma
hist(log2ComBat_res$p.value[,1], main="P-values, 3 settings", 
     breaks=100, xlab="p-value", ylim=c(0,3000), border="red")

# batch handled in limma
hist(log2limmablock_res$p.value[,1], add=T, breaks=100, xlab="p-value" , ylim=c(0,3000), border="blue")

# random ComBat + limma
hist(randomdataComBat_pvalues, add=T, breaks=100, xlab="p-value" , ylim=c(0,3000), border="black")

legend("topright", legend=c("Reald data ComBat adjusted", 
                            "Reald data Limma adjusted", 
                            "Random data ComBat adjusted"),
       text.col=c("red", "blue", "black"))

dev.off()



print(paste( "Figure generated for dataset2. Time spent ", 
             as.integer(round(difftime(Sys.time(),starttime, units="mins"))  ), 
             " minutes", sep="") )
      