#!/usr/bin/env Rscript

#########################################
## Methylation Pre-processing script ####
#########################################

library("optparse")
# get options

option_list = list(
  make_option(c("-f", "--folder"), type="character", default=".", help="folder with idat files [default= %default]", metavar="character"),
  make_option(c("-t", "--target"), type="character", default=".", help="target file with sample information [default= %default]", metavar="character"),
  make_option(c("-c", "--crossreac"), type="character", default=NULL, help="file with list of cross-reactive probes [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out", help="output directory name [default= %default]", metavar="character"),
  make_option(c("-s", "--snp_filter"), action="store_true", default=FALSE, type="logical", help="filter SNPs-associated probes [default= %default]", metavar="logical"),
  make_option(c("-m", "--multimodal_filter"), action="store_true", default=FALSE, type="logical", help="filter multimodal probes [default= %default]", metavar="logical")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

library(minfi)
library(ENmix)
library(wateRmelon)
require(MASS)
require(broom)

# load files
targets = read.table(opt$target,h=T,sep = "\t")
colnames(targets) = tolower(colnames(targets))
rownames(targets) <- targets$barcode
targets$Basename <- rownames(targets)
head(targets)
RGSet <- read.metharray.exp(base = opt$folder, targets = targets, recursive = T)
print(RGSet)

# write data
outdir    <- opt$out
dir.create(outdir, showWarnings = FALSE)
setwd(outdir)
save(RGSet, file="RGSet.RData")

# Quality control
dir.create("QC/", showWarnings = FALSE)
setwd("QC/")
plotCtrl(RGSet)
setwd("../")
qcReport(RGSet, pdf= "QC/qcReport.pdf")

MSet <- preprocessRaw(RGSet)
save(MSet, file="MSet.RData")
ratioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
gset <- mapToGenome(ratioSet, mergeManifest=T)
save(gset, file="gset.RData")

qc <- getQC(MSet)
rm(MSet)
svg("QC/Meth-unmeth_intensities.svg",h=4,w=4)
par(mfrow=c(1,1),family="Times",las=1)
plotQC(qc) # If U and/or M intensity log medians are <10.5, sample is by default of bad quality
dev.off()

#det p values
detP <- detectionP(RGSet)
pdf("QC/detection_pvalues.pdf",h=4,w=4)
par(mfrow=c(1,1),family="Times",las=1)
barplot(colMeans(detP), las=2, cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
dev.off()

if(sum(colMeans(detP)>0.01)>0 ) print(paste("Warning: sample", names(colMeans(detP))[colMeans(detP)>0.01], "has >0.01 mean detection p-value" ) )

# Surrogate variables and Principal components
# Surrogate variables derived from intensity data for non-negative internal control probes.
# These variables can be modeled in association analysis to adjust for experimental batch effects.
sva<-ctrlsva(RGSet, percvar = 0.9, flag = 1)
save(sva, file="sva90.RData")

K = ncol( sva )
lmsvaFull = lapply(1:K, function(i) lm(sva[,i]~Sentrix_id+Sentrix_position+Center,
                                       data.frame("Sentrix_id"=as.factor(pData(RGSet)$sentrix_id),"Sentrix_position"=as.factor(pData(RGSet)$sentrix_position),"Center"=as.factor(pData(RGSet)$center) ) ) )
lmsvaRed = vector("list",K)
for(i in 1:K){
  lmtmp = lmsvaFull[[i]]
  while(1){
    dttmp = dropterm(lmtmp,test = "F")
    if(max(dttmp$`Pr(F)`,na.rm = T)> (0.05) ) ttmp = rownames(dttmp)[which.max(dttmp$`Pr(F)`)]
    else break
    lmtmp = update(lmtmp, paste(".~. - ",ttmp) )
    print(dttmp)
    print(lmtmp)
  }
  lmsvaRed[[i]] = lmtmp
}

for(i in 1:K) write.csv( tidy( anova(lmsvaFull[[i]]) ) ,file = paste("QC/ANOVAfull_sva",i,".csv",sep=""),quote = F,row.names = F)
for(i in 1:K) write.csv( tidy( anova(lmsvaRed[[i]]) ) ,file = paste("QC/ANOVAreduced_sva",i,".csv",sep=""),quote = F,row.names = F)

#require(beeswarm)
pdf("QC/SVA.pdf",h=3*K,w=3*K)
par(mfrow=c(K,K),family="Times",las=1)
for(i in 1:K){
  for(j in 1:K){
    plot(sva[,j], sva[,i], col=rainbow(20)[as.factor(pData(RGSet)$sentrix_id)],pch=as.numeric(as.factor(pData(RGSet)$sentrix_position)) ,   
         xlab=paste("SV",j),ylab=paste("SV",j),main="Effects of Sentrix id (color) & Sentrix position (shape)")
  }
}
#legend("bottomright",legend = c( levels(as.factor(metaM[,22])),levels(as.factor(metaM[,13]))), col=c(colClustLNEN,1,1),pch=c(16,16,16,16,17) )
#beeswarm( sva[,7], pwcol=as.factor(metaM[as.character(pData(fun)$Ms_IDs),2]) , pch=16 ,xlab="",ylab="SV 7",main="Effect of Center of origin",method = "square",spacing = 2)
#legend("bottomright",legend = levels(as.factor(metaM[,2])), col=1:8,pch=16)
dev.off() ## to check

# Normalization
fun <- preprocessFunnorm(RGSet) # includes NOOB background/dye correction
rm(RGSet)
save(fun, file = "fun.RData")

# sex prediction
predictedSex <- getSex(gset, cutoff = -2)$predictedSex

pdf("Sex_pred.pdf",h=10,w=10)
plotSex(getSex(gset, cutoff = -2),id = pData(gset)$barcode)
dev.off()

pData = pData(fun)
pData$predictedSex = predictedSex

write.csv( data.frame("Sample"=as.character(pData(gset)$sample_id), "Prediction"=predictedSex,
                      "Clinical_data"=as.character(pData(gset)$sex),
                      row.names = pData(gset)$barcode), file = "Sex_pred.csv",quote = F)

print("Samples with mismatching sex and predicted sex:")
print( pData(gset)[which( as.character(pData(gset)$sex)!=predictedSex ),] )
rm(gset)

# Age Prediction
predicted_age <- agep(minfi::getBeta(fun))
pData$predictedAge = predicted_age
write.csv( data.frame("Sample"=as.character(pData$sample_id),
                      "Clinical_data"=as.character(pData$age),
                      "Prediction"=predicted_age,row.names = pData$barcode), file = "Age_pred.csv",quote = F)


pdf("Age_pred.pdf")
par(family="Times",las=1)
plot( pData(fun)$age, predicted_age,pch=16,xlim=c(0,120),ylim=c(0,140), xlab="Age", ylab="Predicted age")
abline(lm(pData(fun)$age~predicted_age[,1]), col="red")
ct = cor.test(as.numeric(pData(fun)$age), as.numeric(predicted_age))
text( 0,10 , paste0("r=",ct$estimate) ,pos = 4 )
text( 0,0 , paste0("p=",ct$p.value) ,pos = 4 )
dev.off()


# Smoking "prediction"
cg05575921 <- minfi::getBeta(fun)["cg05575921", ] # within AHRR locus
pdf("Smokig_prediction.pdf",h=6,w=6)
par(mfrow=c(1,1), mar=c(8,4,4,4),family="NimbusSan",las=2)
stripchart(cg05575921~pData(fun)$smoking_status,
             method="jitter",#group.names=c("Current","Former","Never","Passive"),
             pch=16,cex=1,col=c(4,2,3,1),ylab= "Normalized Beta values",
           main="cg05575921", vertical=TRUE,cex.axis=1.2,cex.lab=1)
boxplot(cg05575921,ylab= "Normalized Beta values", main="cg05575921 (smoking related)")
dev.off()

#########################
### Filtering samples ###
#########################
# remove duplicates
pData2 <- pData[order(pData(fun)$sample_id), ]
dim(pData2)
head(pData2)

dups       = pData2$sample_id[duplicated(pData2$sample_id)]
if( length(dups)>0 ){
    whdups     = lapply(dups, function(x) which(pData2$sample_id==x) )
    whdups2rem = sapply( 1:length(dups) , function(i) rbinom( 1 , 1 ,  prob = 1/length(whdups[[1]])  ) )+1
    torem = sapply( 1:length(whdups) , function(i) whdups[[i]][whdups2rem[i]] )
    pData2 <- pData2[-torem,]
}
fun2 <- fun[, rownames(pData2)]

save(fun2, file = "fun_nodup.RData")

########################
### Filtering Probes ###
########################
# remove XY probes
data(Locations)
XY <- rownames(Locations[Locations$chr%in%c("chrY","chrX"),])
fun2<- fun2[ ! featureNames(fun2) %in% XY, ]
rm(XY)
print(fun2)

# remove cross-reactive probes
if(!is.null(opt$crossreac)){
    print("Remove cross-reactive probes")
    Cross_reactive <- read.csv(opt$crossreac,header=F)$V1
    fun2<- fun2[ ! featureNames(fun2) %in% Cross_reactive, ] 
    rm(Cross_reactive)
    print(fun2)
}else{
    print("No file with cross-reactive probes supplied; to remove cross-reactive probes, use the -c option")
}

# remove SNP-containing probes
if(as.logical(opt$snp_filter)){
    print("Remove SNP-associated probes")
    fun2 = dropLociWithSnps(fun2,maf = 0.05) # alternative?
    print(fun2)
}else{
    print("Do not remove SNP-associated probes")
}

# remove any probes that have failed in one or more samples
detP2 = detP[rownames(fun2),colnames(fun2)]
rm(detP)
failed = which(rowSums(detP2 < 0.01) != ncol(fun2) )
fun2 <- fun2[ -failed ,]
rm(detP2)
rm(failed)
print(fun2)
save(fun2, file="fun2.RData")

# remove multimodal CpGs
if(as.logical(opt$multimodal_filter)){
    print("Remove multimodal beta-values probes")
    nmode<- nmode.mc(minfi::getBeta(fun2), minN = 3, modedist=0.2, nCores = 1)
    nmode.hi <- nmode[nmode>2]
    fun2<- fun2[ ! featureNames(fun2) %in% names(nmode.hi), ]

    pData(fun2)<- pData2
}else{
     print("Do not remove multimodal beta-values probes")
}
fun_filtered = fun2
rm(fun2)
save(fun_filtered, file="fun_filtered.RData")

#  Principal component regression analysis plots
# select variables of interest
head(pData(fun_filtered))
pData2 <- data.frame(pData(fun_filtered))
pData2$sentrix_id = as.factor(pData2$sentrix_id)

lev = apply(pData2,2,function(x) table(x))
pData2 = pData2[,(sapply(lev,max,na.rm=T)>1)&(sapply(lev,nrow)>1) ]


svg("PCA.svg",h=6,w=6)
pcrplot(minfi::getBeta(fun_filtered), cov= pData2, npc=20)
dev.off()

# write session info
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
