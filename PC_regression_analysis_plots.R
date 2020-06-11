################
# Principal component regression analysis plots
###############
# select variables of interest
head(pData(fun_filtered))
pData2 <- data.frame(pData(fun_filtered))
pData2$sentrix_id = as.factor(pData2$sentrix_id)

lev = apply(pData2,2,function(x) table(x))
pData2 = pData2[,(sapply(lev,max,na.rm=T)>1)&(sapply(lev,nrow)>1) ]

#svg("PCA.svg",h=6,w=6)
pcrplot(minfi::getBeta(fun_filtered), cov= pData2, npc=20)
#dev.off()

# write session info
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
