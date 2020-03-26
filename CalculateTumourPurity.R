# Function to calculate tumour purity from 44 LUMP probes
# Aran, Dvir, Marina Sirota, and Atul J. Butte. "Systematic pan-cancer analysis of tumour purity." Nature communications 6 (2015): 8971.

# Required function:
CalculateLump <- function(BetaTable) {
  LumpProbes <- c("cg00240653","cg00450164","cg00880290","cg00933696","cg01138020","cg02026204",
                  "cg02053964","cg02167021","cg02997560","cg03431741","cg03436397","cg03841065", 
                  "cg04915566","cg05199874","cg05305434","cg05769344","cg05798125","cg07002058",
                  "cg07598052","cg07641284","cg08854008","cg09302355","cg09606470","cg10511890",
                  "cg10559416","cg13030790","cg13912307","cg14076977","cg14913777","cg17518965",
                  "cg19466818","cg20170223","cg20695297","cg21164509","cg21376733","cg22331159",
                  "cg23114964","cg23553480","cg24796554","cg25384897","cg25574765","cg26427109",
                  "cg26842802","cg27215100")
  GetRowNumber <- function (row_name) {
    RowNumbers <- which(rownames(BetaTable)==row_name)
    return(RowNumbers)
  }
  LumpRows <- as.numeric(sapply(LumpProbes, GetRowNumber))
  BetaLump <- BetaTable[LumpRows,]
  BetaLump <- na.omit(BetaLump)
  Mean <- apply(BetaLump,2,mean)
  LumpCalc <- function(Mean) {
    result1 <- min((Mean/0.85),1)
    return(result1)
  }
  Purity <- sapply(Mean, LumpCalc)
  FinalTable <- rbind(BetaLump, Mean, Purity)
  rownames(FinalTable)[c(nrow(FinalTable)-1,nrow(FinalTable))] <- c("Mean", "Purity")
  return(FinalTable)
}

# Takes as input a dataframe of beta values (beta) for samples of interest generated from HM450k or EPIC methylation arrays
Result <- CalculateLump(beta)



