library('MASS') # required for ginv
library(multtest)
library(gplots)
library(compiler) #required for cmpfun
library("scatterplot3d")

# source("http://www.zzlab.net/GAPIT/emma.txt")
source("GAPIT.library.R")
source("gapit_functions.txt")

myY  <- read.csv("merged_blups.gapit.fixed.csv", head = TRUE)

myG <- read.delim("NAM.imputed.ann.fixed.hapmap.hmp.txt", head = FALSE)

# includes MAT,DTH,ADF,NDF,WSC
myCV <- read.csv("height.fixed.csv", head=TRUE)
myCV <- na.omit(myCV)

myY <- myY[,c(1,2,5,12,15,16,20,35)]
colnames(myY) <- c("Taxa", "ADF_adv", "aNDF_adv", "DryMatter_adv", "Lignin_adv", "Lignin_NDF_adv", "NDF_adv", "WSC_adv")

common <- intersect(myCV$Taxa, myY$Taxa)
myCV <- myCV[myCV$Taxa %in% common,]
myY <- myY[myY$Taxa %in% common,]
myCV <- myCV[order(common),]
myY <- myY[order(common),]
print(head(myY))
row.names(myY) <- NULL
row.names(myCV) <- NULL
myY <- head(myY, n=1203)
myCV <- head(myCV, n=1203)

myCV_WSC <- myCV[,c(1,2,3,6)]
myCV_NDF <- myCV[,c(1,2,3,5)]
myCV_ADF <- myCV[,c(1,2,3,4)]
myCV_NDF_height <- myCV[,c(1,2,3,5,7)]

myY_DryMatter <- myY[,c(1,4)]
myY_NDF <- myY[,c(1,7)]

colnames(myY_DryMatter) <- c("Taxa","Dry_covNDF_height")

myGAPIT <- GAPIT(
                 Y=myY_DryMatter,
                 G=myG,
                 CV=myCV_NDF_height,
                 Model.selection = TRUE,
                 model=c("BLINK"),
                 SNP.MAF=0.01,
                 SNP.FDR=0.05,
                 ncpus=8
                 )

colnames(myY_DryMatter) <- c("Taxa","Dry_covNDF")
print(head(myY_DryMatter))
print(head(myCV_NDF))

# Dry + cov NDF
# myGAPIT <- GAPIT(
#                  Y=myY_DryMatter,
#                  G=myG,
#                  CV=myCV_NDF,
#                  Model.selection = TRUE,
#                  model=c("BLINK"),
#                  SNP.MAF=0.01,
#                  SNP.FDR=0.05,
#                  ncpus=8
#                  )
# 
# print(head(myCV_ADF))
# print(head(myY_DryMatter))
# colnames(myY_DryMatter) <- c("Taxa","Dry_covADF")
# 
# # Dry + cov ADF
# myGAPIT <- GAPIT(
#                  Y=myY_DryMatter,
#                  G=myG,
#                  CV=myCV_ADF,
#                  Model.selection = TRUE,
#                  model=c("BLINK"),
#                  SNP.MAF=0.01,
#                  SNP.FDR=0.05,
#                  ncpus=8
#                  )
# 
# 
# colnames(myY_DryMatter) <- c("Taxa","Dry_covWSC")
# 
# 
# # Dry + cov WSC
# myGAPIT <- GAPIT(
#                  Y=myY_DryMatter,
#                  G=myG,
#                  CV=myCV_WSC,
#                  Model.selection = TRUE,
#                  model=c("BLINK"),
#                  SNP.MAF=0.01,
#                  SNP.FDR=0.05,
#                  ncpus=8
#                  )
# 
# colnames(myY_NDF) <- c("Taxa", "NDF_covADF")
# 
# # NDF + cov ADF
# myGAPIT <- GAPIT(
#                  Y=myY_NDF,
#                  G=myG,
#                  CV=myCV_ADF,
#                  Model.selection = TRUE,
#                  model=c("BLINK"),
#                  SNP.MAF=0.01,
#                  SNP.FDR=0.05,
#                  ncpus=8
#                  )
# 
# 
