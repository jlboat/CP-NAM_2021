library('MASS') # required for ginv
library(multtest)
library(gplots)
library(compiler) #required for cmpfun
library("scatterplot3d")

# source("http://www.zzlab.net/GAPIT/emma.txt")
source("GAPIT.library.R")
source("gapit_functions.txt")

myY  <- read.csv("merged_blups.gapit.fixed.csv", head = TRUE)
# myY <- read.csv("Kathleen_phenotypes_no_mixed.common.csv", head=T)
# colnames(myY) <- c("Taxa","Phenotype")
# myY$Phenotype <- as.numeric(myY$Phenotype)

myG <- read.delim("NAM.imputed.ann.fixed.hapmap.hmp.txt", head = FALSE)

#taxa <- c("rs#", "alleles", "chrom", "pos", "strand", "assembly#", 
#          "center", "protLSID", "assayLSID", "panelLSID", "QCcode")
#taxa <- c(taxa, myG[1,t(myG[1,]) %in% myY$Taxa])
#myG <- myG[,as.character(t(myG[1,])) %in% taxa]

#Step 2: Run GAPIT
# myGAPIT <- GAPIT(
#                  Y=myY,
#                  G=myG,
#                  Model.selection = TRUE,
#                  model=c("MLM","BLINK"),
#                  SNP.MAF=0.01,
#                  SNP.FDR=0.05,
#                  ncpus=8
#                  )
myCV <- read.csv("Maturity_DTH_covariate.final.csv", head=TRUE)
myY <- myY[,c(1,2,5,12,15,16,20,35)]
colnames(myY) <- c("Taxa", "ADF_cov", "aNDF_cov", "DryMatter_cov", "Lignin_cov", "Lignin_NDF_cov", "NDF_cov", "WSC_cov")

common <- intersect(myCV$Taxa, myY$Taxa)
myCV <- myCV[myCV$Taxa %in% common,]
print(head(myCV))
myY <- myY[myY$Taxa %in% common,]
print(head(myY))
row.names(myY) <- NULL
row.names(myCV) <- NULL

myGAPIT <- GAPIT(
                 Y=myY,
                 G=myG,
                 CV=myCV,
                 Model.selection = TRUE,
                 model=c("BLINK"),
                 SNP.MAF=0.01,
                 SNP.FDR=0.05,
                 ncpus=8
                 )

