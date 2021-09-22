library(fastDummies)
library(heritability)

# Example
# data(LDV)
# data(K_atwell)

# out1 <- marker_h2(data.vector=LDV$LDV,geno.vector=LDV$genotype, covariates=LDV[,4:8],K=K_atwell)

# Gemma formatted input
phenotypes <- read.table("CP-NAM_NIR_ALL.gemma.tsv", header=TRUE)
phenotypes <- phenotypes[order(phenotypes$TAXA.1),]
phenotypes <- phenotypes[,c(1,6:40)]
kinship <- read.table("./MV/output/plink.nam_stringent.sXX.txt")
pheno_taxa_order <- phenotypes$TAXA
accession_order <- read.table("./MV/plink.nam_stringent.fam")$V1
rownames(kinship) <- accession_order
colnames(kinship) <- accession_order
kinship <- as.matrix(kinship)
overlap <- intersect(pheno_taxa_order, accession_order)
df <- phenotypes[phenotypes$TAXA %in% overlap,]
# c1.rle <- rle(as.vector(df$Rep))
# repeated_genotypes <- paste0(rep(c1.rle$values, times = c1.rle$lengths), "_",
#        unlist(lapply(c1.rle$lengths, seq_len)))
# covariates <- dummyVars(" ~ .", df)
# covariates <- data.frame(predict(covariates, newdata = df))
# rownames(covariates) <- repeated_genotypes
covariates <- dummy_cols(df$Rep)[,2:4] # one for each rep
for (col in colnames(df)[2:length(colnames(df))]){
    print(paste0("Phenotype: ", col))
    print("Narrow-sense heritability")
    out1 <- marker_h2(data.vector=df[,col], geno.vector=df$TAXA, 
                      covariates=covariates, K=kinship)
    cat(paste0(col, ",narrow,", out1$va, ",",
                 out1$ve, ",", out1$h2, ",",
                 out1$conf.int1[1], ",", out1$conf.int1[2], "\n"))
    # print(out1)
    print("Broad-sense heritability")
    out2 <- repeatability(data.vector=df[,col], geno.vector=df$TAXA, 
                      covariates.frame=covariates[,1:2])
    cat(paste0(col, ",broad,", out2$repeatability, ",",
                 out2$gen.variance, ",", out2$res.variance, ",",
                 out2$conf.int[1], ",", out2$conf.int[2], "\n"))
    # print(out2)
}

# df <- read.csv("Heritability_phenotype_matrix.csv", row.names=1)
# kinship <- read.table("plink.ama_tan_BLUP.sXX.txt")
# blup_pi <- read.table("BLUP.txt")
# pi_names <- read.table("pi_order.txt")
# colnames(kinship) <- pi_names$V1
# rownames(kinship) <- pi_names$V1
# kinship <- as.matrix(kinship)
# overlap <- intersect(blup_pi$V1, pi_names$V1)
# kinship <- kinship[overlap, overlap]
# df <- df[df$PI %in% overlap,]
## overlap <- unique(df$PI)[is.element(unique(df$PI), rownames(K))]
## kinship <- kinship[ %in% overlap, pi_names %]
# out1 <- marker_h2(data.vector=df$cp_AMA, geno.vector=df$PI, covariates=df[6], K=kinship)
