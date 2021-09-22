library(qtl2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
    stop("At least one argument must be supplied (Population)\n", call.=FALSE)
} else if (length(args)==1) {
    pop <- as.character(args[1])
}

cores <- 6
num_phenotypes <- 4

covar <- read.csv("../maturity_DTH.csv")
rownames(covar) <- covar$ID
covar <- covar[,c(2,3)]
biomass <- read_cross2(paste0("rqtl2_", pop, ".yaml"))
biomass$pheno <- biomass$pheno[,c(1,2,6,7)]
map <- readRDS(paste0(pop, ".map.Rdata"))
biomass$gmap <- map
pr <- readRDS(paste0(pop, ".genoprob.pr.Rdata"))
# kinship <- readRDS(paste0(pop, ".kinship.pr.Rdata"))
# kinship_loco <- calc_kinship(pr, "loco", cores=4)
operm <- scan1perm(pr, biomass$pheno, addcovar=covar, cores=cores, n_perm=1000)
# operm <- readRDS(paste0(pop, ".operm.Rdata"))
summary(operm, alpha=c(0.05))
saveRDS(operm, paste0(pop, ".operm.hk.Rdata"))
out <- readRDS(paste0(pop, ".out.Rdata"))
out_pg <- readRDS(paste0(pop, ".lmm.Rdata"))

# lod_peaks <- find_peaks(out, map, threshold=as.numeric(summary(operm)[1]), peakdrop=1.8, drop=1.5)
lod_peaks <- find_peaks(out, map, threshold=3, peakdrop=1.8, drop=1.5)
print(lod_peaks)
# png(paste0(pop, ".lodpeaks.png"))
# plot(lod_peaks)
# dev.off()


# peak_Mbp <- max(out_pg, biomass$pmap)$pos
# variants <- query_variants(2, peak_Mbp - 1, peak_Mbp + 1)
# out_snps <- scan1snps(pr, biomass$pmap, biomass$pheno, 
#                       kinship[[as.character(lod_peaks$chr[1])]], 
#                       query_func=query_variants,
#                       chr=lod_peaks$chr[1], start=peak_Mbp-1, end=peak_Mbp+1, keep_all_snps=TRUE)
# png(paste0(pop, ".snp_assoc.png"))
# par(mar=c(4.1, 4.1, 0.6, 0.6))
# plot(out_snps$lod, out_snps$snpinfo)
# dev.off()

# Plot operm
for(i in 1:num_phenotypes) {
    png(paste0(pop,".", i, ".operm.png"))
    color <- c("slateblue")
    par(mar=c(4.1, 4.1, 1.6, 1.1))
    ymx <- maxlod(operm)
    plot(out, map, lodcolumn=i, col=color[1], type="l", pch=20, altcol="orange")
    legend("topleft", lwd=2, col=color, c("H-K"), bg="gray90", lty=c(1,1,2))
    abline(h=summary(operm)[1], col="red")
    abline(h=3, col="red", lty=2)
}
dev.off()

# snpinfo <- create_snpinfo(biomass)
# snpinfo <- index_snps(biomass$pmap, snpinfo)
# # perform a SNP scan
# out_snp <- scan1snps(pr, biomass$pmap, biomass$pheno[,"pericarp"], snpinfo=snpinfo)
# png(paste0(pop, ".snp_scan.png"))
# # plot the LOD scores
# plot(out_snp$lod, snpinfo, altcol="green3")
# dev.off()
write("permute", stderr())
