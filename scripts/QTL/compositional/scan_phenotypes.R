library(qtl2)
# library(qtlcharts)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (Population)\n", call.=FALSE)
} else if (length(args)==1) {
  pop <- as.character(args[1])
}

cores <- 16
num_phenotypes <- 4 # 8 with maturity

maturity <- read.csv("../maturity_DTH.csv")
rownames(maturity) <- maturity$ID
maturity <- maturity[,c(2,3)]
biomass <- read_cross2(paste0("rqtl2_", pop, ".yaml"))
biomass$pheno <- biomass$pheno[,c(1,2,6,7)]
map <- est_map(cross=biomass, 
                     map_function="haldane", 
                     cores=cores, 
                     save_rf=TRUE)
# map_plot <- iplotMap(map)
# htmlwidgets::saveWidget(map_plot, file=paste0(pop, ".iplotMap.html"))
biomass$gmap <- map

biomass$gmap <- insert_pseudomarkers(biomass$gmap, step=1)
saveRDS(biomass$gmap, paste0(pop, ".map.Rdata"))
pr <- calc_genoprob(biomass, biomass$gmap, cores=cores)
saveRDS(pr, paste0(pop, ".genoprob.pr.Rdata"))
apr <- genoprob_to_alleleprob(pr, cores=cores)
saveRDS(apr, paste0(pop, ".alleleProb.apr.Rdata"))

kinship <- calc_kinship(pr, cores=cores)
saveRDS(kinship, paste0(pop, ".kinship.pr.Rdata"))

out <- scan1(pr, biomass$pheno, addcovar=maturity, cores=cores)
saveRDS(out, paste0(pop, ".out.Rdata"))

# Process
peaks <- find_peaks(out, biomass$gmap, threshold=3, drop=1.5)
saveRDS(peaks, paste0(pop, ".peaks.Rdata"))
# bayes_int(out, map, lodcolumn=1, chr=3, prob=0.95)

out_pg <- scan1(pr, biomass$pheno, addcovar=maturity, kinship, cores=cores)
saveRDS(out_pg, paste0(pop, ".lmm.Rdata"))

kinship_loco <- calc_kinship(pr, "loco", cores=cores)
saveRDS(kinship_loco, paste0(pop, ".loco_kinship.Rdata"))

out_pg_loco <- scan1(pr, biomass$pheno, addcovar=maturity, kinship_loco, cores=cores)
saveRDS(out_pg_loco, paste0(pop, ".loco.Rdata"))

# Estimate heritability
herit <- est_herit(biomass$pheno, kinship)
saveRDS(herit, paste0(pop, ".herit.Rdata"))

# Plot H-K and LMM
for(i in 1:num_phenotypes) {
    png(paste0(pop, ".", i, ".H-K_LMM_LOCO.png"))
    color <- c("slateblue", "violetred", "green3")
    par(mar=c(4.1, 4.1, 1.6, 1.1))
    out_lod <- maxlod(out, lodcolumn=i)
    out_pg_lod <- maxlod(out_pg, lodcolumn=i)
    out_pg_loco_lod <- maxlod(out_pg_loco, lodcolumn=i)
    print(out_lod)
    print(out_pg_lod)
    print(out_pg_loco_lod)
    ymx <- max(out_lod, 
               out_pg_lod, 
               out_pg_loco_lod)
    plot(out, biomass$gmap, lodcolumn=i, col=color[1], 
         main=colnames(biomass$pheno)[i],
         ylim=c(0, ymx*1.02), type="l", pch=20)
    plot(out_pg, biomass$gmap, lodcolumn=i, 
         col=color[2], add=TRUE, type="l", pch=20)
    plot(out_pg_loco, biomass$gmap, lodcolumn=i, 
         col=color[3], add=TRUE, type="l", lty=2, pch=20)
    legend("topleft", lwd=2, col=color, c("H-K", "LMM", "LOCO"), bg="gray90", lty=c(1,1,2))
}
dev.off()

write("scan done", stderr())
