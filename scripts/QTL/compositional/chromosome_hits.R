library(qtl2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
    stop("At least one argument must be supplied (Population)\n", call.=FALSE)
} else if (length(args)==1) {
    pop <- as.character(args[1])
}

# out_pg <- readRDS(paste0(pop, ".operm.hk.Rdata"))
out_pg <- readRDS(paste0(pop, ".out.Rdata"))
map <- readRDS(paste0(pop, ".map.Rdata"))
# operm <- readRDS(paste0(pop, ".operm.Rdata"))
# as.numeric(summary(operm)[1])  # permutation threshold
lod_peaks <- find_peaks(out_pg, map, threshold=3, peakdrop=1.8, drop=1.5)

snps <- c()
for (row in 1:nrow(lod_peaks)){
    chr <- lod_peaks[row, "chr"]
    interval <- c(lod_peaks[row, "ci_lo"], lod_peaks[row, "ci_hi"])
    snps <- c(snps, paste(find_marker(map, chr, interval=interval), collapse=";"))
}
lod_peaks$snps <- snps
write.csv(lod_peaks, file=paste0(pop, ".hits.csv"), row.names=F, quote=F)
write("chrom hits", stderr())
