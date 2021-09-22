library(qtl2)

cores <- 6
all_pops <- c("Chinese_Amber","Leoti",
              "PI_229841","PI_297130",
              "PI_297155","PI_329311",
              "PI_506069","PI_508366",
              "PI_510757","PI_655972","Rio")
phenotypes <- c("height","biomass",
                "DryWeight","WetWeight")
phenotypes_labels <- c("Height","Biomass",
                "DryWeight","WetWeight", "DTH")
num_phenotypes <- length(phenotypes)
colors <- rainbow(num_phenotypes + 1)
max_lods <- c(15,15,15,15)

for(i in 1:length(all_pops)) {
    pop <- all_pops[i]
    png(paste0("All_phenotypes.", all_pops[i], ".png"))
    #par(oma = c(5,4,4,0) + 0.3,
    #    mar = c(0,0,1,1) + 0.3)
    #    # mfrow=c(4,2))
    for (j in 1:num_phenotypes){
        pheno <- phenotypes[j]
        out <- readRDS(paste0("./", pop, "/", pop, ".out.Rdata"))
        map <- readRDS(paste0("./", pop, "/", pop, ".map.Rdata"))
        color <- colors[j]
        if (j==1){
            plot(out, map, lodcolumn=j, xlab="",
                  ylim=c(0, 15), col=color, type="l")
        } else {
            plot(out, map, lodcolumn=j, xlab="",
                  ylim=c(0, 15), col=color, type="l", add=TRUE)
        }
#         if ((i == 10) | (i == 11)){
#             plot(out, map, lodcolumn=j,  
#                  ylim=c(0, max_lods[j]), col=color, type="l")
#         } else if (i/2 == 0){
#             plot(out, map, lodcolumn=j, axes=F, xaxt='n', yaxt='n', 
#                  ylim=c(0, max_lods[j]), col=color, type="l")
#         } else {
#             plot(out, map, lodcolumn=j, xaxt='n', 
#                  ylim=c(0, max_lods[j]),
#                  col=color, type="l")
#         }
        abline(h=3, col="red", lty=2)
        # plot(out, map, lodcolumn=i, col=colors[j], add=TRUE, type="l")
    }
    out <- readRDS(paste0("./", pop, "/", pop, ".out.DTH.Rdata"))
    map <- readRDS(paste0("./", pop, "/", pop, ".map.DTH.Rdata"))
    color <- colors[num_phenotypes+1]
    plot(out, map, lodcolumn=1, xlab="",
                  ylim=c(0, 15), col=color, type="l", add=TRUE)
    legend(x=1, y=15, xpd=NA, ncol=1, lwd=2,
        col=colors, phenotypes_labels, bg="gray90")
    mtext(phenotypes[j], outer = TRUE, cex = 1.5)
    mtext("Chromosome",side=1,line=2.3,outer=TRUE,cex=1.3)
    dev.off()

}
    #mtext("LOD",side=2,line=1.5,outer=TRUE,cex=1.3,las=0)
