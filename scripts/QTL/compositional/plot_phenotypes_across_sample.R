library(qtl2)

cores <- 6
all_pops <- c("Chinese_Amber","Leoti",
              "PI_229841","PI_297130",
              "PI_297155","PI_329311",
              "PI_506069","PI_508366",
              "PI_510757","PI_655972","Rio")
all_phenotypes <- c("ADF", "AD.ICP", "Adj_CP", "aNDF", "aNDFom", "Ash", "Ca", "Cl",
"Crude.protein", "DCAD", "Dry.Matter", "EE.Fat", "K", "Lignin",
"Lignin_.NDF", "Mg", "Moisture", "Na", "NDF", "NDICP_w.oNa2SO3",
"NEG_OARDC", "NEL3x_ADF", "NEL3x_OARDC", "NEM_ADF", "NEM_OARDC",
"NFC", "P", "RFV", "S", "SP.CP", "Starch", "TDN_ADF", "TDN_OARDC",
"WSC_Sugar")
phenotypes <- c("Lignin_.NDF","WSC_Sugar",
                "Dry.Matter","DCAD","ADF",
                "TDN_ADF","Adj_CP")
phenotypes_labels <- c("Lignin (NDF)","WSC",
                "Dry.Matter","DCAD","ADF",
                "TDN (ADF)","Adj.CP")
num_phenotypes <- length(phenotypes)
colors <- rainbow(num_phenotypes)
max_lods <- c(15,15,15,4,15,4,15,15)

for(i in 1:length(all_pops)) {
    pop <- all_pops[i]
    png(paste0("All_phenotypes.", all_pops[i], ".png"))
    #par(oma = c(5,4,4,0) + 0.3,
    #    mar = c(0,0,1,1) + 0.3)
    #    # mfrow=c(4,2))
    first <- 0
    for (j in match(phenotypes, all_phenotypes)){
#        pheno <- all_phenotypes[j]
        out <- readRDS(paste0("./", pop, "/", pop, ".out.Rdata"))
        map <- readRDS(paste0("./", pop, "/", pop, ".map.Rdata"))
        color <- colors[match(j, match(phenotypes, all_phenotypes))]
        if (first==0){
            plot(out, map, lodcolumn=j, xlab="",
                  ylim=c(0, 15), col=color, type="l")
            first <- 1
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
    legend(x=1, y=15, xpd=NA, ncol=1, lwd=2,
        col=colors, phenotypes_labels, bg="gray90")
    mtext(all_phenotypes[j], outer = TRUE, cex = 1.5)
    mtext("Chromosome",side=1,line=2.3,outer=TRUE,cex=1.3)
    dev.off()

}
    #mtext("LOD",side=2,line=1.5,outer=TRUE,cex=1.3,las=0)
