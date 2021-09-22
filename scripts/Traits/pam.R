library(cluster)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
    stop("At least one argument must be supplied (Population).n", call.=FALSE)
} else if (length(args)==1) {
    population <- as.character(args[1])
}

if (population == "Chinese_Amber"){
    # Looks rough
    # Recommend selecting from three clusters separately
    k_count <- 3
} else if (population == "PI_297130"){
    # The selection looks a bit skewed, but no major cluster differences after filter
    # Recommend seleting from two clusters separately
    k_count <- 2
} else if (population == "PI_329311"){
    # The selection looks rough here, but no major cluster differences after filter
    # Recommend seleting from two clusters separately
    k_count <- 2
} else if (population == "PI_506069"){
    # Filter and keep as is
    k_count <- 4
} else if (population == "Rio"){
    # Rio appears decent
    # Recommend seleting from two clusters separately
    # NO LONGER TWO CLUSTERS AFTER FILTERING 250
    k_count <- 0
} else {
    k_count <- 0
}

# Get kinship
kinship <- readRDS(paste0(population,".kinship.pr.Rdata"))
print(paste0("Individuals in kinship: ", dim(kinship)[1]))

grassl_num <- grep("Grassl", row.names(kinship))
alt_num    <- grep(population, row.names(kinship))
all_parents <- c(grassl_num, alt_num)
parent_count <- length(all_parents)
print(paste0("Number of parents: ", parent_count))

# Filter out parents
#filtered_kinship <- kinship[(parent_count + 1):dim(kinship)[1],
#                            (parent_count + 1):dim(kinship)[1]]
filtered_kinship <- kinship[-all_parents, -all_parents]
print(paste0("Number of parents removed: ", dim(kinship)[1] - dim(filtered_kinship)[1]))

if (population == "PI_510757"){
    indvs <- c("15.2_050", "15.2_067", "15.2_L082", "15.2_084", "15.2_089")
    alt_num <- match(indvs, row.names(kinship))
    print(alt_num)
    filtered_kinship <- filtered_kinship[-alt_num, -alt_num]
} else if (population == "Rio"){
    indvs <- as.character(read.table("Rio.250.txt", header=F)$V1)
    indvs <- indvs[!is.na(indvs)]
    print(indvs)
    alt_num <- match(indvs, row.names(kinship))
    alt_num <- alt_num[!is.na(alt_num)]
    print(alt_num)
    filtered_kinship <- filtered_kinship[alt_num, alt_num]
}

# Distance calculation
dissim <- 1 - filtered_kinship
# dissim <- max(kinship) - kinship

# Classical MDS
Res <- cmdscale(dissim)

# Heatmap
png(paste0("./images/", 
           population, 
           "/",
           population,
           ".heatmap.png"))
heatmap(dissim)
dev.off()

# PAM clustering of genotypes
if ((k_count == 3) | (k_count == 2)){
    indivs <- c()
    kmeans_pop <- kmeans(dissim, k_count)
    for (i in 1:k_count){
        i_indivs <- kmeans_pop$cluster[kmeans_pop$cluster == i]
        # print(i_indivs)
        indivs_per_cluster <- round((110/dim(dissim)[1]) * length(i_indivs))
        print(paste0("Individuals to select: ", indivs_per_cluster, 
                     " in cluster ", i,
                     " compared to total " , dim(dissim)[1]))
        alt_num <- c()
        for (indv in names(i_indivs)){
            alt_num <- c(alt_num, grep(indv, row.names(kinship)))
        }
        print(alt_num)
        indiv_filtered_kinship <- filtered_kinship[alt_num, alt_num]
        print(dim(indiv_filtered_kinship))
        filtered_dissim <- 1 - indiv_filtered_kinship
        clust.result <- pam(x=filtered_dissim, 
                            k=indivs_per_cluster, 
                            diss=TRUE)
        indivs <- c(indivs, clust.result$medoids)
    }
    print(paste0("Number of medoids: ", length(indivs)))
} else {
    clust.result <- pam(x=filtered_kinship, k=110, diss=TRUE)
    indivs <- clust.result$medoids
    print(paste0("Number of medoids: ", length(indivs)))
}

# Write 110 medoids
write.csv(indivs, paste0(population, ".medoids.txt"), 
          row.names=FALSE, quote=FALSE)

# Colors
color <- c()
for ( i in row.names(dissim)){
    if (i %in% indivs){
        color <- c(color, "red")
    } else { color <- c(color, "black") }
}

print(paste0("Medoids in red: ", sum(str_count(color, pattern="red"))))

# Plot MDS components -- colored medoids
png(paste0("./images/", 
           population, 
           "/",
           population,
           ".MDS.png"))
plot(Res, col=color, pch=19, main="Multidimensional Scaling", cex=0.4, xlab="",ylab="")
dev.off()

if (k_count != 0){
    kmeans_pop <- kmeans(dissim, k_count)
    png(paste0("./images/",
               population,
               "/",
               population,
               ".KMeans.MDS.png"))
    plot(Res, 
         col=c("red","blue","green","grey")[kmeans_pop$cluster], 
         pch=19, main="Multidimensional Scaling", cex=0.4, xlab="",ylab="")
    dev.off()

    # Missing data after converting H to -
    missing_data <- read.csv(paste0(population, ".missing_after_conversion.csv"), 
                               header=TRUE)
    colnames(missing_data) <- c("row_num","INDV", "F_MISS")
    missing_data$CLUSTER <- kmeans_pop$cluster[missing_data$INDV]
    png(paste0("./images/",
               population,
               "/",
               population,
               ".noH.boxplot.png"))
    boxplot(F_MISS ~ CLUSTER, data=missing_data, 
            col=c("red","blue","green","grey"),
            main="Missing after H -> missing")
    dev.off()

    # Missing data before converting H to -
     missing_data <- read.table(paste0(population, ".Marker.Indiv.MafFiltered.imiss"), header=TRUE)
    missing_data$CLUSTER <- kmeans_pop$cluster[missing_data$INDV]
    png(paste0("./images/",
               population,
               "/",
               population,
               ".boxplot.png"))
    boxplot(F_MISS ~ CLUSTER, data=missing_data, 
            col=c("red","blue","green","grey"),
            main="Missing before H -> missing")
    dev.off()
}

pdf(paste0("./images/", 
           population, 
           "/",
           population,
           ".MDS.pdf"))
plot(Res, col=color, pch=19, main="Multidimensional Scaling", cex=0.4, xlab="",ylab="")
dev.off()

