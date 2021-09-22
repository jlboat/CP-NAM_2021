library(lme4)
library(data.table)

df <- read.csv("NIR_model_data.csv")

for (i in colnames(df)[4:37]){
    model_run <- lmer(paste0(i, " ~ (1|TAXA)"), data=df)
    random_effects <- setDT(ranef(model_run)$TAXA, keep.rownames=TRUE)[]
    colnames(random_effects) <- c("Taxa",i)
    write.csv(random_effects, file=paste0(i, "_BLUPs.csv"), quote=F, row.names=F)
}
