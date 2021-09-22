library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(viridis)
library(purrr)
library(GGally)
library(cluster)
library(stringr)
library(mashr)
library(CDBNgenomics)

source("https://raw.githubusercontent.com/Alice-MacQueen/CDBNgenomics/master/R/handle_mash_results.R")
source("https://raw.githubusercontent.com/Alice-MacQueen/CDBNgenomics/6b00f48eb1c6eec848f11416d7a5fd752cd778bd/R/cdbn_kinship.R")

mash_plot_manhattan_by_condition_custom <- function (m, cond = NA, saveoutput = FALSE, suffix = "", thresh = 0.05) 
{
    # .data$Num_Sig_Conditions
    num_sig_in_cond <- c()
    if (is.na(cond)[1]) {
        cond <- get_colnames(m = m)
    }
    log10bf_df <- get_log10bf(m = m) %>% as.data.frame() %>% 
        rownames_to_column(var = "value") %>% mutate(value = as.integer(.data$value)) %>% 
        as_tibble() %>% left_join(get_marker_df(m = m)) %>% dplyr::rename(log10BayesFactor = .data$V1) %>% 
        dplyr::select(-.data$value)
    write.csv(log10bf_df, file="log10bf.csv")
    ggman_df <- get_n_significant_conditions(m = m, thresh = thresh, 
        conditions = cond) %>% enframe(name = "Marker", value = "Num_Sig_Conditions") %>% 
        separate(.data$Marker, into = c("Chr", "Pos"), remove = FALSE, 
            sep = "_", extra = "merge", convert = TRUE) %>% left_join(log10bf_df, 
        by = "Marker") %>% arrange(.data$Chr, .data$Pos)
    log10BF <- expression(paste("log"[10], plain("(Bayes Factor)")))
    ggmanobject <- ggplot(data = ggman_df, aes(x = .data$Pos, 
        y = .data$log10BayesFactor)) + CDBNgenomics::theme_oeco + 
        geom_point(aes(color = .data$Num_Sig_Conditions, fill = .data$Num_Sig_Conditions, 
            shape = as.factor(.data$Chr)), color="black") + facet_wrap(~.data$Chr, 
        nrow = 1, scales = "free_x", strip.position = "bottom") + 
        scale_color_viridis(option = "B") + scale_fill_viridis(option = "B") + 
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
            panel.background = element_rect(fill = NA)) + labs(x = "Chromosome", 
        y = log10BF) + scale_x_continuous(expand = c(0.25, 0.25)) + 
        scale_shape_manual(values = rep(c(21, 22), 9), guide = FALSE)
    if (saveoutput == TRUE) {
        if (!(str_sub(suffix, end = 1) %in% c("", "_"))) {
            suffix <- paste0("_", suffix)
        }
        if (str_sub(suffix, end = 1) %in% c("")) {
            suffix <- paste0("_", get_date_filename())
        }
        save_plot(paste0("Manhattan_mash", suffix, ".png"), plot = ggmanobject, 
            base_aspect_ratio = 2.4, base_height = 3)
    }
    return(list(ggman_df = ggman_df, ggmanobject = ggmanobject))
}

mash_output <- readRDS("m2.rds")
nbycond <- mash_plot_sig_by_condition(m = mash_output)
png("SigSNPsPerCondition.png")
nbycond$ggobject
dev.off()

mashhattan <- mash_plot_manhattan_by_condition_custom(m = mash_output)
png("mashhattan.png")
mashhattan$ggmanobject
dev.off()

effects <- mash_plot_effects(m = mash_output, n = 1)
png("mash_effects.png")
effects$ggobject +
    scale_x_discrete(labels = str_sub(as_vector(effects$effect_df$mn),
                                      start = 6)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


pairwise_plot <- mash_plot_pairwise_sharing(m = mash_output,
                                            reorder = TRUE) 
png("mash_corr.png")
pairwise_plot$gg_corr
dev.off()
