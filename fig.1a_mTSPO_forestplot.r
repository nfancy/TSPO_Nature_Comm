library(forestplot)
library(dplyr)
library(greekLetters)

load("get_tspo_logFC.RData")
M1 <- readRDS("M1.rds")

mTSPO_table <- tibble::tibble(mean = mTSPO$logFC,
                              lower = mTSPO$CI.L,
                              upper = mTSPO$CI.R,
                              study = stringr::str_extract(M1$studlab, "E.GEOD.[0-9]+|E.MTAB.[0-9]+"),
                              treatment = gsub("E.GEOD.[0-9]+.eset.|E.MTAB.[0-9]+.eset.", "", M1$studlab),
                              n_control = as.character(M1$n.c),
                              n_stim = as.character(M1$n.e),
                              logFC = as.character(round(mTSPO$logFC, 2)),
                              CI = sprintf("[%s, %s]",
                                           as.character(round(mTSPO$CI.L, 2)),
                                           as.character(round(mTSPO$CI.R, 2))))



mTSPO_table <- mTSPO_table %>%
  mutate(study = gsub(".", "-", study, fixed = T)) %>%
  mutate(treatment = case_when(treatment == "IFNa" ~ paste0("IFN", greek_vector["alpha"]),
                               treatment == "IFNg" ~ paste0("IFN", greek_vector["gamma"]),
                               treatment == "LPS.IFNg" ~ paste0("LPS+IFN", greek_vector["gamma"]),
                               TRUE ~ treatment),
         time = c("2 hrs", "4 hrs", "24 hrs", "18 hrs", "2.5 hrs",
                  "2.5 hrs", "24 hrs", "24 hrs", "24 hrs", "24 hrs" )) %>%
  dplyr::relocate(time, .after = treatment) %>%
  arrange(treatment)


treatment_dose <- read.delim("treatment_dose.txt")
treatment_dose <- treatment_dose %>%
  mutate(treatment = case_when(treatment == "IFNa" ~ paste0("IFN", greek_vector["alpha"]),
                               treatment == "IFNg" ~ paste0("IFN", greek_vector["gamma"]),
                               treatment == "LPS+IFNg" ~ paste0("LPS+IFN", greek_vector["gamma"]),
                               TRUE ~ treatment))

mTSPO_table <- left_join(mTSPO_table, treatment_dose)
mTSPO_table$treatment_label <- sprintf("%s, %s, %s",
                                       mTSPO_table$treatment, 
                                       mTSPO_table$time,
                                       mTSPO_table$Concentration)



mTSPO_meta_logFC <- topTable(fit, coef = 2, Inf, confint = T)["Tspo", ]


summary <- tibble::tibble(mean = mTSPO_meta_logFC$logFC,
                          lower = mTSPO_meta_logFC$CI.L,
                          upper = mTSPO_meta_logFC$CI.R,
                          study = "Summary",
                          n_control = as.character(sum(M1$n.c)),
                          n_stim = as.character(sum(M1$n.e)),
                          logFC = as.character(round(mTSPO_meta_logFC$logFC, 2)),
                          CI = sprintf("[%s, %s]", 
                                       as.character(round(mTSPO_meta_logFC$CI.L, 2)), 
                                       as.character(round(mTSPO_meta_logFC$CI.R, 2))),
                          summary = TRUE)


header <- tibble::tibble(study = c("", "Study"),
                         treatment_label = c("", "Treatment"),
                         #time = c("", "Time"),
                         n_control = c("N", "(Control)"),
                         n_stim = c("N", "(Stim)"),
                         logFC = c("","logFC"),
                         CI = c("", "95%-CI"),
                         summary = TRUE)

empty_row <- tibble::tibble(mean = NA_real_)


footer <- tibble::tibble(study = sprintf("meta p = %s",
                                         as.numeric(format(M1$pval.Q, format = "e", digits = 2))),
                         summary = FALSE)

output_df <- dplyr::bind_rows(header,
                              mTSPO_table,
                              empty_row,
                              summary,
                              footer)



output_df %>% 
  forestplot(labeltext = c(study, treatment_label, n_control, n_stim, logFC, CI), 
             is.summary = summary,
             graph.pos = 6,
             #clip = c(-1.5, 2.5), 
             hrzl_lines = list("3" = gpar(lty = 2), 
                               "13" = gpar(lwd = 1, columns = c(1:4, 6), col = "#000044")),
             xlog = FALSE, 
             col = fpColors(box = "black",
                            line = "black",
                            summary = "navy"),
             txt_gp = fpTxtGp(ticks = gpar(cex = 1)),
             ci.vertices = TRUE,
             boxsize = 0.25)


png("forest.mTSPO.pro.png", width = 12, height = 5, units = "in", res = 300)
output_df %>% 
  forestplot(labeltext = c(study, treatment_label, n_control, n_stim, logFC, CI), 
             is.summary = summary,
             graph.pos = 6,
             #clip = c(-1.5, 2.5), 
             hrzl_lines = list("3" = gpar(lty = 2), 
                               "13" = gpar(lwd = 1, columns = c(1:4, 6), col = "#000044")),
             xlog = FALSE, 
             col = fpColors(box = "black",
                            line = "black",
                            summary = "navy"),
             txt_gp = fpTxtGp(ticks = gpar(cex = 1)),
             ci.vertices = TRUE,
             boxsize = 0.25)
dev.off()

mTSPO <- mTSPO %>%
  mutate(treatment = gsub("E.GEOD.[0-9]+.eset.|E.MTAB.[0-9]+.eset.", "", dataset_id) %>%
           gsub(".", "+", ., fixed = T),
         dataset_id = stringr::str_extract(dataset_id, "E.GEOD.[0-9]+|E.MTAB.[0-9]+") %>%
           gsub(".", "-", ., fixed = T),
         time = c("2 hrs", "4 hrs", "24 hrs", "18 hrs", "2.5 hrs",
                  "2.5 hrs", "24 hrs", "24 hrs", "24 hrs", "24 hrs" ),
         treatment = case_when(treatment == "IFNa" ~ paste0("IFN", greek_vector["alpha"]),
                                      treatment == "IFNg" ~ paste0("IFN", greek_vector["gamma"]),
                                      treatment == "LPS+IFNg" ~ paste0("LPS+IFN", greek_vector["gamma"]),
                                      TRUE ~ treatment)) %>%
  dplyr::select(c(dataset_id, treatment, time, everything()))

colnames(treatment_dose)[1] <- "dataset_id"

mTSPO <- left_join(mTSPO, treatment_dose)


mTSPO_meta_logFC$dataset_id <- "Meta-analysis"
mTSPO_meta_logFC$treatment <- NA
mTSPO_meta_logFC$time <- NA
mTSPO_meta_logFC$Concentration <- NA


meta_analysis_res <- rbind(mTSPO, mTSPO_meta_logFC[, c(26:28, 17:20, 22:25, 29)])
meta_analysis_res$P.Value[meta_analysis_res$dataset_id == "Meta-analysis"] <- M1$pval.Q
meta_analysis_res$adj.P.Val[meta_analysis_res$dataset_id == "Meta-analysis"] <- M1$pval.Q

meta_analysis_res <- meta_analysis_res %>%
  dplyr::rename(`Dataset ID` = dataset_id,
                Treatment = treatment,
                Time = time,
                `Gene Name` = Gene_Name) %>%
  arrange(Treatment)



write.table(meta_analysis_res[, c(1:3, 12, 4:11)], file = "mTSPO_meta_analysis_res.tsv", sep = "\t",
            row.names = F, quote = F)
