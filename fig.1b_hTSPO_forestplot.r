library(forestplot)
library(dplyr)
library(greekLetters)

load("get_tspo_logFC.RData")
M1 <- readRDS("M1.rds")


hTSPO <- left_join(hTSPO %>%
                     dplyr::rename(esetName = dataset_id), 
                   dataList[, c(1,4)] )

hTSPO <- hTSPO %>%
  mutate(dataset_id = stringr::str_extract(esetName, "E.GEOD.[0-9]+|E.MTAB.[0-9]+"),
         dataset_id = gsub(".", "-", dataset_id, fixed = T),
         treatment = gsub("E.GEOD.[0-9]+.eset.|E.MTAB.[0-9]+.eset.", "", esetName),
         treatment = gsub(".[0-9]{1,2}[a-z]{1,2}", "", treatment),
         treatment = gsub(".", "+", treatment, fixed = T),
         treatment = case_when(treatment == "SA" ~ "S. aureus",
                               dataset_id == "E-GEOD-24897" & treatment == "LPS" ~ "P. gingivalis LPS",
                               treatment == "fimbriae" ~ "P. gingivalis fimbriae",
                               treatment == "live" ~ "P. gingivalis live",
                               treatment == "GBY" ~ "Glucan (GBY)",
                               treatment == "WT" ~ "L. pneumophila WT",
                               treatment == "ankB" ~ "L. pneumophila ankB",
                               treatment == "T2SS" ~ "L. pneumophila T2SS",
                               treatment == "IFNg" ~ paste0("IFN", greek_vector["gamma"]),
                               treatment == "IFNa" ~ paste0("IFN", greek_vector["alpha"]),
                               treatment == "LPS+IFNg" ~ paste0("LPS+IFN", greek_vector["gamma"]),
                               TRUE ~ treatment
         )) %>%
  mutate(Time = case_when(Time == "8hrs" ~ "8 hrs",
                   TRUE ~ Time)) %>%
  dplyr::select(c( dataset_id, treatment, Time, everything())) 

treatment_dose <- read.delim("treatment_dose.txt")
treatment_dose <- treatment_dose %>%
  mutate(treatment = case_when(treatment == "IFNg" ~ paste0("IFN", greek_vector["gamma"]),
                               treatment == "IFNa" ~ paste0("IFN", greek_vector["alpha"]),
                               treatment == "LPS+IFNg" ~ paste0("LPS+IFN", greek_vector["gamma"]),
                               TRUE ~ treatment))

hTSPO <- left_join(hTSPO, treatment_dose)
hTSPO$treatment_label <- sprintf("%s, %s, %s",
                           hTSPO$treatment, 
                           hTSPO$Time,
                           hTSPO$Concentration)



hTSPO_table <- tibble::tibble(mean = hTSPO$logFC,
                              lower = hTSPO$CI.L,
                              upper = hTSPO$CI.R,
                              study = hTSPO$dataset_id,
                              #treatment = hTSPO$treatment,
                              #time = hTSPO$Time,
                              treatment = hTSPO$treatment_label,
                              n_control = as.character(M1$n.c),
                              n_stim = as.character(M1$n.e),
                              logFC = as.character(round(hTSPO$logFC, 2)),
                              CI = sprintf("[%s, %s]",
                                           as.character(round(hTSPO$CI.L, 2)),
                                           as.character(round(hTSPO$CI.R, 2)))
                              )

hTSPO_table <- hTSPO_table  %>%
  arrange(treatment)

hTSPO_meta_logFC <- topTable(fit, coef = 2, Inf, confint = T)["TSPO", ]

summary <- tibble::tibble(mean = hTSPO_meta_logFC$logFC,
                          lower = hTSPO_meta_logFC$CI.L,
                          upper = hTSPO_meta_logFC$CI.R,
                          study = "Summary",
                          n_control = as.character(sum(M1$n.c)),
                          n_stim = as.character(sum(M1$n.e)),
                          logFC = as.character(round(hTSPO_meta_logFC$logFC, 2)),
                          CI = sprintf("[%s, %s]", 
                                       as.character(round(hTSPO_meta_logFC$CI.L, 2)), 
                                       as.character(round(hTSPO_meta_logFC$CI.R, 2))),
                          summary = TRUE)


header <- tibble::tibble(study = c("", "Study"),
                         treatment = c("", "Treatment"),
                         #time = c("", "Time"),
                         n_control = c("N", "(Control)"),
                         n_stim = c("N", "(Stim)"),
                         logFC = c("","logFC"),
                         CI = c("", "95%-CI"),
                         summary = TRUE)

empty_row <- tibble::tibble(mean = NA_real_)


footer <- tibble::tibble(study = sprintf("meta p = %s", round(M1$pval.Q, 2)),
                         summary = FALSE)

output_df <- dplyr::bind_rows(header,
                       hTSPO_table,
                       empty_row,
                       summary,
                       footer)



output_df %>% 
  forestplot(labeltext = c(study, treatment, n_control, n_stim, logFC, CI), 
             is.summary = summary,
             graph.pos = 6,
             #clip = c(0.1, 2.5), 
             hrzl_lines = list("3" = gpar(lty = 2), 
                               "45" = gpar(lwd = 1, columns = c(1:4, 6), col = "#000044")),
             xlog = FALSE, 
             col = fpColors(box = "black",
                            line = "black",
                            summary = "navy"),
             txt_gp = fpTxtGp(ticks = gpar(cex = 1)),
             ci.vertices = TRUE,
             boxsize = 0.25)


png("forest.hTSPO.pro.png", width = 12, height = 10, units = "in", res = 300)
output_df %>% 
  forestplot(labeltext = c(study, treatment, n_control, n_stim, logFC, CI), 
             is.summary = summary,
             graph.pos = 6,
             #clip = c(0.1, 2.5), 
             hrzl_lines = list("3" = gpar(lty = 2), 
                               "45" = gpar(lwd = 1, columns = c(1:4, 6), col = "#000044")),
             xlog = FALSE, 
             col = fpColors(box = "black",
                            line = "black",
                            summary = "navy"),
             txt_gp = fpTxtGp(ticks = gpar(cex = 1)),
             ci.vertices = TRUE,
             boxsize = 0.25)
dev.off()

hTSPO_meta_logFC$dataset_id <- "Meta-analysis"
hTSPO_meta_logFC$treatment <- NA
hTSPO_meta_logFC$Time <- NA
hTSPO_meta_logFC$esetName <- NA
hTSPO_meta_logFC$Concentration <- NA
hTSPO_meta_logFC$treatment_label <- NA

meta_analysis_res <- rbind(hTSPO, hTSPO_meta_logFC[, c(25:27, 16:19, 21:24, 28:30)])
meta_analysis_res$P.Value[meta_analysis_res$dataset_id == "Meta-analysis"] <- M1$pval.Q
meta_analysis_res$adj.P.Val[meta_analysis_res$dataset_id == "Meta-analysis"] <- M1$pval.Q

meta_analysis_res <- meta_analysis_res %>%
  mutate(treatment = gsub("E.GEOD.[0-9]+.eset.|E.MTAB.[0-9]+.eset.", "", esetName),
         treatment = gsub(".[0-9]{1,2}[a-z]{1,2}", "", treatment),
         treatment = gsub(".", "+", treatment, fixed = T),
         treatment = case_when(treatment == "SA" ~ "S. aureus",
                               dataset_id == "E-GEOD-24897" & treatment == "LPS" ~ "P. gingivalis LPS",
                               treatment == "fimbriae" ~ "P. gingivalis fimbriae",
                               treatment == "live" ~ "P. gingivalis live",
                               treatment == "GBY" ~ "Glucan (GBY)",
                               treatment == "WT" ~ "L. pneumophila WT",
                               treatment == "ankB" ~ "L. pneumophila ankB",
                               treatment == "T2SS" ~ "L. pneumophila T2SS",
                               TRUE ~ treatment)) %>%
  dplyr::rename(`Dataset ID` = dataset_id,
                Treatment = treatment,
                `Gene Name` = Gene_Name) %>%
  arrange(Treatment)

write.table(meta_analysis_res[ , c(1:3, 13, 4:11)], file = "hTSPO_meta_analysis_res.tsv", sep = "\t",
            row.names = F, quote = F)
