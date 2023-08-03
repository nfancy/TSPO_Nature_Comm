library(ggplot2)
library(dplyr)
library(DESeq2)
library(patchwork)


load("meta.analysis.tspo/E-GEOD-19315.eset.RData")
dim(E.GEOD.19315.eset)
col_idx <- which(pData(E.GEOD.19315.eset)$Treatment %in% c("untreated", "lipopolysaccharide"))
row_idx <- which(fData(E.GEOD.19315.eset)$V4 == "ENSG00000100300")
exprs(E.GEOD.19315.eset)[row_idx, col_idx]


tmp <- data.frame(dataset_id = "E.GEOD.19315",
                  simulation = pData(E.GEOD.19315.eset)$Treatment[col_idx],
                  duration = "6h",
                  species = "Human",
                  exprs = exprs(E.GEOD.19315.eset)[row_idx, col_idx],
                  logFC = "log2FC=-0.11",
                  padj = "p.adj=0.759")
dt <- tmp

load("multi_species_tspo_expression/E-MTAB-7572_human/E-MTAB-7572.dds.rdata")
dim(dds)
plotCounts(dds, gene = "TSPO", intgroup = "Treatment",
           normalized = T, transform = F, returnData = T)

col_idx <- which(colData(dds)$Treatment %in% c("none_24", "lipopolysaccharide_24"))
row_idx <- which(rownames(dds) == "TSPO")
counts(dds, normalized = T)[row_idx, col_idx]

tmp <- data.frame(dataset_id = "E.MTAB.7572",
                  simulation = colData(dds)$Treatment[col_idx],
                  duration = "24h",
                  species = "Human",
                  exprs = log2(counts(dds, normalized = T)[row_idx, col_idx]),
                  logFC = "log2FC=-0.65",
                  padj = "p.adj=2.81e-11")

dt <- rbind(dt, tmp)


load("meta.analysis.tspo.mouse/E-GEOD-15610.eset.RData")
dim(E.GEOD.15610.eset)
col_idx <- which(pData(E.GEOD.15610.eset)$Treatment %in% c("untreated", "LPSstimulated"))
row_idx <- which(fData(E.GEOD.15610.eset)$Gene_Name == "Tspo")
exprs(E.GEOD.15610.eset)[row_idx, col_idx]


tmp <- data.frame(dataset_id = "E.GEOD.15610",
                  simulation = pData(E.GEOD.15610.eset)$Treatment[col_idx],
                  duration = "4h",
                  species = "Mouse",
                  exprs = exprs(E.GEOD.15610.eset)[row_idx, col_idx],
                  logFC = "log2FC=0.69",
                  padj = "p.adj=0.0417")

dt <- rbind(dt, tmp)

load("meta.analysis.tspo.mouse/E-GEOD-53986.eset.RData")
dim(E.GEOD.53986.eset)
col_idx <- which(pData(E.GEOD.53986.eset)$Treatment %in% c("untreated", "LPS"))
row_idx <- which(fData(E.GEOD.53986.eset)$Gene_Name == "Tspo")
exprs(E.GEOD.53986.eset)[row_idx, col_idx]


tmp <- data.frame(dataset_id = "E.GEOD.53986",
                  simulation = pData(E.GEOD.53986.eset)$Treatment[col_idx],
                  duration = "24h",
                  species = "Mouse",
                  exprs = exprs(E.GEOD.53986.eset)[row_idx, col_idx],
                  logFC = "log2FC=2.04",
                  padj = "p.adj=9.17e-06")

dt <- rbind(dt, tmp)


load("multi_species_tspo_expression/E-MEXP-3469_Rat_macrophage_LPS_prashant/E.MEXP.3469.eset.rdata")
dim(E.MEXP.3469.eset)
col_idx <- which(pData(E.MEXP.3469.eset)$Treatment %in% 
                   c("none.0hour", "lipopolysaccharide.4hour"))
row_idx <- which(rownames(E.MEXP.3469.eset) == "Tspo")
exprs(E.MEXP.3469.eset)[row_idx, col_idx]


tmp <- data.frame(dataset_id = "E.MEXP.3469",
                  simulation = pData(E.MEXP.3469.eset)$Treatment[col_idx],
                  duration = "4h",
                  species = "Rat",
                  exprs = exprs(E.MEXP.3469.eset)[row_idx, col_idx],
                  logFC = "log2FC=0.21",
                  padj = "p.adj=0.0373")

dt <- rbind(dt, tmp)

col_idx <- which(pData(E.MEXP.3469.eset)$Treatment %in% 
                   c("none.0hour", "lipopolysaccharide.8hour"))
row_idx <- which(rownames(E.MEXP.3469.eset) == "Tspo")
exprs(E.MEXP.3469.eset)[row_idx, col_idx]

tmp <- data.frame(dataset_id = "E.MEXP.3469",
                  simulation = pData(E.MEXP.3469.eset)$Treatment[col_idx],
                  duration = "8h",
                  species = "Rat",
                  exprs = exprs(E.MEXP.3469.eset)[row_idx, col_idx],
                  logFC = "log2FC=1.23",
                  padj = "p.adj=1.92e-08")

dt <- rbind(dt, tmp)


load("multi_species_tspo_expression/E-GEOD-6353_rat/E.GEOD.6353.eset.rdata")
dim(E.GEOD.6353.eset)
col_idx <- which(pData(E.GEOD.6353.eset)$Treatment %in% 
                   c("Control", "IFNg"))
row_idx <- which(fData(E.GEOD.6353.eset)$Composite.Element.Name == "J05122_at")
exprs(E.GEOD.6353.eset)[row_idx, col_idx]


tmp <- data.frame(dataset_id = "E.GEOD.6353",
                  simulation = pData(E.GEOD.6353.eset)$Treatment[col_idx],
                  duration = "16h",
                  species = "Rat",
                  exprs = exprs(E.GEOD.6353.eset)[row_idx, col_idx],
                  logFC = "log2FC=1.53",
                  padj = "p.adj=0.0015")

dt <- rbind(dt, tmp)


load("multi_species_tspo_expression/E-GEOD-59263_rabbit/E.GEOD.59263.4h.eset.rdata")
dim(E.GEOD.59263.4h.eset)
col_idx <- which(pData(E.GEOD.59263.4h.eset)$FactorValue..STIMULATION. %in% 
                   c("control", "LPS"))
row_idx <- which(fData(E.GEOD.59263.4h.eset)$Reporter.Identifier == "A_04_P038882")
exprs(E.GEOD.59263.4h.eset)[row_idx, col_idx]


tmp <- data.frame(dataset_id = "E.GEOD.59263",
                  simulation = pData(E.GEOD.59263.4h.eset)$FactorValue..STIMULATION.[col_idx],
                  duration = "4h",
                  species = "Rabbit",
                  exprs = exprs(E.GEOD.59263.4h.eset)[row_idx, col_idx],
                  logFC = "log2FC=-0.94",
                  padj = "p.adj=0.139")

dt <- rbind(dt, tmp)

load("multi_species_tspo_expression/E-GEOD-59263_rabbit/E.GEOD.59263.24h.eset.rdata")
dim(E.GEOD.59263.24h.eset)
col_idx <- which(pData(E.GEOD.59263.24h.eset)$FactorValue..STIMULATION. %in% 
                   c("control", "LPS"))
row_idx <- which(fData(E.GEOD.59263.24h.eset)$Reporter.Identifier == "A_04_P038882")
exprs(E.GEOD.59263.24h.eset)[row_idx, col_idx]


tmp <- data.frame(dataset_id = "E.GEOD.59263",
                  simulation = pData(E.GEOD.59263.24h.eset)$FactorValue..STIMULATION.[col_idx],
                  duration = "24h",
                  species = "Rabbit",
                  exprs = exprs(E.GEOD.59263.24h.eset)[row_idx, col_idx],
                  logFC = "log2FC=-0.07",
                  padj = "p.adj=0.9318")

dt <- rbind(dt, tmp)

load("multi_species_tspo_expression/GSE71037_sheep/GSE71037_dds.rda")
dim(dds)
plotCounts(dds, gene = "PBR", intgroup = "treatment",
           normalized = T, transform = F, returnData = T)

col_idx <- which(colData(dds)$treatment %in% c("Saline", "LPS"))
row_idx <- which(rownames(dds) == "PBR")
counts(dds, normalized = T)[row_idx, col_idx]

tmp <- data.frame(dataset_id = "GSE71037",
                  simulation = colData(dds)$treatment[col_idx],
                  duration = "6h",
                  species = "Sheep",
                  exprs = log2(counts(dds, normalized = T)[row_idx, col_idx]),
                  logFC = "log2FC=0.08",
                  padj = "p.adj=0.9826")

dt <- rbind(dt, tmp)


load("multi_species_tspo_expression/E-GEOD-30956_pig/E-GEOD-30956_eset.rdata")
dim(E.GEOD.30956.eset)
col_idx <- which(pData(E.GEOD.30956.eset)$FactorValue..TIMEPOINT. %in% 
                   c("0h", "7h"))
row_idx <- which(fData(E.GEOD.30956.eset)$Composite.Element.Name == "Ssc.5936.1.A1_a_at")
exprs(E.GEOD.30956.eset)[row_idx, col_idx]


tmp <- data.frame(dataset_id = "E.GEOD.30956",
                  simulation = pData(E.GEOD.30956.eset)$FactorValue..TIMEPOINT.[col_idx],
                  duration = "7h",
                  species = "Pig",
                  exprs = exprs(E.GEOD.30956.eset)[row_idx, col_idx],
                  logFC = "log2FC=0.12",
                  padj = "p.adj=0.9993")

dt <- rbind(dt, tmp)

col_idx <- which(pData(E.GEOD.30956.eset)$FactorValue..TIMEPOINT. %in% 
                   c("0h", "24h"))
row_idx <- which(fData(E.GEOD.30956.eset)$Composite.Element.Name == "Ssc.5936.1.A1_a_at")
exprs(E.GEOD.30956.eset)[row_idx, col_idx]


tmp <- data.frame(dataset_id = "E.GEOD.30956",
                  simulation = pData(E.GEOD.30956.eset)$FactorValue..TIMEPOINT.[col_idx],
                  duration = "24h",
                  species = "Pig",
                  exprs = exprs(E.GEOD.30956.eset)[row_idx, col_idx],
                  logFC = "log2FC=0.04",
                  padj = "p.adj=1")

dt <- rbind(dt, tmp)

TSPO_expression <- dt %>%
  mutate(species = factor(species, 
                          levels = c("Human", "Pig", "Sheep", "Rabbit", "Mouse", "Rat" )),
         species_duration = paste(species, duration,
                                  sep = " "),
         species_duration = factor(species_duration, 
                                   levels = c("Human 6h",
                                              "Human 24h",
                                              "Pig 2h",
                                              "Pig 7h",
                                              "Pig 24h",
                                              "Cow 2h",
                                              "Sheep 6h",
                                              "Rabbit 4h",
                                              "Rabbit 24h",
                                              "Mouse 4h",
                                              "Mouse 24h",
                                              "Rat 2h",
                                              "Rat 4h",
                                              "Rat 8h",
                                              "Rat 16h"
                                              
                                   )),
         treatment = case_when(simulation %in% c("untreated", "none_24", "none.0hour",
                                                 "Control", "control", "Saline", "0h") ~ "Naive",
                               TRUE ~ "Stimulated"),
         group = paste(species_duration, treatment, sep = " "),
         
         treatment_duration = paste(treatment, duration, sep = " "),
         duration = factor(duration, levels = c("4h", "6h", "7h", "8h", "16h", "24h"))
  ) %>%
  group_by(group) %>%
  mutate(Mean = mean(exprs)) %>%
  ungroup %>%
  group_by(species_duration) %>%
  mutate(exprs_scaled = exprs/unique(Mean[treatment == "Naive"])) %>%
  ungroup() %>%
  mutate(label_sig = gsub("p.adj=", "", padj) %>% as.numeric()) %>%
  mutate_at(., "label_sig", ~case_when(. <= 0.001 ~ "***",
                                       . <= 0.01 ~ "**",
                                       . <= 0.05 ~ "*",
                                       . > 0.05 ~ " ")
  )


ann_text <- TSPO_expression %>%
  select(c(species, duration, treatment, exprs, exprs_scaled, label_sig)) %>%
  group_by(species, duration) %>%
  mutate(y_val = max(exprs_scaled)+ 0.05) %>%
  ungroup() %>%
  select(-c(exprs, exprs_scaled)) %>%
  distinct() %>%
  group_by(species) %>%
  mutate(x = as.numeric(factor(duration))) %>%
  ungroup() %>%
  mutate(xmin = ifelse(label_sig == " ", NA, x - 0.2),
         xmax = ifelse(label_sig == " ", NA, x + 0.2),
         group1 = "Naive",
         group2 = "Stimulated")


palette_choice <- paletteer::paletteer_d("ggsci::nrc_npg")


p1 <- ggplot(TSPO_expression)+
  geom_boxplot(aes(x = duration, y = exprs_scaled, fill = treatment),
               color = "black", outlier.colour = "white", width = 1)+
  geom_point(aes(x = duration, y = exprs_scaled, fill = treatment),
             position = position_dodge(1),
             shape = 21, size = 2, show.legend = F)+
  ylab("Normalized relative expression")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  xlab("")+
  scale_color_manual(name = NULL, values = palette_choice[c(2,1)], 
                     aesthetics = c("colour", "fill")) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 16),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        strip.text = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black", size = 18),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_blank())  +
  facet_grid(~species, scales = "free_x", space = "free") +
  scale_y_continuous(limits = c(0.5, 1.5), 
                     expand = c(0, 0)) +
  ggpubr::stat_pvalue_manual(data = ann_text, label = "label_sig", y.position = "y_val", tip.length = 0.01) 


ggsave("multi_species_tspo_expression/TSPO expression in multiple species box plot.png", 
       units = "in", width = 16, height = 4, dpi = 300, bg = "white")
