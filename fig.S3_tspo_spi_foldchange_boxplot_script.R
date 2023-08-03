library(ggplot2)
library(dplyr)
library(DESeq2)
library(greekLetters)


load("multi_species_tspo_expression/GSE36952_GSE66593_Human_macrophage_RNASeq/GSE36952_dds.RData")

dt <- read.delim("multi_species_tspo_expression/GSE36952_GSE66593_Human_macrophage_RNASeq/GSE36952.foldchnage.IFNgvsCT.txt")
dt %>%
  filter(GeneName %in% c("TSPO", "SPI1"))


col_idx <- which(colData(dds)$Treatment %in% c("None", "IFNg"))
row_idx <- which(rownames(dds) %in% c("TSPO"))
counts(dds, normalized = T)[row_idx, col_idx]

tspo <- data.frame(gene_name = "TSPO",
                  simulation = colData(dds)$Treatment[col_idx],
                  duration = "24h",
                  species = "Human",
                  exprs = log2(counts(dds, normalized = T)[row_idx, col_idx]),
                  logFC = "logFC=-0.41",
                  padj = "p.adj=1")

row_idx <- which(rownames(dds) %in% c("SPI1"))
counts(dds, normalized = T)[row_idx, col_idx]

spi1 <- data.frame(gene_name = "SPI1(PU.1)",
                   simulation = colData(dds)$Treatment[col_idx],
                   duration = "24h",
                   species = "Human",
                   exprs = log2(counts(dds, normalized = T)[row_idx, col_idx]),
                   logFC = "logFC=0.1",
                   padj = "p.adj=1")

dt <- rbind(tspo, spi1)

dt <- dt %>% 
  mutate(treatment = case_when(simulation %in% "None" ~ "Baseline",
                              simulation %in% "IFNg" ~ paste0("IFN", greek_vector["gamma"]))
         )
  

ann_text <- dt %>%
  mutate(label = sprintf("%s \n %s", logFC, padj)) %>%
  group_by(gene_name) %>%
  mutate(y_val = max(exprs) + 1) %>%
  ungroup() %>%
  select(-exprs) %>%
  distinct()


palette_choice <- paletteer::paletteer_d("ggsci::nrc_npg")


ggplot(dt)+
  geom_boxplot(aes(x = gene_name, y = exprs, fill = treatment),
               color = "black", outlier.colour = "white")+
  geom_point(aes(x = gene_name, y = exprs, fill = treatment),
             position = position_jitterdodge(0.1),
             shape = 21, size = 3, show.legend = F)+
  ylab("Normalized expression")+
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
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "transparent"))  
  # geom_text(data = ann_text, 
  #           aes(x = gene_name, y = y_val, label = label), 
  #           color = "black", show.legend = F, size = 6)  +
  # scale_y_continuous(limits = c(NA, max(dt$exprs+2))) 

ggsave("manuscript_fig_doc/TSPO_SPI1_boxplot_IFNg.png", 
       units = "in", width = 6, height = 4, dpi = 300, bg = "transparent")


####

GSE38371 <- read.table("multi_species_tspo_expression/GSE38371_mouse_IFNg_BMDM/GSM940702_IFNg_gene_exp.diff.txt", header = T)


GSE38371 %>%
  filter(test_id %in% c("ENSMUSG00000041736"))

mtspo <- data.frame(dataset_id = "GSE38371",
                   simulation = c("Baseline", "IFNg"),
                   duration = "4h",
                   species = "Mouse",
                   exprs = log2(c(GSE38371 %>%
                     filter(test_id %in% c("ENSMUSG00000041736")) %>%
                     pull(value_1),
                     GSE38371 %>%
                       filter(test_id %in% c("ENSMUSG00000041736")) %>%
                       pull(value_2))),
                   logFC = "logFC=1.28",
                   padj = "p.adj=0.1982")

col_idx <- which(colData(dds)$Treatment %in% c("None", "IFNg"))
row_idx <- which(rownames(dds) %in% c("TSPO"))
counts(dds, normalized = T)[row_idx, col_idx]

htspo <- data.frame(dataset_id = "GSE36952",
                   simulation = colData(dds)$Treatment[col_idx],
                   duration = "24h",
                   species = "Human",
                   exprs = log2(counts(dds, normalized = T)[row_idx, col_idx]),
                   logFC = "logFC=-0.41",
                   padj = "p.adj=1")

dt <- rbind(htspo, mtspo)
dt <- dt %>% 
  mutate(treatment = case_when(simulation %in% c("None", "Baseline") ~ "Baseline",
                               simulation %in% "IFNg" ~ paste0("IFN", greek_vector["gamma"]))
  )



ann_text <- dt %>%
  mutate(label = sprintf("%s \n %s", logFC, padj)) %>%
  group_by(dataset_id) %>%
  mutate(y_val = max(exprs) + 1) %>%
  ungroup() %>%
  select(-c(exprs, simulation)) %>%
  distinct()


palette_choice <- paletteer::paletteer_d("ggsci::nrc_npg")


ggplot(dt)+
  geom_boxplot(aes(x = species, y = exprs, fill = treatment),
               color = "black", outlier.colour = "white")+
  geom_point(aes(x = species, y = exprs, fill = treatment),
             position = position_jitterdodge(0.1),
             shape = 21, size = 3, show.legend = F)+
  ylab("Normalized expression")+
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
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "transparent"))  
  # geom_text(data = ann_text, 
  #           aes(x = species, y = y_val, label = label), 
  #           color = "black", show.legend = F, size = 6)  +
  # scale_y_continuous(limits = c(NA, max(dt$exprs+2))) 

ggsave("manuscript_fig_doc/TSPO_foldchange_IFNg_human_mouse_boxplot.png", 
       units = "in", width = 6, height = 4, dpi = 300, bg = "transparent")


