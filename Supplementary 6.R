library(ggplot2)
library(ggridges)
library(readxl)
library(tidyverse)
library(vapeplot)

GSEAanalysis_up <- read_excel("GSEAanalysis.xlsx", 
                           sheet = "aUP")

GSEA_data_up <- as.data.frame(GSEAanalysis_up)

ggplot(GSEA_data_up, aes(x = logFC, y = fct_rev(fct_reorder(GSEA, logFC*p_val, .fun = sum)), group = GSEA, fill = p_val)) +
  geom_density_ridges_gradient(scale = 0.8, rel_min_height = 0.01, jittered_points = TRUE,
                               position = position_points_jitter(width = 0.05, height = 0),
                               point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
  scale_fill_gradient(low = "#8795E8", high = "#FF6AD5", na.value = NA) +
  labs(title = 'Enriched KEGG defence-related gene sets')+
  theme(text = element_text(size=17))+
  ylab("")+
  geom_vline(xintercept = 0, linetype = "twodash")

