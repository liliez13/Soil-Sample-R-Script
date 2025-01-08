install.packages('ggpubr')
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(cowplot)
library(permute)
library(lattice)
library(vegan)
library(tidyr)
library(tidyverse)

ASV_highph_n <- read.csv("~/ASVtable_highph.csv", header=TRUE, row.names=1)
dist_matrix_highph <- avgdist(ASV_highph_n,sample=4000, dmethod = "bray")
soil_nmds_highph <- metaMDS(dist_matrix_highph) %>% scores() %>% as_tibble(rownames = "Sample")

ASV_lowph_n <- read.csv("~/ASVtable_lowph.csv", header=TRUE, row.names=1)
dist_matrix_lowph <- avgdist(ASV_lowph_n,sample=4000, dmethod = "bray")
soil_nmds_lowph <- metaMDS(dist_matrix_lowph) %>% scores() %>% as_tibble(rownames = "Sample")

ASV <- read.csv("~/ASVtable_rhizosbulk_forR_noCorM.csv", header=TRUE, row.names=1)
metadata <- read.csv("~/Lily_Rhizo_Metadata_modified.csv")
rhizosbulk_dis <- avgdist(ASV,sample=4000, dmethod = "bray")
nmds <- metaMDS(rhizosbulk_dis) %>% scores() %>% as_tibble(rownames = "Sample")

metadata_nmds <- inner_join(metadata, nmds)
metadata_nmds_high <- inner_join(metadata, soil_nmds_highph)
metadata_nmds_low <- inner_join(metadata, soil_nmds_lowph)

stress_value <- metaMDS(rhizosbulk_dis)$stress

stress_value_high <- metaMDS(dist_matrix_highph)$stress

stress_value_low <- metaMDS(dist_matrix_lowph)$stress


p1 <- metadata_nmds %>% ggplot(aes(x = NMDS1, y = NMDS2, color = Depth)) + 
  geom_point(shape = 19, size = 3) + 
  geom_text(aes(x = -0.37, y = 0.45, label = paste("Stress =", round(stress_value, 2))), size = 5, color = "black") +  # Add label at the center
  coord_fixed() + 
  scale_color_manual(name = NULL, 
                     values = c("violet", "blue", "orange"), 
                     breaks = c("Bulk", "Proximal", "Rhizosphere"), 
                     labels = c("Bulk", "Proximal", "Rhizosphere")) +
  labs(x = "", y = "Dimension 2", title = "All Samples (n=120)") + 
  theme_classic() + 
  theme(
    legend.text = element_text(size = 16),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "top",
    plot.title = element_text(size = 16)  # Adjust title size here
  )


p2 <- metadata_nmds_high %>% ggplot(aes(x = NMDS1, y = NMDS2, color = Depth)) + 
  geom_point(shape = 19, size = 3) + 
  geom_text(aes(x = -0.33, y = 0.45, label = paste("Stress =", round(stress_value_high, 2))), size = 5, color = "black") +  # Add label at the center
  coord_fixed() + 
  scale_color_manual(name = NULL, 
                     values = c("violet", "blue", "orange"), 
                     breaks = c("Bulk", "Proximal", "Rhizosphere"), 
                     labels = c("Bulk", "Proximal", "Rhizosphere")) +
  labs(x = "", y = "Dimension 2", title = "pH>5 (n=88)") + 
  theme_classic() + 
  theme(
    legend.text = element_text(size = 16),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "top",
    plot.title = element_text(size = 16)  # Adjust title size here
  )

p3 <- metadata_nmds_low %>% ggplot(aes(x = NMDS1, y = NMDS2, color = Depth)) + 
  geom_point(shape = 19, size = 3) + 
  geom_text(aes(x = -0.39, y = 0.45, label = paste("Stress =", round(stress_value_low, 2))), size = 5, color = "black") +  # Add label at the center
  coord_fixed() + 
  scale_color_manual(name = NULL, 
                     values = c("violet", "blue", "orange"), 
                     breaks = c("Bulk", "Proximal", "Rhizosphere"), 
                     labels = c("Bulk", "Proximal", "Rhizosphere")) +
  labs(x = "Dimension 1", y = "Dimension 2", title = "pH<5 (n=32)") + 
  theme_classic() + 
  theme(
    legend.text = element_text(size = 16),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "top",
    plot.title = element_text(size = 16)  # Adjust title size here
  )

multi_plot1<- ggarrange(p1, p2, p3, 
                        labels = c("A", "B", "C"), 
                        ncol = 1, nrow = 3, 
                        common.legend = T, 
                        font.label = list(size = 25),
                        heights = c(1, 1, 1))
multi_plot1

p1

