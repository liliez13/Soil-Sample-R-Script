library(readxl)
library(ggtext)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)

ASV <- read.csv("~/ASVtable_rhizosbulk_forR_noCorM.csv", header=TRUE) %>%
  pivot_longer(-Sample, names_to = "ASV", values_to = "count")

metadata <- read.csv("~/Lily_Rhizo_Metadata_modified.csv") %>%
  select(Sample, Grassland, Depth, Date, Year, SiteDepth, Depth2)

taxonomy <- read_excel("~/ASV_taxonomy_rhizosbulk_noCorM.xlsx")

ASV_rel_abund <- inner_join(metadata, ASV, by="Sample") %>% inner_join(., taxonomy, by="ASV") %>%
  group_by(Sample) %>% mutate(rel_abund = count/sum(count)) %>%
  ungroup() %>% select(-count) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
               names_to = "Level", values_to = "Taxon") %>%
  mutate(SiteDepth = factor(SiteDepth, levels=c("Ricketts Bulk", "Nescopeck Bulk", "Ricketts Proximal", 
                                                "Nescopeck Proximal", "Ricketts Rhizosphere", "Nescopeck Rhizosphere")))

taxon_rel_abund <- ASV_rel_abund %>%
  filter(Level=="Phylum") %>% 
  group_by(SiteDepth, Sample, Taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
  group_by(SiteDepth, Taxon) %>%
  summarize(mean_rel_abund =100*mean(rel_abund), .groups = "drop") %>%
  mutate(Taxon =str_replace(Taxon,"NA", "Unclassified"), Taxon = str_replace(Taxon, "^(\\S*)$","*\\1*"))


taxon_pool <- taxon_rel_abund %>% group_by(Taxon) %>% summarize(pool=max(mean_rel_abund) < 4, .groups="drop")

inner_join(taxon_rel_abund, taxon_pool, by="Taxon") %>% mutate(Taxon=if_else(pool, "Other", Taxon)) %>%
  filter(!is.na(SiteDepth)) %>%
  group_by(SiteDepth, Taxon) %>%
  summarize(mean_rel_abund =sum(mean_rel_abund), .groups = "drop") %>%
  ggplot(aes(x=SiteDepth, y=mean_rel_abund, fill=Taxon)) + 
  geom_col() + 
  scale_fill_manual(name=NULL, breaks=c("*Acidobacteriota*", "*Actinobacteriota*", "*Bacteroidota*", 
                                        "*Planctomycetota*", "*Proteobacteria*", 
                                        "*Verrucomicrobiota*", "Other"),
                    values=c(brewer.pal(6, "Paired"), "gray")) + 
  scale_y_continuous(expand=c(0,0)) +
  labs(x=NULL, y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(legend.text= element_markdown(),
        legend.key.size = unit(10, "pt"))

taxon_rel_abund_order <- ASV_rel_abund %>%
  filter(Level=="Order") %>% 
  group_by(SiteDepth, Sample, Taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
  group_by(SiteDepth, Taxon) %>%
  summarize(mean_rel_abund =100*mean(rel_abund), .groups = "drop") %>%
  mutate(Taxon =str_replace(Taxon,"NA", "Unclassified"), Taxon = str_replace(Taxon, "^(\\S*)$","*\\1*"))

taxon_pool_order <- taxon_rel_abund_order %>% group_by(Taxon) %>% summarize(pool=max(mean_rel_abund) < 1.75, .groups="drop")

inner_join(taxon_rel_abund_order, taxon_pool_order, by = "Taxon") %>%
  mutate(
    Taxon = if_else(pool, "Other", Taxon),
    Taxon = factor(
      Taxon,
      levels = c(
        "*Acetobacterales*",
        "*Acidobacteriales*", "*Burkholderiales*", "*Chitinophagales*", "*Chthoniobacterales*", "*Elsterales*",
        "*Gemmatales*", "*Isosphaerales*", "*Pedosphaerales*", "*Pirellulales*", "*Planctomycetales*", "*Polyangiales*", "*Pseudomonadales*", "*Tepidisphaerales*",
        "*Rhizobiales*","*Sphingobacteriales*", "*Vicinamibacterales*", "Subgroup 2", "*Unclassified*", "Other"
      )
    )
  ) %>%
  filter(!is.na(SiteDepth)) %>%
  group_by(SiteDepth, Taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund), .groups = "drop") %>%
  ggplot(aes(x = SiteDepth, y = mean_rel_abund, fill = Taxon)) +
  geom_col() +
  scale_fill_manual(
    name = NULL,
    breaks = c(
      "*Acetobacterales*",
      "*Acidobacteriales*", "*Burkholderiales*", "*Chitinophagales*", "*Chthoniobacterales*", "*Elsterales*",
      "*Gemmatales*", "*Isosphaerales*", "*Pedosphaerales*", "*Pirellulales*", "*Planctomycetales*", "*Polyangiales*", "*Pseudomonadales*", "*Tepidisphaerales*",
      "*Rhizobiales*","*Sphingobacteriales*", "*Vicinamibacterales*", "Subgroup 2", "*Unclassified*", "Other"
    ),
    values = c(brewer.pal(12, "Paired"), "blue", "green", "cornflowerblue", "lightgreen", "yellow", "pink", "black", "gray")
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Mean Relative Abundance (%)") +
  theme_classic() +
  theme(
    legend.text = element_markdown(),
    legend.key.size = unit(10, "pt")
  )

library(ggplot2)
library(ggthemes)
library(ggpubr)
library(patchwork)
library(gridExtra)
library(stringr)

P1 <- inner_join(taxon_rel_abund, taxon_pool, by = "Taxon") %>%
  mutate(
    Taxon = if_else(pool, "Other", Taxon)
  ) %>%
  filter(!is.na(SiteDepth)) %>%
  group_by(SiteDepth, Taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund), .groups = "drop") %>%
  ggplot(aes(x = SiteDepth, y = mean_rel_abund, fill = Taxon)) + 
  geom_col() + 
  scale_fill_manual(
    name = NULL, 
    breaks = c("*Acidobacteriota*", "*Actinobacteriota*", "*Bacteroidota*", 
               "*Planctomycetota*", "*Proteobacteria*", 
               "*Verrucomicrobiota*", "Other"),
    values = c(brewer.pal(6, "Paired"), "gray")
  ) + 
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Mean Relative Abundance (%)") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8))+
  theme_classic() +
  theme(
    legend.text = element_markdown(),
    legend.key.size = unit(15, "pt")
  ) +
  guides(fill = guide_legend(title = NULL)) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),  # Adjust the size parameter as needed
        axis.line = element_line(size = 1.5))
P1

P2 <- inner_join(taxon_rel_abund_order, taxon_pool_order, by = "Taxon") %>%
  mutate(
    Taxon = if_else(pool, "Other", Taxon),
    Taxon = factor(
      Taxon,
      levels = c(
        "*Acetobacterales*", "*Acidobacteriales*", "*Burkholderiales*", "*Chitinophagales*", "*Chthoniobacterales*", "*Elsterales*",
        "*Gemmatales*", "*Isosphaerales*", "*Pedosphaerales*", "*Pirellulales*", "*Planctomycetales*", "*Polyangiales*", "*Pseudomonadales*", "*Tepidisphaerales*",
        "*Rhizobiales*", "*Sphingobacteriales*", "*Vicinamibacterales*", "Subgroup 2", "*Unclassified*", "Other"
      )
    )
  ) %>%
  filter(!is.na(SiteDepth)) %>%
  group_by(SiteDepth, Taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund), .groups = "drop") %>%
  ggplot(aes(x = SiteDepth, y = mean_rel_abund, fill = Taxon)) +
  geom_col() +
  scale_fill_manual(
    name = NULL,
    breaks = c(
      "*Acetobacterales*", "*Acidobacteriales*", "*Burkholderiales*", "*Chitinophagales*", "*Chthoniobacterales*", "*Elsterales*",
      "*Gemmatales*", "*Isosphaerales*", "*Pedosphaerales*", "*Pirellulales*", "*Planctomycetales*", "*Polyangiales*", "*Pseudomonadales*", "*Tepidisphaerales*",
      "*Rhizobiales*", "*Sphingobacteriales*", "*Vicinamibacterales*", "Subgroup 2", "*Unclassified*", "Other"
    ),
    values = c(brewer.pal(12, "Paired"), "blue", "green", "cornflowerblue", "lightgreen", "yellow", "pink", "black", "gray")
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Mean Relative Abundance (%)") + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8))+
  theme_classic() +
  theme(
    legend.text = element_markdown(),
    legend.key.size = unit(15, "pt")
  ) +
  guides(fill = guide_legend(title = NULL)) +
  theme(axis.text = element_text(size = 15), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),  # Adjust the size parameter as needed
        axis.line = element_line(size = 1.5)) 

multi_plot <- ggarrange(P1, P2, 
                        labels = c("A", "B"), 
                        ncol = 1, nrow = 2, 
                        common.legend = FALSE,
                        font.label = list(size = 25))

multi_plot

##By Site

ASV_rel_abund_site <- inner_join(metadata, ASV, by="Sample") %>% inner_join(., taxonomy, by="ASV") %>%
  group_by(Sample) %>% mutate(rel_abund = count/sum(count)) %>%
  ungroup() %>% select(-count) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
               names_to = "Level", values_to = "Taxon") %>%
  mutate(Depth2 = factor(Depth2, levels=c("Ricketts Site 1 Rhizosphere", "Ricketts Site 1 Proximal", "Ricketts Site 1 Bulk", "Ricketts Site 2 Rhizosphere", "Ricketts Site 2 Proximal", 
                                          "Ricketts Site 2 Bulk", "Nescopeck Site 1 Rhizosphere", "Nescopeck Site 1 Proximal", "Nescopeck Site 1 Bulk",
                                                "Nescopeck Site 3 Rhizosphere", "Nescopeck Site 3 Proximal", "Nescopeck Site 3 Bulk")))

taxon_rel_abund_site <- ASV_rel_abund_site %>%
  filter(Level=="Phylum") %>% 
  group_by(Depth2, Sample, Taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
  group_by(Depth2, Taxon) %>%
  summarize(mean_rel_abund =100*mean(rel_abund), .groups = "drop") %>%
  mutate(Taxon =str_replace(Taxon,"NA", "Unclassified"), Taxon = str_replace(Taxon, "^(\\S*)$","*\\1*"))


taxon_pool_site <- taxon_rel_abund_site %>% group_by(Taxon) %>% summarize(pool=max(mean_rel_abund) < 4, .groups="drop")

inner_join(taxon_rel_abund_site, taxon_pool_site, by="Taxon") %>% mutate(Taxon=if_else(pool, "Other", Taxon)) %>%
  filter(!is.na(Depth2)) %>%
  group_by(Depth2, Taxon) %>%
  summarize(mean_rel_abund =sum(mean_rel_abund), .groups = "drop") %>%
  ggplot(aes(x=Depth2, y=mean_rel_abund, fill=Taxon)) + 
  geom_col() + 
  scale_fill_manual(name=NULL, breaks=c("*Acidobacteriota*", "*Actinobacteriota*", "*Bacteroidota*", "*Myxococcota*",
                                        "*Planctomycetota*", "*Proteobacteria*", 
                                        "*Verrucomicrobiota*", "Other"),
                    values=c(brewer.pal(7, "Paired"), "gray")) + 
  scale_y_continuous(expand=c(0,0)) +
  labs(x=NULL, y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(legend.text= element_markdown(),
        legend.key.size = unit(10, "pt"))

taxon_rel_abund_order_site <- ASV_rel_abund_site %>%
  filter(Level=="Order") %>% 
  group_by(Depth2, Sample, Taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
  group_by(Depth2, Taxon) %>%
  summarize(mean_rel_abund =100*mean(rel_abund), .groups = "drop") %>%
  mutate(Taxon =str_replace(Taxon,"NA", "Unclassified"), Taxon = str_replace(Taxon, "^(\\S*)$","*\\1*"))

taxon_pool_order_site <- taxon_rel_abund_order_site %>% group_by(Taxon) %>% summarize(pool=max(mean_rel_abund) < 2.75, .groups="drop")

inner_join(taxon_rel_abund_order_site, taxon_pool_order_site, by = "Taxon") %>%
  mutate(
    Taxon = if_else(pool, "Other", Taxon),
    Taxon = factor(
      Taxon,
      levels = c(
        "*Acetobacterales*",
        "*Acidobacteriales*", "*Burkholderiales*", "*Chitinophagales*", "*Chthoniobacterales*", "*Elsterales*",
        "*Gemmatales*", "*Isosphaerales*", "*Pedosphaerales*", "*Pirellulales*", "*Planctomycetales*", "*Polyangiales*", "*Pseudomonadales*", "*Tepidisphaerales*",
        "*Rhizobiales*","*Sphingobacteriales*", "*Vicinamibacterales*", "*WD260*", "Subgroup 2", "*Unclassified*", "Other"
      )
    )
  ) %>%
  filter(!is.na(Depth2)) %>%
  group_by(Depth2, Taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund), .groups = "drop") %>%
  ggplot(aes(x = Depth2, y = mean_rel_abund, fill = Taxon)) +
  geom_col() +
  scale_fill_manual(
    name = NULL,
    breaks = c(
      "*Acetobacterales*",
      "*Acidobacteriales*", "*Burkholderiales*", "*Chitinophagales*", "*Chthoniobacterales*", "*Elsterales*",
      "*Gemmatales*", "*Isosphaerales*", "*Pedosphaerales*", "*Pirellulales*", "*Planctomycetales*", "*Polyangiales*", "*Pseudomonadales*", "*Tepidisphaerales*",
      "*Rhizobiales*","*Sphingobacteriales*", "*Vicinamibacterales*", "*WD260*", "Subgroup 2", "*Unclassified*", "Other"
    ),
    values = c(brewer.pal(12, "Paired"), "blue", "green", "cornflowerblue", "lightgreen", "yellow", "pink", "gold", "black", "gray")
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Mean Relative Abundance (%)") +
  theme_classic() +
  theme(
    legend.text = element_markdown(),
    legend.key.size = unit(10, "pt")
  )

library(ggplot2)
library(ggthemes)
library(ggpubr)
library(patchwork)
library(gridExtra)
library(stringr)

P3 <- inner_join(taxon_rel_abund_site, taxon_pool_site, by = "Taxon") %>%
  mutate(
    Taxon = if_else(pool, "Other", Taxon)
  ) %>%
  filter(!is.na(Depth2)) %>%
  group_by(Depth2, Taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund), .groups = "drop") %>%
  ggplot(aes(x = Depth2, y = mean_rel_abund, fill = Taxon)) + 
  geom_col() + 
  scale_fill_manual(
    name = NULL, 
    breaks = c("*Acidobacteriota*", "*Actinobacteriota*", "*Bacteroidota*", "*Myxococcota*",
               "*Planctomycetota*", "*Proteobacteria*", 
               "*Verrucomicrobiota*", "Other"),
    values = c(brewer.pal(7, "Paired"), "gray")
  ) + 
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Mean Relative Abundance (%)") + 
  scale_x_discrete(limits = c("Ricketts Site 1 Rhizosphere", "Ricketts Site 1 Proximal", "Ricketts Site 1 Bulk", "trt99",
                              "Ricketts Site 2 Rhizosphere", "Ricketts Site 2 Proximal", "Ricketts Site 2 Bulk", "trt98",
                              "Nescopeck Site 1 Rhizosphere", "Nescopeck Site 1 Proximal", "Nescopeck Site 1 Bulk", "trt96",
                              "Nescopeck Site 3 Rhizosphere", "Nescopeck Site 3 Proximal", "Nescopeck Site 3 Bulk"),
                   breaks = c("Ricketts Site 1 Rhizosphere", "Ricketts Site 1 Proximal", "Ricketts Site 1 Bulk",       NA,
                              "Ricketts Site 2 Rhizosphere", "Ricketts Site 2 Proximal", "Ricketts Site 2 Bulk",       NA,
                              "Nescopeck Site 1 Rhizosphere", "Nescopeck Site 1 Proximal", "Nescopeck Site 1 Bulk",       NA,
                              "Nescopeck Site 3 Rhizosphere", "Nescopeck Site 3 Proximal", "Nescopeck Site 3 Bulk"),
                   labels = function(x) str_wrap(x, width = 6))+
  theme_classic() +
  theme(
    legend.text = element_markdown(),
    legend.key.size = unit(15, "pt")
  ) +
  guides(fill = guide_legend(title = NULL)) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(size = 6),
        axis.title.y = element_text(size = 20),  # Adjust the size parameter as needed
        axis.line = element_line(linewidth = 1.5))

P4 <- inner_join(taxon_rel_abund_order_site, taxon_pool_order_site, by = "Taxon") %>%
  mutate(
    Taxon = if_else(pool, "Other", Taxon),
    Taxon = factor(
      Taxon,
      levels = c(
        "*Acetobacterales*",
        "*Acidobacteriales*", "*Burkholderiales*", "*Chitinophagales*", "*Chthoniobacterales*", "*Elsterales*",
        "*Gemmatales*", "*Isosphaerales*", "*Pedosphaerales*", "*Pirellulales*", "*Planctomycetales*", "*Polyangiales*", "*Pseudomonadales*", "*Tepidisphaerales*",
        "*Rhizobiales*","*Sphingobacteriales*", "*Vicinamibacterales*", "*WD260*", "Subgroup 2", "*Unclassified*", "Other"
      )
    )
  ) %>%
  filter(!is.na(Depth2)) %>%
  group_by(Depth2, Taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund), .groups = "drop") %>%
  ggplot(aes(x = Depth2, y = mean_rel_abund, fill = Taxon)) +
  geom_col() +
  scale_fill_manual(
    name = NULL,
    breaks = c(
      "*Acetobacterales*",
      "*Acidobacteriales*", "*Burkholderiales*", "*Chitinophagales*", "*Chthoniobacterales*", "*Elsterales*",
      "*Gemmatales*", "*Isosphaerales*", "*Pedosphaerales*", "*Pirellulales*", "*Planctomycetales*", "*Polyangiales*", "*Pseudomonadales*", "*Tepidisphaerales*",
      "*Rhizobiales*","*Sphingobacteriales*", "*Vicinamibacterales*", "*WD260*", "Subgroup 2", "*Unclassified*", "Other"
    ),
    values = c(brewer.pal(12, "Paired"), "blue", "green", "cornflowerblue", "lightgreen", "yellow", "pink","gold", "black", "gray")
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Mean Relative Abundance (%)") + 
  scale_x_discrete(limits = c("Ricketts Site 1 Rhizosphere", "Ricketts Site 1 Proximal", "Ricketts Site 1 Bulk", "trt99",
                              "Ricketts Site 2 Rhizosphere", "Ricketts Site 2 Proximal", "Ricketts Site 2 Bulk", "trt98",
                              "Nescopeck Site 1 Rhizosphere", "Nescopeck Site 1 Proximal", "Nescopeck Site 1 Bulk", "trt96",
                              "Nescopeck Site 3 Rhizosphere", "Nescopeck Site 3 Proximal", "Nescopeck Site 3 Bulk"),
                   breaks = c("Ricketts Site 1 Rhizosphere", "Ricketts Site 1 Proximal", "Ricketts Site 1 Bulk",       NA,
                              "Ricketts Site 2 Rhizosphere", "Ricketts Site 2 Proximal", "Ricketts Site 2 Bulk",       NA,
                              "Nescopeck Site 1 Rhizosphere", "Nescopeck Site 1 Proximal", "Nescopeck Site 1 Bulk",       NA,
                              "Nescopeck Site 3 Rhizosphere", "Nescopeck Site 3 Proximal", "Nescopeck Site 3 Bulk"),
                   labels = function(x) str_wrap(x, width = 6))+
  theme_classic() +
  theme(
    legend.text = element_markdown(),
    legend.key.size = unit(15, "pt")
  ) +
  guides(fill = guide_legend(title = NULL)) +
  theme(axis.text = element_text(size = 15), 
        axis.text.x = element_text(size = 6),
        axis.title.y = element_text(size = 20),  # Adjust the size parameter as needed
        axis.line = element_line(linewidth = 1.5))
P4

multi_plot <- ggarrange(P3, P4, 
                        labels = c("A", "B"), 
                        ncol = 1, nrow = 2, 
                        common.legend = FALSE,
                        font.label = list(size = 25))

multi_plot





