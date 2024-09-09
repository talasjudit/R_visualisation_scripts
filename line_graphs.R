library(tidyverse)
library(RColorBrewer)
library(viridis)
burst_peak <- read.csv("Burst_peak.csv")

burst_peak <- burst_peak %>%
  mutate(genotype = recode(genotype,
                           `WT x WT` = "WT",
                           `WT x dn4mt1 ko #6` = "dn4mt1_6_KO",
                           `WT x dn4mt1 ko #21` = "dn4mt1_21_KO",
                           `WT x dnmt3b ko D5`="dnmt3b_KO",
                           `WT x cmta ko #64`="cmta_KO",
                           .default = genotype))


burst_peak_long <- burst_peak %>%
  pivot_longer(!genotype,names_to = "days_since_first_mature_sporophyte", values_to = "percent")

burst_peak_long<- burst_peak_long %>%
  mutate(days_since_first_mature_sporophyte = stringr::str_replace(days_since_first_mature_sporophyte, "X", "")) %>% 
  mutate(days_since_first_mature_sporophyte = as.numeric(days_since_first_mature_sporophyte)) %>%
  mutate(genotype = factor(genotype, levels= c("WT","dn4mt1_6_KO","dn4mt1_21_KO","dnmt3b_KO","cmta_KO")))

colour_mapping <- c("WT" = rgb(122/255,151/255,243/255), 
                   "dn4mt1_6_KO" = rgb(18/255,163/255,125/255), 
                   "dn4mt1_21_KO" = rgb(18/255,163/255,125/255), 
                   "dnmt3b_KO" = rgb(168/255,55/255,23/255), 
                   "cmta_KO" = rgb(168/255,55/255,23/255))

readable_labels <- c("WT" = "Wild Type", 
                     "dn4mt1_6_KO" = "dn4mt1 KO #6", 
                     "dn4mt1_21_KO" = "dn4mt1 KO #21", 
                     "dnmt3b_KO" = "dnmt3b KO", 
                     "cmta_KO" = "cmta KO")



burst_peak_filter <- burst_peak_long %>%
  filter(genotype %in% c("WT","dn4mt1_6_KO","dn4mt1_21_KO","dnmt3b_KO","cmta_KO"))

burst_peak_filter %>%
  ggplot(aes(x = days_since_first_mature_sporophyte, y = percent,group = genotype, color = genotype)) +
  geom_ribbon(aes(ymin = 0, ymax = percent, fill = genotype), alpha = 0.5) +
  scale_color_manual(values = colour_mapping, name = "sperm genotype", labels = readable_labels) +
  scale_fill_manual(values = colour_mapping, name = "sperm genotype", labels = readable_labels) + 
  theme_classic() +
  scale_x_continuous(limits = c(0, 26), breaks = c(0, 5, 10, 15, 20, 25)) +
  scale_y_continuous(limits = c(0,36), breaks = c(0, 5, 10, 15, 20, 25,30,35)) +
  labs(x = "days since first mature sporophyte", 
       y = " % burst sporophytes",colour = "sperm genotype") +
  theme(legend.position = "none")  
