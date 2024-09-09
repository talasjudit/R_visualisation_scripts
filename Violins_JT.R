if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
if(!require(rlang)) install.packages("rlang"); library(rlang)
if(!require(patchwork)) install.packages("patchwork"); library(patchwork)
#if(!require(scales)) install.packages("scales"); library(scales)
if(!require(data.table)) install.packages("data.table"); library(data.table)
if(!require(ggpmisc)) install.packages("ggpmisc"); library(ggpmisc)
if(!require(ggpubr)) install.packages("ggpubr"); library(ggpubr)
if(!require(RColorBrewer)) install.packages("RColorBrewer"); library(RColorBrewer)


widetolong <- function(input){
  pivot_longer(input, cols=c("Leaf","Root","Seedling","Tapetum","Meiocyte","Microspore","Sperm_cell","Sperm_nuc"), names_to="Tissue",values_to = "RPKM")
}



outliers <- function(x) {
  
  Q1 <- quantile(x, probs=.25)
  Q3 <- quantile(x, probs=.75)
  iqr = Q3-Q1
  
  upper_limit = Q3 + (iqr*1.5)
  lower_limit = Q1 - (iqr*1.5)
  
  x > upper_limit | x < lower_limit
}

remove_outliers <- function(df, cols = names(df)) {
  for (col in cols) {
    df <- df[!outliers(df[[col]]),]
  }
  df
}

colourCount = 8
getPalette = colorRampPalette(brewer.pal(8, "Paired"))

violinplot <- function(feature,
                       graph_titles,
                       ylabel){
  ggplot2::ggplot(feature, aes(x=Tissue, y=RPKM, fill=Tissue)) + 
    geom_violin(scale="width") + geom_boxplot(outlier.shape =NA,lwd=0.9,width=0.3) +
    scale_fill_manual(values = getPalette(colourCount)) + theme_classic()+
    labs(title=graph_titles, y = expression("24nt siRNA normalised to 21-nt miRNA "(log[2](RPKM+1)))) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 14, color = "black"),
          axis.text = element_text(size = 12, color = "black"),
          axis.line = element_line(size = 1, color = "black"),
          axis.ticks = element_line(size = 1, color = "black"),
          legend.text = element_text(size = 12, color = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30,hjust = 1),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
}

violinplot2 <- function(feature,
                        graph_titles){
  ggplot2::ggplot(feature, aes(x=Tissue, y=RPKM, fill=Tissue)) + 
    geom_violin(scale="width") + geom_boxplot(outlier.shape =NA,lwd=0.9,width=0.3) +
    scale_fill_brewer(palette=3,direction = 1)+theme_classic()+
    labs(title=graph_titles) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 14, color = "black"),
          axis.text = element_text(size = 12, color = "black"),
          axis.line = element_line(size = 1, color = "black"),
          axis.ticks = element_line(size = 1, color = "black"),
          legend.text = element_text(size = 12, color = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30,hjust = 1),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
}



data <- read_tsv("./unfiltered/log_RPKMs/cRdDMs/merged/cRdDM_WT_ALL_merged_24_log2_RPKM_miRNA_norm.tsv")


data_outliers <- remove_outliers(data, c("Leaf","Root","Seedling","Tapetum","Meiocyte","Microspore","Sperm_cell","Sperm_nuc"))
data_filtered <- widetolong(data_outliers)

#Features <- map(Features,function(i){i  %>% widetolong() %>%return()})

data_filtered$Tissue <- as.factor(data_filtered$Tissue)
data_filtered$Tissue <- factor(data_filtered$Tissue, levels = c("Leaf","Root","Seedling","Tapetum","Meiocyte","Microspore","Sperm_cell","Sperm_nuc"))


ggdata <- violinplot(data_filtered, graph_titles = "cRdDMs")


