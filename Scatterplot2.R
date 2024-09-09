if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
if(!require(rlang)) install.packages("rlang"); library(rlang)
if(!require(patchwork)) install.packages("patchwork"); library(patchwork)
#if(!require(scales)) install.packages("scales"); library(scales)
if(!require(data.table)) install.packages("data.table"); library(data.table)
if(!require(ggpmisc)) install.packages("ggpmisc"); library(ggpmisc)
if(!require(ggpubr)) install.packages("ggpubr"); library(ggpubr)


MC_rep1 <- read_tsv('./24v25/bowtie/WT_meiocyte_rep1_24_25_log10RPKM1.tsv')
MC_rep2 <- read_tsv('./24v25/bowtie/WT_meiocyte_rep2_24_25_log10RPKM1.tsv')
MS_rep1 <- read_tsv('./24v25/bowtie/WT_microspore_rep1_24_25_log10RPKM1.tsv')
MS_rep2 <- read_tsv('./24v25/bowtie/WT_microspore_rep2_24_25_log10RPKM1.tsv')
TAP_rep1 <- read_tsv('./24v25/bowtie/WT_tapetum_rep1_24_25_log10RPKM1.tsv')
TAP_rep2 <- read_tsv('./24v25/bowtie/WT_tapetum_rep2_24_25_log10RPKM1.tsv')



scatterplot <- function(input,
                        column_1,
                        column_2,
                        graph_titles=NULL){
  ggplot2::ggplot(data=input, aes(x={{column_1}}, y={{column_2}})) +
  geom_point(size= 0.1, alpha=1) +
  scale_fill_brewer(palette="Paired")+theme_classic()+
  theme(axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.line = element_line(size = 1, color = "black"),
        axis.ticks = element_line(size = 1, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5)) +
           labs(title=graph_titles, x=expression("24nt siRNA "(log[10](RPKM+1))), y=expression("25nt siRNA "(log[10](RPKM+1))))  +
    scale_x_continuous(limits=c(0,4, breaks=seq(0,0.5,4))) +
    scale_y_continuous(limits=c(0,4, breaks=seq(0,0.5,4)))
}

scatterplot <- function(input,
                              column_1,
                              column_2,
                              graph_titles=NULL){
  ggplot2::ggplot(data=input, aes(x={{column_1}}, y={{column_2}})) +
    geom_point(size= 0.1, alpha=1) +
    scale_fill_brewer(palette="Paired")+theme_classic()+
    theme(axis.title = element_text(size = 12, color = "black"),
          axis.text = element_text(size = 12, color = "black"),
          axis.line = element_line(size = 1, color = "black"),
          axis.ticks = element_line(size = 1, color = "black"),
          legend.title = element_blank(),
          legend.text = element_text(size = 12, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    coord_cartesian(ylim = c(0,4), xlim = c(0,4)) +
    labs(title=graph_titles, x=expression("Rep 1 24nt siRNA "(log[10](RPKM+1))), y=expression("Rep 2 24nt siRNA "(log[10](RPKM+1)))) +
    stat_cor(method = "pearson", aes(label = ..r.label..))
}
    # +
    #scale_x_continuous(limits=c(0,15, breaks=seq(0,20,5))) +
    #scale_y_continuous(limits=c(0,15, breaks=seq(0,20,5)))

# + geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95)

Meiocyte_rep1_gg <- scatterplot(MC_rep1,RPKM_24,RPKM_25, graph_titles = "Meiocyte rep1")
Meiocyte_rep2_gg <- scatterplot(MC_rep2,RPKM_24,RPKM_25, graph_titles = "Meiocyte rep2")
Microspore_rep1_gg <- scatterplot(MS_rep1,RPKM_24,RPKM_25, graph_titles = "Microspore rep1")
Microspore_rep2_gg <- scatterplot(MS_rep2,RPKM_24,RPKM_25, graph_titles = "Microspore rep2")
Tapetum_rep1_gg <- scatterplot(TAP_rep1,RPKM_24,RPKM_25, graph_titles = "Tapetum rep1")
Tapetum_rep2_gg <- scatterplot(TAP_rep2,RPKM_24,RPKM_25, graph_titles = "Tapetum rep2")

(Microspore_rep1_gg/Microspore_rep2_gg) | (Meiocyte_rep1_gg/Meiocyte_rep2_gg) | (Tapetum_rep1_gg/Tapetum_rep2_gg)

(Tapetum_rep1_gg | Tapetum_rep2_gg ) / (Meiocyte_rep1_gg | Meiocyte_rep2_gg) / (Microspore_rep1_gg | Microspore_rep2_gg)


## merged between reps

TAP_24_merge <- read_tsv("WT_tapetum_log10RPKH1_merged.tsv")
MC_24_merge <- read_tsv("WT_meiocyte_log10RPKH1_merged.tsv")
MS_24_merge <- read_tsv("./rep_comparison/Jinchengs_method/WT_microspore_log10_RPKM_merged.tsv")

TAP_merge_gg <- scatterplot(TAP_24_merge, Rep1_RPKM, Rep2_RPKM, graph_titles = "Tapetum")
MC_merge_gg <- scatterplot(MC_24_merge, Rep1_RPKM, Rep2_RPKM, graph_titles = "Meiocyte")
MS_merge_gg <- scatterplot(MS_24_merge, Rep1_RPKM, Rep2_RPKM, graph_titles = "Microspore")

TAP_merge_gg + annotate("text", x = 3.5, y = 0, label = "n  == 8513",parse = TRUE) | 
MC_merge_gg + annotate("text", x = 3.5, y = 0, label = "n  == 4843",parse = TRUE)| 
MS_merge_gg + annotate("text", x = 3.5, y = 0, label = "n  == 11748",parse = TRUE)


