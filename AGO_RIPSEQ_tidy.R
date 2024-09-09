# Functions to prepare raw data into format for plotting 
# ago_01 Calculates average of two reps x,y by position
# ago_02: 1. reformat positions relative to sRNA binding site 
#         2. create new column with sum count
#         3. normalisation e.g. A/(A+C+G+T)




ago_01 <- function(x,y){
  rbindlist(list(x,y))[,lapply(.SD,mean), list(pos)] }



ago_02 <- function(x){
  x %>%
    mutate(
      positions = pos - 60, pos=NULL) %>% 
    mutate(
      sums= A+G+C+T) %>%
    print() %>%
    mutate_at(vars(A,G,C,T), ~ (./sums)) %>%
    select(-sums)
}


#options for agoplot

line_colours <- c("#4B0055","#007094","#00BE7D","#FDE333")
graph_titles <- c("A", "G", "C", "T")
y_axis_labels <- c("ALL", "21", "22", "23", "24", "25")

#plot functions

agoplot_all_scales <- function(input,
                    which_column,
                    line_colours="black",
                    graph_titles=NULL,
                    y_axis_labels=NULL){
  ggplot2::ggplot(data=input, aes(x=positions, y={{which_column}}, group=1)) +
    geom_line(size=0.1, colour= line_colours) +
    theme(axis.line = element_line(size = 0.1, linetype = "solid"), 
          axis.title = element_text(size = 20) ,
          axis.title.y = element_text(angle = 0, vjust= 0.5),
          plot.title = element_text(size = 20, hjust = 0.5),
          panel.background = element_rect(fill = NA)) +
    labs(title = graph_titles, y = y_axis_labels , x=NULL) +
    scale_y_continuous(limits = c(0,1), labels = scales::number_format(accuracy = 0.01)) + geom_vline(xintercept = 24, linetype = "dotted", color = "blue", size = 1.5)
}

  

lines_24 <- c(-48,-24,24,48)
lines_23 <- c(-46,-23,23,46)
lines_22 <- c(-44,-22,22,44)
lines_21 <- c(-42,-21,21,42)
lines_19 <- c(-57,-38,-19,19,38,57)

agoplot_no_y_text <-  function(input,
                               which_column,
                               line_colours="black",
                               graph_titles=NULL,
                               y_axis_labels=NULL){
  ggplot2::ggplot(data=input, aes(x=positions, y={{which_column}}, group=1)) +
    geom_line(size=0.1, colour= line_colours) +
    stat_peaks(col = "black", geom = "text", size = 2, vjust = -0.5, ignore_threshold = 0.2 ) +
    stat_valleys(col = "red", geom = "text", size = 2, vjust = 0.5, ignore_threshold = 0.99) +
    theme(axis.line = element_line(size = 0.1, linetype = "solid"), 
          axis.title = element_text(size = 20) ,
          axis.title.y = element_text(angle = 0, vjust= 0.5),
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.text.y = element_blank(),
          panel.background = element_rect(fill = NA)) +
    labs(title = graph_titles, y = y_axis_labels , x=NULL) +
    scale_y_continuous(limits = c(0,1), labels = scales::number_format(accuracy = 1)) #+
    #geom_vline(xintercept = 0,linetype = "dotted",  color = "black", size = 1) #+
   # geom_vline(xintercept = lines_19,linetype = "dotted",  color = "black", size = 1) #+
    #geom_vline(xintercept = lines_21,linetype = "dotted",  color = "blue", size = 1) #+
    #geom_vline(xintercept = lines_22,linetype = "dotted",  color = "orange", size = 1) #+
    #geom_vline(xintercept = lines_23,linetype = "dotted",  color = "green", size = 1) #+
    #geom_vline(xintercept = lines_24,linetype = "dotted",  color = "purple", size = 1)
    }

#import raw data 

AGO_all_r1 <- read_tsv("AGO9_Flowerbud_RIP_rep1_sRNA_0_mismatch_21_25_mapped_sorted.bam_all_60bp.txt")
AGO_all_r2 <- read_tsv("AGO9_Flowerbud_RIP_rep2_sRNA_0_mismatch_21_25_mapped_sorted.bam_all_60bp.txt")

AGO_21_r1 <- read_tsv("AGO9_Flowerbud_RIP_rep1_sRNA_0_mismatch_21_25_mapped_sorted.bam_21_60bp.txt")
AGO_21_r2 <- read_tsv("AGO9_Flowerbud_RIP_rep2_sRNA_0_mismatch_21_25_mapped_sorted.bam_21_60bp.txt")

AGO_22_r1 <- read_tsv("AGO9_Flowerbud_RIP_rep1_sRNA_0_mismatch_21_25_mapped_sorted.bam_22_60bp.txt")
AGO_22_r2 <- read_tsv("AGO9_Flowerbud_RIP_rep2_sRNA_0_mismatch_21_25_mapped_sorted.bam_22_60bp.txt")

AGO_23_r1 <- read_tsv("AGO9_Flowerbud_RIP_rep1_sRNA_0_mismatch_21_25_mapped_sorted.bam_23_60bp.txt")
AGO_23_r2 <- read_tsv("AGO9_Flowerbud_RIP_rep2_sRNA_0_mismatch_21_25_mapped_sorted.bam_23_60bp.txt")

AGO_24_r1 <- read_tsv("AGO9_Flowerbud_RIP_rep1_sRNA_0_mismatch_21_25_mapped_sorted.bam_24_60bp.txt")
AGO_24_r2 <- read_tsv("AGO9_Flowerbud_RIP_rep2_sRNA_0_mismatch_21_25_mapped_sorted.bam_24_60bp.txt")

AGO_25_r1 <- read_tsv("AGO9_Flowerbud_RIP_rep1_sRNA_0_mismatch_21_25_mapped_sorted.bam_25_60bp.txt")
AGO_25_r2 <- read_tsv("AGO9_Flowerbud_RIP_rep2_sRNA_0_mismatch_21_25_mapped_sorted.bam_25_60bp.txt")



  

#Optional: merge reps
AGO_all <- ago_01(AGO_all_r1, AGO_all_r2)
AGO_21 <- ago_01(AGO_21_r1, AGO_21_r2)
AGO_22 <- ago_01(AGO_22_r1, AGO_22_r2)
AGO_23 <- ago_01(AGO_23_r1, AGO_23_r2)
AGO_24 <- ago_01(AGO_24_r1, AGO_24_r2)
AGO_25 <- ago_01(AGO_25_r1, AGO_25_r2)


#make list of data for map()

AGO_list <- list(AGO_all= AGO_all,
                 AGO_21 = AGO_21,
                 AGO_22 = AGO_22,
                 AGO_23 = AGO_23,
                 AGO_24 = AGO_24,
                 AGO_25 = AGO_25)


#wrangle data using functions defined above

AGO_list <- map(AGO_list, as.tibble)

#ago_sums <- function(x){
#  x %>%
#    mutate(
#      positions = pos - 60, pos=NULL) %>% 
#    mutate(
#      sums= A+G+C+T) #%>%
#    #summarise_at(c("A","G","C","T","sums"), mean)
#  }


#AGO_sums <- map(AGO_list,ago_sums)




AGO_list <- map(AGO_list, ago_02)



AGO_gg_A_plots <- map(AGO_list, agoplot_no_y_text,A,line_colours[1])
AGO_gg_G_plots <- map(AGO_list, agoplot_no_y_text,G,line_colours[2])
AGO_gg_C_plots <- map(AGO_list, agoplot_no_y_text,C,line_colours[3])
AGO_gg_T_plots <- map(AGO_list, agoplot_no_y_text,T,"red")


# Plot together


A_patch <- wrap_plots(AGO_gg_A_plots, ncol=1)
G_patch <- wrap_plots(AGO_gg_G_plots, ncol=1)
C_patch <- wrap_plots(AGO_gg_C_plots, ncol=1)
T_patch <- wrap_plots(AGO_gg_T_plots, ncol=1)

A_patch| G_patch | C_patch | T_patch

