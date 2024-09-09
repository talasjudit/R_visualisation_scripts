library(tidyverse)
locus_clusters <- read_tsv("top_MetGenes.tsv")
locus_clusters_mat <- data.matrix(locus_clusters)

install.packages("pheatmap")
library(pheatmap)

pheatmap(locus_clusters_mat,
         main = "MetGenes",
         fontsize = 12,
         color=colorRampPalette(c("#154c79", "white", "red"))(50),
         border_color = NA)


out <-pheatmap(locus_clusters_mat,
         main = "MetGenes rep2",
         fontsize = 12,
         color=colorRampPalette(c("#154c79", "white", "red"))(50)
)

heatmap.clust <- cbind(locus_clusters_mat, 
                      cluster = cutree(out$tree_row, k=4))

install.packages('xlsx')
library('xlsx')

write.xlsx(heatmap.clust, file="Hyper_TEs_4_clusters.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

d <- density(SC_rep1_clusters$SC)
plot(d)




