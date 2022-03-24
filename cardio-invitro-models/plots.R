library(ggplot2)
library(dplyr)
library(readr)

dge_table <- read.csv("cardio-invitro-models/celltype_groups_H9c2Undiff-RA.csv", sep = ",")

colnames(dge_table)[1] <- "gene.symbol"

dge_table %>% mutate(diffexpressed = case_when(log2FoldChange >= 1.25 & padj < 0.05 ~ "up",
                                               log2FoldChange <= -1.25 & padj < 0.05~ "down",
                                               padj >= 0.05 ~ "ns",
                                               log2FoldChange > -1.25 & log2FoldChange < 1.25 & padj < 0.05 ~ "nc")) %>% ggplot(aes(log2FoldChange, -log10(padj))) + 
  geom_point(aes(colour = diffexpressed), alpha = 0.2) + geom_hline(yintercept = -log10(0.05)) + 
  scale_colour_manual(name = "DiffExpr", values = c("blue", "grey", "#363636", "red"),
                      labels = c("down", "nc", "ns", "up")) +
  geom_vline(xintercept = c(1.25, -1.25)) +
  ggrepel::geom_label_repel(data = subset(dge_table, gene.symbol %in% ny_genes), 
                            aes(label = gene.symbol),
                            fontface = "italic",
                            segment.size = 0.8,
                            box.padding = unit(1, "lines"), 
                            point.padding = unit(0.3, "lines"),
                            max.overlaps = 20) + 
  theme_bw() +
  # scale_y_continuous(limits = c(min(-log10(dge_table$padj)), 200)) +
  scale_x_continuous(breaks = c(-5, 0, 5, 10, 15, 20)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted") +
  theme(legend.position = "bottom")

ny_genes <- c("Ank2", "Gata4", "Ppp1r9b", "Actn1","Actn2","Acta1", "Tnnt1", "Tnnt2","Obscn", "Ppp1r12b","Casq2", "Trdn" ,"Ryr2")
