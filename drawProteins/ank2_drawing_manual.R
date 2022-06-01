# ============= #
# load all chunks up to line 61 in "ank2_drawing.Rmd" before continuing with this script

# add MBD region
ank2[286,] <- c("REGION", "Membrane binding domain", as.numeric(10), as.numeric(852), as.numeric(842), "Q01484","ANK2_HUMAN", as.integer(9606), as.integer(1))

# rownames(ank2[286,]) <- "featuresTemp.285" # doesn't work

ank2 %>% mutate(description = case_when(str_detect(description, "ANK *") ~ "ANK repeat",
                                        str_detect(description, "Repeat A*") ~ "Repeat A",
                                        TRUE ~ description)) -> ank2

# defining description factors
ank2$description <- factor(ank2$description, levels = c("Membrane binding domain","ANK repeat", "Interaction with SPTBN1", "ZU5 1", "ZU5 2", "UPA domain", "Death 1", "Repeat-rich region", "Repeat A", "Disordered", "Death 2"))

# type convert
ank2 -> ank2_copy

type.convert(ank2) -> ank2

# manual protein drawing
p + geom_rect(data = ank2[ank2$type == "CHAIN",], aes(xmin = begin, xmax = end,
                                                       ymin = order-0.025, ymax = order+0.025,
                                                      ), fill = "grey") + 
#  geom_rect(aes(xmin = 10, xmax = 852, ymin = 1-0.05, ymax = 1+0.05, fill = "MBD"), colour = "black") +
  geom_rect(data = ank2[ank2$type == "REGION",], aes(xmin = begin, xmax = end,
                                                     ymin = order-0.05, ymax = order+0.05, fill = description), colour = "black") +
  geom_rect(data = ank2[ank2$type == "DOMAIN",], aes(xmin = begin, xmax = end,
                                                        ymin = order-0.04, ymax = order+0.04,
                                                        fill = description), colour = "black") +
  geom_rect(data = ank2[ank2$type == "REPEAT",], aes(xmin = begin, xmax = end,
                                                     ymin = order-0.03, ymax = order+0.03, fill = description), colour = "black") +
  scale_fill_manual("", values = hue, labels = c("Membrane binding domain", "ANK repeat", "Interaction with SPTBN1", "ZU5 1", "ZU5 2", "UPA domain", "Death 1", "Repeat-rich region", "Repeat A", "Disordered", "Death 2"), breaks = c("Membrane binding domain","ANK repeat", "Interaction with SPTBN1", "ZU5 1", "ZU5 2", "UPA domain", "Death 1", "Repeat-rich region", "Repeat A", "Disordered", "Death 2")) +
  scale_x_continuous(limits = c(0, 4025)) +
  # S646F
  geom_segment(aes(x = 646, xend = 646, y = 1.03, yend = 1.15), size = 0.75) +
  geom_point(aes(x = 646, y = 1.16), size=3, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1.5) +
  # R990Q
  geom_segment(aes(x = 990, xend = 990, y = 1.03, yend = 1.15), size = 0.75) +
  geom_point(aes(x = 990, y = 1.16), size=3, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1.5) +
  # R990* - ASD
  geom_segment(aes(x = 990, xend = 990, y = 0.97, yend = 0.85), size = 0.75) +
  geom_point(aes(x = 990, y = 0.84), size=3, color="darkblue", fill=alpha("cyan", 0.3), alpha=0.7, shape=21, stroke=1.5) +
  # Q1283H
  geom_segment(aes(x = 1283, xend = 1283, y = 0.97, yend = 0.85), size = 0.75) +
  geom_point(aes(x = 1283, y = 0.84), size=3, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1.5) +   
  # T1404I - CV
  geom_segment(aes(x = 1404, xend = 1404, y = 1.03, yend = 1.15), size = 0.75) +
  geom_point(aes(x = 1404, y = 1.16), size=3, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1.5) +
  # E1458G
  geom_segment(aes(x = 1458, xend = 1458, y = 0.97, yend = 0.85), size = 0.75) +
  geom_point(aes(x = 1458, y = 0.84), size=3, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1.5) +
  # W1535R - CV
  geom_segment(aes(x = 1535, xend = 1535, y = 1.03, yend = 1.15), size = 0.75) +
  geom_point(aes(x = 1535, y = 1.16), size=3, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1.5) +
  # Q1589Kfs - ASD
  geom_segment(aes(x = 1589, xend = 1589, y = 0.99, yend = 0.85), size = 0.75) +
  geom_point(aes(x = 1589, y = 0.84), size=3, color="darkblue", fill=alpha("cyan", 0.3), alpha=0.7, shape=21, stroke=1.5) +
  # E1813K - CV
  geom_segment(aes(x = 1813, xend = 1813, y = 1.03, yend = 1.15), size = 0.75) +
  geom_point(aes(x = 1813, y = 1.16), size=3, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1.5) +
  # P1843S - ASD
  geom_segment(aes(x = 1843, xend = 1843, y = 0.97, yend = 0.85), size = 0.75) +
  geom_point(aes(x = 1843, y = 0.84), size=3, color="darkblue", fill=alpha("cyan", 0.3), alpha=0.7, shape=21, stroke=1.5) +
  # M1988T - CV
  geom_segment(aes(x = 1988, xend = 1988, y = 1.03, yend = 1.15), size = 0.75) +
  geom_point(aes(x = 1988, y = 1.16), size=3, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1.5) +
  # R2608fs - ASD
  geom_segment(aes(x = 2608, xend = 2608, y = 1.03, yend = 1.15), size = 0.75) +
  geom_point(aes(x = 2608, y = 1.16), size=3, color="darkblue", fill=alpha("cyan", 0.3), alpha=0.7, shape=21, stroke=1.5) +
  # E3429V - ASD
  geom_segment(aes(x = 3429, xend = 3429, y = 0.97, yend = 0.85), size = 0.75) +
  geom_point(aes(x = 3429, y = 0.84), size=3, color="darkblue", fill=alpha("cyan", 0.3), alpha=0.7, shape=21, stroke=1.5) +
  # I3437T - CV
  geom_segment(aes(x = 3437, xend = 3437, y = 1.03, yend = 1.15), size = 0.75) +
  geom_point(aes(x = 3437, y = 1.16), size=3, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1.5) +
  # V3634D - CV
  geom_segment(aes(x = 3634, xend = 3634, y = 1.03, yend = 1.15), size = 0.75) +
  geom_point(aes(x = 3634, y = 1.16), size=3, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1.5) +
  # T3744N - CV
  geom_segment(aes(x = 3744, xend = 3744, y = 0.99, yend = 0.85), size = 0.75) +
  geom_point(aes(x = 3744, y = 0.84), size=3, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1.5) +
  # R3906 - CV
  geom_segment(aes(x = 3906, xend = 3906, y = 1.01, yend = 1.15), size = 0.75) +
  geom_point(aes(x = 3906, y = 1.16), size=3, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1.5) +
  theme_classic() +
  theme(legend.position = "bottom", text = element_text(size = 20), panel.grid.major=element_blank()) +
  theme(axis.ticks = element_blank(),
        axis.text.y = element_blank()) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.line = element_blank()) -> ank2_plot

hue <- c("#a061df","#7ecf57","#c85194","#c9b34e","#43305a","#86caaf","#d3583b","#d1a29e","#4b613d","#727cd0","#d3583b")

ggsave("ank2_drawprotein.png",
       ank2_plot,
       height = 7,
       width = 12,
       dpi = 300,
       units = "in")
