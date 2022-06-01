library(ggplot2)
library(dplyr)
library(readr)

expr_wga <- read.table("human_gestational_expr/human_gest_Cardiac_expr_rld.csv", sep = ",", header = TRUE)

colnames(expr_wga)[1] <- "Gene"

expr_wga %>%  mutate(Gene = case_when(
  Gene == "ENSG00000145362" ~ "ANK2",
  Gene == "ENSG00000152661" ~ "GJA1",
  Gene == "ENSG00000198626" ~ "RYR2",
  Gene == "ENSG00000143632" ~ "ACTA1",
  Gene == "ENSG00000077522" ~ "ACTN2",
  Gene == "ENSG00000072110" ~ "ACTN1",
  Gene == "ENSG00000108819" ~ "PPP1R9B",
  Gene == "ENSG00000151067" ~ "CACNA1C",
  Gene == "ENSG00000066027" ~ "PPP2R5A",
  TRUE ~ NA_character_)) -> expr_wga
na.omit(expr_wga) -> expr_wga

expr_wga %>%pivot_longer(!Gene, 
                         names_to = "Sample", 
                         values_to = "log2.expr") %>% 
  mutate(wGA = case_when(Sample == "ERR2704715" ~ "16",
                         Sample == "ERR2704716" ~ "16",
                         Sample == "ERR2704714" ~ "14", 
                         Sample == "ERR2704720" ~ "9",
                         Sample == "ERR2704719" ~ "9",
                         Sample == "ERR2704713" ~ "12",
                         Sample == "ERR2704717" ~ "16",
                         Sample == "ERR2704718" ~ "9",
                         Sample == "ERR2704712" ~ "12")) -> expr_wga_long

expr_wga_long %>%  filter(wGA != "14") -> expr_wga_long

expr_wga_long$wGA <- factor(expr_wga_long$wGA, levels = c("9", "12", "16"))

expr_wga_long %>%  ggplot(aes(x = wGA, y = log2.expr)) +
  geom_boxplot() + 
  facet_wrap(~Gene) + 
  theme_bw() +
  xlab("Gestational Age (weeks)") +
  ylab("log2(expression)") +
  theme(text = element_text(size = 20), strip.text = element_text(face = "italic")) +
  scale_y_continuous(breaks = c(8,10,12,14), limits = c(7,15))-> human_cardiac_expr
#  labs(caption = "rlog normalized with DESEq2. Raw data source: E-MTAB-7031")

ggsave("human_cardiac_expr.png",
       human_cardiac_expr,
       height = 8,
       width = 6,
       units = "in",
       dpi = 300)

mouse_expr <- read.table("mouse_devl_expr/subset_df.csv", header = T, sep = ",")

mouse_expr <- mouse_expr[,-c(2,3)]

mouse_expr %>% pivot_longer(!Gene, names_to = "age", values_to = "log2.expr") -> mouse_expr_log

mouse_expr_log %>% mutate(age = case_when(
  age == "P1.1" ~ "P1",
  age == "P1.2" ~ "P1",
  age == "P1" ~ "P1",
  age == "P23" ~ "P23",
  age == "P23.1" ~ "P23",
  age == "P23.2" ~ "P23",
  age == "P4" ~ "P4",
  age == "P4.1" ~ "P4",
  age == "P4.2" ~ "P4",
  age == "P9" ~ "P9",
  age == "P9.1" ~ "P9",
  age == "P9.2" ~ "P9"
)) -> mouse_expr_log

mouse_expr_log$age <- factor(mouse_expr_log$age, levels = c("P1", "P4", "P9", "P23"))

mouse_expr_log %>%  ggplot(aes(x = age, y = log2(log2.expr))) +
  geom_boxplot() + 
  facet_wrap(~Gene) + 
  theme_bw() +
  xlab("Postnatal days (P)") +
  ylab("log2(expression)") +
  theme(text = element_text(size = 20), strip.text = element_text(face = "italic")) +
  scale_y_continuous(breaks = c(8,10,12,14), limits = c(7,15))-> mouse_expr_plot

ggsave("mouse_expr_plot.png",
       mouse_expr_plot,
       height = 8,
       width = 6,
       units = "in",
       dpi = 300)

