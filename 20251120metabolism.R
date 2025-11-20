#-----------HFDvsND代谢组KEGG-------
#1.读入文件
library(readr)
kegg <- read_csv("enrichment_kegg_HFD_qian_vs_ND_qian_Up.csv")
#2. 只留 Metabolism，按 p 值升序，取前 10
library(tidyverse)
library(forcats)   # for fct_reorder

top10 <- kegg %>% 
  filter(Classification_level1 == "Metabolism") %>% 
  arrange(p_value) %>% 
  slice_head(n = 10) %>% 
  mutate(Term = str_wrap(Term, width = 38),
         logP = -log10(p_value))

#3.画图：用 fct_reorder 把 Term 按 logP 降序排好
library(showtext)
#install.packages('showtext')
p<- ggplot(top10, aes(x = fct_reorder(Term, logP, .desc = F), y = logP)) +
  geom_col(fill = "#FF0000", width = 0.65) +
  coord_flip() +
  labs(title = "HFD vs ND",
       x = NULL,
       y = expression(-log[10](italic(P)))) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "Arial", face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(colour = "black", size = 1),
    axis.line = element_line(colour = "black", size = 0.5),
    axis.ticks = element_line(colour = "black", size = 0.5),
    axis.ticks.length = unit(4, "pt"),     # 短横线
    axis.text.x = element_text(family = "Arial", size = 10),
    axis.text.y = element_text(family = "Arial", size = 10, hjust = 0)
  )
ggsave("Metabolism_Top10_bar_Arial.png", plot = p,
       width = 6, height = 6, dpi = 600, bg = "white")
