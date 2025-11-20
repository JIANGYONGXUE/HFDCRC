setwd("~/Desktop/002/2025/NC /DATA/00metaanalysis")
group2 <- read_excel("20weeks/group2.xlsx")
test2 <- read_excel("20weeks/test2.xlsx")
library(RColorBrewer)
#install.packages('tidyverse')
library(tidyverse)
library(ggplot2)
library(patchwork)
# 0. 统一列名 -------------------------------------------------------------
df_abu3 <- eggerspecies %>% 
  pivot_longer(-1, names_to = "Sample", values_to = "abundance") %>% 
  rename(taxon = 1) %>% 
  left_join(group1, by = "Sample")

# 1. 分组颜色 -------------------------------------------------------------
group_cols <- c("ND"        = "#E41A1C",
                "HFD"         = "#377EB8", 
                "ND-AOMDSS" = "#FF7F00",
                "HFD-AOMDSS"  = "#4DAF4A")
group_cols1 <- c("ND"        = "#E41A1C",
                "HFD"         = "#377EB8")
# 2. 顶部颜色条 -----------------------------------------------------------
p_bar3 <- ggplot(df_abu3 %>% distinct(Sample, Group),
                aes(x = reorder(Sample, as.numeric(factor(Sample))),
                    y = 1, fill = Group)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_manual(values = group_cols1) +
  theme_void() +
  theme(legend.position = "none")
target_order <- names(df_abu3)[-1] %>%           # 原列顺序
  {.[order(!str_detect(., "^HFD"))]} 
# 3. 主热图（不挑菌，完整3行） -------------------------------------------
p_heat3 <- ggplot(df_abu3,
                 aes(x = reorder(Sample, as.numeric(factor(Sample))),
                     y = taxon,               # 就3个菌，全显示
                     fill = abundance)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_gradientn(colours = c("#2166AC", "white", "#B2182B"),
                       name = "abundance") +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),   # 显示菌名
        panel.grid  = element_blank())

# 4. 拼图 -----------------------------------------------------------------

p_bar3 / p_heat3 +
  plot_layout(heights = c(1, 4))

# 5. 保存（可选） ---------------------------------------------------------
ggsave("10weeks_elenta_taxa_heatmap.pdf", width = 6, height = 4.5, dpi = 300)

# 0. 生成目标顺序（ND 在前，HFD 在后）
target_order <- names(eggerspecies)[-1] %>%  
  {.[order(!str_detect(., "^HFD"))]}

# 1. 整理数据：统一列名小写，锁定顺序
df_abu3 <- eggerspecies %>% 
  pivot_longer(-1, names_to = "Sample", values_to = "abundance") %>% 
  rename(taxon = 1) %>% 
  left_join(group2, by = "Sample") %>% 
  mutate(sample = factor(sample, levels = target_order))
