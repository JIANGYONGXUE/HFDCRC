install.packages("Seurat")
setwd('~/Desktop/002/2025/NC /DATA/scRNAseq')
DimPlot(Tcells_pbmc1, reduction = 'tsne', split.by = 'celltype')
Tcells_pbmc1
## 1. 重新算邻接图（使用现有 2000 个高变基因）
DefaultAssay(Tcells_pbmc1) <- "RNA"
Tcells_pbmc1 <- FindNeighbors(Tcells_pbmc1, dims = 1:30, verbose = FALSE)
## 2. 聚类——分辨率 0.8 起步，想更细就往上加
Tcells_pbmc1 <- FindClusters(Tcells_pbmc1, resolution = 0.8, verbose = FALSE)
## 3. 重新跑 UMAP（基于新的聚类结果）
Tcells_pbmc1 <- RunUMAP(Tcells_pbmc1, dims = 1:30, verbose = FALSE)
Tcells_pbmc1 <- RunTSNE(Tcells_pbmc1, dims = 1:30)
## 4. 可视化
DimPlot(Tcells_pbmc1, label = TRUE, pt.size = 0.8, split.by = 'orig.ident')
## 5. 快速看每个 cluster 多少细胞
table(Idents(Tcells_pbmc1))

FeaturePlot(Tcells_pbmc1, reduction = 'umap', features = c('Cd8a','Cd4'))
#看细胞占比
library(dplyr)
prop.table(
  table(Idents(Tcells_pbmc1), Tcells_pbmc1$orig.ident),
  margin = 2                      # 2 = 按列分组求占比，即每组内总和 100%
) %>%
  as.data.frame() %>%
  rename(Cluster = Var1, Group = Var2, Proportion = Freq)
#对叠图
library(dplyr)
library(ggplot2)

# 1. 提取 meta.data 并转成 tibble
plot_dat <- Tcells_pbmc1@meta.data %>%               # 或者 as.data.frame(Tcells_pbmc1[[]])
  as_tibble() %>%                                    # 确保是 tibble，dplyr 才能用
  group_by(orig.ident, cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(pct = n / sum(n))

# 2. 堆叠柱状图
ggplot(plot_dat,
       aes(x = orig.ident, y = pct, fill = factor(cluster))) +
  geom_col(width = 0.7, color = "black", size = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Group",
       y = "Percentage within group",
       fill = "Cluster") +
  theme_minimal(base_size = 14)

## 1. 计算所有 cluster 的 marker 基因（Wilcoxon，默认只报正向）
markers_all <- FindAllMarkers(
  Tcells_pbmc1,
  only.pos = TRUE,      # 只要上调基因
  min.pct = 0.25,       # 在 cluster 里表达比例 ≥ 25%
  logfc.threshold = 0.25,
  verbose = FALSE
)

## 2. 每个 cluster 只保留 Top10（按 avg_log2FC 排序）
top10_per_cluster <- markers_all %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC) %>%   # 取 FC 最大的 10 条
  arrange(cluster, -avg_log2FC)

#cd4_core   <- c("Cd3d","Cd4")                       # 锁定
#cd4_exclude<- c("Cd8a","Cd19","Ms4a1","Cd14","Ncr1") # 排除
#naive      <- c("Lef1","Tcf7","Sell","Ccr7","Il7r")
#th1        <- c("Tbx21","Ifng","Cxcr3","Gzma")
#th17       <- c("Rorc","Il17a","Il17f","Ccr6","Il23r")
#treg       <- c("Foxp3","Il2ra","Ctla4","Ikzf2")
#tfh        <- c("Bcl6","Cxcr5","Pdcd1","Cxcl13")
#trm        <- c("Cd69","Itgae","Itga1")
#tem        <- c("Ccr5","Gzmk","S1pr1","Cd44")
#tscm       <- c("Slamf6","Tcf7","Cd27"
#all_markers <- c(cd4_core, naive, th1, th17, treg, tfh, trm, tem, tscm)

FeaturePlot(Tcells_pbmc1, features = c('Cd4','Cd8a',
                                       'Lef1','Tcf7','Sell','Ccr7','Il7r',
                                       "Tbx21","Ifng","Cxcr3","Gzma",
                                       "Rorc","Il17a","Il17f","Ccr6","Il23r",
                                       "Foxp3","Il2ra","Ctla4","Ikzf2",
                                       "Bcl6","Cxcr5","Pdcd1","Cxcl13",
                                       "Cd69","Itgae","Itga1",
                                       "Ccr5","Gzmk","S1pr1","Cd44",
                                       "Slamf6","Tcf7","Cd27"))

FeaturePlot(Tcells_pbmc1, features = c("Slamf6","Tcf7","Cd27"
                                       ))
FeaturePlot(Tcells_pbmc1, features = c('Cd4','Cd8a'))

annot2 <- c("CD4_Th1", "CD4_Tm", "CD8_Trm", "CD8_Tem", "Treg",
           "Th17", "CD4_Tm", "other cell", "Cycling_T", "other cell")
DoHeatmap(Tcells_pbmc1 , features = top10_per_cluster$gene, size = 3)
names(annot2)<-levels(Tcells_pbmc1)
Tcells_pbmc1 <- RenameIdents(Tcells_pbmc1, annot2)
DimPlot(Tcells_pbmc1, reduction = 'umap',
        label = TRUE, pt.size = 1) + NoLegend()
Tcells_pbmc1$celltype2 <- Idents(Tcells_pbmc1)
library(scales)
# 1. 按 sample + celltype1 统计
plot_dat <- Tcells_pbmc1@meta.data %>%
  count(orig.ident, celltype2) %>%        # 先数个数
  group_by(orig.ident) %>%                # 每个样本内部分比
  mutate(pct = n / sum(n))

# 2. 画图
ggplot(plot_dat,
       aes(x = orig.ident, y = pct, fill = celltype2)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.15) +
  scale_y_continuous(labels = label_percent(scale = 100)) +
  labs(x = "Sample", y = "Percentage within sample", fill = "Cell type") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
####去除掉other cell再画对叠图

library(dplyr)
library(ggplot2)
library(scales)

# 1. 过滤掉 other cell
plot_dat <- Tcells_pbmc1@meta.data %>%
  filter(celltype2 != "other cell") %>%   # 关键一步
  count(orig.ident, celltype2) %>%
  group_by(orig.ident) %>%
  mutate(pct = n / sum(n))

# 2. 重新堆叠
ggplot(plot_dat,
       aes(x = orig.ident, y = pct, fill = celltype2)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.15) +
  scale_y_continuous(labels = label_percent(scale = 100)) +
  labs(x = "Sample", y = "Percentage within sample", fill = "Cell type") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


table(Tcells_pbmc1$celltype2)
Tcells_pbmc1
##再去掉CD8
plot_dat1 <- Tcells_pbmc1@meta.data %>%
  filter(!celltype2 %in% c('CD8_Trm','CD8_Tem','other cell','Cycling_T')) %>%   # ✅ 正确排除
  count(orig.ident, celltype2) %>%
  group_by(orig.ident) %>%
  mutate(pct = n / sum(n))

##只统计CD4的几个亚群
ggplot(plot_dat1,
       aes(x = orig.ident, y = pct, fill = celltype2)) +
  geom_col(width = 0.9, color = "black", linewidth = 0.15) +
  scale_y_continuous(labels = label_percent(scale = 100)) +
  labs(x = "Sample", y = "Percentage within sample", fill = "Cell type") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    panel.grid    = element_blank(),          # 去掉网格线
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1)  # 加外边框
  )

##重新统计Tcell的TSNE图，分为CD4和CD8，CD4的再细分定义为Th1、Tm、Trg、Th17几种。

Idents(Tcells_pbmc1) <- 'celltype2'
# 1. 留下非目标细胞（=! 反向选择）
keep_cells <- WhichCells(Tcells_pbmc1,
                         expression = !celltype2 %in% c('CD8_Trm','CD8_Tem','other cell','Cycling_T'))
# 2. 生成干净对象
CD4Tcells_clean <- subset(Tcells_pbmc1, cells = keep_cells)
# 3. （可选）重新算 t-SNE——用原高变基因即可
#    如果原对象已做过 PCA/Harmony，可直接用原坐标；这里示范重跑。
CD4Tcells_clean <- RunTSNE(CD4Tcells_clean,
                        dims = 1:30,        # 与你之前一致
                        reduction = 'pca',  # 或 'harmony' 如果你有
                        check.seed = 42)    # 复现用
# 4. 画图
DimPlot(CD4Tcells_clean,
        reduction = 'tsne',
        group.by = 'celltype2',
        label = FALSE,
        pt.size = 1.5)

keep_cells1 <- WhichCells(Tcells_pbmc1,
                         expression = !celltype2 %in% c('other cell','Cycling_T'))
# 2. 生成干净对象
Tcells_clean <- subset(Tcells_pbmc1, cells = keep_cells1)
Tcells_clean <- RunTSNE(Tcells_clean,
                        dims = 1:30,        # 与你之前一致
                        reduction = 'pca',  # 或 'harmony' 如果你有
                        check.seed = 42)    # 复现用
DimPlot(Tcells_clean,
        reduction = 'tsne',
        group.by = 'celltype2',
        label = TRUE,
        pt.size = 1.5)

library(ggforce)   # 先安装一次 install.packages("ggforce")
library(ggplot2)

# 1. 提取 t-SNE 坐标和分组信息
tsne_dat <- FetchData(Tcells_clean, vars = c("tSNE_1", "tSNE_2", "celltype2"))

# 2. 手动画散点 + 虚线椭圆
ggplot(tsne_dat,
       aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(aes(color = celltype2), size = 1.5) +
  geom_mark_ellipse(aes(color = celltype2, group = celltype2),
                    alpha = 0,          # 无填充
                    linewidth = 0.5,    # 线粗
                    linetype = 2) +     # 虚线
  labs(x = "t-SNE 1", y = "t-SNE 2", color = "") +
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank())

install.packages("concaveman")   # 第一次装一次
library(concaveman)
library(ggplot2)
dat <- FetchData(Tcells_clean, vars = c("tSNE_1", "tSNE_2", "celltype2"))
# 1. 重新命名列（与全局 aes 一致）
hull_ls <- lapply(unique(dat$celltype2), function(ct){
  df <- dat[dat$celltype2 == ct, ]
  coords <- as.matrix(df[, c("tSNE_1", "tSNE_2")])
  hull <- concaveman(coords, concavity = 1.5)
  hull <- as.data.frame(hull)
  colnames(hull) <- c("tSNE_1", "tSNE_2")   # 关键！
  hull$celltype2 <- ct
  hull
})
hull_dat <- do.call(rbind, hull_ls)

# 2. 画图（无需再改 aes）
ggplot(dat, aes(tSNE_1, tSNE_2)) +
  geom_point(aes(color = celltype2), size = 1.5, show.legend = TRUE) +
  geom_polygon(data = hull_dat,
               aes(color = celltype2, group = celltype2),
               fill = NA,
               linewidth = 0.6,
               linetype = 2) +
  labs(x = "t-SNE 1", y = "t-SNE 2", color = "") +
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank())

DimPlot(Tcells_clean,
        reduction = 'tsne',
        group.by = 'celltype2',
        label = FALSE,
        pt.size = 1.5) +
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA))  # 仅留外框

FeaturePlot(Tcells_clean,features = c('Cd4'),reduction = 'tsne')


