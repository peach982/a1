library(reshape2)
library(ggplot2)

# 读取数据：第一列是样本名（行名）
df <- read.csv(
  "/data/liuqinglin/ST/spatial/homo_result_heatmap/homo_results.csv",
  row.names = 1, check.names = FALSE, stringsAsFactors = FALSE
)

mat <- as.matrix(df)

# 固定列顺序（cell type）
order_cell <- c(
  "Bcell","EEC","Enterocyte","Goblet","Myeloid",
  "Paneth","Plasma","Stem","TA","T_NKcell"
)
mat <- mat[, order_cell, drop = FALSE]

# ============== 转置矩阵 ==============
mat_transpose <- t(mat)
# =====================================

# 转成长表（转置后：行名→celltype，列名→sample）
data <- melt(mat_transpose)
colnames(data) <- c("celltype", "sample", "value")

# 固定显示顺序
data$celltype <- factor(data$celltype, levels = rownames(mat_transpose))
data$sample   <- factor(data$sample,   levels = colnames(mat_transpose))

# 自定义颜色
custom_colors <- c(
  "#21529b",  # 深蓝
  "#7fa6d8",  
  "#ffffff",  # 白
  "#e6a19a",  
  "#c13933"   # 深红
)
# 固定色阶范围
fill_limits <- c(-1, 2)

# 绘图：转置后 sample 在 x 轴，celltype 在 y 轴（横向）
p_transpose <- ggplot(data, aes(x = sample, y = celltype, fill = value)) +
  geom_tile(color = "white", alpha = 0.85) +
  scale_fill_gradientn(
    colours  = custom_colors,
    limits   = fill_limits,
    oob      = scales::squish,
    na.value = "grey90",
    name     = "Score"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 12, angle = 90, hjust = 1),
    axis.text.y  = element_text(size = 12),
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    panel.background = element_blank(),
    legend.position = "right"
  )

# 保存：宽度 > 高度 ⇒ 视觉横向
ggsave("homoresult_heatmap_transposed.pdf",
       p_transpose, width = 9, height = 8)