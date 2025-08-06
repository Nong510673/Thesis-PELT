# ============================================================
# 0.  (可選) 清理環境
# ============================================================
# rm(list = ls())

# ============================================================
# 1. 套件
# ============================================================
library(phyloseq)
library(ape)
library(vegan)
library(ggplot2)
library(readxl)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(scales)
library(grid)        # for unit()
set.seed(123)

# ============================================================
# 2. 讀取 OTU / metadata / 樹
# ============================================================
otu        <- read.csv("E:/2025_UASB_MOTHER/KTU/KTU.table.summary.csv",
                       row.names = 1, check.names = FALSE)
sample_df  <- read_excel("E:/2025_UASB_MOTHER/KTU/sample_list.xlsx")
tree       <- read.tree("E:/2025_UASB_MOTHER/pruned_tree.nwk")

# ============================================================
# 3. 保留 Anl_count == "on" 的樣本
# ============================================================
sample_df <- sample_df %>%
  filter(Anl_count == "on") %>%
  mutate(Day = as.numeric(Day))

# ============================================================
# 4. 組裝 phyloseq 物件
# ============================================================
otu_mat        <- as.matrix(otu[, sample_df$Sample_ID, drop = FALSE])
otu_table_ps   <- otu_table(otu_mat, taxa_are_rows = TRUE)
sample_meta    <- sample_df %>% select(-Anl_count) %>% as.data.frame()
rownames(sample_meta) <- sample_meta$Sample_ID
sample_meta$Sample_ID <- NULL
sample_data_ps <- sample_data(sample_meta)
physeq         <- phyloseq(otu_table_ps, sample_data_ps, phy_tree(tree))

# ============================================================
# 5. Weighted UniFrac 距離
# ============================================================
wuni_dist <- UniFrac(physeq, weighted = TRUE, normalized = TRUE)

# ============================================================
# 6. PCoA (前 3 軸)
# ============================================================
pcoa3_res <- cmdscale(wuni_dist, eig = TRUE, k = 3)
pcoa3_df  <- as.data.frame(pcoa3_res$points)
colnames(pcoa3_df) <- c("PC1", "PC2", "PC3")
pcoa3_df$Sample_ID <- rownames(pcoa3_df)

eig_pct <- round(100 * pcoa3_res$eig / sum(pcoa3_res$eig), 2)
message(sprintf("PC1–PC3 variance explained (%%): %.2f, %.2f, %.2f",
                eig_pct[1], eig_pct[2], eig_pct[3]))

# ============================================================
# 7. 合併 metadata ＋ 分段
# ============================================================
pcoa3_df <- pcoa3_df %>%
  left_join(sample_df, by = "Sample_ID") %>%
  arrange(Day) %>%
  mutate(
    seg = cut(
      Day,
      breaks = c(-Inf, 200, 345, 597, 870, Inf),
      labels = c("0–200d", "200–345d", "345–597d",
                 "597–870d", ">870d")
    )
  )

# -- 關鍵日期 (三角形顯示) ----------------------
special_days <- c(345, 870)
special_df <- pcoa3_df %>% filter(Day %in% special_days)
normal_df  <- pcoa3_df %>% filter(!Day %in% special_days)

# ============================================================
# 8. 互動 3D PCoA (Plotly) → 匯出 HTML
# ============================================================
fig <- plot_ly() %>%
  ## 普通圓點
  add_trace(
    data = normal_df,
    x = ~PC1, y = ~PC2, z = ~PC3,
    type = "scatter3d", mode = "markers+text",
    textposition = "top center",
    color = ~seg,
    colors = c("grey60", "#00A600", "#005AB5", "purple", "red"),
    marker = list(size = 5, symbol = "circle",
                  line = list(width = 1, color = "black")),
    text = ~Day,
    textfont = list(size = 10, family = "Arial"),
    hoverinfo = "text",
    hovertext = ~paste0(
      "Sample: ", Sample_ID, "<br>",
      "Day: ", Day, "<br>",
      "PC1: ", round(PC1, 2), "<br>",
      "PC2: ", round(PC2, 2), "<br>",
      "PC3: ", round(PC3, 2)
    ),
    name = "Normal"
  ) %>%
  ## 關鍵三角點
  add_trace(
    data = special_df,
    x = ~PC1, y = ~PC2, z = ~PC3,
    type = "scatter3d", mode = "markers+text",
    textposition = "top center",
    color = ~seg,
    colors = c("grey60", "#00A600", "#005AB5", "purple", "red"),
    marker = list(size = 14, symbol = "triangle-up",
                  line = list(width = 1.5, color = "black")),
    text = ~Day,
    textfont = list(size = 12, family = "Arial"),
    hoverinfo = "text",
    hovertext = ~paste0(
      "Sample: ", Sample_ID, "<br>",
      "Day: ", Day, "<br>",
      "PC1: ", round(PC1, 2), "<br>",
      "PC2: ", round(PC2, 2), "<br>",
      "PC3: ", round(PC3, 2)
    ),
    name = "Key day"
  ) %>%
  layout(
    title = list(text = "3D PCoA (PC1–PC3)",
                 font = list(size = 22, family = "Arial Bold")),
    scene = list(
      xaxis = list(
        title = list(text = paste0("<b>PC1 (", eig_pct[1], "%)</b>")),
        showgrid = FALSE, zeroline = FALSE,
        ticks = "outside", ticklen = 5, tickwidth = 2,
        linecolor = "black", linewidth = 2,
        tickfont = list(size = 15, family = "Arial Bold")
      ),
      yaxis = list(
        title = list(text = paste0("<b>PC2 (", eig_pct[2], "%)</b>")),
        showgrid = FALSE, zeroline = FALSE,
        ticks = "outside", ticklen = 5, tickwidth = 2,
        linecolor = "black", linewidth = 2,
        tickfont = list(size = 15, family = "Arial Bold")
      ),
      zaxis = list(
        title = list(text = paste0("<b>PC3 (", eig_pct[3], "%)</b>")),
        showgrid = FALSE, zeroline = FALSE,
        ticks = "outside", ticklen = 5, tickwidth = 2,
        linecolor = "black", linewidth = 2,
        tickfont = list(size = 15, family = "Arial Bold")
      ),
      bgcolor = "rgba(0,0,0,0)"
    ),
    legend = list(title = list(text = "Time segment"))
  )

saveWidget(fig,
           file = "C:/Users/User/Desktop/PCoA3D_interactive.html",
           selfcontained = TRUE)

# ============================================================
# 9. 灰階背景層 for ggplot
# ============================================================
background_layer_2 <- list(
  annotate("rect", xmin =   0, xmax = 120,  fill = "grey80", ymin = -Inf, ymax = Inf, alpha = 0.6),
  annotate("rect", xmin = 120, xmax = 200,  fill = "grey90", ymin = -Inf, ymax = Inf, alpha = 0.6),
  annotate("rect", xmin = 200, xmax = 345,  fill = "grey70", ymin = -Inf, ymax = Inf, alpha = 0.6),
  annotate("rect", xmin = 345, xmax = 395,  fill = "grey80", ymin = -Inf, ymax = Inf, alpha = 0.6),
  annotate("rect", xmin = 395, xmax = 596,  fill = "grey90", ymin = -Inf, ymax = Inf, alpha = 0.6),
  annotate("rect", xmin = 597, xmax = 870,  fill = "grey70", ymin = -Inf, ymax = Inf, alpha = 0.6),
  annotate("rect", xmin = 870, xmax = 1000, fill = "grey80", ymin = -Inf, ymax = Inf, alpha = 0.6)
)

# ============================================================
# 10. PC 軸 vs Day 散布圖 (只顯示，不匯出)
# ============================================================
base_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid        = element_blank(),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.title.x      = element_text(size = 14, face = "bold"),
    axis.title.y      = element_text(size = 14, face = "bold"),
    axis.text         = element_text(size = 12, face = "bold"),  # 數值粗體
    axis.ticks        = element_line(linewidth = 0.8),           # 刻度線粗
    axis.ticks.length = unit(5, "pt"),                           # 刻度長
    legend.position   = "none"
  )

axis_scales <- scale_x_continuous(
  name   = "Operation Time (Day)",
  limits = c(0, 1000),
  breaks = seq(0, 1000, 100),
  expand = expansion(mult = c(0, 0))
)

# ---- PC1 ----
p_pc1 <- ggplot(pcoa3_df, aes(Day, PC1)) +
  background_layer_2 +
  geom_point(size = 2) +
  axis_scales +
  scale_y_continuous(name = "PCoA 1 Coordinate") +
  base_theme
print(p_pc1)

# ---- PC2 ----
p_pc2 <- ggplot(pcoa3_df, aes(Day, PC2)) +
  background_layer_2 +
  geom_point(size = 2) +
  axis_scales +
  scale_y_continuous(name = "PCoA 2 Coordinate") +
  base_theme
print(p_pc2)

# ---- PC3 ----
p_pc3 <- ggplot(pcoa3_df, aes(Day, PC3)) +
  background_layer_2 +
  geom_point(size = 2) +
  axis_scales +
  scale_y_continuous(name = "PCoA 3 Coordinate") +
  base_theme
print(p_pc3)


# ============================================================
# 11. 匯出 PCoA 軸對應數值變化 (CSV)
# ============================================================
write.csv(
  pcoa3_df %>%
    select(Sample_ID, Day, PC1, PC2, PC3, seg),
  file = "C:/Users/User/Desktop/PCoA3_axes_with_day.csv",
  row.names = FALSE
)
