###############################################################################
# 稀釋曲線與多樣性指標計算（KTU 表）                               R 4.x
# ---------------------------------------------------------------------------
# 1. 讀取 KTU 計數表 → 稀釋曲線 (Chao1 / Shannon)
# 2. 將所有樣本隨機稀釋至 8 000 reads，計算 Observed、Chao1、Shannon、
#    Simpson 與 Pielou Evenness，並匯出 CSV
###############################################################################

# ─────────────────────────────────────────────────────────────────────────────
# 0) 套件
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
})

# ─────────────────────────────────────────────────────────────────────────────
# 1) 讀取 KTU 計數表
# ─────────────────────────────────────────────────────────────────────────────
file_path <- "E:/2025_Batch/KTU/KTU.table.summary.csv"   # ★自行修改
otu_table <- read.csv(file_path, row.names = 1, check.names = FALSE)

# 將 read 數四捨五入 → 整數
otu_table <- round(otu_table)

# ─────────────────────────────────────────────────────────────────────────────
# 2) 稀釋曲線 (Chao1、Shannon) —— 抽樣深度 0–15 000，每 1 000 一點
# ─────────────────────────────────────────────────────────────────────────────
depths    <- seq(0, 15000, by = 1000)
num_iter  <- 10   # 每深度重複抽樣次數

compute_rarefaction <- function(counts, depths, num_iter = 10) {
  chao1 <- shannon <- numeric(length(depths))
  tot   <- sum(counts)

  for (i in seq_along(depths)) {
    d <- depths[i]
    if (d == 0) {
      chao1[i]   <- 0
      shannon[i] <- 0
    } else if (d >= tot) {
      chao1[i]   <- estimateR(counts)[2]
      shannon[i] <- diversity(counts, index = "shannon")
    } else {
      chao1_iter <- shannon_iter <- numeric(num_iter)
      for (j in seq_len(num_iter)) {
        subs <- rrarefy(counts, sample = d)
        chao1_iter[j]   <- estimateR(subs)[2]
        shannon_iter[j] <- diversity(subs, index = "shannon")
      }
      chao1[i]   <- mean(chao1_iter)
      shannon[i] <- mean(shannon_iter)
    }
  }
  list(chao1 = chao1, shannon = shannon)
}

# 計算並整理成長格式
raref_list <- lapply(colnames(otu_table), function(smp) {
  out <- compute_rarefaction(otu_table[[smp]], depths, num_iter)
  rbind(
    data.frame(Depth = depths, Sample = smp,
               Index = "Chao1",  Value = out$chao1),
    data.frame(Depth = depths, Sample = smp,
               Index = "Shannon", Value = out$shannon)
  )
})
raref_df <- bind_rows(raref_list)

# ─────────────────────────────────────────────────────────────────────────────
# 3) 稀釋曲線繪圖
# ─────────────────────────────────────────────────────────────────────────────
plot_theme <- theme(
  panel.background = element_blank(),
  panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1.2),
  axis.title       = element_text(size = 14, face = "bold"),
  axis.text        = element_text(size = 12, face = "bold"),
  axis.ticks       = element_line(colour = "black", linewidth = 1.2),
  legend.position  = "none"
)

# Chao1
p_chao1 <- ggplot(filter(raref_df, Index == "Chao1"),
                  aes(Depth, Value, colour = Sample)) +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(limits = c(0, 15000), breaks = seq(0, 15000, 1000)) +
  scale_y_continuous(limits = c(0, 550),   breaks = seq(0, 550,    50))  +
  labs(x = "Depth", y = "Chao1") + plot_theme
print(p_chao1)

# Shannon
p_shannon <- ggplot(filter(raref_df, Index == "Shannon"),
                    aes(Depth, Value, colour = Sample)) +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(limits = c(0, 15000), breaks = seq(0, 15000, 1000)) +
  scale_y_continuous(limits = c(0, 6),     breaks = seq(0, 6,      1))  +
  labs(x = "Depth", y = "Shannon") + plot_theme
print(p_shannon)

# ─────────────────────────────────────────────────────────────────────────────
# 4) 稀釋至 8 000 reads → 多樣性指標
# ─────────────────────────────────────────────────────────────────────────────
target_depth <- 8000
sample_reads <- colSums(otu_table)
keep_idx     <- sample_reads >= target_depth

if (any(!keep_idx)) {
  warning("下列樣本 reads < ", target_depth, "，已略過： ",
          paste(colnames(otu_table)[!keep_idx], collapse = ", "))
}

# 稀釋（樣本為列 → 需轉置）
rare_mat <- t(rrarefy(t(otu_table[, keep_idx, drop = FALSE]),
                      sample = target_depth))

div_df <- data.frame(
  Sample    = colnames(rare_mat),
  Observed  = colSums(rare_mat > 0),
  Chao1     = apply(rare_mat, 2, function(x) estimateR(x)[2]),
  Shannon   = apply(rare_mat, 2, function(x) diversity(x, "shannon")),
  Simpson   = apply(rare_mat, 2, function(x) diversity(x, "simpson"))
)
div_df$Evenness <- with(div_df, Shannon / log(Observed))

# 匯出
out_csv <- "E:/2025_Batch/KTU/diversity_rarefied_8000.csv"  # ★自行修改
write.csv(div_df, out_csv, row.names = FALSE)
message("✓ 已匯出多樣性指標： ", out_csv)
