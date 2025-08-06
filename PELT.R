# ==============================================================================
# Change-Point Analysis for COD Removal Efficiency in UASB Reactor
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. 安裝與載入所需套件（如未安裝請先 run install.packages）
# ------------------------------------------------------------------------------
required_pkgs <- c("readxl", "dplyr", "ggplot2", "changepoint", "ecp", "cpm", "strucchange", "segmented")
invisible(lapply(required_pkgs, function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

# ------------------------------------------------------------------------------
# 1. 參數設定
# ------------------------------------------------------------------------------
data_path   <- "C:/Users/User/Desktop/UASB_DATA_Final_final.xlsx"
sheet_name  <- "工作表1"

# 各方法參數
pelt_method    <- "PELT"
pelt_teststat  <- "Normal"
pelt_penalty   <- "MBIC"
pelt_penvalue  <- NULL
pelt_minseglen <- 16

ecp_siglvl <- 0.05
ecp_R      <- 199

cpm_type   <- "Student"
sc_h       <- 0.1
seg_npsi   <- 2

# ------------------------------------------------------------------------------
# 2. 讀取與整理資料
# ------------------------------------------------------------------------------
df <- read_excel(data_path, sheet = sheet_name) %>%
  rename(COD_raw = `COD Removal effeciency (%)`) %>%
  mutate(Day = as.numeric(Day), COD = as.numeric(COD_raw))

if (any(is.na(df$Day))) stop("Day 欄含 NA")
if (any(is.na(df$COD))) warning("COD 欄含 NA 將略過")

# ------------------------------------------------------------------------------
# 3. 各方法變點偵測
# ------------------------------------------------------------------------------
# 3.1 PELT
pelt_args <- list(method = pelt_method, test.stat = pelt_teststat,
                  penalty = pelt_penalty, minseglen = pelt_minseglen)
if (identical(pelt_penalty, "Manual") && !is.null(pelt_penvalue)) {
  pelt_args$pen.value <- pelt_penvalue
}
cp_pelt  <- do.call(cpt.meanvar, c(list(df$COD), pelt_args))
idx_pelt <- cpts(cp_pelt)

# 3.2 e.divisive
ecp_out  <- e.divisive(as.matrix(df$COD), sig.lvl = ecp_siglvl, R = ecp_R)
idx_ecp  <- setdiff(ecp_out$estimates, c(1, length(df$COD)))

# 3.3 cpm
cpm_out  <- detectChangePointBatch(df$COD, cpmType = cpm_type)
idx_cpm  <- cpm_out$changePoint

# 3.4 strucchange
bp       <- breakpoints(COD ~ Day, data = df, h = sc_h)
idx_struc <- na.omit(bp$breakpoints)

# 3.5 segmented
lm0      <- lm(COD ~ Day, data = df)
seg_out  <- segmented(lm0, seg.Z = ~Day, npsi = seg_npsi)
idx_seg  <- round(seg_out$psi[, "Est."])

# ------------------------------------------------------------------------------
# 4. 索引對應時間點
# ------------------------------------------------------------------------------
map_day <- function(idx) df$Day[idx]
day_pelt   <- map_day(idx_pelt)
day_ecp    <- map_day(idx_ecp)
day_cpm    <- map_day(idx_cpm)
day_struc  <- map_day(idx_struc)
day_seg    <- map_day(idx_seg)

# ------------------------------------------------------------------------------
# 5. 匯總結果表
# ------------------------------------------------------------------------------
library(tibble)
seg_df <- bind_rows(
  tibble(method = "PELT",        cp_idx = idx_pelt,  cp_day = day_pelt),
  tibble(method = "e.divisive",  cp_idx = idx_ecp,   cp_day = day_ecp),
  tibble(method = "cpm",         cp_idx = idx_cpm,   cp_day = day_cpm),
  tibble(method = "strucchange", cp_idx = idx_struc, cp_day = day_struc),
  tibble(method = "segmented",   cp_idx = idx_seg,   cp_day = day_seg)
)

# ------------------------------------------------------------------------------
# 5.5 印出變點時間點
# ------------------------------------------------------------------------------
cat("\n===== Detected Change Points (By Method) =====\n")
seg_df %>%
  arrange(method, cp_day) %>%
  group_by(method) %>%
  summarise(change_points = paste(cp_day, collapse = ", ")) %>%
  mutate(output = paste0(method, ": ", change_points)) %>%
  pull(output) %>%
  cat(sep = "\n")

# ------------------------------------------------------------------------------
# 6. 繪圖（各方法變點以虛線表示，X 軸每 100 天一格）
# ------------------------------------------------------------------------------
ggplot(df, aes(x = Day, y = COD)) +
  geom_line(color = "gray60") +
  geom_point(size = 1) +
  geom_vline(
    data = seg_df,
    aes(xintercept = cp_day, color = method),
    linetype = "dashed", size = 0.8
  ) +
  facet_wrap(~ method, ncol = 1, scales = "free_x") +
  scale_x_continuous(breaks = seq(0, max(df$Day, na.rm = TRUE), by = 100)) +
  labs(
    x     = "Operation Day",
    y     = "COD Removal Efficiency (%)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
