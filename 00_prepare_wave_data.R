# Prepare the Wave 2-priority data set used in the sensitivity analysis.
# For participants examined in both waves, the Wave 2 phenotype is retained;
# Wave 1 is used when no Wave 2 record is available.

df <- read.csv("df_clean.csv", fileEncoding = "UTF-8", stringsAsFactors = FALSE)
cat("=== 原始数据 ===\n")
cat("总行数：", nrow(df), "\n")

required_cols <- c("id_child", "id_parent", "survey_child", "survey_parent")
if (!all(required_cols %in% colnames(df))) {
  stop("缺少必要列：", paste(setdiff(required_cols, colnames(df)), collapse=", "))
}

# ====== 阶段1：剔除跨波次错配 ======
cat("\n=== 阶段1：按波次过滤 ===\n")
df_wave_matched <- df[df$survey_child == df$survey_parent, ]
n_removed <- nrow(df) - nrow(df_wave_matched)
cat(sprintf("过滤前：%d，过滤后：%d，剔除跨波次：%d (%.1f%%)\n",
            nrow(df), nrow(df_wave_matched), n_removed,
            n_removed / nrow(df) * 100))

df_mismatch <- df[df$survey_child != df$survey_parent, ]
write.csv(df_mismatch, "result_wave_mismatch_removed.csv", row.names = FALSE)

# ====== 阶段1.5：诊断重复数量 ======
pair_key_wm    <- paste(df_wave_matched$id_child, df_wave_matched$id_parent, sep="_")
pair_counts_wm <- table(pair_key_wm)
cat("\n=== 阶段1后：重复分布 ===\n")
print(table(pair_counts_wm))

dup_keys_wm <- names(pair_counts_wm[pair_counts_wm > 1])
cat("仍重复的关系数：", length(dup_keys_wm), "\n")

if (max(pair_counts_wm) > 2) {
  cat("警告：存在出现次数>2的关系，请检查！\n")
}

# ====== 阶段2：以Wave2为准合并 ======
cat("\n=== 阶段2：同波次重复合并（Wave2优先策略）===\n")
cat("策略：有Wave2用Wave2；只有Wave1用Wave1；逆转（Wave1阳→Wave2阴）判为阴性\n")

tori_cols  <- grep("^tori", colnames(df), value = TRUE)
age_cols   <- grep("^age",  colnames(df), value = TRUE)
other_cols <- setdiff(colnames(df_wave_matched),
                      c(tori_cols, age_cols, "survey_child", "survey_parent"))

cat("骨隆起合并列：", paste(tori_cols, collapse=", "), "\n")
cat("年龄合并列：",   paste(age_cols,  collapse=", "), "\n")

# 确认Wave1/Wave2的实际字符串
cat("\nsurvey_child实际取值（确认Wave标签）：\n")
print(table(df_wave_matched$survey_child, useNA="always"))
# ⚠️ 根据上面输出确认，默认如下：
WAVE1_LABEL <- "1次"
WAVE2_LABEL <- "2次"
POS_LABEL   <- "1. あり"
NEG_LABEL   <- "0. なし"

# Wave2优先合并函数
merge_tori_wave2 <- function(vals, surveys,
                              wave2_label = WAVE2_LABEL,
                              pos_label   = POS_LABEL,
                              neg_label   = NEG_LABEL) {
  if (all(is.na(vals))) return(NA_character_)

  w2_idx <- which(surveys == wave2_label)

  # 有Wave2：直接用Wave2的值
  if (length(w2_idx) > 0) {
    w2_val <- vals[w2_idx[1]]
    if (!is.na(w2_val)) return(w2_val)
  }

  # 没有Wave2（只有Wave1）：用Wave1的值
  w1_val <- vals[!is.na(vals)][1]
  return(w1_val)
}

safe_max_numeric <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

df_unique <- df_wave_matched[!(pair_key_wm %in% dup_keys_wm), ]
df_dup    <- df_wave_matched[pair_key_wm %in% dup_keys_wm, ]
key_dup   <- paste(df_dup$id_child, df_dup$id_parent, sep="_")

# 统计各类情况
n_w2_pos    <- 0  # Wave2阳（含新发）
n_w2_neg    <- 0  # Wave2阴
n_reversal  <- 0  # 逆转（Wave1阳→Wave2阴）
n_incident  <- 0  # 新发（Wave1阴→Wave2阳）
n_w1_only   <- 0  # 只有Wave1

merged_list <- vector("list", length(unique(key_dup)))
names(merged_list) <- unique(key_dup)

for (k in unique(key_dup)) {
  rows       <- df_dup[key_dup == k, ]
  merged_row <- rows[1, ]

  # 骨隆起列：Wave2优先
  for (col in tori_cols) {
    vals    <- rows[[col]]
    surveys <- rows$survey_child

    w1_idx <- which(surveys == WAVE1_LABEL)
    w2_idx <- which(surveys == WAVE2_LABEL)
    w1_val <- if (length(w1_idx) > 0) vals[w1_idx[1]] else NA
    w2_val <- if (length(w2_idx) > 0) vals[w2_idx[1]] else NA

    # 统计（只在tori_child列统计，避免重复计数）
    if (col == tori_cols[1]) {
      if (!is.na(w2_val)) {
        if (w2_val == POS_LABEL) n_w2_pos <- n_w2_pos + 1
        else n_w2_neg <- n_w2_neg + 1
        if (!is.na(w1_val)) {
          if (w1_val == POS_LABEL && w2_val != POS_LABEL)
            n_reversal <- n_reversal + 1
          if (w1_val != POS_LABEL && w2_val == POS_LABEL)
            n_incident <- n_incident + 1
        }
      } else {
        n_w1_only <- n_w1_only + 1
      }
    }

    merged_row[[col]] <- merge_tori_wave2(vals, surveys)
  }

  # 年龄列：取较大值（等价于Wave2的年龄）
  for (col in age_cols) {
    merged_row[[col]] <- safe_max_numeric(rows[[col]])
  }

  # survey列标记为merged
  merged_row$survey_child  <- "merged"
  merged_row$survey_parent <- "merged"

  # 其余列一致性检查
  for (col in other_cols) {
    vals <- unique(rows[[col]])
    if (length(vals) > 1) {
      cat(sprintf("警告：关系%s 在列%s 不一致：%s\n",
                  k, col, paste(vals, collapse=" vs ")))
    }
  }

  merged_list[[k]] <- merged_row
}

df_merged <- do.call(rbind, merged_list)
rownames(df_merged) <- NULL

# 拼回完整数据集
df_final <- rbind(df_unique, df_merged)
rownames(df_final) <- NULL

cat("\n=== 最终汇总 ===\n")
cat("原始df_clean.csv行数：              ", nrow(df), "\n")
cat("剔除跨波次错配后：                   ", nrow(df_wave_matched), "\n")
cat("最终去重后（分析用数据集）：          ", nrow(df_final), "\n")
cat("\n--- 合并策略统计（重复亲子对）---\n")
cat(sprintf("  Wave2阳性（含新发）：%d\n",  n_w2_pos))
cat(sprintf("  Wave2阴性：          %d\n",  n_w2_neg))
cat(sprintf("  其中逆转（W1阳→W2阴）：%d（此部分OR合并会高估为阳性）\n", n_reversal))
cat(sprintf("  其中新发（W1阴→W2阳）：%d（此部分两策略一致，均为阳性）\n", n_incident))
cat(sprintf("  只有Wave1记录：      %d\n",  n_w1_only))

write.csv(df_final, "df_sensitivity_wave2.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\n已导出：df_sensitivity_wave2.csv\n")
