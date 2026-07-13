# Multivariable logistic regression and stratified analyses.
# Input: df_final_corrected.csv

# ① 读取数据
df <- read.csv("df_final_corrected.csv",
               fileEncoding = "UTF-8", stringsAsFactors = FALSE)
cat("总行数：", nrow(df), "\n")

# ② 强制重建所有binary变量
df$tori_child_bin  <- as.integer(df$tori_child  == "1. あり")
df$tori_parent_bin <- as.integer(df$tori_parent == "1. あり")
df$sex_child_bin   <- as.integer(df$sex_child   == "F")   # F=1, M=0（ref=男性）
df$relation_bin    <- as.integer(df$relation    == "母")  # 母=1, 父=0

df$pair_type <- paste(
  ifelse(df$relation  == "母", "Mother", "Father"),
  ifelse(df$sex_child == "F",  "Daughter", "Son"),
  sep = "-"
)

# ③ survey wave编码
# 先打印实际取值确认
cat("\nsurvey_child实际取值：\n")
print(table(df$survey_child, useNA = "always"))

# 【需确认】根据上方输出修改这一行的匹配值
# 常见情况：
#   "1次調査" / "2次調査"        → c("2次調査")
#   "1次" / "2次" / "merged"    → c("2次", "merged")
#   "一次" / "二次"              → c("二次")
WAVE2_VALUES <- c("2次調査")  # ← 根据实际取值修改这里

df$survey_child_bin <- as.integer(df$survey_child %in% WAVE2_VALUES)
cat("survey_child_bin编码结果（0=Wave1, 1=Wave2）：\n")
print(table(df$survey_child, df$survey_child_bin, useNA = "always"))

# ④ 编码确认
cat("\n--- 编码确认 ---\n")
cat("tori_child_bin:\n")
print(table(df$tori_child, df$tori_child_bin, useNA = "always"))
cat("tori_parent_bin:\n")
print(table(df$tori_parent, df$tori_parent_bin, useNA = "always"))

# ⑤ 完整案例筛选
vars_needed <- c("tori_child_bin", "tori_parent_bin",
                 "age_child", "age_parent", "sex_child_bin", "survey_child_bin")
na_counts <- sapply(df[vars_needed], function(x) sum(is.na(x)))
cat("\n--- 各变量NA数量 ---\n"); print(na_counts)

df_complete <- df[complete.cases(df[vars_needed]), ]
cat(sprintf("\n完整案例：%d / %d（%.1f%%）\n",
            nrow(df_complete), nrow(df),
            nrow(df_complete) / nrow(df) * 100))

# ⑥ 四个嵌套模型
# M0：粗OR
# M1：+ 子代年龄 + 子代性别
# M2：+ 亲代年龄（主模型，论文用）
# M3：+ survey wave（用于判断是否需要纳入，不作为主模型）
cat("\n--- 构建模型 ---\n")
glm_m0 <- glm(tori_child_bin ~ tori_parent_bin,
              data = df_complete, family = binomial())
glm_m1 <- glm(tori_child_bin ~ tori_parent_bin + age_child + sex_child_bin,
              data = df_complete, family = binomial())
glm_m2 <- glm(tori_child_bin ~ tori_parent_bin + age_child + sex_child_bin + age_parent,
              data = df_complete, family = binomial())
glm_m3 <- glm(tori_child_bin ~ tori_parent_bin + age_child + sex_child_bin + age_parent +
                survey_child_bin,
              data = df_complete, family = binomial())

# 亲代OR变化
or_progression <- sapply(list(glm_m0, glm_m1, glm_m2, glm_m3),
                         function(m) round(exp(coef(m)["tori_parent_bin"]), 3))
cat("\n亲代骨隆起OR变化：\n")
cat(sprintf("  M0 粗OR                    : %.3f\n", or_progression[1]))
cat(sprintf("  M1 + 子代年龄/性别         : %.3f\n", or_progression[2]))
cat(sprintf("  M2 + 亲代年龄（主模型）    : %.3f\n", or_progression[3]))
cat(sprintf("  M3 + survey wave           : %.3f\n", or_progression[4]))

# M3中wave效应（判断用）
cat("\n--- M3 survey wave效应（判断是否纳入）---\n")
coef_m3 <- summary(glm_m3)$coefficients
or_m3   <- exp(coef(glm_m3))
ci_m3   <- exp(confint.default(glm_m3))
cat(sprintf("Survey wave aOR = %.3f (95%% CI: %.3f–%.3f), p = %.4f\n",
            or_m3["survey_child_bin"],
            ci_m3["survey_child_bin", 1],
            ci_m3["survey_child_bin", 2],
            coef_m3["survey_child_bin", 4]))
cat("→ p>0.05且AIC不降低则最终模型用M2；否则改用M3\n")

# ⑦ 主模型（M2）详细结果
cat("\n--- 主模型M2详细结果 ---\n")
coef_m2 <- summary(glm_m2)$coefficients
or_m2   <- exp(coef(glm_m2))
ci_m2   <- exp(confint.default(glm_m2))

label_map <- c(
  "tori_parent_bin" = "Parental oral tori (yes vs no)",
  "age_child"       = "Offspring age (per 1 year)",
  "sex_child_bin"   = "Offspring sex (female vs male, ref=male)",
  "age_parent"      = "Parental age (per 1 year)"
)
result_m2 <- data.frame(
  variable = sapply(rownames(coef_m2)[-1], function(nm)
               ifelse(nm %in% names(label_map), label_map[nm], nm)),
  aOR      = round(or_m2[-1],      3),
  CI_low   = round(ci_m2[-1, 1],  3),
  CI_high  = round(ci_m2[-1, 2],  3),
  p_value  = round(coef_m2[-1, 4], 4),
  stringsAsFactors = FALSE,
  row.names = NULL
)
print(result_m2)
cat(sprintf("\n【核心】亲代骨隆起 aOR = %.3f (95%% CI: %.3f–%.3f), p = %.4f\n",
            result_m2$aOR[1], result_m2$CI_low[1],
            result_m2$CI_high[1], result_m2$p_value[1]))

# ⑧ 四种亲子对分层aOR
cat("\n--- 分层aOR ---\n")
pair_types <- c("Mother-Daughter", "Mother-Son", "Father-Daughter", "Father-Son")
result_stratified_df <- do.call(rbind, lapply(pair_types, function(pt) {
  sub <- df_complete[df_complete$pair_type == pt, ]
  if (nrow(sub) < 30) { cat(sprintf("[%s] 样本量过小，跳过\n", pt)); return(NULL) }
  tryCatch({
    m    <- glm(tori_child_bin ~ tori_parent_bin + age_child + age_parent,
                data = sub, family = binomial())
    or_s <- exp(coef(m))
    ci_s <- exp(confint.default(m))
    p_s  <- summary(m)$coefficients["tori_parent_bin", 4]
    cat(sprintf("[%s] n=%d, aOR=%.3f (%.3f–%.3f), p=%.4f\n",
                pt, nrow(sub), or_s["tori_parent_bin"],
                ci_s["tori_parent_bin",1], ci_s["tori_parent_bin",2], p_s))
    data.frame(pair_type=pt, n=nrow(sub),
               aOR    =round(or_s["tori_parent_bin"],     3),
               CI_low =round(ci_s["tori_parent_bin", 1], 3),
               CI_high=round(ci_s["tori_parent_bin", 2], 3),
               p_value=round(p_s, 4), stringsAsFactors=FALSE)
  }, error=function(e){ cat(sprintf("[%s] 模型失败: %s\n", pt, e$message)); NULL })
}))
rownames(result_stratified_df) <- NULL

# ⑨ 模型拟合优度
nagelkerke_r2 <- function(model) {
  n      <- nrow(model$model)
  log_L0 <- -model$null.deviance / 2
  log_L1 <- -deviance(model) / 2
  r2_cs  <- 1 - exp((2/n) * (log_L0 - log_L1))
  r2_max <- 1 - exp((2/n) * log_L0)
  round(r2_cs / r2_max, 4)
}

model_fit <- data.frame(
  model = c("M0: tori_parent only",
            "M1: + age_child + sex_child",
            "M2: + age_parent",
            "M3: + survey_wave"),
  AIC           = round(c(AIC(glm_m0), AIC(glm_m1), AIC(glm_m2), AIC(glm_m3)), 2),
  BIC           = round(c(BIC(glm_m0), BIC(glm_m1), BIC(glm_m2), BIC(glm_m3)), 2),
  Nagelkerke_R2 = c(nagelkerke_r2(glm_m0), nagelkerke_r2(glm_m1),
                    nagelkerke_r2(glm_m2), nagelkerke_r2(glm_m3)),
  stringsAsFactors = FALSE
)
cat("\n--- 模型拟合比较 ---\n"); print(model_fit)

or_prog_df <- data.frame(
  model     = c("M0: Crude OR", "M1: + age/sex (offspring)",
                "M2: + age (parent)", "M3: + survey wave"),
  OR_parent = unname(or_progression),
  stringsAsFactors = FALSE
)

# ⑩ 导出
write.csv(result_m2,            "result_step5_main_model.csv",      row.names=FALSE, fileEncoding="UTF-8")
write.csv(or_prog_df,           "result_step5_or_progression.csv",   row.names=FALSE, fileEncoding="UTF-8")
write.csv(result_stratified_df, "result_step5_stratified_aOR.csv",   row.names=FALSE, fileEncoding="UTF-8")
write.csv(model_fit,            "result_step5_model_fit.csv",        row.names=FALSE, fileEncoding="UTF-8")

cat("\n=== Step 5 完成 ===\n")
cat("✓ result_step5_main_model.csv\n")
cat("✓ result_step5_or_progression.csv\n")
cat("✓ result_step5_stratified_aOR.csv\n")
cat("✓ result_step5_model_fit.csv\n")
cat("\n【提醒】运行后请确认：\n")
cat("  1. survey_child实际取值是否和WAVE2_VALUES匹配\n")
cat("  2. M3中wave的p值和AIC，决定最终模型用M2还是M3\n")
