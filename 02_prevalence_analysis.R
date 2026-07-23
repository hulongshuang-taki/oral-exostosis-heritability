# Prevalence analyses by age, sex, and parental relationship.
# Input: df_final_corrected.csv

# 1. Read the data
df <- read.csv("df_final_corrected.csv",
               fileEncoding = "UTF-8", stringsAsFactors = FALSE)
cat("总行数：", nrow(df), "\n")

# 2. Create binary variables (rebuild them without relying on existing columns)
df$tori_child_bin  <- as.integer(df$tori_child  == "1. あり")
df$tori_parent_bin <- as.integer(df$tori_parent == "1. あり")
df$sex_child_bin   <- as.integer(df$sex_child   == "F")   # F = 1, M = 0 (reference: male)
df$relation_bin    <- as.integer(df$relation    == "母")  # Mother = 1, father = 0 (reference: father)

# Verify the coding after execution to ensure there are no NAs or unexpected values
cat("\n--- 编码确认 ---\n")
cat("tori_child_bin:\n");  print(table(df$tori_child,  df$tori_child_bin,  useNA="always"))
cat("tori_parent_bin:\n"); print(table(df$tori_parent, df$tori_parent_bin, useNA="always"))
cat("sex_child_bin:\n");   print(table(df$sex_child,   df$sex_child_bin))
cat("relation_bin:\n");    print(table(df$relation,    df$relation_bin))

# ============================================================
# 3. Age-stratified prevalence
# Offspring: ages 20–59 (exclude age 60+ because of the small sample size)
# Parents: ages 40–79 (exclude both extremes because of small sample sizes)
# ============================================================
calc_age_prev <- function(data, age_col, tori_bin_col, id_col,
                          age_min, age_max, breaks, labels, group_name) {
  sub <- data[data[[age_col]] >= age_min & data[[age_col]] < age_max, ]
  cat(sprintf("%s：保留%d行（%d–%d岁），排除%d行\n",
              group_name, nrow(sub), age_min, age_max-1,
              nrow(data) - nrow(sub)))
  sub$age_grp <- cut(sub[[age_col]], breaks = breaks,
                     right = FALSE, labels = labels)
  # Use tapply to calculate n and n_tori at once, avoiding row-order errors from merge
  n_vec     <- tapply(sub[[id_col]],       sub$age_grp, length)
  ntori_vec <- tapply(sub[[tori_bin_col]], sub$age_grp, sum, na.rm = TRUE)
  data.frame(
    age_group  = labels,
    n          = as.integer(n_vec),
    n_tori     = as.integer(ntori_vec),
    prevalence = round(as.integer(ntori_vec) / as.integer(n_vec) * 100, 1),
    group      = group_name,
    stringsAsFactors = FALSE
  )
}

cat("\n--- 年龄分层患病率 ---\n")
age_child <- calc_age_prev(
  df, "age_child", "tori_child_bin", "id_child",
  20, 60, seq(20, 60, 10), c("20s","30s","40s","50s"), "Offspring"
)
age_parent <- calc_age_prev(
  df, "age_parent", "tori_parent_bin", "id_parent",
  40, 80, seq(40, 80, 10), c("40s","50s","60s","70s"), "Parent"
)
age_result_all <- rbind(age_child, age_parent)
print(age_result_all)

# ============================================================
# 4. Compare prevalence by offspring sex and parental relationship (chi-squared tests)
# Offspring: F vs M
# Parents: mother vs father (do not use sex_parent because it is perfectly collinear with relationship)
# ============================================================
cat("\n--- 子代 性别比较 ---\n")
tab_child  <- table(sex      = df$sex_child, tori = df$tori_child_bin)
chisq_child <- chisq.test(tab_child)
print(tab_child); print(chisq_child)

cat("\n--- 亲代 父母比较 ---\n")
tab_parent  <- table(relation = df$relation,  tori = df$tori_parent_bin)
chisq_parent <- chisq.test(tab_parent)
print(tab_parent); print(chisq_parent)

# Format the results as an output table
calc_prev_row <- function(data, filter_col, filter_val, tori_col) {
  sub <- data[data[[filter_col]] == filter_val, ]
  n   <- nrow(sub)
  nt  <- sum(sub[[tori_col]], na.rm = TRUE)
  c(n = n, n_tori = nt, prevalence = round(nt / n * 100, 1))
}

r_F   <- calc_prev_row(df, "sex_child", "F", "tori_child_bin")
r_M   <- calc_prev_row(df, "sex_child", "M", "tori_child_bin")
r_mom <- calc_prev_row(df, "relation",  "母", "tori_parent_bin")
r_dad <- calc_prev_row(df, "relation",  "父", "tori_parent_bin")

sex_comparison <- data.frame(
  group      = c("Offspring_Female", "Offspring_Male",
                 "Parent_Mother",    "Parent_Father"),
  n          = as.integer(c(r_F["n"],   r_M["n"],   r_mom["n"],   r_dad["n"])),
  n_tori     = as.integer(c(r_F["n_tori"], r_M["n_tori"],
                             r_mom["n_tori"], r_dad["n_tori"])),
  prevalence = c(r_F["prevalence"], r_M["prevalence"],
                 r_mom["prevalence"], r_dad["prevalence"]),
  chisq_p    = as.character(c(round(chisq_child$p.value,  4), "",
                               round(chisq_parent$p.value, 4), "")),
  stringsAsFactors = FALSE
)
print(sex_comparison)

# ============================================================
# 5. Compare parent and offspring prevalence (McNemar's test)
# Use binary columns to avoid coding issues in the original character columns
# ============================================================
cat("\n--- McNemar検験（亲代 vs 子代）---\n")
prev_child  <- mean(df$tori_child_bin,  na.rm = TRUE) * 100
prev_parent <- mean(df$tori_parent_bin, na.rm = TRUE) * 100
cat(sprintf("Offspring: %.1f%%  Parent: %.1f%%\n", prev_child, prev_parent))

mcnemar_tab    <- table(offspring = df$tori_child_bin, parent = df$tori_parent_bin)
mcnemar_result <- mcnemar.test(mcnemar_tab)
print(mcnemar_tab); print(mcnemar_result)

# Stratify by offspring sex
for (sx in c("F", "M")) {
  sub <- df[df$sex_child == sx, ]
  cat(sprintf("Sex=%s | Offspring: %.1f%%  Parent: %.1f%%\n", sx,
              mean(sub$tori_child_bin,  na.rm = TRUE) * 100,
              mean(sub$tori_parent_bin, na.rm = TRUE) * 100))
}

prev_comparison <- data.frame(
  group      = c("Offspring", "Parent"),
  prevalence = round(c(prev_child, prev_parent), 1),
  mcnemar_p  = as.character(c(round(mcnemar_result$p.value, 4), "")),
  stringsAsFactors = FALSE
)

# ============================================================
# 6. Univariable logistic regression
# Offspring: age_child and sex_child_bin (F = 1; reference: male)
# Parents: age_parent and relation_bin (mother = 1; reference: father)
# ============================================================
cat("\n--- 单因素logistic回归 ---\n")

fit_glm <- function(formula, data) {
  m  <- glm(formula, data = data, family = binomial(link = "logit"))
  or <- exp(coef(m))
  ci <- exp(confint.default(m))
  p  <- summary(m)$coefficients[, 4]
  list(or = or, ci = ci, p = p)
}

g1 <- fit_glm(tori_child_bin  ~ age_child,    df)
g2 <- fit_glm(tori_child_bin  ~ sex_child_bin, df)
g3 <- fit_glm(tori_parent_bin ~ age_parent,   df)
g4 <- fit_glm(tori_parent_bin ~ relation_bin,  df)

logistic_result <- data.frame(
  group    = c("Offspring", "Offspring", "Parent", "Parent"),
  variable = c("Age (per 1 year)",
               "Sex (Female vs Male, ref=Male)",
               "Age (per 1 year)",
               "Relation (Mother vs Father, ref=Father)"),
  OR      = round(c(g1$or[2], g2$or[2], g3$or[2], g4$or[2]), 3),
  CI_low  = round(c(g1$ci[2,1], g2$ci[2,1], g3$ci[2,1], g4$ci[2,1]), 3),
  CI_high = round(c(g1$ci[2,2], g2$ci[2,2], g3$ci[2,2], g4$ci[2,2]), 3),
  p_value = round(c(g1$p[2],    g2$p[2],    g3$p[2],    g4$p[2]),    4),
  stringsAsFactors = FALSE
)
print(logistic_result)

# ============================================================
# 7. Export CSV files
# ============================================================
write.csv(age_result_all,   "result_step2_age_prevalence.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sex_comparison,   "result_step2_sex_comparison.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(prev_comparison,  "result_step2_parent_vs_child.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(logistic_result,  "result_step2_logistic.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\n=== Step 2 完成 ===\n")
cat("✓ result_step2_age_prevalence.csv\n")
cat("✓ result_step2_sex_comparison.csv\n")
cat("✓ result_step2_parent_vs_child.csv\n")
cat("✓ result_step2_logistic.csv\n")
cat("注：子代性别ref=男性（OR<1=女性风险低）；亲代比较ref=父亲（OR<1=母亲患病率低）\n")
