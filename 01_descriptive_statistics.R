# Descriptive statistics for the parent-offspring sample.
# Input: df_final_corrected.csv
# Output files contain the main descriptive table and supporting summaries.

# 1. Read the corrected data
df <- read.csv("df_final_corrected.csv", fileEncoding = "UTF-8", stringsAsFactors = FALSE)

cat("=== Basic data information ===\n")
cat("Total rows (corrected parent-offspring pairs): ", nrow(df), "\n")
cat("Column names:\n")
print(colnames(df))

# Count merged rows (survey_child_bin = NA indicates records merged across both waves)
n_merged_rows <- sum(is.na(df$survey_child_bin))
n_single_wave <- nrow(df) - n_merged_rows
cat(sprintf("\nSingle-wave records: %d rows; records merged across both waves (survey = NA): %d rows\n",
            n_single_wave, n_merged_rows))

# 2. Confirm that required variables exist; recreate them if needed
if (!"tori_child_bin" %in% colnames(df)) {
  df$tori_child_bin  <- ifelse(df$tori_child  == "有", 1, 0)
}
if (!"tori_parent_bin" %in% colnames(df)) {
  df$tori_parent_bin <- ifelse(df$tori_parent == "有", 1, 0)
}
if (!"pair_type" %in% colnames(df)) {
  df$pair_type <- paste(
    ifelse(df$relation  == "母", "Mother", "Father"),
    ifelse(df$sex_child == "F",  "Daughter", "Son"),
    sep = "-"
  )
}

# ============================================================
# 3. Offspring statistics
# ============================================================
cat("\n=== Offspring ===\n")

n_off          <- nrow(df)
age_off_mean   <- mean(df$age_child, na.rm = TRUE)
age_off_sd     <- sd(df$age_child,   na.rm = TRUE)
n_off_female   <- sum(df$sex_child == "F", na.rm = TRUE)
pct_off_female <- n_off_female / n_off * 100
n_off_tori     <- sum(df$tori_child == "有", na.rm = TRUE)
pct_off_tori   <- n_off_tori / n_off * 100

cat(sprintf("N = %d\n", n_off))
cat(sprintf("Age: %.1f ± %.1f\n", age_off_mean, age_off_sd))
cat(sprintf("Female: %d (%.1f%%)\n", n_off_female, pct_off_female))
cat(sprintf("Tori prevalence: %d (%.1f%%)\n", n_off_tori, pct_off_tori))

# Stratify offspring oral exostosis by sex
tori_by_sex_child <- tapply(
  df$tori_child == "有", df$sex_child,
  function(x) round(mean(x, na.rm = TRUE) * 100, 1)
)
cat("\nOffspring oral exostosis prevalence by sex:\n")
print(tori_by_sex_child)

# ============================================================
# 4. Parent statistics
# ============================================================
cat("\n=== Mothers ===\n")
df_mother    <- df[df$relation == "母", ]
n_mom        <- nrow(df_mother)
age_mom_mean <- mean(df_mother$age_parent, na.rm = TRUE)
age_mom_sd   <- sd(df_mother$age_parent,   na.rm = TRUE)
n_mom_tori   <- sum(df_mother$tori_parent == "有", na.rm = TRUE)
pct_mom_tori <- n_mom_tori / n_mom * 100

cat(sprintf("N = %d\n", n_mom))
cat(sprintf("Age: %.1f ± %.1f\n", age_mom_mean, age_mom_sd))
cat(sprintf("Tori prevalence: %d (%.1f%%)\n", n_mom_tori, pct_mom_tori))

cat("\n=== Fathers ===\n")
df_father    <- df[df$relation == "父", ]
n_dad        <- nrow(df_father)
age_dad_mean <- mean(df_father$age_parent, na.rm = TRUE)
age_dad_sd   <- sd(df_father$age_parent,   na.rm = TRUE)
n_dad_tori   <- sum(df_father$tori_parent == "有", na.rm = TRUE)
pct_dad_tori <- n_dad_tori / n_dad * 100

cat(sprintf("N = %d\n", n_dad))
cat(sprintf("Age: %.1f ± %.1f\n", age_dad_mean, age_dad_sd))
cat(sprintf("Tori prevalence: %d (%.1f%%)\n", n_dad_tori, pct_dad_tori))

# ============================================================
# 5. Parent-offspring pair-type statistics
# ============================================================
cat("\n=== Parent-offspring pair types ===\n")

pair_tab <- table(df$pair_type)
pair_pct <- prop.table(pair_tab) * 100
pair_summary <- data.frame(
  pair_type = names(pair_tab),
  n         = as.integer(pair_tab),
  pct       = round(as.numeric(pair_pct), 1),
  stringsAsFactors = FALSE
)
print(pair_summary)
cat("Total pairs:", nrow(df), "\n")

# ============================================================
# 6. Prevalence by survey wave (single-wave records only; merged rows excluded)
# ============================================================
cat("\n=== Prevalence by wave (single-wave records only; merged rows excluded) ===\n")

df_single <- df[!is.na(df$survey_child_bin), ]
cat(sprintf("Rows used for wave-stratified analysis: %d (%d merged rows excluded)\n",
            nrow(df_single), n_merged_rows))

cat("Distribution of survey_child values in df_single (should contain no NAs):\n")
print(table(df_single$survey_child, useNA = "always"))

# Offspring stratified by survey wave
survey_child_tab <- table(df_single$survey_child, df_single$tori_child)
cat("\nOffspring oral exostosis status by survey wave:\n")
print(survey_child_tab)
print(round(prop.table(survey_child_tab, margin = 1) * 100, 1))
chisq_child <- chisq.test(survey_child_tab)
print(chisq_child)

# Parents stratified by survey wave
survey_parent_tab <- table(df_single$survey_parent, df_single$tori_parent)
cat("\nParent oral exostosis status by survey wave:\n")
print(survey_parent_tab)
print(round(prop.table(survey_parent_tab, margin = 1) * 100, 1))
chisq_parent <- chisq.test(survey_parent_tab)
print(chisq_parent)

# ============================================================
# 7. Export each result to a separate CSV file
# ============================================================

# --- CSV 1: Main descriptive-statistics table ---
table1 <- data.frame(
  group      = c("Offspring", "Mother", "Father"),
  n          = c(n_off, n_mom, n_dad),
  age_mean   = round(c(age_off_mean, age_mom_mean, age_dad_mean), 1),
  age_sd     = round(c(age_off_sd,   age_mom_sd,   age_dad_sd),   1),
  female_n   = c(n_off_female, n_mom, 0),
  female_pct = round(c(pct_off_female, 100.0, 0.0), 1),
  tori_n     = c(n_off_tori, n_mom_tori, n_dad_tori),
  tori_pct   = round(c(pct_off_tori, pct_mom_tori, pct_dad_tori), 1),
  stringsAsFactors = FALSE
)
write.csv(table1, "result_step1_table1.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("\n✓ result_step1_table1.csv\n")

# --- CSV 2: Distribution of parent-offspring pair types ---
write.csv(pair_summary, "result_step1_pair_types.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("✓ result_step1_pair_types.csv\n")

# --- CSV 3: Offspring oral exostosis prevalence by sex ---
tori_sex_df <- data.frame(
  sex        = names(tori_by_sex_child),
  tori_pct   = as.vector(tori_by_sex_child),
  stringsAsFactors = FALSE
)
write.csv(tori_sex_df, "result_step1_tori_by_sex.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("✓ result_step1_tori_by_sex.csv\n")

# --- CSV 4: Offspring prevalence by survey wave ---
survey_child_long <- as.data.frame(survey_child_tab, stringsAsFactors = FALSE)
colnames(survey_child_long) <- c("wave", "tori_status", "count")
survey_child_pct  <- as.data.frame(
  round(prop.table(survey_child_tab, margin = 1) * 100, 1),
  stringsAsFactors = FALSE
)
colnames(survey_child_pct) <- c("wave", "tori_status", "pct")
survey_child_out <- merge(survey_child_long, survey_child_pct, by = c("wave", "tori_status"))
survey_child_out$chisq_p <- ifelse(
  survey_child_out$wave == survey_child_out$wave[1] &
  survey_child_out$tori_status == survey_child_out$tori_status[1],
  round(chisq_child$p.value, 4), ""
)
# Write the p-value in the first row only; leave the remaining rows blank
survey_child_out$chisq_p <- ""
survey_child_out$chisq_p[1] <- round(chisq_child$p.value, 4)

write.csv(survey_child_out, "result_step1_wave_offspring.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("✓ result_step1_wave_offspring.csv\n")

# --- CSV 5: Parent prevalence by survey wave ---
survey_parent_long <- as.data.frame(survey_parent_tab, stringsAsFactors = FALSE)
colnames(survey_parent_long) <- c("wave", "tori_status", "count")
survey_parent_pct  <- as.data.frame(
  round(prop.table(survey_parent_tab, margin = 1) * 100, 1),
  stringsAsFactors = FALSE
)
colnames(survey_parent_pct) <- c("wave", "tori_status", "pct")
survey_parent_out <- merge(survey_parent_long, survey_parent_pct, by = c("wave", "tori_status"))
survey_parent_out$chisq_p <- ""
survey_parent_out$chisq_p[1] <- round(chisq_parent$p.value, 4)

write.csv(survey_parent_out, "result_step1_wave_parent.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("✓ result_step1_wave_parent.csv\n")

cat("\n=== Step 1 completed ===\n")
cat(sprintf("Total sample size based on the corrected data set: N = %d\n", nrow(df)))
cat(sprintf("Of these, %d rows were merged across both waves and excluded from the wave-stratified analysis.\n", n_merged_rows))
