# Heritability and sex-specific transmission analyses.
# The tetrachoric estimate is the primary analysis; the Falconer estimate is
# retained as a sensitivity analysis. Requires the polycor package.

# 1. Read the data and rebuild the binary variables
df <- read.csv("df_final_corrected.csv",
               fileEncoding = "UTF-8", stringsAsFactors = FALSE)
cat("Total rows: ", nrow(df), "\n")

df$tori_child_bin  <- as.integer(df$tori_child  == "1. あり")
df$tori_parent_bin <- as.integer(df$tori_parent == "1. あり")
df$sex_child_bin   <- as.integer(df$sex_child   == "F")  # F=1, M=0
df$relation_bin    <- as.integer(df$relation    == "母") # Mother = 1, father = 0
df$pair_type <- paste(
  ifelse(df$relation  == "母", "Mother", "Father"),
  ifelse(df$sex_child == "F",  "Daughter", "Son"),
  sep = "-"
)

# Verify the coding
cat("\n--- Coding verification ---\n")
cat("tori_child_bin:\n");  print(table(df$tori_child,  df$tori_child_bin,  useNA="always"))
cat("tori_parent_bin:\n"); print(table(df$tori_parent, df$tori_parent_bin, useNA="always"))

# ============================================================
# 2. Falconer liability-threshold model
# Formula: h² = 2 × r_falconer
# r_falconer = (K_R - K_o) × t_p / (z_p × (1 - K_o))
#   K_p = parent prevalence; K_o = offspring prevalence
#   K_R = offspring prevalence among affected parents (conditional prevalence)
#   t_p = qnorm(1 - K_p)，z_p = dnorm(t_p)
# Resample at the pair level for the bootstrap to preserve parent-offspring pairing
# ============================================================
falconer_h2 <- function(parent_bin, child_bin, label = "", n_boot = 1000) {
  valid  <- !is.na(parent_bin) & !is.na(child_bin)
  p_vec  <- parent_bin[valid]
  o_vec  <- child_bin[valid]
  n      <- sum(valid)

  K_p <- mean(p_vec)
  K_o <- mean(o_vec)

  # Boundary check: the Falconer formula fails when prevalence is too extreme
  if (K_p <= 0 || K_p >= 1 || K_o <= 0 || K_o >= 1) {
    cat(sprintf("Warning: %s prevalence is outside the valid range (K_p=%.3f, K_o=%.3f); skipping.\n",
                label, K_p, K_o))
    return(data.frame(group=label, n=n, K_parent=NA, K_offspring=NA,
                      r_falconer=NA, h2=NA, CI_low=NA, CI_high=NA))
  }

  t_p <- qnorm(1 - K_p)
  z_p <- dnorm(t_p)
  K_R <- mean(o_vec[p_vec == 1])  # Offspring prevalence among affected parents

  r  <- (K_R - K_o) * t_p / (z_p * (1 - K_o))
  h2 <- 2 * r

  # Bootstrap by resampling pairs to preserve parent-offspring pairing
  h2_boot <- replicate(n_boot, {
    idx   <- sample(n, n, replace = TRUE)
    p_b   <- p_vec[idx]; o_b <- o_vec[idx]
    K_p_b <- mean(p_b);  K_o_b <- mean(o_b)
    if (K_p_b <= 0 || K_p_b >= 1 ||
        K_o_b <= 0 || K_o_b >= 1 ||
        sum(p_b == 1) == 0) return(NA)
    t_b   <- qnorm(1 - K_p_b)
    z_b   <- dnorm(t_b)
    K_R_b <- mean(o_b[p_b == 1])
    2 * (K_R_b - K_o_b) * t_b / (z_b * (1 - K_o_b))
  })
  h2_boot <- h2_boot[!is.na(h2_boot)]

  data.frame(
    group       = label,
    n           = n,
    K_parent    = round(K_p * 100, 1),
    K_offspring = round(K_o * 100, 1),
    r_falconer  = round(r,  4),
    h2          = round(h2, 4),
    CI_low      = round(quantile(h2_boot, 0.025), 4),
    CI_high     = round(quantile(h2_boot, 0.975), 4),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 3. Tetrachoric-correlation method
# h² = 2 × r_tetrachoric
# CI: propagate the SE returned by polychor() using the delta method
# ============================================================
library(polycor)

tetra_h2 <- function(parent_bin, child_bin, label = "") {
  valid <- !is.na(parent_bin) & !is.na(child_bin)
  p_f   <- factor(parent_bin[valid], levels = c(0, 1))
  o_f   <- factor(child_bin[valid],  levels = c(0, 1))
  n     <- sum(valid)

  tryCatch({
    res     <- polychor(p_f, o_f, std.err = TRUE)
    r       <- res$rho
    se      <- sqrt(res$var[1, 1])
    h2      <- 2 * r
    CI_low  <- 2 * (r - 1.96 * se)
    CI_high <- 2 * (r + 1.96 * se)
    data.frame(group=label, n=n,
               r_tetra=round(r, 4), h2=round(h2, 4),
               CI_low=round(CI_low, 4), CI_high=round(CI_high, 4),
               stringsAsFactors=FALSE)
  }, error = function(e) {
    cat(sprintf("Warning: Tetrachoric analysis failed for %s - %s\n", label, e$message))
    data.frame(group=label, n=n,
               r_tetra=NA, h2=NA, CI_low=NA, CI_high=NA,
               stringsAsFactors=FALSE)
  })
}

# ============================================================
# 4. Z-test function for comparing two h² estimates
# Approximate SE = (CI_high - CI_low) / (2 × 1.96)
# ============================================================
z_test_h2 <- function(h2_a, h2_b) {
  se_a <- (h2_a$CI_high - h2_a$CI_low) / (2 * 1.96)
  se_b <- (h2_b$CI_high - h2_b$CI_low) / (2 * 1.96)
  Z    <- (h2_a$h2 - h2_b$h2) / sqrt(se_a^2 + se_b^2)
  p    <- 2 * pnorm(-abs(Z))
  list(Z = round(Z, 3), p = round(p, 4))
}

# ============================================================
# 5. Define grouping indices
# ============================================================
idx_mother   <- df$relation  == "母"
idx_father   <- df$relation  == "父"
idx_daughter <- df$sex_child == "F"
idx_son      <- df$sex_child == "M"
pair_types   <- c("Mother-Daughter","Mother-Son","Father-Daughter","Father-Son")

# ============================================================
# 6. Run the Falconer analysis
# ============================================================
cat("\n--- Falconer method (1,000 bootstrap iterations; please wait) ---\n")
set.seed(42)  # Global seed to make the bootstrap reproducible across all groups
f_overall  <- falconer_h2(df$tori_parent_bin, df$tori_child_bin, "Overall")
f_mother   <- falconer_h2(df$tori_parent_bin[idx_mother],   df$tori_child_bin[idx_mother],   "Mother")
f_father   <- falconer_h2(df$tori_parent_bin[idx_father],   df$tori_child_bin[idx_father],   "Father")
f_daughter <- falconer_h2(df$tori_parent_bin[idx_daughter], df$tori_child_bin[idx_daughter], "Female offspring")
f_son      <- falconer_h2(df$tori_parent_bin[idx_son],      df$tori_child_bin[idx_son],      "Male offspring")

f_pairs <- lapply(pair_types, function(pt) {
  idx <- df$pair_type == pt
  if (sum(idx) < 50) return(NULL)
  falconer_h2(df$tori_parent_bin[idx], df$tori_child_bin[idx], pt)
})

h2_falconer <- do.call(rbind, c(
  list(f_overall, f_mother, f_father, f_daughter, f_son),
  f_pairs[!sapply(f_pairs, is.null)]
))
cat("\n=== Falconer results ===\n"); print(h2_falconer)

# ============================================================
# 7. Run the tetrachoric analysis
# ============================================================
cat("\n--- Tetrachoric method ---\n")
t_overall  <- tetra_h2(df$tori_parent_bin, df$tori_child_bin, "Overall")
t_mother   <- tetra_h2(df$tori_parent_bin[idx_mother],   df$tori_child_bin[idx_mother],   "Mother")
t_father   <- tetra_h2(df$tori_parent_bin[idx_father],   df$tori_child_bin[idx_father],   "Father")
t_daughter <- tetra_h2(df$tori_parent_bin[idx_daughter], df$tori_child_bin[idx_daughter], "Female offspring")
t_son      <- tetra_h2(df$tori_parent_bin[idx_son],      df$tori_child_bin[idx_son],      "Male offspring")

t_pairs <- lapply(pair_types, function(pt) {
  idx <- df$pair_type == pt
  tetra_h2(df$tori_parent_bin[idx], df$tori_child_bin[idx], pt)
})

h2_tetra <- do.call(rbind, c(
  list(t_overall, t_mother, t_father, t_daughter, t_son),
  t_pairs
))
cat("\n=== Tetrachoric results ===\n"); print(h2_tetra)

# ============================================================
# 8. Z-tests: maternal vs paternal and female vs male offspring
# Use the tetrachoric results (primary method)
# ============================================================
cat("\n--- Z-tests ---\n")
z_parent   <- z_test_h2(t_mother,   t_father)
z_offspring <- z_test_h2(t_daughter, t_son)

cat(sprintf("Maternal (%.4f) vs paternal (%.4f): Z=%.3f, p=%.4f\n",
            t_mother$h2, t_father$h2, z_parent$Z, z_parent$p))
cat(sprintf("Female offspring (%.4f) vs male offspring (%.4f): Z=%.3f, p=%.4f\n",
            t_daughter$h2, t_son$h2, z_offspring$Z, z_offspring$p))

z_result <- data.frame(
  comparison = c("Mother vs Father", "Female vs Male offspring"),
  h2_A       = c(t_mother$h2,   t_daughter$h2),
  h2_B       = c(t_father$h2,   t_son$h2),
  Z          = c(z_parent$Z,    z_offspring$Z),
  p_value    = c(z_parent$p,    z_offspring$p),
  significant = c(z_parent$p < 0.05, z_offspring$p < 0.05),
  stringsAsFactors = FALSE
)
print(z_result)

# ============================================================
# 9. Likelihood-ratio tests for sex-specific transmission
# Model A (base) vs B (+ offspring-sex interaction) vs C (+ parent-sex interaction) vs D (both interactions)
# ============================================================
cat("\n--- Likelihood-ratio tests ---\n")

glm_base     <- glm(tori_child_bin ~ tori_parent_bin + sex_child_bin + relation_bin + age_child,
                    data=df, family=binomial())
glm_int_sex  <- glm(tori_child_bin ~ tori_parent_bin * sex_child_bin + relation_bin + age_child,
                    data=df, family=binomial())
glm_int_rel  <- glm(tori_child_bin ~ tori_parent_bin * relation_bin + sex_child_bin + age_child,
                    data=df, family=binomial())
glm_int_both <- glm(tori_child_bin ~ tori_parent_bin * sex_child_bin +
                      tori_parent_bin * relation_bin + age_child,
                    data=df, family=binomial())

lrt_sex  <- anova(glm_base, glm_int_sex,  test="LRT")
lrt_rel  <- anova(glm_base, glm_int_rel,  test="LRT")
lrt_both <- anova(glm_base, glm_int_both, test="LRT")

p_lrt_sex  <- lrt_sex[2,  "Pr(>Chi)"]
p_lrt_rel  <- lrt_rel[2,  "Pr(>Chi)"]
p_lrt_both <- lrt_both[2, "Pr(>Chi)"]

cat(sprintf("LRT tori_parent × offspring_sex：p = %.4f\n", p_lrt_sex))
cat(sprintf("LRT tori_parent × parent_sex   ：p = %.4f\n", p_lrt_rel))
cat(sprintf("LRT both (joint test)                     : p = %.4f\n", p_lrt_both))

lrt_result <- data.frame(
  test        = c("tori_parent × offspring_sex",
                  "tori_parent × parent_sex (relation)",
                  "Both interactions (joint)"),
  LRT_p       = round(c(p_lrt_sex, p_lrt_rel, p_lrt_both), 4),
  significant = c(p_lrt_sex < 0.05, p_lrt_rel < 0.05, p_lrt_both < 0.05),
  stringsAsFactors = FALSE
)
print(lrt_result)

# ============================================================
# 10. Manuscript Table 3: tetrachoric as primary and Falconer as supplementary
# ============================================================
table3 <- merge(
  h2_tetra[,    c("group","n","r_tetra","h2","CI_low","CI_high")],
  h2_falconer[, c("group","K_parent","K_offspring","h2","CI_low","CI_high")],
  by = "group", suffixes = c("_tetra","_falconer"), all.x = TRUE
)
# Arrange rows in group order
group_order <- c("Overall","Mother","Father","Female offspring","Male offspring",
                 "Mother-Daughter","Mother-Son","Father-Daughter","Father-Son")
table3 <- table3[match(group_order, table3$group), ]
rownames(table3) <- NULL
cat("\n=== Manuscript Table 3 ===\n"); print(table3)

# ============================================================
# 11. Export CSV files
# ============================================================
write.csv(h2_falconer, "result_step4_falconer.csv",    row.names=FALSE, fileEncoding="UTF-8")
write.csv(h2_tetra,    "result_step4_tetrachoric.csv", row.names=FALSE, fileEncoding="UTF-8")
write.csv(z_result,    "result_step4_z_test.csv",      row.names=FALSE, fileEncoding="UTF-8")
write.csv(lrt_result,  "result_step4_lrt.csv",         row.names=FALSE, fileEncoding="UTF-8")
write.csv(table3,      "result_step4_table3.csv",      row.names=FALSE, fileEncoding="UTF-8")

cat("\n=== Step 4 completed ===\n")
cat("✓ result_step4_falconer.csv\n")
cat("✓ result_step4_tetrachoric.csv\n")
cat("✓ result_step4_z_test.csv\n")
cat("✓ result_step4_lrt.csv\n")
cat("✓ result_step4_table3.csv\n")
cat("Note: Tetrachoric correlation is the primary method and underlies the Z-tests and Table 3; Falconer is included as a supplementary comparison.\n")
