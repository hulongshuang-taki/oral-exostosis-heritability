# Prepare the Wave 2-priority data set used in the sensitivity analysis.
# For participants examined in both waves, the Wave 2 phenotype is retained;
# Wave 1 is used when no Wave 2 record is available.

df <- read.csv("df_clean.csv", fileEncoding = "UTF-8", stringsAsFactors = FALSE)
cat("=== Raw data ===\n")
cat("Total rows: ", nrow(df), "\n")

required_cols <- c("id_child", "id_parent", "survey_child", "survey_parent")
if (!all(required_cols %in% colnames(df))) {
  stop("Missing required columns: ", paste(setdiff(required_cols, colnames(df)), collapse=", "))
}

# ====== Stage 1: Remove cross-wave mismatches ======
cat("\n=== Stage 1: Filter by survey wave ===\n")
df_wave_matched <- df[df$survey_child == df$survey_parent, ]
n_removed <- nrow(df) - nrow(df_wave_matched)
cat(sprintf("Before filtering: %d; after filtering: %d; cross-wave records removed: %d (%.1f%%)\n",
            nrow(df), nrow(df_wave_matched), n_removed,
            n_removed / nrow(df) * 100))

df_mismatch <- df[df$survey_child != df$survey_parent, ]
write.csv(df_mismatch, "result_wave_mismatch_removed.csv", row.names = FALSE)

# ====== Stage 1.5: Diagnose duplicate counts ======
pair_key_wm    <- paste(df_wave_matched$id_child, df_wave_matched$id_parent, sep="_")
pair_counts_wm <- table(pair_key_wm)
cat("\n=== Duplicate distribution after Stage 1 ===\n")
print(table(pair_counts_wm))

dup_keys_wm <- names(pair_counts_wm[pair_counts_wm > 1])
cat("Number of relationships still duplicated: ", length(dup_keys_wm), "\n")

if (max(pair_counts_wm) > 2) {
  cat("Warning: Some relationships occur more than twice; please review them.\n")
}

# ====== Stage 2: Merge records using Wave 2 as the reference ======
cat("\n=== Stage 2: Merge same-wave duplicates (Wave 2 priority) ===\n")
cat("Strategy: use Wave 2 when available; otherwise use Wave 1; classify reversals (Wave 1 positive -> Wave 2 negative) as negative.\n")

tori_cols  <- grep("^tori", colnames(df), value = TRUE)
age_cols   <- grep("^age",  colnames(df), value = TRUE)
other_cols <- setdiff(colnames(df_wave_matched),
                      c(tori_cols, age_cols, "survey_child", "survey_parent"))

cat("Oral exostosis columns to merge: ", paste(tori_cols, collapse=", "), "\n")
cat("Age columns to merge: ",   paste(age_cols,  collapse=", "), "\n")

# Confirm the actual strings used for Wave 1 and Wave 2
cat("\nActual survey_child values (verify wave labels):\n")
print(table(df_wave_matched$survey_child, useNA="always"))
# Confirm these values using the output above; defaults are shown below:
WAVE1_LABEL <- "1次"
WAVE2_LABEL <- "2次"
POS_LABEL   <- "1. あり"
NEG_LABEL   <- "0. なし"

# Merge function that prioritizes Wave 2
merge_tori_wave2 <- function(vals, surveys,
                              wave2_label = WAVE2_LABEL,
                              pos_label   = POS_LABEL,
                              neg_label   = NEG_LABEL) {
  if (all(is.na(vals))) return(NA_character_)

  w2_idx <- which(surveys == wave2_label)

  # If Wave 2 is available, use its value directly
  if (length(w2_idx) > 0) {
    w2_val <- vals[w2_idx[1]]
    if (!is.na(w2_val)) return(w2_val)
  }

  # If Wave 2 is unavailable (Wave 1 only), use the Wave 1 value
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

# Count each type of record
n_w2_pos    <- 0  # Wave 2 positive (including incident cases)
n_w2_neg    <- 0  # Wave 2 negative
n_reversal  <- 0  # Reversal (Wave 1 positive -> Wave 2 negative)
n_incident  <- 0  # Incident case (Wave 1 negative -> Wave 2 positive)
n_w1_only   <- 0  # Wave 1 only

merged_list <- vector("list", length(unique(key_dup)))
names(merged_list) <- unique(key_dup)

for (k in unique(key_dup)) {
  rows       <- df_dup[key_dup == k, ]
  merged_row <- rows[1, ]

  # Oral exostosis columns: prioritize Wave 2
  for (col in tori_cols) {
    vals    <- rows[[col]]
    surveys <- rows$survey_child

    w1_idx <- which(surveys == WAVE1_LABEL)
    w2_idx <- which(surveys == WAVE2_LABEL)
    w1_val <- if (length(w1_idx) > 0) vals[w1_idx[1]] else NA
    w2_val <- if (length(w2_idx) > 0) vals[w2_idx[1]] else NA

    # Count using tori_child only to avoid double counting
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

  # Age columns: use the larger value (equivalent to the Wave 2 age)
  for (col in age_cols) {
    merged_row[[col]] <- safe_max_numeric(rows[[col]])
  }

  # Mark the survey column as merged
  merged_row$survey_child  <- "merged"
  merged_row$survey_parent <- "merged"

  # Check consistency across the remaining columns
  for (col in other_cols) {
    vals <- unique(rows[[col]])
    if (length(vals) > 1) {
      cat(sprintf("Warning: Relationship %s has inconsistent values in column %s: %s\n",
                  k, col, paste(vals, collapse=" vs ")))
    }
  }

  merged_list[[k]] <- merged_row
}

df_merged <- do.call(rbind, merged_list)
rownames(df_merged) <- NULL

# Reassemble the complete data set
df_final <- rbind(df_unique, df_merged)
rownames(df_final) <- NULL

cat("\n=== Final summary ===\n")
cat("Rows in the original df_clean.csv:              ", nrow(df), "\n")
cat("Rows after removing cross-wave mismatches:       ", nrow(df_wave_matched), "\n")
cat("Rows after final deduplication (analysis data):  ", nrow(df_final), "\n")
cat("\n--- Merge-strategy statistics for duplicate parent-offspring pairs ---\n")
cat(sprintf("  Wave 2 positive (including incident cases): %d\n", n_w2_pos))
cat(sprintf("  Wave 2 negative:                            %d\n", n_w2_neg))
cat(sprintf("  Reversals (W1 positive -> W2 negative):     %d (OR merging would overestimate these as positive)\n", n_reversal))
cat(sprintf("  Incident cases (W1 negative -> W2 positive): %d (both strategies classify these as positive)\n", n_incident))
cat(sprintf("  Wave 1 records only:                        %d\n", n_w1_only))

write.csv(df_final, "df_sensitivity_wave2.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\nExported: df_sensitivity_wave2.csv\n")
