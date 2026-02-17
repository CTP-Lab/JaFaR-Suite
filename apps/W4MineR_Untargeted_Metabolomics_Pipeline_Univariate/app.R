# =====================================================================
# Untargeted Metabolomics — Univariate v3.9.18 (KW switch) — ENGLISH UI
# Auto-converted from app_v3_9_9_kw_switch.txt with comprehensive FR→EN replacements.
# =====================================================================


# app_v3_9_9_kw_switch.R
# Untargeted Metabolomics — Univariate (W4M or MZmine CSV)
# This version adds a test_mode switch: "auto" (ANOVA/Welch/KW) or "kw_only" (force Kruskal–Wallis + Dunn)

suppressPackageStartupMessages({
  library(shiny); library(shinyWidgets); library(DT)
  library(readr); library(dplyr); library(tidyr); library(stringr); library(tibble)
  library(purrr); library(ggplot2); library(scales); library(ggrepel); library(cowplot)
  ggplot2::theme_set(ggplot2::theme_bw(base_size = 12))
  library(rstatix)  # games_howell_test, dunn_test
  library(car)      # leveneTest
})

# ---------- helpers ----------

# Flexible TSV/CSV reader
read_tabflex <- function(path) {
  stopifnot(length(path) == 1, !is.na(path))
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv","tab","tabular","txt")) {
    suppressMessages(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE, guess_max = 100000))
  } else {
    suppressMessages(readr::read_csv(path, show_col_types = FALSE, progress = FALSE, guess_max = 100000))
  }
}

# Detect a column by name among candidates (case-insensitive)
detect_col <- function(df, candidates) {
  nm <- names(df)
  for (x in candidates) {
    hit <- nm[tolower(nm) == tolower(x)]
    if (length(hit)) return(hit[[1]])
  }
  for (x in candidates) {
    hit <- nm[grepl(paste0("^", x, "$"), nm, ignore.case = TRUE)]
    if (length(hit)) return(hit[[1]])
  }
  NA_character_
}

# Robust canonicalization of sample IDs (handles paths, extensions, Unicode spaces, diacritics, punctuation)
canon_id <- function(x) {
  y <- as.character(x)
  y <- gsub("\\\\", "/", y, fixed = TRUE)     # normalize slashes
  y <- basename(y)                            # drop path
  y <- sub("\\.[A-Za-z0-9]+$", "", y)         # drop simple extension
  # normalize NBSP/NNBSP to regular space (literal non-breaking spaces)
  y <- gsub(" | ", " ", y, perl = TRUE)
  # transliterate diacritics when possible
  y2 <- tryCatch(iconv(y, to = "ASCII//TRANSLIT"), error = function(e) y)
  y[!is.na(y2)] <- y2[!is.na(y2)]
  y <- trimws(y)
  # replace whitespace, dots, dashes & punctuation by underscore
  y <- gsub("[[:space:].-]+", "_", y)
  y <- gsub("[^A-Za-z0-9_]+", "_", y)
  y <- gsub("_+", "_", y)
  y <- gsub("^_+|_+$", "", y)
  tolower(y)
}

# PQN (basic)
do_pqn <- function(X) {
  ref <- apply(X, 1, median, na.rm = TRUE)
  ratio <- sweep(X, 1, ref, "/")
  f <- apply(ratio, 2, median, na.rm = TRUE)
  X / rep(f, each = nrow(X))
}

# Robust QRILC extractor (handles different imputeLCMD versions)
safe_qrilc <- function(log2M, tune.sigma = 0.3) {
  if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
    stop("Package 'imputeLCMD' is required for QRILC.")
  }
  res <- imputeLCMD::impute.QRILC(as.matrix(log2M), tune.sigma = tune.sigma)
  if (is.matrix(res) || is.data.frame(res)) {
    imp <- as.matrix(res)
  } else if (is.list(res)) {
    if (!is.null(res$imputed) && (is.matrix(res$imputed) || is.data.frame(res$imputed))) {
      imp <- as.matrix(res$imputed)
    } else if (!is.null(res$imputedData) && (is.matrix(res$imputedData) || is.data.frame(res$imputedData))) {
      imp <- as.matrix(res$imputedData)
    } else if (length(res) >= 1 && (is.matrix(res[[1]]) || is.data.frame(res[[1]]))) {
      imp <- as.matrix(res[[1]])
    } else {
      stop(paste0("impute.QRILC returned unsupported object of class: ", paste(class(res), collapse = "/")))
    }
  } else {
    stop(paste0("impute.QRILC returned unsupported class: ", paste(class(res), collapse = "/")))
  }
  if (!is.null(rownames(log2M))) rownames(imp) <- rownames(log2M)
  if (!is.null(colnames(log2M))) colnames(imp) <- colnames(log2M)
  imp
}

# ---- AUTO imputation helpers ----
infer_mnar_flags <- function(log2X, min_missing = 0.05, max_missing = 0.95, trunc_score_thr = 3) {
  n <- nrow(log2X); flags <- rep(FALSE, n)
  for (i in seq_len(n)) {
    v <- as.numeric(log2X[i, ])
    mis <- is.na(v); m <- mean(mis)
    if (all(is.na(v)) || sum(!is.na(v)) < 6L) { flags[i] <- FALSE; next }
    w <- v[!is.na(v)]
    q10 <- tryCatch(as.numeric(stats::quantile(w, probs = 0.10, na.rm = TRUE)), error = function(e) NA_real_)
    q50 <- tryCatch(as.numeric(stats::quantile(w, probs = 0.50, na.rm = TRUE)), error = function(e) NA_real_)
    mn  <- suppressWarnings(min(w, na.rm = TRUE))
    if (!is.finite(q10) || !is.finite(q50) || !is.finite(mn)) { flags[i] <- FALSE; next }
    denom <- max(q10 - mn, .Machine$double.eps)
    trunc_score <- (q50 - q10) / denom
    flags[i] <- (m >= min_missing && m <= max_missing && trunc_score >= trunc_score_thr)
  }
  flags
}


choose_auto_method <- function(log2X,
                               pres_rate_mat = NULL,
                               onoff_high = 0.80, onoff_low = 0.20,
                               mean_mrate_cutoff = 0.02,
                               prop_onoff_cutoff = 0.05,
                               prop_mnar_cutoff = 0.20,
                               mnar_max_missing = 0.95,
                               mnar_min_missing = 0.05,
                               trunc_score_thr = 3) {
  mrate <- rowMeans(is.na(log2X))
  mean_mrate <- mean(mrate, na.rm = TRUE)
  if (!is.finite(mean_mrate)) mean_mrate <- 0

  # ON/OFF-like proportion (group-aware): features with high presence in ≥1 group and low presence in ≥1 other group
  prop_onoff <- NA_real_
  onoff_like <- NULL
  if (!is.null(pres_rate_mat)) {
    pr_max <- apply(pres_rate_mat, 1, max, na.rm = TRUE)
    pr_min <- apply(pres_rate_mat, 1, min, na.rm = TRUE)
    onoff_like <- is.finite(pr_max) & is.finite(pr_min) & (pr_max >= onoff_high) & (pr_min <= onoff_low)
    prop_onoff <- mean(onoff_like, na.rm = TRUE)
    if (!is.finite(prop_onoff)) prop_onoff <- 0
  }

  # Preset 1: LC–MS safe (v3.9.18)
  # Priority: (1) ON/OFF structure, (2) MNAR/left-censoring, then (3) allow missForest only for truly rare holes.
  if (is.finite(prop_onoff) && prop_onoff >= prop_onoff_cutoff) {
    return(list(
      method = "QRILC",
      reason = sprintf("Presence/absence patterns detected (ON/OFF-like=%.1f%% ≥ %.1f%%): using left-censored imputation (QRILC).",
                       100*prop_onoff, 100*prop_onoff_cutoff),
      mean_mrate = mean_mrate, prop_onoff = prop_onoff, prop_mnar = NA_real_,
      onoff_like = onoff_like
    ))
  }

  flags <- infer_mnar_flags(log2X, min_missing = mnar_min_missing, max_missing = mnar_max_missing, trunc_score_thr = trunc_score_thr)
  prop_mnar <- mean(flags, na.rm = TRUE)
  if (!is.finite(prop_mnar)) prop_mnar <- 0

  if (prop_mnar >= prop_mnar_cutoff) {
    return(list(
      method = "QRILC",
      reason = sprintf("Left-censoring likely for %.1f%% of features (≥ %.1f%%): using QRILC.",
                       100*prop_mnar, 100*prop_mnar_cutoff),
      flags = flags,
      mean_mrate = mean_mrate, prop_onoff = prop_onoff, prop_mnar = prop_mnar,
      onoff_like = onoff_like
    ))
  }

  if (mean_mrate < mean_mrate_cutoff) {
    return(list(
      method = "missForest",
      reason = sprintf("Very low missingness (mean per-feature missing=%.1f%% < %.1f%%) and no ON/OFF or MNAR signal: missForest is fine for rare technical holes.",
                       100*mean_mrate, 100*mean_mrate_cutoff),
      mean_mrate = mean_mrate, prop_onoff = prop_onoff, prop_mnar = prop_mnar,
      onoff_like = onoff_like, flags = flags
    ))
  }

  list(
    method = "missForest",
    reason = sprintf("MNAR-left-censoring not dominant (%.1f%% < %.1f%%) and ON/OFF-like below threshold: using missForest.",
                     100*prop_mnar, 100*prop_mnar_cutoff),
    flags = flags,
    mean_mrate = mean_mrate, prop_onoff = prop_onoff, prop_mnar = prop_mnar,
    onoff_like = onoff_like
  )
}


# ---------- pipeline ----------

# ---------- pipeline ----------
run_pipeline <- function(dataMatrix, sampleMetadata, variableMetadata,
                         outdir = "R_RESULTS",
                         alpha = 0.05, fc = 2.0,
                         up_only = FALSE, aggregator = "intersection", seed = 123,
                         center_param = "geometric", center_nonparam = "median",
                         do_auto_pqn = TRUE, pqn_cv_threshold = 0.20,
                         impute_method = "QRILC",
                         on_step = NULL, presence_thr = 0.80,
                         onoff_high = 0.80, onoff_low = 0.20,
                         test_mode = "auto") {   # <-- NEW

  set.seed(as.integer(seed))

# --- Sanitize ON/OFF thresholds ---
onoff_high <- max(0, min(1, as.numeric(onoff_high)))
onoff_low  <- max(0, min(1, as.numeric(onoff_low)))
if (is.na(onoff_high)) onoff_high <- 0.80
if (is.na(onoff_low))  onoff_low  <- 0.20
if (onoff_low > onoff_high) { tmp <- onoff_low; onoff_low <- onoff_high; onoff_high <- tmp }


  # === [1/14] Read & alignment ===
  if (!is.null(on_step)) on_step(1, "Read & alignment")

  dm <- read_tabflex(dataMatrix)
  sm <- read_tabflex(sampleMetadata)
  vm <- read_tabflex(variableMetadata)

  # Sanity checks & shapes
  if (nrow(dm) == 0) stop("dataMatrix has 0 rows.")
  if (ncol(dm) < 2) stop("dataMatrix must have >= 2 sample columns.")
  if (nrow(sm) < 2) stop("sampleMetadata must have >= 2 rows.")

  # Identify columns
  sid_col <- detect_col(sm, c("sample_name","sample","id","Sample","Sample_name"))
  cls_col <- detect_col(sm, c("class","Class","group","Group"))
  typ_col <- detect_col(sm, c("sampleType","sample_type","type","Type"))
  if (is.na(sid_col) || is.na(cls_col)) stop("sampleMetadata must contain sample_name and class columns")

  # Filter to biological samples only (strict; no fallback)
  if (!is.na(typ_col)) {
    type_vec <- tolower(sm[[typ_col]])
    bio_syn <- c("bio","biological","biologique","sample","samples","study","study sample","study_sample")
    sm <- sm[type_vec %in% bio_syn, , drop = FALSE]
    if (nrow(sm) < 2) stop("After filtering sampleType to biological samples (BIO), fewer than 2 samples remain.")
  }

  # dataMatrix row id
  rid_col <- detect_col(dm, c("name","feature_id","id","row_id","compound","feature"))
  if (is.na(rid_col)) stop("dataMatrix must contain a row id column (e.g., 'name')")
  row_ids <- dm[[rid_col]]
  dm <- dm %>% dplyr::select(-all_of(rid_col))
  dm <- as.data.frame(dm)
  rownames(dm) <- row_ids

  # Create output dir early (for debug files)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # Robust alignment of samples via canonical IDs
  samp_ids  <- sm[[sid_col]] %>% as.character()
  dm_cols   <- colnames(dm)
  samp_ids_n <- canon_id(samp_ids)
  dm_cols_n  <- canon_id(dm_cols)

  # Check duplicates after normalization
  dup_sm <- samp_ids[duplicated(samp_ids_n) | duplicated(samp_ids_n, fromLast = TRUE)]
  dup_dm <- dm_cols[duplicated(dm_cols_n)  | duplicated(dm_cols_n,  fromLast = TRUE)]
  if (length(dup_sm) || length(dup_dm)) {
    utils::write.table(unique(dup_sm), file.path(outdir, "DEBUG_duplicates_sampleMetadata.txt"),
                       row.names = FALSE, col.names = FALSE, quote = FALSE)
    utils::write.table(unique(dup_dm), file.path(outdir, "DEBUG_duplicates_dataMatrix.txt"),
                       row.names = FALSE, col.names = FALSE, quote = FALSE)
    stop("Duplicate IDs after normalization. See DEBUG_duplicates_* in output directory.")
  }

  # Map normalized SM id -> original SM id
  sm_map <- setNames(samp_ids, samp_ids_n)
  matched_sm <- sm_map[dm_cols_n]
  has_match  <- !is.na(matched_sm)

  # Debug preview
  preview <- data.frame(
    dataMatrix_col        = dm_cols,
    dataMatrix_col_norm   = dm_cols_n,
    sampleMetadata_id     = matched_sm,
    sampleMetadata_id_norm= names(sm_map)[match(matched_sm, sm_map)],
    has_match             = has_match,
    stringsAsFactors = FALSE
  )
  utils::write.table(preview, file.path(outdir, "DEBUG_matching_preview.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  common_idx <- which(has_match)
  if (length(common_idx) < 2) {
    only_dm <- dm_cols[!has_match]
    only_sm <- samp_ids[!samp_ids_n %in% dm_cols_n]
    utils::write.table(only_dm, file.path(outdir, "DEBUG_unmatched_dataMatrix_cols.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    utils::write.table(only_sm, file.path(outdir, "DEBUG_unmatched_sampleMetadata_ids.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    stop("No or too few matched sample IDs (need >=2). See DEBUG_matching_preview.tsv and DEBUG_unmatched_* in output directory.")
  }

  # Keep matched columns and order as in sampleMetadata
  dm <- dm[, common_idx, drop = FALSE]
  colnames(dm) <- as.character(matched_sm[common_idx])          # overwrite with matched SM ids
  sm <- sm[ match(colnames(dm), sm[[sid_col]]), , drop = FALSE] # reorder SM rows to DM columns

  # Align variableMetadata to dataMatrix rows if possible
  vm_rid_col <- detect_col(vm, c("name","feature_id","id","row_id","compound","feature"))
  if (!is.na(vm_rid_col)) {
    ord <- match(rownames(dm), vm[[vm_rid_col]])
    if (!any(is.na(ord))) {
      vm <- vm[ord, , drop = FALSE]
    } else {
      # fallback: leave vm as-is
    }
  } else if (nrow(vm) != nrow(dm)) {
    vm <- tibble::tibble(name = rownames(dm))
  }

  # === [2/14] Diagnostics des missings ===
  if (!is.null(on_step)) on_step(2, "Diagnostics des missings")

  X <- as.matrix(dm)
  storage.mode(X) <- "numeric"
  X[!is.finite(X) | X <= 0] <- NA_real_

  # Sanitize presence threshold from UI (0..1)
  presence_thr <- max(min(as.numeric(presence_thr), 1), 0)

  # Per-group presence filter (ANY group)
  n_features_before <- nrow(X)
  pres_mat <- !is.na(X)
  cls_vec <- as.character(sm[[cls_col]])
  classes_now <- unique(cls_vec)
  grp_indices <- lapply(classes_now, function(cl) which(cls_vec == cl))
  pres_rate_mat <- sapply(grp_indices, function(idx) {
    if (length(idx) == 0) return(rep(NA_real_, nrow(pres_mat)))
    rowMeans(pres_mat[, idx, drop = FALSE], na.rm = TRUE)
  })
  colnames(pres_rate_mat) <- classes_now
  keep_pg <- apply(pres_rate_mat >= presence_thr, 1, function(z) any(z, na.rm = TRUE))
  keep_pg[is.na(keep_pg)] <- FALSE
  if (!all(keep_pg)) {
    dropped <- rownames(X)[!keep_pg]
    if (length(dropped) > 0) {
      utils::write.table(dropped, file.path(outdir, "Filtered_out_perGroup_rule_features.txt"),
                         row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    X <- X[keep_pg, , drop = FALSE]
    pres_mat <- pres_mat[keep_pg, , drop = FALSE]
    if (exists("pres_rate_mat")) pres_rate_mat <- pres_rate_mat[keep_pg, , drop = FALSE]
  }
  n_features_after <- nrow(X)

  # === [3/14] Imputation ===
  if (!is.null(on_step)) on_step(3, "Imputation")

  log2X <- suppressWarnings(log2(X))
  log2X[!is.finite(log2X)] <- NA_real_
  impute_used <- NULL
  auto_reason <- NULL
  auto_counts <- NULL

  do_MinDet <- function(M, q = 0.01, delta = 0.1) {
    mins <- apply(M, 1, function(v) suppressWarnings(quantile(v, probs = q, na.rm = TRUE)))
    out <- M
    na_idx <- which(is.na(out), arr.ind = TRUE)
    if (nrow(na_idx) > 0) out[na_idx] <- mins[na_idx[,1]] - delta
    out
  }

  if (toupper(impute_method) == "AUTO") {
    auto <- choose_auto_method(log2X, pres_rate_mat = pres_rate_mat,
                              onoff_high = onoff_high, onoff_low = onoff_low,
                              mean_mrate_cutoff = 0.02, prop_onoff_cutoff = 0.05,
                              prop_mnar_cutoff = 0.20, mnar_max_missing = 0.95)
    chosen <- toupper(auto$method)
    auto_reason <- auto$reason
    mrate <- rowMeans(is.na(log2X))
    heavy <- which(is.finite(mrate) & mrate >= 0.80)
    if (chosen == "QRILC") {
      qr_ok <- FALSE
      try({
        imp_log2 <- safe_qrilc(log2X, tune.sigma = 0.3)
        qr_ok <- TRUE
      }, silent = TRUE)
      if (!qr_ok) imp_log2 <- do_MinDet(log2X, q = 0.01, delta = 0.1)
      if (length(heavy)) imp_log2[heavy, ] <- do_MinDet(log2X[heavy, , drop = FALSE], q = 0.01, delta = 0.1)
    } else if (chosen == "MISSFOREST") {
      if (requireNamespace("missForest", quietly = TRUE)) {
        imp_log2 <- missForest::missForest(as.matrix(log2X))$ximp
      } else if (requireNamespace("impute", quietly = TRUE)) {
        imp_log2 <- impute::impute.knn(as.matrix(log2X))$data
        auto_reason <- paste0(auto_reason, " (missForest unavailable → fallback to kNN)")
      } else {
        imp_log2 <- do_MinDet(log2X, q = 0.01, delta = 0.1)
        auto_reason <- paste0(auto_reason, " (missForest & kNN unavailable → fallback to MinDet)")
      }
      if (length(heavy)) imp_log2[heavy, ] <- do_MinDet(log2X[heavy, , drop = FALSE], q = 0.01, delta = 0.1)
    } else {
      stop(paste0("Unknown AUTO method: ", chosen))
    }
    impute_used <- paste0("AUTO→", chosen)
    if (!is.null(auto$flags)) {
      auto_counts <- list(
        features_total = nrow(log2X),
        mean_missing = if (!is.null(auto$mean_mrate)) auto$mean_mrate else mean(rowMeans(is.na(log2X)), na.rm = TRUE),
        prop_ONOFF_like = if (!is.null(auto$prop_onoff)) auto$prop_onoff else NA_real_,
        features_ONOFF_like = if (!is.null(auto$onoff_like)) sum(auto$onoff_like, na.rm = TRUE) else NA_integer_,
        prop_MNAR_like = if (!is.null(auto$prop_mnar)) auto$prop_mnar else NA_real_,
        features_MNAR_like = if (!is.null(auto$flags)) sum(auto$flags, na.rm = TRUE) else NA_integer_,
        features_heavy_missing = length(heavy)
      )
    }
  } else if (toupper(impute_method) == "QRILC") {
    qr_ok <- FALSE
    try({
      imp_log2 <- safe_qrilc(log2X, tune.sigma = 0.3)
      if (!all(dim(imp_log2) == dim(log2X))) stop("QRILC returned wrong dimensions")
      impute_used <- sprintf("QRILC (tune.sigma = %.1f)", 0.3)
      qr_ok <- TRUE
    }, silent = TRUE)
    if (!qr_ok) {
      impute_used <- sprintf("MinDet (q = %.2f)", 0.01)
      imp_log2 <- do_MinDet(log2X, q = 0.01, delta = 0.1)
    }
  } else if (toupper(impute_method) == "MINDET") {
    impute_used <- sprintf("MinDet (q = %.2f)", 0.01)
    imp_log2 <- do_MinDet(log2X, q = 0.01, delta = 0.1)
  } else if (toupper(impute_method) == "KNN") {
    if (!requireNamespace("impute", quietly = TRUE)) stop("Package 'impute' is required for kNN imputation")
    impute_used <- "kNN"
    imp_log2 <- impute::impute.knn(as.matrix(log2X))$data
  } else if (toupper(impute_method) == "MISSFOREST") {
    if (!requireNamespace("missForest", quietly = TRUE)) stop("Package 'missForest' is required for missForest imputation")
    impute_used <- "missForest"
    imp_log2 <- missForest::missForest(as.matrix(log2X))$ximp
  } else {
    stop("Unknown impute_method")
  }

  RAW <- 2^imp_log2
  LOG <- suppressWarnings(log10(RAW))

  # === [4/14] Info imputation ===
  if (!is.null(on_step)) on_step(4, "Info imputation")

  # === [5/14] Auto-PQN ===
  if (!is.null(on_step)) on_step(5, "PQN auto (option)")
  pqn_applied <- FALSE; cv_medians <- NA_real_
  med_by_sample <- apply(RAW, 2, median, na.rm = TRUE)
  cv_medians <- stats::sd(med_by_sample, na.rm = TRUE) / mean(med_by_sample, na.rm = TRUE)
  if (isTRUE(do_auto_pqn) && is.finite(cv_medians) && cv_medians >= pqn_cv_threshold) {
    RAW <- do_pqn(RAW)
    LOG <- suppressWarnings(log10(RAW))
    pqn_applied <- TRUE
  }

  # === [6/14] Groups ===
  if (!is.null(on_step)) on_step(6, "Groups")
  sm <- sm %>% mutate(`_sid` = .data[[sid_col]], `_cls` = .data[[cls_col]])
  classes <- unique(sm$`_cls`)
  if (length(classes) < 2) stop(sprintf("Need at least 2 classes (found: %d). Check sampleMetadata$class.", length(classes)))

  # === [7/14] Hypothèses ===
  if (!is.null(on_step)) on_step(7, if (identical(test_mode,"kw_only")) "Route forcée: Kruskal–Wallis (KW)" else "Vérif. hypothèses (Shapiro & Levene)")
  route <- rep(NA_character_, nrow(LOG)); names(route) <- rownames(LOG)

  if (identical(test_mode, "kw_only")) {
    route[] <- "KW"
  } else {
    for (i in seq_len(nrow(LOG))) {
      df <- tibble(value = as.numeric(LOG[i, ]), sample = colnames(LOG)) %>%
        left_join(sm %>% select(`_sid`,`_cls`) %>% distinct(), by = c("sample" = "_sid")) %>%
        mutate(group = factor(`_cls`))
      if (n_distinct(df$group) < 2) { route[i] <- "KW"; next }
      shaps <- df %>% group_by(group) %>% summarize(p = tryCatch(stats::shapiro.test(value)$p.value, error = function(e) NA_real_), .groups="drop")
      sh_ok <- all(shaps$p > alpha, na.rm = TRUE)
      lev_p <- tryCatch(car::leveneTest(value ~ group, center = median, data = df)[["Pr(>F)"]][1], error = function(e) NA_real_)
      if (isTRUE(sh_ok) && is.finite(lev_p) && lev_p > alpha) route[i] <- "ANOVA"
      else if (isTRUE(sh_ok)) route[i] <- "WELCH"
      else route[i] <- "KW"
    }
  }

  # === [8/14] Test global + FDR ===
  if (!is.null(on_step)) on_step(8, "Test global + FDR")
  mk_long <- function(mat_row) {
    vals <- as.numeric(mat_row)
    tibble(value = vals, sample = colnames(LOG)) %>%
      left_join(sm %>% select(`_sid`,`_cls`) %>% distinct(), by = c("sample" = "_sid")) %>%
      mutate(group = factor(`_cls`)) %>% select(value, group)
  }

  p_global <- rep(NA_real_, nrow(RAW))
  kw_epsilon2 <- rep(NA_real_, nrow(RAW))

  
if (identical(test_mode, "kw_only")) {
  for (i in seq_len(nrow(RAW))) {
    df <- mk_long(RAW[i, ]); if (n_distinct(df$group) < 2) next
    kt <- kruskal.test(value ~ group, data = df)
    p_global[i] <- as.numeric(kt$p.value)
    H <- as.numeric(kt$statistic)
    k <- nlevels(df$group)
    n <- nrow(df)
    if (is.finite(H) && is.finite(k) && is.finite(n) && n > k) {
      eps2 <- (H - k + 1) / (n - k)
      kw_epsilon2[i] <- max(0, min(1, as.numeric(eps2)))
    }
  }
} else {
    which_anova <- which(route == "ANOVA"); which_welch <- which(route == "WELCH"); which_kw <- which(route == "KW")
    if (length(which_anova)) {
      p_global[which_anova] <- sapply(which_anova, function(i) {
        df <- mk_long(LOG[i, ]); if (n_distinct(df$group) < 2) return(NA_real_)
        fit <- stats::aov(value ~ group, data = df)
        as.numeric(summary(fit)[[1]][["Pr(>F)"]][1])
      })
    }
    if (length(which_welch)) {
      p_global[which_welch] <- sapply(which_welch, function(i) {
        df <- mk_long(LOG[i, ]); if (n_distinct(df$group) < 2) return(NA_real_)
        as.numeric(stats::oneway.test(value ~ group, data = df, var.equal = FALSE)$p.value)
      })
    }
    
if (length(which_kw)) {
  for (i in which_kw) {
    df <- mk_long(RAW[i, ]); if (n_distinct(df$group) < 2) next
    kt <- kruskal.test(value ~ group, data = df)
    p_global[i] <- as.numeric(kt$p.value)
    H <- as.numeric(kt$statistic)
    k <- nlevels(df$group)
    n <- nrow(df)
    if (is.finite(H) && is.finite(k) && is.finite(n) && n > k) {
      eps2 <- (H - k + 1) / (n - k)
      kw_epsilon2[i] <- max(0, min(1, as.numeric(eps2)))
    }
  }
}}
  q_global <- p.adjust(p_global, method = "BH")

  # === [9/14] Post-hoc ===
  if (!is.null(on_step)) on_step(9, "Post-hoc")

  candidates <- which(is.finite(q_global) & q_global <= max(0.2, alpha))
  classes_vec <- sort(unique(sm$`_cls`))
  CLRS <- setNames(scales::hue_pal()(length(classes_vec)), classes_vec)
  pairs <- combn(classes_vec, 2, simplify = FALSE)
  post_q <- lapply(pairs, function(pr) rep(NA_real_, nrow(RAW)))
  names(post_q) <- sapply(pairs, function(pr) paste(pr, collapse = "_vs_"))

# Cliff's delta (pairwise effect size; robust ranking beyond q-values)
cliffs_delta <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 1 || length(y) < 1) return(NA_real_)
  mean(sign(outer(x, y, "-")))
}
idx_by_class <- lapply(classes_vec, function(cl) which(sm$`_cls` == cl))
names(idx_by_class) <- classes_vec
post_cliffs <- lapply(pairs, function(pr) rep(NA_real_, nrow(RAW)))
names(post_cliffs) <- sapply(pairs, function(pr) paste(pr, collapse = "_vs_"))


  for (ii in candidates) {
    dfA <- mk_long(if (identical(test_mode,"kw_only") || route[ii] == "KW") RAW[ii, ] else LOG[ii, ])
    if (n_distinct(dfA$group) < 2) next

# Pairwise Cliff's delta computed on RAW (useful even when values are detection-limited)
for (pi in seq_along(pairs)) {
  pr <- pairs[[pi]]
  xa <- RAW[ii, idx_by_class[[pr[1]]]]
  xb <- RAW[ii, idx_by_class[[pr[2]]]]
  post_cliffs[[pi]][ii] <- cliffs_delta(xa, xb)
}

    if (identical(test_mode, "kw_only")) {
      dn <- try(rstatix::dunn_test(dfA, value ~ group, p.adjust.method = "none"), silent = TRUE)
      if (!inherits(dn, "try-error") && is.data.frame(dn)) {
        pcol <- intersect(c("p","p.value","p.adj"), names(dn)); pcol <- if (length(pcol)) pcol[1] else NA_character_
        for (pi in seq_along(pairs)) {
          pr <- pairs[[pi]]
          row <- dn %>% dplyr::filter((group1 == pr[1] & group2 == pr[2]) | (group1 == pr[2] & group2 == pr[1])) %>% dplyr::slice(1)
          if (nrow(row) && !is.na(pcol)) {
            val <- suppressWarnings(as.numeric(row[[pcol]][1]))
            if (length(val) == 1 && is.finite(val)) post_q[[pi]][ii] <- val
          }
        }
      }
    
} else if (route[ii] == "ANOVA") {
  # Pairwise pooled-variance t-tests (unadjusted p-values). FDR control is applied later across features per contrast.
  pw <- try(stats::pairwise.t.test(dfA$value, dfA$group, p.adjust.method = "none", pool.sd = TRUE), silent = TRUE)
  if (!inherits(pw, "try-error") && is.list(pw) && !is.null(pw$p.value)) {
    pmat <- pw$p.value
    for (pi in seq_along(pairs)) {
      pr <- pairs[[pi]]
      a <- pr[1]; b <- pr[2]
      val <- NA_real_
      if (!is.null(rownames(pmat)) && !is.null(colnames(pmat))) {
        if (a %in% rownames(pmat) && b %in% colnames(pmat)) val <- pmat[a, b]
        if (!is.finite(val) && b %in% rownames(pmat) && a %in% colnames(pmat)) val <- pmat[b, a]
      }
      if (length(val) == 1 && is.finite(val)) post_q[[pi]][ii] <- as.numeric(val)
    }
  }
} else if (route[ii] == "WELCH") {
  # Pairwise Welch t-tests (unadjusted p-values). FDR control is applied later across features per contrast.
  pw <- try(stats::pairwise.t.test(dfA$value, dfA$group, p.adjust.method = "none", pool.sd = FALSE), silent = TRUE)
  if (!inherits(pw, "try-error") && is.list(pw) && !is.null(pw$p.value)) {
    pmat <- pw$p.value
    for (pi in seq_along(pairs)) {
      pr <- pairs[[pi]]
      a <- pr[1]; b <- pr[2]
      val <- NA_real_
      if (!is.null(rownames(pmat)) && !is.null(colnames(pmat))) {
        if (a %in% rownames(pmat) && b %in% colnames(pmat)) val <- pmat[a, b]
        if (!is.finite(val) && b %in% rownames(pmat) && a %in% colnames(pmat)) val <- pmat[b, a]
      }
      if (length(val) == 1 && is.finite(val)) post_q[[pi]][ii] <- as.numeric(val)
    }
  }

      gh <- try(rstatix::games_howell_test(dfA, value ~ group), silent = TRUE)
      if (!inherits(gh, "try-error") && is.data.frame(gh)) {
        pcol <- intersect(c("p.adj","p.value","p"), names(gh)); pcol <- if (length(pcol)) pcol[1] else NA_character_
        for (pi in seq_along(pairs)) {
          pr <- pairs[[pi]]
          row <- gh %>% dplyr::filter((group1 == pr[1] & group2 == pr[2]) | (group1 == pr[2] & group2 == pr[1])) %>% dplyr::slice(1)
          if (nrow(row) && !is.na(pcol)) {
            val <- suppressWarnings(as.numeric(row[[pcol]][1]))
            if (length(val) == 1 && is.finite(val)) post_q[[pi]][ii] <- val
          }
        }
      }
    } else { # KW in auto-route
      dn <- try(rstatix::dunn_test(dfA, value ~ group, p.adjust.method = "none"), silent = TRUE)
      if (!inherits(dn, "try-error") && is.data.frame(dn)) {
        pcol <- intersect(c("p","p.value","p.adj"), names(dn)); pcol <- if (length(pcol)) pcol[1] else NA_character_
        for (pi in seq_along(pairs)) {
          pr <- pairs[[pi]]
          row <- dn %>% dplyr::filter((group1 == pr[1] & group2 == pr[2]) | (group1 == pr[2] & group2 == pr[1])) %>% dplyr::slice(1)
          if (nrow(row) && !is.na(pcol)) {
            val <- suppressWarnings(as.numeric(row[[pcol]][1]))
            if (length(val) == 1 && is.finite(val)) post_q[[pi]][ii] <- val
          }
        }
      }
    }
  }
  # BH inter-features per contrast
  for (pi in seq_along(pairs)) {
    post_q[[pi]] <- p.adjust(post_q[[pi]], method = "BH")
  }

  # === [10/14] OVR log2FC & agrégation ===
  if (!is.null(on_step)) on_step(10, "OVR log2FC et agrégation")

  geo <- function(v) exp(mean(log(pmax(v, .Machine$double.eps))))

  classes_vec <- unique(sm$`_cls`)
  if (!exists("CLRS")) CLRS <- setNames(scales::hue_pal()(length(classes_vec)), classes_vec)
  LFC  <- matrix(NA_real_, nrow(RAW), length(classes_vec), dimnames = list(rownames(RAW), classes_vec))
  QOVR <- matrix(NA_real_, nrow(RAW), length(classes_vec), dimnames = list(rownames(RAW), classes_vec))
  for (cl in classes_vec) {
    idx_g <- which(sm$`_cls` == cl); idx_r <- which(sm$`_cls` != cl)
    if (length(idx_g) < 1 || length(idx_r) < 1) next
    for (i in seq_len(nrow(RAW))) {
      if (!identical(test_mode,"kw_only") && route[i] %in% c("ANOVA","WELCH")) {
        cen_g <- if (center_param == "arithmetic") mean(RAW[i, idx_g], na.rm = TRUE) else geo(RAW[i, idx_g])
        cen_r <- if (center_param == "arithmetic") mean(RAW[i, idx_r], na.rm = TRUE) else geo(RAW[i, idx_r])
      } else {
        cen_g <- if (center_nonparam == "geometric") geo(RAW[i, idx_g]) else if (center_nonparam == "arithmetic") mean(RAW[i, idx_g], na.rm = TRUE) else median(RAW[i, idx_g], na.rm = TRUE)
        cen_r <- if (center_nonparam == "geometric") geo(RAW[i, idx_r]) else if (center_nonparam == "arithmetic") mean(RAW[i, idx_r], na.rm = TRUE) else median(RAW[i, idx_r], na.rm = TRUE)
      }
      LFC[i, cl] <- log2(cen_g / cen_r)
    }
    # Aggregate post-hoc q across pairs for class vs REST
    cols <- which(sapply(pairs, function(pr) cl %in% pr))
    vv <- rep(NA_real_, nrow(RAW))
    if (length(cols)) {
      agg_fun <- if (tolower(aggregator) == "intersection") pmax else pmin
      qagg <- do.call(agg_fun, c(lapply(cols, function(ci) post_q[[ci]]), list(na.rm = TRUE)))
      vv <- qagg
    }
    QOVR[, cl] <- vv
  }

  # UNIV_hit by class
  lfc_thr <- log2(fc)
  vmenr <- vm
  vm_rid_col <- detect_col(vm, c("name","feature_id","id","row_id","compound","feature"))
  if (!is.na(vm_rid_col)) {
    ord <- match(rownames(RAW), vm[[vm_rid_col]])
    if (!any(is.na(ord))) vmenr <- vm[ord, , drop = FALSE]
  }
  if (nrow(vmenr) != nrow(RAW)) {
    vmenr <- tibble::tibble(name = rownames(RAW))
  }

  for (cl in classes_vec) {
    univ <- (q_global <= alpha) & (QOVR[, cl] <= alpha) & (is.finite(LFC[, cl]) & (if (isTRUE(up_only)) LFC[, cl] >= lfc_thr else abs(LFC[, cl]) >= lfc_thr))
    vmenr[[paste0("UNIV_hit_", cl)]] <- as.integer(ifelse(is.na(univ), 0L, univ))
    vmenr[[paste0("log2FC_OVR_final_", cl)]] <- LFC[, cl]
    vmenr[[paste0("q_posthoc_agg_OVR_", cl)]] <- QOVR[, cl]
  }
  vmenr$route <- if (identical(test_mode,"kw_only")) "KW" else route
  vmenr$q_global <- q_global

# Effect size (global KW): epsilon-squared (computed for KW features; NA otherwise)
vmenr$kw_epsilon2 <- kw_epsilon2

# Pairwise post-hoc q-values and effect sizes (Cliff's delta)
if (exists("post_q") && length(post_q)) {
  for (nm in names(post_q)) {
    vmenr[[paste0("q_pairwise_", nm)]] <- post_q[[nm]]
  }
}
if (exists("post_cliffs") && length(post_cliffs)) {
  for (nm in names(post_cliffs)) {
    vmenr[[paste0("cliffs_", nm)]] <- post_cliffs[[nm]]
  }
}


# --- Presence rates & missingness (helps interpret ON/OFF features) ---
if (exists("pres_rate_mat")) {
  for (cl in colnames(pres_rate_mat)) {
    vmenr[[paste0("presence_rate_", cl)]] <- pres_rate_mat[, cl]
  }
}
if (exists("pres_mat")) {
  vmenr$missing_rate <- 1 - rowMeans(pres_mat)
}
# Simple ON/OFF-like flag: at least one class >=80% present and at least one class <=20% present
if (exists("pres_rate_mat") && ncol(pres_rate_mat) >= 2) {
  pr_max <- apply(pres_rate_mat, 1, max, na.rm = TRUE)
  pr_min <- apply(pres_rate_mat, 1, min, na.rm = TRUE)
  vmenr$ONOFF_like <- as.integer(pr_max >= onoff_high & pr_min <= onoff_low)
  # Optional: which classes are max/min (useful in tri+ group designs)
  vmenr$ONOFF_max_class <- colnames(pres_rate_mat)[apply(pres_rate_mat, 1, which.max)]
  vmenr$ONOFF_min_class <- colnames(pres_rate_mat)[apply(pres_rate_mat, 1, which.min)]
} else {
  vmenr$ONOFF_like <- 0L
}
  # Overall
  hit_cols_any <- grep("^UNIV_hit_", names(vmenr), value = TRUE)
  if (length(hit_cols_any)) {
    hits_any <- as.matrix(vmenr[, hit_cols_any, drop = FALSE])
    storage.mode(hits_any) <- "numeric"
    vmenr$UNIV_hit <- as.integer(rowSums(hits_any, na.rm = TRUE) >= 1L)
  } else if (nrow(vmenr) > 0) { vmenr$UNIV_hit <- 0L } else { vmenr$UNIV_hit <- integer(0) }

  # === [11/14] Tables & Plots ===
  if (!is.null(on_step)) on_step(11, "Tables & Plots")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  readr::write_tsv(as_tibble(vmenr), file.path(outdir, "Post_Univariate_variableMetadata_enriched.tsv"))
  readr::write_tsv(as_tibble(tibble::rownames_to_column(as.data.frame(RAW), var = "name")), file.path(outdir, "Post_Univariate_dataMatrix_imputed.tsv"))
  readr::write_tsv(as_tibble(sm), file.path(outdir, "Post_Univariate_sampleMetadata_enriched.tsv"))

  # Barplots
  hit_cols <- grep("^UNIV_hit_", names(vmenr), value = TRUE)
  if (length(hit_cols)) {
    hits_mat <- as.matrix(vmenr[, hit_cols, drop = FALSE]); storage.mode(hits_mat) <- "numeric"
    classes_plot <- sub("^UNIV_hit_", "", hit_cols)

    sig_counts <- colSums(hits_mat, na.rm = TRUE)
    df_sig <- tibble::tibble(class = classes_plot, n = as.numeric(sig_counts)) %>% dplyr::arrange(dplyr::desc(n))
    readr::write_tsv(df_sig, file.path(outdir, "Counts_SignificantFeatures_per_class.tsv"))
    sub_lab <- sprintf("α = %s; FC ≥ %s; OVR aggregator = %s; up_only = %s",
                      signif(alpha, 3), format(fc, trim = TRUE), aggregator, if (isTRUE(up_only)) "TRUE" else "FALSE")
    
    if (nrow(df_sig) > 0) {
      p_sig <- ggplot2::ggplot(df_sig, ggplot2::aes(x = reorder(class, n), y = n, label = n, fill = class)) +
        ggplot2::geom_col() + ggplot2::geom_text(hjust = -0.15) + ggplot2::coord_flip() +
        ggplot2::expand_limits(y = max(df_sig$n, na.rm = TRUE) * 1.08) +
        ggplot2::scale_fill_manual(values = CLRS)
      p_sig <- p_sig + ggplot2::labs(title = "Significant features per class", subtitle = sub_lab, x = "Number of features") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5),
                       axis.title.y = ggplot2::element_blank())
      ggplot2::ggsave(filename = file.path(outdir, "Barplot_SignificantFeatures_per_class.png"),
                      plot = p_sig, width = 7, height = 4, dpi = 150)
    }

    row_sums <- rowSums(hits_mat, na.rm = TRUE)
    unique_idx <- which(row_sums == 1)
    if (length(unique_idx)) {
      which_class <- apply(hits_mat[unique_idx, , drop = FALSE], 1, function(v) classes_plot[which(v == 1)[1]])
      df_uniq <- as.data.frame(table(which_class), stringsAsFactors = FALSE); colnames(df_uniq) <- c("class", "n")
      df_uniq <- df_uniq %>% dplyr::arrange(dplyr::desc(n))
    } else df_uniq <- tibble::tibble(class = classes_plot, n = 0L)
    readr::write_tsv(df_uniq, file.path(outdir, "Counts_SignificantFeatures_only_one_class.tsv"))
    
    if (nrow(df_uniq) > 0) {
      p_uni <- ggplot2::ggplot(df_uniq, ggplot2::aes(x = reorder(class, n), y = n, label = n, fill = class)) +
        ggplot2::geom_col() + ggplot2::geom_text(hjust = -0.15) + ggplot2::coord_flip() +
        ggplot2::expand_limits(y = max(df_uniq$n, na.rm = TRUE) * 1.08) +
        ggplot2::scale_fill_manual(values = CLRS)
      p_uni <- p_uni + ggplot2::labs(title = "Significant features in exactly one class", subtitle = sub_lab, x = "Number of features") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5),
                       axis.title.y = ggplot2::element_blank())
      ggplot2::ggsave(filename = file.path(outdir, "Barplot_SignificantFeatures_only_one_class.png"),
                      plot = p_uni, width = 7, height = 4, dpi = 150)
    }
  }

  # === [12/14] README ===
  if (!is.null(on_step)) on_step(12, "README")
  readme <- c(
    paste0("Version: v3.9.18 (KW switch + effect sizes + ON/OFF thresholds)"),
    paste0("Univariate test mode: ", test_mode),
    paste0("Per-group presence >= ", sprintf("%.0f%%", 100*presence_thr), " (ANY class): removed ", (n_features_before - n_features_after), " of ", n_features_before, " features"),
    paste0("Imputation: ", if (!is.null(impute_used)) impute_used else impute_method),

if (!is.null(auto_reason)) paste0("AUTO decision: ", auto_reason) else NULL,
if (!is.null(auto_counts)) paste0("AUTO metrics: mean_missing=", sprintf("%.3f", auto_counts$mean_missing),
                                 "; prop_ONOFF_like=", ifelse(is.finite(auto_counts$prop_ONOFF_like), sprintf("%.3f", auto_counts$prop_ONOFF_like), "NA"),
                                 "; prop_MNAR_like=", ifelse(is.finite(auto_counts$prop_MNAR_like), sprintf("%.3f", auto_counts$prop_MNAR_like), "NA"),
                                 "; cutoffs: mean_mrate<0.02→missForest; ONOFF>=0.05→QRILC; MNAR>=0.20→QRILC; mnar_max_missing=0.95") else NULL,
    paste0("Auto-PQN applied: ", ifelse(isTRUE(pqn_applied), "YES", "NO"),
           " (CV(medians)=", sprintf('%.3f', cv_medians), ", thr=", pqn_cv_threshold, ")"),
    paste0("Alpha=", alpha, "; FC=", fc, "; |log2FC|>=", log2(fc)),
    paste0("ON/OFF thresholds: high=", sprintf("%.2f", onoff_high), "; low=", sprintf("%.2f", onoff_low), " (ONOFF_like computed from presence_rate_<class>)"),
    paste0("Effect sizes in variableMetadata: kw_epsilon2 (Kruskal–Wallis epsilon-squared); cliffs_<A_vs_B> (Cliff's delta per pairwise contrast); q_pairwise_<A_vs_B> (BH across features per contrast)"),
    "Post-hoc p-values (pre-FDR): KW uses Dunn (unadjusted p); ANOVA-route uses pairwise pooled-variance t-tests (unadjusted p); WELCH-route uses pairwise Welch t-tests (unadjusted p). Then, BH is applied across features separately for each contrast to yield q_pairwise_<A_vs_B>.",
    paste0("Aggregator OVR: ", aggregator),
    paste0("Univariate stats: global test depends on test_mode (kw_only forces Kruskal-Wallis; auto chooses ANOVA/Welch/KW per feature). q_global = BH-adjusted p-values across features"),
    paste0("Post-hoc: KW uses Dunn with raw p-values (p.adjust.method=none); ANOVA uses TukeyHSD adjusted p-values; Welch uses Games-Howell adjusted p-values. Then BH correction is applied across all features separately for each pairwise contrast to obtain post-hoc q-values"),
    paste0("Interpretation helpers in variableMetadata: presence_rate_<class>, missing_rate, ONOFF_like, ONOFF_max_class, ONOFF_min_class"),
    paste0("BIO filter: strict; sampleType in {bio, biological, sample, samples, study, study sample, study_sample}"),
    paste0("Matching: robust normalization (basename, extension drop, Unicode spaces, diacritics transliteration, punctuation to underscore)")
  )
  writeLines(readme, con = file.path(outdir, "README_results.txt"))

  # === [13/14] Final export ===
  if (!is.null(on_step)) on_step(13, "Final export")

  # === [14/14] Finished ===
  if (!is.null(on_step)) on_step(14, "Finished")

  invisible(list(outdir = outdir))
}

# --- Utilities for UI helpers ---
replace_unicode_spaces <- function(x) {
  for (ch in c("\u00A0","\u202F","\u2009","\u200A","\u200B","\u2060")) x <- gsub(ch, " ", x, fixed = TRUE)
  x
}

norm_names <- function(x) {
  x <- as.character(x)
  x <- replace_unicode_spaces(x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  x <- tolower(x)
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

strip_suffixes <- function(x) {
  y <- x
  y <- gsub("\\.(mzml|mzxml|raw)\\s*(peak\\s*(area|height))?$", "", y, ignore.case = TRUE)
  y <- gsub("\\.(mzml|mzxml|raw)$", "", y, ignore.case = TRUE)
  y <- gsub("\\s*(peak\\s*area|peak\\s*height)\\s*$", "", y, ignore.case = TRUE)
  y <- gsub("\\s*area\\s*$", "", y, ignore.case = TRUE)
  y <- gsub("\\s*height\\s*$", "", y, ignore.case = TRUE)
  y <- gsub("_+", "_", y)
  y <- gsub("^_|_$", "", y)
  y
}

canonicalize_sampleType <- function(x, class_vec = NULL) {
  z <- tolower(as.character(x))
  z <- dplyr::case_when(
    z %in% c("qc","pool","qcs","qualitycontrol","quality_control") ~ "QC",
    z %in% c("blank","blanc","bkc","proc","procedural_blank") ~ "BLANK",
    z %in% c("bio","biological","sample","samples") ~ "BIO",
    TRUE ~ NA_character_
  )
  if (!is.null(class_vec)) {
    z <- ifelse(is.na(z) & grepl("(?i)^(QC|POOL)$", class_vec), "QC", z)
    z <- ifelse(is.na(z) & grepl("(?i)^(BLANK|BKC|BLANC|PROC)$", class_vec), "BLANK", z)
    z[is.na(z)] <- "BIO"
  } else {
    z[is.na(z)] <- "BIO"
  }
  z
}

read_delim_auto <- function(path) {
  if (grepl("\\.(tsv|txt|tab|tabular)$", path, ignore.case = TRUE)) {
    readr::read_tsv(path, show_col_types = FALSE, guess_max = 200000)
  } else {
    readr::read_csv(path, show_col_types = FALSE, guess_max = 200000)
  }
}

detect_key_cols <- function(df) {
  nm <- colnames(df); nn <- norm_names(nm)
  pick_exact <- function(cands) {
    cands_n <- norm_names(cands)
    for (cn in cands_n) {
      hit <- which(nn == cn)
      if (length(hit) > 0) return(hit[1])
    }
    NA_integer_
  }
  id_idx <- pick_exact(c("row ID","# featureID","featureID","ID","row_ID","rowid","feature_id","name"))
  mz_idx <- pick_exact(c("row m/z","m/z","mz","row mz","mz_centroid","mzmed"))
  rt_idx <- pick_exact(c("row retention time (sec)","row retention time (min)","retention time","rt","RT [min]","rt_min","rt_sec","row_retention_time","rtmed"))
  list(id=id_idx, mz=mz_idx, rt=rt_idx)
}

canonicalize_sm_headers <- function(SM) {
  cn <- replace_unicode_spaces(colnames(SM))
  colnames(SM) <- trimws(cn)
  colnames(SM) <- make.unique(colnames(SM), sep = "_")
  nn <- norm_names(colnames(SM))
  cand_sample <- c("sample_name","sample","file","filename","file_name","raw_file","rawfile","name","sampleid","sample_id","samplemetadata")
  cand_class  <- c("class","group","classe","condition","biogroup","category","group_name","groupid","group_id","class_name")
  cand_type   <- c("sampletype","type","sample_type")
  pick <- function(cands) which(nn %in% norm_names(cands))
  iS <- pick(cand_sample); if (length(iS)) colnames(SM)[iS[1]] <- "sample_name"
  iC <- pick(cand_class);  if (length(iC)) colnames(SM)[iC[1]]  <- "class"
  iT <- pick(cand_type);   if (length(iT)) colnames(SM)[iT[1]]   <- "sampleType"
  SM
}

# ---- MZmine CSV builder ----
# --------- Preprocess W4M triplet (read + normalize headers + write temp TSV) ---------
preprocess_w4m_triplet <- function(dm_path, sm_path, vm_path) {
  stopifnot(file.exists(dm_path), file.exists(sm_path), file.exists(vm_path))
  DM <- read_delim_auto(dm_path)
  SM <- read_delim_auto(sm_path)
  VM <- read_delim_auto(vm_path)
  # Ensure first column is 'name' for DM and VM (W4M convention)
  if (ncol(DM) >= 1) colnames(DM)[1] <- "name"
  if (ncol(VM) >= 1) colnames(VM)[1] <- "name"
  # Canonicalize sampleMetadata headers (sample_name, class, sampleType, etc.)
  if (exists("canonicalize_sm_headers", mode = "function")) {
    SM <- canonicalize_sm_headers(SM)
  }
  # Write to a temp folder to have clean TSVs
  d <- file.path(tempdir(), paste0("w4m_", as.integer(Sys.time())))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dm_p <- file.path(d, "dataMatrix.tsv")
  sm_p <- file.path(d, "sampleMetadata.tsv")
  vm_p <- file.path(d, "variableMetadata.tsv")
  readr::write_tsv(as.data.frame(DM), dm_p)
  readr::write_tsv(as.data.frame(SM), sm_p)
  readr::write_tsv(as.data.frame(VM), vm_p)
  list(dm = dm_p, sm = sm_p, vm = vm_p, dir = d)
}
# ---------------- UI ----------------
ui <- fluidPage(
  titlePanel("Untargeted Metabolomics — Univariate"),
  sidebarLayout(
    sidebarPanel(
        fileInput("dm", "dataMatrix (TSV/CSV)", accept=c('.tsv','.csv','.txt','.tab','.tabular')),
        fileInput("sm", "sampleMetadata (TSV/CSV)", accept=c('.tsv','.csv','.txt','.tab','.tabular')),
        fileInput("vm", "variableMetadata (TSV/CSV)", accept=c('.tsv','.csv','.txt','.tab','.tabular')),
      textInput("outdir", "Output folder", value="R_RESULTS"),
      tags$hr(),
      # --- NEW: test mode switch ---
      selectInput("test_mode", "Univariate test mode",
                  choices = c("Auto (ANOVA/Welch/KW)" = "auto", "KW only" = "kw_only"),
                  selected = "auto"),
      tags$hr(),
      sliderInput("presence_thr", "Presence threshold (by class)", min=0.5, max=1.0, value=0.80, step=0.05),
      sliderInput("onoff_high", "ON/OFF high presence threshold", min=0, max=1, value=0.80, step=0.05),
      sliderInput("onoff_low",  "ON/OFF low presence threshold",  min=0, max=1, value=0.20, step=0.05),
      numericInput("alpha", "alpha (BH)", value=0.05, min=0.0001, step=0.01),
      numericInput("fc", "FC (|log2FC| >= log2(FC))", value=2.0, min=1.0, step=0.1),
      checkboxInput("up_only", "Up only (one-sided)", FALSE),
      selectInput("aggregator","OVR aggregator", choices=c("intersection","min"), selected="intersection"),
      selectInput("center_param","Parametric center", choices=c("geometric","mean","median"), selected="geometric"),
      selectInput("center_nonparam","Non-parametric center", choices = c("median","mean","geometric"), selected = "median"),
      checkboxInput("do_auto_pqn","Auto PQN if high median CV", TRUE),
      numericInput("pqn_cv_threshold","Median CV threshold for PQN", value=0.20, min=0, step=0.05),
      selectInput("impute_method","Imputation method", choices = c("QRILC","kNN","MinDet","missForest","AUTO"), selected = "AUTO"),
      numericInput("seed","seed", value=123, min=1, step=1),
      actionButton("run","Run analysis", class="btn-primary")
    ),
    mainPanel(
      verbatimTextOutput("log"),
      DTOutput("preview")
    )
  )
)

# ---------------- Server ----------------
server <- function(input, output, session) {
  logs <- reactiveVal(character(0))
  add_log <- function(x) logs(c(logs(), paste0(format(Sys.time(), "%H:%M:%S"), " - ", x)))
  append_log <- add_log
  output$log <- renderText(paste(logs(), collapse = "\n"))

  observeEvent(input$run, {
    req(input$outdir)
    outdir <- input$outdir
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    withProgress(message = "Pipeline", value = 0, {
      on_step <- function(step, msg) {
        incProgress(1/14, detail = paste0("[", step, "/14] ", msg))
        add_log(paste0("[", step, "/14] ", msg))
      }


        req(input$dm, input$sm, input$vm)
        tmpw <- preprocess_w4m_triplet(input$dm$datapath, input$sm$datapath, input$vm$datapath)
        dm_path <- tmpw$dm; sm_path <- tmpw$sm; vm_path <- tmpw$vm

      res <- tryCatch({
        run_pipeline(
          dataMatrix = dm_path,
          sampleMetadata = sm_path,
          variableMetadata = vm_path,
          outdir = outdir,
          alpha = input$alpha, fc = input$fc,
          up_only = isTRUE(input$up_only), aggregator = input$aggregator, seed = input$seed,
          center_param = input$center_param, center_nonparam = input$center_nonparam,
          do_auto_pqn = isTRUE(input$do_auto_pqn), pqn_cv_threshold = input$pqn_cv_threshold,
          impute_method = input$impute_method,
          on_step = on_step, presence_thr = input$presence_thr,
          onoff_high = input$onoff_high, onoff_low = input$onoff_low,
          test_mode = input$test_mode          # <-- pass new param
        )
      }, error = function(e) e)

      if (inherits(res, "error")) {
        add_log(paste0("ERROR: ", res$message))
        showNotification(res$message, type="error", duration=10)
        return(invisible(NULL))
      }

      add_log("Pipeline done.")
    })
  })
}

shinyApp(ui, server)
