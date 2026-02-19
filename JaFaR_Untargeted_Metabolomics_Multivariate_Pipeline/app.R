# =====================================================================
# JaFaR_Untargeted_Metabolomics_Multivariate_Pipeline
# =====================================================================


# ---- Utilities ----
safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

maybe_log2 <- function(DM, already_log = TRUE, offset = NULL) {
  if (already_log) return(DM)
  if (is.null(offset)) {
    min_pos <- suppressWarnings(min(DM[DM > 0], na.rm = TRUE))
    if (!is.finite(min_pos)) min_pos <- 1
    offset <- 0.01 * min_pos
  }
  log2(DM + offset)
}

check_classes <- function(SM, class_col, groups_needed = NULL) {
  stopifnot(class_col %in% colnames(SM))
  cl <- as.character(SM[[class_col]])
  if (!is.null(groups_needed)) {
    miss <- setdiff(groups_needed, unique(cl))
    if (length(miss)) stop("Missing groups in SM: ", paste(miss, collapse = ", "))
  }
  cl
}

# ---- Dunn provider (returns data.frame columns: Comparison, p) ----
.get_dunn_fun <- function() {
  if (requireNamespace("FSA", quietly = TRUE)) {
    return(function(x, g) {
      out <- FSA::dunnTest(x ~ g, method = "none")
      data.frame(Comparison = out$res$Comparison, p = out$res$P.unadj, stringsAsFactors = FALSE)
    })
  } else if (requireNamespace("DescTools", quietly = TRUE)) {
    return(function(x, g) {
      out <- DescTools::DunnTest(x ~ g, method = "none")
      tab <- as.data.frame(out$P)
      comps <- c(); ps <- c()
      rn <- rownames(tab); cn <- colnames(tab)
      for (i in seq_along(rn)) for (j in seq_along(cn)) {
        if (j <= i) next
        if (is.finite(tab[i,j])) {
          comps <- c(comps, paste(rn[i], "-", cn[j]))
          ps    <- c(ps,    tab[i,j])
        }
      }
      data.frame(Comparison = comps, p = ps, stringsAsFactors = FALSE)
    })
  } else if (requireNamespace("PMCMRplus", quietly = TRUE)) {
    return(function(x, g) {
      out <- PMCMRplus::kwAllPairsDunnTest(x ~ g, p.adjust.method = "none")
      tab <- as.data.frame(out$p.value)
      comps <- c(); ps <- c()
      rn <- rownames(tab); cn <- colnames(tab)
      for (i in seq_along(rn)) for (j in seq_along(cn)) {
        if (j <= i) next
        if (is.finite(tab[i,j])) {
          comps <- c(comps, paste(rn[i], "-", cn[j]))
          ps    <- c(ps,    tab[i,j])
        }
      }
      data.frame(Comparison = comps, p = ps, stringsAsFactors = FALSE)
    })
  } else {
    stop("No Dunn-capable package available. Install one of: FSA, DescTools, PMCMRplus.")
  }
}

# ---- Public API: Pairwise ----
compute_pairwise_wilcoxon <- function(DM, SM, class_col = "class",
                                      groupA, groupB,
                                      already_log = TRUE, offset = NULL) {
  cl <- check_classes(SM, class_col, groups_needed = c(groupA, groupB))
  DM <- maybe_log2(DM, already_log, offset)
  idxA <- which(cl == groupA); idxB <- which(cl == groupB)
  stopifnot(length(idxA) > 0, length(idxB) > 0)
  nfeat <- nrow(DM)
  pvals <- numeric(nfeat); l2fc <- numeric(nfeat)
  for (i in seq_len(nfeat)) {
    xa <- as.numeric(DM[i, idxA]); xa <- xa[is.finite(xa)]
    xb <- as.numeric(DM[i, idxB]); xb <- xb[is.finite(xb)]
    if (length(xa) >= 1 && length(xb) >= 1) {
      pvals[i] <- tryCatch(stats::wilcox.test(xa, xb, exact = FALSE, correct = FALSE)$p.value,
                           error = function(e) NA_real_)
      l2fc[i]  <- safe_median(xa) - safe_median(xb)
    } else { pvals[i] <- NA_real_; l2fc[i] <- NA_real_ }
  }
  qvals <- stats::p.adjust(pvals, method = "BH")
  data.frame(
    name = rownames(DM),
    log2FC = l2fc,
    p = pvals,
    q = qvals,
    stringsAsFactors = FALSE
  )
}

# ---- Public API: OVR ----
compute_OVR_KW_Dunn <- function(DM, SM, class_col = "class",
                                target_group,
                                already_log = TRUE, offset = NULL,
                                agg_mode = c("STRICT","UNION")) {
  agg_mode <- match.arg(agg_mode)
  cl <- check_classes(SM, class_col)
  DM <- maybe_log2(DM, already_log, offset)
  classes <- unique(cl)
  if (!(target_group %in% classes)) stop("target_group not in classes")
  others <- setdiff(classes, target_group)
  if (length(others) < 1) stop("OVR requires >=2 groups total")
  dunn_fun <- .get_dunn_fun()

  # pattern keys for C vs X regardless of order "C - X" or "X - C"
  comp_keys <- paste0("^(", target_group, " - (", paste(others, collapse="|"), ")|(", paste(others, collapse="|"), ") - ", target_group, ")$")

  nfeat <- nrow(DM)
  p_kw   <- numeric(nfeat)
  p_aggW <- numeric(nfeat)
  l2fc   <- numeric(nfeat)

  for (i in seq_len(nfeat)) {
    xi <- as.numeric(DM[i, ])
    gi <- factor(cl, levels = classes)
    # KW
    p_kw[i] <- tryCatch(stats::kruskal.test(xi ~ gi)$p.value, error = function(e) NA_real_)

    # Dunn C vs rest (3 p raw) -> Holm within-feature
    p_cr <- rep(NA_real_, length(others))
    if (is.finite(p_kw[i])) {
      dtab <- tryCatch(dunn_fun(xi, gi), error = function(e) NULL)
      if (!is.null(dtab) && nrow(dtab)) {
        dtab$Comparison <- gsub("\\s+", " ", trimws(dtab$Comparison))
        # extract each C vs X
        k <- 1
        for (oth in others) {
          pat <- paste0("^(", target_group, " - ", oth, "|", oth, " - ", target_group, ")$")
          hit <- grepl(pat, dtab$Comparison, perl = TRUE)
          if (any(hit)) p_cr[k] <- dtab$p[which(hit)[1]]
          k <- k + 1
        }
      }
    }
    p_within <- stats::p.adjust(p_cr, method = "BH")
    p_agg <- if (agg_mode == "STRICT") max(p_within, na.rm = TRUE) else min(p_within, na.rm = TRUE)
    if (!is.finite(p_agg)) p_agg <- NA_real_
    p_aggW[i] <- p_agg

    # Effect = diff of medians (log2)
    idxC <- which(cl == target_group); idxR <- which(cl != target_group)
    l2fc[i] <- safe_median(DM[i, idxC]) - safe_median(DM[i, idxR])
  }

  q_kw   <- stats::p.adjust(p_kw,   method = "BH")
  q_aggW <- stats::p.adjust(p_aggW, method = "BH")

  data.frame(
    feature_id = rownames(DM),
    p_global   = p_kw,   q_global = q_kw,
    p_posthoc_agg_OVR = p_aggW, q_posthoc_agg_OVR = q_aggW,
    log2FC_OVR_final  = l2fc,
    group = target_group, agg_mode = agg_mode,
    stringsAsFactors = FALSE
  )
}

# ---- Convenience: enrich VM with OVR columns for each class ----
enrich_VM_with_OVR <- function(VM, DM, SM, class_col = "class", classes = NULL,
                               already_log = TRUE, mode = c("STRICT","UNION")) {
  mode <- match.arg(mode)
  cl_all <- as.character(SM[[class_col]])
  if (is.null(classes)) classes <- unique(cl_all)
  for (cls in classes) {
    res <- stop("OVR stats recomputation disabled; using VM columns only")
    m <- match(rownames(DM), res$feature_id)
    # Always (re)store q_global (shared across modes)
    VM[["q_global"]] <- res$q_global[m]
    # Mode-specific columns
    suf <- if (mode == "UNION") "_UNION" else ""
    VM[[paste0("q_posthoc_agg_OVR_", cls, suf)]]  <- res$q_posthoc_agg_OVR[m]
    VM[[paste0("log2FC_OVR_final_", cls, suf)]]   <- res$log2FC_OVR_final[m]
  }
  VM
}

# ---- Override: pairwise_volcano_stats() used by the app ----
if (exists("pairwise_volcano_stats", mode = "function")) {
  rm(pairwise_volcano_stats)
}
pairwise_volcano_stats <- function(Xlog10, cls, g1, g2) {
  idx1 <- which(cls == g1); idx2 <- which(cls == g2)
  # effect = diff of medians on log10 -> convert to log2 by dividing log10 by log10(2)
  log2FC <- apply(Xlog10, 1, function(v) {
    med1 <- safe_median(v[idx1]); med2 <- safe_median(v[idx2])
    (med1 - med2) / log10(2)
  })
  pvals <- rep(NA_real_, nrow(Xlog10))
  for (i in seq_len(nrow(Xlog10))) {
    v1 <- Xlog10[i, idx1]; v2 <- Xlog10[i, idx2]
    if (sum(is.finite(v1)) >= 1 && sum(is.finite(v2)) >= 1) {
      pvals[i] <- tryCatch(stats::wilcox.test(v1, v2, exact = FALSE, correct = FALSE)$p.value, error = function(e) NA_real_)
    } else pvals[i] <- NA_real_
  }
  tibble::tibble(name = rownames(Xlog10), log2FC = log2FC, p = pvals, q = stats::p.adjust(pvals, method = "BH"))
}

# ======================= end stats_patch_OVR_pairwise_u10g.R =======================
# ===== END: Integrated stats patch =====

# Post-Run Diagnostics & Visuals - Shiny GUI v1.4.27f (EN) - S-plot fix

# --- v1.4.21 changes ---
# --- v1.4.22 changes ---
# --- v1.4.23 changes ---
# --- v1.4.24 changes ---
# - Added GUI control 'Number of cross-validation segments' passed as crossvalI to ropls::opls.
# - Runtime check clamps k <= number of samples (>=2) with a log message.

# - Heatmaps UNIV_hit: added fontsize controls; titles now explicitly state 'z-score' or 'log2FC'.

# - Use biosigner method name 'randomforest' (not 'rf').
# - Signature extraction accepts both 'A' and 'AS' labels.

# - Biosigner OVR: added per-method (RF/PLSDA/SVM) box/violin plots alongside the union plot.
# - Enriched dataMatrix: removed per-class mean_log10_* columns as requested.
# - Random Forest: method is already auto-enabled if package 'randomForest' is installed; logs clarify when missing.

suppressPackageStartupMessages({
  library(shiny)
# ==== DPI / PNG helpers (injected) ============================================
`%||%` <- function(a,b) if (!is.null(a) && !is.na(a)) a else b

.get_export_dpi <- function() {
  dpi <- getOption("PNG_EXPORT_DPI", 300L)
  dpi <- suppressWarnings(as.integer(dpi))
  if (!is.finite(dpi) || dpi < 300L) dpi <- 300L
  dpi
}

# Save a ggplot to exact pixel size at (>=)300 ppi without changing geometry
save_gg_fixedpx <- function(plot, file, width_px, height_px, dpi = .get_export_dpi(), bg = "white") {
  width_px  <- as.integer(width_px)
  height_px <- as.integer(height_px)
  dpi <- as.integer(dpi)
  if (!requireNamespace("ragg", quietly = TRUE)) {
    stop("Package 'ragg' is required for fixed-pixel saves.", call. = FALSE)
  }
  dev <- ragg::agg_png(filename = file, width = width_px, height = height_px, units = "px", res = dpi, background = bg)
  print(plot)
  grDevices::dev.off()
  invisible(file)
}

# One png() to rule them all (pixel-true, >=300 ppi), routing by filename
png <- function(filename, width, height, res = NULL, units = "px", ...) {
  dpi <- .get_export_dpi()

  bn <- basename(filename)
  sel <- function(k, w_fallback=width, h_fallback=height){
    w_def <- getOption("FIG_DEFAULT_W", w_fallback)
    h_def <- getOption("FIG_DEFAULT_H", h_fallback)
    c(getOption(paste0("FIG_", k, "_W"), w_def),
      getOption(paste0("FIG_", k, "_H"), h_def))
  }

  if (grepl("Splot",      bn, ignore.case=TRUE)) { wh <- sel("SPLOT",   2500,1600)
  } else if (grepl("Volcano",    bn, ignore.case=TRUE)) { wh <- sel("VOLCANO", 1400,1100)
  } else if (grepl("PCA_",       bn, ignore.case=TRUE)) { wh <- sel("PCA",     1400,1100)
  } else if (grepl("PLSDA|OPLSDA",bn, ignore.case=TRUE) && grepl("scores", bn, ignore.case=TRUE)) {
       wh <- sel("PLS",     1400,1100)
  } else if (grepl("permutation", bn, ignore.case=TRUE)) { wh <- sel("PERM",    1400,1100)
  } else if (grepl("confusion|cm",bn, ignore.case=TRUE)) { wh <- sel("CM",      1200, 900)
  } else if (grepl("pred|perf",   bn, ignore.case=TRUE)) { wh <- sel("PRED",    1400,1100)
  } else if (grepl("Venn",        bn, ignore.case=TRUE)) { wh <- sel("VENN",    1400,1100)
  } else if (grepl("Biosigner",    bn, ignore.case=TRUE)) { wh <- sel("BIOFEAT", 1400,1100)
  } else { wh <- sel("DEFAULT",   1400,1100) }

  width  <- as.integer(wh[1]); height <- as.integer(wh[2])

  if (requireNamespace("ragg", quietly=TRUE)) {
    return(ragg::agg_png(filename=filename, width=width, height=height, units="px", res=dpi, background="white"))
  } else {
    if (is.null(res) || is.na(res) || res < 300) res <- dpi
    return(grDevices::png(filename=filename, width=width, height=height, units="px", res=res, type="cairo-png"))
  }
}



# Shadow ggplot2::ggsave with pixel-true sizing based on filename
.orig_ggsave <- getFromNamespace("ggsave","ggplot2")
ggsave <- function(filename, plot = ggplot2::last_plot(), device = NULL, path = NULL,
                   scale = 1, width = NA, height = NA, units = c("in", "cm", "mm", "px"),
                   dpi = .get_export_dpi(), limitsize = TRUE, bg = NULL, ...) {
  units <- units[1]
  dpi <- .get_export_dpi()
  bn <- basename(filename)
  pick <- function(k, w_def, h_def) {
    w0 <- getOption("FIG_DEFAULT_W", w_def)
    h0 <- getOption("FIG_DEFAULT_H", h_def)
    c(getOption(paste0("FIG_", k, "_W"), w0),
      getOption(paste0("FIG_", k, "_H"), h0))
  }
  # Decide desired pixel geometry from filename if possible
  if (grepl("Volcano", bn, ignore.case=TRUE)) {
    wh <- pick("VOLCANO", 1400,1100)
  } else if (grepl("VIP", bn, ignore.case=TRUE)) {
    wh <- pick("VIP", 1400,1100)
  } else if (grepl("UNIV_hit", bn, ignore.case=TRUE) && grepl("heatmap|hm", bn, ignore.case=TRUE)) {
    wh <- pick("UNIVHM", 1400,1100)
  } else if (grepl("Biosigner", bn, ignore.case=TRUE)) {
    wh <- pick("BIOFEAT", 1400,1100)
  } else if (grepl("Venn", bn, ignore.case=TRUE)) {
    wh <- pick("VENN", 1400,1100)
  } else if (grepl("PCA_", bn, ignore.case=TRUE)) {
    wh <- pick("PCA", 1400,1100)
  } else if (grepl("PLSDA|OPLSDA", bn, ignore.case=TRUE) && grepl("scores", bn, ignore.case=TRUE)) {
    wh <- pick("PLS", 1400,1100)
  } else if (grepl("permutation", bn, ignore.case=TRUE)) {
    wh <- pick("PERM", 1400,1100)
  } else {
    wh <- pick("DEFAULT", 1400,1100)
  }
  # Force inches to match desired pixel geometry at our DPI
  width  <- as.numeric(wh[1]) / dpi
  height <- as.numeric(wh[2]) / dpi
  units <- "in"
  dpi <- dpi
  .orig_ggsave(filename = filename, plot = plot, device = device, path = path,
               scale = scale, width = width, height = height, units = units,
               dpi = dpi, limitsize = limitsize, bg = bg, ...)
}
# =============================================================================
# --- helper: save pheatmap object with pixel geometry driven by FIG_* options ---
save_pheatmap_px <- function(ph, file, key = "UNIVHM") {
  wh <- c(getOption(paste0("FIG_", key, "_W"), getOption("FIG_DEFAULT_W", 1400L)),
          getOption(paste0("FIG_", key, "_H"), getOption("FIG_DEFAULT_H", 1100L)))
  png(filename = file, width = as.integer(wh[1]), height = as.integer(wh[2]), res = .get_export_dpi())
  grid::grid.newpage(); grid::grid.draw(ph$gtable); grDevices::dev.off()
}



# ===== DPI & PNG WRAPPERS (injected) =====
.ragg_ok <- function(){ requireNamespace("ragg", quietly = TRUE) }
.get_export_dpi <- function(){ dpi <- getOption("PNG_EXPORT_DPI", 300); as.integer(max(300, suppressWarnings(as.numeric(dpi)))) }
.open_png_fixedpx <- function(file, width_px, height_px, dpi){
  width_px  <- as.integer(max(50, width_px))
  height_px <- as.integer(max(50, height_px))
  dpi <- as.integer(max(300, dpi))
  maxdim <- getOption("ragg.max_dim", 50000)
  options(ragg.max_dim = maxdim)
  width_px  <- min(width_px,  maxdim)
  height_px <- min(height_px, maxdim)
  if (.ragg_ok()) {
    ragg::agg_png(filename = file, width = width_px, height = height_px, units = "px", res = dpi, background = "white")
  } else {
    png(filename = file, width = width_px, height = height_px, units = "px", res = dpi, bg = "white", type = "cairo-png")
  }
}
.with_png_fixedpx <- function(file, width_px, height_px, dpi, code){
  .open_png_fixedpx(file, width_px, height_px, dpi)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}
save_gg_fixedpx <- function(plot, file, width_px, height_px, dpi = .get_export_dpi()){
  .with_png_fixedpx(file, width_px, height_px, dpi, print(plot))
}
save_grid_fixedpx <- function(grob, file, width_px, height_px, dpi = .get_export_dpi()){
  .with_png_fixedpx(file, width_px, height_px, dpi, grid::grid.draw(grob))
}
save_pheatmap_fixedpx <- function(pheat, file, width_px, height_px, dpi = .get_export_dpi()){
  if (!inherits(pheat, "pheatmap")) stop("save_pheatmap_fixedpx: object is not a pheatmap result")
  .with_png_fixedpx(file, width_px, height_px, dpi, grid::grid.draw(pheat$gtable))
}
# --- Monkey-patch ggsave & png to keep pixel geometry but write desired DPI metadata ---
.gg_to_px_dims <- function(width, height, units = "in", dpi = NA){
  units <- match.arg(as.character(units), c("px","in","cm","mm"))
  if (units == "px"){
    wpx <- as.integer(width); hpx <- as.integer(height)
  } else {
    if (is.na(dpi) || !is.finite(dpi)) dpi <- 96
    if (units == "cm"){ width <- width/2.54; height <- height/2.54 }
    if (units == "mm"){ width <- width/25.4; height <- height/25.4 }
    wpx <- as.integer(round(width  * dpi))
    hpx <- as.integer(round(height * dpi))
  }
  list(w=wpx, h=hpx)
}
# Shadow grDevices::png locally in the app environment
# =========================================

# ===== BIOSIGNER PNG/TSV COHERENCE (u10g) — single-file module =====
suppressPackageStartupMessages({
  library(ggplot2); library(readr); library(dplyr); library(tidyr); library(stringr); library(grid)
})
.norm_tier <- function(x){ x <- toupper(trimws(as.character(x))); x[x=="AS"]<-"S"; ok<-c("S","A","B","C","D","E"); x[!x%in%ok]<-NA_character_; x }

.parse_three <- function(s){
  # Prefer quoted tokens (they preserve duplicates and order). Fallback to bare tokens.
  q <- regmatches(s, gregexpr('\"([A-Za-z]+)\"', s, perl = TRUE))[[1]]
  if (length(q) >= 3) {
    return(gsub('^"|"$', '', q[1:3]))
  } else {
    bare <- unlist(strsplit(s, "[^A-Za-z]+", perl = TRUE))
    bare <- bare[nzchar(bare)]
    return(bare[1:3])
  }
}
parse_biosigner_model_u10g <- function(path){
  txt <- readLines(path, warn = FALSE, encoding = "UTF-8")
  i <- which(grepl("^\\s*plsda\\s+randomforest\\s+svm\\s*$", tolower(txt))); if (!length(i)) stop("Header not found in ", path); i <- i[1]
  block <- character(0)
  for (k in (i+1):length(txt)){
    z <- txt[k]; if (!nzchar(trimws(z))) break
    if (grepl("^\\s*accuracy\\s*:", tolower(z))) break
    block <- c(block, z)
  }
  block <- block[nzchar(trimws(block))]; if (!length(block)) stop("No data lines after header in ", path)
  rows <- lapply(block, function(z){
    m <- regexec("^\\s*(\\S+)\\s+(.*)$", z, perl=TRUE); g <- regmatches(z, m)[[1]]; if (length(g) < 3) return(NULL)
    fid <- g[2]; rest <- g[3]; tt <- .parse_three(rest)
    data.frame(feature_id=fid, plsda=tt[1], randomforest=tt[2], svm=tt[3], stringsAsFactors=FALSE)
  })
  df <- dplyr::bind_rows(rows)
  df <- df %>% dplyr::mutate(dplyr::across(c(plsda,randomforest,svm), .norm_tier))
  df <- df %>% dplyr::filter(!(is.na(plsda) & is.na(randomforest) & is.na(svm)))
  if (!nrow(df)) stop("All rows NA in ", path)
  ids <- df$feature_id
  ord <- if (all(grepl("^\\d+$", ids))) order(as.integer(ids), method="radix") else order(ids, method="radix")
  df[ord,,drop=FALSE]
}
plot_biosigner_heatmap_u10g <- function(df_long, out_png, title="Tiers of the selected features"){
  lv_tier <- c("S","A","B","C","D","E"); lv_model <- c("PLS-DA","Random Forest","SVM")
  dfm <- df_long %>%
    dplyr::rename(`PLS-DA`=plsda, `Random Forest`=randomforest, `SVM`=svm) %>%
    tidyr::pivot_longer(cols=c(`PLS-DA`,`Random Forest`,`SVM`), names_to="model", values_to="tier") %>%
    dplyr::mutate(feature_id = factor(feature_id, levels = unique(df_long$feature_id)),
           model = factor(model, levels = lv_model),
           tier = factor(.norm_tier(tier), levels = lv_tier))
  pal <- c(S="#1F9D55", A="#8BC34A", B="#DCE775", C="#FFE082", D="#FF8A65", E="#D84343")
  g <- ggplot2::ggplot(dfm, ggplot2::aes(x=model, y=feature_id, fill=tier)) +
    ggplot2::geom_tile(color="black", size=0.6, na.rm=FALSE) +
    ggplot2::scale_fill_manual(values=pal, drop=FALSE, na.value="grey90") +
    ggplot2::scale_y_discrete(breaks=levels(dfm$feature_id), labels=levels(dfm$feature_id), expand=ggplot2::expansion(mult=c(0,0))) +
    ggplot2::scale_x_discrete(expand=ggplot2::expansion(mult=c(0,0))) +
    ggplot2::labs(title=title, x=NULL, y=NULL, fill=NULL) +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold",hjust=0,size=20,margin=ggplot2::margin(b=10)),
          axis.text.x=ggplot2::element_text(angle=90,vjust=0.5,hjust=1,face="bold"),
          axis.text.y=ggplot2::element_text(face="bold"),
          panel.grid=ggplot2::element_blank(),
          legend.position="right", legend.box="vertical", legend.spacing.y=grid::unit(6,"pt"))
  ggsave(out_png, g, width=12, height=9, dpi=150, units="in")
}
biosigner_export <- function(path_model_txt, out_prefix){
  df <- parse_biosigner_model_u10g(path_model_txt)
  readr::write_tsv(df, paste0(out_prefix, "_model.tsv"))
  plot_biosigner_heatmap_u10g(df, paste0(out_prefix, "_plot.png"))
  invisible(list(tsv=paste0(out_prefix,"_model.tsv"), png=paste0(out_prefix,"_plot.png"), n=nrow(df)))
}
# ===== END BIOSIGNER u10g =====
; library(shinyWidgets); library(shinyFiles)


library(readr); library(dplyr); library(tidyr); library(tibble); library(stringr); library(purrr)
  library(ggplot2); library(ggrepel); library(cowplot); library(cowplot); library(patchwork)
  library(pheatmap); library(ggVennDiagram); library(ropls)

# --- OVR Volcano (VM-only) helper ---
`%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a
make_volcano_OVR_from_VM <- function(VM, class_name, lfc_thr, q_thr) {
  log2fc_col <- paste0("log2FC_OVR_final_", class_name)
  q_col      <- paste0("q_posthoc_agg_OVR_", class_name)
  if (!all(c(log2fc_col, q_col) %in% colnames(VM))) {
    stop("Missing required VM columns for OVR volcano: ", log2fc_col, " or ", q_col)
  }
  nm <- if ("name" %in% colnames(VM)) VM$name else if ("feature_id" %in% colnames(VM)) VM$feature_id else rownames(VM)
  df <- data.frame(
    name   = nm,
    log2FC = suppressWarnings(as.numeric(VM[[log2fc_col]])),
    q      = suppressWarnings(as.numeric(VM[[q_col]]))
  )
  df <- df[is.finite(df$log2FC) & is.finite(df$q), , drop = FALSE]
  df$neglog10q <- -log10(pmax(df$q, .Machine$double.xmin))
  df$sign_lfc <- abs(df$log2FC) >= lfc_thr
  df$sign_q   <- df$q <= q_thr
  df
}

; library(biosigner)
  library(viridisLite)

# ---- Helper: normalize row ID to numeric-only (e.g., "X64" -> "64") ----
normalize_row_id <- function(x) {
  sub("^\\D+", "", as.character(x))
}


})


# ==== BIOSIGNER SYNC U10 (regex-free, safe) ====
left_trim <- function(s) { while (nchar(s) > 0 && substr(s, 1, 1) == " ") s <- substr(s, 2, nchar(s)); s }
right_trim <- function(s) { while (nchar(s) > 0 && substr(s, nchar(s), nchar(s)) == " ") s <- substr(s, 1, nchar(s) - 1); s }
trim <- function(s) right_trim(left_trim(s))
collapse_spaces <- function(s) { if (!nzchar(s)) return(s); repeat { s2 <- gsub("  ", " ", s, fixed = TRUE); if (identical(s2, s)) break; s <- s2 }; s }
split_spaces <- function(s) { s1 <- collapse_spaces(trim(s)); toks <- character(0); cur <- ""; for (i in seq_len(nchar(s1))) { ch <- substr(s1, i, i); if (ch == " ") { if (nzchar(cur)) { toks <- c(toks, cur); cur <- "" } } else { cur <- paste0(cur, ch) } }; if (nzchar(cur)) toks <- c(toks, cur); toks }
is_digit_char <- function(ch) ch %in% as.character(0:9)
read_leading_integer <- function(s) { if (!nzchar(s)) return(NA_integer_); i <- 1; while (i <= nchar(s) && substr(s, i, i) == " ") i <- i + 1; if (i > nchar(s)) return(NA_integer_); j <- i; while (j <= nchar(s) && is_digit_char(substr(s, j, j))) j <- j + 1; if (j == i) return(NA_integer_); as.integer(substr(s, i, j - 1)) }

.bio_levels  <- function() c("S","A","B","C","D","E")
.bio_norm    <- function(x){ x <- toupper(trim(as.character(x))); x[!(x %in% .bio_levels())] <- NA_character_; factor(x, levels = .bio_levels(), ordered = TRUE) }
.bio_palette <- function() c("#1B9E77","#7CAE61","#D9EF8B","#FEE08B","#FC8D62","#D73027")

bio_parse_model_txt_safe <- function(path) {
  txt <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (!length(txt)) stop("File vide: ", path)
  has_sig <- which(vapply(txt, function(z) grepl("Significant features", z, fixed = TRUE), logical(1)))[1]
  if (is.na(has_sig)) stop("Pas de bloc 'Significant features' dans: ", path)
  norm_line <- function(z) collapse_spaces(tolower(trim(z)))
  norm <- vapply(txt, norm_line, character(1))
  hdr <- which(norm == "plsda randomforest svm"); hdr <- hdr[hdr > has_sig][1]
  if (is.na(hdr)) stop("Pas d'entête 'plsda randomforest svm' après le bloc, dans: ", path)
  acc <- which(vapply(txt, function(z) grepl("Accuracy", z, fixed = TRUE), logical(1)))
  acc <- acc[acc > hdr]; end <- if (length(acc)) acc[1] - 1L else length(txt)
  lines <- txt[seq.int(hdr + 1L, end)]
  keep <- vapply(lines, function(z) !is.na(read_leading_integer(z)), logical(1))
  lines <- lines[keep]
  if (!length(lines)) stop("Aucune ligne de features parsable dans: ", path)
  parse_row <- function(z) {
    z1 <- gsub('"', "", z,  fixed = TRUE)
    z1 <- gsub("'", "", z1, fixed = TRUE)
    toks <- split_spaces(z1)
    if (length(toks) < 4) return(NULL)
    if (is.na(suppressWarnings(as.integer(toks[1])))) return(NULL)
    data.frame(feature = toks[1],
               plsda = toupper(substr(toks[2], 1, 1)),
               randomforest = toupper(substr(toks[3], 1, 1)),
               svm = toupper(substr(toks[4], 1, 1)),
               stringsAsFactors = FALSE)
  }
  rows <- do.call(rbind, lapply(lines, parse_row))
  if (is.null(rows) || !nrow(rows)) stop("Impossible de parser les lignes de features: ", path)
  suppressWarnings(fnum <- as.numeric(rows$feature))
  ord <- if (all(!is.na(fnum))) order(fnum, rows$feature) else order(rows$feature)
  rows[ord, , drop = FALSE]
}

bio_build_matrix <- function(df){
  tab <- data.frame(feature = as.character(df$feature),
                    plsda = .bio_norm(df$plsda),
                    randomforest = .bio_norm(df$randomforest),
                    svm = .bio_norm(df$svm),
                    stringsAsFactors = FALSE)
  suppressWarnings(fnum <- as.numeric(tab$feature))
  ord <- if (all(!is.na(fnum))) order(fnum, tab$feature) else order(tab$feature)
  tab <- tab[ord, , drop = FALSE]; rownames(tab) <- tab$feature; tab$feature <- NULL
  as.matrix(tab)
}

bio_plot_matrix <- function(mat, out_png, title = "Tiers of the selected features"){
  lvl <- .bio_levels()
  z   <- apply(mat, 2, function(col) as.integer(factor(col, levels = lvl, ordered = TRUE)))
  pal <- .bio_palette()
  z_rev <- z[nrow(z):1, , drop = FALSE]; row_labels_rev <- rev(rownames(mat))
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE); par(mar = c(5, 6, 3, 6))
  png(out_png, width = 1100, height = 850, res = 120)
  image(x = seq_len(ncol(z_rev)), y = seq_len(nrow(z_rev)), z = t(z_rev),
        axes = FALSE, xlab = "", ylab = "", main = title,
        col = pal, useRaster = TRUE, zlim = c(1, length(lvl)))
  box(lwd = 3)
  axis(1, at = seq_len(ncol(z_rev)), labels = colnames(mat), las = 2, tick = FALSE, line = -0.5)
  axis(2, at = seq_len(nrow(z_rev)), labels = row_labels_rev, las = 2, tick = FALSE, line = -0.5)
  usr <- par("usr"); par(xpd = NA)
  legend(x = usr[2] + 0.5, y = mean(usr[3:4]), legend = lvl, fill = pal, bty = "n", y.intersp = 1.2, x.intersp = 0.8, cex = 1.1)
  dev.off()
}

biosigner_sync <- function(outdir) {
  files <- list.files(outdir, full.names = TRUE)
  if (!length(files)) return(invisible(NULL))
  b <- basename(files)
  ok <- startsWith(b, "Biosigner_OVR_") & endsWith(b, "_model.txt")
  files <- files[ok]
  if (!length(files)) return(invisible(NULL))
  for (path in files) {
    bn <- basename(path)
    scope <- substr(bn, nchar("Biosigner_OVR_") + 1L, nchar(bn) - nchar("_model.txt"))
    df  <- bio_parse_model_txt_safe(path)
    mat <- bio_build_matrix(df)
    out_tsv <- file.path(outdir, paste0("Biosigner_OVR_", scope, "_model.tsv"))
    write.table(data.frame(feature = rownames(mat),
                           plsda = as.character(mat[, "plsda"]),
                           randomforest = as.character(mat[, "randomforest"]),
                           svm = as.character(mat[, "svm"]), stringsAsFactors = FALSE),
                file = out_tsv, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    out_png <- file.path(outdir, paste0("Biosigner_OVR_", scope, "_plot.png"))
    bio_plot_matrix(mat, out_png, title = "Tiers of the selected features")
  }
  invisible(TRUE)
}
# ==== end helpers U10 ====


# Center titles globally for ggplot
theme_set(theme_bw(base_size = 12))
theme_update(plot.title = element_text(hjust = 0.5, face = "bold"))

# ---------- Helpers ----------
read_tabflex <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv","tab","tabular")) suppressMessages(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE))
  else suppressMessages(readr::read_delim(path, delim = ifelse(ext=="csv", ",", "\t"), show_col_types = FALSE, progress = FALSE, guess_max = 100000))
}

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

# ---- Build a single MZmine-style quant CSV from DM/SM/VM ----
make_mzmine_quant_enriched_csv <- function(dm_enr, sm_enr, vm_enr, out_csv, sid_col = NULL) {
  # Reuse detect_col() defined elsewhere in this file
  feat_vm <- detect_col(vm_enr, c("name","feature_id","id","row ID","row_ID","row_id","feature")); if (is.na(feat_vm)) feat_vm <- "name"
  feat_dm <- detect_col(dm_enr, c("name","feature_id","id","row ID","row_ID","row_id","feature")); if (is.na(feat_dm)) feat_dm <- "name"

  vm_enr[[feat_vm]] <- as.character(vm_enr[[feat_vm]])
  dm_enr[[feat_dm]] <- as.character(dm_enr[[feat_dm]])

  # m/z & RT (always output RT in minutes like MZmine; convert seconds -> minutes when needed)
  mz_col <- detect_col(vm_enr, c("row m/z","m/z","mz","mzmed","mz_mean","mean_mz","Mass"))
  rt_col <- detect_col(vm_enr, c("row retention time","retention time","retention_time","rt","RT","rtmed","rt_mean","retention_time_min"))

  # Header block (row ID, row m/z, row retention time)
  left <- data.frame(`row ID` = dm_enr[[feat_dm]], check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.na(mz_col)) {
    left[["row m/z"]] <- vm_enr[[mz_col]][match(left[["row ID"]], vm_enr[[feat_vm]])]
  }
  if (!is.na(rt_col)) {
    rt_vec <- vm_enr[[rt_col]][match(left[["row ID"]], vm_enr[[feat_vm]])]
    rt_num <- suppressWarnings(as.numeric(rt_vec))
    # Assume minutes only if column name explicitly mentions 'min'
    col_says_min <- grepl("min", rt_col, ignore.case = TRUE)
    if (!col_says_min) rt_num <- rt_num / 60
    left[["row retention time"]] <- round(rt_num, 8)
  }

  # Sample intensity columns
  sample_cols <- setdiff(colnames(dm_enr), feat_dm)
  if (is.null(sid_col)) sid_col <- detect_col(sm_enr, c("sample_name","Sample","Sample_name","sample","id","ID","filename","file"))
  if (!is.na(sid_col) && all(as.character(sm_enr[[sid_col]]) %in% sample_cols)) {
    sample_cols <- as.character(sm_enr[[sid_col]])
  }

  # Rename to "Sample Peak area" (MZmine style)
  rename_one <- function(sid) sprintf("%s Peak area", sid)
  new_names <- setNames(vapply(sample_cols, rename_one, character(1)), sample_cols)
  samp_df <- dm_enr[, c(feat_dm, sample_cols), drop = FALSE]
  colnames(samp_df)[match(sample_cols, colnames(samp_df))] <- unname(new_names)
  samp_df <- samp_df[, setdiff(colnames(samp_df), feat_dm), drop = FALSE]

  # Extra variables (all VM enriched columns excluding keys/mz/rt)
  vm_extra <- vm_enr[match(left[["row ID"]], vm_enr[[feat_vm]]), , drop = FALSE]
  vm_extra <- vm_extra[, setdiff(colnames(vm_extra), unique(c(feat_vm, mz_col, rt_col))), drop = FALSE]

  # Bind and write
  out_df <- cbind(left, samp_df, vm_extra, stringsAsFactors = FALSE)

  # Ensure row ID is numeric only (e.g., 'X64' -> 64)
  if ("row ID" %in% names(out_df)) {
    rid <- as.character(out_df[["row ID"]])
    rid_num <- suppressWarnings(as.integer(gsub("[^0-9]+", "", rid)))
    out_df[["row ID"]] <- rid_num
  }

  # Enforce filename suffix once: _post_univariate_and_multivariate.csv
  if (!grepl("(?i)_post_univariate_and_multivariate\\.csv$", out_csv, perl = TRUE)) {
    out_csv2 <- sub("(?i)\\.csv$", "_post_univariate_and_multivariate.csv", out_csv, perl = TRUE)
  } else {
    out_csv2 <- out_csv
  }

  readr::write_csv(out_df, out_csv2, na = "")
}



as_matrix_numeric <- function(df) { M <- as.matrix(df); storage.mode(M) <- "numeric"; M }

# --- Robust join helper to make keys character ---
ensure_char_key <- function(df, key = "name", source_col = NULL, .before = 1) {
  if (!is.data.frame(df)) stop("ensure_char_key: df must be a data.frame")
  if (!(key %in% names(df))) {
    if (is.null(source_col)) stop("ensure_char_key: key is missing and source_col is NULL")
    if (!(source_col %in% names(df))) stop(sprintf("ensure_char_key: source_col '%s' not found", source_col))
    df <- tibble::add_column(df, !!key := as.character(df[[source_col]]), .before = .before)
  } else {
    df[[key]] <- as.character(df[[key]])
  }
  df
}



pareto_scale <- function(M) { mu <- apply(M, 1, mean, na.rm = TRUE); sdv <- apply(M, 1, sd, na.rm = TRUE); sdv[!is.finite(sdv) | sdv == 0] <- 1; sweep(sweep(M, 1, mu, "-"), 1, sqrt(sdv), "/") }
unitvar_scale <- function(M) { mu <- apply(M, 1, mean, na.rm = TRUE); sdv <- apply(M, 1, sd, na.rm = TRUE); sdv[!is.finite(sdv) | sdv == 0] <- 1; sweep(sweep(M, 1, mu, "-"), 1, sdv, "/") }
ensure_dir <- function(path) { dir.create(path, showWarnings = FALSE, recursive = TRUE); path }
safe_write_gg <- function(plot, file, width=7, height=5, dpi=150) {
  bn <- basename(file)
  pick <- function(k, w_def, h_def) {
    w0 <- getOption("FIG_DEFAULT_W", w_def)
    h0 <- getOption("FIG_DEFAULT_H", h_def)
    c(getOption(paste0("FIG_", k, "_W"), w0),
      getOption(paste0("FIG_", k, "_H"), h0))
  }
  if (grepl("Volcano", bn, ignore.case=TRUE)) {
    wh <- pick("VOLCANO", 1400,1100)
  } else if (grepl("VIP", bn, ignore.case=TRUE)) {
    wh <- pick("VIP", 1400,1100)
  } else if (grepl("UNIV_hit", bn, ignore.case=TRUE) && grepl("heatmap|hm", bn, ignore.case=TRUE)) {
    wh <- pick("UNIVHM", 1400,1100)
  } else if (grepl("Biosigner", bn, ignore.case=TRUE)) {
    wh <- pick("BIOFEAT", 1400,1100)
  } else if (grepl("Venn", bn, ignore.case=TRUE)) {
    wh <- pick("VENN", 1400,1100)
  } else {
    wh <- pick("DEFAULT", 1400,1100)
  }
  save_gg_fixedpx(plot, file, width_px=as.integer(wh[1]), height_px=as.integer(wh[2]), dpi=.get_export_dpi())
}
stubify <- function(x) { x <- gsub("[^A-Za-z0-9]+", "_", x); x <- gsub("_+", "_", x); x <- gsub("^_|_$", "", x); x }

# Add a centered title to base graphics (ropls) plots
add_top_title <- function(txt, line = 0.7) { mtext(txt, side = 3, line = line, cex = 1.2, font = 2) }

# Diverging palette factory
get_div_palette <- function(name = "Blue-White-Red", n = 100) {
  switch(name,
    "RdBu" = colorRampPalette(c("#2166AC","white","#B2182B"))(n),
    "RdYlBu" = colorRampPalette(c("#2c7bb6","#ffffbf","#d7191c"))(n),
    "PRGn" = colorRampPalette(c("#276419","white","#762a83"))(n),
    "PuOr" = colorRampPalette(c("#2b8cbe","white","#d95f0e"))(n),
    "Blue-White-Red" = colorRampPalette(c("#2b8cbe","white","#d7301f"))(n),
    "Viridis (sequential)" = viridisLite::viridis(n),
    "Plasma (sequential)" = viridisLite::plasma(n),
    colorRampPalette(c("#2b8cbe","white","#d7301f"))(n)
  )
}

# Per-class log2FC (class mean vs rest mean) from raw intensities
log2fc_by_class <- function(Xraw, classes) {
  Xlog2 <- log(pmax(Xraw, .Machine$double.eps), base = 2)
  lv <- levels(classes)
  out <- matrix(NA_real_, nrow = nrow(Xraw), ncol = length(lv),
                dimnames = list(rownames(Xraw), lv))
  for (i in seq_along(lv)) {
    cl <- lv[i]
    idx1 <- which(classes == cl); idx2 <- which(classes != cl)
    out[, i] <- rowMeans(Xlog2[, idx1, drop = FALSE], na.rm = TRUE) -
                rowMeans(Xlog2[, idx2, drop = FALSE], na.rm = TRUE)
  }
  out
}

# Volcano utility
make_volcano <- function(df, lfc_col, p_col, title, lfc_thr=1, q_thr=0.05, top_label=10, pal = "Blue-White-Red") {
  dfx <- df %>% dplyr::mutate(
    sig = dplyr::case_when(
      is.finite(.data[[p_col]]) & .data[[p_col]] <= q_thr & is.finite(.data[[lfc_col]]) & .data[[lfc_col]] >= lfc_thr ~ "Up",
      is.finite(.data[[p_col]]) & .data[[p_col]] <= q_thr & is.finite(.data[[lfc_col]]) & .data[[lfc_col]] <= -lfc_thr ~ "Down",
      TRUE ~ "NS"
    )
  )
  ord <- order(dfx[[p_col]], -abs(dfx[[lfc_col]]), na.last = NA)
  lab_ids <- head(ord, top_label)
  dfx$label <- ifelse(seq_len(nrow(dfx)) %in% lab_ids, dfx$name, NA_character_)
  palv <- get_div_palette(pal, 3)
  
p <- ggplot(dfx, aes(x = .data[[lfc_col]], y = -log10(.data[[p_col]]), color = sig)) +
    geom_point(alpha = 0.8, size = 1.6, shape = 16) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
    geom_hline(yintercept = -log10(q_thr), linetype = "dashed") +
    scale_color_manual(values = c("Down" = palv[1], "Up" = palv[3], "NS" = "grey70")) +
    labs(x = "log2FC", y = "-log10(q)", title = title, color = NULL) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.title = element_blank())
p <- p + guides(color = ggplot2::guide_legend(override.aes = list(shape = 16)))
if (top_label > 0) {
    p <- p + ggrepel::geom_text_repel(aes(label = label, show.legend = FALSE, inherit.aes = FALSE, color = "black"),
           min.segment.length = 0, max.overlaps = Inf, size = 3, na.rm = TRUE)
  }
  p
}

pairwise_volcano_stats <- function(Xlog10, cls, g1, g2) {
  idx1 <- which(cls == g1); idx2 <- which(cls == g2)
  lfc <- apply(Xlog10, 1, function(v) mean(v[idx1], na.rm = TRUE) - mean(v[idx2], na.rm = TRUE))
  pvals <- rep(NA_real_, nrow(Xlog10))
  for (i in seq_len(nrow(Xlog10))) {
    v1 <- Xlog10[i, idx1]; v2 <- Xlog10[i, idx2]
    if (sum(is.finite(v1)) >= 2 && sum(is.finite(v2)) >= 2) {
      pv <- tryCatch(stats::t.test(v1, v2, var.equal = FALSE)$p.value, error = function(e) NA_real_)
      if (!is.finite(pv) || is.na(pv)) pv <- tryCatch(stats::wilcox.test(v1, v2)$p.value, error = function(e) NA_real_)
      pvals[i] <- pv
    } else pvals[i] <- NA_real_
  }
  tibble::tibble(name = rownames(Xlog10), log2FC = lfc / log10(2), p = pvals)
}

export_scores_loadings <- function(model, prefix, outdir) {
  sn <- slotNames(model)
  for (sl in sn[grepl("ScoreMN$", sn)]) {
    mat <- tryCatch(slot(model, sl), error = function(e) NULL)
    if (is.matrix(mat) && nrow(mat) > 0 && ncol(mat) > 0) {
      df <- as.data.frame(mat); df <- tibble::rownames_to_column(df, var = "sample")
      readr::write_tsv(df, file.path(outdir, paste0(prefix, "_", sl, ".tsv")))
    }
  }
  for (sl in sn[grepl("LoadingMN$", sn)]) {
    mat <- tryCatch(slot(model, sl), error = function(e) NULL)
    if (is.matrix(mat) && nrow(mat) > 0 && ncol(mat) > 0) {
      df <- as.data.frame(mat); df <- tibble::rownames_to_column(df, var = "feature")
      readr::write_tsv(df, file.path(outdir, paste0(prefix, "_", sl, ".tsv")))
    }
  }
  if ("vipVn" %in% sn) {
    vip <- tryCatch(slot(model, "vipVn"), error = function(e) NULL)
    if (is.numeric(vip) && length(vip) > 0) {
      df <- tibble::tibble(feature = names(vip), VIP = as.numeric(vip))
      readr::write_tsv(df, file.path(outdir, paste0(prefix, "_VIP.tsv")))
    }
  }
}



export_predictions <- function(model, prefix, outdir, X, y_true, append_log = function(...) {}) {
  # Prefer cross-validated predictions if available
  yPred <- tryCatch(slot(model, "yPredMN"), error = function(e) NULL)
  used_cv <- TRUE
  if (is.null(yPred) || (is.matrix(yPred) && nrow(yPred) == 0)) {
    used_cv <- FALSE
    append_log(paste0(prefix, ": yPredMN not found - using predict() (not CV)"))
    yPred <- tryCatch(predict(model, X), error = function(e) NULL)
  }
  if (is.null(yPred)) { append_log(paste0(prefix, ": predictions unavailable")); return(invisible(NULL)) }

  # Turn predictions into a data.frame and extract predicted class
  if (is.matrix(yPred)) {
    if (!is.numeric(yPred)) suppressWarnings(storage.mode(yPred) <- "numeric")
    idx <- max.col(yPred, ties.method = "first")
    pred_class <- colnames(yPred)[idx]
    all_na <- apply(yPred, 1, function(r) all(!is.finite(r)))
    pred_class[all_na] <- NA_character_
    df <- as.data.frame(yPred, stringsAsFactors = FALSE)
  } else {
    pred_class <- as.character(yPred)
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE);
df <- data.frame(pred_score = pred_class, stringsAsFactors = FALSE)
  }

  # Ensure rownames of df reflect sample names from X (avoid default '1..n')
  rn <- rownames(df)
  rn_is_seq <- is.null(rn) || identical(rn, as.character(seq_len(nrow(df)))) || identical(rn, seq_len(nrow(df)))
  if (rn_is_seq && !is.null(rownames(X))) rownames(df) <- rownames(X)[seq_len(nrow(df))]

  out <- tibble::rownames_to_column(df, var = "sample")

  # Align true labels by sample name
  names_y <- names(y_true)
  if (is.null(names_y) && !is.null(rownames(X))) names_y <- rownames(X)
  out$true <- as.character(y_true)[match(out$sample, names_y)]
  out$predicted <- pred_class
  readr::write_tsv(out, file.path(outdir, paste0(prefix, "_Predictions.tsv")))

  # Confusion matrix on complete cases
  keep <- stats::complete.cases(out$true, out$predicted)
  if (!any(keep)) { 
    append_log(paste0(prefix, ": no complete cases for confusion matrix"))
    return(invisible(out))
  }
  y_levels <- unique(as.character(na.omit(y_true)))
  tab <- table(factor(out$true[keep], levels = y_levels), factor(out$predicted[keep], levels = y_levels))
  cm_df <- as.data.frame.matrix(tab)
  cm_df <- tibble::rownames_to_column(cm_df, var = "True")
  acc <- if (sum(tab) > 0) sum(diag(tab)) / sum(tab) else NA_real_
  cm_file <- file.path(outdir, paste0(prefix, if (used_cv) "_ConfusionMatrix_CV.tsv" else "_ConfusionMatrix.tsv"))
  readr::write_tsv(cm_df, cm_file)
  metrics <- tibble::tibble(metric = c("accuracy","n","cv"), value = c(acc, sum(tab), as.numeric(used_cv)))
  readr::write_tsv(metrics, file.path(outdir, paste0(prefix, "_Metrics.tsv")))
  invisible(out)
}







splot_export <- function(
  model,
  Xmv,                 # scaled/centered matrix used by ropls (row = features)
  X_unscaled,          # same rows/cols as Xmv but NOT z-scaled (log10 OK)
  prefix, outdir,
  topN = 15,
  title = "S-plot (OPLS-DA)",
  corr_thr = 0.5,
  cov_q = 0.95,
  color_by_vip = TRUE,
  vip_thr = 1.0,
  lfc_by_feature = NULL,   # named numeric by feature (e.g. vm_enr[[paste0("log2FC_OVR_final_", cl)]], names = feature ids)
  append_log = function(...) {}
) {
  if (!is.matrix(Xmv)) Xmv <- as.matrix(Xmv)
  if (!is.matrix(X_unscaled)) X_unscaled <- as.matrix(X_unscaled)
  if (is.null(rownames(Xmv))) rownames(Xmv) <- paste0("feat_", seq_len(nrow(Xmv)))
  if (is.null(rownames(X_unscaled))) rownames(X_unscaled) <- rownames(Xmv)
  X_unscaled <- X_unscaled[rownames(Xmv), , drop = FALSE]

  # Prefer computing cov/corr against predictive scores t1 on the UN-SCALED matrix
  t1 <- NULL
  for (cand in c("scoreMN","xScoreMN")) {
    t1 <- tryCatch(slot(model, cand)[, 1], error = function(e) NULL)
    if (!is.null(t1)) break
  }
  if (!is.null(t1)) {
    corv <- apply(X_unscaled, 1, function(v) suppressWarnings(stats::cor(v, t1, use = "pairwise.complete.obs")))
    covv <- apply(X_unscaled, 1, function(v) suppressWarnings(stats::cov(v, t1, use = "pairwise.complete.obs")))
  } else {
    # fallback to loadings if t1 not available
    p1  <- tryCatch(slot(model, "loadingMN")[, 1, drop = TRUE],    error = function(e) NULL)
    pc1 <- tryCatch(slot(model, "loadingCorMN")[, 1, drop = TRUE], error = function(e) NULL)
    if (is.null(p1) || is.null(pc1)) { append_log(paste0(prefix, ": cannot extract scores/loadings — skip S-plot")); return(invisible(NULL)) }
    if (!is.null(names(p1))  && all(rownames(Xmv) %in% names(p1)))  p1  <- p1[rownames(Xmv)]
    if (!is.null(names(pc1)) && all(rownames(Xmv) %in% names(pc1))) pc1 <- pc1[rownames(Xmv)]
    covv <- as.numeric(p1); corv <- as.numeric(pc1)
  }

  # VIP
  vip <- tryCatch(slot(model, "vipVn"), error = function(e) NULL)
  if (!is.null(vip)) {
    if (is.null(names(vip)) && length(vip) == nrow(Xmv)) names(vip) <- rownames(Xmv)
    vip <- vip[rownames(Xmv)]
    if (is.null(vip) || length(vip) != nrow(Xmv)) vip <- rep(NA_real_, nrow(Xmv))
  } else vip <- rep(NA_real_, nrow(Xmv))

  df <- data.frame(
    feature = rownames(Xmv),
    cov  = as.numeric(covv),
    corr = as.numeric(corv),
    VIP  = as.numeric(vip),
    stringsAsFactors = FALSE
  )

  # Strict orientation using provided log2FC vector (named by feature)
  if (!is.null(lfc_by_feature)) {
    lfc <- as.numeric(lfc_by_feature[df$feature])
    ok <- is.finite(lfc)
    if (sum(ok) >= 3) {
      s <- suppressWarnings(sign(stats::cor(lfc[ok], df$corr[ok], use = "pairwise.complete.obs")))
      if (is.finite(s) && s < 0) { df$cov <- -df$cov; df$corr <- -df$corr }
    }
  }

  # thresholds / candidates
  df$abs_prod <- abs(df$cov * df$corr)
  finite_mask <- is.finite(df$cov) & is.finite(df$corr)
  if (!any(finite_mask)) finite_mask <- rep(TRUE, nrow(df))
  cov_thr_val  <- stats::quantile(abs(df$cov[finite_mask]), probs = min(max(cov_q, 0), 0.999), na.rm = TRUE)
  corr_thr_val <- min(max(corr_thr, 0), 1)
  vip_ok <- is.finite(df$VIP) & (df$VIP >= vip_thr)
  if (!any(is.finite(df$VIP))) vip_ok <- rep(TRUE, nrow(df))
  thr_ok <- (abs(df$corr) >= corr_thr_val) & (abs(df$cov) >= cov_thr_val)
  df$VIP_ok <- vip_ok; df$thr_ok <- thr_ok; df$candidate <- vip_ok & thr_ok
  if (!is.null(lfc_by_feature)) {
    l2 <- as.numeric(lfc_by_feature[df$feature])
    agree <- sign(l2) == sign(df$corr)
    agree[!is.finite(agree)] <- FALSE
    df$candidate_strict <- df$candidate & agree
  } else {
    df$candidate_strict <- df$candidate
  }

  # exports
  readr::write_tsv(df,   file.path(outdir, paste0(prefix, "_Splot_data.tsv")))
  cand <- df[df$candidate, , drop = FALSE][order(df$abs_prod[df$candidate], decreasing = TRUE), , drop = FALSE]
  readr::write_tsv(cand, file.path(outdir, paste0(prefix, "_Splot_candidates.tsv")))
  ord <- order(df$abs_prod, decreasing = TRUE)
  top <- if (is.finite(topN) && topN > 0) df[head(ord, topN), , drop = FALSE] else df[0, , drop = FALSE]
  readr::write_tsv(top,  file.path(outdir, paste0(prefix, "_Splot_top", if (is.finite(topN) && topN > 0) topN else 0, ".tsv")))

  # ----- PLOT (ggplot2) -----
  suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(cowplot) })

  # xlim based on quantile but DO NOT winsorize points; coord_cartesian will clip.
    # x-limits: start from q99.5% but ALWAYS include all candidates
  x_q <- stats::quantile(abs(df$cov[finite_mask]), probs = 0.995, na.rm = TRUE)
  if (!is.finite(x_q) || x_q <= 0) x_q <- max(abs(df$cov[finite_mask]), na.rm = TRUE)
  if (nrow(cand)) {
    x_cand <- max(abs(cand$cov), na.rm = TRUE)
    if (is.finite(x_cand) && x_cand > x_q) x_q <- x_cand * 1.02  # small padding
  }
  x_lim <- c(-x_q, x_q)

  # Base scatter
  df$col_group <- "All"
  if (isTRUE(color_by_vip) && any(is.finite(df$VIP))) df$col_group[df$VIP_ok] <- "VIP>=thr"

  p <- ggplot(df, aes(x = cov, y = corr)) +
    geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.2, color = "grey60") +
    geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.2, color = "grey60") +
    geom_point(aes(shape = ifelse(candidate, "candidate", "all"),
                   color = col_group), size = 1.2, alpha = 0.9, stroke = 0.2, show.legend = FALSE) +
    scale_shape_manual(
      name = "Candidatee",
      breaks = c("all", "candidate"),
      labels = c("Feature", "Candidatees"),
      values = c("all" = 16, "candidate" = 21)
    ) +
    scale_color_manual(
      name = "VIP status",
      breaks = c("All", "VIP>=thr"),
      labels = c("Other features", sprintf("VIP ≥ %.2f", vip_thr)),
      values = c("All" = "grey40", "VIP>=thr" = "red")
    ) +
    geom_hline(yintercept = c(-corr_thr_val, corr_thr_val), linetype = "dashed", linewidth = 0.3, color = "grey50") +
    geom_vline(xintercept = c(-cov_thr_val, cov_thr_val), linetype = "dashed", linewidth = 0.3, color = "grey50") +
    labs(x = "Covariance p[1] (x)", y = "Correlation p(corr)[1] (y)", title = title) +
    coord_cartesian(xlim = x_lim, ylim = c(-1.02, 1.02), clip = "on") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")

  # Candidate rings
  if (nrow(cand)) {
    p <- p + ggplot2::geom_point(data = cand, aes(x = cov, y = corr), show.legend = FALSE,
                                 shape = 21, fill = "orange", color = "black", size = 2.2, stroke = 0.4,
                                 inherit.aes = FALSE)
  }

  # Labels for top points
  if (nrow(top)) {
    p <- p + ggrepel::geom_text_repel(
      data = top,
      aes(x = cov, y = corr, label = feature, show.legend = FALSE, inherit.aes = FALSE, color = "black"),
      max.overlaps = Inf, size = 3,
      box.padding = 0.25, point.padding = 0.15, min.segment.length = 0,
      segment.color = "grey50", segment.size = 0.3
    )
  }


  # Custom legend panel (right side) built with cowplot
  leg_df <- data.frame(
    x = c(0.05, 0.05, 0.05),
    y = c(0.82, 0.74, 0.60),
    label = c("Other features", sprintf("VIP \u2265 %.2f", vip_thr), "Candidatees"),
    shape = c(16, 16, 21),
    col = c("grey40", "red", "black"),
    fill = c(NA, NA, "orange"),
    stroke = c(0.3, 0.3, 0.4),
    stringsAsFactors = FALSE
  )

  
  # --- LEGEND: publication-ready (wide spacing, compact panel) ---
  legend_box <- ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = "white", colour = "grey70", alpha = 0.96, linewidth = 0.30) +
    annotate("text", x = 0.10, y = 0.95, label = "Legend",
             hjust = 0, vjust = 1, size = 4.8, fontface = "bold") +
    annotate("text", x = 0.10, y = 0.83, label = sprintf("|corr| ≥ %.2f", corr_thr_val),
             hjust = 0, vjust = 1, size = 3.5) +
    annotate("text", x = 0.10, y = 0.76, label = sprintf("|cov| ≥ q%d", as.integer(100*cov_q)),
             hjust = 0, vjust = 1, size = 3.5) +
    annotate("text", x = 0.10, y = 0.69, label = sprintf("Candidatees: %d", nrow(cand)),
             hjust = 0, vjust = 1, size = 3.5) +
    geom_point(aes(x = 0.12, y = 0.56), shape = 16, size = 3.2, colour = "grey40") +
    annotate("text", x = 0.19, y = 0.56, label = "Other features",
             hjust = 0, vjust = 0.5, size = 3.5) +
    geom_point(aes(x = 0.12, y = 0.43), shape = 16, size = 3.2, colour = "red") +
    annotate("text", x = 0.19, y = 0.43, label = sprintf("VIP ≥ %.2f", vip_thr),
             hjust = 0, vjust = 0.5, size = 3.5) +
    geom_point(aes(x = 0.12, y = 0.30), shape = 21, size = 3.2, fill = "orange",
               colour = "black", stroke = 0.48) +
    annotate("text", x = 0.19, y = 0.30, label = "Candidatees",
             hjust = 0, vjust = 0.5, size = 3.5) +
    xlim(0, 1) + ylim(0, 1) + theme_void() +
    theme(plot.margin = margin(4, 8, 6, 8))

  # Compose with legend inside (right, vertically centered)
  composed <- cowplot::ggdraw() +
    cowplot::draw_plot(p) +
    cowplot::draw_plot(legend_box, x = 0.84, y = 0.32, width = 0.155, height = 0.22)

  print(composed)
  grDevices::dev.off()
  # -- Compose plot + legend panel and export PNG --
  composed <- cowplot::ggdraw() +
    cowplot::draw_plot(p) +
    cowplot::draw_plot(legend_box, x = 0.80, y = 0.33, width = 0.22, height = 0.32)

  png(file.path(outdir, paste0(prefix, "_Splot.png")), width = getOption("FIG_SPLOT_W", 2500L), height = getOption("FIG_SPLOT_H", 1600L), res = .get_export_dpi(), units = "px")
  print(composed)
  grDevices::dev.off()


  invisible(df)
}


export_validation_metrics <- function(model, prefix, outdir, append_log = function(...) {}) {
  # Safely extract metrics from ropls model
  safe_slot <- function(obj, name) { tryCatch(slot(obj, name), error = function(e) NULL) }
  summ <- safe_slot(model, "summaryDF")
  r2cum <- tryCatch(as.numeric(summ[1, grep("R2", colnames(summ))[1]]), error = function(e) NA_real_)
  q2cum <- tryCatch(as.numeric(summ[1, grep("Q2", colnames(summ))[1]]), error = function(e) NA_real_)
  an <- tryCatch(anova(model), error = function(e) NULL)
  p_cvanova <- tryCatch(as.numeric(an[1, ncol(an)]), error = function(e) NA_real_)
  # Permutations
  permMN <- safe_slot(model, "permMN"); permCorVn <- safe_slot(model, "permCorVn")
  R2int <- NA_real_; Q2int <- NA_real_; R2slope <- NA_real_; Q2slope <- NA_real_
  if (!is.null(permMN) && !is.null(permCorVn)) {
    colR2 <- grep("R2", colnames(permMN), value = TRUE)[1]
    colQ2 <- grep("Q2", colnames(permMN), value = TRUE)[1]
    dfp <- data.frame(cor = as.numeric(permCorVn),
                      R2 = as.numeric(permMN[, colR2]),
                      Q2 = as.numeric(permMN[, colQ2]))
    fitR2 <- tryCatch(lm(R2 ~ cor, data = dfp), error = function(e) NULL)
    fitQ2 <- tryCatch(lm(Q2 ~ cor, data = dfp), error = function(e) NULL)
    if (!is.null(fitR2)) { R2int <- coef(fitR2)[1]; R2slope <- coef(fitR2)[2] }
    if (!is.null(fitQ2)) { Q2int <- coef(fitQ2)[1]; Q2slope <- coef(fitQ2)[2] }
  }
  # Write TSV
  out <- tibble::tibble(
    metric = c("R2_cum","Q2_cum","CVANOVA_p","Perm_R2_intercept","Perm_Q2_intercept","Perm_R2_slope","Perm_Q2_slope"),
    value  = c(r2cum, q2cum, p_cvanova, R2int, Q2int, R2slope, Q2slope)
  )
  readr::write_tsv(out, file.path(outdir, paste0(prefix, "_ModelValidation.tsv")))
  invisible(out)
}
export_perm_intercepts <- function(model, prefix, outdir, append_log = function(...) {}) {
  sn <- slotNames(model)
  permMN <- tryCatch(if ("permMN" %in% sn) slot(model, "permMN") else NULL, error = function(e) NULL)
  permCorVn <- tryCatch(if ("permCorVn" %in% sn) slot(model, "permCorVn") else NULL, error = function(e) NULL)
  summ <- tryCatch(slot(model, "summaryDF"), error = function(e) NULL)
  if (is.null(permMN) || is.null(permCorVn) || is.null(summ) || nrow(permMN) == 0) { append_log(paste0(prefix, ": no permutation info")); return(invisible(NULL)) }
  colR2 <- grep("R2", colnames(permMN), value = TRUE)[1]
  colQ2 <- grep("Q2", colnames(permMN), value = TRUE)[1]
  if (is.na(colR2) || is.na(colQ2)) { append_log(paste0(prefix, ": R2/Q2 columns missing in permMN")); return(invisible(NULL)) }
  df <- data.frame(cor = as.numeric(permCorVn), R2 = as.numeric(permMN[, colR2]), Q2 = as.numeric(permMN[, colQ2]))
  fitR2 <- tryCatch(lm(R2 ~ cor, data = df), error = function(e) NULL)
  fitQ2 <- tryCatch(lm(Q2 ~ cor, data = df), error = function(e) NULL)
  R2orig <- suppressWarnings(as.numeric(summ[1, grep("R2", colnames(summ))[1]]))
  Q2orig <- suppressWarnings(as.numeric(summ[1, grep("Q2", colnames(summ))[1]]))
  out <- tibble::tibble(
    metric = c("R2", "Q2"),
    intercept = c(if (!is.null(fitR2)) coef(fitR2)[1] else NA_real_, if (!is.null(fitQ2)) coef(fitQ2)[1] else NA_real_),
    slope = c(if (!is.null(fitR2)) coef(fitR2)[2] else NA_real_, if (!is.null(fitQ2)) coef(fitQ2)[2] else NA_real_),
    p_empirical = c(mean(df$R2 >= R2orig, na.rm = TRUE), mean(df$Q2 >= Q2orig, na.rm = TRUE)),
    original = c(R2orig, Q2orig)
  )
  readr::write_tsv(out, file.path(outdir, paste0(prefix, "_PermIntercepts.tsv")))
}


# --- helper: simple vertical repel for labels ---
repel_vertical <- function(x, y, min_dy = 0.03, lo=-0.98, hi=0.98) {
  o <- order(y)
  yy <- y[o]
  # forward pass
  for (i in 2:length(yy)) {
    if (!is.finite(yy[i-1]) || !is.finite(yy[i])) next
    if (yy[i] - yy[i-1] < min_dy) yy[i] <- yy[i-1] + min_dy
  }
  # backward pass
  for (i in (length(yy)-1):1) {
    if (!is.finite(yy[i+1]) || !is.finite(yy[i])) next
    if (yy[i+1] - yy[i] < min_dy) yy[i] <- yy[i+1] - min_dy
  }
  yy <- pmin(pmax(yy, lo), hi)
  res <- numeric(length(y)); res[o] <- yy; res
}
write_html_report <- function(outdir, params) {
  images <- list.files(outdir, pattern = "\\.png$", full.names = FALSE)
  html <- c(
    "<!doctype html><html><head><meta charset='utf-8'><title>Post-Run Report</title>",
    "<style>body{font-family:sans-serif;margin:24px;} table{border-collapse:collapse;} td,th{border:1px solid #ccc;padding:6px;} img{max-width:100%;height:auto;border:1px solid #eee;margin:8px 0;} h2{margin-top:28px;}</style>",
    "</head><body>",
    "<h1>Post-Run Report</h1>",
    "<h2>Parameters</h2>",
    "<table><tbody>",
    paste0("<tr><th>", names(params), "</th><td>", unlist(params), "</td></tr>", collapse = ""),
    "</tbody></table>",
    "<h2>Figures</h2>"
  )
  for (img in images) html <- c(html, sprintf("<div><h3>%s</h3><img src='%s'></div>", img, img))
  html <- c(html, "</body></html>")
  writeLines(html, file.path(outdir, "PostRun_Report.html"))
}

# Robust signature extraction (biosigner versions differ)
extract_biosigner_signatures <- function(bobj, levels = c("S","AS")) {
  # Try exported helpers
  try_list <- list(
    function(z) getFromNamespace("significantFeatures","biosigner")(z, levels),
    function(z) getFromNamespace("signifFeatures","biosigner")(z, levels)
  )
  for (f in try_list) {
    res <- tryCatch(f(bobj), error = function(e) e)
    if (!inherits(res, "error")) return(res)
  }
  # Try direct slots
  sig <- NULL
  sn <- methods::slotNames(bobj)
  cand <- intersect(sn, c("signaturesLs","signatures","signLs","sign"))
  for (s in cand) {
    obj <- tryCatch(methods::slot(bobj, s), error = function(e) NULL)
    if (!is.null(obj)) sig <- obj
  }
  if (is.null(sig)) return(NULL)
  # Flatten list structure to named list of character vectors
  flatten_chars <- function(x) {
    out <- list()
    if (is.character(x)) return(list(unknown = x))
    if (is.list(x)) {
      for (nm in names(x)) {
        sub <- flatten_chars(x[[nm]])
        for (k in names(sub)) out[[paste(nm, k, sep = "_")]] <- sub[[k]]
      }
    }
    out
  }
  res <- flatten_chars(sig)
  # Keep unique char vectors
  res <- res[vapply(res, function(v) is.character(v) && length(v) > 0, logical(1))]
  if (!length(res)) return(NULL)
  res
}

# VIP selection helper
select_vip_features <- function(vipVn, mode = c("topN","threshold","cumprop"), topN = 15, thr = 1.0, cum_prop = 0.8, minN = 5) {
  mode <- match.arg(mode); v <- sort(vipVn, decreasing = TRUE)
  feats <- switch(mode,
    "topN" = names(v)[seq_len(min(length(v), topN))],
    "threshold" = names(v)[which(v >= thr)],
    "cumprop" = { cs <- cumsum(v); pr <- cs / max(cs); names(v)[which(pr <= cum_prop)] }
  )
  if (length(feats) < minN) feats <- names(v)[seq_len(min(minN, length(v)))]
  list(features = feats, vip_sorted = v)
}

# VIP dotplot + class log2FC heatmap (returns selected features)
vip_dot_heatmap <- function(vipVn, select_mode, topN, thr, cum_prop, vm, name_col, classes, Xraw, out_png, title = "VIP top features", pal = "Blue-White-Red", hm_ratio = getOption("vip_hm_ratio", 0.9)) {
  if (length(vipVn) == 0) return(invisible(character(0)))
  sel <- select_vip_features(vipVn, mode = select_mode, topN = topN, thr = thr, cum_prop = cum_prop, minN = 5)
  feats <- sel$features; vip_sorted <- sel$vip_sorted
  vip_df <- tibble::tibble(feature = factor(feats, levels = rev(feats)), VIP = as.numeric(vip_sorted[feats]))

  # Per-class log2FC matrix
  cls_levels <- levels(classes)
  if (!is.null(vm) && !is.na(name_col) && all(paste0("log2FC_OVR_final_", cls_levels) %in% names(vm))) {
    mat <- as.matrix(vm[match(feats, vm[[name_col]]), paste0("log2FC_OVR_final_", cls_levels), drop = FALSE])
    rownames(mat) <- feats; colnames(mat) <- cls_levels
  } else {
    mat <- log2fc_by_class(Xraw, classes)[feats, , drop = FALSE]
  }

  ylev <- levels(vip_df$feature)
  p1 <- ggplot(vip_df, aes(x = VIP, y = feature)) +
    geom_point(size = 2) +
    scale_y_discrete(limits = ylev, drop = FALSE, expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(clip = "off") + labs(x = "VIP score", y = NULL) +
    theme(axis.text.y = element_text(size = 9), axis.ticks.length = grid::unit(0, "pt"),
          plot.margin = margin(t = 6, r = 8, b = 14, l = 14))

  palv <- get_div_palette(pal, 100)
  hdf <- as.data.frame(mat) %>% tibble::rownames_to_column("feature") %>% tidyr::pivot_longer(-feature, names_to = "class", values_to = "log2FC")
  hdf$feature <- factor(hdf$feature, levels = ylev)
  p2 <- ggplot(hdf, aes(x = class, y = feature, fill = log2FC)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_y_discrete(limits = ylev, drop = FALSE, expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(clip = "off") + scale_fill_gradientn(colors = palv) +
    labs(x = NULL, y = NULL, fill = "log2FC") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.ticks.length = grid::unit(0, "pt"),
          plot.margin = margin(t = 6, r = 14, b = 14, l = 2))

  final <- (p1 + p2 + patchwork::plot_layout(widths = c(1.8, hm_ratio))) + patchwork::plot_annotation(title = title)
  ggsave(filename = out_png, plot = final, width = 10, height = 7.0, dpi = 280)
  readr::write_tsv(tibble::tibble(feature = feats, VIP = as.numeric(vip_sorted[feats])), sub("\\.png$", "_selected.tsv", out_png))
  return(invisible(feats))
}

# Improved Venn (centered title, wider margins)
venn_from_lists <- function(lst, file_png, width=8.5, height=6.5, dpi=220, title = "Venn") {
  p <- ggVennDiagram::ggVennDiagram(lst) +
    ggtitle(title) +
    coord_cartesian(clip = "off") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.margin = margin(20, 50, 20, 50))
  ggsave(file_png, plot = p, width = width, height = height, dpi = dpi)
}

# Per-feature distribution plots (box/violin + jitter)
plot_feature_distributions <- function(Xraw, classes, features, out_png, title = "Signature features", type = c("box","violin"), ncol = 4) {
  type <- match.arg(type)
  feats <- intersect(features, rownames(Xraw))
  if (!length(feats)) return(invisible(NULL))
  M <- log10(pmax(Xraw[feats, , drop = FALSE], .Machine$double.eps))
  df <- as.data.frame(t(M))
  df$sample <- rownames(df)
  df$class <- as.character(classes)[match(df$sample, colnames(M))]
  long <- tidyr::pivot_longer(df, cols = all_of(feats), names_to = "feature", values_to = "log10int")
  g <- ggplot(long, aes(x = class, y = log10int)) +
    { if (type == "violin") ggplot2::geom_violin(trim = FALSE, alpha = 0.9) else ggplot2::geom_boxplot(outlier.shape = NA) } +
    geom_jitter(width = 0.1, alpha = 0.5, size = 0.8) +
    facet_wrap(~feature, scales = "free_y", ncol = ncol) +
    labs(x = NULL, y = "log10(intensity)", title = title)
  n_panels <- length(unique(long$feature))
  w <- max(6, 3 * ncol); h <- max(4, 2 + 2.2 * ceiling(n_panels / ncol))
  ggsave(out_png, g, width = w, height = h, dpi = 220)
}

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("JaFaR - Untargeted Metabolomics Multivariate Pipeline"),
  tags$style(".folder-display{padding:.4em .6em;border:1px solid #ccc;border-radius:4px;background:#fff; white-space:nowrap; overflow:hidden; text-overflow:ellipsis;} #time_status { color:#444; } .preset-row button{margin-right:6px;} .small-note{color:#777;font-size:12px;}"),
  sidebarLayout(
    sidebarPanel(width = 4,
      h4("PNG export — DPI/Anti-crash"),
      numericInput("png_dpi", "PNG DPI (min 300)

", value = 420, min = 300, step = 25),
      numericInput("ragg_max_dim", "ragg.max_dim (anti-crash)", value = 50000, min = 10000, step = 5000),
      helpText("DPI does not change the geometry (pixels). It only writes pixel density inside the PNG file."),
      tags$hr(),
      fileInput("dm", "Post_Univariate_dataMatrix_imputed.tsv", accept = c(".tsv",".tab",".tabular",".csv")),
      fileInput("sm", "Post_Univariate_sampleMetadata_enriched.tsv", accept = c(".tsv",".tab",".tabular",".csv")),
      fileInput("vm", "Post_Univariate_variableMetadata_enriched.tsv", accept = c(".tsv",".tab",".tabular",".csv")),
      shinyDirButton("outdir_btn", "Output folder...", "Choose", multiple = FALSE),
      uiOutput("outdir_display"),
      hr(),
      div(class="preset-row",
        actionButton("preset_w4m", "W4M settings", class = "btn btn-success"),
        actionButton("preset_ropls", "ropls-only", class = "btn btn-secondary")
      ),
      br(),
      checkboxInput("do_log10", "Log10-transform for multivariate & pairwise", value = TRUE),
      selectInput("scaling", "Scaling (multivariate)", choices = c("pareto","unitVar","none"), selected = "unitVar"),
      selectInput("ropls_scaleC", "ropls::opls scaleC (if no pre-scaling)", choices = c("none","center","pareto","standard"), selected = "none"),
      selectInput("pls_plot_type", "Plot type (PLS/OPLS-DA)", choices = c("x-score","xy-score"), selected = "x-score"),
      numericInput("splot_topN", "S-plot labels (top N)", value = 20, min = 0, step = 1),
      sliderInput("splot_corr_thr", "|corr| threshold", min = 0, max = 1, value = 0.7, step = 0.05),
      sliderInput("splot_cov_q", "|cov| threshold (quantile)", min = 0.5, max = 0.999, value = 0.99, step = 0.01),
      checkboxInput("splot_color_vip", "Color by VIP if available", value = TRUE),
      numericInput("splot_vip_thr", "VIP threshold (highlight)", value = 1.0, min = 0, step = 0.1),

      hr(),
      h5("S-plot - significance filters"),
      checkboxInput("splot_use_sig_filters", "Enable (q_global <= ... AND q_posthoc_agg_OVR_Class <= ... AND |log2FC_OVR_final_Class| >= ...)", value = TRUE),
      numericInput("splot_sig_q_global", "q_global <=", value = 0.05, min = 0, max = 1, step = 0.001),
      numericInput("splot_sig_q_ovr", "q_posthoc_agg_OVR <=", value = 0.05, min = 0, max = 1, step = 0.001),
      numericInput("splot_sig_abs_l2fc", "|log2FC_OVR_final| >=", value = 1.0, min = 0, step = 0.1),
      div(class = "small-note", "These filters apply per class (OVR) using the columns q_global, q_posthoc_agg_OVR_<class>, log2FC_OVR_final_<class> from variableMetadata."),
          selectInput("pca_predI", "PCA - predictive components", choices = c("Auto (NA)" = "NA", "1","2","3","4","5","6","7","8","9","10"), selected = "NA"),
      selectInput("pls_predI", "PLS-DA - predictive components", choices = c("Auto (NA)" = "NA", "1","2","3","4","5","6","7","8","9","10"), selected = "NA"),
      numericInput("orthoI", "Orthogonal components (OPLS-DA)", value = 1, min = 0, step = 1),
      numericInput("permI", "Permutations (PLS-DA/OPLS-DA)", value = 1000, min = 0, step = 10),
      numericInput("cv_segments", "Number of cross-validation segments", value = 6, min = 2, step = 1),
      helpText("Must be less than or equal to the number of replicates"),
      numericInput("opls_predI", "OPLS-DA - predictive components", value = 1, min = 1, step = 1),
      helpText("1) PCA and PLS(-DA): NA can be selected to get a suggestion of the optimal number of predictive components; 2) OPLS(-DA) modeling: select 1 predictive component"),
      numericInput("q_thr", "q threshold (BH)", value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput("lfc_thr", "|log2FC| threshold (volcano)", value = 1.0, min = 0, step = 0.1),
      numericInput("volcano_top_label", "Volcano: top labels (N)", value = 10, min = 0, step = 1),
      selectInput("div_palette", "Diverging palette (volcano & VIP heatmap)", choices = c("Blue-White-Red","RdBu","RdYlBu","PRGn","PuOr","Viridis (sequential)","Plasma (sequential)"), selected = "Blue-White-Red"),
      numericInput("biosigner_boot", "Biosigner bootstraps", value = 50, min = 20, step = 10),
      radioButtons("biosigner_sigset", "Biosigner signature set", choices = c("S", "S+A"), selected = "S+A", inline = TRUE),
      selectInput("biosigner_plot_type", "Biosigner per-feature plots", choices = c("box","violin"), selected = "box"),
      numericInput("biosigner_plot_n", "Max signature features per class", value = 12, min = 1, step = 1),
      checkboxInput("do_pairwise_volcano", "Pairwise volcanoes (A vs B)", value = TRUE),
      checkboxInput("do_ovr_volcano", "One-vs-Rest volcanoes (per class)", value = TRUE),
      hr(),
      h4("Figure sizes (px)"),
      # --- Volcano & S-plot ---
      numericInput("px_volcano_w", "Volcano width (px)", value = 3000, min = 800, step = 50),
      numericInput("px_volcano_h", "Volcano height (px)", value = 3000, min = 600, step = 50),
      # --- PCA / PLS-DA / OPLS(-DA) base plots ---
      numericInput("px_pca_w",     "PCA width (px)",      value = 3000, min = 800, step = 50),
      numericInput("px_pca_h",     "PCA height (px)",     value = 3000, min = 600, step = 50),
      numericInput("px_pls_w",     "PLS/OPLS width (px)", value = 3000, min = 800, step = 50),
      numericInput("px_pls_h",     "PLS/OPLS height (px)",value = 3000, min = 600, step = 50),
      numericInput("px_perm_w",    "Permutations width (px)", value = 3000, min = 800, step = 50),
      numericInput("px_perm_h",    "Permutations height (px)", value = 3000, min = 600, step = 50),
      numericInput("px_cm_w",      "Confusion matrix width (px)", value = 3000, min = 600, step = 50),
      numericInput("px_cm_h",      "Confusion matrix height (px)", value = 3000, min = 500, step = 50),
      numericInput("px_pred_w",    "Predictions/Perf. width (px)", value = 3000, min = 800, step = 50),
      numericInput("px_pred_h",    "Predictions/Perf. height (px)", value = 3000, min = 600, step = 50),
      # --- Heatmaps / VIP / Biosigner ---
      numericInput("px_vip_w",     "VIP heatmap width (px)", value = 3000, min = 800, step = 50),
      numericInput("px_vip_h",     "VIP heatmap height (px)", value = 3000, min = 600, step = 50),
      numericInput("px_univhm_w",  "UNIV_hit heatmap width (px)", value = 3000, min = 800, step = 50),
      numericInput("px_univhm_h",  "UNIV_hit heatmap height (px)", value = 3000, min = 600, step = 50),
      numericInput("px_biofeat_w", "Biosigner per-feature width (px)", value = 3000, min = 800, step = 50),
      numericInput("px_biofeat_h", "Biosigner per-feature height (px)", value = 3000, min = 600, step = 50),
      # --- Fallback (others) ---
      numericInput("px_default_w", "Default width (px)", value = 3000, min = 600, step = 50),
      numericInput("px_default_h", "Default height (px)", value = 3000, min = 400, step = 50),
      numericInput("px_venn_w", "Venn width (px)", value = 3000, min = 800, step = 50),
      numericInput("px_venn_h", "Venn height (px)", value = 3000, min = 600, step = 50),
      numericInput("px_splot_w", "S-plot width (px)", value = 3000, min = 600, step = 50),
      numericInput("px_splot_h", "S-plot height (px)", value = 3000, min = 400, step = 50),
      numericInput("px_volcano_w", "Volcano width (px)", value = 3000, min = 800, step = 50),
      numericInput("px_volcano_h", "Volcano height (px)", value = 3000, min = 600, step = 50),

      hr(),
      h4("UNIV_hit heatmaps"),
      numericInput("hm_univ_top", "Top N features by variance (fallback)", value = 200, min = 2, step = 10),
      numericInput("hm_univ_top_per_class", "Top K per class (most up-regulated)", value = 20, min = 0, step = 1),
      checkboxInput("hm_univ_scale", "z-score by feature", value = TRUE),
      checkboxInput("hm_univ_cluster_rows", "Cluster rows (features)", value = TRUE),
      checkboxInput("hm_univ_cluster_cols", "Cluster columns (samples/classes)", value = TRUE),
      numericInput("hm_univ_font_row", "Row label size (features)", value = 6, min = 4, step = 1),
      numericInput("hm_univ_font_col", "Col label size (samples/classes)", value = 10, min = 6, step = 1),
      sliderInput("hm_univ_w", "Width (inches)", min = 4, max = 16, value = 10, step = 0.5),
      sliderInput("hm_univ_h", "Height (inches)", min = 4, max = 16, value = 10, step = 0.5),
      selectInput("hm_palette", "Heatmap palette", choices = c("Blue-White-Red","RdBu","RdYlBu","PRGn","PuOr","Viridis (sequential)","Plasma (sequential)"), selected = "Blue-White-Red"),
      uiOutput("hm_ann_cols_ui"),
      hr(),
      h4("VIP options"),
      sliderInput("vip_hm_ratio", "VIP heatmap width ratio", min = 0.6, max = 2.5, value = 0.9, step = 0.1),
      radioButtons("vip_mode", "VIP selection", inline = FALSE,
                   choices = c("Top N" = "topN", "VIP >= threshold" = "threshold", "Cumulative VIP <= %" = "cumprop"), selected = "topN"),
      numericInput("vip_topN", "Top N VIP (if Top N)", value = 20, min = 5, step = 1),
      numericInput("vip_thr", "VIP threshold (if threshold)", value = 1.0, min = 0, step = 0.1),
      sliderInput("vip_cum", "Cumulative VIP % (if cumulative)", min = 0.5, max = 0.99, value = 0.8, step = 0.01),
      div(class = "small-note", "Tip: VIP >= 1 is common; 80% cumulative helps compact the list."),
      actionButton("run", "Run", class = "btn btn-primary btn-block")
    ) ,
    mainPanel(width = 8,
      shinyWidgets::progressBar(id = "pb", value = 0, total = 18, title = "Idle...", display_pct = TRUE, striped = TRUE, status = "info", size = "sm"),
      textOutput("step_status"),
      textOutput("time_status"),
      verbatimTextOutput("log")
    )

  )
)

# ---------- Server ----------
server <- function(input, output, session) {

# -- Sync UI -> options for DPI and figure sizes (single source of truth) --
observe({
  options(
    PNG_EXPORT_DPI = max(300L, as.integer(input$png_dpi %||% 300L)),

    FIG_VOLCANO_W = as.integer(input$px_volcano_w %||% 1400L),
    FIG_VOLCANO_H = as.integer(input$px_volcano_h %||% 1100L),

    FIG_SPLOT_W   = as.integer(input$px_splot_w   %||% 2500L),
    FIG_SPLOT_H   = as.integer(input$px_splot_h   %||% 1600L),

    FIG_PCA_W     = as.integer(input$px_pca_w     %||% 1400L),
    FIG_PCA_H     = as.integer(input$px_pca_h     %||% 1100L),

    FIG_PLS_W     = as.integer(input$px_pls_w     %||% 1400L),
    FIG_PLS_H     = as.integer(input$px_pls_h     %||% 1100L),

    FIG_PERM_W    = as.integer(input$px_perm_w    %||% 1400L),
    FIG_PERM_H    = as.integer(input$px_perm_h    %||% 1100L),

    FIG_CM_W      = as.integer(input$px_cm_w      %||% 1200L),
    FIG_CM_H      = as.integer(input$px_cm_h      %||%  900L),

    FIG_PRED_W    = as.integer(input$px_pred_w    %||% 1400L),
    FIG_PRED_H    = as.integer(input$px_pred_h    %||% 1100L),

    FIG_VIP_W     = as.integer(input$px_vip_w     %||% 1400L),
    FIG_VIP_H     = as.integer(input$px_vip_h     %||% 1100L),

    FIG_UNIVHM_W  = as.integer(input$px_univhm_w  %||% 1400L),
    FIG_UNIVHM_H  = as.integer(input$px_univhm_h  %||% 1100L),

    FIG_BIOFEAT_W = as.integer(input$px_biofeat_w %||% 1400L),
    FIG_BIOFEAT_H = as.integer(input$px_biofeat_h %||% 1100L),

    FIG_VENN_W    = as.integer(input$px_venn_w    %||% 1400L),
    FIG_VENN_H    = as.integer(input$px_venn_h    %||% 1100L),

    FIG_DEFAULT_W = as.integer(input$px_default_w %||% 1400L),
    FIG_DEFAULT_H = as.integer(input$px_default_h %||% 1100L)
  )
})






# -- Sync UI -> options for DPI and figure sizes --



  
  

  # --- Biosigner coherence auto-harmonizer (u10d) ---
  

  volumes <- tryCatch(shinyFiles::getVolumes()(), error = function(e) NULL, warning = function(w) NULL)
  if (is.null(volumes) || length(volumes) == 0) {
    drv <- paste0(LETTERS, ":/"); drv <- drv[file.exists(drv) | dir.exists(drv)]
    if (length(drv) > 0) { names(drv) <- gsub("/$", "", drv); volumes <- c(Home = normalizePath("~", winslash = "/", mustWork = FALSE), drv) }
    else { volumes <- c(Home = normalizePath("~", winslash = "/", mustWork = FALSE)) }
  }

  rv <- reactiveValues(
    log = "", t0 = NULL, outdir = normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    sm_cols = NULL, vip_pls_feats = character(0), vip_ovr_feats = list(),
    biosigner_sigs = list(), pred_plsda = NULL, pred_ovr = list(),
    splot_ovr_union = NULL,
    splot_by_class = list()
  )
  output$outdir_display <- renderUI(tagList(tags$label("Output folder"), div(class="folder-display", rv$outdir)))
  shinyFiles::shinyDirChoose(input, "outdir_btn", roots = volumes, session = session)
  observeEvent(input$outdir_btn, {
    sel <- try(shinyFiles::parseDirPath(volumes, input$outdir_btn), silent = FALSE)
    if (!inherits(sel, "try-error") && length(sel) > 0 && nchar(sel[1]) > 0) {
      rv$outdir <- as.character(sel[1])
      output$outdir_display <- renderUI(tagList(tags$label("Output folder"), div(class="folder-display", rv$outdir)))
    }
  })

  # Update annotation column choices when sampleMetadata chosen
  observeEvent(input$sm, {
    if (!is.null(input$sm$datapath) && file.exists(input$sm$datapath)) {
      tmp <- try(read_tabflex(input$sm$datapath), silent = FALSE)
      if (!inherits(tmp, "try-error")) rv$sm_cols <- names(tmp)
    }
  })
  output$hm_ann_cols_ui <- renderUI({
    cols <- rv$sm_cols; if (is.null(cols)) cols <- character(0)
    # Do not allow re-adding the class/group columns (prevents duplicated legends)
    cols <- setdiff(cols, c("sample_name","Sample","Sample_name","sample","id","ID","class","Class","group","Group"))
    selectizeInput("hm_ann_cols", "Column annotations (from sampleMetadata)", choices = cols, multiple = TRUE,
                   options = list(placeholder = 'Select 0..n columns (batch, day, etc.)'))
  })

  # Presets
  observeEvent(input$preset_w4m, {
    updateCheckboxInput(session, "do_log10", value = TRUE)
    updateSelectInput(session, "scaling", selected = "pareto")
    updateSelectInput(session, "ropls_scaleC", selected = "none")
    updateSelectInput(session, "pls_plot_type", selected = "x-score")
    updateSelectInput(session, "pca_predI", selected = "NA")
    updateSelectInput(session, "pls_predI", selected = "NA")
    updateNumericInput(session, "orthoI", value = 1)
    updateNumericInput(session, "opls_predI", value = 1)
    updateNumericInput(session, "permI", value = 1000)
    updateNumericInput(session, "q_thr", value = 0.05)
    updateNumericInput(session, "lfc_thr", value = 1.0)
    updateNumericInput(session, "volcano_top_label", value = 10)
    updateSelectInput(session, "div_palette", selected = "Blue-White-Red")
    updateNumericInput(session, "biosigner_boot", value = 50)
    updateRadioButtons(session, "biosigner_sigset", selected = "S+A")
    updateSliderInput(session, "vip_hm_ratio", value = 0.9)
    updateRadioButtons(session, "vip_mode", selected = "topN")
    updateNumericInput(session, "vip_topN", value = 15)
    updateNumericInput(session, "vip_thr", value = 1.0)
    updateSliderInput(session, "vip_cum", value = 0.8)
  })
  observeEvent(input$preset_ropls, {
    updateCheckboxInput(session, "do_log10", value = TRUE)
    updateSelectInput(session, "scaling", selected = "none")
    updateSelectInput(session, "ropls_scaleC", selected = "pareto")
    updateSelectInput(session, "pls_plot_type", selected = "x-score")
    updateSelectInput(session, "pca_predI", selected = "NA")
    updateSelectInput(session, "pls_predI", selected = "NA")
    updateNumericInput(session, "orthoI", value = 1)
    updateNumericInput(session, "opls_predI", value = 1)
    updateNumericInput(session, "permI", value = 1000)
    updateNumericInput(session, "q_thr", value = 0.05)
    updateNumericInput(session, "lfc_thr", value = 1.0)
    updateNumericInput(session, "volcano_top_label", value = 10)
    updateSelectInput(session, "div_palette", selected = "Blue-White-Red")
    updateNumericInput(session, "biosigner_boot", value = 50)
    updateRadioButtons(session, "biosigner_sigset", selected = "S+A")
    updateSliderInput(session, "vip_hm_ratio", value = 0.9)
    updateRadioButtons(session, "vip_mode", selected = "topN")
    updateNumericInput(session, "vip_topN", value = 15)
    updateNumericInput(session, "vip_thr", value = 1.0)
    updateSliderInput(session, "vip_cum", value = 0.8)
  })

  t_total <- 18
  append_log <- function(txt) { stamp <- format(Sys.time(), "%H:%M:%S"); rv$log <- paste0(rv$log, sprintf("%s - %s\n", stamp, txt)); output$log <- renderText(rv$log) }

  # Utility: Enrich and export the three tables
  enrich_and_export_tables <- function(Xraw, classes, vm, sm, rid_col, sid_col, outdir,
                                     vip_pls_feats, vip_ovr_feats, biosigner_sigs,
                                     pred_plsda, pred_ovr, splot_by_class = NULL) {
  # --- Type safety on keys ---
  name_col_vm <- detect_col(vm, c("name","feature_id","feature","id")); if (is.na(name_col_vm)) name_col_vm <- "name"
  vm_enr <- vm
  vm_enr <- ensure_char_key(vm_enr, key = "name", source_col = name_col_vm)
  vm_enr$name <- as.character(vm_enr$name)

  if (!(sid_col %in% names(sm))) {
    sid_guess <- detect_col(sm, c("sample_name","Sample","Sample_name","sample","id","ID","filename","file"))
    if (!is.na(sid_guess)) sid_col <- sid_guess
  }
  sm[[sid_col]] <- as.character(sm[[sid_col]])

  if (is.null(rownames(Xraw))) stop("Xraw must have rownames equal to feature identifiers")

  # --- log2FC per class ---
  l2fc_mat <- log2fc_by_class(Xraw, classes)
  l2fc_df  <- as.data.frame(l2fc_mat)
  l2fc_df  <- tibble::rownames_to_column(l2fc_df, "name")
  l2fc_df$name <- as.character(l2fc_df$name)
  vm_enr <- dplyr::left_join(vm_enr, l2fc_df %>% dplyr::rename_with(~paste0("calc_log2FC_OVR_", .x), -name), by = "name")

  # --- VIP flags ---
  vm_enr$VIP_PLSDA_selected <- FALSE
  if (length(vip_pls_feats)) vm_enr$VIP_PLSDA_selected <- vm_enr$name %in% vip_pls_feats
  if (length(vip_ovr_feats)) {
    for (cl in names(vip_ovr_feats)) {
      coln <- paste0("VIP_OVR_", stubify(cl), "_selected")
      vm_enr[[coln]] <- vm_enr$name %in% vip_ovr_feats[[cl]]
    }
  }

  # --- S-plot annotations ---
  if (!is.null(splot_by_class) && length(splot_by_class)) {
    name_col_vm2 <- detect_col(vm_enr, c("name","feature_id","id")); if (is.na(name_col_vm2)) name_col_vm2 <- "name"
    for (cl in names(splot_by_class)) {
      dfsp <- splot_by_class[[cl]]
      if (is.null(dfsp) || !is.data.frame(dfsp) || !all(c("feature","cov","corr") %in% names(dfsp))) next
      keep_cols <- intersect(c("feature","cov","corr","VIP","thr_ok","VIP_ok","candidate","candidate_strict"), names(dfsp))
      dfsp <- dfsp[, keep_cols, drop = FALSE]
      colnames(dfsp)[colnames(dfsp) == "feature"] <- name_col_vm2
      other_cols <- setdiff(colnames(dfsp), name_col_vm2)
      pref <- paste0("Splot_OVR_", stubify(cl), "_")
      colnames(dfsp)[match(other_cols, colnames(dfsp))] <- paste0(pref, other_cols)
      vm_enr[[name_col_vm2]] <- as.character(vm_enr[[name_col_vm2]])
      dfsp[[name_col_vm2]]   <- as.character(dfsp[[name_col_vm2]])
      vm_enr <- dplyr::left_join(vm_enr, dfsp, by = setNames(name_col_vm2, name_col_vm2))
    }
  }

  # --- Sample metadata enrichment ---
  sm_enr <- sm
  if (!(sid_col %in% names(sm_enr))) sm_enr[[sid_col]] <- sm_enr[[detect_col(sm_enr, c("sample_name","Sample","Sample_name","sample","id","ID","filename","file"))]]
  sm_enr[[sid_col]] <- as.character(sm_enr[[sid_col]])
  if (!is.null(pred_plsda) && is.data.frame(pred_plsda) && "sample" %in% names(pred_plsda)) {
    sm_enr <- sm_enr %>% dplyr::left_join(pred_plsda %>% dplyr::select(sample, predicted) %>% dplyr::rename(PLSDA_predicted = predicted),
                                          by = setNames("sample", sid_col))
  }
  if (length(pred_ovr)) {
    for (cl in names(pred_ovr)) {
      dfp <- pred_ovr[[cl]]
      if (is.data.frame(dfp) && "sample" %in% names(dfp)) {
        coln <- paste0("OPLSDA_OVR_", stubify(cl), "_predicted")
        sm_enr <- sm_enr %>% dplyr::left_join(dfp %>% dplyr::select(sample, predicted) %>% dplyr::rename(!!coln := predicted),
                                              by = setNames("sample", sid_col))
      }
    }
  }

  # --- dataMatrix enriched ---
  dm_enr <- data.frame(name = as.character(rownames(Xraw)), stringsAsFactors = FALSE)
  dm_enr <- cbind(dm_enr, as.data.frame(Xraw[, , drop = FALSE], check.names = FALSE))

  # --- Exports ---
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  readr::write_tsv(vm_enr, file.path(outdir, "PostRun_variableMetadata_enriched_plus.tsv"))
  readr::write_tsv(sm_enr, file.path(outdir, "PostRun_sampleMetadata_enriched_plus.tsv"))
  readr::write_tsv(dm_enr, file.path(outdir, "PostRun_dataMatrix_enriched_plus.tsv"))

  # Export combined MZmine quant enriched CSV
  try({
    out_csv <- file.path(outdir, "MZmine_quant_enriched_post_univariate_and_multivariate.csv")
    make_mzmine_quant_enriched_csv(dm_enr = dm_enr, sm_enr = sm_enr, vm_enr = vm_enr, out_csv = out_csv)
  }, silent = TRUE)


  readr::write_tsv(dm_enr, file.path(outdir, "dataMatrix_enriched.tsv"))
  readr::write_tsv(sm_enr, file.path(outdir, "sampleMetadata_enriched.tsv"))
  readr::write_tsv(vm_enr, file.path(outdir, "variableMetadata_enriched.tsv"))

  # --- BIOSIGNER TXT↔PNG COHERENCE (u6) ---
  # Rebuild PNG and a normalized TSV *from the pretty-printed *_model.txt*
  # so that text and plot always match 1:1.
  try({
    # Helpers inlined here to avoid sourcing issues
    .bio_levels <- function() c('S','A','B','C','D','E')
    .bio_norm <- function(x) { x <- toupper(trimws(as.character(x))); x[!x %in% .bio_levels()] <- NA_character_; factor(x, levels = .bio_levels(), ordered = TRUE) }
    .bio_palette <- function() c('#1B9E77','#7CAE61','#D9EF8B','#FEE08B','#FC8D62','#D73027')
    .bio_parse_model_txt <- function(path) {
      stopifnot(file.exists(path))
      txt <- readLines(path, warn = FALSE, encoding = 'UTF-8')
      start <- grep('^\\s*Significant features', txt)[1]
      if (is.na(start)) stop('Block \"Significant features\" not found in: ', path)
      hdr <- grep('^\\s*plsda\\s+randomforest\\s+svm\\s*$', txt, ignore.case = TRUE); hdr <- hdr[hdr > start][1]
      if (is.na(hdr)) stop('Header line (plsda randomforest svm) not found after \"Significant features\" block.')
      end <- grep('^\\s*Accuracy\\s*:', txt); end <- end[end > hdr]; end <- if (length(end)) end[1] - 1L else length(txt)
      lines <- txt[seq(hdr + 1L, end)]; lines <- lines[grepl('^\\s*\\d+', lines)]; if (!length(lines)) stop('No feature rows detected.')
      parse_row <- function(z) {
        m <- regexec('^\\s*(\\d+)\\s+\"?([A-Za-z])\"?\\s+\"?([A-Za-z])\"?\\s+\"?([A-Za-z])\"?\\s*$', z); g <- regmatches(z, m)[[1]]
        if (length(g) == 5) {
          data.frame(feature=g[2], plsda=toupper(g[3]), randomforest=toupper(g[4]), svm=toupper(g[5]), stringsAsFactors=FALSE)
        } else {
          toks <- gsub('[\\\"\\\']', '', z); toks <- unlist(strsplit(toks, '\\\\s+')); toks <- toks[nzchar(toks)]
          if (length(toks) >= 4 && grepl('^\\\\d+$', toks[1])) data.frame(feature=toks[1], plsda=toupper(substr(toks[2],1,1)), randomforest=toupper(substr(toks[3],1,1)), svm=toupper(substr(toks[4],1,1)), stringsAsFactors=FALSE) else NULL
        }
      }
      rows <- do.call(rbind, lapply(lines, parse_row)); if (is.null(rows) || !nrow(rows)) stop('Failed to parse feature rows: ', path)
      suppressWarnings(fnum <- as.numeric(rows$feature)); ord <- if (all(!is.na(fnum))) order(fnum, rows$feature) else order(rows$feature)
      rows <- rows[ord, , drop = FALSE]; rownames(rows) <- NULL; rows
    }
    .bio_build_matrix <- function(df) {
      tab <- data.frame(feature=as.character(df$feature),
                        plsda=.bio_norm(df$plsda),
                        randomforest=.bio_norm(df$randomforest),
                        svm=.bio_norm(df$svm),
                        stringsAsFactors=FALSE)
      suppressWarnings({ fnum <- suppressWarnings(as.numeric(tab$feature)) })
      ord <- if (all(!is.na(fnum))) order(fnum, tab$feature) else order(tab$feature)
      tab <- tab[ord, , drop = FALSE]; rownames(tab) <- tab$feature; tab$feature <- NULL; as.matrix(tab)
    }
    .bio_plot <- function(mat, path_png=NULL, title='Tiers of the selected features') {
      lvl <- .bio_levels(); z <- apply(mat, 2, function(col) as.integer(factor(col, levels=lvl, ordered=TRUE)))
      row_labels <- rownames(mat); col_labels <- colnames(mat); pal <- .bio_palette()
      z_rev <- z[nrow(z):1, , drop=FALSE]; row_labels_rev <- rev(row_labels)
      op <- par(no.readonly=TRUE); on.exit(par(op), add=TRUE); par(mar=c(5,6,3,6))
      if (!is.null(path_png)) png(path_png, width=1100, height=850, res=120)
      image(x=seq_len(ncol(z_rev)), y=seq_len(nrow(z_rev)), z=t(z_rev), axes=FALSE, xlab='', ylab='', main=title, col=pal, useRaster=TRUE, zlim=c(1, length(lvl)))
      box(lwd=3); axis(1, at=seq_len(ncol(z_rev)), labels=col_labels, las=2, tick=FALSE, line=-0.5); axis(2, at=seq_len(nrow(z_rev)), labels=row_labels_rev, las=2, tick=FALSE, line=-0.5)
      legend_x <- par('usr')[2] + 0.5; legend_y <- mean(par('usr')[3:4]); par(xpd=NA)
      legend(x=legend_x, y=legend_y, legend=lvl, fill=pal, bty='n', title=NULL, y.intersp=1.2, x.intersp=0.8, cex=1.1)
      if (!is.null(path_png)) dev.off(); invisible(path_png)
    }

    files <- list.files(outdir, pattern = '^Biosigner_OVR_internal_.*_model\\.txt$', full.names = TRUE)
    for (mtxt in files) {
      scope <- sub('^Biosigner_OVR_internal_(.+?)_model\\.txt$', '\\\\1', basename(mtxt))
      df_bio <- .bio_parse_model_txt(mtxt)
      mat    <- .bio_build_matrix(df_bio)
      out_png <- file.path(outdir, sprintf('Biosigner_OVR_internal_%s_plot.png', scope))
      out_tsv <- file.path(outdir, sprintf('Biosigner_OVR_internal_%s_model.tsv', scope))
      readr::write_tsv(data.frame(feature=rownames(mat),
                                  plsda=as.character(mat[, 'plsda']),
                                  randomforest=as.character(mat[, 'randomforest']),
                                  svm=as.character(mat[, 'svm']),
                                  stringsAsFactors=FALSE),
                       out_tsv)
      .bio_plot(mat, out_png, title='Tiers of the selected features')
    }
  }, silent = TRUE)

  # Late-write (U10): regenerate PNG/TSV from *_model.txt to ensure coherence
  try(later::later(function(){ try(biosigner_sync(outdir), silent = TRUE) }, delay = 0), silent = TRUE)
}



  observeEvent(input$run, {
    req(input$dm, input$sm, input$vm)
    dir.create(rv$outdir, showWarnings = FALSE, recursive = TRUE)
    rv$log <- ""; append_log("Starting post-run...")
    shinyWidgets::updateProgressBar(session, "pb", value = 0, total = t_total, title = "Starting", status = "info")
    output$step_status <- renderText("Ready."); output$time_status <- renderText(""); rv$t0 <- Sys.time()

    tryCatch({
      # [1/18] Read
      shinyWidgets::updateProgressBar(session, "pb", value = 1, total = t_total, title = sprintf("[1/%d] Reading tables", t_total), status = "info")
      output$step_status <- renderText(sprintf("[1/%d] Reading tables", t_total)); append_log("Reading the 3 tables")
      dm <- read_tabflex(input$dm$datapath); sm <- read_tabflex(input$sm$datapath); vm <- read_tabflex(input$vm$datapath)

      rid_col <- detect_col(dm, c("name","feature_id","id","row_id","compound","feature"))
      if (is.na(rid_col)) stop("dataMatrix must contain a row id column (e.g., 'name')")
      row_ids <- dm[[rid_col]]; dm <- dm %>% dplyr::select(-all_of(rid_col)); dm <- as.data.frame(dm); rownames(dm) <- row_ids

      sid_col <- detect_col(sm, c("sample_name","sample","id","Sample","Sample_name"))
      cls_col <- detect_col(sm, c("class","Class","group","Group"))
      if (is.na(sid_col) || is.na(cls_col)) stop("sampleMetadata must contain sample_name and class columns")

      common <- intersect(colnames(dm), as.character(sm[[sid_col]]))
      if (length(common) < 3) stop("Too few matched samples between dataMatrix and sampleMetadata")
      dm <- dm[, common, drop = FALSE]; sm <- sm[match(common, sm[[sid_col]]), , drop = FALSE]

      classes <- factor(as.character(sm[[cls_col]])); cls_levels <- levels(classes)
      W4M_palette <- function(n) scales::hue_pal()(n); CLRS <- setNames(W4M_palette(length(cls_levels)), cls_levels)

      Xraw <- as_matrix_numeric(dm); Xraw[!is.finite(Xraw) | Xraw < 0] <- NA_real_
      # ---- OVR enrichment (STRICT) injected patch ----
      # Requires: stats_patch_OVR_pairwise_u10g.R (already sourced at top)
      # Builds log2 matrix from Xraw, then enriches vm with q_global, q_posthoc_agg_OVR_<class>, log2FC_OVR_final_<class>
      try({
        DM_log2 <- suppressWarnings(log2(Xraw)); DM_log2[!is.finite(DM_log2)] <- NA_real_
        classes <- factor(as.character(sm[[cls_col]])); cls_levels <- levels(classes)
        vm <- enrich_VM_with_OVR(vm, DM = DM_log2, SM = sm, class_col = cls_col, classes = cls_levels,
                                 already_log = TRUE, mode = "STRICT")
        vm <- enrich_VM_with_OVR(vm, DM = DM_log2, SM = sm, class_col = cls_col, classes = cls_levels,
                                 already_log = TRUE, mode = "UNION")
        append_log("OVR STRICT+UNION enrichment done (KW + Dunn (within-feature BH) -> BH across).")
      }, silent = TRUE)
      # ---- end injected patch ----


      # [2/18] Prep multivariate
      shinyWidgets::updateProgressBar(session, "pb", value = 2, total = t_total, title = sprintf("[2/%d] Multivariate prep", t_total), status = "info")
      output$step_status <- renderText(sprintf("[2/%d] Multivariate prep", t_total)); append_log("Preparing log10 & scaling for multivariate")
      Xlog10 <- suppressWarnings(log10(Xraw)); Xlog10[!is.finite(Xlog10)] <- NA_real_
      Xmv <- if (isTRUE(input$do_log10)) Xlog10 else Xraw
      if (input$scaling == "pareto") Xmv <- pareto_scale(Xmv) else if (input$scaling == "unitVar") Xmv <- unitvar_scale(Xmv)
      if (input$ropls_scaleC != "none" && input$scaling != "none") append_log("WARNING: double scaling (pre-scaling + ropls::scaleC)")
      Xmv[!is.finite(Xmv)] <- 0; Xmv_t <- t(Xmv)

      # Predictive components inputs (support NA for PCA/PLS-DA)
      pred_pca <- if (identical(input$pca_predI, "NA")) NA else as.integer(input$pca_predI)
      if (!is.na(pred_pca) && pred_pca < 1) pred_pca <- 1

      pred_pls <- if (identical(input$pls_predI, "NA")) NA else as.integer(input$pls_predI)
      if (!is.na(pred_pls) && pred_pls < 1) pred_pls <- 1

      pred_opls <- as.integer(input$opls_predI)
      if (!is.finite(pred_opls) || pred_opls < 1) pred_opls <- 1
      if (pred_opls != 1) append_log(paste0("NOTE: OPLS-DA with ", pred_opls, " predictive components (guideline is 1)"))
    

      # CV segments sanity check
      cv_k <- as.integer(input$cv_segments)
      if (!is.finite(cv_k) || cv_k < 2) cv_k <- 2
      if (cv_k > nrow(Xmv_t)) { append_log(paste0("CV segments reduced from ", input$cv_segments, " to ", nrow(Xmv_t), " (samples)")); cv_k <- nrow(Xmv_t) }


      # [3/18] PCA
      shinyWidgets::updateProgressBar(session, "pb", value = 3, total = t_total, title = sprintf("[3/%d] PCA (ropls)", t_total), status = "primary")
      output$step_status <- renderText(sprintf("[3/%d] PCA (ropls)", t_total)); append_log("PCA (ropls) with 95% ellipses")
      set.seed(123); pca <- ropls::opls(Xmv_t, predI = pred_pca, scaleC = input$ropls_scaleC, crossvalI = cv_k, info.txtC = "none")
      oldpal <- palette(); palette(unname(CLRS))
      png(file.path(rv$outdir, "PCA_scores_ropls.png"), width = 1800, height = 1400, res = 200)
      par(mar = c(5,5,2,2)); plot(pca, typeVc = "x-score", parAsColFcVn = classes, parEllipsesL = TRUE); add_top_title("PCA (ropls)"); dev.off(); palette(oldpal)
      sink(file.path(rv$outdir, "PCA_model_ropls.txt")); print(pca); sink(NULL)
      export_scores_loadings(pca, "PCA", rv$outdir)

      # [4/18] PLS-DA
      shinyWidgets::updateProgressBar(session, "pb", value = 4, total = t_total, title = sprintf("[4/%d] PLS-DA (ropls)", t_total), status = "primary")
      output$step_status <- renderText(sprintf("[4/%d] PLS-DA (ropls)", t_total)); append_log("PLS-DA (ropls) with 95% ellipses + permutations")
      set.seed(123); y <- classes
      plsda <- ropls::opls(Xmv_t, y, predI = pred_pls, orthoI = 0, scaleC = input$ropls_scaleC, permI = input$permI, crossvalI = cv_k, info.txtC = "none")
      oldpal <- palette(); palette(unname(CLRS))
      png(file.path(rv$outdir, "PLSDA_scores_ropls.png"), width = 1800, height = 1400, res = 200)
      par(mar = c(5,5,2,2)); plot(plsda, typeVc = input$pls_plot_type, parAsColFcVn = classes, parEllipsesL = TRUE); add_top_title("PLS-DA (ropls)"); dev.off(); palette(oldpal)
      sink(file.path(rv$outdir, "PLSDA_model_ropls.txt")); print(plsda); sink(NULL)
      try(export_scores_loadings(plsda, "PLSDA", rv$outdir), silent = FALSE)
      rv$pred_plsda <- try(export_predictions(plsda, "PLSDA", rv$outdir, Xmv_t, classes, append_log), silent = FALSE)
      try(export_perm_intercepts(plsda, "PLSDA", rv$outdir, append_log), silent = FALSE)

      # [5/18] VIP (PLS-DA)
      shinyWidgets::updateProgressBar(session, "pb", value = 5, total = t_total, title = sprintf("[5/%d] VIP (PLS-DA)", t_total), status = "primary")
      output$step_status <- renderText(sprintf("[5/%d] VIP (PLS-DA)", t_total))
      name_col_vm <- detect_col(vm, c("name","feature_id","id")); if (is.na(name_col_vm)) name_col_vm <- "name"
      vip_pls <- tryCatch(slot(plsda, "vipVn"), error = function(e) NULL)
      options(vip_hm_ratio = input$vip_hm_ratio)
      if (!is.null(vip_pls) && length(vip_pls) > 0) {
        rv$vip_pls_feats <- vip_dot_heatmap(vip_pls, input$vip_mode, input$vip_topN, input$vip_thr, input$vip_cum, vm, name_col_vm, classes, Xraw,
                        file.path(rv$outdir, "VIP_PLSDA_dot_heatmap.png"),
                        title = "VIP (PLS-DA) - top features",
                        pal = input$div_palette)
      }

      # [6/18] PLS-DA permutations
      shinyWidgets::updateProgressBar(session, "pb", value = 6, total = t_total, title = sprintf("[6/%d] PLS-DA permutations", t_total), status = "primary")
      output$step_status <- renderText(sprintf("[6/%d] PLS-DA permutations", t_total))
      if (!is.null(input$permI) && input$permI > 0) {
        png(file.path(rv$outdir, "PLSDA_permutations_ropls.png"), width = 1800, height = 1400, res = 200)
        par(mar = c(5,5,2,2)); plot(plsda, typeVc = "permutation"); add_top_title("PLS-DA permutations"); dev.off()
      } else append_log("Permutations disabled (permI = 0)")

      # [7-10/18] OPLS-DA (binary) or OVR per class
      if (length(cls_levels) == 2) {
        shinyWidgets::updateProgressBar(session, "pb", value = 7, total = t_total, title = sprintf("[7/%d] OPLS-DA (ropls)", t_total), status = "primary")
        output$step_status <- renderText(sprintf("[7/%d] OPLS-DA (ropls)", t_total)); append_log("OPLS-DA (ropls) + permutations")
        set.seed(123)
        oplsda <- ropls::opls(Xmv_t, y, predI = pred_pls, orthoI = input$orthoI, scaleC = input$ropls_scaleC, permI = input$permI, crossvalI = cv_k, info.txtC = "none")
        oldpal <- palette(); palette(unname(CLRS))
        png(file.path(rv$outdir, "OPLSDA_scores_ropls.png"), width = 1800, height = 1400, res = 200)
        par(mar = c(5,5,2,2)); plot(oplsda, typeVc = input$pls_plot_type, parAsColFcVn = classes, parEllipsesL = TRUE); add_top_title("OPLS-DA (ropls)"); dev.off(); palette(oldpal)

      # S-plot for OPLS-DA (binary)
      try(splot_export(oplsda, Xmv, Xlog10, "OPLSDA", rv$outdir, topN = input$splot_topN, title = "S-plot (OPLS-DA)",
                      corr_thr = input$splot_corr_thr, cov_q = input$splot_cov_q,
                      color_by_vip = input$splot_color_vip, vip_thr = input$splot_vip_thr), silent = FALSE)
        sink(file.path(rv$outdir, "OPLSDA_model_ropls.txt")); print(oplsda); sink(NULL)
        try(export_scores_loadings(oplsda, "OPLSDA", rv$outdir), silent = FALSE)
        rv$pred_oplsda <- try(export_predictions(oplsda, "OPLSDA", rv$outdir, Xmv_t, classes, append_log), silent = FALSE)
        try(export_perm_intercepts(oplsda, "OPLSDA", rv$outdir, append_log), silent = FALSE)
        if (!is.null(input$permI) && input$permI > 0) {
          png(file.path(rv$outdir, "OPLSDA_permutations_ropls.png"), width = 1800, height = 1400, res = 200)
          par(mar = c(5,5,2,2)); plot(oplsda, typeVc = "permutation"); add_top_title("OPLS-DA permutations"); dev.off()
        }
        # VIP for OPLSDA
        vip_opls <- tryCatch(slot(oplsda, "vipVn"), error = function(e) NULL)
        options(vip_hm_ratio = input$vip_hm_ratio)
        if (!is.null(vip_opls) && length(vip_opls) > 0) {
          vip_dot_heatmap(vip_opls, input$vip_mode, input$vip_topN, input$vip_thr, input$vip_cum, vm, name_col_vm, classes, Xraw,
                          file.path(rv$outdir, "VIP_OPLSDA_dot_heatmap.png"),
                          title = "VIP (OPLS-DA) - top features",
                          pal = input$div_palette)
        }
      } else {
        i <- 7
        for (cl in cls_levels) {
          i <- i + 1; shinyWidgets::updateProgressBar(session, "pb", value = i, total = t_total, title = sprintf("[%d/%d] OPLS-DA OVR: %s", i, t_total, cl), status = "primary")
          output$step_status <- renderText(sprintf("[%d/%d] OPLS-DA OVR: %s", i, t_total, cl))
          y_bin <- factor(ifelse(classes == cl, cl, paste0("not_", cl)), levels = c(cl, paste0("not_", cl)))
          append_log(paste0("OPLS-DA OVR for ", cl))
          set.seed(123)
          opls_ovr <- ropls::opls(Xmv_t, y_bin, predI = pred_opls, orthoI = input$orthoI, scaleC = input$ropls_scaleC, permI = input$permI, crossvalI = cv_k, info.txtC = "none")
          oldpal <- palette(); palette(c(unname(CLRS[cl]), "#999999"))
          png(file.path(rv$outdir, paste0("OPLSDA_OVR_", stubify(cl), "_scores_ropls.png")), width = 1800, height = 1400, res = 200)
          par(mar = c(5,5,2,2)); plot(opls_ovr, typeVc = input$pls_plot_type, parAsColFcVn = y_bin, parEllipsesL = TRUE); add_top_title(paste0("OPLS-DA OVR - ", cl)); dev.off(); palette(oldpal)
          sink(file.path(rv$outdir, paste0("OPLSDA_OVR_", stubify(cl), "_model_ropls.txt"))); print(opls_ovr); sink(NULL)

# S-plot for OPLS-DA OVR (per class)
# Compute strict orientation from CURRENT run (no dependency on vm_enr):
lfc_vec <- tryCatch({
  tmp_mat <- log2fc_by_class(Xraw, classes)    # rows = features, cols = class levels
  v <- tmp_mat[, cl, drop = TRUE]
  names(v) <- rownames(tmp_mat)
  v
}, error = function(e) NULL)
dir.create(rv$outdir, recursive = TRUE, showWarnings = FALSE);
res_splot_ovr <- try({

  m <- splot_export(opls_ovr, Xmv, Xlog10, paste0("OPLSDA_OVR_", stubify(cl)), rv$outdir,

                    topN = input$splot_topN, title = paste0("S-plot (OPLS-DA OVR) - ", cl),

                    corr_thr = input$splot_corr_thr, cov_q = input$splot_cov_q,

                    color_by_vip = input$splot_color_vip, vip_thr = input$splot_vip_thr,

                    lfc_by_feature = lfc_vec)

  append_log("S-plot OK")

  m

}, silent = FALSE)

if (!inherits(res_splot_ovr, "try-error") && is.data.frame(res_splot_ovr)) {
  cand_tmp <- res_splot_ovr[which(res_splot_ovr$candidate), , drop = FALSE]
  if (nrow(cand_tmp) > 0) {
    cand_tmp$class_OVR <- cl
    rv$splot_ovr_union <- dplyr::bind_rows(rv$splot_ovr_union, cand_tmp)
  }
}

          
          # --- S-plot: filtres stricts (q_global & q_posthoc_agg_OVR & |log2FC_OVR_final|) ---
          if (!inherits(res_splot_ovr, "try-error") && isTRUE(input$splot_use_sig_filters)) {
            name_col_vm <- detect_col(vm, c("name","feature_id","id")); if (is.na(name_col_vm)) name_col_vm <- "name"
            # Construire noms de colonnes spécifiques à la classe
            qovr_col  <- paste0("q_posthoc_agg_OVR_", cl)
            l2fc_col  <- paste0("log2FC_OVR_final_", cl)
            # Sélection VM minimale pour joindre
            vm_join <- try(dplyr::select(vm, !!rlang::sym(name_col_vm), q_global, !!rlang::sym(qovr_col), !!rlang::sym(l2fc_col)), silent = FALSE)
            if (!inherits(vm_join, "try-error")) {
              vm_join <- dplyr::rename(vm_join,
                                       feature = !!rlang::sym(name_col_vm),
                                       q_posthoc = !!rlang::sym(qovr_col),
                                       log2FC_OVR_final = !!rlang::sym(l2fc_col))
              vm_join$feature <- as.character(vm_join$feature); res_splot_ovr$feature <- as.character(res_splot_ovr$feature);
              df_annot <- dplyr::left_join(res_splot_ovr, vm_join, by = "feature")
              df_annot$candidate_strict <- with(df_annot,
                candidate & is.finite(q_global) & is.finite(q_posthoc) & is.finite(log2FC_OVR_final) &
                (q_global <= input$splot_sig_q_global) &
                (q_posthoc <= input$splot_sig_q_ovr) &
                (abs(log2FC_OVR_final) >= input$splot_sig_abs_l2fc)
              )
              # Export annoté
              readr::write_tsv(df_annot, file.path(rv$outdir, paste0("OPLSDA_OVR_", stubify(cl), "_Splot_data_annotated.tsv")))
              # Export candidats stricts
              cand_strict <- df_annot[df_annot$candidate_strict & !is.na(df_annot$candidate_strict), , drop = FALSE]
              cand_strict <- cand_strict[order(cand_strict$abs_prod, decreasing = TRUE), , drop = FALSE]
              readr::write_tsv(cand_strict, file.path(rv$outdir, paste0("OPLSDA_OVR_", stubify(cl), "_Splot_candidates_STRICT.tsv")))
              # Stockage pour enrichissement VM ultérieur
              rv$splot_by_class[[cl]] <- df_annot
            } else {
              rv$splot_by_class[[cl]] <- res_splot_ovr
            }
          } else {
            # Stockage simple si pas de filtre strict
            rv$splot_by_class[[cl]] <- res_splot_ovr
          }
try(export_scores_loadings(opls_ovr, paste0("OPLSDA_OVR_", stubify(cl)), rv$outdir), silent = FALSE)
          # store predictions for SM enrichment
          rv$pred_ovr[[cl]] <- try(export_predictions(opls_ovr, paste0("OPLSDA_OVR_", stubify(cl)), rv$outdir, Xmv_t, y_bin, append_log), silent = FALSE)
          try(export_perm_intercepts(opls_ovr, paste0("OPLSDA_OVR_", stubify(cl)), rv$outdir, append_log), silent = FALSE)
          if (!is.null(input$permI) && input$permI > 0) {
            png(file.path(rv$outdir, paste0("OPLSDA_OVR_", stubify(cl), "_permutations_ropls.png")), width = 1800, height = 1400, res = 200)
            par(mar = c(5,5,2,2)); plot(opls_ovr, typeVc = "permutation"); add_top_title(paste0("OPLS-DA OVR permutations - ", cl)); dev.off()
          }
          vip_ovr <- tryCatch(slot(opls_ovr, "vipVn"), error = function(e) NULL)
          options(vip_hm_ratio = input$vip_hm_ratio)
          if (!is.null(vip_ovr) && length(vip_ovr) > 0) {
            rv$vip_ovr_feats[[cl]] <- vip_dot_heatmap(vip_ovr, input$vip_mode, input$vip_topN, input$vip_thr, input$vip_cum, vm, name_col_vm, classes, Xraw,
                            file.path(rv$outdir, paste0("VIP_OPLSDA_OVR_", stubify(cl), "_dot_heatmap.png")),
                            title = paste0("VIP (OPLS-DA OVR) - ", cl),
                            pal = input$div_palette)
          }
        }
      }

      # [11/18] UNIV_hit heatmap(s)
      shinyWidgets::updateProgressBar(session, "pb", value = 11, total = t_total, title = sprintf("[11/%d] UNIV_hit heatmaps", t_total), status = "primary")
      output$step_status <- renderText(sprintf("[11/%d] UNIV_hit heatmaps", t_total))
      univ_cols <- grep("^UNIV_hit", names(vm), value = TRUE)
      if (!("UNIV_hit" %in% names(vm)) && length(univ_cols) > 0) vm$UNIV_hit <- as.integer(rowSums(as.matrix(vm[, univ_cols, drop = FALSE]), na.rm = TRUE) >= 1L)
      name_col <- detect_col(vm, c("name","feature_id","id")); if (is.na(name_col)) name_col <- "name"
      hit_idx <- if ("UNIV_hit" %in% names(vm)) which(vm$UNIV_hit == 1L) else integer(0)

      select_top_red_union <- function(hit_idx, K, Xraw, classes, vm, name_col) {
        if (length(hit_idx) < 1 || K <= 0) return(character(0))
        l2_all <- log2fc_by_class(Xraw, classes);
          keys <- as.character(vm[[name_col]][hit_idx]);
          keys <- keys[keys %in% rownames(l2_all)];
          if (!length(keys)) return(character(0));
          l2 <- l2_all[keys, , drop = FALSE]
        keep <- c()
        for (cl in colnames(l2)) {
          ord <- order(l2[, cl], decreasing = TRUE)
          sel <- head(rownames(l2)[ord], min(K, nrow(l2)))
          keep <- unique(c(keep, sel))
        }
        keep[keep %in% rownames(Xraw)]
      }

      if (length(hit_idx) >= 2) {
        feats_union <- character(0)
        if (!is.null(input$hm_univ_top_per_class) && input$hm_univ_top_per_class > 0) feats_union <- select_top_red_union(hit_idx, input$hm_univ_top_per_class, Xraw, classes, vm, name_col)

        Mhit_sel <- if (length(feats_union)) Xraw[feats_union, , drop = FALSE] else {
          Mhit0 <- Xraw[hit_idx, , drop = FALSE]; vars <- apply(Mhit0, 1, stats::var, na.rm = TRUE)
          ord <- order(vars, decreasing = TRUE)[seq_len(min(length(vars), input$hm_univ_top))]; Mhit0[ord, , drop = FALSE]
        }

        Mz <- if (isTRUE(input$hm_univ_scale)) t(scale(t(Mhit_sel))) else Mhit_sel
        Mz[!is.finite(Mz)] <- 0

        # Column annotations: build 'Class' + optional user-selected, deduplicate
        ann_col <- data.frame(Class = as.character(classes), stringsAsFactors = FALSE); rownames(ann_col) <- colnames(Mz)
        if (!is.null(input$hm_ann_cols) && length(input$hm_ann_cols) > 0) {
          add_cols <- intersect(input$hm_ann_cols, names(sm))
          add_cols <- setdiff(add_cols, c("class","Class","group","Group"))
          if (length(add_cols) > 0) ann_col <- cbind(ann_col, sm[match(rownames(ann_col), sm[[sid_col]]), add_cols, drop = FALSE])
        }
        keep <- !duplicated(tolower(colnames(ann_col))); ann_col <- ann_col[, keep, drop = FALSE]

        pal_hm <- get_div_palette(input$hm_palette, 100)
        ph_z <- pheatmap::pheatmap(Mz, annotation_col = ann_col,
                           fontsize_row = input$hm_univ_font_row,
                           fontsize_col = input$hm_univ_font_col,
                           filename = NA,
                           main = if (length(feats_union)) sprintf("UNIV_hit z-score - top %d per class (most up-regulated)", input$hm_univ_top_per_class) else "UNIV_hit z-score (top variance)",
                           width = getOption("FIG_UNIVHM_W", input$hm_univ_w * .get_export_dpi()) / .get_export_dpi(), height = getOption("FIG_UNIVHM_H", input$hm_univ_h * .get_export_dpi()) / .get_export_dpi(), silent = FALSE,
                           cluster_rows = isTRUE(input$hm_univ_cluster_rows),
                           cluster_cols = isTRUE(input$hm_univ_cluster_cols),
                           color = pal_hm)
save_pheatmap_px(ph_z, file.path(rv$outdir, if (length(feats_union)) "Heatmap_UNIV_hit_TopKPerClass.png" else "Heatmap_UNIV_hit.png"), key = "UNIVHM")


        # Log2FC heatmap
        l2fc_hit <- log2fc_by_class(Xraw, classes)[rownames(Mhit_sel), , drop = FALSE]
        pal_hm2 <- get_div_palette(input$hm_palette, 100)
        ph_fc <- pheatmap::pheatmap(l2fc_hit,
                           fontsize_row = input$hm_univ_font_row,
                           fontsize_col = input$hm_univ_font_col,
                           filename = NA,
                           main = if (length(feats_union)) sprintf("UNIV_hit log2FC - top %d per class", input$hm_univ_top_per_class) else "UNIV_hit log2FC (top variance)",
                           width = getOption("FIG_UNIVHM_W", input$hm_univ_w * .get_export_dpi()) / .get_export_dpi(), height = getOption("FIG_UNIVHM_H", input$hm_univ_h * .get_export_dpi()) / .get_export_dpi(), silent = FALSE,
                           cluster_rows = isTRUE(input$hm_univ_cluster_rows),
                           cluster_cols = TRUE,
                           color = pal_hm2)
save_pheatmap_px(ph_fc, file.path(rv$outdir, if (length(feats_union)) "Heatmap_UNIV_hit_TopKPerClass_Log2FC.png" else "Heatmap_UNIV_hit_Log2FC.png"), key = "UNIVHM")

      } else append_log("No UNIV_hit detected - heatmap skipped")

      # [12/18] Venn UNIV_hit
      shinyWidgets::updateProgressBar(session, "pb", value = 12, total = t_total, title = sprintf("[12/%d] Venn UNIV_hit", t_total), status = "primary")
      output$step_status <- renderText(sprintf("[12/%d] Venn UNIV_hit", t_total))
      class_levels <- cls_levels; sets <- list()
      for (cl in class_levels) {
        coln <- paste0("UNIV_hit_", cl)
        if (coln %in% names(vm)) sets[[cl]] <- as.character(vm[[name_col]][which(vm[[coln]] == 1L)])
      }
      if (length(sets) >= 2) venn_from_lists(sets, file_png = file.path(rv$outdir, "Venn_UNIV_hit.png"), title = "Venn - UNIV_hit by class")
      else append_log("Not enough classes with UNIV_hit_* to draw a Venn")

      
      

# [13/18] OVR volcanoes
shinyWidgets::updateProgressBar(session, "pb", value = 13, total = t_total,
  title = sprintf("[13/%d] OVR volcanoes", t_total), status = "primary")
output$step_status <- renderText(sprintf("[13/%d] OVR volcanoes", t_total))
if (isTRUE(input$do_ovr_volcano)) {
  gate_q_global <- 0.05  # rigor gate (applied only if q_global exists)
  name_col <- detect_col(vm, c("name","feature_id","id"))
  if (is.na(name_col)) name_col <- "name"
  for (cl in cls_levels) {
    lfc_col <- paste0("log2FC_OVR_final_", cl)
    q_col   <- paste0("q_posthoc_agg_OVR_", cl)
    if (!all(c(lfc_col, q_col) %in% names(vm))) {
      append_log(paste0("OVR volcano (", cl, "): missing ", lfc_col, " or ", q_col, " in VM")); next
    }
    dfv <- vm %>% dplyr::select(all_of(name_col), all_of(lfc_col), all_of(q_col), dplyr::any_of("q_global"))
    names(dfv)[1] <- "name"
    if ("q_global" %in% names(dfv)) dfv <- dfv %>% dplyr::filter(q_global <= gate_q_global)
    if (nrow(dfv) == 0) { append_log(paste0("OVR volcano (", cl, "): no points pass q_global <= ", gate_q_global)); next }
    names(dfv)[names(dfv) == lfc_col] <- "log2FC"
    names(dfv)[names(dfv) == q_col]   <- "q"
    p <- make_volcano(dfv, "log2FC", "q",
          paste0("Volcano OVR (q_global≤", gate_q_global, "): ", cl),
          lfc_thr = input$lfc_thr, q_thr = input$q_thr,
          top_label = input$volcano_top_label, pal = input$div_palette)
    safe_write_gg(p,
      file.path(rv$outdir, paste0("Volcano_OVR_", stubify(cl), ".png")),
      width = 7, height = 7)
  }
}
# [14/18] Pairwise volcanoes
      shinyWidgets::updateProgressBar(session, "pb", value = 14, total = t_total, title = sprintf("[14/%d] Pairwise volcanoes", t_total), status = "primary")
      output$step_status <- renderText(sprintf("[14/%d] Pairwise volcanoes", t_total))
      if (isTRUE(input$do_pairwise_volcano)) {
        Xv <- if (isTRUE(input$do_log10)) Xlog10 else log10(pmax(Xraw, .Machine$double.eps))
        pairs <- combn(cls_levels, 2, simplify = FALSE)
        for (pr in pairs) {
          g1 <- pr[1]; g2 <- pr[2]
          stats_df <- pairwise_volcano_stats(Xv, as.character(classes), g1, g2)
          stats_df$q <- p.adjust(stats_df$p, method = "BH")
          p <- make_volcano(stats_df, "log2FC", "q", paste0("Volcano: ", g1, " vs ", g2),
                            lfc_thr = input$lfc_thr, q_thr = input$q_thr, top_label = input$volcano_top_label, pal = input$div_palette)
          safe_write_gg(p, file.path(rv$outdir, paste0("Volcano_", stubify(g1), "_vs_", stubify(g2), ".png")), width = 7, height = 5, dpi = 150)
          readr::write_tsv(stats_df, file.path(rv$outdir, paste0("VolcanoStats_", stubify(g1), "_vs_", stubify(g2), ".tsv")))
        }
      }
# [15/18] Biosigner (+ per-feature plots incl. OVR)
      shinyWidgets::updateProgressBar(session, "pb", value = 15, total = t_total, title = sprintf("[15/%d] Biosigner", t_total), status = "primary")
      output$step_status <- renderText(sprintf("[15/%d] Biosigner", t_total))

      biosign_safe <- function(Xdf, yfac, methods) {
        attempt <- function(m) tryCatch(biosigner::biosign(Xdf, y = yfac, bootI = input$biosigner_boot, method = m), error = function(e) e)
        mth <- methods; res <- attempt(mth)
        if (inherits(res, "error") && grepl("not available", conditionMessage(res))) {
          msg <- conditionMessage(res)
          if (grepl("\\brf\\b", msg)) mth <- setdiff(mth, "rf")
          if (grepl("\\bsvm\\b", msg)) mth <- setdiff(mth, "svm")
          if (grepl("\\bplsda\\b", msg)) mth <- setdiff(mth, "plsda")
          if (length(mth) > 0) res <- attempt(mth)
        }
        if (inherits(res, "error")) {
          for (m in methods) { r2 <- attempt(m); if (!inherits(r2, "error")) return(r2) }
        }
        res
      }

      mth <- c(); if (requireNamespace("pls", quietly = TRUE)) mth <- c(mth, "plsda") else append_log("Biosigner: package 'pls' missing - plsda disabled")
      if (requireNamespace("randomForest", quietly = TRUE)) mth <- c(mth, "randomforest") else append_log("Biosigner: package 'randomForest' missing - rf disabled")
      if (requireNamespace("e1071", quietly = TRUE)) mth <- c(mth, "svm") else append_log("Biosigner: package 'e1071' missing - svm disabled")
      if (!length(mth)) { append_log("Biosigner: no method available - skipped") } else {
        set.seed(123)
        Xbio_df <- as.data.frame(Xmv_t, stringsAsFactors = FALSE)
        if (length(cls_levels) == 2) {
          append_log(paste0("Biosigner binary (", paste(mth, collapse = ","), ")"))
          ybio <- factor(as.character(classes))
          bres <- biosign_safe(Xbio_df, ybio, mth)
          sink(file.path(rv$outdir, "Biosigner_model.txt")); print(bres); sink(NULL)
          ## [u10g] concordant export
          try(biosigner_export(file.path(rv$outdir, "Biosigner_model.txt"), file.path(rv$outdir, "Biosigner")), silent = TRUE)
          ## [u10f] concordant export
          try(biosigner_export(file.path(rv$outdir, "Biosigner_model.txt"), file.path(rv$outdir, "Biosigner")), silent = TRUE)
          ## [u10f] old PNG disabled: ## [u10g] old biosigner PNG disabled
          sign_levels <- if (identical(input$biosigner_sigset, "S+A")) c("S","AS") else c("S")
          sign_list <- extract_biosigner_signatures(bres, sign_levels)
          if (!is.null(sign_list) && length(sign_list)) {
            for (m in names(sign_list)) readr::write_tsv(tibble::tibble(name = sign_list[[m]]), file.path(rv$outdir, paste0("Biosigner_signatures_", m, ".tsv")))
            dfb <- tibble::tibble(method = names(sign_list), n = vapply(sign_list, length, integer(1)))
            p_b <- ggplot(dfb, aes(x = method, y = n)) + geom_col() + labs(title = "Biosigner - number of signature features")
            safe_write_gg(p_b, file.path(rv$outdir, "Biosigner_counts.png"), width = 6, height = 4, dpi = 150)

            feat_union <- unique(unlist(sign_list))
            if (length(feat_union)) {
              feat_union <- head(feat_union, input$biosigner_plot_n)
              plot_feature_distributions(Xraw, classes, feat_union,
                                         file.path(rv$outdir, paste0("Biosigner_Signatures_", input$biosigner_plot_type, "plots.png")),
                                         title = "Biosigner - signature features", type = input$biosigner_plot_type)

                # Also create per-method plots (RF / PLSDA / SVM) if desired
                for (m in names(sign_list)) {
                  feats_m <- head(sign_list[[m]], input$biosigner_plot_n)
                  if (length(feats_m)) {
                    plot_feature_distributions(
                      Xraw, classes, feats_m,
                      file.path(rv$outdir, paste0(prefix, "_Signatures_", toupper(m), "_", input$biosigner_plot_type, "plots.png")),
                      title = paste0("Biosigner OVR - ", cl, " - ", toupper(m)),
                      type = input$biosigner_plot_type
                    )
                  }
                }
            } else append_log("Biosigner: no signature features to plot")
          } else append_log("Biosigner: no signatures extracted")
        } else {
          append_log(paste0("Biosigner OVR (methods: ", paste(mth, collapse = ","), ")"))
          all_counts <- list(); rv$biosigner_sigs <- list()
          for (cl in cls_levels) {
            y_bin <- factor(ifelse(classes == cl, cl, paste0("not_", cl)), levels = c(cl, paste0("not_", cl)))
            prefix <- paste0("Biosigner_OVR_internal_", stubify(cl))
            append_log(paste0("Biosigner OVR for ", cl))
            bres <- biosign_safe(Xbio_df, y_bin, mth)
            if (inherits(bres, "error")) { append_log(paste0("Biosigner OVR ", cl, " - ERROR: ", conditionMessage(bres))); bres <- NULL }
            if (is.null(bres)) next
            sink(file.path(rv$outdir, paste0(prefix, "_model.txt"))); print(bres); sink(NULL)
            ## [u10g] concordant export
            try(biosigner_export(file.path(rv$outdir, paste0(prefix, "_model.txt")), file.path(rv$outdir, prefix)), silent = TRUE)
            ## [u10f] concordant export
            try(biosigner_export(file.path(rv$outdir, paste0(prefix, "_model.txt")), file.path(rv$outdir, prefix)), silent = TRUE)
            ## [u10f] old PNG disabled; using biosigner_export below
# ## [u10f] old PNG disabled: ## [u10g] old biosigner PNG disabled
            sign_levels <- if (identical(input$biosigner_sigset, "S+A")) c("S","AS") else c("S")
            sign_list <- extract_biosigner_signatures(bres, sign_levels)
            if (!is.null(sign_list) && length(sign_list)) {
              # persist signatures for enrichment
              rv$biosigner_sigs[[cl]] <- sign_list
              for (m in names(sign_list)) readr::write_tsv(tibble::tibble(name = sign_list[[m]]), file.path(rv$outdir, paste0(prefix, "_signatures_", m, ".tsv")))
              all_counts[[cl]] <- tibble::tibble(class = cl, method = names(sign_list), n = vapply(sign_list, length, integer(1)))

              feat_union <- unique(unlist(sign_list))
              if (length(feat_union)) {
                feat_union <- head(feat_union, input$biosigner_plot_n)
                plot_feature_distributions(Xraw, classes, feat_union,
                                           file.path(rv$outdir, paste0(prefix, "_Signatures_", input$biosigner_plot_type, "plots.png")),
                                           title = paste0("Biosigner OVR - ", cl, " (signatures)"), type = input$biosigner_plot_type)
              } else append_log(paste0("Biosigner OVR ", cl, ": no signature features to plot"))
            } else append_log(paste0("Biosigner OVR ", cl, ": no signatures extracted"))
          }
          if (length(all_counts)) {
            valid <- purrr::compact(purrr::keep(all_counts, ~ is.data.frame(.x) && all(c("class","method","n") %in% names(.x)) && nrow(.x) > 0))
            if (length(valid)) {
              dfb <- dplyr::bind_rows(valid)
              if (nrow(dfb) > 0 && all(c("class","method","n") %in% names(dfb))) {
                p_b <- ggplot(dfb, aes(x = class, y = n, fill = method)) +
                  geom_col(position = "dodge") + labs(title = "Biosigner OVR - number of signature features")
                safe_write_gg(p_b, file.path(rv$outdir, "Biosigner_OVR_internal_counts.png"), width = 8, height = 4, dpi = 150)
                readr::write_tsv(dfb, file.path(rv$outdir, "Biosigner_OVR_internal_counts.tsv"))
              } else append_log("Biosigner OVR - no counts to export")
            } else append_log("Biosigner OVR - no signature results")
          }
        }
      }

      # [16/18] HTML report
      shinyWidgets::updateProgressBar(session, "pb", value = 16, total = t_total, title = sprintf("[16/%d] HTML report", t_total), status = "primary")
      output$step_status <- renderText(sprintf("[16/%d] HTML report", t_total))
      params <- list(
        log10 = input$do_log10, scaling = input$scaling, ropls_scaleC = input$ropls_scaleC,
        plot_type = input$pls_plot_type, pca_comp = input$pca_comp, pls_comp = input$pls_comp,
        orthoI = input$orthoI, permI = input$permI, q_thr = input$q_thr, lfc_thr = input$lfc_thr,
        biosigner_boot = input$biosigner_boot, biosigner_plot_type = input$biosigner_plot_type, biosigner_plot_n = input$biosigner_plot_n,
        vip_mode = input$vip_mode, vip_topN = input$vip_topN, vip_thr = input$vip_thr, vip_cum = input$vip_cum,
        volcano_top_label = input$volcano_top_label, div_palette = input$div_palette,
        hm_univ_top = input$hm_univ_top, hm_univ_top_per_class = input$hm_univ_top_per_class, hm_univ_scale = input$hm_univ_scale,
        hm_univ_cluster_rows = input$hm_univ_cluster_rows, hm_univ_cluster_cols = input$hm_univ_cluster_cols,
        hm_univ_w = input$hm_univ_w, hm_univ_h = input$hm_univ_h, hm_palette = input$hm_palette,
        hm_ann_cols = paste(input$hm_ann_cols, collapse = ",")
      )
      write_html_report(rv$outdir, params)

      # [17/18] Export enriched tables (dm/sm/vm)
      shinyWidgets::updateProgressBar(session, "pb", value = 17, total = t_total, title = sprintf("[17/%d] Export enriched tables", t_total), status = "primary")
      output$step_status <- renderText(sprintf("[17/%d] Export enriched tables", t_total)); append_log("Exporting enriched DM/SM/VM tables")
      enrich_and_export_tables(
        Xraw = Xraw, classes = classes, vm = vm, sm = sm,
        rid_col = rid_col, sid_col = sid_col, outdir = rv$outdir,
        vip_pls_feats = rv$vip_pls_feats, vip_ovr_feats = rv$vip_ovr_feats,
        biosigner_sigs = rv$biosigner_sigs,
        pred_plsda = if (is.list(rv$pred_plsda)) rv$pred_plsda[[1]] else rv$pred_plsda,
        pred_ovr = lapply(rv$pred_ovr, function(x) if (is.list(x)) x[[1]] else x),
        splot_by_class = rv$splot_by_class
      )

      # [18/18] Done
      shinyWidgets::updateProgressBar(session, "pb", value = t_total, total = t_total, title = sprintf("[%d/%d] Done", t_total, t_total), status = "success")
      output$step_status <- renderText(sprintf("[%d/%d] Done", t_total, t_total))
      elapsed <- as.numeric(difftime(Sys.time(), rv$t0, units = "secs"))
      output$time_status <- renderText(sprintf("Elapsed time: %0.0fs", elapsed))
      append_log("Finished.")
    }, error = function(e) {
      shinyWidgets::updateProgressBar(session, "pb", value = 0, total = t_total, title = "Error", status = "danger")
      append_log(paste("ERROR:", conditionMessage(e)))
    })
  })
}

shinyApp(ui, server)
