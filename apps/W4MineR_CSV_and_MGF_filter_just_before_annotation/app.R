
# CSV and MGF filter (just before annotation)
# Author: ChatGPT for Jason
# Date: 2025-10-22
# Patch: English-only UI + revised Flow sentence

options(shiny.maxRequestSize = 8 * 1024^3)  # 8 GB

suppressPackageStartupMessages({
  library(shiny)
  library(readr)
  library(dplyr)
  library(tools)
  if (!requireNamespace("shinyFiles", quietly = TRUE)) {
    stop("Package 'shinyFiles' is required. Please install it with install.packages('shinyFiles').")
  }
})

# ---- Helpers ----
norm_id <- function(x, force_numeric = TRUE) {
  x <- as.character(x); x <- trimws(x)
  if (force_numeric) {
    y <- gsub("\\D+", "", x)
    y[nchar(y) == 0] <- x[nchar(y) == 0]
    return(y)
  } else x
}

detect_id_column <- function(df) {
  nms <- names(df)
  hits <- which(
    grepl("^row[ _]?id$", nms, ignore.case = TRUE) |
    grepl("^feature[ _]?id$", nms, ignore.case = TRUE) |
    tolower(nms) %in% c("row id","rowid","feature id","feature_id","id","name")
  )
  if (length(hits) >= 1) return(nms[hits[1]])
  names(df)[1]
}

read_feature_set_from_csv <- function(path, force_numeric = TRUE) {
  mz <- suppressMessages(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  idcol <- detect_id_column(mz)
  ids <- norm_id(mz[[idcol]], force_numeric = force_numeric)
  unique(ids)
}

mgf_split_blocks <- function(lines) {
  L <- trimws(lines)
  begins <- which(L == "BEGIN IONS"); ends <- which(L == "END IONS")
  if (length(begins) != length(ends)) stop(sprintf("MGF parse error: %d BEGIN vs %d END.", length(begins), length(ends)))
  if (!length(begins)) return(list())
  blocks <- vector("list", length(begins))
  e_idx <- 1
  for (i in seq_along(begins)) {
    b <- begins[i]
    while (e_idx <= length(ends) && ends[e_idx] < b) e_idx <- e_idx + 1
    if (e_idx > length(ends)) stop("MGF parse error: unmatched BEGIN IONS.")
    e <- ends[e_idx]
    blocks[[i]] <- lines[b:e]
    e_idx <- e_idx + 1
  }
  blocks
}

mgf_block_feature_id <- function(block_lines, force_numeric = TRUE) {
  idx <- grep("^[Ff][Ee][Aa][Tt][Uu][Rr][Ee]_?[Ii][Dd]\\s*=", block_lines)
  if (!length(idx)) return(NA_character_)
  ln <- block_lines[idx[1]]
  val <- sub("^[Ff][Ee][Aa][Tt][Uu][Rr][Ee]_?[Ii][Dd]\\s*=\\s*", "", ln)
  norm_id(val, force_numeric = force_numeric)
}

mgf_filter_by_ids <- function(mgf_path, keep_ids, out_path, force_numeric = TRUE) {
  lines <- readLines(mgf_path, warn = FALSE)
  blocks <- mgf_split_blocks(lines)
  if (!length(blocks)) stop("No MGF blocks found.")
  ids <- vapply(blocks, mgf_block_feature_id, FUN.VALUE = character(1), force_numeric = force_numeric)
  keep <- ids %in% keep_ids
  kept_blocks <- blocks[keep]
  if (!length(kept_blocks)) {
    writeLines(character(0), con = out_path)
  } else {
    out_lines <- unlist(lapply(kept_blocks, function(bl) c(bl, "")), use.names = FALSE)
    writeLines(out_lines, con = out_path)
  }
  list(total_blocks = length(blocks), kept_blocks = sum(keep),
       dropped_blocks = sum(!keep),
       blocks_missing_feature_id = sum(is.na(ids) | ids == ""))
}

mgf_read_feature_ids <- function(mgf_path, force_numeric = TRUE) {
  lines <- readLines(mgf_path, warn = FALSE)
  blocks <- mgf_split_blocks(lines)
  if (!length(blocks)) return(character(0))
  ids <- vapply(blocks, mgf_block_feature_id, FUN.VALUE = character(1), force_numeric = force_numeric)
  ids <- ids[!is.na(ids) & ids != ""]
  unique(ids)
}

csv_filter_by_ids <- function(csv_path, out_path, keep_ids, force_numeric = TRUE) {
  df <- suppressMessages(readr::read_csv(csv_path, show_col_types = FALSE, progress = FALSE))
  idcol <- detect_id_column(df)
  orig_n <- nrow(df)
  df[[idcol]] <- norm_id(df[[idcol]], force_numeric = force_numeric)
  df2 <- dplyr::filter(df, .data[[idcol]] %in% keep_ids)
  readr::write_csv(df2, out_path)
  list(total_rows = orig_n, kept_rows = nrow(df2), id_column = idcol)
}

# ---- UI ----
ui <- fluidPage(
  titlePanel("CSV and MGF filter (just before annotation)"),
  sidebarLayout(
    sidebarPanel(
      tags$h4("Inputs"),
      fileInput("mzmine_csv", "CSV inputed after statistical analysis", accept = c(".csv")),
      fileInput("w4m_csv", "CSV not inputed before statistical analysis", accept = c(".csv")),
      fileInput("mgf_file", "MGF (e.g., filtered_STRICT_numeric.mgf)", accept = c(".mgf")),
      checkboxInput("force_numeric", 'Force numeric IDs (e.g., "X64" → "64")', value = TRUE),
      tags$hr(),
      tags$h4("Outputs (single destination folder)"),
      shinyFiles::shinyDirButton("outdir", "Browse... (output folder)", "Choose a folder"),
      verbatimTextOutput("outdir_text"),
      div(style="margin-top:10px"),
      textInput("w4m_name", "CSV not inputed before statistical analysis_FILTRATED",
                value = "CSV_not_inputed_before_statistical_analysis_FILTRATED"),
      textInput("mzmine_name", "CSV inputed after statistical analysis_FILTRATED",
                value = "CSV_inputed_after_statistical_analysis_FILTRATED"),
      textInput("mgf_name", "MGF final", value = "mgf_final.mgf"),
      verbatimTextOutput("paths_preview"),
      tags$hr(),
      checkboxInput("refilter_after_mgf", "After MGF filtering: restrict both CSVs to the features actually present in the filtered MGF (overwrite CSVs)", value = TRUE),
      tags$hr(),
      actionButton("run", "Run filtering", class = "btn btn-primary")
    ),
    mainPanel(
      tags$h4("Log"),
      verbatimTextOutput("log"),
      tags$hr(),
      tags$h4("Summary"),
      tableOutput("summary_tbl"),
      tags$hr(),
      tags$p(em("Flow: 1) IDs from 'CSV inputed after statistical analysis' → filter CSV not inputed before statistical analysis & MGF. 2) (optional) re-read the filtered MGF and restrict both CSVs to its IDs."))
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  volumes <- tryCatch(shinyFiles::getVolumes()(), error = function(e) NULL)
  if (is.null(volumes)) {
    volumes <- c(Home = "~", Root = "/")
    if (.Platform$OS.type == "windows") volumes <- c(volumes, "C:" = "C:\\", "D:" = "D:\\")
  }
  shinyFiles::shinyDirChoose(input, "outdir", roots = volumes, session = session, allowDirCreate = TRUE, defaultRoot = names(volumes)[1])

  out_dir <- reactiveVal(NULL)
  observeEvent(input$outdir, {
    parsed <- shinyFiles::parseDirPath(volumes, input$outdir)
    if (length(parsed)) out_dir(as.character(parsed))
  })

  output$outdir_text <- renderText({
    od <- out_dir()
    if (is.null(od)) "No folder selected."
    else paste0("Output folder -> ", od)
  })

  output$paths_preview <- renderText({
    od <- out_dir()
    if (is.null(od)) return("Path preview: (select a folder)")
    f <- function(nm, ext) {
      if (!nzchar(nm)) nm <- paste0("file.", ext)
      if (!grepl(paste0("\\.", ext, "$"), tolower(nm))) nm <- paste0(nm, ".", ext)
      file.path(od, nm)
    }
    paste0("W4M CSV -> ", f(input$w4m_name, "csv"), "\n",
           "MZmine CSV -> ", f(input$mzmine_name, "csv"), "\n",
           "MGF -> ", f(input$mgf_name, "mgf"))
  })

  log_lines <- reactiveVal(character(0))
  append_log <- function(...) {
    txt <- paste(...)
    log_lines(c(log_lines(), paste(format(Sys.time(), "%H:%M:%S"), "-", txt)))
  }
  output$log <- renderText(paste(log_lines(), collapse = "\n"))

  summary_state <- reactiveVal(
    data.frame(Item = character(0), Before = integer(0), After = integer(0), stringsAsFactors = FALSE)
  )
  output$summary_tbl <- renderTable(summary_state(), striped = TRUE, bordered = TRUE, spacing = "m", align = "lrr")

  observeEvent(input$run, {
    req(input$mzmine_csv, input$w4m_csv, input$mgf_file)
    force_num <- isTRUE(input$force_numeric)
    log_lines(character(0))
    summary_items <- list()

    tryCatch({
      od <- out_dir()
      if (is.null(od) || !nzchar(od)) stop("No output folder. Click 'Browse... (output folder)'.")

      fpath <- function(nm, ext) {
        if (!nzchar(nm)) nm <- paste0("file.", ext)
        if (!grepl(paste0("\\.", ext, "$"), tolower(nm))) nm <- paste0(nm, ".", ext)
        file.path(od, nm)
      }
      w4m_out    <- fpath(input$w4m_name, "csv")
      mzmine_out <- fpath(input$mzmine_name, "csv")
      mgf_out    <- fpath(input$mgf_name, "mgf")

      append_log("Reading IDs from 'CSV inputed after statistical analysis'...")
      keep_ids <- read_feature_set_from_csv(input$mzmine_csv$datapath, force_numeric = force_num)
      append_log(sprintf("IDs to keep (from that CSV): %d", length(keep_ids)))

      append_log("Filtering W4M CSV (pass 1, by those IDs)...")
      w4m_res <- csv_filter_by_ids(input$w4m_csv$datapath, w4m_out, keep_ids, force_numeric = force_num)
      append_log(sprintf("W4M CSV: %d -> %d rows kept (ID column: %s).",
                         w4m_res$total_rows, w4m_res$kept_rows, w4m_res$id_column))
      summary_items[["W4M CSV (pass1)"]] <- c(w4m_res$total_rows, w4m_res$kept_rows)

      append_log("Writing a filtered copy of 'CSV inputed after statistical analysis' (pass 1)...")
      mzmine_res1 <- csv_filter_by_ids(input$mzmine_csv$datapath, mzmine_out, keep_ids, force_numeric = force_num)
      append_log(sprintf("'CSV inputed after' : %d -> %d rows kept.",
                         mzmine_res1$total_rows, mzmine_res1$kept_rows))
      summary_items[["CSV inputed after (pass1)"]] <- c(mzmine_res1$total_rows, mzmine_res1$kept_rows)

      append_log("Filtering MGF by those IDs...")
      mgf_res <- mgf_filter_by_ids(input$mgf_file$datapath, keep_ids, mgf_out, force_numeric = force_num)
      append_log(sprintf("MGF: %d -> %d blocks kept (%d without FEATURE_ID, %d dropped).",
                         mgf_res$total_blocks, mgf_res$kept_blocks,
                         mgf_res$blocks_missing_feature_id, mgf_res$dropped_blocks))
      summary_items[["MGF (pass1, blocks)"]] <- c(mgf_res$total_blocks, mgf_res$kept_blocks)

      if (isTRUE(input$refilter_after_mgf)) {
        append_log("Option ON: read filtered MGF and restrict both CSVs to IDs actually present in the MGF...")
        mgf_ids <- mgf_read_feature_ids(mgf_out, force_numeric = force_num)
        append_log(sprintf("IDs present in filtered MGF: %d", length(mgf_ids)))

        append_log("Re-filter W4M CSV (pass 2 -> overwrite)...")
        w4m_res2 <- csv_filter_by_ids(input$w4m_csv$datapath, w4m_out, mgf_ids, force_numeric = force_num)
        append_log(sprintf("W4M CSV (pass2 via MGF): %d rows.", w4m_res2$kept_rows))
        summary_items[["W4M CSV (pass2 via MGF)"]] <- c(NA_integer_, w4m_res2$kept_rows)

        append_log("Re-filter 'CSV inputed after' (pass 2 -> overwrite)...")
        mzmine_res2 <- csv_filter_by_ids(input$mzmine_csv$datapath, mzmine_out, mgf_ids, force_numeric = force_num)
        append_log(sprintf("'CSV inputed after' (pass2 via MGF): %d rows.", mzmine_res2$kept_rows))
        summary_items[["CSV inputed after (pass2 via MGF)"]] <- c(NA_integer_, mzmine_res2$kept_rows)
      }

      if (length(summary_items)) {
        items <- names(summary_items)
        Before <- vapply(summary_items, function(v) ifelse(length(v)>=1, v[1], NA_integer_), numeric(1))
        After  <- vapply(summary_items, function(v) ifelse(length(v)>=2, v[2], NA_integer_), numeric(1))
        sm <- data.frame(Item = items, Before = Before, After = After, check.names = FALSE)
        summary_state(sm)
      }

      append_log(sprintf("Done. Files written to %s", od))
    }, error = function(e) {
      append_log(paste("ERROR:", conditionMessage(e)))
      showNotification(conditionMessage(e), type = "error", duration = NULL)
    })
  })
}

shinyApp(ui, server)
