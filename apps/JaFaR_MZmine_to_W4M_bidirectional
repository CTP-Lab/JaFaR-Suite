# =====================================================================
# JaFaR_MZmine_to_W4M_bidirectional
# =====================================================================


options(shiny.maxRequestSize = 1024 * 1024^2) # ~1 GB

suppressPackageStartupMessages({
  library(shiny)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(data.table)
  suppressWarnings(try(library(readxl), silent = TRUE))
  suppressWarnings(try(library(zip),    silent = TRUE))
  suppressWarnings(try(library(DT),     silent = TRUE))
})

# ---------- Helpers I/O ----------
read_any_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv")) {
    readr::read_csv(path, show_col_types = FALSE)
  } else if (ext %in% c("tsv", "tab", "tabular", "txt")) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else if (ext %in% c("xlsx")) {
    if (requireNamespace("readxl", quietly = TRUE)) {
      readxl::read_xlsx(path)
    } else {
      stop("XLSX file requires the 'readxl' package. Install it or provide TSV/CSV.")
    }
  } else {
    readr::read_delim(path, delim = "\t", show_col_types = FALSE, progress = FALSE)
  }
}

# Normalisation : lower + trim + underscores/points -> spaces + compact spaces
.norm <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[_.]+", " ", x)
  x <- gsub("\\s+", " ", x)
  x
}

# Find first matching column name among candidates (case/underscore/space insensitive)
find_col <- function(df, candidates) {
  nms <- names(df)
  nms_n <- .norm(nms)
  cand_n <- .norm(candidates)
  idx <- match(cand_n, nms_n)
  if (all(is.na(idx))) return(NA_character_)
  nms[na.omit(idx)[1]]
}

standardize_vm_sm_dm <- function(vm, sm, dm) {
  # W4M aliases in first column
  if (!is.null(names(vm)[1]) && identical(names(vm)[1], "variableMetadata")) names(vm)[1] <- "name"
  if (!is.null(names(dm)[1]) && identical(names(dm)[1], "dataMatrix")) names(dm)[1] <- "name"
  if (!is.null(names(sm)[1]) && identical(names(sm)[1], "sampleMetadata")) names(sm)[1] <- "sample_name"

  # m/z
  mz_col <- find_col(vm, c("mz", "m/z", "row m/z", "row mz"))
  if (is.na(mz_col)) stop("variableMetadata must contain a 'mz' column (aliases: m/z, row m/z).")
  names(vm)[names(vm) == mz_col] <- "mz"

  # RT
  rt_col <- find_col(vm, c("rt", "retention time", "row retention time (sec)", "row retention time", "retention time (min)", "rt [min]", "row retention time (min)"))
  if (is.na(rt_col)) stop("variableMetadata must contain a 'rt' column (aliases: retention time, rt [min], etc.).")
  names(vm)[names(vm) == rt_col] <- "rt"

  # name
  nm_col <- find_col(vm, c("name"))
  if (is.na(nm_col)) {
    if (!is.null(names(vm)[1])) names(vm)[1] <- "name" else stop("variableMetadata must contain a 'name' column as first column.")
  } else {
    names(vm)[names(vm) == nm_col] <- "name"
  }

  # sample_name
  sm_name_col <- find_col(sm, c("sample_name"))
  if (is.na(sm_name_col)) stop("sampleMetadata must contain 'sample_name' (alias 'sampleMetadata' automatically renamed if first column).")
  names(sm)[names(sm) == sm_name_col] <- "sample_name"

  # name in DM
  dm_name_col <- find_col(dm, c("name"))
  if (is.na(dm_name_col)) {
    if (!is.null(names(dm)[1])) names(dm)[1] <- "name" else stop("dataMatrix must contain a 'name' column as first column.")
  } else {
    names(dm)[names(dm) == dm_name_col] <- "name"
  }

  vm <- vm %>% mutate(name = as.character(name),
                      mz   = suppressWarnings(as.numeric(mz)),
                      rt   = suppressWarnings(as.numeric(rt)))
  list(vm = vm, sm = sm, dm = dm)
}

rt_to_minutes <- function(rt) {
  if (all(is.na(rt))) return(rt)
  med <- stats::median(rt, na.rm = TRUE)
  if (is.finite(med) && med > 100) rt/60 else rt
}

extract_numeric_id <- function(x) {
  x <- as.character(x %||% NA_character_)
  out <- stringr::str_extract(x, "(?<!\\d)\\d+(?!\\d)")
  suppressWarnings(as.integer(out))
}

`%||%` <- function(a,b) if (!is.null(a)) a else b

detect_id_in_block <- function(block_lines) {
  if (length(block_lines) == 0) return(list(id = NA_integer_, lines = block_lines))
  key_lines <- block_lines[grepl("^(TITLE|SCANS|FEATURE[_ ]?ID|FEATUREFULLID|FEATURE_FULL_ID|SPECTRUMID|FEATUREID)\\s*[:=]", block_lines, ignore.case = TRUE)]
  cat_line <- paste(block_lines, collapse = " ")
  cand <- NA_character_
  if (length(key_lines) > 0) {
    tmp <- stringr::str_extract(paste(key_lines, collapse = " "), "(?<!\\d)\\d+(?!\\d)")
    cand <- tmp
  }
  if (is.na(cand)) {
    tmp2 <- stringr::str_extract(cat_line, "(?<!\\d)\\d+(?!\\d)")
    cand <- tmp2
  }
  id <- suppressWarnings(as.integer(cand))
  list(id = id, lines = block_lines)
}

stream_filter_mgf <- function(in_path, out_path, keep_ids) {
  con_in  <- file(in_path, open = "r", blocking = TRUE)
  on.exit(close(con_in), add = TRUE)
  con_out <- file(out_path, open = "w")
  on.exit(close(con_out), add = TRUE)

  in_block <- FALSE
  block <- character()
  total_blocks <- 0L
  kept_blocks  <- 0L
  seen_ids <- integer()

  while (TRUE) {
    ln <- readLines(con_in, n = 1, warn = FALSE)
    if (length(ln) == 0) {
      if (in_block && length(block) > 0) {
        total_blocks <- total_blocks + 1L
        det <- detect_id_in_block(block)
        if (!is.na(det$id) && det$id %in% keep_ids) {
          writeLines(det$lines, con_out)
          writeLines("", con_out)
          kept_blocks <- kept_blocks + 1L
          seen_ids <- unique(c(seen_ids, det$id))
        }
      }
      break
    }
    if (identical(ln, "BEGIN IONS")) {
      in_block <- TRUE
      block <- "BEGIN IONS"
    } else if (identical(ln, "END IONS")) {
      block <- c(block, "END IONS")
      in_block <- FALSE
      total_blocks <- total_blocks + 1L
      det <- detect_id_in_block(block)
      if (!is.na(det$id) && det$id %in% keep_ids) {
        writeLines(det$lines, con_out)
        writeLines("", con_out)
        kept_blocks <- kept_blocks + 1L
        seen_ids <- unique(c(seen_ids, det$id))
      }
      block <- character()
    } else {
      if (in_block) block <- c(block, ln)
    }
  }
  list(total_blocks = total_blocks, kept_blocks = kept_blocks, seen_ids = seen_ids)
}

make_zip_safe <- function(zip_path, files) {
  ok <- FALSE
  if (requireNamespace("zip", quietly = TRUE)) {
    try({
      zip::zipr(zipfile = zip_path, files = files, include_directories = FALSE)
      ok <- TRUE
    }, silent = TRUE)
  }
  if (!ok) {
    z <- try(utils::zip(zipfile = zip_path, files = files, flags = "-r9Xq"), silent = TRUE)
    if (!inherits(z, "try-error") && file.exists(zip_path)) ok <- TRUE
  }
  if (!ok) return(NULL)
  zip_path
}

ui <- fluidPage(
  tags$head(tags$title("JaFaR - MZmine to W4M - bidirectional")),
  h2("JaFaR - MZmine to W4M - bidirectional"),
  tabsetPanel(
    tabPanel("MZmine → W4M",
      br(),
      fluidRow(
        column(6,
          fileInput("mzmine_quant", "MZmine quant (CSV/TSV/TXT/XLSX/TABULAR)", accept = c(".csv",".tsv",".txt",".tab",".tabular",".xlsx")),
          fileInput("sm_in", "sampleMetadata (TSV/TXT/CSV/XLSX/TABULAR)", accept = c(".tsv",".txt",".csv",".xlsx",".tab",".tabular"))
        ),
        column(6,
          checkboxInput("include_meta_cols", "Include non-intensity columns from quant in variableMetadata (prefix meta_)", value = FALSE),
          actionButton("build_w4m", "Build W4M tables")
        )
      ),
      hr(),
      h4("Previews"),
      tableOutput("preview_vm1"),
      tableOutput("preview_dm1"),
      tableOutput("preview_sm1"),
      hr(),
      h4("Downloads"),
      downloadButton("dl_vm1", "variableMetadata.tsv"),
      downloadButton("dl_dm1", "dataMatrix.tsv"),
      downloadButton("dl_sm1", "sampleMetadata.tsv"),
      downloadButton("dl_zip1", "ZIP (3 files)")
    ),
    tabPanel("W4M → MZmine",
      br(),
      fluidRow(
        column(6,
          fileInput("vm_in", "variableMetadata (TSV/TXT/CSV/XLSX/TABULAR)", accept = c(".tsv",".txt",".csv",".xlsx",".tab",".tabular")),
          fileInput("dm_in", "dataMatrix (TSV/TXT/CSV/XLSX/TABULAR)", accept = c(".tsv",".txt",".csv",".xlsx",".tab",".tabular")),
          fileInput("sm_in2","sampleMetadata (TSV/TXT/CSV/XLSX/TABULAR)", accept = c(".tsv",".txt",".csv",".xlsx",".tab",".tabular"))
        ),
        column(6,
          fileInput("mgf_in", "Original MGF (very large files supported)", accept = c(".mgf",".txt")),
          checkboxInput("extract_numeric", "Extract numeric ID from vm$name (e.g., X14427 → 14427)", value = TRUE),
          actionButton("run_filter", "Filter MGF + rebuild quant")
        )
      ),
      hr(),
      h4("Diagnostics"),
      verbatimTextOutput("diag_counts"),
      h5("Expected IDs (from vm$name)"),
      DT::dataTableOutput("ids_expected_tbl"),
      h5("IDs detected in MGF"),
      DT::dataTableOutput("ids_seen_tbl"),
      h5("Missing IDs (expected but not found)"),
      DT::dataTableOutput("ids_missing_tbl"),
      h5("Extra IDs (found in MGF but absent from W4M)"),
      DT::dataTableOutput("ids_extra_tbl"),
      hr(),
      h4("Downloads"),
      downloadButton("dl_mgf2", "filtered_STRICT_numeric.mgf"),
      downloadButton("dl_quant2", "quant_from_W4M.csv"),
      downloadButton("dl_sm2", "sampleMetadata.tsv"),
      downloadButton("dl_zip2", "ZIP (MGF + quant + sampleMetadata)")
    )
  )
)

server <- function(input, output, session) {

  # ===== Tab 1 : MZmine -> W4M =====
  mzmine_quant <- reactive({
    req(input$mzmine_quant)
    # Keep MZmine headers as-is
    read_any_table(input$mzmine_quant$datapath)
  })

  sm1 <- reactive({
    req(input$sm_in)
    read_any_table(input$sm_in$datapath)
  })

  vm1_dm1_sm1 <- eventReactive(input$build_w4m, {
    q <- mzmine_quant()
    sm <- sm1()

    id_col <- find_col(q, c("row id","# featureid","featureid","id"))
    mz_col <- find_col(q, c("row m/z","m/z","mz","row mz","row_m/z","row_mz"))
    rt_col <- find_col(q, c("row retention time (sec)","row retention time","retention time","rt","rt [min]","row retention time (min)","retention_time"))

    if (is.na(id_col) || is.na(mz_col) || is.na(rt_col)) {
      stop("Missing essential columns in MZmine quant (ID/mz/rt).")
    }

    q <- q %>%
      rename(feature_id = !!id_col,
             mz = !!mz_col,
             rt_raw = !!rt_col) %>%
      mutate(name = as.character(feature_id))

    # Heuristic: if median RT < 30, RT is likely in minutes → convert to seconds for W4M VM? (Keep as-is; DM uses intensities)
    rt_med <- stats::median(q$rt_raw, na.rm = TRUE)
    if (is.finite(rt_med) && rt_med < 30) {
      q <- q %>% mutate(rt = rt_raw * 60)
    } else {
      q <- q %>% mutate(rt = rt_raw)
    }

    sm_name_col <- find_col(sm, c("sample_name"))
    if (is.na(sm_name_col)) {
      if (identical(names(sm)[1], "sampleMetadata")) names(sm)[1] <- "sample_name"
      sm_name_col <- find_col(sm, c("sample_name"))
    }
    if (is.na(sm_name_col)) stop("sampleMetadata must contain 'sample_name'.")
    names(sm)[names(sm) == sm_name_col] <- "sample_name"

    sample_names <- sm$sample_name %>% as.character()
    intensity_cols <- character()
    for (s in sample_names) {
      cand <- c(s, paste0(s, ".mzML Peak area"))
      hit <- find_col(q, cand)
      if (is.na(hit)) stop(sprintf("Intensity column not found for sample '%s' (tried: %s).", s, paste(cand, collapse=" / ")))
      intensity_cols <- c(intensity_cols, hit)
    }

    dm <- q %>%
      select(name, all_of(intensity_cols)) %>%
      { setNames(., c("name", sample_names)) }

    keep_vm <- c("name","mz","rt")
    vm <- q %>% select(any_of(keep_vm))

    if (isTRUE(input$include_meta_cols)) {
      non_intensity <- setdiff(names(q), c("feature_id","mz","rt_raw","rt","name", intensity_cols))
      if (length(non_intensity) > 0) {
        vm <- q %>% select(any_of(keep_vm), any_of(non_intensity)) %>%
          rename_with(~ifelse(.x %in% non_intensity, paste0("meta_", .x), .x), all_of(non_intensity))
      }
    }

    list(vm = vm, dm = dm, sm = sm)
  }, ignoreInit = TRUE)

  output$preview_vm1 <- renderTable({ head(vm1_dm1_sm1()$vm, 10) })
  output$preview_dm1 <- renderTable({ head(vm1_dm1_sm1()$dm, 10) })
  output$preview_sm1 <- renderTable({ head(vm1_dm1_sm1()$sm, 10) })

  output$dl_vm1 <- downloadHandler(
    filename = function() "variableMetadata.tsv",
    content  = function(file) readr::write_tsv(vm1_dm1_sm1()$vm, file)
  )
  output$dl_dm1 <- downloadHandler(
    filename = function() "dataMatrix.tsv",
    content  = function(file) readr::write_tsv(vm1_dm1_sm1()$dm, file)
  )
  output$dl_sm1 <- downloadHandler(
    filename = function() "sampleMetadata.tsv",
    content  = function(file) readr::write_tsv(vm1_dm1_sm1()$sm, file)
  )
  output$dl_zip1 <- downloadHandler(
    filename = function() "MZmine_to_W4M.zip",
    content = function(file) {
      tmpdir <- tempfile("w4m1_"); dir.create(tmpdir)
      f_vm <- file.path(tmpdir, "variableMetadata.tsv"); readr::write_tsv(vm1_dm1_sm1()$vm, f_vm)
      f_dm <- file.path(tmpdir, "dataMatrix.tsv");       readr::write_tsv(vm1_dm1_sm1()$dm, f_dm)
      f_sm <- file.path(tmpdir, "sampleMetadata.tsv");   readr::write_tsv(vm1_dm1_sm1()$sm, f_sm)
      z <- make_zip_safe(file, c(f_vm,f_dm,f_sm))
      if (is.null(z)) stop("Could not create a ZIP on this system. Please download files individually.")
    }
  )

  # ===== Tab 2 : W4M -> MZmine =====
  vm2 <- reactive({ req(input$vm_in); read_any_table(input$vm_in$datapath) })
  dm2 <- reactive({ req(input$dm_in); read_any_table(input$dm_in$datapath) })
  sm2 <- reactive({ req(input$sm_in2); read_any_table(input$sm_in2$datapath) })

  vals <- reactiveValues(
    vm = NULL, dm = NULL, sm = NULL,
    ids_expected = NULL, ids_seen = NULL, ids_missing = NULL, ids_extra = NULL,
    diag = NULL, mgf_out_path = NULL, quant_out_path = NULL, sm_out_path = NULL
  )

  observeEvent(input$run_filter, {
    req(input$mgf_in)
    vm <- vm2(); dm <- dm2(); sm <- sm2()

    std <- standardize_vm_sm_dm(vm, sm, dm)
    vm <- std$vm; dm <- std$dm; sm <- std$sm

    if (!all(c("name","mz","rt") %in% names(vm))) stop("variableMetadata must contain: name, mz, rt")
    if (!"sample_name" %in% names(sm)) stop("sampleMetadata must contain: sample_name")
    if (!"name" %in% names(dm)) stop("dataMatrix must contain: name")

    # Expected numeric IDs
    ids_expected <- extract_numeric_id(vm$name)
    ids_expected <- ids_expected[!is.na(ids_expected)] %>% unique() %>% sort()

    # MGF Streaming filter
    mgf_in  <- input$mgf_in$datapath
    mgf_out <- file.path(tempdir(), "filtered_STRICT_numeric.mgf")
    diag <- stream_filter_mgf(mgf_in, mgf_out, keep_ids = ids_expected)
    ids_seen <- sort(unique(diag$seen_ids))

    ids_missing <- setdiff(ids_expected, ids_seen)
    ids_extra   <- setdiff(ids_seen, ids_expected)

    # MZmine-like quant with numeric "row ID"
    vm_rid <- extract_numeric_id(vm$name)
    dm_rid <- extract_numeric_id(dm$name)

    vm_exp <- vm %>%
      mutate(`row ID` = vm_rid,
             `row m/z` = mz,
             `row retention time` = rt_to_minutes(rt)) %>%
      select(`row ID`,`row m/z`,`row retention time`)

    sm_names <- sm$sample_name %>% as.character()
    missing_dm_cols <- setdiff(sm_names, names(dm))
    if (length(missing_dm_cols) > 0) {
      stop(sprintf("Missing intensity columns in dataMatrix for: %s", paste(missing_dm_cols, collapse=", ")))
    }

    q2 <- dm %>%
      transmute(`row ID` = dm_rid, across(all_of(sm_names))) %>%
      distinct(`row ID`, .keep_all = TRUE)

    names(q2) <- c("row ID", paste0(sm_names, ".mzML Peak area"))

    quant_out <- vm_exp %>% dplyr::left_join(q2, by = "row ID")

    out_quant <- file.path(tempdir(), "quant_from_W4M.csv")
    readr::write_csv(quant_out, out_quant)

    out_sm <- file.path(tempdir(), "sampleMetadata.tsv")
    readr::write_tsv(sm, out_sm)

    zipfile <- file.path(tempdir(), "W4M_to_MZmine_STREAM.zip")
    z <- make_zip_safe(zipfile, c(mgf_out, out_quant, out_sm))

    vals$vm <- vm; vals$dm <- dm; vals$sm <- sm
    vals$ids_expected <- ids_expected
    vals$ids_seen <- ids_seen
    vals$ids_missing <- ids_missing
    vals$ids_extra <- ids_extra
    vals$diag <- diag
    vals$mgf_out_path <- mgf_out
    vals$quant_out_path <- out_quant
    vals$sm_out_path <- out_sm
    vals$zip2 <- if (!is.null(z)) z else NULL
  })

  output$diag_counts <- renderText({
    req(vals$diag)
    sprintf("MGF blocks read: %d | kept: %d | unique IDs detected: %d",
            vals$diag$total_blocks, vals$diag$kept_blocks, length(vals$ids_seen %||% integer()))
  })

  output$ids_expected_tbl <- DT::renderDataTable({
    req(vals$ids_expected)
    DT::datatable(data.frame(id = vals$ids_expected), options = list(pageLength = 10))
  })
  output$ids_seen_tbl <- DT::renderDataTable({
    req(vals$ids_seen)
    DT::datatable(data.frame(id = vals$ids_seen), options = list(pageLength = 10))
  })
  output$ids_missing_tbl <- DT::renderDataTable({
    req(vals$ids_missing)
    DT::datatable(data.frame(id = vals$ids_missing), options = list(pageLength = 10))
  })
  output$ids_extra_tbl <- DT::renderDataTable({
    req(vals$ids_extra)
    DT::datatable(data.frame(id = vals$ids_extra), options = list(pageLength = 10))
  })

  output$dl_mgf2 <- downloadHandler(
    filename = function() "filtered_STRICT_numeric.mgf",
    content  = function(file) file.copy(vals$mgf_out_path, file, overwrite = TRUE)
  )
  output$dl_quant2 <- downloadHandler(
    filename = function() "quant_from_W4M.csv",
    content  = function(file) file.copy(vals$quant_out_path, file, overwrite = TRUE)
  )
  output$dl_sm2 <- downloadHandler(
    filename = function() "sampleMetadata.tsv",
    content  = function(file) file.copy(vals$sm_out_path, file, overwrite = TRUE)
  )
  output$dl_zip2 <- downloadHandler(
    filename = function() "W4M_to_MZmine_STREAM.zip",
    content  = function(file) {
      if (is.null(vals$zip2)) stop("Could not create a ZIP on this system. Please download files individually.")
      file.copy(vals$zip2, file, overwrite = TRUE)
    }
  )
}

if (interactive()) shinyApp(ui, server)
