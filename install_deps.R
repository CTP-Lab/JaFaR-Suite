# install_deps.R
# Auto-detect + install dependencies for JaFaR-Suite Shiny apps

detect_packages <- function(apps_dir = "apps") {
  if (!dir.exists(apps_dir)) {
    warning("apps_dir not found: ", apps_dir)
    return(character(0))
  }

  app_files <- list.files(apps_dir, pattern = "app\\.R$", recursive = TRUE, full.names = TRUE)
  if (length(app_files) == 0) {
    warning("No app.R files found under: ", apps_dir)
    return(character(0))
  }

  pkgs <- character(0)

  # Regex patterns
  re_lib  <- "library\\s*\\(\\s*([A-Za-z0-9._]+)\\s*\\)"
  re_req  <- "require\\s*\\(\\s*([A-Za-z0-9._]+)\\s*(,|\\))"
  re_reqn <- "requireNamespace\\s*\\(\\s*[\"']([A-Za-z0-9._]+)[\"']\\s*,\\s*quietly\\s*=\\s*TRUE\\s*\\)"
  re_ns   <- "([A-Za-z0-9._]+)\\s*::\\s*[A-Za-z0-9._]+"

  for (f in app_files) {
    x <- tryCatch(readLines(f, warn = FALSE, encoding = "UTF-8"), error = function(e) character(0))
    if (length(x) == 0) next

    # Strip simple comments to reduce false hits
    x2 <- sub("#.*$", "", x)

    m1 <- regmatches(x2, gregexpr(re_lib,  x2, perl = TRUE))
    m2 <- regmatches(x2, gregexpr(re_req,  x2, perl = TRUE))
    m3 <- regmatches(x2, gregexpr(re_reqn, x2, perl = TRUE))
    m4 <- regmatches(x2, gregexpr(re_ns,   x2, perl = TRUE))

    # Extract capture groups
    extract_group <- function(matches, pattern, group = 1) {
      out <- character(0)
      for (vec in matches) {
        if (length(vec) == 0) next
        for (s in vec) {
          g <- regexec(pattern, s, perl = TRUE)
          r <- regmatches(s, g)[[1]]
          if (length(r) >= group + 1) out <- c(out, r[group + 1])
        }
      }
      out
    }

    pkgs <- c(
      pkgs,
      extract_group(m1, re_lib,  1),
      extract_group(m2, re_req,  1),
      extract_group(m3, re_reqn, 1),
      extract_group(m4, re_ns,   1)
    )
  }

  pkgs <- unique(pkgs[nzchar(pkgs)])

  # Remove base/recommended packages that do not need installation
  base_pkgs <- c(
    "base","compiler","datasets","graphics","grDevices","grid","methods","parallel",
    "splines","stats","stats4","tcltk","tools","utils"
  )

  setdiff(pkgs, base_pkgs)
}

install_deps <- function(pkgs = NULL,
                         apps_dir = "apps",
                         repos = "https://cloud.r-project.org",
                         try_bioc = TRUE,
                         quiet = FALSE) {

  if (is.null(pkgs)) {
    pkgs <- detect_packages(apps_dir = apps_dir)
  }

  # Always ensure shiny is available (core requirement)
  pkgs <- unique(c("shiny", pkgs))
  pkgs <- pkgs[nzchar(pkgs)]

  is_installed <- function(p) requireNamespace(p, quietly = TRUE)

  missing <- pkgs[!vapply(pkgs, is_installed, FUN.VALUE = logical(1))]

  if (length(missing) == 0) {
    if (!quiet) message("All dependencies are already installed.")
    return(invisible(character(0)))
  }

  if (!quiet) {
    message("Installing missing packages (CRAN first):\n- ", paste(missing, collapse = "\n- "))
  }

  # Try CRAN install
  tryCatch(
    install.packages(missing, repos = repos),
    error = function(e) {
      warning("CRAN install encountered an error: ", conditionMessage(e))
    }
  )

  # Re-check what is still missing
  still_missing <- missing[!vapply(missing, is_installed, FUN.VALUE = logical(1))]

  # Optionally try Bioconductor for remaining packages
  if (try_bioc && length(still_missing) > 0) {
    if (!is_installed("BiocManager")) {
      if (!quiet) message("Installing BiocManager (required to try Bioconductor packages)...")
      install.packages("BiocManager", repos = repos)
    }
    if (is_installed("BiocManager")) {
      if (!quiet) {
        message("Trying Bioconductor install for remaining packages:\n- ",
                paste(still_missing, collapse = "\n- "))
      }
      tryCatch(
        BiocManager::install(still_missing, ask = FALSE, update = FALSE),
        error = function(e) {
          warning("Bioconductor install encountered an error: ", conditionMessage(e))
        }
      )
    }
  }

  # Final check
  final_missing <- pkgs[!vapply(pkgs, is_installed, FUN.VALUE = logical(1))]

  if (length(final_missing) > 0) {
    warning(
      "Some packages could not be installed automatically:\n- ",
      paste(final_missing, collapse = "\n- "),
      "\nCheck whether they are on CRAN/Bioconductor, or require system libraries."
    )
  } else {
    if (!quiet) message("Done. All dependencies installed.")
  }

  invisible(final_missing)
}

list_deps <- function(apps_dir = "apps") {
  detect_packages(apps_dir = apps_dir)
}
