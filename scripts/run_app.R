run_app <- function(app_dir) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required. Please install it first.", call. = FALSE)
  }

  if (!dir.exists(app_dir)) {
    stop(paste0("App directory not found: ", app_dir), call. = FALSE)
  }

  shiny::runApp(appDir = app_dir, launch.browser = TRUE)
}
