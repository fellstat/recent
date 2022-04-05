#' Run Shiny UI
#' @details
#' Runs a Shiny application
#' @export
shiny_rita <- function() {
  app_dir <- system.file("shiny_ui", package = "recent")
  if (app_dir == "") {
    stop("Could not find example directory. Try re-installing `recent`.", call. = FALSE)
  }
  shiny::runApp(app_dir, display.mode = "normal")
}
