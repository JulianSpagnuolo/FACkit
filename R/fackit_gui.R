#' FACkit GUI
#'
#' @export
#'
#'
fackit <- function() {

  appDir <- system.file("fackit_gui", "app.R", package = "FACkit")
  if (appDir == "") {
    stop("Could not find fackit_gui directory. Try re-installing `FACkit`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}
