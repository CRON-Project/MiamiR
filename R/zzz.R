
.onAttach <- function(libname, pkgname) {

  manhattan_plot_diagram <- "

    +-------------------------------------+
    |                MiamiR               |
    |                                     |
    |                                     |
    |   *           *                     |
    |   *           *         *           |
    |   *    *      *         *           |
    | * *    *      *    *    *     *     |
    | * *    *      *    *    *     *     |
    | * *    *      *    *    *     *     |
    | * * * * * * * * * * * * * * * * * * |
    | * * * * * * * * * * * * * * * * * * |
    | * * * * * * * * * * * * * * * * * * |
    +-------------------------------------+

  "
  packageStartupMessage(manhattan_plot_diagram)
  packageStartupMessage("Thank you for using MiamiR!")
  packageStartupMessage("Citation greatly appreciated :)")
  packageStartupMessage("Visit the MiamiR Github page (https://github.com/CRON-Project/MiamiR)")
  packageStartupMessage("For requests, feedback, and info.")
  packageStartupMessage("Access vignette (https://cron-project.github.io/MiamiR/) for full tutorial.")
  packageStartupMessage("Companion web app available at https://cron-project.shinyapps.io/MiamiR_Web_App_Pub/ .")
  packageStartupMessage("Zenodo repository (https://zenodo.org/records/18141343) contains archived and latest releases.")
  packageStartupMessage("\n")

  # Run would you like METASOFT too?
  if (!interactive()) return(invisible(NULL))

  final_jar <- file.path(tools::R_user_dir("MiamiR", "data"), "java", "Metasoft.jar")

  if (!file.exists(final_jar)) {

    ans <- readline("METASOFT not installed. Install now? [y/N] ")

    if (tolower(substr(ans, 1, 1)) %in% "y") {
      install_metasoft_simple()
    } else {
      packageStartupMessage("Skipping METASOFT install. You can run install_metasoft_simple() later.")
    }

  } else {

    # Set for immediate use (optional but useful)
    assign("METASOFT_jar_location", final_jar, envir = .GlobalEnv)

    packageStartupMessage(
      paste0(
        "METASOFT already installed.\n",
        "Path: ", final_jar, "\n",
        "To reinstall/overwrite, run: install_metasoft_simple()\n"
      )
    )

  }

}
