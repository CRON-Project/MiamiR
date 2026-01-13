
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
  packageStartupMessage("Visit the MiamiR Github page (https://github.com/CRON-Project/MiamiR) for requests, feedback, and info.")
  packageStartupMessage("Access vignette (https://cron-project.github.io/MiamiR/) for full tutorial.")
  packageStartupMessage("Companion web app available at https://cron-project.shinyapps.io/MiamiR_Web_App_Pub/ .")
  packageStartupMessage("Zenodo repository (https://zenodo.org/records/18141343) contains archived and latest releases.")
  packageStartupMessage("\n")

}

