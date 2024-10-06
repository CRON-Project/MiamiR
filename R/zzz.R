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
  packageStartupMessage("Citation not required, but greatly appreciated :)")
  packageStartupMessage("Visit Github (https://github.com/CRON-Project/MiamiR) for requests, feedback, and info.")
  packageStartupMessage("Visit www.theCRONproject.com for more data visualisation examples.")
  packageStartupMessage("\n")

    }


#Need non-CRAN ggmanh first

.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  if (!requireNamespace("ggmanh", quietly = TRUE)) {
    message("ggmanh package not found. Installing it from Bioconductor...")
    BiocManager::install("ggmanh")
  }
}

.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("ggmanh", quietly = TRUE)) {
    message("ggmanh package not found. Installing it from Bioconductor...")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("ggmanh")
  }
}
