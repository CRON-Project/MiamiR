.onAttach <- function(libname, pkgname) {
  manhattan_plot_diagram <- "
    ┌───────────────────────────────────────────────────┐
    │                                                   │
    │                      MiamiR                       │
    │                                                   │
    │       █          █                                │
    │       █          █        █                       │
    │       █   █      █        █                       │
    │   █   █   █      █    █   █                  █    │
    │   █   █   █      █    █   █    █             █    │
    │   █   █   █      █    █   █    █   █         █    │
    │ █ █ █ █ █ █      █ █  █ █ █ █  █   █   █     █    │
    │ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ │
    │ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ │
    └───────────────────────────────────────────────────┘
  "

  packageStartupMessage(manhattan_plot_diagram)
  packageStartupMessage("Thank you for using MiamiR! Citation not required, but greatly appreciated :)")
  packageStartupMessage("Visit Github for requests, feedback, and info, and the www.theCRONproject.com for more data visualisation examples.")
}
