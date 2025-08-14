  library(shiny)

  library(ggplot2)
  library(vroom)
  library(dplyr)
  library(plotly)
  library(htmlwidgets)
  library(shinyBS)
  library(mailR)
  library(googledrive)



  # Fallback operator: use 'a' if not NULL, otherwise 'b'
  `%||%` <- function(a, b) {
    if (!is.null(a)) a else b
  }



  # Authenticate once (automatically)
  drive_auth(path = "C:/Users/callumon/Miami_Package_R/MiamiR/miamir-storage-361ee7c11e56.json")


  if (!googledrive::drive_has_token()) {
    googledrive::drive_auth(path = "C:/Users/callumon/Miami_Package_R/MiamiR/miamir-storage-361ee7c11e56.json")
  }


  # Increase file upload limit to 10 GB
  options(shiny.maxRequestSize = 10 * 1024^3)
  options(shiny.http.timeout = 300)  # Timeout increased to 5 minutes
  options(shiny.maxRequestSize = 10 * 1024^3)  # 30 MB


  # Helper function to generate dynamic UI inputs, specifically handling list inputs
  generate_inputs_with_defaults <- function(func, saved_values = list()) {
    arg_list <- formals(func)

    inputs <- lapply(names(arg_list), function(arg_name) {
      if (arg_name %in% c("Data", "Top_Data", "Bottom_Data")) return(NULL)

      arg_value <- arg_list[[arg_name]]
      if (is.symbol(arg_value) || missing(arg_value)) return(NULL)
      if (is.call(arg_value)) arg_value <- eval(arg_value)

      # Use saved value if present
      value_to_use <- saved_values[[arg_name]] %||% arg_value

      if (is.numeric(value_to_use) && length(value_to_use) > 1) {
        textInput(arg_name, label = arg_name, value = paste(value_to_use, collapse = ","))
      } else if (is.character(value_to_use) && length(value_to_use) > 1) {
        textInput(arg_name, label = arg_name, value = paste(value_to_use, collapse = ","))
      } else if (is.numeric(value_to_use)) {
        numericInput(arg_name, label = arg_name, value = value_to_use)
      } else if (is.character(value_to_use)) {
        textInput(arg_name, label = arg_name, value = value_to_use)
      } else if (is.logical(value_to_use)) {
        checkboxInput(arg_name, label = arg_name, value = value_to_use)
      } else {
        textInput(arg_name, label = arg_name, value = as.character(value_to_use))
      }
    })

    return(Filter(Negate(is.null), inputs))
  }





  generate_combined_inputs_with_defaults <- function(primary_func, inherited_func, saved_values = list()) {
    # Preserve all args, including NULL, by using pairlist merging
    args_combined <- as.list(c(formals(inherited_func), formals(primary_func)))
    args_combined <- args_combined[!names(args_combined) %in% c("Data", "Top_Data", "Bottom_Data")]

    inputs <- lapply(names(args_combined), function(arg_name) {
      val <- args_combined[[arg_name]]

      # Safely evaluate each default value
      arg_value <- tryCatch({
        if (is.symbol(val)) {
          if (deparse(val) == "NULL") NULL else return(NULL)  # Accept NULL, skip ...
        } else if (is.call(val)) {
          eval(val)
        } else {
          val
        }
      }, error = function(e) NULL)

      # Use saved value if available
      value_to_use <- saved_values[[arg_name]] %||% arg_value


      # UI rendering
      if (is.null(value_to_use)) {
        return(textInput(arg_name, label = arg_name, value = "", placeholder = "NULL (default)"))
      } else if (is.numeric(value_to_use) && length(value_to_use) > 1) {
        textInput(arg_name, label = arg_name, value = paste(value_to_use, collapse = ","))
      } else if (is.character(value_to_use) && length(value_to_use) > 1) {
        textInput(arg_name, label = arg_name, value = paste(value_to_use, collapse = ","))
      } else if (is.numeric(value_to_use)) {
        numericInput(arg_name, label = arg_name, value = value_to_use)
      } else if (is.character(value_to_use)) {
        textInput(arg_name, label = arg_name, value = value_to_use)
      } else if (is.logical(value_to_use)) {
        checkboxInput(arg_name, label = arg_name, value = value_to_use)
      } else {
        textInput(arg_name, label = arg_name, value = as.character(value_to_use))
      }
    })

    Filter(Negate(is.null), inputs)
  }





  generate_combined_inputs <- function(primary_func, inherited_func) {
    # Combine arguments from both
    args_primary <- formals(primary_func)
    args_inherited <- formals(inherited_func)

    # Remove duplicates from inherited that exist in primary
    args_combined <- modifyList(as.list(args_inherited), as.list(args_primary))

    # Remove "Data" or any disallowed ones
    args_combined <- args_combined[!names(args_combined) %in% c("Data", "Top_Data", "Bottom_Data")]

    # Render inputs
    inputs <- lapply(names(args_combined), function(arg_name) {
      arg_value <- args_combined[[arg_name]]
      if (missing(arg_value) || is.symbol(arg_value)) return(NULL)
      if (is.call(arg_value)) arg_value <- eval(arg_value)

      if (is.numeric(arg_value) && length(arg_value) > 1) {
        textInput(arg_name, label = arg_name, value = paste(arg_value, collapse = ","))
      } else if (is.character(arg_value) && length(arg_value) > 1) {
        textInput(arg_name, label = arg_name, value = paste(arg_value, collapse = ","))
      } else if (is.numeric(arg_value)) {
        numericInput(arg_name, label = arg_name, value = arg_value)
      } else if (is.character(arg_value)) {
        textInput(arg_name, label = arg_name, value = arg_value)
      } else if (is.logical(arg_value)) {
        checkboxInput(arg_name, label = arg_name, value = arg_value)
      } else {
        textInput(arg_name, label = arg_name, value = as.character(arg_value))
      }
    })

    # Download options
    file_type_input <- textInput("File_Type", label = "File Type", value = "jpg")
    width_input <- numericInput("Width", label = "Plot Width (in)", value = 30)
    height_input <- numericInput("Height", label = "Plot Height (in)", value = 15)
    dpi_input <- numericInput("Quality", label = "Plot DPI", value = 300)

    return(c(
      list(file_type_input, width_input, height_input, dpi_input),
      Filter(Negate(is.null), inputs)
    ))
  }


  # Define UI for the Shiny app
  # Define UI for the Shiny app
  # Define UI for the Shiny app
  # Define UI for the Shiny app
  # Define UI for the Shiny app

  ui <- fluidPage(
    uiOutput("mainUI"),     # This will hold your full tabbed app — only shown after "Continue"
    uiOutput("landingUI")   # This will show first, as the overlay screen
  )






  server <- function(input, output, session) {

    batch_settings <- reactiveValues()


    # Place at the top of server()
    batch_pdfs <- reactiveVal(list())

    batch_regional_indices <- reactiveValues()
    batch_regional_plots <- reactiveValues()


    show_landing <- reactiveVal(TRUE)

    show_menu <- reactiveVal(FALSE)
    selectedTab <- reactiveVal("Home")

    observeEvent(input$goMenu, {
      show_menu(TRUE)
    })



    output$functionSelectionUI <- renderUI({
      tagList(
        checkboxGroupInput("selectedFunctions", "Select Functions to Run:",
                           choices = c("Single_Plot", "Regional_Plot", "Miami_Plot", "Forest_Plot", "Annotate_Data", "Model_Munge"),
                           selected = input$selectedFunctions),
        lapply(input$selectedFunctions, function(fn) {
          actionButton(paste0("edit_", fn), paste("⚙️ Settings for", fn), class = "btn-info", style = "margin-bottom: 10px;")
        })
      )
    })



    observe({
      lapply(c("Single_Plot", "Regional_Plot", "Miami_Plot", "Forest_Plot", "Annotate_Data", "Model_Munge"), function(fn) {
        observeEvent(input[[paste0("edit_", fn)]], {
          showModal(modalDialog(
            title = paste("Configure Settings for", fn),
            size = "l",
            easyClose = FALSE,
            footer = tagList(
              actionButton(paste0("save_", fn), "Save Settings", class = "btn-primary"),
              modalButton("Cancel")
            ),
            uiOutput(paste0("batch_settings_ui_", fn))
          ))
        })
      })
    })



    lapply(c("Single_Plot", "Regional_Plot", "Miami_Plot", "Forest_Plot", "Annotate_Data", "Model_Munge"), function(fn) {
      output[[paste0("batch_settings_ui_", fn)]] <- renderUI({
        saved <- batch_settings[[fn]] %||% list()

        content <- if (fn == "Regional_Plot") {
          generate_combined_inputs_with_defaults(.Regional_Plot_original, .Single_Plot_original, saved_values = saved)
        } else if (fn == "Single_Plot") {
          generate_inputs_with_defaults(.Single_Plot_original, saved_values = saved)
        } else {
          generate_inputs_with_defaults(get(fn), saved_values = saved)
        }


        tagList(
          actionButton(paste0("restore_", fn), "Restore Defaults", class = "btn-danger", style = "margin-bottom: 15px;"),
          content
        )
      })




    })



    lapply(c("Single_Plot", "Regional_Plot", "Miami_Plot", "Forest_Plot", "Annotate_Data", "Model_Munge"), function(fn) {
      observeEvent(input[[paste0("save_", fn)]], {
        args <- if (fn == "Regional_Plot") {
          modifyList(as.list(formals(.Single_Plot_original)), as.list(formals(.Regional_Plot_original)))
        }else if (fn == "Single_Plot") {
          formals(.Single_Plot_original)
        }
        else {
          formals(get(fn))
        }

        user_values <- list()

        for (arg in names(args)) {
          val <- input[[arg]]
          if (!is.null(val)) {
            val <- if (is.character(val) && val == "") NULL else val

            if (!is.null(val) && is.character(val) && grepl(",", val)) {
              val_split <- strsplit(val, ",")[[1]]
              val_split <- trimws(val_split)
              user_values[[arg]] <- if (all(val_split == "")) NULL else val_split
            } else {
              user_values[[arg]] <- val
            }
          }
        }


        batch_settings[[fn]] <- user_values
        removeModal()

        # Optional: Save to user profile if signed in
        if (!user_session$is_guest) {
          saveRDS(batch_settings, file = paste0("user_settings_", user_session$username, ".rds"))
        }
      })
    })


    lapply(c("Single_Plot", "Regional_Plot", "Miami_Plot", "Forest_Plot", "Annotate_Data", "Model_Munge"), function(fn) {
      observeEvent(input[[paste0("restore_", fn)]], {
        batch_settings[[fn]] <- list()  # Wipe saved settings

        # Optional: update user file on disk if signed in
        if (!user_session$is_guest) {
          saveRDS(batch_settings, file = paste0("user_settings_", user_session$username, ".rds"))
        }

        showModal(modalDialog(
          title = "✅ Defaults Restored",
          paste0(fn, " settings have been reset to defaults."),
          easyClose = TRUE
        ))
      })
    })




    output$dashboardMenuUI <- renderUI({
      fluidPage(
        tags$style(HTML("
        .tile-grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(350px, 1fr));
    gap: 40px;
    padding: 60px 80px;
    justify-items: center;
  }

  .tile-button {
    background: linear-gradient(135deg, #4e54c8, #8f94fb);
    color: white;
    font-size: 2.5em;
    font-weight: 600;
    padding: 80px 60px;
    border: none;
    border-radius: 32px;
    text-align: center;
    transition: all 0.3s ease;
    box-shadow: 0 8px 16px rgba(0,0,0,0.25);
    width: 100%;
    max-width: 500px;
    word-wrap: break-word;
    overflow-wrap: break-word;
    line-height: 1.3;
  }

  .tile-button:hover {
    background: linear-gradient(135deg, #8f94fb, #4e54c8);
    cursor: pointer;
    transform: scale(1.08);
  }


      ")),
        h2("Select a Module", style = "text-align:center; padding-top: 20px;"),
        div(class = "tile-grid",
            actionButton("goHome", "Information", class = "tile-button"),
            actionButton("goSingle", "Single_Plot()", class = "tile-button"),
            actionButton("goRegional", "Regional_Plot()", class = "tile-button"),
            actionButton("goMiami", "Miami_Plot()", class = "tile-button"),
            actionButton("goForest", "Forest_Plot()", class = "tile-button"),
            actionButton("goAnnotate", "Annotate_Data()", class = "tile-button"),
            actionButton("goModel", "Model_Munge()", class = "tile-button"),
            actionButton("goBatch", "Batch Mode", class = "tile-button")
        )
      )
    })



    # User authentication reactive storage
    credentials <- reactiveVal(data.frame(username = character(), password = character(), stringsAsFactors = FALSE))
    user_session <- reactiveValues(logged_in = FALSE, username = NULL, is_guest = FALSE)

    # Load stored credentials from file if exists
    if (file.exists("user_credentials.rds")) {
      credentials(readRDS("user_credentials.rds"))
    }

    observeEvent(input$goHome, {
      selectedTab("Home")
      show_menu(FALSE)
    })
    observeEvent(input$goSingle, {
      selectedTab("Single_Plot()")
      show_menu(FALSE)
    })
    observeEvent(input$goRegional, {
      selectedTab("Regional_Plot()")
      show_menu(FALSE)
    })
    observeEvent(input$goMiami, {
      selectedTab("Miami_Plot()")
      show_menu(FALSE)
    })
    observeEvent(input$goForest, {
      selectedTab("Forest_Plot()")
      show_menu(FALSE)
    })
    observeEvent(input$goAnnotate, {
      selectedTab("Annotate_Data()")
      show_menu(FALSE)
    })
    observeEvent(input$goModel, {
      selectedTab("Model_Munge()")
      show_menu(FALSE)
    })
    observeEvent(input$goBatch, {
      selectedTab("Automated Batch Mode")
      show_menu(FALSE)
    })


    observeEvent(input$auth_submit, {
      if (input$auth_choice == "Continue as Guest") {
        user_session$logged_in <- TRUE
        show_landing(FALSE)
        show_menu(TRUE)

        user_session$is_guest <- TRUE


      } else if (input$auth_choice == "Login") {
        users <- credentials()
        match <- users %>%
          filter(username == input$auth_username, password == input$auth_password)

        if (nrow(match) == 1) {
          user_session$logged_in <- TRUE
          user_session$is_guest <- FALSE
          user_session$username <- input$auth_username
          show_landing(FALSE)
        } else {
          showModal(modalDialog("❌ Invalid login credentials", easyClose = TRUE))
        }


        if (!user_session$is_guest && file.exists(paste0("user_settings_", input$auth_username, ".rds"))) {
          batch_settings_loaded <- readRDS(paste0("user_settings_", input$auth_username, ".rds"))
          for (name in names(batch_settings_loaded)) {
            batch_settings[[name]] <- batch_settings_loaded[[name]]
          }
        }


      } else if (input$auth_choice == "Register") {
        users <- credentials()
        if (input$auth_username %in% users$username) {
          showModal(modalDialog("❌ Username already exists", easyClose = TRUE))
        } else {
          new_user <- data.frame(username = input$auth_username,
                                 password = input$auth_password,
                                 stringsAsFactors = FALSE)
          all_users <- rbind(users, new_user)
          credentials(all_users)
          saveRDS(all_users, "user_credentials.rds")  # Persist between sessions
          showModal(modalDialog("✅ Registered successfully. You can now login.", easyClose = TRUE))
        }
      }
    })


    output$landingUI <- renderUI({
      req(show_landing())

      fluidPage(
        tags$head(
          tags$style(HTML("


    .login-container {
      display: flex;
      flex-direction: column;
      align-items: center;
      justify-content: flex-start;
      height: 100vh;
      padding-top: 80px;
      background: linear-gradient(135deg, #1e3c72, #2a5298);
      color: white;
      font-family: 'Segoe UI', sans-serif;
      text-align: center;
      padding: 20px;
    }

    .login-container h2 {
      font-size: 3em;
      margin-bottom: 30px;
    }

    .login-container label,
    .login-container input,
    .login-container .btn {
      font-size: 1.5em;
    }

    .login-container .radio {
      margin-bottom: 30px;
    }

    .login-container input[type='text'],
    .login-container input[type='password'] {
      width: 300px;
      height: 45px;
      border-radius: 8px;
      border: none;
      padding-left: 15px;
      margin-bottom: 20px;
    }

   .login-container {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: flex-start;
    height: 100vh;
    padding-top: 220px;  /* ⬅️ Increase to move everything down */
    background: linear-gradient(135deg, #1e3c72, #2a5298);
    color: white;
    font-family: 'Segoe UI', sans-serif;
    text-align: center;
    padding: 20px;
  }

    .login-container .btn:hover {
      background-color: #e67e22;
    }

    .login-container .radio label {
      font-size: 1.2em;
      margin-right: 20px;
    }

  .manhattan-wrapper {
    display: flex;
    justify-content: center;
    margin-top: 80px;
    width: 100%;
  }

  .login-inner {
    margin-top: 100px; /* ⬅️ This moves the whole block down */
    display: flex;
    flex-direction: column;
    align-items: center;
  }


  .manhattan-preview {
    position: relative;
    height: 200px;
    width: 160px;
    border-top: 2px dashed white;
  }
  .chrom {
    position: absolute;
    width: 10px;
    border-radius: 2px;
  }

  .up {
    bottom: 50%;                /* Anchored at midpoint */
    transform-origin: bottom;
  }

  .down {
    top: 50%;                   /* Anchored at midpoint */
    transform-origin: top;
  }

  .h50 { height: 50px; }
  .h70 { height: 70px; }
  .h90 { height: 90px; }
  .h110 { height: 110px; }
  .h130 { height: 130px; }

  /* Colors: up green, down blue */
  .up {
    background-color: #2ecc71;
  }

  .down {
    background-color: #3498db;
  }

  "))

        ),
        div(class = "login-container",
            div(class = "login-inner",
                tags$h2("Welcome to the MiamiR App"),
                radioButtons("auth_choice", label = NULL,
                             choices = c("Login", "Register", "Continue as Guest"),
                             selected = "Login",
                             inline = TRUE),
                conditionalPanel(
                  condition = "input.auth_choice != 'Continue as Guest'",
                  textInput("auth_username", "Username"),
                  passwordInput("auth_password", "Password")
                ),
                actionButton("auth_submit", "Continue", class = "btn"),

                div(class = "manhattan-wrapper",
                    div(class = "manhattan-preview",
                        div(class = "chrom up h90", style = "left: 0px;"),
                        div(class = "chrom down h110", style = "left: 15px;"),
                        div(class = "chrom up h70", style = "left: 30px;"),
                        div(class = "chrom down h90", style = "left: 45px;"),
                        div(class = "chrom up h50", style = "left: 60px;"),
                        div(class = "chrom down h130", style = "left: 75px;"),
                        div(class = "chrom up h110", style = "left: 90px;"),
                        div(class = "chrom down h70", style = "left: 105px;"),
                        div(class = "chrom up h130", style = "left: 120px;"),
                        div(class = "chrom down h50", style = "left: 135px;")
                    )
                )
            )
        )

      )


    })




    observeEvent(input$continueBtn, {
      show_landing(FALSE)
    })

    output$mainUI <- renderUI({
      req(user_session$logged_in)
      req(!show_landing())

      if (show_menu()) {
        return(uiOutput("dashboardMenuUI"))
      }

      fluidPage(
        tags$head(
          tags$style(HTML("body { background: linear-gradient(135deg, #1e3c72, #2a5298); color: white; }")),
          tags$style(HTML(".shiny-input-container { color: black; }"))
        ),

        fluidRow(
          column(9, titlePanel("MiamiR User Web App", windowTitle = "MiamiR Web App")),
          column(3, align = "right", br(), actionButton("goMenu", "Menu", class = "btn-warning", style = "margin-top: 10px;"))
        ),

        switch(selectedTab(),

               "Home" = tabPanel("Home",
                                 mainPanel(
                                   h2("Welcome to the MiamiR package Shiny App"),
                                   p("This application allows you to process data and generate different types of plots, including Manhattan, Miami and Forest Plots."),
                                   p("To get started, select one of the tabs above:"),
                                   tags$ul(
                                     tags$li("Single_Plot() for generating Manhattan Plots"),
                                     tags$li("Miami_Plot() for generating Miami Plots"),
                                     tags$li("Forest_Plot() for generating Forest Plots"),
                                     tags$li("Annotate_Data() for annotating data with RS codes for index SNPs"),
                                     tags$li("Model_Munge() for creating and munging LM and GLM models ahead of plotting")
                                   ),
                                   p("Click on a tab to get started!"),
                                   tags$img(src = "https://i.redd.it/l5sfjf60tclc1.gif", width = "500px", height = "auto")
                                 )
               ),

               "Single_Plot()" = tabPanel("Single_Plot()",
                                          style = "background: linear-gradient(135deg, #1e3c72, #2a5298); color: white;",
                                          sidebarLayout(
                                            sidebarPanel(
                                              fileInput("file1", "Upload Any File Type (up to 10GB)", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                              fluidRow(
                                                column(6, uiOutput("generateBtn1")),
                                                column(6, uiOutput("downloadBtn1"))
                                              ),
                                              tags$br(), tags$br(),
                                              uiOutput("dynamic_inputs_1")
                                            ),
                                            mainPanel(
                                              h4("Data Preview:"),
                                              tableOutput("dataPreview1"),
                                              h4("Rendered Plot:"),
                                              uiOutput("fullscreenButton1"),
                                              tags$div(
                                                id = "plotProgressContainer",
                                                style = "background-color: white; padding: 30px; text-align: center; font-family: 'Courier New', monospace;",
                                                tags$pre(id = "plotProgressStatus", " ")  # leave blank initially
                                              ),

                                              uiOutput("interactivePlot1_ui")
                                            )
                                          )
               ),

               "Regional_Plot()" = tabPanel("Regional_Plot()",
                                            style = "background: linear-gradient(135deg, #1e3c72, #2a5298); color: white;",
                                            sidebarLayout(
                                              sidebarPanel(
                                                fileInput("file6", "Upload Data (up to 10GB)", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                                numericInput("Chromosome", "Chromosome", value = 1, min = 1),
                                                fluidRow(
                                                  column(6, actionButton("submit6", "Generate Regional Plot")),
                                                  column(6, downloadButton("downloadPlot6", "Download All Plots"))
                                                ),
                                                tags$br(), tags$br(),
                                                uiOutput("dynamic_inputs_6")
                                              ),
                                              mainPanel(
                                                h4("Data Preview:"),
                                                tableOutput("dataPreview6"),
                                                h4("Rendered Regional Plots:"),
                                                uiOutput("renderedPlot6")
                                              )
                                            )
               ),

               "Miami_Plot()" = tabPanel("Miami_Plot()",
                                         style = "background: linear-gradient(135deg, #1e3c72, #2a5298); color: white;",
                                         sidebarLayout(
                                           sidebarPanel(
                                             fileInput("topData", "Upload Top Data (up to 10GB)", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                             fileInput("bottomData", "Upload Bottom Data (up to 10GB)", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                             fluidRow(
                                               column(6, actionButton("submit2", "Generate Miami Plot")),
                                               column(6, downloadButton("downloadPlot2", "Download Miami Plot"))
                                             ),
                                             tags$br(),
                                             uiOutput("dynamic_inputs_2")
                                           ),
                                           mainPanel(
                                             h4("Top Data Preview:"),
                                             tableOutput("topDataPreview"),
                                             h4("Bottom Data Preview:"),
                                             tableOutput("bottomDataPreview"),
                                             h4("Rendered Miami Plot:"),
                                             plotOutput("renderedPlot2")
                                           )
                                         )
               ),

               "Forest_Plot()" = tabPanel("Forest_Plot()",
                                          style = "background: linear-gradient(135deg, #1e3c72, #2a5298); color: white;",
                                          sidebarLayout(
                                            sidebarPanel(
                                              numericInput("numDatasets", "How many datasets do you want to upload?", value = 1, min = 1, max = 10),
                                              uiOutput("dynamic_uploads"),
                                              fluidRow(
                                                column(6, actionButton("submit3", "Generate Forest Plot")),
                                                column(6, downloadButton("downloadPlot3", "Download Forest Plot"))
                                              ),
                                              tags$br(), tags$br(),
                                              uiOutput("dynamic_inputs_3")
                                            ),
                                            mainPanel(
                                              h4("Data Preview:"),
                                              uiOutput("dataPreview3"),
                                              h4("Rendered Forest Plot:"),
                                              plotOutput("renderedPlot3")
                                            )
                                          )
               ),

               "Annotate_Data()" = tabPanel("Annotate_Data()",
                                            style = "background: linear-gradient(135deg, #1e3c72, #2a5298); color: white;",
                                            sidebarLayout(
                                              sidebarPanel(
                                                fileInput("file4", "Upload Data File (up to 10GB)", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                                fluidRow(
                                                  column(6, actionButton("submit4", "Annotate Data")),
                                                  column(6, downloadButton("downloadData", "Download Annotated Data"))
                                                ),
                                                tags$br(), tags$br(),
                                                uiOutput("dynamic_inputs_4")
                                              ),
                                              mainPanel(
                                                h4("Data Preview:"),
                                                tableOutput("dataPreview4"),
                                                h4("Annotated SNPs:"),
                                                tableOutput("annotatedSNPs")
                                              )
                                            )
               ),

               "Model_Munge()" = tabPanel("Model_Munge()",
                                          style = "background: linear-gradient(135deg, #1e3c72, #2a5298); color: white;",
                                          sidebarLayout(
                                            sidebarPanel(
                                              fileInput("file5", "Upload Data File (up to 10GB)", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                              textAreaInput("modelCode", "Write Your Model Code (e.g., lm(Var1 ~ Var2 + Var3))", value = "lm(Var1 ~ Var2 + Var3)", rows = 3),
                                              fluidRow(
                                                column(6, actionButton("submit5", "Run Model")),
                                                column(6, actionButton("mungeButton", "Munge Data"))
                                              ),
                                              tags$br(), tags$br(),
                                              uiOutput("dynamic_inputs_5"),
                                              downloadButton("downloadMunge", "Download Munged Data")
                                            ),
                                            mainPanel(
                                              h4("Data Preview:"),
                                              tableOutput("dataPreview5"),
                                              h4("Model Summary:"),
                                              verbatimTextOutput("modelSummary"),
                                              h4("Munged Data:"),
                                              tableOutput("mungedData")
                                            )
                                          )
               ),

               "Automated Batch Mode" = tabPanel("Automated Batch Mode",
                                                 style = "background: linear-gradient(135deg, #1e3c72, #2a5298); color: white;",
                                                 sidebarLayout(
                                                   sidebarPanel(
                                                     fileInput("batchData", "Upload Top Data (for Single_Plot, Regional_Plot, Miami_Plot)", accept = ".csv", multiple = TRUE),
                                                     conditionalPanel(
                                                       condition = "input.selectedFunctions.includes('Miami_Plot')",
                                                       fileInput("bottomBatchData", "Upload Bottom Data (only for Miami_Plot)", accept = ".csv")
                                                     ),
                                                     uiOutput("functionSelectionUI")

                                                     ,
                                                     actionButton("runBatch", "Run Selected Functions", class = "btn-success"),
                                                     checkboxInput("sendEmail", "Email me the results?", value = FALSE),
                                                     conditionalPanel(
                                                       condition = "input.sendEmail == true",
                                                       textInput("userEmail", "Enter your email address:", placeholder = "your@email.com")
                                                     ),
                                                     tags$br(),
                                                     downloadButton("downloadBatchPDF", "Download All Batch Plots as PDF", class = "btn-primary"),
                                                     tags$br(),  # Add spacing before this section

                                                     tags$div(style = "margin-top: 15px;",  # Push everything down a bit
                                                              fluidRow(
                                                                column(8,
                                                                       textInput("manualEmailInput", "Enter your email address:",
                                                                                 placeholder = "name@example.com", width = "100%")
                                                                ),
                                                                column(4,
                                                                       tags$br(),  # ⬅️ Aligns the button with input box
                                                                       actionButton("manualResendEmail", "📧 Resend Email",
                                                                                    class = "btn-info", style = "width: 100%;")
                                                                )
                                                              )
                                                     )
                                                   ),
                                                   mainPanel(
                                                     h4("Batch Execution Results:"),
                                                     uiOutput("batchResults")
                                                   )
                                                 )
               )
        ),
        tags$script(HTML("
    Shiny.addCustomMessageHandler('plot_progress', function(message) {
      const totalBars = 30;
      const filled = Math.round((message.pct / 100) * totalBars);
      const empty = totalBars - filled;

      const bar = '[' + '='.repeat(filled) + ' '.repeat(empty) + '] ' +
                  message.pct + '% (' + message.msg + ')';

      const el = document.getElementById('plotProgressStatus');
      if (el) {
        el.innerText = bar;
      }
    });
  "))


      )
    })


    # library(blastula)


    #  email_credentials <- blastula::creds_anonymous(
    #    host = "smtp.mail.yahoo.com",
    #    port = 465,
    #    username = "callum.oneill@rocketmail.com",
    #    password = "vyfnwlbnfqckumuh",
    #    use_ssl = TRUE
    #  )

    #email_credentials <- blastula::creds(
    #  host = "smtp.mail.yahoo.com",
    #  port = 465,
    #  user = "callum.oneill@rocketmail.com",
    #  pass = "vyfnwlbnfqckumuh",  # paste the 16-char app password here (no spaces)
    #  use_ssl = TRUE
    #  )




    batch_pdfs <- reactiveVal(list())

    # Helper function to build Single_Plot arguments for batch
    build_single_plot_args_batch <- function(batch_top_data, file_name = NULL) {
      args <- formals(Single_Plot)
      args <- lapply(args, function(x) if (is.call(x)) eval(x) else x)

      args$Data <- batch_top_data

      if (is.null(args$Title) || args$Title == "") {
        if (!is.null(file_name)) {
          args$Title <- tools::file_path_sans_ext(file_name)
        } else {
          args$Title <- "Auto Single Plot"
        }
      }

      if (is.null(args$Interactive) || args$Interactive == "") {
        args$Interactive <- FALSE
      }

      return(args)
    }

    # Helper function to build Regional_Plot arguments for batch
    safe_eval <- function(x) {
      if (missing(x)) return(NULL)
      if (is.symbol(x)) return(NULL)
      if (is.call(x)) return(eval(x))
      return(x)
    }

    build_regional_plot_args_batch <- function(batch_top_data) {
      # Combine defaults from both
      defaults_single <- formals(Single_Plot)
      defaults_regional <- formals(Regional_Plot)

      all_defaults <- modifyList(as.list(defaults_single), as.list(defaults_regional))
      all_defaults <- lapply(all_defaults, function(x) if (is.call(x)) eval(x) else x)

      # Combine user-saved settings from both tabs!
      saved_single <- batch_settings$Single_Plot %||% list()
      saved_regional <- batch_settings$Regional_Plot %||% list()
      all_saved <- modifyList(saved_single, saved_regional)

      args <- modifyList(all_defaults, all_saved)

      args$Data <- batch_top_data
      #  args$Chromosome <- args$Chromosome %||% 1
      args$Condense_Scale <- args$Condense_Scale %||% TRUE

      # Clean comma-separated strings
      for (arg in names(args)) {
        val <- args[[arg]]
        if (is.character(val) && length(val) == 1 && grepl(",", val)) {
          val_split <- trimws(strsplit(val, ",")[[1]])
          args[[arg]] <- if (all(!is.na(suppressWarnings(as.numeric(val_split))))) as.numeric(val_split) else val_split
        }
        if (is.character(val) && val %in% c("TRUE", "FALSE")) {
          args[[arg]] <- as.logical(val)
        }
      }

      return(args)
    }





    # Safe fallback operator
    `%||%` <- function(a, b) if (!is.null(a)) a else b

    build_single_plot_args_batch <- function(batch_top_data, file_name = NULL) {
      # 1. Start with function defaults
      defaults <- formals(Single_Plot)
      defaults <- lapply(defaults, function(x) if (is.call(x)) eval(x) else x)

      # 2. Merge in saved user settings (if any)
      saved <- batch_settings$Single_Plot %||% list()
      args <- modifyList(defaults, saved)

      # 3. Required overrides
      args$Data <- batch_top_data
      args$Title <- args$Title %||% tools::file_path_sans_ext(file_name %||% "Auto Single Plot")
      args$Interactive <- args$Interactive %||% FALSE

      # 4. Convert character lists (if needed)
      for (arg in c("Chromosome_Colours", "Break_Point", "Chromosome_Labels")) {
        if (!is.null(args[[arg]]) && is.character(args[[arg]])) {
          if (length(args[[arg]]) == 1 && grepl(",", args[[arg]])) {
            args[[arg]] <- trimws(strsplit(args[[arg]], ",")[[1]])
            if (arg == "Break_Point") args[[arg]] <- as.numeric(args[[arg]])
          }
        }
      }

      return(args)
    }


    build_regional_plot_args_batch <- function(batch_top_data) {
      args <- batch_settings$Regional_Plot %||% list()
      args$Data <- batch_top_data
      #  args$Chromosome <- args$Chromosome %||% 1
      args$Condense_Scale <- args$Condense_Scale %||% TRUE
      return(args)
    }

    build_miami_plot_args_batch <- function(top_data, bottom_data) {
      args <- batch_settings$Miami_Plot %||% list()
      args$Top_Data <- top_data
      args$Bottom_Data <- bottom_data
      return(args)
    }

    build_annotate_data_args_batch <- function(data) {
      args <- batch_settings$Annotate_Data %||% list()
      args$Data <- data
      return(args)
    }

    build_forest_plot_args_batch <- function(dataset_names) {
      args <- batch_settings$Forest_Plot %||% list()
      args$Data_Sets <- dataset_names
      return(args)
    }

    build_model_munge_args_batch <- function(model_obj_name = "Model") {
      args <- batch_settings$Model_Munge %||% list()
      args$Model_Object <- model_obj_name
      return(args)
    }



    # # Helper function to build Miami_Plot arguments for batch
    # build_miami_plot_args_batch <- function(batch_data) {
    #   midpoint <- floor(nrow(batch_data) / 2)
    #
    #   args <- formals(Miami_Plot)
    #   args <- lapply(args, function(x) if (is.call(x)) eval(x) else x)
    #
    #   args$Top_Data <- batch_data[1:midpoint, ]
    #   args$Bottom_Data <- batch_data[(midpoint+1):nrow(batch_data), ]
    #
    #   return(args)
    # }
    #
    # # Helper function to build Annotate_Data arguments for batch
    # build_annotate_data_args_batch <- function(batch_data) {
    #   args <- formals(Annotate_Data)
    #   args <- lapply(args, function(x) if (is.call(x)) eval(x) else x)
    #
    #   args$Data <- batch_data
    #
    #   return(args)
    # }


    render_interactive_plot <- function(plot_obj, width = NULL, height = NULL) {
      ggplotly(plot_obj, tooltip = "text") %>%
        layout(
          hoverlabel = list(font = list(size = 35)),
          autosize = is.null(width),  # FALSE only when width is explicitly set
          width = width,
          height = height,
          margin = list(t = 100, b = 50, l = 50, r = 150)
        ) %>%
        config(displayModeBar = TRUE, responsive = TRUE) %>%
        onRender("
  function(el, x) {
    let hideTimeout = null;

    if (!document.getElementById('snpMenu')) {
      var menu = document.createElement('div');
      menu.id = 'snpMenu';
      menu.style.position = 'absolute';
      menu.style.display = 'none';
      menu.style.zIndex = 10000;
      menu.style.background = 'rgba(255, 255, 255, 0.9)';
      menu.style.border = '1px solid rgba(204, 204, 204, 1)';
      menu.style.borderRadius = '5px';
      menu.style.boxShadow = '1px 1px 1px rgba(0, 0, 0, 0.2)';
      menu.style.fontFamily = 'Segoe UI, Helvetica, Arial, sans-serif';
      menu.style.fontSize = '12px';
      menu.style.color = 'rgb(68, 68, 68)';
      menu.style.padding = '6px 8px';
      menu.style.maxWidth = '150px';
      menu.style.pointerEvents = 'auto';
      menu.style.transform = 'translate(-100%, -100%)';

      menu.innerHTML = '<div style=\\\"font-weight:bold; margin-bottom:5px;\\\">View SNP in:</div>' +
                       '<div style=\\\"margin-bottom:4px;\\\"><button id=\\\"dbsnpBtn\\\" style=\\\"width:100%;\\\">dbSNP</button></div>' +
                       '<div style=\\\"margin-bottom:4px;\\\"><button id=\\\"otBtn\\\" style=\\\"width:100%;\\\">Open Targets</button></div>' +
                       '<div><button id=\\\"varsomeBtn\\\" style=\\\"width:100%;\\\">VarSome</button></div>';

      document.body.appendChild(menu);

      menu.addEventListener('mouseenter', function() {
        clearTimeout(hideTimeout);
        menu.setAttribute('data-hovering', 'true');
      });

      menu.addEventListener('mouseleave', function() {
        menu.setAttribute('data-hovering', 'false');
        menu.style.display = 'none';
      });
    }

    el.on('plotly_click', function(data) {
      const text = data.points[0].text.replace(/\\n/g, ' ').replace(/\\s+/g, ' ').trim();
      const snp_id = text.match(/rs\\d+/);
      const chrom = text.match(/CHR:\\s?(\\w+)/)?.[1];
      const pos = text.match(/POS:\\s?(\\d+)/)?.[1];
      const allele_line = text.match(/REF:\\s?([A-Za-z0-9]+)\\s?ALT:\\s?([A-Za-z0-9]+)/);
      const ref = allele_line?.[1];
      const alt = allele_line?.[2];

      const pointX = data.event.pageX;
      const pointY = data.event.pageY;
      const menu = document.getElementById('snpMenu');

      menu.style.left = pointX + 'px';
      menu.style.top = pointY + 'px';
      menu.style.display = 'block';

      document.getElementById('dbsnpBtn').onclick = function() {
        if (snp_id) window.open('https://www.ncbi.nlm.nih.gov/snp/' + snp_id[0], '_blank');
        menu.style.display = 'none';
      };
      document.getElementById('otBtn').onclick = function() {
        if (chrom && pos && ref && alt) {
          window.open('https://genetics.opentargets.org/variant/' + chrom + '_' + pos + '_' + ref + '_' + alt + '/associations', '_blank');
          menu.style.display = 'none';
        } else alert('Missing info for Open Targets link');
      };
      document.getElementById('varsomeBtn').onclick = function() {
        if (snp_id) {
          const gb = document.getElementById('Genome_Build')?.value.toLowerCase() || '';
          const build = gb.includes('38') ? 'hg38' : 'hg19';
          window.open('https://varsome.com/variant/' + build + '/' + snp_id[0], '_blank');
          menu.style.display = 'none';
        }
      };
    });

    el.on('plotly_unhover', function(data) {
      hideTimeout = setTimeout(function() {
        const menu = document.getElementById('snpMenu');
        if (menu.getAttribute('data-hovering') !== 'true') {
          menu.style.display = 'none';
        }
      }, 200);
    });

    document.addEventListener('click', function(e) {
      const menu = document.getElementById('snpMenu');
      if (!e.target.closest('#snpMenu') && !e.target.closest('.plotly')) {
        menu.style.display = 'none';
      }
    });

    console.log('🧪 onRender triggered in batch');
  }
  ")






    }



  run_with_counter <- function(func, args = list(), session = NULL) {
    exprs <- as.list(body(func))[-1]
    total <- length(exprs)
    env <- new.env(parent = environment(func))

    # Inject arguments and defaults
    defaults <- formals(func)
    for (name in names(defaults)) {
      assign(name, args[[name]] %||% eval(defaults[[name]], envir = env), envir = env)
    }

    result <- NULL
    for (i in seq_along(exprs)) {
      pct <- round(100 * i / total)
      bar_length <- 30
      filled_len <- round(bar_length * pct / 100)
      bar <- paste0("[", strrep("=", filled_len), strrep(" ", bar_length - filled_len), "] ", sprintf("%3d%%", pct))

      expr_text <- paste(deparse(exprs[[i]]), collapse = " ")
      is_vroom_line <- grepl("vroom", expr_text)

      if (!is_vroom_line) cat("\r", bar)

      # Evaluate each line and only assign to result if it's the last line
      if (i == length(exprs)) {
        result <- eval(exprs[[i]], envir = env)
      } else {
        # Suppress outputs for intermediate lines
        suppressMessages(suppressWarnings(
          capture.output(eval(exprs[[i]], envir = env), type = "output")
        ))
      }

      if (!is.null(session)) {
        session$sendCustomMessage("plot_progress", list(
          pct = pct,
          msg =  "Progress" # paste("Line", i, "of", total)
        ))
      }

      flush.console()
    }

    cat("\r[", strrep("=", 30), "] 100%\n", sep = "")
    return(result)  # ✅ return only once
  }






  # Insert this inside your observeEvent(input$runBatch, {...})
  # Assumes your Single_Plot, Regional_Plot, Annotate_Data etc. are defined
  # Insert this inside your observeEvent(input$runBatch, {...})
  # Assumes your Single_Plot, Regional_Plot, Annotate_Data etc. are defined

  # Insert this inside your observeEvent(input$runBatch, {...})
  # Assumes your Single_Plot, Regional_Plot, Annotate_Data etc. are defined
  # Insert this inside your observeEvent(input$runBatch, {...})
  # Assumes your Single_Plot, Regional_Plot, Annotate_Data etc. are defined
  # Insert this inside your observeEvent(input$runBatch, {...})
  # Assumes your Single_Plot, Regional_Plot, Annotate_Data etc. are defined

  # Insert this inside your observeEvent(input$runBatch, {...})
  # Assumes your Single_Plot, Regional_Plot, Annotate_Data etc. are defined

  # Insert this inside your observeEvent(input$runBatch, {...})
  # Assumes your Single_Plot, Regional_Plot, Annotate_Data etc. are defined

  # Insert this inside your observeEvent(input$runBatch, {...})
  # Assumes your Single_Plot, Regional_Plot, Annotate_Data etc. are defined

  # Insert this inside your observeEvent(input$runBatch, {...})
  # Assumes your Single_Plot, Regional_Plot, Annotate_Data etc. are defined

  # Insert this inside your observeEvent(input$runBatch, {...})
  # Assumes your Single_Plot, Regional_Plot, Annotate_Data etc. are defined

  # Safe wrapper to assign image output to avoid scoping issues
  assign_static_image <- function(id, path) {
    output[[id]] <- renderImage({
      list(src = path, contentType = "image/jpeg", width = "100%", height = "auto")
    }, deleteFile = FALSE)
  }

  # Inside observeEvent(input$runBatch, { ... })
  observeEvent(input$runBatch, {
    req(input$batchData)
    req(input$selectedFunctions)

    # Reset batch PDFs
    batch_pdfs(list())

    batch_top_files <- input$batchData
    bottom_file <- input$bottomBatchData

    all_names <- batch_top_files$name
    all_temp_paths <- vector("character", length(batch_top_files$datapath))

    for (i in seq_along(batch_top_files$datapath)) {
      data <- vroom::vroom(batch_top_files$datapath[i], show_col_types = FALSE)
      attr(data, "filename") <- all_names[i]
      temp_path <- tempfile(pattern = paste0("datafile_", tools::file_path_sans_ext(all_names[i])), fileext = ".rds")
      saveRDS(data, temp_path)
      all_temp_paths[i] <- temp_path
    }

    tab_panels <- lapply(seq_along(all_temp_paths), function(i) {
      this_data <- readRDS(all_temp_paths[[i]])
      file_name <- attr(this_data, "filename")
      panel_list <- list()

      # --- 🔹 Cover page PDF with filename ---
      title_pdf <- tempfile(fileext = ".pdf")
      grDevices::pdf(title_pdf, width = 11, height = 8.5)
      par(mar = c(0, 0, 0, 0))
      plot.new()
      text(0.5, 0.5, file_name, cex = 3, font = 2)
      dev.off()
      dataset_pdf_paths <- list(title_pdf)

      # --- 🔹 Single Plot ---
      # --- 🔹 Single Plot ---
      # --- 🔹 Single Plot ---
      # --- 🔹 Single Plot ---
      if ("Single_Plot" %in% input$selectedFunctions) {
        single_id <- paste0("single_plot_", i)
        plotly_id <- paste0("interactive_batch_plot_", i)
        args <- build_single_plot_args_batch(this_data, file_name = file_name)

        if (isTRUE(args$Interactive)) {
          single_plot <- run_with_counter(Single_Plot, args = args, session = session)
          output[[plotly_id]] <- renderPlotly({
            render_interactive_plot(single_plot)
          })

          panel_list <- append(panel_list, list(
            tabPanel("Single Plot (Interactive)",
                     tags$div(
                       style = "transform: scale(0.35); transform-origin: top left;",
                       plotlyOutput(plotly_id,
                                    width = paste0(30 * 96, "px"),
                                    height = paste0(15 * 96, "px"))
                     )
            )
          ))
        } else {
          static_img <- tempfile(fileext = ".jpg")
          plot <- run_with_counter(.Single_Plot_original, args = args, session = session)


          ggsave(static_img, plot = plot, width = 30, height = 15, units = "in", dpi = 300)
          single_pdf <- tempfile(fileext = ".pdf")
          ggsave(single_pdf, plot = plot, width = 30, height = 15, units = "in")
          dataset_pdf_paths <- c(dataset_pdf_paths, single_pdf)

          output[[single_id]] <- renderImage({
            list(src = static_img, contentType = "image/jpeg", width = "100%", height = "auto")
          }, deleteFile = FALSE)

          panel_list <- append(panel_list, list(tabPanel("Single Plot", imageOutput(single_id))))
        }
      }




      if ("Regional_Plot" %in% input$selectedFunctions) {
        regional_id <- paste0("batch_regional_", i)
        args <- build_regional_plot_args_batch(this_data)
        plots <- do.call(Regional_Plot, args)

        is_interactive <- isTRUE(args$Interactive)

        batch_regional_plots[[regional_id]] <- plots
        batch_regional_indices[[regional_id]] <- 1

        # Static version: save plots for PDF
        if (!is_interactive) {
          for (plot_name in names(plots)) {
            plot <- plots[[plot_name]]
            regional_pdf <- tempfile(fileext = ".pdf")
            height <- attr(plot, "dynamic_height", exact = TRUE)
            if (is.null(height)) height <- 25
            ggsave(regional_pdf, plot = plot, width = 30, height = height, units = "in")
            dataset_pdf_paths <- c(dataset_pdf_paths, regional_pdf)
          }
        }

        output[[paste0("batch_regional_ui_", regional_id)]] <- renderUI({
          idx <- batch_regional_indices[[regional_id]]
          plot_names <- names(batch_regional_plots[[regional_id]])
          plot <- batch_regional_plots[[regional_id]][[plot_names[idx]]]
          height <- attr(plot, "dynamic_height", exact = TRUE)
          if (is.null(height)) height <- 25

          # Output block for either static or interactive
          plot_output <- if (is_interactive) {
            tags$div(
              style = "transform: scale(0.35); transform-origin: top left;",
              plotlyOutput(paste0("batch_regional_interactive_", regional_id),
                           width = paste0(30 * 96, "px"),
                           height = paste0(height * 96, "px"))
            )
          } else {
            imageOutput(paste0("batch_regional_img_", regional_id))
          }

          tagList(
            tags$br(),
            fluidRow(
              column(12, align = "center",
                     div(style = "display: inline-block;",
                         actionButton(paste0("batch_regional_prev_", regional_id), "<< Previous")),
                     div(style = "display: inline-block; margin: 0 20px;",
                         tags$strong(paste(idx, "of", length(plot_names)))),
                     div(style = "display: inline-block;",
                         actionButton(paste0("batch_regional_next_", regional_id), "Next >>"))
              )
            ),
            tags$br(),
            plot_output
          )
        })

        # ⬅️ Render STATIC IMAGE
        output[[paste0("batch_regional_img_", regional_id)]] <- renderImage({
          idx <- batch_regional_indices[[regional_id]]
          plot_names <- names(batch_regional_plots[[regional_id]])
          plot <- batch_regional_plots[[regional_id]][[plot_names[idx]]]

          temp_path <- tempfile(fileext = ".jpg")
          height <- attr(plot, "dynamic_height", exact = TRUE)
          if (is.null(height)) height <- 25
          ggsave(temp_path, plot = plot, width = 30, height = height, units = "in", dpi = 300)

          list(src = temp_path, contentType = "image/jpeg", width = "100%", height = "auto")
        }, deleteFile = FALSE)

        # ⬅️ Render INTERACTIVE SUBPLOT (main + gene)
        output[[paste0("batch_regional_interactive_", regional_id)]] <- renderPlotly({
          idx <- batch_regional_indices[[regional_id]]
          plot_names <- names(batch_regional_plots[[regional_id]])
          current_plot <- batch_regional_plots[[regional_id]][[plot_names[idx]]]

          p_main <- attr(current_plot, "main_ggplot")
          p_genes <- attr(current_plot, "gene_track_plot")
          height <- attr(current_plot, "dynamic_height")
          if (is.null(height)) height <- 25

          fig1 <- ggplotly(p_main, tooltip = "text")
          fig2 <- ggplotly(p_genes, tooltip = "text")

          render_interactive_plot(
            subplot(fig1, fig2, nrows = 2, shareX = TRUE, heights = c(0.5, 0.5)),
            width = 30 * 96,
            height = height * 96
          )
        })

        # Pagination
        observeEvent(input[[paste0("batch_regional_prev_", regional_id)]], {
          idx <- batch_regional_indices[[regional_id]]
          new_idx <- ifelse(idx > 1, idx - 1, length(batch_regional_plots[[regional_id]]))
          batch_regional_indices[[regional_id]] <- new_idx
        })

        observeEvent(input[[paste0("batch_regional_next_", regional_id)]], {
          idx <- batch_regional_indices[[regional_id]]
          new_idx <- ifelse(idx < length(batch_regional_plots[[regional_id]]), idx + 1, 1)
          batch_regional_indices[[regional_id]] <- new_idx
        })

        # Add tab
        panel_list <- append(panel_list, list(tabPanel("Regional Plot", uiOutput(paste0("batch_regional_ui_", regional_id)))))
      }


      # --- 🔹 Miami Plot ---
      if ("Miami_Plot" %in% input$selectedFunctions && !is.null(bottom_file)) {
        bottom_data <- vroom::vroom(bottom_file$datapath, show_col_types = FALSE)
        # Get default args
        args <- formals(Miami_Plot)
        args <- lapply(args, function(x) if (is.call(x)) eval(x) else x)

        # Merge saved settings if any
        user_args <- batch_settings[["Miami_Plot"]]
        if (!is.null(user_args)) {
          for (arg in names(user_args)) {
            args[[arg]] <- user_args[[arg]]
          }
        }

        # Inject required data
        args$Top_Data <- this_data
        args$Bottom_Data <- bottom_data

        # Run plot
        plot <- do.call(Miami_Plot, args)


        static_img <- tempfile(fileext = ".jpg")
        ggsave(static_img, plot = plot, width = 30, height = 15, units = "in", dpi = 300)

        miami_pdf <- tempfile(fileext = ".pdf")
        ggsave(miami_pdf, plot = plot, width = 30, height = 15, units = "in")
        dataset_pdf_paths <- c(dataset_pdf_paths, miami_pdf)

        miami_id <- paste0("miami_plot_", i)
        output[[miami_id]] <- renderImage({
          list(src = static_img, contentType = "image/jpeg", width = "100%", height = "auto")
        }, deleteFile = FALSE)

        panel_list <- append(panel_list, list(tabPanel("Miami Plot", imageOutput(miami_id))))
      }

      # Save this file's PDF group to global batch_pdfs
      old_pdfs <- batch_pdfs()
      batch_pdfs(c(old_pdfs, list(dataset_pdf_paths)))

      # Return tabPanel for this file
      tabPanel(
        title = file_name,
        do.call(tabsetPanel, c(id = paste0("file_tabs_", i), panel_list))
      )
    })

    # Render UI with all file panels
    output$batchResults <- renderUI({
      do.call(tabsetPanel, tab_panels)  # 🔥 Show all results as separate tabs
    })

    folder_id <-  as_id("https://drive.google.com/drive/u/1/folders/184gvqSKKbfVy_fKKflWHVSDmiObdqIZL")

    # ✅ Combine all batch PDFs
    pdf_paths <- unname(unlist(batch_pdfs()))
    combined_path <- tempfile(fileext = ".pdf")

    if (length(pdf_paths) > 0) {
      pdftools::pdf_combine(pdf_paths, output = combined_path)

      # ✅ Upload to Google Drive
      unique_name <- paste0("Batch_Plots_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
      drive_file <- drive_upload(media = combined_path, path = folder_id, name = unique_name, type = "application/pdf")

      # ✅ Make public
      drive_share(drive_file, role = "reader", type = "anyone")
      public_link <- paste0("https://drive.google.com/uc?id=", drive_file$id, "&export=download")

      # ✅ Save for use in download and email
      options(latest_combined_pdf = combined_path)
      options(latest_batch_link = public_link)
    }


    # ✅ Send email if checked
    if (isTRUE(input$sendEmail)) {
      email <- input$userEmail
      public_link <- getOption("latest_batch_link", default = NULL)

      if (!is.null(public_link) && !is.null(email) && nzchar(email)) {
        tryCatch({
          library(mailR)
          send.mail(
            from = "callum.oneill@rocketmail.com",
            to = email,
            subject = "Your MiamiR Batch Results",
            body = paste0(
              "Hello,\n\nYour MiamiR batch results are ready.\n",
              "Download them here (valid for 24 hours):\n\n",
              public_link,
              "\n\nThank you."
            ),
            smtp = list(
              host.name = "smtp.mail.yahoo.com", port = 465,
              user.name = "callum.oneill@rocketmail.com",
              passwd = "vyfnwlbnfqckumuh", ssl = TRUE
            ),
            authenticate = TRUE,
            send = TRUE
          )
        }, error = function(e) {
          showModal(modalDialog("❌ Failed to send email: ", e$message, easyClose = TRUE))
        })
      }
    }
  })



  showModal(modalDialog(
    title = "✅ Batch Run Complete",
    "Your batch run is finished. You can now download the combined PDF or send it via email using the fields below.",
    easyClose = TRUE,
    footer = modalButton("OK")
  ))



  # Add this download handler in server()
  output$downloadBatchPDF <- downloadHandler(
    filename = function() {
      paste0("Batch_Plots_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      combined <- getOption("latest_combined_pdf", default = NULL)
      if (!is.null(combined) && file.exists(combined)) {
        file.copy(combined, file)
      } else {
        stop("❌ PDF not available. Please re-run the batch.")
      }
    }
  )





  observeEvent(input$manualResendEmail, {
    email <- input$userEmail
    public_link <- getOption("latest_batch_link", default = NULL)

    if (!nzchar(email)) {
      showModal(modalDialog("❌ Please enter an email address.", easyClose = TRUE))
      return()
    }

    if (is.null(public_link)) {
      showModal(modalDialog("❌ No batch PDF found to send. Please run a batch first.", easyClose = TRUE))
      return()
    }

    tryCatch({
      send.mail(
        from = "callum.oneill@rocketmail.com",
        to = email,
        subject = "Your MiamiR Batch Results",
        body = paste0(
          "Hello,\n\nHere is the latest batch results from the MiamiR app.\n",
          "Download your combined PDF here:\n\n",
          public_link,
          "\n\nThank you!"
        ),
        smtp = list(
          host.name = "smtp.mail.yahoo.com", port = 465,
          user.name = "callum.oneill@rocketmail.com",
          passwd = "vyfnwlbnfqckumuh", ssl = TRUE
        ),
        authenticate = TRUE,
        send = TRUE
      )

      showModal(modalDialog("✅ Email successfully sent!", easyClose = TRUE))
    }, error = function(e) {
      showModal(modalDialog("❌ Failed to send email: ", e$message, easyClose = TRUE))
    })
  })


  ### Single_Plot Logic ###
  # Load the uploaded data for Single_Plot
  # Global variable to store the path of the generated plot
  plot_width <- reactive({
    req(input$Width)
    input$Width * 96
  })

  plot_height <- reactive({
    req(input$Height)
    input$Height * 96
  })

  plot_file_1 <- reactiveVal(NULL)

  ### Single_Plot Logic ###
  # Load the uploaded data for Single_Plot

  user_inputs_cache <- reactiveVal(NULL)




  observe({
    req(input$file1)

    output$dataPreview1 <- renderTable({
      vroom(input$file1$datapath, n_max = 10)  # Preview first 10 rows
    })

    output$generateBtn1 <- renderUI({
      req(input$file1)
      actionButton("submit1", "Generate Manhattan Plot", class = "btn-primary")
    })

    output$downloadBtn1 <- renderUI({
      req(plot_file_1())  # Only show after plot file exists
      downloadButton("downloadPlot1", "Download Manhattan Plot")
    })

    output$fullscreenButton1 <- renderUI({
      req(input$submit1)  # Show only after "Generate" is clicked
      actionButton("fullscreenBtn1", "View Fullscreen", class = "btn-primary")
    })


    # Assign the uploaded data to the 'Data' argument
    full_data_1 <- reactive({
      req(input$file1)
      vroom(input$file1$datapath)  # Load the full dataset
    })

    # Dynamically generate UI inputs for Single_Plot
    output$dynamic_inputs_1 <- renderUI({
      generate_inputs_with_plot_settings(Single_Plot)
    })


    # Collect user inputs and dynamically pass them to Single_Plot, including Data
    # Collect user inputs and dynamically pass them to Single_Plot, including Data
    reactive_inputs_1 <- reactive({
      args <- formals(Single_Plot)

      user_inputs <- lapply(names(args), function(arg_name) {
        user_input <- input[[arg_name]]
        if (is.null(user_input) || user_input == "") return(NULL)
        if (grepl(",", user_input)) {
          input_value <- strsplit(user_input, ",")[[1]]
          return(trimws(input_value))
        } else {
          return(user_input)
        }
      })

      names(user_inputs) <- names(args)
      user_inputs$Data <- full_data_1()  # Assign the uploaded data to 'Data'

      # 🔥 Inject default title if not provided
      if (is.null(user_inputs$Title)) {
        if (!is.null(input$file1$name)) {
          user_inputs$Title <- tools::file_path_sans_ext(input$file1$name)
        }
      }

      return(user_inputs)
    })


    observeEvent(input$fullscreenBtn1, {
      showModal(modalDialog(
        easyClose = TRUE,
        footer = NULL,
        size = "l",

        tags$style(HTML("
          .modal-dialog {
            position: fixed;
            top: 0; left: 0;
            width: 100vw;
            height: 100vh;
            margin: 0;
            padding: 0;
            z-index: 1050;
          }
          .modal-content {
            height: 100vh;
            border: none;
            border-radius: 0;
            background-color: #fff;
            display: flex;
            flex-direction: column;
          }
          .modal-body {
            flex-grow: 1;
            overflow-y: hidden;
            display: flex;
            flex-direction: column;
          }
          #interactivePlot1_fullscreen {
            flex-grow: 1;
          }
        ")),

        # Exit button row
        fluidRow(
          column(4, align = "left", actionButton("exitFullscreen1", "Exit Fullscreen", class = "btn-danger")),
          column(4),
          column(4)
        ),
        tags$br(),

        # Fullscreen plot: same plotly code, just bigger!
        tags$div(
          style = "transform: scale(0.55); transform-origin: top left;",
          if (isTRUE(user_inputs_cache()$Interactive)) {
            plotlyOutput("interactivePlot1_fullscreen",
                         width = paste0(plot_width(), "px"),
                         height = paste0(plot_height(), "px"))
          } else {
            imageOutput("staticPlot1_fullscreen",
                        width = paste0(plot_width(), "px"),
                        height = paste0(plot_height(), "px"))
          }
        )
      ))
    })



    observeEvent(input$exitFullscreen1, {
      removeModal()
    })


    output$interactivePlot1_fullscreen <- renderPlotly({
      req(user_inputs_cache())
      req(user_inputs_cache()$Interactive)

      user_args <- isolate(reactive_inputs_1())
      user_args$Interactive <- TRUE
      plot_obj <- run_with_counter(Single_Plot, args = user_args, session = session)

      render_interactive_plot(plot_obj, width = input$Width * 96, height = input$Height * 96)
    })






    observeEvent(input$submit1, {

      output$plotStatusUI <- renderUI({
        req(input$submit1)  # Only show after "Generate Manhattan Plot" is clicked

        tags$div(
          id = "plotProgressContainer",
          style = "background-color: white; padding: 30px; text-align: center; font-family: 'Courier New', monospace;",
          tags$pre(id = "plotProgressStatus", "[                              ] 0% (Initializing...)")
        )
      })






      user_args <- isolate(reactive_inputs_1())
      user_inputs_cache(user_args)

      temp_file <- tempfile(fileext = paste0(".", input$File_Type))
      plot_file_1(temp_file)

      # Always generate the plot
      plot_obj <- run_with_counter(Single_Plot, args = user_args, session = session)

      # Save static version
      ggsave(temp_file, plot = plot_obj,
             width = input$Width, height = input$Height,
             units = "in", dpi = input$Quality)

      if (isTRUE(user_args$Interactive)) {
        output$interactivePlot1 <- renderPlotly({
          req(input$submit1)
          req(user_inputs_cache())
          req(user_inputs_cache()$Interactive)

          user_args <- isolate(reactive_inputs_1())
          user_args$Interactive <- TRUE
          plot_obj <- run_with_counter(Single_Plot, args = user_args, session = session)

          render_interactive_plot(plot_obj)
        })


        output$interactivePlot1_ui <- renderUI({
          tags$div(
            style = "transform: scale(0.35); transform-origin: top left;",
            plotlyOutput("interactivePlot1",
                         width = paste0(input$Width * 96, "px"),
                         height = paste0(input$Height * 96, "px"))
          )
        })
      } else {
        output$staticPlot1 <- renderImage({
          list(
            src = plot_file_1(),
            contentType = ifelse(input$File_Type == "png", "image/png", "image/jpeg"),
            width = "100%",   # Let the browser handle sizing
            height = "auto",
            alt = "Static plot"
          )
        }, deleteFile = FALSE)


        output$interactivePlot1_ui <- renderUI({
          imageOutput("staticPlot1",
                      width = "100%",
                      height = "auto")
        })
      }
    })



    output$interactivePlot1_ui <- renderUI({
      tags$div(
        style = "padding: 100px; text-align: center; font-size: 18px; color: #333;",
        "",
        tags$br(),
        tags$div(
          id = "plotProgressStatus",
          style = "font-weight: bold; font-size: 20px; color: #007bff;"
        )
      )
    })


    output$interactivePlot1 <- renderPlotly({
      req(input$submit1)
      req(user_inputs_cache())
      req(user_inputs_cache()$Interactive)

      user_args <- isolate(reactive_inputs_1())
      user_args$Interactive <- TRUE

      plot_obj <- run_with_counter(Single_Plot, args = user_args, session = session)

      ggplotly(plot_obj, tooltip = "text") %>%
        layout(
          autosize = TRUE,
          margin = list(t = 50, b = 50, l = 50, r = 50)
        ) %>%
        config(displayModeBar = FALSE)
    })



    output$staticPlot1_fullscreen <- renderImage({
      req(plot_file_1())

      list(
        src = plot_file_1(),
        contentType = ifelse(input$File_Type == "png", "image/png", "image/jpeg"),
        width = plot_width(),
        height = plot_height(),
        alt = "Fullscreen static plot"
      )
    }, deleteFile = FALSE)


    # output$interactivePlot1_ui <- renderUI({
    #   req(input$Width, input$Height)
    #   req(user_inputs_cache())
    #
    #   if (isTRUE(user_inputs_cache()$Interactive)) {
    #     plotlyOutput("interactivePlot1",
    #                  width = paste0(input$Width * 96, "px"),
    #                  height = paste0(input$Height * 96, "px"))
    #   } else {
    #     imageOutput("staticPlot1", width = paste0(input$Width * 96, "px"),
    #                 height = paste0(input$Height * 96, "px"))
    #   }
    # })


    #
    #       output$staticPlot1 <- renderImage({
    #         req(plot_file_1())
    #
    #         list(
    #           src = plot_file_1(),
    #           contentType = 'image/jpeg',
    #           width = input$Width * 96,
    #           height = input$Height * 96,
    #           alt = "Static plot"
    #         )
    #       }, deleteFile = FALSE)


    # Download handler for Single_Plot
    output$downloadPlot1 <- downloadHandler(
      filename = function() { paste0("Single_Plot_", Sys.Date(), ".", input$File_Type) },
      content = function(file) {
        if (!is.null(plot_file_1()) && input$File_Type != "pdf") {
          file.copy(plot_file_1(), file)
        } else {
          Plot_Outcome <- do.call(Single_Plot, user_inputs_cache())

          if (input$File_Type == "pdf") {
            # PDF-specific handling (no dpi here, default is vector)
            ggsave(file, plot = Plot_Outcome,
                   width = input$Width, height = input$Height,
                   units = "in")
          } else {
            ggsave(file, plot = Plot_Outcome,
                   width = input$Width, height = input$Height,
                   units = "in", dpi = input$Quality)
          }
        }

      }
    )
  })


  ###REIGONAL PLOT LOGICS


  plot_file_6 <- reactiveVal(NULL)
  regional_plot_paths <- reactiveVal(NULL)

  observe({
    req(input$file6)

    output$dataPreview6 <- renderTable({
      vroom(input$file6$datapath, n_max = 10)
    })

    output$dynamic_inputs_6 <- renderUI({
      generate_combined_inputs(Regional_Plot, Single_Plot)
    })

  })

  Regional_Plots_Results <- reactiveVal(NULL)

  observeEvent(input$submit6, {

    full_data_6 <- vroom(input$file6$datapath)

    fn_args_primary <- formals(Regional_Plot)
    fn_args_inherited <- formals(Single_Plot)

    all_args <- modifyList(as.list(fn_args_inherited), as.list(fn_args_primary))  # Inherit properly

    user_args <- lapply(names(all_args), function(arg_name) {
      if (arg_name %in% c("Data")) return(NULL)
      value <- input[[arg_name]]
      if (is.null(value) || value == "") return(NULL)
      if (grepl(",", value)) {
        return(trimws(strsplit(value, ",")[[1]]))
      } else {
        return(value)
      }
    })

    names(user_args) <- names(all_args)
    user_args$Data <- full_data_6

    # ✅ Generate Plots (this gives you named list of patchwork or ggplot objects)
    Regional_Plots <- do.call(Regional_Plot, user_args)
    Regional_Plots_Results(Regional_Plots)  # ✅ Store actual plot objects 🔥

    # Save Static Plots
    plot_dir <- file.path(tempdir(), paste0("regional_plots_", as.integer(Sys.time())))
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    saved_files <- list()
    for (plot_name in names(Regional_Plots)) {
      plot <- Regional_Plots[[plot_name]]
      safe_name <- gsub("[^A-Za-z0-9]", "_", plot_name)
      save_path <- file.path(plot_dir, paste0(safe_name, ".jpg"))

      height_to_use <- attr(plot, "dynamic_height")
      if (is.null(height_to_use)) height_to_use <- 25

      if (inherits(plot, "gg")) {
        ggsave(save_path, plot = plot, width = 30, height = height_to_use, units = "in", dpi = 300)
        saved_files[[plot_name]] <- save_path
      } else {
        # It's a plotly object. You CANNOT save it with ggsave.
        # You might optionally save it as HTML, or just skip.
        message(paste0("Skipping saving for interactive plot: ", plot_name))
      }
    }


    regional_plot_paths(saved_files)
  })


  # Track index of current plot
  current_plot_index <- reactiveVal(1)

  observeEvent(input$prevPlot, {
    new_index <- current_plot_index() - 1
    if (new_index < 1) new_index <- length(regional_plot_paths())
    current_plot_index(new_index)
  })

  observeEvent(input$nextPlot, {
    new_index <- current_plot_index() + 1
    if (new_index > length(regional_plot_paths())) new_index <- 1
    current_plot_index(new_index)
  })


  output$renderedPlot6 <- renderUI({
    req(Regional_Plots_Results())

    idx <- current_plot_index()  # ✅ NOT isolated
    plot_names <- names(Regional_Plots_Results())
    current_plot <- Regional_Plots_Results()[[plot_names[idx]]]

    if (isTRUE(input$Interactive)) {
      p_interactive <- attr(current_plot, "interactive_panel", exact = TRUE)
      if (is.null(p_interactive)) p_interactive <- current_plot

      height_in <- attr(current_plot, "dynamic_height", exact = TRUE)
      if (is.null(height_in)) height_in <- 25
      height_px <- height_in * 96

      tagList(
        tags$br(),
        fluidRow(
          column(4, align = "left", tagList(
            actionButton("prevPlot", "<< Previous"),
            actionButton("fullscreenBtn", "View Fullscreen", class = "btn-primary")
          )),
          column(4, align = "center", tags$h5(paste(idx, "of", length(plot_names)))),
          column(4, align = "right", actionButton("nextPlot", "Next >>"))
        ),
        tags$br(),
        tags$div(
          style = "transform: scale(0.35); transform-origin: top left;",
          plotlyOutput("interactivePlot6", height = paste0(height_px, "px"))
        )
      )
    }
    else {
      # --- STATIC HANDLING ---
      req(regional_plot_paths())  # Static version requires saved images
      plot_names <- names(regional_plot_paths())
      idx <- current_plot_index()

      tagList(
        tags$br(),
        fluidRow(
          column(4, align = "left", tagList(
            actionButton("prevPlot", "<< Previous"),
            actionButton("fullscreenBtn", "View Fullscreen", class = "btn-primary")
          )),
          column(4, align = "center", tags$h5(paste(idx, "of", length(plot_names)))),
          column(4, align = "right", actionButton("nextPlot", "Next >>"))
        ),
        tags$br(),

        imageOutput("regional_img")
      )
    }
  })




  output$regional_img <- renderImage({
    req(regional_plot_paths())
    plot_names <- names(regional_plot_paths())
    idx <- current_plot_index()

    list(
      src = regional_plot_paths()[[plot_names[idx]]],
      contentType = "image/jpeg",
      width = "100%",
      height = "auto",
      alt = "Static regional plot"
    )
  }, deleteFile = FALSE)



  output$inlineGeneTrackPlot6 <- renderPlot({
    req(Regional_Plots_Results())
    plot_names <- names(Regional_Plots_Results())
    idx <- current_plot_index()
    full_plot <- Regional_Plots_Results()[[plot_names[idx]]]

    p_gene <- attr(full_plot, "gene_track_panel", exact = TRUE)
    if (is.null(p_gene)) return(NULL)

    print(p_gene)
  })






  # 🔥 Place this directly below output$renderedPlot6
  output$interactivePlot6 <- renderPlotly({
    req(Regional_Plots_Results(), current_plot_index())

    plots <- Regional_Plots_Results()
    idx <- current_plot_index()
    plot_names <- names(plots)
    current_plot <- plots[[plot_names[idx]]]

    p_main <- attr(current_plot, "main_ggplot")
    p_genes <- attr(current_plot, "gene_track_plot")
    height_in <- attr(current_plot, "dynamic_height")
    if (is.null(height_in)) height_in <- 25

    fig1 <- ggplotly(p_main, tooltip = "text") %>%
      layout(margin = list(t = 50, b = 50, l = 50, r = 50)) %>%
      config(displayModeBar = TRUE, responsive = TRUE)

    fig2 <- ggplotly(p_genes, tooltip = "text") %>%
      layout(margin = list(t = 100, b = 50, l = 50, r = 50)) %>%
      config(displayModeBar = TRUE, responsive = TRUE)

    interactive_full <- subplot(
      fig1, fig2,
      nrows = 2,
      shareY = TRUE,
      shareX = FALSE,
      titleY = TRUE,
      margin = 0,
      heights = c(0.5, 0.5)
    ) %>%
      layout(
        hovermode = "closest",
        margin = list(t = 50, b = 50, l = 50, r = 50)
      ) %>%
      config(displayModeBar = TRUE, responsive = TRUE)

    render_interactive_plot(interactive_full, width = 30 * 96, height = height_in * 96)
  })




  output$modalCounter <- renderUI({
    req(regional_plot_paths())
    tags$h4(paste(current_plot_index(), "of", length(regional_plot_paths())))
  })

  observeEvent(input$prevPlotModal, {
    current_plot_index(ifelse(current_plot_index() > 1,
                              current_plot_index() - 1,
                              length(regional_plot_paths())))
  })

  observeEvent(input$nextPlotModal, {
    current_plot_index(ifelse(current_plot_index() < length(regional_plot_paths()),
                              current_plot_index() + 1,
                              1))
  })


  observeEvent(input$exitFullscreen6, {
    removeModal()
  })



  output$modalImage <- renderImage({
    req(regional_plot_paths())
    plot_names <- names(regional_plot_paths())
    idx <- current_plot_index()

    list(
      src = regional_plot_paths()[[plot_names[idx]]],
      contentType = "image/jpeg",
      width = "100%",
      height = "auto",
      alt = "Fullscreen static plot"
    )
  }, deleteFile = FALSE)

  observeEvent(input$exitFullscreen, {
    removeModal()
  })


  observeEvent(input$fullscreenBtn, {

    showModal(modalDialog(
      easyClose = TRUE,
      footer = NULL,
      size = "l",

      tags$style(HTML("
          .modal-dialog {
            position: fixed;
            top: 0; left: 0;
            width: 100vw;
            height: 100vh;
            margin: 0;
            padding: 0;
            z-index: 1050;
          }
          .modal-content {
            height: 100vh;
            border: none;
            border-radius: 0;
            background-color: #fff;
          }
          .modal-body {
            height: calc(100vh - 80px);
            overflow-y: auto;
          }
          .btn { margin-right: 8px; }  /* Optional spacing between buttons */
        ")),

      fluidRow(
        column(4, align = "left", tagList(
          actionButton("prevPlotModal", "<< Previous"),
          actionButton("exitFullscreen", "Exit Fullscreen", class = "btn-danger")
        )),
        column(4, align = "center", uiOutput("modalCounter")),
        column(4, align = "right", actionButton("nextPlotModal", "Next >>"))
      ),
      tags$br(),
      uiOutput("modalPlot6_ui")
    ))
  })





  output$modalPlot6_ui <- renderUI({
    req(Regional_Plots_Results())
    req(regional_plot_paths())

    idx <- current_plot_index()
    plot_names <- names(Regional_Plots_Results())
    full_plot <- Regional_Plots_Results()[[plot_names[idx]]]

    if (isTRUE(input$Interactive)) {
      # Extract dynamic dimensions
      height_in <- attr(full_plot, "dynamic_height")
      if (is.null(height_in)) height_in <- 25
      height_px <- height_in * 96

      p_gene <- attr(full_plot, "gene_track_panel", exact = TRUE)
      plot_width <- attr(p_gene, "plot_width", exact = TRUE)
      if (is.null(plot_width)) plot_width <- 30
      plot_width_px <- plot_width * 96

      tagList(
        tags$div(
          style = "transform: scale(0.55); transform-origin: top left;",
          plotlyOutput("modalInteractivePlot6", width = 30 * 96, height = height_px)
        ),
        tags$div(
          style = "margin-top: 40px;",
          plotOutput("modalGeneTrackPlot6",
                     width = paste0(plot_width_px, "px"),
                     height = paste0(height_px, "px"))
        )
      )
    } else {
      imageOutput("modalImage", width = "100%", height = "auto")
    }
  })





  output$modalGeneTrackPlot6 <- renderPlot({
    req(Regional_Plots_Results())
    plot_names <- names(Regional_Plots_Results())
    idx <- current_plot_index()
    full_plot <- Regional_Plots_Results()[[plot_names[idx]]]

    p_gene <- attr(full_plot, "gene_track_panel", exact = TRUE)
    if (is.null(p_gene)) return(NULL)

    plot_width <- attr(p_gene, "plot_width", exact = TRUE)
    if (is.null(plot_width)) plot_width <- 30
    plot_height <- attr(full_plot, "dynamic_height", exact = TRUE)
    if (is.null(plot_height)) plot_height <- 25

    print(p_gene)
  })





  output$downloadPlot6 <- downloadHandler(
    filename = function() {
      if (input$File_Type == "pdf") {
        paste0("Regional_Plots_", Sys.Date(), ".pdf")
      } else {
        paste0("Regional_Plots_", Sys.Date(), ".zip")
      }
    },
    content = function(file) {
      plots <- Regional_Plots_Results()
      req(plots)

      temp_dir <- tempfile()
      dir.create(temp_dir)
      saved_files <- c()
      height_list <- list()

      # Save plots as JPGs
      for (pname in names(plots)) {
        plot <- plots[[pname]]
        safe_name <- gsub("[^A-Za-z0-9]", "_", pname)
        save_path <- file.path(temp_dir, paste0(safe_name, ".jpg"))

        height_in <- attr(plot, "dynamic_height", exact = TRUE)
        if (is.null(height_in)) height_in <- 25  # Default

        width_in <- 30

        if (inherits(plot, "gg")) {
          ggsave(save_path, plot = plot, width = width_in, height = height_in, units = "in", dpi = input$Quality)
          saved_files <- c(saved_files, save_path)
          height_list[[save_path]] <- height_in  # record each plot's height
        } else {
          warning(paste("Skipping interactive plot", pname))
        }
      }

      if (input$File_Type == "pdf") {
        # --- New method: Create a multi-page PDF with different page sizes ---
        pdf_paths <- c()

        for (img_path in saved_files) {
          height_in <- height_list[[img_path]]
          width_in <- 30

          # Create a temporary PDF for each image with correct height
          single_pdf <- tempfile(fileext = ".pdf")

          grDevices::cairo_pdf(single_pdf, width = width_in, height = height_in)
          img <- jpeg::readJPEG(img_path)

          # --- Fixed Margin Version ---
          page_width <- width_in
          page_height <- height_in

          image_width <- unit((page_width - .5) / page_width, "npc")    # leave 1 inch each side
          image_height <- unit((page_height - .5) / page_height, "npc")

          grid::grid.raster(
            img,
            width = image_width,
            height = image_height,
            just = "center"
          )

          dev.off()

          pdf_paths <- c(pdf_paths, single_pdf)
        }



        # Merge all single PDFs into the final PDF
        pdftools::pdf_combine(pdf_paths, output = file)

      } else {
        # Zip the JPG files
        zip::zipr(zipfile = file, files = saved_files)
      }
    }
  )




  ### Miami_Plot Logic ###
  # Load and preview the uploaded data for Miami_Plot
  # Define reactive expressions for top_data and bottom_data globally so they are accessible
  top_data <- reactive({
    req(input$topData)  # Load the full top dataset after it's uploaded
    vroom(input$topData$datapath)
  })

  bottom_data <- reactive({
    req(input$bottomData)  # Load the full bottom dataset after it's uploaded
    vroom(input$bottomData$datapath)
  })

  # Load and preview the uploaded data for Miami_Plot
  # Global variable to store the path of the generated Miami plot
  plot_file_2 <- reactiveVal(NULL)

  ### Miami_Plot Logic ###
  # Load and preview the uploaded data for Miami_Plot
  observe({
    req(input$topData)  # Only require topData for previewing it

    output$topDataPreview <- renderTable({
      vroom(input$topData$datapath, n_max = 10)  # Preview top data after it's uploaded
    })
  })

  observe({
    req(input$bottomData)  # Only require bottomData for previewing it

    output$bottomDataPreview <- renderTable({
      vroom(input$bottomData$datapath, n_max = 10)  # Preview bottom data after it's uploaded
    })
  })

  # Dynamically generate UI inputs for Miami_Plot
  output$dynamic_inputs_2 <- renderUI({
    req(input$topData, input$bottomData)

    inputs <- generate_inputs_with_defaults(Miami_Plot)
    saved <- batch_settings$Miami_Plot

    if (!is.null(saved)) {
      for (i in seq_along(inputs)) {
        this_input <- inputs[[i]]
        this_name <- this_input$name
        if (!is.null(this_name) && this_name %in% names(saved)) {
          inputs[[i]]$value <- saved[[this_name]]
        }
      }
    }

    tagList(inputs)
  })


  # Collect user inputs and dynamically pass them to Miami_Plot
  reactive_inputs_2 <- reactive({
    args <- formals(Miami_Plot)

    user_inputs <- lapply(names(args), function(arg_name) {
      user_input <- input[[arg_name]]
      if (is.null(user_input) || user_input == "") return(NULL)
      if (grepl(",", user_input)) {
        input_value <- strsplit(user_input, ",")[[1]]
        return(trimws(input_value))
      } else {
        return(user_input)
      }
    })

    names(user_inputs) <- names(args)
    user_inputs$Top_Data <- top_data()  # Assign top data
    user_inputs$Bottom_Data <- bottom_data()  # Assign bottom data
    return(user_inputs)
  })

  # Save and display the Miami_Plot plot as an image
  output$renderedPlot2 <- renderImage({
    req(input$submit2)

    # Generate a temporary file to store the saved plot
    temp_file <- tempfile(fileext = paste0(".", input$File_Type))
    plot_file_2(temp_file)  # Store the temp file globally

    # Collect user inputs and dynamically pass them to Miami_Plot
    user_args <- reactive_inputs_2()

    # Create the plot and save it to the temporary file
    Plot_Outcome <- do.call(Miami_Plot, user_args)
    ggplot2::ggsave(temp_file, plot = Plot_Outcome, width = input$Width,
                    height = input$Height, units = "in", dpi = input$Quality)

    # Return the image to be displayed in the UI
    list(src = temp_file, contentType = paste0("image/", input$File_Type),
         width = "100%", height = "auto")
  }, deleteFile = FALSE)  # Do not delete the file after rendering

  # Download handler for Miami_Plot
  output$downloadPlot2 <- downloadHandler(
    filename = function() { paste0("Miami_Plot_", Sys.Date(), ".", input$File_Type) },
    content = function(file) {
      # If the plot file already exists, use it
      if (!is.null(plot_file_2())) {
        file.copy(plot_file_2(), file)
      } else {
        # Otherwise, regenerate the plot
        user_args <- reactive_inputs_2()
        Plot_Outcome <- do.call(Miami_Plot, user_args)
        ggplot2::ggsave(file, plot = Plot_Outcome, width = input$Width,
                        height = input$Height, units = "in", dpi = input$Quality)
      }
    }
  )


  ### Forest_Plot Logic ###
  # Dynamically generate the file upload inputs for Forest_Plot based on the number of datasets
  # Global variable to store the path of the generated Forest plot
  # Global variable to store the path of the generated Forest plot
  plot_file_3 <- reactiveVal(NULL)

  ### Forest_Plot Logic ###
  # Dynamically generate the file upload inputs for Forest_Plot based on the number of datasets
  output$dynamic_uploads <- renderUI({
    req(input$numDatasets)

    # Create file input elements based on the number of datasets, using unique IDs for Forest_Plot
    lapply(1:input$numDatasets, function(i) {
      fileInput(paste0("forest_file", i), paste("Upload Dataset", i, "(up to 10GB)"),
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
    })
  })

  # Dynamically generate data previews for the uploaded datasets, using unique IDs for Forest_Plot
  output$dataPreview3 <- renderUI({
    req(input$numDatasets)

    # Create table output elements for each dataset preview
    lapply(1:input$numDatasets, function(i) {
      tableOutput(paste0("forestDataPreview", i))
    })
  })

  # Load and preview each dataset based on user input, using unique IDs for Forest_Plot
  #  observe({
  #    lapply(1:input$numDatasets, function(i) {
  #      local({
  #        datasetInput <- paste0("forest_file", i)  # Updated input IDs for Forest_Plot
  #        previewOutput <- paste0("forestDataPreview", i)  # Updated preview IDs for Forest_Plot
  #
  #       output[[previewOutput]] <- renderTable({
  #        req(input[[datasetInput]])  # Only show preview for Forest_Plot uploads
  #       vroom(input[[datasetInput]]$datapath, n_max = 10)  # Preview first 10 rows of the uploaded data
  #     })
  #    })
  #    })
  #  })


  observe({
    req(input$numDatasets)  # Prevent this from running with NULL

    lapply(1:input$numDatasets, function(i) {
      local({
        datasetInput <- paste0("forest_file", i)
        previewOutput <- paste0("forestDataPreview", i)

        output[[previewOutput]] <- renderTable({
          req(input[[datasetInput]])
          vroom(input[[datasetInput]]$datapath, n_max = 10)
        })
      })
    })
  })


  # Load full datasets into the environment and pass them as a list to Forest_Plot
  full_datasets <- reactive({
    req(input$numDatasets)

    dataset_names <- lapply(1:input$numDatasets, function(i) {
      datasetInput <- paste0("forest_file", i)  # Updated to use Forest_Plot file inputs
      req(input[[datasetInput]])  # Ensure this dataset was uploaded for Forest_PPlot
      data <- vroom(input[[datasetInput]]$datapath)  # Load full dataset

      # Create a temporary variable name
      dataset_name <- paste0("forest_dataset", i)  # Updated variable names for Forest_Plot
      assign(dataset_name, data, envir = .GlobalEnv)

      return(dataset_name)
    })

    return(unlist(dataset_names))  # Return the names of the datasets
  })

  # Dynamically generate UI inputs for Forest_Plot (excluding Data_Sets)
  output$dynamic_inputs_3 <- renderUI({
    req(input$numDatasets)

    # Ensure all datasets have been uploaded before generating UI inputs
    all_files_uploaded <- all(sapply(1:input$numDatasets, function(i) {
      !is.null(input[[paste0("forest_file", i)]])  # Check for Forest_Plot file uploads
    }))

    req(all_files_uploaded)  # Only proceed if all files are uploaded

    generate_inputs_without_datasets(Forest_Plot)  # Generate inputs for Forest_Plot
  })

  # Helper function to generate dynamic inputs excluding 'Data_Sets'
  generate_inputs_without_datasets <- function(func) {
    arg_list <- formals(func)

    inputs <- lapply(names(arg_list), function(arg_name) {
      if (arg_name == "Data_Sets") return(NULL)  # Exclude 'Data_Sets' argument from UI
      arg_value <- arg_list[[arg_name]]
      if (is.call(arg_value)) arg_value <- eval(arg_value)

      if (is.numeric(arg_value) && length(arg_value) > 1) {
        textInput(arg_name, label = arg_name, value = paste(arg_value, collapse = ","))
      } else if (is.character(arg_value) && length(arg_value) > 1) {
        textInput(arg_name, label = arg_name, value = paste(arg_value, collapse = ","))
      } else if (is.numeric(arg_value)) {
        numericInput(arg_name, label = arg_name, value = arg_value)
      } else if (is.character(arg_value)) {
        textInput(arg_name, label = arg_name, value = arg_value)
      } else if (is.logical(arg_value)) {
        checkboxInput(arg_name, label = arg_name, value = arg_value)
      } else {
        textInput(arg_name, label = arg_name, value = as.character(arg_value))
      }
    })

    return(Filter(Negate(is.null), inputs))  # Return all inputs excluding Data_Sets
  }


  generate_inputs_with_plot_settings <- function(func, saved_values = list()) {
    tagList(
      textInput("File_Type", "File Type", value = "jpg"),
      numericInput("Width", "Plot Width (in)", value = 30),
      numericInput("Height", "Plot Height (in)", value = 15),
      numericInput("Quality", "Plot DPI", value = 300),
      generate_inputs_with_defaults(func, saved_values)
    )
  }


  # Collect user inputs and dynamically pass them to Forest_Plot along with dataset names
  reactive_inputs_3 <- reactive({
    req(input$numDatasets)
    args <- formals(Forest_Plot)

    user_inputs <- lapply(names(args), function(arg_name) {
      if (arg_name == "Data_Sets") return(NULL)  # Skip Data_Sets argument
      user_input <- input[[arg_name]]
      if (is.null(user_input) || user_input == "") return(NULL)
      if (grepl(",", user_input)) {
        input_value <- strsplit(user_input, ",")[[1]]
        input_value <- trimws(input_value)
        if (all(!is.na(as.numeric(input_value)))) return(as.numeric(input_value))
        return(as.character(input_value))
      } else {
        return(user_input)
      }
    })

    names(user_inputs) <- names(args)
    user_inputs$Data_Sets <- full_datasets()  # Assign the list of dataset names to Data_Sets
    return(user_inputs)
  })

  # Save and display the Forest_Plot as an image
  output$renderedPlot3 <- renderImage({
    req(input$numDatasets)
    req(input$submit3)

    # Generate a temporary file to store the saved plot
    temp_file <- tempfile(fileext = paste0(".", input$File_Type))
    plot_file_3(temp_file)  # Store the temp file globally

    # Collect user inputs and dynamically pass them to Forest_Plot
    user_args <- reactive_inputs_3()

    # Create the plot and save it to the temporary file
    Plot_Outcome <- do.call(Forest_Plot, user_args)
    ggplot2::ggsave(temp_file, plot = Plot_Outcome, width = input$Width,
                    height = input$Height, units = "in", dpi = input$Quality)

    # Return the image to be displayed in the UI
    list(src = temp_file, contentType = paste0("image/", input$File_Type),
         width = "100%", height = "auto")
  }, deleteFile = FALSE)  # Do not delete the file after rendering

  # Download handler for Forest_Plot
  output$downloadPlot3 <- downloadHandler(
    filename = function() { paste0("Forest_Plot_", Sys.Date(), ".", input$File_Type) },
    content = function(file) {
      # If the plot file already exists, use it
      if (!is.null(plot_file_3())) {
        file.copy(plot_file_3(), file)
      } else {
        # Otherwise, regenerate the plot
        user_args <- reactive_inputs_3()
        Plot_Outcome <- do.call(Forest_Plot, user_args)
        ggplot2::ggsave(file, plot = Plot_Outcome, width = input$Width,
                        height = input$Height, units = "in", dpi = input$Quality)
      }
    }
  )


  ### Annotate_Data Logic ###
  ### Annotate_Data Logic ###
  # Global variable to store the path of the generated annotated data
  annotated_file <- reactiveVal(NULL)

  ### Annotate_Data Logic ###
  ### Annotate_Data Logic ###
  # Observe file upload and generate dynamic UI
  observe({
    req(input$file4)

    # Preview the uploaded data
    output$dataPreview4 <- renderTable({
      vroom(input$file4$datapath, n_max = 10)  # Preview first 10 rows of the original dataset
    })

    full_data_4 <- reactive({
      req(input$file4)
      vroom(input$file4$datapath)  # Load the full dataset
    })

    # Dynamically generate UI inputs for Annotate_Data
    output$dynamic_inputs_4 <- renderUI({
      req(full_data_4())  # Ensure data is loaded before generating inputs
      generate_inputs_with_defaults(Annotate_Data)  # Generate dynamic inputs for Annotate_Data
    })

    # Collect user inputs and dynamically pass them to Annotate_Data, including Data
    reactive_inputs_4 <- reactive({
      args <- formals(Annotate_Data)

      user_inputs <- lapply(names(args), function(arg_name) {
        user_input <- input[[arg_name]]
        if (is.null(user_input) || user_input == "") return(NULL)
        if (grepl(",", user_input)) {
          input_value <- strsplit(user_input, ",")[[1]]
          return(trimws(input_value))
        } else {
          return(user_input)
        }
      })

      names(user_inputs) <- names(args)
      user_inputs$Data <- full_data_4()  # Assign the uploaded data to 'Data'
      return(user_inputs)
    })

    # Reactive expression for filtered data (containing "rs" in Lab column)
    filtered_data <- reactive({
      req(input$submit4)  # Ensure the button was clicked to annotate
      user_args <- reactive_inputs_4()
      Annotated_Data <- do.call(Annotate_Data, user_args)

      # Filter rows where 'Lab' contains "rs"
      Annotated_Data %>% filter(grepl("rs", Lab, ignore.case = TRUE))
    })

    # Display filtered data after annotation
    output$annotatedSNPs <- renderTable({
      req(filtered_data())  # Only proceed if the filtered data exists
      head(filtered_data(), 100)  # Show first 100 rows that match the condition
    })

    # Store the result in a temporary file for download
    observeEvent(input$submit4, {
      temp_file <- tempfile(fileext = ".csv")
      write.csv(filtered_data(), temp_file, row.names = FALSE)
      annotated_file(temp_file)  # Store the temp file path in reactive value
    })

    # Download handler for Annotate_Data
    output$downloadData <- downloadHandler(
      filename = function() { paste0("Annotated_Data_", Sys.Date(), ".csv") },
      content = function(file) {
        # If the annotated file already exists, use it
        if (!is.null(annotated_file())) {
          file.copy(annotated_file(), file)
        } else {
          # Otherwise, regenerate the filtered annotated data and save it
          write.csv(filtered_data(), file, row.names = FALSE)
        }
      }
    )
  })


  ### Model_Munge Logic ###
  observe({
    req(input$file5)

    # Preview the uploaded data
    output$dataPreview5 <- renderTable({
      vroom(input$file5$datapath, n_max = 10)  # Preview first 10 rows of the original dataset
    })

    full_data_5 <- reactive({
      req(input$file5)
      vroom(input$file5$datapath)  # Load the full dataset
    })

    # Dynamically generate UI inputs for Model_Munge
    output$dynamic_inputs_5 <- renderUI({
      generate_inputs_with_defaults(Model_Munge)  # Generate dynamic inputs for Model_Munge
    })

    # Collect user inputs and dynamically pass them to Model_Munge
    reactive_inputs_5 <- reactive({
      args <- formals(Model_Munge)

      user_inputs <- lapply(names(args), function(arg_name) {
        user_input <- input[[arg_name]]
        if (is.null(user_input) || user_input == "") return(NULL)
        if (grepl(",", user_input)) {
          input_value <- strsplit(user_input, ",")[[1]]
          return(trimws(input_value))
        } else {
          return(user_input)
        }
      })

      names(user_inputs) <- names(args)
      return(user_inputs)
    })

    # When the user clicks "Run Model"
    observeEvent(input$submit5, {
      req(input$modelCode)

      # Dynamically evaluate the model code input
      model_expr <- parse(text = input$modelCode)

      # Evaluate the model formula using the uploaded dataset
      Model <- eval(model_expr, envir = full_data_5())  # Create the model object

      # Assign the model object to "Model" in the global environment
      assign("Model", Model, envir = .GlobalEnv)

      # Display the model summary
      output$modelSummary <- renderPrint({
        summary(Model)
      })
    })

    # When the user clicks "Munge"
    observeEvent(input$mungeButton, {
      output$mungedData <- renderTable({
        # Use the model object name "Model"
        req(exists("Model", envir = .GlobalEnv))

        # Prepare inputs for the Model_Munge function
        user_args <- reactive_inputs_5()  # Collect the inputs dynamically
        user_args$Model_Object <- "Model"  # Pass the model object name as a string

        # Call the Model_Munge function and preview the munged data
        Munged_Data <- do.call(Model_Munge, user_args)
        return(head(Munged_Data, 10))  # Preview first 10 rows of munged data
      })
    })

    # Download handler for munged data
    output$downloadMunge <- downloadHandler(
      filename = function() { paste0("Munged_Data_", Sys.Date(), ".csv") },
      content = function(file) {
        # Ensure the Model object exists before munging
        req(exists("Model", envir = .GlobalEnv))

        # Prepare inputs for the Model_Munge function
        user_args <- reactive_inputs_5()  # Collect the inputs dynamically
        user_args$Model_Object <- "Model"  # Pass the model object name as a string

        # Call the Model_Munge function and save the munged data to a CSV file
        Munged_Data <- do.call(Model_Munge, user_args)
        write.csv(Munged_Data, file, row.names = FALSE)
      }
    )
  })
  }

  # Run the Shiny app
  shinyApp(ui = ui, server = server)
