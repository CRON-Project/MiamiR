  library(shiny)
  library(ggplot2)
  library(vroom)
  library(dplyr)

  # Increase file upload limit to 10 GB
  options(shiny.maxRequestSize = 10 * 1024^3)
  options(shiny.http.timeout = 300)  # Timeout increased to 5 minutes
  options(shiny.maxRequestSize = 10 * 1024^3)  # 30 MB


  # Helper function to generate dynamic UI inputs, specifically handling list inputs
  generate_inputs <- function(func) {
    arg_list <- formals(func)

    inputs <- lapply(names(arg_list), function(arg_name) {
      if (arg_name == "Data" | arg_name == "Top_Data" | arg_name == "Bottom_Data") return(NULL)  # Skip the Data argument

      arg_value <- arg_list[[arg_name]]

      # Automatically handle vector or list arguments as comma-separated inputs
      if (is.call(arg_value)) arg_value <- eval(arg_value)

      if (is.numeric(arg_value) && length(arg_value) > 1) {
        # Convert numeric vector to comma-separated string
        textInput(arg_name, label = arg_name, value = paste(arg_value, collapse = ","))

      } else if (is.character(arg_value) && length(arg_value) > 1) {
        # Convert character vector to comma-separated string
        textInput(arg_name, label = arg_name, value = paste(arg_value, collapse = ","))

      } else if (is.numeric(arg_value)) {
        # Single numeric input
        numericInput(arg_name, label = arg_name, value = arg_value)

      } else if (is.character(arg_value)) {
        # Single text input
        textInput(arg_name, label = arg_name, value = arg_value)

      } else if (is.logical(arg_value)) {
        # Checkbox input for booleans
        checkboxInput(arg_name, label = arg_name, value = arg_value)

      } else {
        # Fallback for unknown types
        textInput(arg_name, label = arg_name, value = as.character(arg_value))
      }
    })

    # Return the inputs, filtering out any NULL values (i.e., excluding 'Data')
    return(Filter(Negate(is.null), inputs))
  }


  # Define UI for the Shiny app
  # Define UI for the Shiny app
  # Define UI for the Shiny app
  # Define UI for the Shiny app
  # Define UI for the Shiny app
  ui <- fluidPage(
    titlePanel("MiamiR User Web App"),

    tabsetPanel(

      tabPanel("Home",  # New Home Tab
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
              #   tags$img(src = "2024-10-11-MiamiR", width = "500px", height = "auto")  # Add the gif stored in the www folder

                )



      ),


      tabPanel("Single_Plot()",  # First Tab for Single_Plot
               sidebarLayout(
                 sidebarPanel(
                   fileInput("file1", "Upload Any File Type (up to 10GB)",
                             accept = c("text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".csv")),

                   # Place Generate and Download buttons right after the file upload
                   fluidRow(
                     column(6, actionButton("submit1", "Generate Manhattan Plot")),
                     column(6, downloadButton("downloadPlot1", "Download Manhattan Plot"))
                   ),
                   tags$br(), tags$br(),  # Add some spacing below the buttons

                   uiOutput("dynamic_inputs_1")  # Dynamically generated inputs for first plot
                 ),

                 mainPanel(
                   h4("Data Preview:"),
                   tableOutput("dataPreview1"),  # Data preview for first plot

                   h4("Rendered Plot:"),
                   plotOutput("renderedPlot1")  # Directly render the ggplot object for first plot
                 )
               )
      ),

      tabPanel("Miami_Plot()",  # Second Tab for Miami_Plot
               sidebarLayout(
                 sidebarPanel(
                   fileInput("topData", "Upload Top Data (up to 10GB)",
                             accept = c("text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".csv")),
                   fileInput("bottomData", "Upload Bottom Data (up to 10GB)",
                             accept = c("text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".csv")),

                   # Wrap Generate and Download buttons in a fluidRow with two columns for side-by-side layout
                   fluidRow(
                     column(6, actionButton("submit2", "Generate Miami Plot")),
                     column(6, downloadButton("downloadPlot2", "Download Miami Plot"))
                   ),
                   tags$br(),  # Add some spacing below the buttons

                   uiOutput("dynamic_inputs_2")  # Dynamically generated inputs for Miami_Plot
                 ),

                 mainPanel(
                   h4("Top Data Preview:"),
                   tableOutput("topDataPreview"),  # Data preview for Top_Data

                   h4("Bottom Data Preview:"),
                   tableOutput("bottomDataPreview"),  # Data preview for Bottom_Data

                   h4("Rendered Miami Plot:") ,
                  plotOutput("renderedPlot2")  # Directly render the ggplot object for Miami_Plot
                 )
               )
      ),

      tabPanel("Forest_Plot()",  # Third Tab for Forest_Plot
               sidebarLayout(
                 sidebarPanel(
                   numericInput("numDatasets", "How many datasets do you want to upload?",
                                value = 1, min = 1, max = 10),

                   uiOutput("dynamic_uploads"),  # Dynamically generated file inputs for datasets

                   # Place Generate and Download buttons right after the file uploads
                   fluidRow(
                     column(6, actionButton("submit3", "Generate Forest Plot")),
                     column(6, downloadButton("downloadPlot3", "Download Forest Plot"))
                   ),
                   tags$br(), tags$br(),  # Add some spacing below the buttons

                   uiOutput("dynamic_inputs_3")  # Dynamically generated inputs for Forest_Plot
                 ),

                 mainPanel(
                   h4("Data Preview:"),
                   uiOutput("dataPreview3"),  # Data preview for Forest_Plot (dynamically generated)

                   h4("Rendered Forest Plot:"),
                   plotOutput("renderedPlot3")  # Directly render the ggplot object for Forest_Plot
                 )
               )
      ),

      tabPanel("Annotate_Data()",  # New Tab for Annotate_Data
               sidebarLayout(
                 sidebarPanel(
                   fileInput("file4", "Upload Data File (up to 10GB)",
                             accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),

                   # Place a Download button after file upload
                   fluidRow(
                     column(6, actionButton("submit4", "Annotate Data")),
                     column(6, downloadButton("downloadData", "Download Annotated Data"))
                   ),
                   tags$br(), tags$br(),  # Add some spacing below the buttons

                   uiOutput("dynamic_inputs_4")  # Dynamically generated inputs for Annotate_Data
                 ),

                 mainPanel(
                   h4("Data Preview:"),
                   tableOutput("dataPreview4"),  # Data preview for Annotate_Data

                   h4("Annotated SNPs:"),
                   tableOutput("annotatedSNPs")  # Display the filtered rows where 'Lab' is not blank
                 )
               )
      ),


      tabPanel("Model_Munge()",
               sidebarLayout(
                 sidebarPanel(
                   fileInput("file5", "Upload Data File (up to 10GB)",
                             accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),

                   textAreaInput("modelCode", "Write Your Model Code (e.g., lm(Var1 ~ Var2 + Var3))",
                                 value = "lm(Var1 ~ Var2 + Var3)", rows = 3),

                   # Buttons to run the model and munge
                   fluidRow(
                     column(6, actionButton("submit5", "Run Model")),
                     column(6, actionButton("mungeButton", "Munge Data"))
                   ),
                   tags$br(), tags$br(),

                   uiOutput("dynamic_inputs_5"),  # Dynamically generated inputs for Model_Munge

                   # Download button for the munged data
                   downloadButton("downloadMunge", "Download Munged Data")
                 ),

                 mainPanel(
                   h4("Data Preview:"),
                   tableOutput("dataPreview5"),  # Data preview for uploaded dataset

                   h4("Model Summary:"),
                   verbatimTextOutput("modelSummary"),  # Display the model summary

                   h4("Munged Data:"),
                   tableOutput("mungedData")  # Display the munged data
                 )
               )
      )
    )
  )

  server <- function(input, output, session) {

    ### Single_Plot Logic ###
    # Load the uploaded data for Single_Plot
    # Global variable to store the path of the generated plot
    plot_file_1 <- reactiveVal(NULL)

    ### Single_Plot Logic ###
    # Load the uploaded data for Single_Plot
    observe({
      req(input$file1)

      output$dataPreview1 <- renderTable({
        vroom(input$file1$datapath, n_max = 10)  # Preview first 10 rows
      })

      # Assign the uploaded data to the 'Data' argument
      full_data_1 <- reactive({
        req(input$file1)
        vroom(input$file1$datapath)  # Load the full dataset
      })

      # Dynamically generate UI inputs for Single_Plot
      output$dynamic_inputs_1 <- renderUI({
        generate_inputs(Single_Plot)  # Generate dynamic inputs for Single_Plot
      })

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
        return(user_inputs)
      })

      # Save and display the Single_Plot plot as an image
      output$renderedPlot1 <- renderImage({
        req(input$submit1)

        # Generate a temporary file to store the saved plot
        temp_file <- tempfile(fileext = paste0(".", input$File_Type))
        plot_file_1(temp_file)  # Store the temp file globally

        # Collect user inputs and dynamically pass them to Single_Plot
        user_args <- reactive_inputs_1()

        # Create the plot and save it to the temporary file
        Plot_Outcome <- do.call(Single_Plot, user_args)
        ggplot2::ggsave(temp_file, plot = Plot_Outcome, width = input$Width,
                        height = input$Height, units = "in", dpi = input$Quality)

        # Return the image to be displayed in the UI
        list(src = temp_file, contentType = paste0("image/", input$File_Type),
             width = "100%", height = "auto")
      }, deleteFile = FALSE)  # Do not delete the file after rendering

      # Download handler for Single_Plot
      output$downloadPlot1 <- downloadHandler(
        filename = function() { paste0("Single_Plot_", Sys.Date(), ".", input$File_Type) },
        content = function(file) {
          # If the plot file already exists, use it
          if (!is.null(plot_file_1())) {
            file.copy(plot_file_1(), file)
          } else {
            # Otherwise, regenerate the plot
            user_args <- reactive_inputs_1()
            Plot_Outcome <- do.call(Single_Plot, user_args)
            ggplot2::ggsave(file, plot = Plot_Outcome, width = input$Width,
                            height = input$Height, units = "in", dpi = input$Quality)
          }
        }
      )
    })


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
      req(input$topData, input$bottomData)  # Wait for both datasets to be uploaded
      generate_inputs(Miami_Plot)  # Generate dynamic inputs for Miami_Plot only after data is uploaded
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
    observe({
      lapply(1:input$numDatasets, function(i) {
        local({
          datasetInput <- paste0("forest_file", i)  # Updated input IDs for Forest_Plot
          previewOutput <- paste0("forestDataPreview", i)  # Updated preview IDs for Forest_Plot

          output[[previewOutput]] <- renderTable({
            req(input[[datasetInput]])  # Only show preview for Forest_Plot uploads
            vroom(input[[datasetInput]]$datapath, n_max = 10)  # Preview first 10 rows of the uploaded data
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
        generate_inputs(Annotate_Data)  # Generate dynamic inputs for Annotate_Data
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
        generate_inputs(Model_Munge)  # Generate dynamic inputs for Model_Munge
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
