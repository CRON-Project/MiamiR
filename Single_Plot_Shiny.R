  library(shiny)

  library(ggplot2)
  library(vroom)
  library(dplyr)
  library(plotly)
  library(htmlwidgets)

  # Increase file upload limit to 10 GB
  options(shiny.maxRequestSize = 10 * 1024^3)
  options(shiny.http.timeout = 300)  # Timeout increased to 5 minutes
  options(shiny.maxRequestSize = 10 * 1024^3)  # 30 MB


  # Helper function to generate dynamic UI inputs, specifically handling list inputs
  generate_inputs <- function(func) {
    arg_list <- formals(func)

    inputs <- lapply(names(arg_list), function(arg_name) {
      # Skip known data inputs and arguments with no defaults
      if (arg_name %in% c("Data", "Top_Data", "Bottom_Data")) return(NULL)

      arg_value <- arg_list[[arg_name]]

      # Skip missing (no default) arguments
      if (missing(arg_value) || is.symbol(arg_value)) return(NULL)

      # Handle vector or list inputs
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

    # Add the download config inputs
    file_type_input <- textInput("File_Type", label = "File Type", value = "jpg")
    width_input <- numericInput("Width", label = "Plot Width (in)", value = 30)
    height_input <- numericInput("Height", label = "Plot Height (in)", value = 15)
    dpi_input <- numericInput("Quality", label = "Plot DPI", value = 300)

    return(c(
      list(file_type_input, width_input, height_input, dpi_input),
      Filter(Negate(is.null), inputs)
    ))
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


      tabPanel("Automated Batch Mode",
               sidebarLayout(
                 sidebarPanel(
                   fileInput("batchData", "Upload Data (applies to all selected functions)", accept = ".csv"),

                   checkboxGroupInput("selectedFunctions", "Select Functions to Run:",
                                      choices = c("Single_Plot", "Regional_Plot", "Miami_Plot", "Forest_Plot", "Annotate_Data", "Model_Munge"),
                                      selected = NULL),

                   actionButton("runBatch", "Run Selected Functions", class = "btn-success")
                 ),
                 mainPanel(
                   h4("Batch Execution Results:"),
                   uiOutput("batchResults")
                 )
               )
      ),



      tabPanel("Single_Plot()",  # First Tab for Single_Plot
               sidebarLayout(
                 sidebarPanel(
                   fileInput("file1", "Upload Any File Type (up to 10GB)",
                             accept = c("text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".csv")),

                   fluidRow(
                     column(6, uiOutput("generateBtn1")),
                     column(6, uiOutput("downloadBtn1"))
                   ),
                   tags$br(), tags$br(),

                   uiOutput("dynamic_inputs_1"),  # Dynamically generated inputs

                   #   textInput("File_Type", label = "File Type", value = "jpg")  # Add this line
                 ),
                 mainPanel(
                   h4("Data Preview:"),
                   tableOutput("dataPreview1"),

                   h4("Rendered Plot:"),
                   uiOutput("fullscreenButton1"),  # <- This stays above the plot
                   uiOutput("plotStatusUI"),
                   uiOutput("interactivePlot1_ui")
                 )

               )
      ),


      tabPanel("Regional_Plot()",
               sidebarLayout(
                 sidebarPanel(
                   fileInput("file6", "Upload Data (up to 10GB)",
                             accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),

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
    ),
    tags$script(HTML("
    Shiny.addCustomMessageHandler('plot_progress', function(message) {
      const el = document.getElementById('plotProgressStatus');
      if (el) {
        el.textContent = 'Progress: ' + message.pct + '% (' + message.msg + ') based on script lines';
      }
    });
  "))
  )

  server <- function(input, output, session) {


    # Helper function to build Single_Plot arguments for batch
    build_single_plot_args_batch <- function(batch_data) {
      args <- formals(Single_Plot)
      args <- lapply(args, function(x) if (is.call(x)) eval(x) else x)

      args$Data <- batch_data

      # Custom fixes matching tab logic
      if (is.null(args$Title) || args$Title == "") {
        if (is.null(args$Title) || args$Title == "") {
          if (!is.null(input$batchData$name)) {
            args$Title <- tools::file_path_sans_ext(input$batchData$name)
          } else {
            args$Title <- "Auto Single Plot"
          }
        }

      }
      if (is.null(args$Interactive) || args$Interactive == "") {
        args$Interactive <- FALSE
      }

      return(args)
    }

    # Helper function to build Regional_Plot arguments for batch
    build_regional_plot_args_batch <- function(batch_data) {
      args_primary <- formals(Regional_Plot)
      args_inherited <- formals(Single_Plot)

      args_combined <- modifyList(lapply(args_inherited, function(x) if (is.call(x)) eval(x) else x),
                                  lapply(args_primary, function(x) if (is.call(x)) eval(x) else x))

      args_combined$Data <- batch_data

      # 🔥 Regional_Plot-specific fixes
      if (is.null(args_combined$Chromosome) || args_combined$Chromosome == "") {
        args_combined$Chromosome <- 1
      }

      # 🔥 Single_Plot-specific inherited fixes
      if (is.null(args_combined$Title) || args_combined$Title == "") {
        if (!is.null(input$batchData$name)) {
          args_combined$Title <- tools::file_path_sans_ext(input$batchData$name)
        } else {
          args_combined$Title <- "Auto Regional Plot"
        }
      }
      if (is.null(args_combined$Interactive) || args_combined$Interactive == "") {
        args_combined$Interactive <- FALSE
      }

      # 🧹 Important: Remove any arguments not actually expected by Regional_Plot
      expected_args <- names(formals(Regional_Plot))
      all_args <- args_combined

      all_args <- c(all_args, list(... = args_combined[setdiff(names(args_combined), expected_args)]))

      # Keep only expected top-level Regional_Plot args and ... (extra passed separately)
      args_final <- all_args[expected_args]

      args_final$... <- all_args$...

      return(args_final)
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

      menu.innerHTML = `
        <div style='font-weight:bold; margin-bottom:5px;'>View SNP in:</div>
        <div style='margin-bottom:4px;'><button id='dbsnpBtn' style='width:100%;'>dbSNP</button></div>
        <div style='margin-bottom:4px;'><button id='otBtn' style='width:100%;'>Open Targets</button></div>
        <div><button id='varsomeBtn' style='width:100%;'>VarSome</button></div>
      `;

      document.body.appendChild(menu);

      // Track hover into menu
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


setTimeout(() => {
  const isFullscreen = window.innerWidth === screen.width && window.innerHeight === screen.height;

  menu.style.setProperty('font-size', isFullscreen ? '20px' : '12px', 'important');
  menu.style.setProperty('padding', isFullscreen ? '16px 18px' : '6px 8px', 'important');
  menu.style.setProperty('max-width', isFullscreen ? '300px' : '150px', 'important');

  const buttons = menu.querySelectorAll('button');
  buttons.forEach(btn => {
    btn.style.setProperty('font-size', isFullscreen ? '18px' : '12px', 'important');
    btn.style.setProperty('padding', isFullscreen ? '10px 12px' : '4px 6px', 'important');
  });
}, 50); // Slight delay to ensure fullscreen is applied


      document.getElementById('dbsnpBtn').onclick = function() {
        if (snp_id) window.open('https://www.ncbi.nlm.nih.gov/snp/' + snp_id[0], '_blank');
        menu.style.display = 'none';
      };
      document.getElementById('otBtn').onclick = function() {
        if (chrom && pos && ref && alt) {
          window.open(`https://genetics.opentargets.org/variant/${chrom}_${pos}_${ref}_${alt}/associations`, '_blank');
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

    // Also hide when clicking elsewhere
    document.addEventListener('click', function(e) {
      const menu = document.getElementById('snpMenu');
      if (!e.target.closest('#snpMenu') && !e.target.closest('.plotly')) {
        menu.style.display = 'none';
      }
    });
  }
")





    }


    run_with_counter <- function(func, args = list(), session = NULL) {
      exprs <- as.list(body(func))[-1]
      total <- length(exprs)
      env <- new.env(parent = environment(func))

      # Merge defaults from function formals
      defaults <- formals(func)
      for (name in names(defaults)) {
        if (name %in% names(args)) {
          assign(name, args[[name]], envir = env)
        } else {
          assign(name, eval(defaults[[name]], envir = env), envir = env)
        }
      }

      result <- NULL
      for (i in seq_along(exprs)) {
        pct <- round(100 * i / total)
        message(sprintf("Progress: %d%% (Line %d of %d)", pct, i, total))

        if (!is.null(session)) {
          session$sendCustomMessage("plot_progress", list(
            pct = pct,
            msg = paste("Line", i, "of", total)
          ))
        }

        result <- eval(exprs[[i]], envir = env)
      }

      return(result)
    }



    observeEvent(input$runBatch, {
      req(input$batchData)
      req(input$selectedFunctions)

      batch_data <- vroom(input$batchData$datapath)

      results <- list()

      ### 🏗 1. Run all selected functions correctly with preprocessing ###
      if ("Single_Plot" %in% input$selectedFunctions) {
        args <- build_single_plot_args_batch(batch_data)
        results$Single_Plot <- run_with_counter(Single_Plot, args = args, session = session)
      }

      if ("Regional_Plot" %in% input$selectedFunctions) {
        args <- build_regional_plot_args_batch(batch_data)
        results$Regional_Plot <- do.call(Regional_Plot, args)
      }

      if ("Miami_Plot" %in% input$selectedFunctions) {
        args <- build_miami_plot_args_batch(batch_data)
        results$Miami_Plot <- do.call(Miami_Plot, args)
      }

      if ("Annotate_Data" %in% input$selectedFunctions) {
        args <- build_annotate_data_args_batch(batch_data)
        results$Annotate_Data <- do.call(Annotate_Data, args)
      }

      ### 🛠 2. Now Create Outputs Dynamically ###

      # Single_Plot output
      if (!is.null(results$Single_Plot)) {
        if (inherits(results$Single_Plot, "plotly") || isTRUE(attr(results$Single_Plot, "Interactive"))) {
          output$autoSingle <- renderPlotly({
            results$Single_Plot
          })
        } else {
          output$autoSingle <- renderImage({
            temp_file <- tempfile(fileext = ".jpg")
            ggsave(temp_file, plot = results$Single_Plot,
                   width = 30, height = 15, units = "in", dpi = 300)

            list(
              src = temp_file,
              contentType = "image/jpeg",
              width = "100%",
              height = "auto",
              alt = "Static Manhattan plot"
            )
          }, deleteFile = FALSE)
        }
      }

      # Regional_Plot output
      if (!is.null(results$Regional_Plot)) {
        first_plot <- results$Regional_Plot[[1]]

        if (!is.null(attr(first_plot, "interactive_panel"))) {
          output$autoRegional <- renderPlotly({
            attr(first_plot, "interactive_panel")
          })
        } else {
          output$autoRegional <- renderImage({
            temp_file <- tempfile(fileext = ".jpg")
            ggsave(temp_file, plot = first_plot,
                   width = 30, height = 15, units = "in", dpi = 300)

            list(
              src = temp_file,
              contentType = "image/jpeg",
              width = "100%",
              height = "auto",
              alt = "Static Regional plot"
            )
          }, deleteFile = FALSE)
        }
      }

      # Miami_Plot output (always static)
      if (!is.null(results$Miami_Plot)) {
        output$autoMiami <- renderImage({
          temp_file <- tempfile(fileext = ".jpg")
          ggsave(temp_file, plot = results$Miami_Plot,
                 width = 30, height = 15, units = "in", dpi = 300)

          list(
            src = temp_file,
            contentType = "image/jpeg",
            width = "100%",
            height = "auto",
            alt = "Static Miami plot"
          )
        }, deleteFile = FALSE)
      }

      # Annotate_Data output (always table)
      if (!is.null(results$Annotate_Data)) {
        output$autoAnno <- renderTable({
          results$Annotate_Data %>% filter(grepl("rs", Lab, ignore.case = TRUE))
        })
      }

      ### 🏗 3. Finally create the correct batch UI ###
      output$batchResults <- renderUI({
        tagList(
          if (!is.null(results$Single_Plot)) {
            if (inherits(results$Single_Plot, "plotly") || isTRUE(attr(results$Single_Plot, "Interactive"))) {
              plotlyOutput("autoSingle")
            } else {
              imageOutput("autoSingle")
            }
          },
          if (!is.null(results$Regional_Plot)) {
            first_plot <- results$Regional_Plot[[1]]
            if (!is.null(attr(first_plot, "interactive_panel"))) {
              plotlyOutput("autoRegional")
            } else {
              imageOutput("autoRegional")
            }
          },
          if (!is.null(results$Miami_Plot)) imageOutput("autoMiami"),
          if (!is.null(results$Annotate_Data)) tableOutput("autoAnno")
        )
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
        generate_inputs(Single_Plot)  # Generate dynamic inputs for Single_Plot
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
          tags$div(
            style = "padding: 50px; text-align: center; font-size: 18px; color: #333;",
            "⏳ Rendering plot... Please wait.",
            tags$br(),
            tags$div(
              id = "plotProgressStatus",
              style = "font-weight: bold; font-size: 20px; color: #007bff;"
            )
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
          "Please generate your plot!",
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





      output$modalInteractivePlot6 <- renderPlotly({
        req(Regional_Plots_Results())
        plot_names <- names(Regional_Plots_Results())

        idx <- current_plot_index()
        full_plot <- Regional_Plots_Results()[[plot_names[idx]]]

        # Extract components
        p_main <- attr(full_plot, "main_ggplot")
        p_genes <- attr(full_plot, "gene_track_plot")
        height_in <- attr(full_plot, "dynamic_height")
        if (is.null(height_in)) height_in <- 25

        fig1 <- ggplotly(p_main, tooltip = "text") %>%
          layout(margin = list(t = 50, b = 50, l = 50, r = 50)) %>%
          config(displayModeBar = TRUE, responsive = TRUE)

        fig2 <- ggplotly(p_genes + theme(legend.position = "none"), tooltip = "text") %>%
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


      output$downloadPlot6 <- downloadHandler(
        filename = function() { paste0("Regional_Plots_", Sys.Date(), ".zip") },
        content = function(file) {
          saved_files <- Regional_Plots_Results()
          zip::zipr(zipfile = file, files = saved_files)
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
