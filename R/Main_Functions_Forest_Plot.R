
#' Annotated, Customised Forest Plots
#'
#' @param Names These are a list of manual data set names to be applied to the legend text/key on the main plot - this also specifies the order down the page; defaults to NULL
#' @param Data_Set_Colours These are a list of manual data set colours allocated within the main plot - this also specifies the order down the page; defaults to NULL
#' @param Model_Reference Manually specify where the statistics to be used derive from (either from GWAS or Models); defaults to FALSE - not a model
#' @param Shapes These are a list of manual data set shapes allocated within the main plot - this also specifies the order down the page; defaults to NULL
#' @param Null_Line_Colour Colour of null line; defaults to "red"
#' @param Null_Line_Type Type of null line; defaults to "dashed"
#' @param X_Axis_Title X axis title; defaults to NULL
#' @param X_Axis_Title_Size X axis title size; defaults to 13
#' @param X_Axis_Label Display X axis title; defaults to TRUE
#' @param X_Axis_Separation Continuous value of the specific break separation to show on the X axis; defaults to NULL - leading to auto-calculation
#' @param Strip_Colour Colour of the strips highlighting the groups of covariates or SNPs under a certain category; defaults to NULL
#' @param Strips Display strips highlighting groups of interest; defaults to TRUE
#' @param Pre_Calculated_CIs Indicate that custom CIs are precalculated and available in datasets; defaults to FALSE
#' @param Legend_Title_Size Legend title size; defaults to 13
#' @param Legend_Text_Size Legend text size; defaults to 11
#' @param Legend_Title Title to display next to the legend; defaults to "Data"
#' @param Left_Title Title to display for the information in the grid text module on the left hand side of the plot; defaults to NULL
#' @param Match_Allele_Direction Match the allele directions to avoid effect direction flips; defaults to TRUE
#' @param Match_Allele_Study Data set used as reference point for effect allele direction; defaults to NULL - automatically assigned as first set
#' @param Selected_SNPs A list of SNPs to be shown in the forest plot; defaults to NULL - automatically selected for from the summary statistics data set
#' @param Selected_Covariates A list of covariates to be shown in the forest plot; defaults to NULL - automatically selected for from the munged model data set
#' @param X_Axis_Text_Resolution Number of decimal places to display on X axis text; defaults to NULL
#' @param Legend_On Display the main plot legend; defaults to TRUE
#' @param X_Axis_Text_Size Size of the X axis text labels; defaults to 11
#' @param Null_Buffer Units of space to avoid for X axis text labels around the NULL point; defaults to 0.00
#' @param Data These are the summary statistics (either file path(s), environmental object(s) or a variable(s) containing such contents) to be plotted; defaults to NULL
#' @param Styles These are a list of styles to be applied to the grid text on the right hand side of the main plot - this also specifies the order down the page; defaults to NULL and accepts "underline" or "bold" as deviations
#' @param Missings Simulate shape positions for data points which may be missing from specific data sets plotted; defaults to TRUE
#' @param Null_Line_Thickness Manually specify the thickness of the corresponding Null Line; defaults to 0.5
#' @param Test_Statistic_Choice Manually specify which type of test statistic output is desired, OR or BETA; defaults to NULL - regardless of input type either output may be specified
#' @param Lead_Studies Manually specify which study's lead SNPs should be plotted when Selected_SNPs is NULL; defaults to NULL - leads to index SNPs from any set to be plotted
#' @param Separation_Distance The minimum genomic distance between significant regions of interest/lead SNPs to be defined as distinct peaks, defined by P_threshold value when Selected_SNPs is NULL; defaults to 1e6
#' @param P_threshold Maximum value of P (minimum LOG10P) allowed for a SNP to be classified as initiating detection as the lead SNP of a peak in the function, when Selected_SNPs is NULL; defaults to 5e-8
#' @param Chromosomes Manually specify which chromosomes from the summary statistics to run the function on, when Selected_SNPs is NULL; defaults to NULL - leads to all chromosomes being processed
#' @param Double_Label Display a second level of information under individual grid text layers on the left hand side of the main plot; defaults to FALSE
#' @param Shape_Size Manually specify the size of the Shapes added on main plot; defaults to 4
#' @param CI_Line_Height Manually specify the vertical length of the CI lines added on main plot; defaults to 0.35
#' @param CI_Line_Size Manually specify the thickness of the CI lines added on main plot; defaults to 0.4
#' @param Shape_Thickness Manually specify the thickness of the Shapes added on main plot; defaults to 1
#' @param Grid_Thickness Manually specify the thickness of the grid lines, when Grid_Lines_On is TRUE; defaults to 0.25
#' @param Grid_Lines_On Show grid lines with specified settings around left and right hand side annotations; defaults to TRUE
#' @param SE_Column Manually specify SE of Test Statistic column; defaults to NULL, leading to auto-detection
#' @param OR_Column Manually specify OR Test Statistic column; defaults to NULL, leading to auto-detection
#' @param Chromosome_Column Manually specify chromosome column; defaults to NULL, leading to auto-detection
#' @param BETA_Column Manually specify BETA Test Statistic column; defaults to NULL, leading to auto-detection
#' @param Position_Column Manually specify chromosome genomic position column; defaults to NULL, leading to auto-detection
#' @param SNP_ID_Column Manually specify SNP ID column; defaults to NULL, leading to auto-detection
#' @param PValue_Column Manually specify P Value column; defaults to NULL, leading to auto-detection
#' @param Reference_Allele_Column Manually specify reference allele column; defaults to NULL, leading to auto-detection
#' @param Effect_Allele_Column Manually specify effect allele column; defaults to NULL, leading to auto-detection
#' @param Upper_CI_Column Manually specify upper CI column, when Pre_Calculated_CIs is TRUE; defaults to NULL, leading to auto-detection
#' @param Lower_CI_Column Manually specify lower CI column, when Pre_Calculated_CIs is TRUE; defaults to NULL, leading to auto-detection
#' @param Left_Size Manually specify the size of the grid text on the left hand side of the main plot; defaults to 9
#' @param Right_Size Manually specify the size of the grid text on the right hand side of the main plot; defaults to 9
#' @param Right_Buffer Manually specify the size of space around the grid text on the right hand side of the main plot; defaults to 1
#' @param Left_Buffer Manually specify the size of space around the grid text on the left hand side of the main plot; defaults to 2
#' @param Axis_Buffer Manually specify the size of space around the axis end positions of the main plot; defaults to 0.02
#' @param Left_Fill Manually specify the colour of the grid text background on the left hand side of the main plot; defaults to "white"
#' @param Right_Fill Manually specify the colour of the grid text background on the right hand side of the main plot; defaults to "white"
#' @param Grid_Colour_Left Manually specify the colour of the grid lines on the left hand side of the main plot; defaults to "black"
#' @param Grid_Colour_Right Manually specify the colour of the grid lines on the right hand side of the main plot; defaults to "black"
#' @param CI_Label Allow the display of a column containing CI information in the grid text module on the right hand side of the plot; defaults to TRUE
#' @param BETA_Label Allow the display of a column containing BETA/OR information, depending on Test_Statistic_Choice in the grid text module on the right hand side of the plot; defaults to TRUE
#' @param SE_Label Allow the display of a column containing SE information in the grid text module on the right hand side of the plot; defaults to TRUE
#' @param BETA_SE_Label Allow the display of a column containing BETA & SE information in the grid text module on the right hand side of the plot; defaults to TRUE
#' @param BETA_CI_Label Allow the display of a column containing BETA & CI information in the grid text module on the right hand side of the plot; defaults to TRUE
#' @param P_Label Allow the display of a column containing P-Value information in the grid text module on the right hand side of the plot; defaults to TRUE
#' @param Special_P Format the column containing P-Value information in the grid text module on the right hand side of the plot when P_Label is TRUE and P_Label is in order() to display as * or NS; defaults to FALSE
#' @param CI_Label_Resolution Specify the number of decimal places of the column containing CI information in the grid text module on the right hand side of the plot; defaults to 3
#' @param BETA_Label_Resolution Specify the number of decimal places of the column containing BETA information in the grid text module on the right hand side of the plot; defaults to 3
#' @param SE_Label_Resolution Specify the number of decimal places of the column containing SE information in the grid text module on the right hand side of the plot; defaults to 3
#' @param P_Label_Resolution Specify the number of decimal places of the column containing P-Value information in the grid text module on the right hand side of the plot; defaults to 2
#' @param CI_Label_Title Specify the title of the column containing CI information in the grid text module on the right hand side of the plot; defaults to "CI"
#' @param BETA_Label_Title Specify the title of the column containing BETA information in the grid text module on the right hand side of the plot; defaults to NULL
#' @param SE_Label_Title Specify the title of the column containing SE information in the grid text module on the right hand side of the plot; defaults to "SE"
#' @param BETA_SE_Label_Title Specify the title of the column containing SE information in the grid text module on the right hand side of the plot; defaults to NULL
#' @param BETA_CI_Label_Title Specify the title of the column containing BETA & CI information in the grid text module on the right hand side of the plot; defaults to NULL
#' @param P_Label_Title Specify the title of the column containing P information in the grid text module on the right hand side of the plot; defaults to "P"
#' @param order Specify the order and selection of information containing columns (Labels) in the grid text module on the right hand side of the plot; defaults to c("BETA_SE_Label", "P_Label")
#' @param Verbose Prevent display of progress bar as function is running and instead show key milestone outputs/messages (mainly for debugging purposes); defaults to FALSE
#' @param Lab_Col Allocate the second level of information in the grid text module on the left hand side of the plot, when Double_Label is TRUE; defaults to "Real_ID"
#'
#' @return Image of Forest Plot is allocated to specified object and the resulting ggplot object can then be saved to an image
#' @export
#'
#' @examples Forest_Plot <- Forest_Plot(Data = "Intelligence_Sum_Stats", Chromosomes = 1)
#'


      Forest_Plot <- function(
        Data = NULL,
        Names = NULL,
        Styles = NULL,
        Shapes = NULL,
        Data_Set_Colours = NULL,
        Missings = TRUE,
        Null_Line_Colour = "red",
        Null_Line_Type = "dashed",
        Null_Line_Thickness = 0.5,
        Strip_Colour = "goldenrod",
        Strips = TRUE,
        Legend_Title = "Data",
        Legend_Title_Size = 13.0,
        Legend_Text_Size = 11.0,
        Legend_On = TRUE,
        X_Axis_Label = TRUE,
        X_Axis_Title = NULL,
        X_Axis_Text_Size = 11,
        X_Axis_Title_Size = 13,
        Test_Statistic_Choice = "BETA",
        X_Axis_Separation = NULL,
        X_Axis_Text_Resolution = NULL,
        Lead_Studies = NULL,
        Separation_Distance = 1e6,
        P_threshold = 5e-8,
        Chromosomes = NULL,
        Double_Label = FALSE,
        Shape_Size = 4,
        CI_Line_Height = 0.35,
        CI_Line_Size = 0.4,
        Shape_Thickness = 1,
        Match_Allele_Study = NULL,
        Match_Allele_Direction = TRUE,
        Grid_Thickness = 0.25,
        Grid_Lines_On = TRUE,
        SE_Column = NULL,
        OR_Column = NULL,
        Chromosome_Column = NULL,
        BETA_Column = NULL,
        Position_Column = NULL,
        SNP_ID_Column = NULL,
        PValue_Column = NULL,
        Reference_Allele_Column = NULL,
        Effect_Allele_Column = NULL,
        Upper_CI_Column = NULL,
        Lower_CI_Column = NULL,
        Pre_Calculated_CIs = FALSE,
        Left_Size = 9,
        Right_Size = 9,
        Right_Buffer = 1,
        Left_Buffer = 2,
        Axis_Buffer = 0.02,
        Left_Fill = "white",
        Right_Fill = "white",
        Grid_Colour_Left = "black",
        Grid_Colour_Right = "black",
        CI_Label = TRUE,
        BETA_Label = TRUE,
        SE_Label = TRUE,
        BETA_SE_Label = TRUE,
        BETA_CI_Label = TRUE,
        P_Label = TRUE,
        Special_P = FALSE,
        CI_Label_Resolution = 3,
        BETA_Label_Resolution = 3,
        SE_Label_Resolution = 3,
        P_Label_Resolution = 2,
        Null_Buffer = 0.00,
        Left_Title = NULL,
        CI_Label_Title = "CI",
        BETA_Label_Title = NULL,
        SE_Label_Title = "SE",
        BETA_SE_Label_Title = NULL,
        BETA_CI_Label_Title = NULL,
        P_Label_Title = "P",
        order = c("BETA_SE_Label", "P_Label"),
        Selected_Covariates = NULL,
        Selected_SNPs = NULL,
        Model_Reference = FALSE,
        Verbose = FALSE,
        Lab_Col = "Real_ID"

      )

{



  if (is.null(Data) || length(Data) == 0) {

    stop("Error: 'Data' is NULL or empty. Please provide valid data sets.")

  }

  if(!is.null(Selected_SNPs))
  {

    message("No chromosome filtering required as specific SNPs supplied")

    Chromosomes <- NULL

  }

  Names_Run <- if (exists(".NF_Orig_Names", envir = .GlobalEnv)) {

    message("Retrieving file names")

    get(".NF_Orig_Names", envir = .GlobalEnv)

  } else {

    message("Using custom file names")

    Names  # use the user-supplied Names if present

  }



  call_expr <- Names_Run

  Orig_Names <- call_expr

  if(Grid_Lines_On == FALSE)
  {

    message("Setting plot element arguments for minimalist plot")

    Grid_Colour_Left <- "transparent"
    Grid_Colour_Right <- "transparent"
    Grid_Thickness <- 0
    Top_Grid_Thickness <- 0

  }else{

    Top_Grid_Thickness <- Grid_Thickness

  }

  message("Checking names for path formatting")
  Orig_Names <- gsub('^"|"$', '', Orig_Names)

  if(Model_Reference == T && is.null(Left_Title))
  {
    Left_Title <- "Covariate"
  }

  if(Model_Reference == F && is.null(Left_Title))
  {
    Left_Title <- "SNP"
  }

  if(is.null(Lab_Col))
  {
    Lab_Col <- "ID"
  }

  if(is.null(Data_Set_Colours))
  {
    Data_Set_Colours <- viridis::viridis(length(Orig_Names))
  }

  if(!is.null(Names))
  {
    message("Using Names Provided")

  }else{

    Names <- Orig_Names

  }

  if (any(grepl("[/\\\\]", Names))) {
    Names <- sub(".*[\\\\/]", "", Names)
  }


  #Goes below above as need to match based on STUDY assigned from Name
  if (is.null(Match_Allele_Study) & Match_Allele_Direction == TRUE) {

    Match_Allele_Study <- Names[1]
    message(paste0("Matching allele directions automatically to ", Match_Allele_Study) )

  }

  if (is.null(Shapes)) {

    Shapes <- rep("square", length(Names))

  }

  if (is.null(Styles)) {

    Styles <- rep("normal", length(Names))

  }

  message("Assigning statistic information titles")

  if(is.null(BETA_SE_Label_Title))
  {

    if(Test_Statistic_Choice == "BETA")
    {

      BETA_SE_Label_Title <- "BETA (SE)"

    }else{

      BETA_SE_Label_Title <- "OR (SE)"

    }

  }

  if(is.null(BETA_CI_Label_Title))
  {

    if(Test_Statistic_Choice == "BETA")
    {

      BETA_CI_Label_Title <- "BETA (CI)"

    }else{

      BETA_CI_Label_Title <- "OR (CI)"

    }

  }

  if(is.null(BETA_Label_Title))
  {

    if(Test_Statistic_Choice == "BETA")
    {

      BETA_Label_Title <- "BETA"

    }else{

      BETA_Label_Title <- "OR"

    }

  }

  if (is.null((X_Axis_Title))) {

    message("Assigning default X Axis Title")

    X_Axis_Title <- Test_Statistic_Choice

    X_Axis_Title <- paste0("\n", X_Axis_Title, "\n")
  }

  Combined_Processed_Data <- data.frame()
  Combined_Processed_SNPs <- data.frame()

  message("Loading in data")

  for (i in seq_along(Orig_Names)) {

    Orig_Name <- Orig_Names[i]
    Spec_Name <- Names[i]

    message("Index:", i, " Name:", Names[i], "\n")

    Data <- if (exists(Orig_Names[i], envir = .GlobalEnv)) {

      message("Reading data frame from environment")

      get(Orig_Names[i], envir = .GlobalEnv)

    } else if (file.exists(Orig_Names[i])) {

      message("Reading from file path")

      vroom::vroom(Orig_Names[i], show_col_types = FALSE)
    } else {

      stop("The name or file path '", Orig_Names[i], "' does not exist.")

    }

    if(Model_Reference == TRUE)
    {

      message("Munging model object")

      Data <- Model_Munge(Model_Object = Orig_Names[i], Verbose = TRUE)

    }

    corresponding_name <- Names[i]
    corresponding_shape <- Shapes[i]
    corresponding_style <- Styles[i]

    message("Deducing key column names")

    Chromosome_Column <- detect_chromosome_column(Data, Chromosome_Column)

    if(!is.null(Chromosome_Column))

    {

      Data <- Data %>%
        dplyr::mutate(!!Chromosome_Column := ifelse(.data[[Chromosome_Column]] == "23", "X", .data[[Chromosome_Column]]))


    }

    PValue_Column      <- detect_pvalue_column(Data, PValue_Column)
    Position_Column    <- detect_position_column(Data, Position_Column)
    SNP_ID_Column      <- detect_snp_column(Data, SNP_ID_Column)
    Ref_Allele_Column  <- detect_reference_allele_column(Data, Reference_Allele_Column)
    Alt_Allele_Column  <- detect_effect_allele_column(Data, Effect_Allele_Column)
    BETA_Column      <- detect_beta_column(Data, BETA_Column)
    OR_Column      <- detect_or_column(Data, OR_Column)
    SE_Column      <- detect_se_column(Data, SE_Column)
    Upper_CI_Column      <- detect_upper_ci_column (Data, Upper_CI_Column)
    Lower_CI_Column      <- detect_lower_ci_column (Data, Lower_CI_Column)

    message("Assigning key columns to standardised nomenclature")

    if(!is.null(Chromosome_Column))

    {

      Data$CHROM <- Data[[Chromosome_Column]]

      if(!is.null(Chromosomes))

      {

        message("Filtering for key chromosomes")

        Data <- Data[Data$CHROM %in% Chromosomes, ]

      }

      Data$GENPOS <- Data[[Position_Column]]
      Data$ID <- Data[[SNP_ID_Column]]

      Data$ALLELE0 <- toupper(Data[[Ref_Allele_Column]])
      Data$ALLELE1 <- toupper(Data[[Alt_Allele_Column]])

    }

    Data <- ensure_P(Data)
    Data$P <- Data[[PValue_Column]]

    message("Deducing type of test statistic supplied and performing required OR/BETA calculations")

    Test_Statistic <- NULL

    if (is.null(Test_Statistic)) {

      if (!is.null(OR_Column) && is.null(BETA_Column)) {

        Test_Statistic <- "OR"

      } else if (!is.null(BETA_Column) && is.null(OR_Column)) {

        Test_Statistic <- "BETA"

      } else {

        Test_Statistic <- "BETA" # Default if both are NULL
      }

    }

    if (Test_Statistic == "BETA") {

      # Assign BETA and rename SE
      Data$BETA    <- Data[[BETA_Column]]
      Data$BETA_SE <- Data[[SE_Column]]

      # Calculate OR and OR_SE
      Data$OR    <- exp(Data$BETA)
      Data$OR_SE <- Data$OR * Data$BETA_SE

    }

    if (Test_Statistic == "OR") {

      # Assign OR and rename SE
      Data$OR    <- Data[[OR_Column]]
      Data$OR_SE <- Data[[SE_Column]]

      # Calculate BETA and BETA_SE
      Data$BETA    <- log(Data$OR)
      Data$BETA_SE <- Data$OR_SE / Data$OR

    }

    if(Pre_Calculated_CIs == T)
    {

      Data$LL <- Data[[Upper_CI_Column]]
      Data$UL <- Data[[Lower_CI_Column]]

    }

    if(Model_Reference == T)
    {

      Data$ID <- Data$Covariate

    }

    if(Model_Reference == T)
    {

      if(!is.null(Selected_Covariates)){
        message("Filtering for key Covariates (manual)")
        Selected_SNPs <- Selected_Covariates
        message(Selected_SNPs)
        Data <- Data[Data$group %in% Selected_SNPs,]

      }else{

        message("No covariates specified - retaining all (automatic)")
        message("Automatically generated order:")
        Selected_Covariates <- Data$group[!duplicated(Data$group)]
        message(Selected_Covariates)
        Selected_SNPs <- Selected_Covariates

      }

    }

    Data$STUDY <- corresponding_name
    Data$Shape <- corresponding_shape
    Data$Style <- corresponding_style
    Data$Real_ID <- Data$ID

    if(is.null(Selected_SNPs) & Model_Reference == F)
    {

      find_peaks <- function(data, Separation_Distance, P_threshold) {

        # Group by chromosome (or single group if only one CHROM)
        chrom_list <- if (length(unique(data$CHROM)) == 1) list(data) else base::split(data, data$CHROM)

        results <- base::lapply(chrom_list, function(chrom_data) {

          # Filter to only SNPs below genome-wide threshold - MAXIMUM P Value for significant hit/minimum LOG10P
          chrom_data <- chrom_data[chrom_data$P < P_threshold, , drop = FALSE]

          # If none pass threshold, return nothing
          if (nrow(chrom_data) == 0) return(chrom_data)

          # Sort by smallest P first
          chrom_data <- chrom_data[base::order(chrom_data$P), ]

          selected <- chrom_data[0, , drop = FALSE]

          for (i in base::seq_len(nrow(chrom_data))) {

            snp <- chrom_data[i, , drop = FALSE]

            if (nrow(selected) == 0) {

              selected <- snp

            } else {

              # Only keep if this SNP is at least Separation_Distance away from all previously selected - basically loop through peaks.
              too_close <- base::any(base::abs(snp$GENPOS - selected$GENPOS) < Separation_Distance)

              if (!too_close) {

                selected <- dplyr::bind_rows(selected, snp)

              }

            }

          }

          selected
        })

        dplyr::bind_rows(results)
      }

      message("Calling find_peaks algorithm based on user supplied significance and distinct peak distance threholds")

      peaks <- find_peaks(Data, Separation_Distance = Separation_Distance, P_threshold = P_threshold)

      message("Lead SNPs of defined peaks ascertained")

      found_snps <- peaks %>% dplyr::select(ID, STUDY, CHROM)

      message("Concatenating SNPs of interest")

      Combined_Processed_SNPs <- rbind(Combined_Processed_SNPs, found_snps)

    }

    message("Cleaning processed data")

    if(Model_Reference == F & Double_Label == T)

    {

      Data <- Data %>% dplyr::select(ID, ALLELE0, ALLELE1, CHROM, GENPOS, BETA, OR, BETA_SE, OR_SE,  P, STUDY, Shape, Style, !!Lab_Col,Real_ID)

    }

    if(Model_Reference == F & Double_Label == F)

    {

      Data <- Data %>% dplyr::select(ID, ALLELE0, ALLELE1, CHROM, GENPOS, BETA, OR, BETA_SE, OR_SE,  P, STUDY, Shape, Style,Real_ID)

    }

    if(Model_Reference == T)
    {

      Data <- Data %>% dplyr::select(ID,  OR, BETA, BETA_SE, OR_SE, P, STUDY, Shape, Style, group, Reference)

      Data$Backup_ID <- Data$ID

    }

    if(Model_Reference == F)
    {

      if(Double_Label == T)
      {

       Data$Backup_ID <- Data[[Lab_Col]]

      }else{

        Data$Backup_ID <- Data$ID

      }

      message("Unifying coordinates")

      Data$Real_ID <- Data$ID
      Data$ID <- stringi::stri_c("chr", Data$CHROM, ":", Data$GENPOS, ":", Data$ALLELE0, ":", Data$ALLELE1)
      Data$COORD_Norm <- stringi::stri_c("chr", Data$CHROM, ":", Data$GENPOS, ":", Data$ALLELE0, ":", Data$ALLELE1)
      Data$COORD_Alt <- stringi::stri_c("chr", Data$CHROM, ":", Data$GENPOS, ":", Data$ALLELE1, ":", Data$ALLELE0)

    }

    if(Model_Reference == F)
    {

      #Reset all before next iteration as datasets vary and don't want saved allocation

      Chromosome_Column <- NULL
      PValue_Column     <- NULL
      Position_Column   <- NULL
      SNP_ID_Column     <- NULL
      Ref_Allele_Column <- NULL
      Alt_Allele_Column <- NULL
      BETA_Column       <- NULL
      OR_Column         <- NULL
      SE_Column         <- NULL
      Upper_CI_Column   <- NULL
      Lower_CI_Column   <- NULL
      Test_Statistic    <- NULL

    }

    Combined_Processed_Data <- rbind(Combined_Processed_Data, Data)

  }

  if(Model_Reference == FALSE)

    {

    ALL_STUDIES <- unique(Combined_Processed_Data$STUDY)

    if(!is.null(Lead_Studies))
    {

      messafe("Filtering auto-generated SNPs of interest to particular Studies")

      Combined_Processed_SNPs <- Combined_Processed_SNPs[Combined_Processed_SNPs$STUDY %in% Lead_Studies, ]

    }

    if(is.null(Selected_SNPs))
    {

      Selected_SNPs <- Combined_Processed_SNPs$ID

    }

    if (length(Selected_SNPs) == 0) {

      Combined_Processed_Data <- Combined_Processed_Data

    } else {

      Combined_Processed_Data <- Combined_Processed_Data[Combined_Processed_Data$COORD_Norm %in% Selected_SNPs | Combined_Processed_Data$COORD_Alt %in% Selected_SNPs |  Combined_Processed_Data$ID %in% Selected_SNPs | Combined_Processed_Data$Backup_ID %in% Selected_SNPs , ] # will select even if there is a flip

    }

  }

  if (isTRUE(Model_Reference) && any(Data$Reference != "None", na.rm = TRUE)) {

    message("Formatting model derived processed data for plotting")

    # Get unique references and their corresponding group values where Reference is not "None"
    unique_references_with_group <- Combined_Processed_Data %>%
      dplyr::filter(Reference != "None") %>%
      dplyr::distinct(Reference, group)

    # Create blank rows for each unique Reference value with NA for numeric columns,
    # set the ID to the Reference value, and the group to the corresponding group value
    blank_rows <- data.frame(ID = unique_references_with_group$Reference,
                             OR = 1,
                             BETA = 0,
                             SE = 0,
                             P = 1,
                             STUDY = "",
                             Shape = "square",
                             Style = "normal",
                             group = unique_references_with_group$group,
                             Reference = unique_references_with_group$Reference,
                             stringsAsFactors = FALSE)

    # Combine the original dataframe with the new blank rows
    Combined_Processed_Data <- dplyr::bind_rows(Combined_Processed_Data, blank_rows)

    Combined_Processed_Data$Left_Plot_Value <- ifelse(
      Combined_Processed_Data$ID == Combined_Processed_Data$Reference,
      paste0(Combined_Processed_Data$group, " (Ref: ", Combined_Processed_Data$Reference, ")"),
      Combined_Processed_Data$ID
    )

  }else{

    message("Assigning left sided plot metrics")

    Combined_Processed_Data$Left_Plot_Value <- Combined_Processed_Data$ID

  }

  if(Missings == T & Model_Reference == T)

  {

    message("Simulating missing values from model")

    # Baseline: variables in Model
    baseline_ids <- Combined_Processed_Data %>%
      dplyr::filter(STUDY == "Model", !grepl("\\(Ref:", Left_Plot_Value)) %>%
      dplyr::pull(Left_Plot_Value) %>%
      unique()

    # For each study, only expect baseline_ids
    full_combination <- expand.grid(
      Left_Plot_Value = baseline_ids,
      STUDY = unique(Combined_Processed_Data$STUDY),
      stringsAsFactors = FALSE
    )

    missing_rows <- full_combination %>%
      dplyr::anti_join(Combined_Processed_Data, by = c("Left_Plot_Value", "STUDY"))

    # Add the missing rows with BETA, SE, and P set to 0
    missing_rows <- missing_rows %>%
      dplyr::mutate(
        BETA = 0,
        OR = 1,
        SE = 0,
        P = 1
      )

    #because of blank refs to be skipped
    missing_rows <- missing_rows[!is.na(missing_rows$Left_Plot_Value) & missing_rows$Left_Plot_Value != "" &
                                   !is.na(missing_rows$STUDY) & missing_rows$STUDY != "", ]

    missing_rows <- missing_rows %>%
      dplyr::mutate(
        Shape = "cross",
        Style = "normal"
      )

    # Combine the original dataset with the missing rows
    Combined_Processed_Data <-  dplyr::bind_rows(Combined_Processed_Data, missing_rows)

    Combined_Processed_Data <- dplyr::mutate(Combined_Processed_Data,
                                             Shape = dplyr::case_when(
                                               BETA == 0 & BETA_SE == 0 & P == 1 ~ "cross",
                                               TRUE ~ Shape  # Keeps existing values
                                             )
    )


  }

  if(Missings == T & Model_Reference == F)

  {

    message("Simulating missing values from summary statistics")

    # Extract all unique IDs and STUDY values
    unique_ids <- unique(Combined_Processed_Data$Real_ID)

    unique_studies <- unique(Combined_Processed_Data$STUDY)

    unique_studies <- ALL_STUDIES

    # Create a full combination of all unique IDs and STUDY values
    full_combination <- expand.grid(Real_ID = unique_ids, STUDY = unique_studies)

    # Identify missing combinations (rows not in the original dataset)
    missing_rows <- full_combination %>%
      dplyr::anti_join(Combined_Processed_Data, by = c("Real_ID", "STUDY"))

    # Add the missing rows with BETA, SE, and P set to 0
    missing_rows <- missing_rows %>%
      dplyr::mutate(
        BETA = 0,
        OR = 1,
        SE = 0,
        P = 1
      )

    missing_rows <- missing_rows %>%
      dplyr::mutate(
        Shape = "cross",
        Style = "normal"
      )

    # Combine the original dataset with the missing rows
    Combined_Processed_Data <-  dplyr::bind_rows(Combined_Processed_Data, missing_rows)

    Combined_Processed_Data <- dplyr::mutate(Combined_Processed_Data,
                                             Shape = dplyr::case_when(
                                               BETA == 0 & BETA_SE == 0 & P == 1 ~ "cross",
                                               TRUE ~ Shape  # Keeps existing values
                                             )
    )

  }

  if(Model_Reference == F)
  {

    message(paste0("Processing of matching allele directions to: ", Match_Allele_Study))

    Match_Allele_Study_Clean <- Match_Allele_Study

    reference <- Combined_Processed_Data %>%
      dplyr::filter(!(BETA == 0 & P == 1 & BETA_SE == 0)) %>%   # Exclude placeholder rows
      dplyr::mutate(STUDY_Clean = STUDY) %>%
      dplyr::filter(STUDY_Clean == Match_Allele_Study_Clean) %>%
      dplyr::select(ID, ALLELE0, ALLELE1, COORD_Norm, COORD_Alt, Real_ID) %>%
      dplyr::rename(Ref_ALLELE0 = ALLELE0, Ref_ALLELE1 = ALLELE1)

    reference <- reference %>%
      dplyr::rename_with(~ paste0("REF_", .), everything())

    message("Aligning reference alleles with data...")

    Combined_Processed_Data2 <- Combined_Processed_Data %>%
      dplyr::left_join(reference, by = c("COORD_Norm" = "REF_COORD_Norm"))

    message("Finding alleles that failed to match at first attempt...")

    unmatched <- Combined_Processed_Data2 %>%
      dplyr::filter(is.na(REF_ID)) %>%      dplyr::select(names(Combined_Processed_Data))  # Retain only original columns for the fallback join not joined ones too.

    message("Trying to align reference to failures...")

    fallback <- unmatched %>%
      dplyr::left_join(reference, by = c("COORD_Norm" = "REF_COORD_Alt"))

    message("Removing failures from initial success...and binding second attempt of failures...")

    Combined_Processed_Data2 <- Combined_Processed_Data2 %>%
      dplyr::filter(!is.na(REF_ID)) %>%  # Keep rows matched on COORD_Norm
      dplyr::bind_rows(fallback)  # Add rows matched on COORD_Alt

    Combined_Processed_Data <- Combined_Processed_Data2

    message("Reversing BETA values of alleles where REF and ALT match reference study REF & ALT but are flipped...")

    if(Match_Allele_Direction == TRUE )

    {

      Combined_Processed_Data <- Combined_Processed_Data %>%
        dplyr::mutate(
          BETA_Flipped = !is.na(REF_Ref_ALLELE0) & !is.na(REF_Ref_ALLELE1) &
            ALLELE0 == REF_Ref_ALLELE1 & ALLELE1 == REF_Ref_ALLELE0,
          BETA = dplyr::if_else(BETA_Flipped, BETA * -1, BETA),
          OR = dplyr::if_else(BETA_Flipped & !is.na(OR), 1 / OR, OR),
          No_Match = (is.na(REF_COORD_Norm) & is.na(REF_COORD_Alt))
        )

    }

    # Collect the Real_IDs from rows where Study == Match_Allele_Study
    real_ids_in_match <- Combined_Processed_Data %>%
      dplyr::filter(STUDY == Match_Allele_Study  & !is.na(ID) &
                      P != 1) %>%
      dplyr::pull(Real_ID) %>%
      unique()

    Combined_Processed_Data <- Combined_Processed_Data %>%
      dplyr::mutate(
        use_alt_alleles = STUDY != Match_Allele_Study &
          !is.na(ID) &
          P != 1 &
          !(Real_ID %in% real_ids_in_match),

        COORD_Uni = ifelse(
          use_alt_alleles,
          stringi::stri_c("chr", CHROM, ":", GENPOS, ":", ALLELE0, ":", ALLELE1),
          stringi::stri_c("chr", CHROM, ":", GENPOS, ":", REF_Ref_ALLELE0, ":", REF_Ref_ALLELE1)
        )
      )

    Combined_Processed_Data <- Combined_Processed_Data %>%
      dplyr::group_by(Real_ID) %>%
      dplyr::mutate(COORD_Uni = ifelse(is.na(ID) & P == 1,
                                       dplyr::first(na.omit(COORD_Uni)),
                                       COORD_Uni)) %>%
      dplyr::ungroup()


  }else{

    message("Model References don't require allele adjustment!")

  }

  if(Model_Reference == F)
  {

    Combined_Processed_Data$ID <- Combined_Processed_Data$COORD_Uni

    Combined_Processed_Data$RS <- Combined_Processed_Data$COORD_Uni

    Combined_Processed_Data$Left_Plot_Value <- Combined_Processed_Data$COORD_Uni

    message("Making coordinates look fancy")

    Combined_Processed_Data <- Combined_Processed_Data %>%
      dplyr::mutate(
        Left_Plot_Value = stringr::str_replace(
          Left_Plot_Value,
          "^(chr\\d+):(\\d+):([A-Za-z]+):([A-Za-z]+)$",
          "\\1:\\2(\\3>\\4)"
        )
      )

    if(Double_Label == T)
    {

      Combined_Processed_Data$Left_Plot_Value <- paste0(
        Combined_Processed_Data$Left_Plot_Value,
        "<br><br>(",  # Adds an extra line break to create spacing
        Combined_Processed_Data[[Lab_Col]], ")"
      )

    }

  }

  res <- Combined_Processed_Data

  blank_row <- as.data.frame(lapply(res, function(x) NA))

  res <- dplyr::bind_rows(blank_row, res)

  blank_row <- as.data.frame(lapply(res, function(x) NA))
  res <- dplyr::bind_rows(blank_row, res)

  blank_row <- as.data.frame(lapply(res, function(x) NA))
  res <- dplyr::bind_rows(blank_row, res)

  message("Add necessary dummy rows for spacing")

  res$RS[1] <- "rs99999999"
  res$RS[2] <- "-a-aaarModel"
  res$RS[3] <- "----a"

  if(Pre_Calculated_CIs == F)
  {

    message("Calculating CIs")

    if(Test_Statistic_Choice == "BETA")

    {
      res$LL <- res$BETA - 1.96*(res$BETA_SE)
      res$UL <- res$BETA + 1.96*(res$BETA_SE)

    }else{

      #Assigned OR to BETA earlier.
      res$LL <- res$OR * exp(-1.96 * res$OR_SE)
      res$UL <- res$OR * exp(1.96 * res$OR_SE)

    }

    res$UL <- as.numeric(res$UL)
    res$LL <- as.numeric(res$LL)

  }

  message("Set args for special P notation")

  spaced_star <- paste0(
    "<span style='color:transparent;'>A</span>",
    "*",
    "<span style='color:transparent;'>A</span>"
  )

  res$Special_P <- spaced_star
  res$Special_P[res$P >= 0.05] <- "NS"

  res$P <- sprintf(paste0("%.", P_Label_Resolution, "e"), res$P)

  res$UL <- as.numeric(res$UL)
  res$LL <- as.numeric(res$LL)

  message("Ordering data set and SNPs/covars")

  # Create a vector of letters to use as prefixes
  prefixes <- LETTERS[1:length(Names)]
  postfixes <- LETTERS[1:length(Selected_SNPs)]

  # Modify the res data frame by adding the letter prefix to each dataset
  res <- res %>%
    dplyr::mutate(Study = dplyr::case_when(
      STUDY %in% Names ~ paste(prefixes[match(STUDY, Names)], STUDY, sep = "-"),
      TRUE ~ STUDY  # Default if none match
    ))

  if(Double_Label == T)
  {

    res <- res %>%
      dplyr::mutate(Backup_Single = stringr::str_extract(Backup_ID, "rs\\d+"))

  }else{

    res <- res %>%
      dplyr::mutate(Backup_Single = Backup_ID)

  }


  message("Ordering")

  if(Model_Reference == F)
  {

    res <- res %>%
      dplyr::mutate(
        RS = dplyr::case_when(
          RS %in% Selected_SNPs | Backup_Single %in% Selected_SNPs | Real_ID %in% Selected_SNPs ~ paste(
            postfixes[match(
              ifelse(
                RS %in% Selected_SNPs,
                RS,
                ifelse(Backup_Single %in% Selected_SNPs, Backup_Single, Real_ID)
              ),
              Selected_SNPs
            )],
            RS,
            sep = "-"
          ),
          TRUE ~ RS
        )
      )

  }else{

    res <- res %>%
      dplyr::mutate(RS = dplyr::case_when(
        group %in% Selected_SNPs ~ paste(postfixes[match(group, Selected_SNPs)], RS, sep = "-"),
        TRUE ~ RS  # Default if none match
      ))

  }

  if(Model_Reference == TRUE)
  {

  message("Adjusting Model version for ref values")

  res <- res %>%
    dplyr::mutate(
      RS = ifelse(
        grepl("Ref:", Left_Plot_Value),
        sub("-", "-11A", RS),
        RS
      )
    )

  }

  res$BETA <- as.numeric(res$BETA)

  res$RS <- paste0(res$RS,res$Study)

  #Adjust back these placeholders
  res$RS[res$RS ==  "rs99999999NA"  ] <- "-aaa-rs99999999"
  res$RS[res$RS ==   "-a-aaarModelNA" ] <- "-a-aaarModel"
  res$RS[res$RS ==  "----aNA"  ] <- "---a"
  res$RS[res$RS ==   "zzzNA" ] <- "zzz"

  unique_study <- unique(na.omit(res$Study))

  res <- res %>%
    dplyr::arrange(desc(RS)) %>%  # Order by RS in descending order
    dplyr::mutate(Overall_Row_Number = dplyr::row_number())  # Create a new column with row numbers

  #In model ref need to remove dummy " for ref cols
  unique_study <- unique_study[unique_study != ""]

  res2 <- res %>%
    dplyr::arrange(RS)

  res2 <- res2 %>% dplyr::arrange(desc(dplyr::row_number()))

  res$Left_Plot_Value[res$RS == "-aaa-rs99999999"] <- ""
  res$Left_Plot_Value[res$RS == "---a"] <- ""
  res$Left_Plot_Value[res$RS == "-a-aaarModel"] <- Left_Title
  res$Plot_Value[res$RS == "-a-aaarModel"] <- res$Overall_Row_Number[res$RS == "-a-aaarModel"]

  message("Assigning required value to test statistic for plotting")

  if(Test_Statistic_Choice == "BETA")
  {

    res$SE <- res$BETA_SE
    res$BETA <- res$BETA

  }else{

    res$SE <- res$OR_SE
    res$BETA <- res$OR

  }

  MAX <- (.Machine$double.xmin) #may need in future

  res$outer <- res$Left_Plot_Value

  message("Creating buffers")

  stat_buf <- Left_Buffer

  invisible_z <- paste(rep("z", stat_buf), collapse = "")
  transparent_z_span <- paste0("<span style='color:transparent;'>", invisible_z, "</span>")

  # Invisible spacer
  spacer_unit <- "Z"
  invisible_padding <- paste(rep(transparent_z_span, stat_buf), collapse = "")

  num_z <- 2

  # Create transparent spacer
  zzz_span <- sprintf(
    "<span style='color:transparent;'>%s</span>",
    strrep("Z", num_z)
  )

  zzz_span <- transparent_z_span

  invis_bracket <- "<span style='color:transparent;'>)</span>"

  # Model row index
  model_row <- which(res$inner == "-a-aaarModel")[1]

  # Check if any non-model row has a ')' before the first <br>
  any_has_bracket <- any(sapply(res$outer[res$inner != "-a-aaarModel"], function(x) {
    if (grepl("<br", x, fixed = TRUE)) {

      stringr::str_detect(strsplit(x, "<br>", fixed = TRUE)[[1]][1], stringr::fixed(")"))

    } else {

      stringr::str_detect(x, stringr::fixed(")"))

    }

  }))

  # Check if model row does NOT have a ')' before first <br>
  model_lacks_bracket <- {
    x <- res$outer[model_row]
    if (grepl("<br", x, fixed = TRUE)) {

      !stringr::str_detect(strsplit(x, "<br>", fixed = TRUE)[[1]][1], stringr::fixed(")"))

    } else {

      !stringr::str_detect(x, stringr::fixed(")"))

    }
  }

  # Now build transformed version of res$outer
  res$outer <- vapply(seq_along(res$outer), function(i) {

    x <- res$outer[i]

    # Always start with leading transparent zzz
    out <- zzz_span

    if (grepl("<br", x, fixed = TRUE)) {

      parts <- strsplit(x, "<br>", fixed = TRUE)[[1]]

      # Model row may need invisible bracket in the first part
      if (i == model_row && any_has_bracket && model_lacks_bracket) {

        parts[1] <- paste0(parts[1], invis_bracket)

      }

      # Add zzz after every part *before* a <br> (i.e. all but the last part)
      if (length(parts) > 1) {

        for (j in 1:(length(parts) - 1)) {
          parts[j] <- paste0(parts[j], zzz_span)

        }
      }

      out <- paste0(out, paste(parts, collapse = "<br>"), zzz_span)
    } else {

      # No <br>: optionally add bracket to model row if needed
      if (i == model_row && any_has_bracket && model_lacks_bracket) {

        x <- paste0(x, invis_bracket)

      }
      out <- paste0(out, x, zzz_span)

    }

    out
  }, character(1))

  message("Defining styling function")

  style_label <- function(label, style) {

    if (is.na(style) || style == "normal") {

      return(label)

    } else if (style == "bold") {

      return(paste0("<b>", label, "</b>"))

    } else if (style == "underline") {

      # Unicode combining underline (works!)
      return(paste0(paste0(strsplit(label, "")[[1]], collapse = "\u0332"), "\u0332"))

    } else {

      return(label)

    }
  }

  message("Assigning final plot values")

  res$inner <- res$RS

  res <- dplyr::mutate(res, num = dplyr::row_number())

  if(Test_Statistic_Choice == "BETA")
  {

    res$beta <- res$BETA
    res$beta <- as.numeric(res$beta)
    res$se <- res$BETA_SE
    res$se <- as.numeric(res$se)

  }else{

    res$beta <- res$OR
    res$beta <- as.numeric(res$beta)
    res$se <- res$OR_SE
    res$se <- as.numeric(res$se)

  }

  res$p <- as.numeric(res$P)

  if(Special_P == TRUE)
  {

    res$p <- res$Special_P

  }

  res$x = 1
  res$y = 1

  res<- res[!(res$inner %in% c("-aaa-rs99999999", "---a")), ]
  res$inner <- as.factor(res$inner)
  res$outer <- as.factor(res$outer)

  message("Computing CIs")

  if(Test_Statistic_Choice == "BETA")
  {

    res$ci_lower <- res$beta - 1.96 * res$se
    res$ci_upper <- res$beta + 1.96 * res$se

  }else{

    res$ci_lower <- res$beta - exp(1.96 * res$se)
    res$ci_upper <- res$beta + exp(1.96 * res$se)

  }

  message("Rounding Text Labels")

  if(Special_P == FALSE)
  {
    # Scientific notation with 'x10^' style
    res$p_label <- formatC(res$p, format = "e", digits = P_Label_Resolution)

  }

  res$beta <- round(res$beta, digits = BETA_Label_Resolution)
  res$se <- round(res$se, digits = SE_Label_Resolution)
  res$ci_lower <- round(res$ci_lower, digits = CI_Label_Resolution)
  res$ci_upper <- round(res$ci_upper, digits = CI_Label_Resolution)

  message("Formatting Text Labels")

  res$beta[res$inner == "-a-aaarModel"] <- NA

  res$p[res$inner == "-a-aaarModel"] <- NA

  res$beta_label <- sprintf(paste0("%.", BETA_Label_Resolution, "f"), res$beta)

  formatted <- sprintf("%.2f", res$beta[match(levels(res$inner), res$inner)])
  unique_levels <- make.unique(rev(formatted))

  res$beta_label <- factor(res$beta_label, levels = unique_levels)

  # If outer is not a character yet, convert it
  res$outer <- as.character(res$outer)

  # Define buffer size
  stat_buf <- Right_Buffer

  invisible_z <- paste(rep("z", stat_buf), collapse = "")
  transparent_z_span <- paste0("<span style='color:transparent;'>", invisible_z, "</span>")

  # Invisible spacer
  spacer_unit <- "Z"
  invisible_padding <- paste(rep(transparent_z_span, stat_buf), collapse = "")

  # Wrap function
  wrap_with_padding <- function(vec) {
    paste0(invisible_padding, vec, invisible_padding)
  }

  display_labels <- rep("", length(res$inner))

  res$ci_upper[res$inner == "-a-aaarModel"] <- NA
  res$ci_lower[res$inner == "-a-aaarModel"] <- NA

  res$beta_label <- sprintf(paste0("%.", BETA_Label_Resolution, "f"), res$beta)

  res$beta_label[res$inner == "-a-aaarModel"] <- BETA_Label_Title

  res$beta_label[(is.na(res$UL)) & res$inner != "-a-aaarModel" ] <- ""

  res$beta_label[ (is.na(res$ALLELE0) | is.na(res$UL)) & res$inner != "-a-aaarModel" ] <- ""

  if(Test_Statistic_Choice == "BETA")
  {

    # Add invisible minus (visually hidden but layout-preserving)
    res$beta_label <- mapply(function(label, inner) {

      if (label != "" && inner != "-a-aaarModel" && !startsWith(label, "-")) {

        paste0("<span style='color:transparent;'>-</span>", label)
      } else {

        label

      }
    }, res$beta_label, res$inner, USE.NAMES = FALSE)

  }

  # Tag-aware underline version of style_label()
  style_label_safe <- function(label, style) {

    if (is.na(style) || style == "normal") {

      return(label)

    } else if (style == "bold") {

      return(paste0("<b>", label, "</b>"))

    } else if (style == "underline") {

      # Split into HTML and non-HTML segments
      parts <- unlist(strsplit(label, "(<[^>]+>)", perl = TRUE))
      parts <- vapply(parts, function(part) {


        if (grepl("^<", part)) {
          # HTML tag -> leave it untouched
          part

        } else if (nzchar(part)) {

          # Visible text -> add combining underline
          paste0(paste0(strsplit(part, "")[[1]], collapse = "\u0332"), "\u0332")

        } else {

          part

        }

      }, character(1))
      paste0(parts, collapse = "")

    } else {

      return(label)

    }
  }

  # Use the safe version here
  res$beta_label <- setNames(
    mapply(style_label_safe, res$beta_label, res$Style, SIMPLIFY = TRUE),
    res$beta
  )

  res$beta_label <- wrap_with_padding(res$beta_label)

  display_labels_beta <- setNames(res$beta_label, res$inner)

  res$se_label <- sprintf(paste0("%.", SE_Label_Resolution, "f"), res$se)

  res$se_label <- paste0("", res$se_label, "")

  res$se_label[(is.na(res$UL)) & res$inner != "-a-aaarModel" ] <- ""

  res$se_label[is.na(res$ALLELE0) & res$inner != "-a-aaarModel" ] <- ""

  res$se_label[res$inner == "-a-aaarModel"] <- SE_Label_Title


  if ("underline" %in% res$Style) {

    res$se_label <- paste0("\u200B", res$se_label)

  }

  #Build display label mappings
  res$se_label  <- setNames(
    mapply(style_label, res$se_label, res$Style, SIMPLIFY = TRUE),
    res$se
  )

  res$se_label <- wrap_with_padding(res$se_label)

  display_labels_se <- setNames(res$se_label, res$inner)

  if(Special_P == FALSE){

    # Scientific notation
    res$p_label <- sprintf(paste0("%.", P_Label_Resolution, "e"), res$p)

  }else{

    res$p_label <- res$p
  }

  res$p_label[(is.na(res$UL)) & res$inner != "-a-aaarModel" ] <- ""

  res$p_label <- paste0("", res$p_label, "")

  res$p_label[is.na(res$ALLELE0) & res$inner != "-a-aaarModel" ] <- ""

  res$p_label[res$inner == "-a-aaarModel"] <- P_Label_Title

  if ("underline" %in% res$Style) {

    res$p_label <- paste0("\u200B", res$p_label)

  }

  res$p_label  <- setNames(
    mapply(style_label_safe, res$p_label, res$Style, SIMPLIFY = TRUE),
    res$p
  )

  # Overwrite the raw display labels directly
  res$p_label <- wrap_with_padding(res$p_label)

  display_labels_p <- setNames(res$p_label, res$inner)

  # Ensure numeric
  res$beta <- as.numeric(res$beta)
  res$se   <- as.numeric(res$se)

  # Format values
  res$beta_fmt <- sprintf(paste0("%.", BETA_Label_Resolution, "f"), res$beta)

  res$se_fmt <- sprintf(paste0("%.", SE_Label_Resolution, "f"), res$se)

  res$beta_fmt[(is.na(res$UL)) & res$inner != "-a-aaarModel" ] <- ""
  res$se_fmt[(is.na(res$UL)) & res$inner != "-a-aaarModel" ] <- ""

  # Add invisible minus (visually hidden but layout-preserving)
  res$beta_fmt <- mapply(function(label, inner) {

    if (label != "" && inner != "-a-aaarModel" && !startsWith(label, "-")) {

      paste0("<span style='color:transparent;'>-</span>", label)

    } else {

      label

    }
  }, res$beta_fmt, res$inner, USE.NAMES = FALSE)

  # Combine beta (se)
  res$beta_se_label <- paste0(res$beta_fmt, " (", res$se_fmt, ")")

  # Handle model header
  res$beta_se_label[res$inner == "-a-aaarModel"] <- BETA_SE_Label_Title

  # Blank for missing SNPs
  res$beta_se_label[(is.na(res$ALLELE0)) & res$inner != "-a-aaarModel" ] <- ""
  res$beta_se_label[res$beta_se_label == " ()" & res$inner != "-a-aaarModel" ] <- ""

  res$beta_se_label  <- setNames(
    mapply(style_label_safe, res$beta_se_label, res$Style, SIMPLIFY = TRUE),
    res$beta_se_label
  )

  # Overwrite the raw display labels directly
  res$beta_se_label       <- wrap_with_padding(res$beta_se_label)

  # Name the labels for faceting
  display_labels_beta_se <- setNames(res$beta_se_label, res$inner)

  # Ensure numeric
  res$beta     <- as.numeric(res$beta)
  res$ci_lower <- as.numeric(res$ci_lower)
  res$ci_upper <- as.numeric(res$ci_upper)

  # Format values to 2 decimal places
  res$beta_fmt <- sprintf(paste0("%.", BETA_Label_Resolution, "f"), res$beta)

  res$beta_fmt[(is.na(res$UL)) & res$inner != "-a-aaarModel" ] <- ""

  res$ci_lower_fmt <- sprintf(paste0("%.", CI_Label_Resolution, "f"), res$ci_lower)

  res$ci_lower_fmt[(is.na(res$UL)) & res$inner != "-a-aaarModel" ] <- ""

  res$ci_upper_fmt <- sprintf(paste0("%.", CI_Label_Resolution, "f"), res$ci_upper)

  res$ci_upper_fmt[(is.na(res$UL)) & res$inner != "-a-aaarModel" ] <- ""

  # Add invisible minus (visually hidden but layout-preserving)
  res$beta_fmt <- mapply(function(label, inner) {

    if (label != "" && inner != "-a-aaarModel" && !startsWith(label, "-")) {

      paste0("<span style='color:transparent;'>-</span>", label)

    } else {

      label

    }

  }, res$beta_fmt, res$inner, USE.NAMES = FALSE)

  # Add invisible minus (visually hidden but layout-preserving)
  res$ci_lower_fmt <- mapply(function(label, inner) {

    if (label != "" && inner != "-a-aaarModel" && !startsWith(label, "-")) {

      paste0("<span style='color:transparent;'>-</span>", label)

    } else {

      label

    }
  }, res$ci_lower_fmt, res$inner, USE.NAMES = FALSE)

  # Add invisible minus (visually hidden but layout-preserving)
  res$ci_upper_fmt <- mapply(function(label, inner) {

    if (label != "" && inner != "-a-aaarModel" && !startsWith(label, "-")) {

      paste0("<span style='color:transparent;'>-</span>", label)

    } else {

      label
    }
  }, res$ci_upper_fmt, res$inner, USE.NAMES = FALSE)

  # Combine into label
  res$beta_ci_label <- paste0("", res$beta_fmt, " (", res$ci_lower_fmt, ", ", res$ci_upper_fmt, ")")

  # Handle model header
  res$beta_ci_label[res$inner == "-a-aaarModel"] <- BETA_CI_Label_Title

  # Blank for missing SNP rows
  res$beta_ci_label[is.na(res$ALLELE0) & res$inner != "-a-aaarModel"] <- ""

  res$beta_ci_label[res$beta_ci_label == " (, )"] <- ""

  res$beta_ci_label  <- setNames(
    mapply(style_label_safe, res$beta_ci_label, res$Style, SIMPLIFY = TRUE),
    res$beta_ci_label
  )

  res$beta_ci_label       <- wrap_with_padding(res$beta_ci_label)

  # Named vector for use in facet labels
  display_labels_beta_ci <- setNames(res$beta_ci_label, res$inner)

  # Combine into label
  res$ci_label <- paste0(""," (", res$ci_lower_fmt, ", ", res$ci_upper_fmt, ")")

  # Handle model header
  res$ci_label[res$inner == "-a-aaarModel"] <- CI_Label_Title

  res$ci_label[res$ci_label == " (, )"] <- ""

  # Blank for missing SNP rows
  res$ci_label[is.na(res$ALLELE0) & res$inner != "-a-aaarModel"] <- ""

  res$ci_label  <- setNames(
    mapply(style_label_safe, res$ci_label, res$Style, SIMPLIFY = TRUE),
    res$ci_label
  )

  # Overwrite the raw display labels directly
  res$ci_label <- wrap_with_padding(res$ci_label)

  # Named vector for use in facet labels
  display_labels_ci <- setNames(res$ci_label, res$inner)

  # User input: number of transparent padding A's between labels
  n_padding_A <- 10
  invisible_A_pad <- paste0("<span style='color:transparent;'>", strrep("A", n_padding_A), "</span>")

  #Strip HTML and compute visible widths
  strip_html <- function(x) gsub("<[^>]*>", "", x)
  visible_labels <- sapply(res$beta_se_label, strip_html)
  visible_widths <- nchar(visible_labels, type = "width")

  #Find max and second max widths
  sorted_widths <- sort(unique(visible_widths), decreasing = TRUE)
  max_width <- sorted_widths[1]
  second_max <- if (length(sorted_widths) >= 2) sorted_widths[2] else max_width

  #Compute the difference in width
  width_diff <- max_width - second_max

  #Create invisible spacer (e.g., transparent underscores)
  transparent_spacer <- paste(rep("<span style='color:transparent;'>_</span>", width_diff), collapse = "")

  #Append the extra_spacer only where needed
  res$extra_spacers <- ifelse(res$inner == "-a-aaarModel", transparent_spacer, "")

  #Combine all into the desired format
  res$beta_se_ci_label <- paste0(
    res$beta_se_label,
    invisible_A_pad,
    res$extra_spacers,
    res$beta_ci_label

  )

  selected_cols <- order

  strip_html <- function(x) gsub("<[^>]*>", "", x)

  res$beta_se_ci_label[is.na(res$ALLELE0) & res$inner != "-a-aaarModel"] <- ""

  res$beta_se_ci_label  <- wrap_with_padding(res$beta_se_ci_label)

  #Create display vector
  display_labels_beta_se_ci <- setNames(res$beta_se_ci_label, res$inner)

  res$outer <- as.factor(res$outer)

  res$special_row <- res$inner == "-a-aaarModel"

  res <- res[order(!res$special_row, res$inner), ]

  res$special_row <- NULL

  outer_levels <- unique(res$outer)
  res$outer <- factor(res$outer, levels = outer_levels)

  inner_levels <- sort(unique(res$inner))

  inner_levels <- c("-a-aaarModel", setdiff(inner_levels, "-a-aaarModel"))

  res$inner <- factor(res$inner, levels = inner_levels)

  message("Creating left plot")

  left <- ggplot2::ggplot(res, ggplot2::aes(x, y)) +
    ggplot2::theme_void()+
    ggplot2::geom_blank()+
    ggplot2::scale_y_discrete(position = "bottom")+
    ggh4x::facet_nested(
      outer + inner ~ .,
      scales = "free_y",
      space = "free_y",
      strip = ggh4x::strip_nested(size = "variable"),
      labeller = ggplot2::labeller(inner = setNames(display_labels, levels(res$inner)))
    ) +
      ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = "white", linewidth = 0),
      panel.background = ggplot2::element_rect(fill = "white", colour = "white", linewidth = 0),
      panel.grid.major = ggplot2::element_line(colour = "white", linewidth = 0),
      panel.grid.minor = ggplot2::element_line(colour = "white", linewidth = 0),
      axis.line = ggplot2::element_line(colour = "white", linewidth = 0),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "white", linewidth = 0, fill = NA),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_text(size = 36, colour = "navy", angle = 0, hjust = 1),
      strip.text.y.right = ggtext::element_markdown(size = Left_Size, colour = "black", angle = 0, hjust = 1),
      strip.background.y = ggplot2::element_rect(fill = Left_Fill, colour = Grid_Colour_Left, linewidth = Grid_Thickness),
      strip.background.y.right = ggplot2::element_rect(fill = "lightblue", colour = "lightblue", linewidth = 2),
      strip.placement = "outside",
      panel.spacing = ggplot2::unit(0.0, "cm"),
      strip.switch.pad.grid = ggplot2::unit(0, "cm"),
      plot.margin = ggplot2::unit(c(0.0, 0.0, 0.0, 0.0), "cm")
    )

  #Ensure character for label keys
  res$se_label <- as.character(res$se_label)
  res$beta_label <- as.character(res$beta_label)
  res$p_label <- as.character(res$p_label)

  message("Ordering")

  #Assign random number of spaces (0 to 5) per row
  set.seed(123)
  res$outer_sec <- sapply(sample(0:5, nrow(res), replace = TRUE), function(n) paste(rep(" ", n), collapse = ""))

  res$inner <- factor(res$inner, levels = unique(res$inner))

  res$num <- seq_len(nrow(res))

  res$outer_sec <- as.character(res$outer_sec)

  res$num <- factor(res$num, levels = res$num)

  display_labels_outer_sec <- setNames(res$outer_sec, res$num)

  all_labels <- list(
    ci_label = display_labels_ci,
    beta_label = display_labels_beta,
    se_label = display_labels_se,
    beta_se_label = display_labels_beta_se,
    beta_ci_label = display_labels_beta_ci,
    p_label = display_labels_p
  )

  show_flags <- c(
    ci_label = CI_Label,
    beta_label = BETA_Label,
    se_label = SE_Label,
    beta_se_label = BETA_SE_Label,
    beta_ci_label = BETA_CI_Label,
    p_label = P_Label
  )

  selected_labels <- all_labels[show_flags]

  selected_labels$num <- display_labels_outer_sec

  labeller_obj <- do.call(ggplot2::labeller, selected_labels)

  selected_vars <- names(show_flags[show_flags])

  # Map to show_flags names
  name_map <- c(
    P_Label = "p_label",
    BETA_Label = "beta_label",
    SE_Label = "se_label",
    BETA_SE_Label = "beta_se_label",
    BETA_CI_Label = "beta_ci_label",
    CI_Label = "ci_label"
  )

  # Replace with correct names
  order_matched <- name_map[order]

  order <- order_matched

  # Always start with "num"
  facet_vars <- c("num", order)

  facet_formula_str <- paste(facet_vars, collapse = " + ")

  facet_formula <- as.formula(paste(facet_formula_str, "~ ."))

  message("Creating right plot")

  right <- ggplot2::ggplot(res, ggplot2::aes(x, y)) +
    ggplot2::theme_void() +
    ggplot2::geom_blank() +
    ggplot2::scale_y_discrete(position = "bottom") +
    ggh4x::facet_nested(
      facet_formula,
      switch = "y",
      bleed = FALSE,
      scales = "free_y",
      space = "free_y",
      strip = ggh4x::strip_nested(size = "variable"),
      labeller = labeller_obj
    ) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = "white", linewidth = 0),
      panel.background = ggplot2::element_rect(fill = "white", colour = "white", linewidth = 0),
      panel.grid.major = ggplot2::element_line(colour = "white", linewidth = 0),
      panel.grid.minor = ggplot2::element_line(colour = "white", linewidth = 0),
      axis.line = ggplot2::element_line(colour = "white", linewidth = 0),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "white", linewidth = 0, fill = NA),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      strip.text.y.left = ggtext::element_markdown(size = Right_Size, colour = "black",  angle = 0, hjust = 0),
      strip.background.y = ggplot2::element_rect(fill = Right_Fill, colour = Grid_Colour_Right, linewidth = Grid_Thickness),
      strip.background.y.right = ggplot2::element_rect(fill = "lightblue", colour = "lightblue", linewidth = 2),
      strip.placement = "outside",
      panel.spacing = ggplot2::unit(0.0, "cm"),
      strip.switch.pad.grid = ggplot2::unit(0, "cm"),
      plot.margin = ggplot2::unit(c(0.0, 0.0, 0.0, 0.0), "cm")
    )


  if(Strips == TRUE )

  {

  message("Determining strip placement")

  }

  # Ensure outer is a factor in the desired order
  res$outer <- factor(res$outer, levels = unique(res$outer))

  # Identify the outer value for the SNP row
  snp_outer <- res$outer[res$inner == "-a-aaarModel"][1]

  # Get all unique outer levels
  outer_levels <- levels(res$outer)

  # Place SNP outer first, then alternate others
  ordered_outers <- c(snp_outer, setdiff(outer_levels, snp_outer))

  if(is.null(Strip_Colour))
  {

    outer_colors <- rep(c("goldenrod", "white"), length.out = length(ordered_outers))

  }else{

    outer_colors <- rep(c(Strip_Colour, "white"), length.out = length(ordered_outers))

  }

  outer_colors[1] <- "#bebebe"

  # Create lookup table
  outer_fill_map <- setNames(outer_colors, ordered_outers)

  # Apply to df
  res <- res %>%
    dplyr::mutate(
      fills = outer_fill_map[as.character(outer)],
      cols  = outer_fill_map[as.character(outer)]
    )

  res2 <- res[res$inner == "-a-aaarModel", ]

  message("Simulating NULL line")

  if(Test_Statistic_Choice == "OR")

  {

    vline_df <- res |>
      dplyr::distinct(outer, inner) |>
      dplyr::arrange(outer, inner) |>
      dplyr::mutate(
        vline = ifelse(trimws(inner) != "-a-aaarModel", 1, NA_real_)
      )

  }else{

    vline_df <- res |>
      dplyr::distinct(outer, inner) |>
      dplyr::arrange(outer, inner) |>
      dplyr::mutate(
        vline = ifelse(trimws(inner) != "-a-aaarModel", 0, NA_real_)
      )
  }

  message("Creating middle plot")

  mid <- ggplot2::ggplot(res, ggplot2::aes(x = beta, y )) +
    ggplot2::theme_void()+
    ggplot2::facet_grid(
      outer + inner ~ .,
      scales = "free_y",
      space = "free_y",
      labeller = ggplot2::labeller(inner = setNames(display_labels, levels(res$inner)))
    ) + ggplot2::geom_rect(data = res2, fill = "transparent",  colour = "black" ,xmin = -Inf,xmax = Inf,
              ymin = -Inf,ymax = Inf,alpha = 0.3,
              linewidth = Top_Grid_Thickness)+
    ggplot2::geom_point(ggplot2::aes(x=beta, y, colour = STUDY),  shape=res$Shape, size= Shape_Size, stroke = Shape_Thickness) +
    ggplot2::scale_y_discrete(position = "bottom")+
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lower, xmax = ci_upper), height = CI_Line_Height, size = CI_Line_Size, color = "black") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 10, color = "black"),
      axis.title.x = ggplot2::element_text(size = 12,  margin = ggplot2::margin(t = 0)),
      axis.line.x = ggplot2::element_line(colour = "black", size = 0),
      axis.ticks.x = ggplot2::element_line(colour = "black"),
      plot.background = ggplot2::element_rect(fill = "white", colour = "white", linewidth = 0),
      panel.grid.major = ggplot2::element_line(colour = "white", linewidth = 0),
      panel.grid.minor = ggplot2::element_line(colour = "white", linewidth = 0),
      axis.ticks.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "white", linewidth = 0, fill = NA),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_blank(),
      strip.text.y.right = ggplot2::element_blank(),
      strip.placement = "outside",
      panel.spacing = ggplot2::unit(0.0, "cm"),
      strip.switch.pad.grid = ggplot2::unit(0, "cm"),
      plot.margin = ggplot2::unit(c(0.0, 0.0, 0.0, 0.0), "cm")
    )

  if(is.null(Null_Line_Colour))
  {
    Null_Line_Colour <- "red"
  }

  if(is.null(Null_Line_Type))
  {
    Null_Line_Type <- "dotted"
  }

  if(is.null(Null_Line_Thickness))
  {
    Null_Line_Thickness <- 1
  }

  message("Adding null line")

  mid <- mid + ggplot2::geom_vline(
    data = vline_df,
    ggplot2::aes(xintercept = vline),
    colour = Null_Line_Colour, linetype = Null_Line_Type, linewidth = Null_Line_Thickness
  )

  if(Strips == TRUE)

  {
    message("Colouring strips")

    mid <- mid + ggplot2::geom_rect(data = res, fill = res$fills,  colour = res$cols ,xmin = -Inf,xmax = Inf,
                           ymin = -Inf,ymax = Inf,alpha = 0.1,
                           linewidth = 0 )

  }

  message("Creating Key")

  values <- setNames(Data_Set_Colours, Names)
  labels <- setNames(Names, Names)

  mid <- mid + ggplot2::scale_color_manual(
    values = values,
    labels = labels,
    breaks = Names

  )

  if(Legend_On == TRUE)

  {

    message("Creating legend")

    n_groups <- length(unique(res$STUDY))  # or group variable

    mid <- mid + ggplot2::guides(
      color = ggplot2::guide_legend(
        title = Legend_Title,
        override.aes = list(shape = Shapes)
      )
    )

  }

  message("Formatting Legend and X Axis")

  mid <- mid + ggplot2::theme(legend.text = ggplot2::element_text(size = Legend_Text_Size, color ="black"),
                              legend.title = ggplot2::element_text(size = Legend_Title_Size, color ="black"),
                              axis.text.x = ggplot2::element_text(size = X_Axis_Text_Size, color ="black"),
                              axis.title.x = ggplot2::element_text(margin = ggplot2::margin(), color = "black", size = X_Axis_Title_Size ))

  mid <- mid + ggplot2::labs(x=X_Axis_Title, y="YL")


  mid <- mid + ggplot2::theme(
    axis.ticks.length.x  = ggplot2::unit(0.2,"cm"),
    axis.text.x = ggplot2::element_text( size = X_Axis_Text_Size, vjust = -1),
    legend.margin=ggplot2::margin()
  )


  if(X_Axis_Label == T)
  {

    X_Axis_Title <- paste0("\n", X_Axis_Title)
    mid <- mid + ggplot2::xlab(X_Axis_Title)

  }

  if(X_Axis_Label == F)
  {

    mid <- mid + ggplot2::xlab(" ")

  }

  message("Determining X axis format")

  maxcalc <-  max(res$ci_upper, na.rm = T)
  mincalc <-  min(res$ci_lower, na.rm = T)

  # Determine range width
  range_width <- max(abs(mincalc), abs(maxcalc))  # Symmetrical range

  .ensure_two_sides <- function(breaks, null, minlim, maxlim, step = NULL, log_scale = FALSE) {
    if (!length(breaks)) return(breaks)

    has_below <- any(breaks < null)
    has_above <- any(breaks > null)

    if (has_below && has_above) return(breaks)

    if (log_scale) {

      if (is.null(step) || !is.finite(step) || step <= 0) step <- 0.30103  # ~log10(2)
      cand_below <- 10^(log10(null) - step)
      cand_above <- 10^(log10(null) + step)

    } else {

      if (is.null(step) || !is.finite(step) || step <= 0) step <- diff(range(breaks)) / 6
      cand_below <- null - step
      cand_above <- null + step

    }

    if (!has_below && cand_below >= minlim && cand_below < null) {

      breaks <- sort(c(breaks, cand_below))

    }

    if (!has_above && cand_above <= maxlim && cand_above > null) {

      breaks <- sort(c(breaks, cand_above))

    }

    breaks
  }

  mid <- (function(
    mid, mincalc, maxcalc,
    Test_Statistic_Choice = c("BETA","OR"),
    Axis_Buffer = 0.02,
    X_Axis_Separation = NULL,
    X_Axis_Text_Resolution = NULL,
    Null_Buffer = 0.0
  ){

    Test_Statistic_Choice <- match.arg(Test_Statistic_Choice)

    # Helper for nice step rounding
    nice_steps <- function(step_raw) {

      cands <- c(1, 2, 2.5, 5) * 10^floor(log10(step_raw))
      cands[which.min(abs(step_raw - cands))]

    }
    dec_for_step <- function(step) if (step >= 1) 0 else abs(floor(log10(step)))

    #Null center value
    null_center <- if (Test_Statistic_Choice == "BETA") 0 else 1

    #Linear scale
    range_width <- max(abs(mincalc - null_center), abs(maxcalc - null_center))

    #Auto separation if user inputs missing
    if (is.null(X_Axis_Separation) || is.na(X_Axis_Separation)) {
      raw_step <- (2 * range_width) / 6
      X_Axis_Separation <- nice_steps(raw_step)
    }
    if (is.null(X_Axis_Text_Resolution) || is.na(X_Axis_Text_Resolution)) {
      X_Axis_Text_Resolution <- dec_for_step(X_Axis_Separation)
    }

    # Breaks
    min_axis <- null_center - range_width
    max_axis <- null_center + range_width
    breaks <- seq(min_axis, max_axis, by = X_Axis_Separation)

    # Null buffer
    if (Null_Buffer > 0) {
      keep <- (abs(breaks - null_center) >= Null_Buffer) | (abs(breaks - null_center) < .Machine$double.eps)
      breaks <- breaks[keep]
    }

    # Ensure at least one tick each side of center
    breaks <- .ensure_two_sides(breaks, null_center, mincalc, maxcalc, X_Axis_Separation, log_scale = FALSE)

    # Clip to CI range
    breaks <- breaks[breaks >= mincalc & breaks <= maxcalc]

    # Labels
    fmt <- paste0("%.", X_Axis_Text_Resolution, "f")
    labels <- sprintf(fmt, breaks)

    # Axis padding
    buffer <- (maxcalc - mincalc) * Axis_Buffer

    mid <- mid + ggplot2::scale_x_continuous(
      limits = c(mincalc - buffer, maxcalc + buffer),
      breaks = breaks,
      labels = labels,
      expand = c(0, 0)
    )

    mid
  })(mid, mincalc, maxcalc, Test_Statistic_Choice, Axis_Buffer,
     X_Axis_Separation, X_Axis_Text_Resolution, Null_Buffer)


  if(Legend_On == FALSE)
  {

    mid <- mid +  ggplot2::theme(legend.position = "none")

  }

  if(Legend_On == TRUE)
  {

    message("Positioning Legend")

      mid <- mid +
      ggplot2::theme(
        legend.position = "bottom",
        legend.box = "horizontal",

      )

  }

  message("Final overall formatting and combining")

  if ("underline" %in% res$Style) {

    message("Adding underline")

    mid <- mid + ggplot2::theme(axis.line.x = ggplot2::element_line(linewidth  = 0.1) )
    combined <- invisible(left + mid + right + patchwork::plot_layout(widths = c(0.05, 0.9, 0.05)))
    g <- combined

  }
  else if(Grid_Lines_On == FALSE)
  {

    message("Stripping gird lines")

    mid <- mid + ggplot2::theme(axis.line.x = ggplot2::element_line(linewidth  = 0.1) )
    combined <- left + mid + right + patchwork::plot_layout(widths = c(0.05, 0.9, 0.05))
    g <- combined

  }
  else{

    message("Combining")

    #Combine plots with patchwork
    combined <- left + mid + right + patchwork::plot_layout(widths = c(0.05, 0.9, 0.05))

    #Convert to grob
    g <- patchwork::patchworkGrob(combined)

    #Identify panel rows
    mid_panel_rows <- grep("^panel", g$layout$name)
    mid_layout <- g$layout[mid_panel_rows, ]

    #Find the horizontal column span for the mid plot
    mid_cols <- range(mid_layout$l)

    #Find the bottom-most row (i.e., bottom-most facet)
    bottom_row <- max(mid_layout$t)

    #Add green line inside panel, slightly inset to avoid strip overhang
    g <- gtable::gtable_add_grob(
      g,
      grid::segmentsGrob(
        x0 = grid::unit(0.05, "npc"),  # 5% inset from left edge
        x1 = grid::unit(0.95, "npc"),  # 5% inset from right edge
        y0 = grid::unit(0, "npc"),
        y1 = grid::unit(0, "npc"),
        gp = grid::gpar(col = "black", lwd = 0.5)
      ),
      t = bottom_row,
      b = bottom_row,
      l = mid_cols[1],
      r = mid_cols[2]
    )

  }

  # Number of rows in res
  n_rows <- nrow(res)

  # Base + per-row scaling
  total_height <- 0.1 + (n_rows)

  scale_factor <- 1  # adjust as needed
  total_height <- total_height * scale_factor

  # Attach as attribute to your grob
  attr(g, "dynamic_height") <- total_height

  return(invisible(g))

}

if (!exists("use_wrapper")) use_wrapper <- TRUE

if(use_wrapper == TRUE)
{

.Forest_Plot_original <- Forest_Plot

Forest_Plot <- function(Data, ..., session = NULL) {

  user_args <- c(list(Data = Data), list(...))
  orig <- .Forest_Plot_original

  fmls <- formals(orig)
  defaults <- lapply(fmls, function(x) eval(x, envir = parent.frame()))
  defaults[["..."]] <- NULL

  if (!("session" %in% names(fmls))) {

    session <- NULL

  } else {

    defaults$session <- session

  }

  merged <- utils::modifyList(defaults, user_args, keep.null = TRUE)

  call_expr <- match.call(expand.dots = FALSE)$Data
  eval_data <- user_args$Data

  if (is.character(eval_data)) {

    # Data supplied as character vector (file paths)
    Orig_Names <- unname(eval_data)  #

    # If Names not provided, derive: prefer names(), else basename sans ext
    if (is.null(merged$Names)) {

      nm <- names(eval_data)
      if (is.null(nm)) {

        nm <- tools::file_path_sans_ext(basename(eval_data))

      }

      merged$Names <- nm

    }

  } else {

    # Data supplied as symbols or c(symbol, "path", ...)
    parts <- if (is.call(call_expr) && identical(call_expr[[1L]], quote(c))) {

      as.list(call_expr)[-1L]

    } else {

      list(call_expr)

    }

    Orig_Names <- vapply(parts, function(p) {

      if (is.character(p)) p else deparse1(p)

    }, character(1))

    Orig_Names <- gsub('^"|"$', '', Orig_Names)
    if (is.null(merged$Names)) merged$Names <- Orig_Names

  }

  assign(".NF_Orig_Names", Orig_Names, envir = .GlobalEnv)

  ds_labels <- merged$Names

  if (is.null(ds_labels) || all(is.na(ds_labels))) {

    ds_labels <- tools::file_path_sans_ext(basename(Orig_Names))

  }

  ds_labels <- gsub('^"|"$', '', ds_labels)

  message(sprintf(

    "Processing dataset%s: %s",
    if (length(ds_labels) > 1) "s" else "",
    paste(ds_labels, collapse = ", ")

  ))

  on.exit({

    if (exists(".NF_Orig_Names", envir = .GlobalEnv, inherits = FALSE)) {

      rm(".NF_Orig_Names", envir = .GlobalEnv)

    }

  }, add = TRUE)

  default_verbose <- isTRUE(merged$Verbose)

  verbose_mode <- if ("Verbose" %in% names(user_args)) isTRUE(user_args$Verbose) else default_verbose

  merged$Verbose <- verbose_mode

  if (verbose_mode) {

    return(do.call(orig, merged))

  } else {

    merged$Verbose <- FALSE
    return(withCallingHandlers(
      run_with_counter(
        func    = orig,
        args    = merged,
        session = session
      ),
      message = function(m) invokeRestart("muffleMessage"),
      warning = function(w) invokeRestart("muffleWarning")
    ))

  }

}


}
