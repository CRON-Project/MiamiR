#' Title Create forest plots of test statistics of regression model covariates and GWAS summary statistics for SNPs of interest
#'
#' @param Data_Sets These are a list of GWAS summary statistics or munged model test statistics to be plotted; defaults to c("ModelSum", "ModelSum")
#' @param Names These are a list of names for the data sets to be labelled with in the forest plot - this also specifies the order down the page; defaults to c("Vals 1", "Vals 2")
#' @param Data_Set_Colours These are a list of colours for the data sets to be coloured with in the forest plot - this also specifies the order down the page; defaults to c("blue", "darkgreen")
#' @param Chromosome_Columns These are a list of manual chromosome column names for the data sets used in the order specified; defaults to c("chromosome", "chromosome")
#' @param Model_Reference Do you want to use statistics from an R regression model (T/F); defaults to T
#' @param Line_Space A continuous value specifying the boundary between the plot margin and the left and right lines; defaults to 1
#' @param Test_Statistic The type of values being plotted (OR or BETA); defaults to "OR"
#' @param Display_Test_Stat_Se_Column Do you want to display the BETA (SE) raw values on the right hand side (T/F); defaults to F
#' @param Display_Test_Stat_CI_Column Do you want to display the OR (CI - LL to UL) raw values on the right hand side (T/F); defaults to T
#' @param Display_P_Value_Column Do you want to display the P raw values on the right hand side before the test stat column (T/F); defaults to T
#' @param Shapes These are a list of shapes corresponding (in the order of) to the data sets to be plotted; defaults to c("square", "diamond")
#' @param Null_Line_Colour Colour of null line; defaults to "red"
#' @param Null_Line_Type Type of null line; defaults to "dashed"
#' @param X_Axis_Title X axis title; defaults to "BETA"
#' @param X_Axis_Title_Size X axis title; defaults to 20
#' @param X_Axis_Label Do you want an X Axis Label (T/F); defaults to F
#' @param X_Axis_Separation Continuous value of the specific break separation to show on the X axis; defaults to 0.02
#' @param Strip_Colour Colour of the strips highlighting the groups of covariates or SNPs under a certain study group; defaults to "lightblue"
#' @param Strips Do you want above strips to be displayed (T/F); defaults to T
#' @param Pre_Calculated_CIs Have you pre-calculated UL and LL for CIs (T/F); defaults to F
#' @param Legend_Title_Size Legend title size; defaults to 15
#' @param Legend_Text_Size Legend text size; defaults to 12.5
#' @param Legend_Title The name of the title to display next to the legend; defaults to "Study"
#' @param Left_Title The name of the title to display on the left hand side SNPs or covariates; defaults to "SNP"
#' @param P_Value_Title The name of the title to display above P values on right hand side; defaults to "p-value"
#' @param Test_Stat_Se_Title The name of the title to display above raw test statistics values on right hand side; defaults to "BETA (SE)"
#' @param OR_Columns These are a list of manual OR column names for the data sets used in the order specified; defaults to c("OR", "OR")
#' @param Position_Columns These are a list of manual position column names for the data sets used in the order specified; defaults to c("base_pair_location", "base_pair_location")
#' @param SNP_ID_Columns These are a list of manual SNP ID column names for the data sets used in the order specified; defaults to c("variant_id","variant_id")
#' @param Beta_Columns These are a list of manual BETA column names for the data sets used in the order specified; defaults to  c("beta", "beta")
#' @param Standard_Error_Columns These are a list of manual SE column names for the data sets used in the order specified; defaults to c("SE","SE")
#' @param PValue_Columns These are a list of manual P column names for the data sets used in the order specified; defaults to c("P","P")
#' @param Match_Allele_Direction Do you want to match the allele directions to avoid flips (T/F); defaults to T
#' @param Match_Allele_Study If the above is T, which data sets should be the reference point for effect allele direction; defaults to "LbDementia_Sum_Stats"; THIS SHOULD BE THE NAME YOU ASSIGN
#' @param Selected_SNPs A list of SNPs to be shown in the forest plot, automatically selected for from the summary statistics data set; defaults to c("rs2616526", "rs7974838", "rs59867714", "rs9571588", "rs79007041")
#' @param Selected_Covariates A list of covariates to be shown in the forest plot, automatically selected for from the munged model data set; defaults to c()
#' @param Reference_Alleles These are a list of manual reference allele column names for the data sets used in the order specified; defaults to c("other_allele","other_allele")
#' @param Effect_Alleles These are a list of manual effect allele column names for the data sets used in the order specified; defaults to c("effect_allele","effect_allele")
#' @param Upper_CI_Columns These are a list of manual CI UL column names for the data sets used in the order specified; defaults to c("UL","UL")
#' @param Lower_CI_Columns These are a list of manual CI LL column names for the data sets used in the order specified; defaults to c("LL","LL")
#' @param File_Name File name to save plot as; defaults to "Forest_Plot"
#' @param Width Width of saved plot; defaults to 10
#' @param Height Height of saved plot; defaults to 6
#' @param Quality Quality of saved plot (dpi); defaults to 900
#' @param File_Type File type of saved plot; defaults to "jpg"
#' @param X_Axis_Text_Resolution Number of decimal places to display on X axis text; defaults to 1
#' @param Legend_On Do you want to display the legend - TRUE/FALSE; defaults to TRUE
#' @param X_Axis_Text_Size Size of the X axis text labels; defaults to 15
#' @param Null_Buffer Units of space to avoid X axis text labels around the nul point; defaults to 0.1
#'
#' @return Image of Single Forest Plot is saved to the current directory and ggplot object is saved
#' @export
#'
#' @examples
#'
#'
#'
#'Forest_Plot_SNPs_BETA_Peak_Finder <- Forest_Plot(Data_Sets = c("Intelligence_Peaks"),
#'                                                 X_Axis_Separation = 0.05,
#'                                                 File_Name = "Forest_Plot", Width =10, Height = 9, Quality = 900,
#'                                                 File_Type = "jpg")
#'
#'
#'
#'
#'


Forest_Plot <- function(Data_Sets = c(),
                        Names = NULL,
                        Data_Set_Colours = viridis::viridis(length(Data_Sets)),
                        Chromosome_Columns = c(),
                        Left_Spaces = 2, #def 2
                        Right_Spaces = 2,
                        Missings = F,
                        Double_Label = F,
                        P_Stat_Spaces = 3,
                        Model_Reference = FALSE,
                        X_Axis_Text_Size = 125,
                        Line_Space = 0.01,
                        Test_Statistic = NULL,
                        Display_Test_Stat_Se_Column = FALSE,
                        Display_Test_Stat_CI_Column = FALSE,
                        Display_P_Value_Column = TRUE,
                        Shapes = NULL,
                        Null_Line_Colour = "red",
                        Null_Line_Type = "dashed",
                        X_Axis_Title = NULL,
                        Null_Buffer = 0.0,
                        X_Axis_Title_Size = 125,
                        SNP_Stat_Text_Size = 125,  #some non odd tens dont work?
                        X_Axis_Label = TRUE,
                        X_Axis_Separation = NULL,
                        Strip_Colour = "goldenrod",
                        Strips = TRUE,
                        X_Axis_Text_Resolution = 2,
                        P_Value_Resolution = 1,
                        Pre_Calculated_CIs = FALSE,
                        Legend_On = FALSE,
                        Legend_Title_Size = 100,
                        Legend_Text_Size = 100,
                        Legend_Title = "STUDY",
                        Left_Title = "SNP",
                        P_Value_Title = "p-value",
                        Test_Stat_Se_Title = "BETA (SE)",
                        OR_Columns = c(),
                        Position_Columns = c() , SNP_ID_Columns = c(),
                        Beta_Columns = c(), Standard_Error_Columns = c(),
                        PValue_Columns = c(),
                        Match_Allele_Direction = TRUE,
                        Match_Allele_Study = NULL,
                        Selected_SNPs = c(),
                        Selected_Covariates = c(),
                        Reference_Alleles = c(),  Effect_Alleles = c(),
                        Upper_CI_Columns = c(), Lower_CI_Columns = c(),
                        File_Name = "Forest_Plot", Width =10, Height = 10, Quality = 600,
                        File_Type = "jpg"
                        )

  {


  if (is.null(Data_Sets) || length(Data_Sets) == 0) {
    stop("Error: 'Data_Sets' is NULL or empty. Please provide valid data sets.")
  }

  bases <- basename(Data_Sets)
  num_pattern_count <- sum(grepl("^[0-9]+_", bases)) #changes later, only this for now

#  print(Data_Sets)
#  print(num_pattern_count)


  if (num_pattern_count > 2) {

  # Extract common suffix, removing file paths and leading numbers
  dataset_suffixes <- gsub(".*[/\\\\]\\d+_", "", Data_Sets)

  print(dataset_suffixes)
 # print(Data_Sets)

  # Function to find common substrings of a minimum length
  # Extract unique phrases between underscores and group dataset names
  # Function to find unique phrases between underscores and merge similar groups
  # Function to find unique phrases between underscores and merge similar groups
  find_and_merge_strict_groups <- function(strings) {
    groups <- list()

    for (s in strings) {
      # Extract all phrases between underscores, including start and end boundaries
      phrases <- unlist(regmatches(s, gregexpr("(?<=_)([^_]+)(?=_)", s, perl = TRUE)))

      # Add full filename as a phrase to capture whole patterns
      phrases <- c(phrases, s)

      for (phrase in phrases) {
        if (!phrase %in% names(groups)) {
          groups[[phrase]] <- c(s)
        } else {
          groups[[phrase]] <- unique(c(groups[[phrase]], s))
        }
      }
    }

    # Merge groups with overlapping dataset entries but only if they match strictly
    merged_groups <- list()
    for (name in names(groups)) {
      merged <- FALSE
      for (existing_name in names(merged_groups)) {
        if (any(groups[[name]] %in% merged_groups[[existing_name]]) && grepl(paste0("\\b", existing_name, "\\b"), name)) {
          merged_groups[[existing_name]] <- unique(c(merged_groups[[existing_name]], groups[[name]]))
          merged <- TRUE
          break
        }
      }
      if (!merged) {
        merged_groups[[name]] <- groups[[name]]
      }
    }

    return(merged_groups)
  }

  # Find and merge groups with strict matching
  strict_merged_groups <- find_and_merge_strict_groups(dataset_suffixes)

  # Display the strict merged groups
  strict_merged_groups

#  print(strict_merged_groups)

 # z

  # Get unique suffixes and assign colors
  unique_suffixes <- unique(dataset_suffixes)
  colour_palette <- viridis::viridis(length(unique_suffixes))
  suffix_to_color <- setNames(colour_palette, unique_suffixes)

  # Assign colors based on suffix
  Data_Set_Colours <- suffix_to_color[dataset_suffixes]

  # Print to check assignments
#  print(data.frame(Data_Sets, dataset_suffixes, Data_Set_Colours))

  print("Colouring")


}
  # Check if Data_Sets contains file paths
  if (all(file.exists(Data_Sets))) {
    message("Loading datasets from file paths...")

    dataset_names <- c()  # Initialize empty vector to store dataset names

    for (path in Data_Sets) {
      # Extract filename without extension
      dataset_name <- tools::file_path_sans_ext(basename(path))

      message("Processing file: ", path)
      message("Extracted dataset name: ", dataset_name)

      # Read the data
      df <- if (grepl("\\.csv$", path, ignore.case = TRUE)) {
        read.csv(path, stringsAsFactors = FALSE)
      } else if (grepl("\\.rds$", path, ignore.case = TRUE)) {
        readRDS(path)
      } else if (grepl("\\.(txt|tab|tsv)$", path, ignore.case = TRUE)) {
     #   read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,  colClasses = "character")
        vroom::vroom(path, col_types = vroom::cols(
          ALLELE0 = vroom::col_character(),
          ALLELE1 = vroom::col_character()
        ))
        } else {
        stop("Unsupported file format. Supported formats: CSV, RDS, TXT, TAB, TSV.")
      }

      # Store dataset in environment
      assign(dataset_name, df, envir = .GlobalEnv)
      message(paste("Dataset", dataset_name, "loaded into environment."))

      # Append dataset name to vector
      dataset_names <- c(dataset_names, dataset_name)
    }

    # Update Data_Sets with dataset names only
    Data_Sets <- dataset_names

  } else {
    message("Using datasets from the R environment...")
    Data_Sets <- lapply(Data_Sets, function(name) {
      if (!is.character(name)) stop("Dataset name must be a character string.")

      if (exists(name, envir = .GlobalEnv)) {
        return(name)  # Keep dataset name in Data_Sets
      } else {
        stop(paste("Dataset", name, "not found in environment."))
      }
    })

    # Convert list to character vector
    Data_Sets <- unlist(Data_Sets)
  }

  message("Final Data_Sets list: ", paste(Data_Sets, collapse = ", "))

  # Return Data_Sets with correct names
 # return(Data_Sets)






print(Data_Sets)

print(Names)

if(!is.null(Names))
{
  print("Using Names Provided")
}
else{
if (!is.null(Data_Sets)) {
    print("Using Data Set Names UPDATED")




#    Data_Sets <- sub("^[0-9]+_", "", Data_Sets)


  #  print(Data_Sets)


    Names <- (Data_Sets)  # Use names of the Data_Sets vector





  } else {
    print("Auto-naming")
    Names <- paste("Dataset", seq_along(Data_Sets))  # Fallback to generic names
  }

}
  print(Names)


#  print("fail here")

 # Left_Spaces <- strrep(" ", Left_Spaces) - #due to html white out later

 # print(Left_Spaces)



 # Left_Spaces <- strrep("Z", 10)




  #Goes below above as need to match based on STUDY assigned from Name
  if (is.null(Match_Allele_Study)) {

    Match_Allele_Study <- Names[1]
    print(paste0("Matching allele directions automatically to ", Match_Allele_Study) )

    print(Match_Allele_Study)

#    if (num_pattern_count > 2) {
#
 #     Match_Allele_Study <- sub("^[0-9]+_", "", Match_Allele_Study)
  #    }

   # print(Match_Allele_Study)

  }

  if (is.null(Shapes)) {
    Shapes <- rep("square", length(Data_Sets))  # Assign "square" to all
  }

  if (is.null((Test_Statistic))) {

    print("No Test Stat Allocated")

  # Randomly get one dataset from Data_Sets before the loop
  dataset_name <- sample(Data_Sets, 1)

#  print(dataset_name)

#  print("Here")

  Data <- get(dataset_name)



  # Initialize variables
  Beta_Column <- NULL
  OR_Column <- NULL
  Test_Statistic <- NULL

  # Check for BETA columns first
  allowed_names_betas <- c("BETA", "Beta", "beta", "B", "stdBeta")
  for (allowed_name in allowed_names_betas) {
    if (allowed_name %in% colnames(Data)) {
      print(paste0("Detected BETA values: Using ", allowed_name, " as the Test_Statistic column option (automatic)"))
      Beta_Column <- allowed_name
      Test_Statistic <- "BETA"  # Assign Test_Statistic as BETA
      break
    }
  }

  # If no beta column is found, check for OR columns
  if (is.null(Beta_Column)) {
    allowed_names_ors <- c("OR", "or", "Or", "Odds", "Odd Ratio", "odds_ratio")
    for (allowed_name in allowed_names_ors) {
      if (allowed_name %in% colnames(Data)) {
        print(paste0("Detected OR values: Using ", allowed_name, " as the Test_Statistic column option (automatic)"))
        OR_Column <- allowed_name
        Test_Statistic <- "OR"  # Assign Test_Statistic as OR
        break
      }
    }
  }



  }


  if (is.null((X_Axis_Title))) {

    X_Axis_Title <- Test_Statistic

 #   print(X_Axis_Title)

  }





  # Initialize variables back to nothing
  Beta_Column <- NULL
  OR_Column <- NULL
#  Test_Statistic <- NULL

  Combined_Processed_Data <- data.frame()

 # print(Data_Sets)

  for (i in seq_along(Data_Sets)) {

 #   print(Data_Sets)

    # dataset_name <- Data_Sets[i]
    # Data <- dataset_name

    dataset_name <- Data_Sets[i]  # Get the dataset name

#    print(dataset_name)

    Data <- get(dataset_name)  # Load the dataset using get()

#    print(Data)



    corresponding_name <- Names[i]  # Get the corresponding name
    corresponding_shape <- Shapes[i]  # Get the corresponding shape


    if(Model_Reference == F)
    {
    Chromosome_Column <- Chromosome_Columns[i]
    Position_Column <- Position_Columns[i]
    SNP_ID_Column <- SNP_ID_Columns[i]
    }
    PValue_Column <- PValue_Columns[i]
    Standard_Error_Column <- Standard_Error_Columns[i]

 #   print(Test_Statistic)

    if(Test_Statistic == "BETA")
    {
    Beta_Column <- Beta_Columns[i]
    }else{
    OR_Column <- OR_Columns[i]
    }

    if(Model_Reference == F)
    {
    Reference_Allele <- Reference_Alleles[i]
    Effect_Allele <- Effect_Alleles[i]
    Upper_CI_Column <- Upper_CI_Columns[i]
    Lower_CI_Column <- Lower_CI_Columns[i]
    }

    # Print dataset name and corresponding name
    print(paste("Processing dataset:", dataset_name))
    print(paste("Corresponding name:", corresponding_name))


    if(Model_Reference == F)
    {

    if(is.null(Reference_Allele))
    {

    allowed_names_reference_allele <- c("A0", "A2", "REF", "other_allele", "a0", "a2", "ref", "Non_effect_Allele", "Reference", "reference", "allele0", "allele2", "ALLELE0", "ALLELE2")

    for (allowed_name in allowed_names_reference_allele) {
      if (allowed_name %in% colnames(Data)) {
        usable_ref_allele <- colnames(Data)[which(colnames(Data) == allowed_name)]
        break
      }
    }



    print(paste0("Using", " ", allowed_name, " ", "as reference allele (automatic)"))
    Reference_Allele <- allowed_name

    }else{
      print(paste0("Using", " ", Reference_Allele, " ", "as reference allele (manual)"))
    }


    if(is.null(Effect_Allele))
    {

    allowed_names_alternate_allele <- c("A1", "ALT", "a1", "alt", "Alternate", "Effect_Allele", "effect_allele", "alternate", "allele1", "ALLELE1")

    for (allowed_name in allowed_names_alternate_allele) {
      if (allowed_name %in% colnames(Data)) {
        usable_effect_allele <- colnames(Data)[which(colnames(Data) == allowed_name)]
        break
      }
    }



    print(paste0("Using", " ", allowed_name, " ", "as effect allele (automatic)"))
    Effect_Allele <- allowed_name


    }else{
      print(paste0("Using", " ", Effect_Allele, " ", "as effect allele (manual)"))
    }




    if(is.null(Chromosome_Column))
    {

  allowed_names_chromosomes <- c("chromosome", "chrom", "chr", "CHROM", "Chromosome", "CHR", "Chr", "Chrom")

  for (allowed_name in allowed_names_chromosomes) {
    if (allowed_name %in% colnames(Data)) {
      usable_chrom_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as chromosome column (automatic)"))
  Chromosome_Column <- allowed_name


  }else{
    print(paste0("Using", " ", Chromosome_Column, " ", "as chromosome column (manual)"))
  }


    }

if(Test_Statistic == "BETA")
{

    if(is.null(Beta_Column))
    {


  allowed_names_betas <- c("BETA", "Beta", "beta", "B", "stdBeta")



  for (allowed_name in allowed_names_betas) {
    if (allowed_name %in% colnames(Data)) {
      usable_chrom_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as beta column (automatic)"))
  Beta_Column <- allowed_name



  }else{
    print(paste0("Using", " ", Beta_Column, " ", "as beta column (manual)"))
  }

}

    if(Test_Statistic == "OR")
    {

      if(is.null(OR_Column))
      {


        allowed_names_ors <- c("OR", "or", "Or", "Odds", "Odd Ratio", "odds_ratio") #maybe add more



        for (allowed_name in allowed_names_ors) {
          if (allowed_name %in% colnames(Data)) {
            usable_chrom_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
            break
          }
        }


        print(paste0("Using", " ", allowed_name, " ", "as OR column (automatic)"))
        OR_Column <- allowed_name



      }else{
        print(paste0("Using", " ", OR_Column, " ", "as OR column (manual)"))
      }

    }

    if(is.null(Standard_Error_Column))
    {


  allowed_names_ses<- c("SE", "se", "Se", "Standard_Error", "standard_error", "Standard_Error_of_Beta")

  for (allowed_name in allowed_names_ses) {
    if (allowed_name %in% colnames(Data)) {
      usable_chrom_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as standard error column (automatic)"))
  Standard_Error_Column <- allowed_name



    }else{
      print(paste0("Using", " ", Standard_Error_Column, " ", "as standard error column (manual)"))
    }



    if(is.null(PValue_Column))
    {
  allowed_names_pvals <- c("P", "p", "Pvalue", "pvalue", "P-Value", "p-value", "p-Value",
                           "P-VALUE", "logp","LogP", "LOGP", "Logp", "log10p","Log10P",
                           "LOG10P", "Log10p",
                           "log10p", "LOG10P", "-LOG10P", "", "p_value")

  for (allowed_name in allowed_names_pvals) {
    if (allowed_name %in% colnames(Data)) {
      usable_p_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as P-Val column (automatic)"))
  PValue_Column <- allowed_name

    }else{
      print(paste0("Using", " ", PValue_Column, " ", "as P-Val column (manual)"))
    }



  if (!("P" %in% colnames(Data)) & any(colnames(Data) %in% c("logp", "LogP", "LOGP", "Logp",
                                                             "log10p", "Log10P", "LOG10P",
                                                             "Log10p", "-LOG10P"))) {

    print("No P Value in first dataset, Calculating from LOG10P column detected")
    Data$P <- 10^-(as.numeric(Data[[PValue_Column]]))
    PValue_Column <- "P"

  }


    if(Model_Reference == F)
    {

    if(is.null(Position_Column))
    {

  allowed_names_pos <- c("POS", "pos", "Pos", "BP", "BPos", "bpos", "BPOS", "bPos",
                         "Position", "position", "POSITION", "genpos", "GENPOS",
                         "Genpos", "base_pair_location")


  for (allowed_name in allowed_names_pos) {
    if (allowed_name %in% colnames(Data)) {
      usable_pos_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as position column (automatic)"))
  Position_Column <- allowed_name

    }else{
      print(paste0("Using", " ", Position_Column, " ", "as position column (manual)"))
    }


    if(is.null(SNP_ID_Column))
    {


  allowed_names_SNP <- c("ID", "Id", "ID", "RsID", "RsId","RSID", "snp", "SNP", "Snp",
                         "snv" ,"SNV" , "Snv", "RS", "rs", "variant_id" )

  for (allowed_name in allowed_names_SNP) {
    if (allowed_name %in% colnames(Data)) {
      usable_snp_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as SNP column (automatic)"))
  SNP_ID_Column <- allowed_name

    }else{
      print(paste0("Using", " ", SNP_ID_Column, " ", "as SNP column (manual)"))
    }


    if(Pre_Calculated_CIs == T)
    {

    if(is.null(Upper_CI_Column))
    {

      allowed_names_upper_ci <- c("UL", "Upper_CI", "CIU")

      for (allowed_name in allowed_names_upper_ci) {
        if (allowed_name %in% colnames(Data)) {
          usable_chrom_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
          break
        }
      }


      print(paste0("Using", " ", allowed_name, " ", "as upper CI column (automatic)"))
      Upper_CI_Column <- allowed_name


    }else{
      print(paste0("Using", " ", Upper_CI_Column, " ", "as upper CI column (manual)"))
    }


    if(is.null(Lower_CI_Column))
    {

      allowed_names_Lower_ci <- c("UL", "Lower_CI", "CIU")

      for (allowed_name in allowed_names_Lower_ci) {
        if (allowed_name %in% colnames(Data)) {
          usable_chrom_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
          break
        }
      }


      print(paste0("Using", " ", allowed_name, " ", "as lower CI column (automatic)"))
      Lower_CI_Column <- allowed_name


    }else{
      print(paste0("Using", " ", Lower_CI_Column, " ", "as lower CI column (manual)"))
    }


    }else{
      print("No Upper of Lower CIs specified - these will be manually calculated")
    }


  #Manually assign columns for ease of use
  Data$CHROM <- Data[[Chromosome_Column]]
  Data$GENPOS <- Data[[Position_Column]]
  Data$ID <- Data[[SNP_ID_Column]]



    }
  if(Model_Reference == T)
  {
    Data$ID <- Data$Covariate
  }
  Data$P <- Data[[PValue_Column]]
  Data$SE <- Data[[Standard_Error_Column]]
  if(Pre_Calculated_CIs == T)
  {
  Data$LL <- Data[[Lower_CI_Column]]
  Data$UL <- Data[[Upper_CI_Column]]
  }
  if(Test_Statistic == "BETA")
  {
  Data$BETA <- Data[[Beta_Column]]
  }
  if(Test_Statistic == "OR")
  {

  Data$BETA <- Data[[OR_Column]] # keep as BETA in code to save time

  }
  Data$STUDY <- corresponding_name
  Data$Shape <- corresponding_shape
  if(Model_Reference == F)
  {
  Data$ALLELE0 <- Data[[Reference_Allele]]
  Data$ALLELE1 <- Data[[Effect_Allele]]

}

  if(Model_Reference == F)
  {


  print("Filtering for key SNPs")
  print(Selected_SNPs)


#  print(Data)



  }
  if(Model_Reference == T)
  {
    if(!is.null(Selected_Covariates)){
    print("Filtering for key Covariates (manual)")
    Selected_SNPs <- Selected_Covariates
    print(Selected_SNPs)
    Data <- Data[Data$group %in% Selected_SNPs,]
    }else{
    print("No covariates specified - retaining all (automatic)")
    print("Automatically generated order:")
    Selected_Covariates <- Data$group[!duplicated(Data$group)]
    print(Selected_Covariates)
    Selected_SNPs <- Selected_Covariates
    }
  }



#sill fine if for OR as renaming BETA

 # print(Data)

  if(Model_Reference == F)
  {
  Data <- Data %>% dplyr::select(ID, ALLELE0, ALLELE1, CHROM, GENPOS, BETA, SE, P, STUDY, Shape)
  }

  if(Model_Reference == T)
  {
    Data <- Data %>% dplyr::select(ID,  BETA, SE, P, STUDY, Shape, group, Reference)

  }



  if(Model_Reference == F)
  {

  #keep maybe an RS also provided
  Data$Backup_ID <- Data$ID
  #need to remake earlier in case RS also provided
  Data$ID <- stringi::stri_c("chr", Data$CHROM, ":", Data$GENPOS, ":", Data$ALLELE0, ":", Data$ALLELE1)

  Data$COORD_Norm <- stringi::stri_c("chr", Data$CHROM, ":", Data$GENPOS, ":", Data$ALLELE0, ":", Data$ALLELE1)
  Data$COORD_Alt <- stringi::stri_c("chr", Data$CHROM, ":", Data$GENPOS, ":", Data$ALLELE1, ":", Data$ALLELE0)

  }


  if(Model_Reference == F)
  {


    #if you want all of them, like in peak finder plot
  if (length(Selected_SNPs) == 0) {
    Data <- Data
  } else {
  Data <-Data[Data$COORD_Norm %in% Selected_SNPs | Data$COORD_Alt %in% Selected_SNPs |  Data$ID %in% Selected_SNPs | Data$Backup_ID %in% Selected_SNPs , ] # will select even if there is a flip
  }


  }

  Combined_Processed_Data <- rbind(Combined_Processed_Data, Data)



  }





  print(Combined_Processed_Data)



#bef <-  print(nrow(Combined_Processed_Data)) # if missings removed mod end point - that

if(Missings == F)

{

  Combined_Processed_Data <- Combined_Processed_Data[!(Combined_Processed_Data$BETA == 0 & Combined_Processed_Data$SE == 0 & Combined_Processed_Data$P == 1), ]

}

#aft <-   print(nrow(Combined_Processed_Data))




 #  print(Combined_Processed_Data)





if(Model_Reference == T)
{

  # Get unique references and their corresponding group values where Reference is not "None"
  unique_references_with_group <- Combined_Processed_Data %>%
    dplyr::filter(Reference != "None") %>%
    dplyr::distinct(Reference, group)

  # Create blank rows for each unique Reference value with NA for numeric columns,
  # set the ID to the Reference value, and the group to the corresponding group value

  if(Test_Statistic == "OR")
  {

  blank_rows <- data.frame(ID = unique_references_with_group$Reference,
                           BETA = 1,
                           SE = 0,
                           P = 1,
                           STUDY = "",
                           Shape = "square",
                           group = unique_references_with_group$group,
                           Reference = unique_references_with_group$Reference,
                           stringsAsFactors = FALSE)

  }else{


    blank_rows <- data.frame(ID = unique_references_with_group$Reference,
                             BETA = 0,
                             SE = 0,
                             P = 1,
                             STUDY = "",
                             Shape = "square",
                             group = unique_references_with_group$group,
                             Reference = unique_references_with_group$Reference,
                             stringsAsFactors = FALSE)


  }


  # Combine the original dataframe with the new blank rows
  Combined_Processed_Data <- dplyr::bind_rows(Combined_Processed_Data, blank_rows)



  Combined_Processed_Data$Left_Plot_Value <- ifelse(
    Combined_Processed_Data$ID == Combined_Processed_Data$Reference,
    paste0(Combined_Processed_Data$group, " (Ref: ", Combined_Processed_Data$Reference, ")"),
    Combined_Processed_Data$ID
  )



}else{





  Combined_Processed_Data$Left_Plot_Value <- Combined_Processed_Data$ID




}

#  print(Combined_Processed_Data)

#Peak finder says now but keep for other circum?

  if(Missings == T)

  {

  # Extract all unique IDs and STUDY values
  unique_ids <- unique(Combined_Processed_Data$ID)
  unique_studies <- unique(Combined_Processed_Data$STUDY)

  # Create a full combination of all unique IDs and STUDY values
  full_combination <- expand.grid(ID = unique_ids, STUDY = unique_studies)

  # Identify missing combinations (rows not in the original dataset)
  missing_rows <- full_combination %>%
    dplyr::anti_join(Combined_Processed_Data, by = c("ID", "STUDY"))

  # Add the missing rows with BETA, SE, and P set to 0
  missing_rows <- missing_rows %>%
    dplyr::mutate(
      BETA = 0,
      SE = 0,
      P = 1
    )


#  print(missing_rows)

  missing_rows <- missing_rows %>%
    dplyr::mutate(
      ALLELE0 = sub("^[^:]+:[^:]+:([^:]+):.*$", "\\1", ID),
      ALLELE1 = sub("^[^:]+:[^:]+:[^:]+:([^:]+)$", "\\1", ID), # Part after 3rd colon
      CHROM = stringr::str_extract(ID, "(?<=chr)[^:]+"),                # Extract part after 'chr' and before ':'
      GENPOS = stringr::str_extract(ID, "(?<=:)[^:]+(?=:)"),
      Left_Plot_Value = ID,
      Shape = "cross"

      # Extract part after 1st colon and before 2nd
    )


  #print(missing_rows)

#zzz
  missing_rows$CHROM <- as.numeric(missing_rows$CHROM)

  missing_rows$GENPOS <- as.numeric(missing_rows$GENPOS)

#  print(missing_rows)

#  print("fails here")

#  print(nrow(missing_rows))

  # if(nrow(missing_rows) >= 1) #sometimes no missing
  # {
  # Combine the original dataset with the missing rows
  Combined_Processed_Data <-  dplyr::bind_rows(Combined_Processed_Data, missing_rows)

  # }
  # View the updated dataframe
 # print(Combined_Processed_Data)


 #  zzz





  Combined_Processed_Data <- dplyr::mutate(Combined_Processed_Data,
                                           Shape = dplyr::case_when(
                                             BETA == 0 & SE == 0 & P == 1 ~ "cross",
                                             TRUE ~ Shape  # Keeps existing values
                                           )
  )

 # print(Combined_Processed_Data)


}


#   return(Combined_Processed_Data)

  if(Model_Reference == F)
  {

if(Match_Allele_Direction == T)

{

  print(paste0("Processing of matching allele directions to: ", Match_Allele_Study))




  if(Test_Statistic == "BETA")
  {
#   print("Matching study effect allele directions...")
#
#   #  print(Combined_Processed_Data$STUDY)
#    # print(Match_Allele_Study)
#
#
#   # Convert ALLELE0 and ALLELE1 to uppercase to ensure uniformity
#   Combined_Processed_Data$ALLELE0 <- toupper(Combined_Processed_Data$ALLELE0)
#   Combined_Processed_Data$ALLELE1 <- toupper(Combined_Processed_Data$ALLELE1)
#
#
# #  print("Data before BETA flips")
# #  print(Combined_Processed_Data)
#
#
#   # Loop through the data frame and check for matching IDs with swapped alleles
#   for (i in 1:nrow(Combined_Processed_Data)) {
#
#  #   print(Combined_Processed_Data$STUDY[i] == Match_Allele_Study)
#     # Ensure the row with the 'Intelligence_Sum_Stats' in STUDY is the master row
#     if (Combined_Processed_Data$STUDY[i] == Match_Allele_Study) {
#
#       for (j in 1:nrow(Combined_Processed_Data)) {
#         # Skip the master row (Intelligence_Sum_Stats)
#         if (i != j && Combined_Processed_Data$STUDY[j] != Match_Allele_Study) {
#
#      #     print("Hi")
#     #      print(i)
#      #     print(j)
#
#
#           print(Combined_Processed_Data$ID[i])
#           print(Combined_Processed_Data$ID[j])
#                 print(Combined_Processed_Data$ALLELE0[i])
#                       print(Combined_Processed_Data$ALLELE1[j])
#                             print(Combined_Processed_Data$ALLELE1[i])
#                                   print(Combined_Processed_Data$ALLELE0[j])
#
#
#
#          # print("HEEW")
#     #      print(Combined_Processed_Data$ID[i])
#
#           # Check if IDs are the same, and alleles are swapped
#           if (Combined_Processed_Data$ID[i] == Combined_Processed_Data$ID[j] &&
#               Combined_Processed_Data$ALLELE0[i] == Combined_Processed_Data$ALLELE1[j] &&
#               Combined_Processed_Data$ALLELE1[i] == Combined_Processed_Data$ALLELE0[j]) {
#
#
#
#
#             # Multiply the BETA value of the swapped allele SNP by -1 (for non-master row)
#             Combined_Processed_Data$BETA[j] <- Combined_Processed_Data$BETA[j] * -1
#           }
#         }
#       }
#     }
#   }
#


#    print(Combined_Processed_Data)

 #   print(Combined_Processed_Data$STUDY)
#    print(Match_Allele_Study)




    if (num_pattern_count > 2) {

    # Remove number prefix from Match_Allele_Study for comparison
    Match_Allele_Study_Clean <- stringr::str_remove(Match_Allele_Study, "^[0-9]+_")

    }else{
      Match_Allele_Study_Clean <- Match_Allele_Study
    }

    # Ensure STUDY column has no leading/trailing spaces
    Combined_Processed_Data$STUDY <- trimws(Combined_Processed_Data$STUDY)


    if (num_pattern_count > 2) {

    # Extract reference alleles, matching STUDY even if the number prefix differs
    reference <- Combined_Processed_Data %>%
      dplyr::mutate(STUDY_Clean = stringr::str_remove(STUDY, "^[0-9]+_")) %>%  # Remove number prefix
      dplyr::filter(STUDY_Clean == Match_Allele_Study_Clean) %>%  # Compare cleaned names
      dplyr::select(ID, ALLELE0, ALLELE1, COORD_Norm, COORD_Alt) %>%
      dplyr::rename(Ref_ALLELE0 = ALLELE0, Ref_ALLELE1 = ALLELE1)

    }else{

      # Extract reference alleles, matching STUDY even if the number prefix differs
      reference <- Combined_Processed_Data %>%
        dplyr::mutate(STUDY_Clean = STUDY) %>%  # Remove number prefix
        dplyr::filter(STUDY_Clean == Match_Allele_Study_Clean) %>%  # Compare cleaned names
        dplyr::select(ID, ALLELE0, ALLELE1, COORD_Norm, COORD_Alt) %>%
        dplyr::rename(Ref_ALLELE0 = ALLELE0, Ref_ALLELE1 = ALLELE1)

    }





#If from peak finder because of COORD_Uni they will also have ref match - think about raw later


    reference <- reference %>%
      dplyr::rename_with(~ paste0("REF_", .), everything())

 # print(reference)


    #print(Combined_Processed_Data)


    print("Aligning reference alleles with data...")


    Combined_Processed_Data2 <- Combined_Processed_Data %>%
      dplyr::left_join(reference, by = c("COORD_Norm" = "REF_COORD_Norm"))


    #print(sum(is.na(Combined_Processed_Data2$REF_Ref_ALLELE1)))
    #print(sum(is.na(Combined_Processed_Data2$REF_Ref_ALLELE0)))


    #return(Combined_Processed_Data2)

    print("Finding alleles that failed to match at first attempt...")

    # Step 2: Identify rows with no match in COORD_Norm
    unmatched <- Combined_Processed_Data2 %>%
      dplyr::filter(is.na(REF_ID)) %>%
      dplyr::select(names(Combined_Processed_Data))  # Retain only original columns for the fallback join not joined ones too.





    print("Trying to align reference to failures...")

    #print(unmatched)


    # Step 3: Perform the second join using COORD_Alt
    fallback <- unmatched %>%
      dplyr::left_join(reference, by = c("COORD_Norm" = "REF_COORD_Alt"))

    #print(fallback)


    print("Removing failures from initial success...and binding second attempt of failures...")

    # Step 4: Combine matched rows from the first and second join
    Combined_Processed_Data2 <- Combined_Processed_Data2 %>%
      dplyr::filter(!is.na(REF_ID)) %>%  # Keep rows matched on COORD_Norm
      dplyr::bind_rows(fallback)  # Add rows matched on COORD_Alt


    Combined_Processed_Data <- Combined_Processed_Data2


    # Combined_Processed_Data <- Combined_Processed_Data  %>%
    #       dplyr::mutate(
    #           # Flip BETA only if alleles are swapped
    #           BETA = dplyr::if_else(!is.na(Ref_ALLELE0) & !is.na(Ref_ALLELE1) &
    #                                   ALLELE0 == Ref_ALLELE1 & ALLELE1 == Ref_ALLELE0,
    #                                 BETA * -1, BETA)
    #         )


    print("Reversing BETA values of alleles where REF and ALT match reference study REF & ALT but are flipped...")

    Combined_Processed_Data <- Combined_Processed_Data %>%
      dplyr::mutate(
        BETA_Flipped = !is.na(REF_Ref_ALLELE0) & !is.na(REF_Ref_ALLELE1) &
          ALLELE0 == REF_Ref_ALLELE1 & ALLELE1 == REF_Ref_ALLELE0,
        BETA = dplyr::if_else(BETA_Flipped, BETA * -1, BETA),
        No_Match = (is.na(REF_COORD_Norm) & is.na(REF_COORD_Alt))
      )


    print("Unified ID dervies from the reference one as alleles matched to this direction")


    # View only the rows where BETA was flipped
    flipped_rows <- Combined_Processed_Data %>%
      dplyr::filter(BETA_Flipped)

    print(paste0("Number of alleles flipped across any study: ", nrow(flipped_rows)))


    #print(table(is.na(Combined_Processed_Data$REF_COORD_Norm) & is.na(Combined_Processed_Data$REF_COORD_Alt)))

    #unmatched_rows <- Combined_Processed_Data %>%
    # dplyr::filter(is.na(Combined_Processed_Data$REF_COORD_Norm) & is.na(Combined_Processed_Data$REF_COORD_Alt))


    #print(unmatched_rows)



    # View only the rows where BETA was flipped
    unmatched_rows_2 <- Combined_Processed_Data %>%
      dplyr::filter(No_Match)

    #print(unmatched_rows_2)

    print(paste0("Number of alleles which were not present either aligned or flipped from any study in the reference study: ", nrow(unmatched_rows_2)))





    #Make a unified ID showing direction of effect for all based on reference picked for allels flips

    print("Creating unified IDs for Flips, based on reference study direction of effect")

    Combined_Processed_Data$COORD_Uni <- stringi::stri_c("chr", Combined_Processed_Data$CHROM, ":", Combined_Processed_Data$GENPOS, ":", Combined_Processed_Data$REF_Ref_ALLELE0, ":", Combined_Processed_Data$REF_Ref_ALLELE1)


#return(Combined_Processed_Data)


# Print final dataset
#print(Combined_Processed_Data)
   # print(Match_Allele_Study)


   }

#  print(Test_Statistic)





 print(Combined_Processed_Data)
# return(Combined_Processed_Data)



  if(Test_Statistic == "OR")
  {
    print("Matching study effect allele directions...")

    # Convert ALLELE0 and ALLELE1 to uppercase to ensure uniformity
    Combined_Processed_Data$ALLELE0 <- toupper(Combined_Processed_Data$ALLELE0)
    Combined_Processed_Data$ALLELE1 <- toupper(Combined_Processed_Data$ALLELE1)

    # Loop through the data frame and check for matching IDs with swapped alleles
    for (i in 1:nrow(Combined_Processed_Data)) {
      # Ensure the row with the 'Intelligence_Sum_Stats' in STUDY is the master row
      if (Combined_Processed_Data$STUDY[i] == Match_Allele_Study) {

        for (j in 1:nrow(Combined_Processed_Data)) {
          # Skip the master row (Intelligence_Sum_Stats)
          if (i != j && Combined_Processed_Data$STUDY[j] != Match_Allele_Study) {

            # Check if IDs are the same, and alleles are swapped
            if (Combined_Processed_Data$ID[i] == Combined_Processed_Data$ID[j] &&
                Combined_Processed_Data$ALLELE0[i] == Combined_Processed_Data$ALLELE1[j] &&
                Combined_Processed_Data$ALLELE1[i] == Combined_Processed_Data$ALLELE0[j]) {

              # Multiply the BETA value of the swapped allele SNP by -1 (for non-master row)
              Combined_Processed_Data$BETA[j] <- 1 / Combined_Processed_Data$BETA[j]
            }
          }
        }
      }
    }

  }

}

  }else{
    print("Model References don't require allele adjustment!")
  }

   #Want to make sure I am actually plotting what has been adjusted for ! direction wise.



if(Model_Reference == F)
{

 # Combined_Processed_Data$Backup_ID <-    Combined_Processed_Data$ID
 # print( Combined_Processed_Data$Backup_ID)
  Combined_Processed_Data$ID <- Combined_Processed_Data$COORD_Uni

#print(Combined_Processed_Data)


  Combined_Processed_Data$RS <- Combined_Processed_Data$COORD_Uni


  Combined_Processed_Data$Left_Plot_Value <- Combined_Processed_Data$COORD_Uni

#  Combined_Processed_Data$Left_Plot_Value <- paste0( Combined_Processed_Data$Left_Plot_Value, "<br>",
 #                                                    Combined_Processed_Data$Backup_ID)


  if(Double_Label == T)
  {

  # Combined_Processed_Data$Left_Plot_Value <- paste0(
  #   Combined_Processed_Data$Left_Plot_Value, "<br>", "(",
  #   Combined_Processed_Data$Backup_ID, ")"
  # )
    Combined_Processed_Data$Left_Plot_Value <- paste0(
      Combined_Processed_Data$Left_Plot_Value,
      "<br><br>(",  # Adds an extra line break to create spacing
      Combined_Processed_Data$Backup_ID, ")"
    )



  }

  print(Combined_Processed_Data)

}

  res <- Combined_Processed_Data


#  print(res)



  #as long as a few parallels then do

  num_pattern_count <- sum(grepl("^[0-9]+_", res$STUDY))

  # If more than 2, remove the numeric prefix and print the message
  if (num_pattern_count > 2) {
    res$STUDY <- sub("^[0-9]+_", "", res$STUDY)
    print("Detected parallel peak finder run")
  }



#  res$STUDY <- sub("^[0-9]+_", "",   res$STUDY)




  blank_row <- as.data.frame(lapply(res, function(x) NA))

  # Combine the blank row and the original data frame using bind_rows from dplyr



  res <- dplyr::bind_rows(blank_row, res)
  #res <- bind_rows(blank_row, res)

  blank_row <- as.data.frame(lapply(res, function(x) NA))
  res <- dplyr::bind_rows(blank_row, res)

  blank_row <- as.data.frame(lapply(res, function(x) NA))
  res <- dplyr::bind_rows(blank_row, res)


  #Add necessary dummy rows

  res$RS[1] <- "rs99999999"
  res$RS[2] <- "-a-aaarModel"
  res$RS[3] <- "----a"

  #Create CIs


  if(Pre_Calculated_CIs == F)
  {



    if(Test_Statistic == "BETA")

    {

  #    print(res)
   #   print(str(res))

  res$LL <- res$BETA - 1.96*(res$SE)
  res$UL <- res$BETA + 1.96*(res$SE)

    }else{

      #Assigned OR to BETA earlier.
      res$LL <- res$BETA * exp(-1.96 * res$SE)
      res$UL <- res$BETA * exp(1.96 * res$SE)

    }




  res$UL <- as.numeric(res$UL)


  res$LL <- as.numeric(res$LL)


  }



  res$P <- sprintf(paste0("%.", P_Value_Resolution, "e"), res$P)


  #res$P <- sprintf("%.2e", res$P)
  res$UL <- as.numeric(res$UL)
  res$LL <- as.numeric(res$LL)




  #change this and legends later

  # Create a vector of letters to use as prefixes
  prefixes <- LETTERS[1:length(Names)]  # A, B, C, etc.
  postfixes <- LETTERS[1:length(Selected_SNPs)]  # A, B, C, etc.


  # Modify the res data frame by adding the letter prefix to each dataset
  res <- res %>%
    dplyr::mutate(Study = dplyr::case_when(
      STUDY %in% Names ~ paste(prefixes[match(STUDY, Names)], STUDY, sep = "-"),
      TRUE ~ STUDY  # Default if none match
    ))


  print(Selected_SNPs)
  print(res$RS)
  print(res$Backup_ID)
  print(res)
  print("check")
  res <- res %>%
    mutate(Backup_Single = stringr::str_extract(Backup_ID, "rs\\d+"))
  print(res)
  print(res$Backup_Single)


  if(Model_Reference == F)
  {




  # res <- res %>%
  #   dplyr::mutate(RS = dplyr::case_when( #need to allow RS label to work
  #     RS %in% Selected_SNPs %in%  Selected_SNPs ~ paste(postfixes[match(RS, Selected_SNPs)], RS, sep = "-"),
  #     TRUE ~ RS  # Default if none match
  #   ))
  res <- res %>%
    dplyr::mutate(RS = dplyr::case_when(
      RS %in% Selected_SNPs | Backup_Single %in% Selected_SNPs ~ paste(
        postfixes[match(ifelse(RS %in% Selected_SNPs, RS, Backup_Single), Selected_SNPs)],
        RS,
        sep = "-"
      ),
      TRUE ~ RS  # Default if none match
    ))

  }else{
    res <- res %>%
      dplyr::mutate(RS = dplyr::case_when(
        group %in% Selected_SNPs ~ paste(postfixes[match(group, Selected_SNPs)], RS, sep = "-"),
        TRUE ~ RS  # Default if none match
      ))

}


  res <- res %>%
    dplyr::mutate(
      RS = ifelse(
        grepl("Ref:", Left_Plot_Value),  # Check if 'Ref:' appears anywhere in 'Left_Plot_Value'
        sub("-", "-11A", RS),              # Insert 'A' after the first hyphen in 'RS'
        RS                               # Leave 'RS' unchanged otherwise
      )
    )


  print(res)
  print(res$RS)
  print(res$Backup_Single)



  if(Test_Statistic == "BETA")
  {



    maxcalc <-  max(res$UL, na.rm = T) * 1
    mincalc <-  min(res$LL, na.rm = T) * 1

  #  print(maxcalc)
  #  print(mincalc)

  }else{


    maxcalc <-  max(res$UL, na.rm = T) * 1
    mincalc <-  min(res$LL, na.rm = T) / 1
  }


  #Needs to be symmetrical

  if(Test_Statistic == "BETA")
  {
    midmax <- max(abs(mincalc), maxcalc)
  #  midmaxneg <- midmax * -1
  #  midmaxpos <- midmax * 1

    midmaxneg <- mincalc * 1
    midmaxpos <- maxcalc * 1

  }else{

    midmaxneg <- mincalc * 1
    midmaxpos <- maxcalc * 1


  }







  if(Test_Statistic == "BETA")
  {
    #left SNPs names

    mincalcL <- midmaxneg + (0.001 * midmaxneg) * (1 + (Line_Space/10))

    mincalcLFull <- mincalcL * 1 * (1 + (Left_Spaces/100))

    mincalcR <- midmaxpos + (0.001 * midmaxneg * -1) *  (1 + (Line_Space/10))

    mincalcRFull <- mincalcR * 1 * (1 + (Right_Spaces/100))

  }else{

    #maybe add exact =/- here adjusted for log10
    mincalcR <- midmaxpos * (1 + (0.001 * abs(midmaxneg) * (1 + (Line_Space / 10))))
    mincalcL <- midmaxneg * (1 - (0.001 * abs(midmaxneg) * (1 + (Line_Space / 10))))

#     mincalcL <- midmaxneg - ( (0.001 * midmaxneg) * (1 + (Line_Space/10)) )

    mincalcLFull <- mincalcL / 1 / (1 + (Left_Spaces/100))

#    mincalcR <- midmaxpos + (0.001 * midmaxneg)  *  (1 + (Line_Space/10))

    mincalcRFull <- mincalcR * 1 * (1 + (Right_Spaces/100))



  }


  # if(Test_Statistic == "OR")
  # {
  #   # Use log10 scaling to adjust left and right limits more symmetrically
  #   # Add a small proportion based on the log values to create space
  #   mincalcL <- log10(midmaxneg) - (0.01 * log10(midmaxneg)) * (1 + (Line_Space/100))
  #   mincalcLFull <- mincalcL / (1 + (Border_Space_Left/100))
  #
  #   mincalcR <- log10(midmaxpos) + (0.01 * log10(midmaxneg)) * (1 + (Line_Space/100))
  #   mincalcRFull <- mincalcR * (1 + (Border_Space_Right/100))
  # }


  midmaxneg1dp <- round(midmaxneg, 10)
  midmaxneg1dp <- min(res$LL, na.rm = T)

 # midmaxneg1dp <- round(midmaxneg1dp, 2)
  midmaxneg1dp <- floor(midmaxneg1dp * 100) / 100
  #round up/right to always give left space

  midmaxpos1dp <-  round(midmaxpos, 10)
  midmaxpos1dp <- max(res$UL, na.rm = T)

#  print(midmaxpos1dp)


#  midmaxpos1dp <-  round(midmaxpos1dp, 2)
  #then down
  midmaxpos1dp <- ceiling(midmaxpos1dp * 100) / 100



  res$BETA <- as.numeric(res$BETA)

  # Print the updated data frame


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



  # if(Model_Reference == F)
  # {
  #
  # if(length(unique_study) %% 2 == 0)
  # {
  #
  # #scale halfway
  # res <- res %>%
  #   # Reverse the factor levels of RS and assign increasing values starting from 0.5
  #   dplyr::mutate(mid_y = seq(1.5, by = 1, length.out = dplyr::n()))
  #
  # res <- res %>%
  #   # Reverse the factor levels of RS and arrange the data frame by the reversed order
  #   dplyr::arrange(forcats::fct_rev(RS)) %>%
  #   # Assign increasing values starting from 1.5
  #   dplyr::mutate(mid_y = seq(1.5, by = 1, length.out = dplyr::n()))
  #
  # }else{
  #
  #   #scale halfway
  #   res <- res %>%
  #     # Reverse the factor levels of RS and assign increasing values starting from 0.5
  #     dplyr::mutate(mid_y = seq(1, by = 1, length.out = dplyr::n()))
  #
  #   res <- res %>%
  #     # Reverse the factor levels of RS and arrange the data frame by the reversed order
  #     dplyr::arrange(forcats::fct_rev(RS)) %>%
  #     # Assign increasing values starting from 1.5
  #     dplyr::mutate(mid_y = seq(1, by = 1, length.out = dplyr::n()))
  #
  #
  # }
  #
  # }else{

    #In model ref need to remove dummy " for ref cols
    unique_study <- unique_study[unique_study != ""]


      res2 <- res %>%
        dplyr::arrange(RS)




      res2 <- res2 %>% dplyr::arrange(desc(dplyr::row_number()))


      print(res2)
   #   return(res2)

      first_non_na <- min(which(!is.na(res2$Study)))

      # Initialize the row_number column with NA
      res2$row_number <- NA

      # Variable to count the number of NA rows encountered in the Study column after the first non-NA row
      na_count <- 0

      # Start numbering at the first non-NA row in Study and stop after the second NA in Study
      for (i in first_non_na:nrow(res2)) {
        if (is.na(res2$Study[i])) {
          na_count <- na_count + 1
          # Stop the loop when the second NA row is encountered in the Study column
          if (na_count == 2) {
            break
          }
        } else {
          res2$row_number[i] <- i - first_non_na + 1
        }
      }



      # # Initialize the Plot_Value column with NA
      # res2$Plot_Value <- NA
      #
      # # Loop through the dataframe
      # i <- 1
      # while (i <= nrow(res2)) {
      #   # Skip rows with NA in row_number
      #   if (is.na(res2$row_number[i])) {
      #     i <- i + 1
      #     next
      #   }
      #
      #   # Check if the current and next row have the same ID and both have non-NA row_number
      #   if (i < nrow(res2) && res2$ID[i] == res2$ID[i + 1] && !is.na(res2$row_number[i + 1])) {
      #     # If two consecutive rows have the same ID, take the median of their row_number
      #     median_value <- median(c(res2$row_number[i], res2$row_number[i + 1]))
      #     res2$Plot_Value[i] <- median_value
      #     res2$Plot_Value[i + 1] <- median_value
      #     # Skip the next row since it's already handled
      #     i <- i + 2
      #   } else {
      #     # If the row has no consecutive match or next row is NA, set Plot_Value equal to row_number
      #     res2$Plot_Value[i] <- res2$row_number[i]
      #     i <- i + 1
      #   }
      # }
      #
      #


      # Initialize the Plot_Value column with NA
      res2$Plot_Value <- NA

      # Loop through the dataframe
      i <- 1
      while (i <= nrow(res2)) {

        # Skip rows with NA in row_number
        if (is.na(res2$row_number[i])) {
          i <- i + 1
          next
        }

        # Initialize a list to collect indices of consecutive rows with the same ID
        consecutive_indices <- c(i)

        # Check for consecutive rows with the same ID
        while (i < nrow(res2) && res2$ID[i] == res2$ID[i + 1] && !is.na(res2$row_number[i + 1])) {
          consecutive_indices <- c(consecutive_indices, i + 1)
          i <- i + 1
        }

        # Calculate the median of row_number for consecutive rows
        median_value <- median(res2$row_number[consecutive_indices])

        # Assign the median value to Plot_Value for each consecutive row
        res2$Plot_Value[consecutive_indices] <- median_value

        # Move to the next row after the current sequence
        i <- i + 1
      }

      # Preview the updated dataframe
      head(res2)

      # View the updated dataframe
     # print(res2)


      res2$Plot_Value <- res2$Plot_Value

      res <- res2


  #    print(res$Plot_Value)



    if(length(unique_study) %% 2 == 0)
    {

      #scale halfway
      res <- res %>%
        # Reverse the factor levels of RS and assign increasing values starting from 0.5
        dplyr::mutate(mid_y = seq(1.5, by = 1, length.out = dplyr::n()))

      res <- res %>%
        # Reverse the factor levels of RS and arrange the data frame by the reversed order
        dplyr::arrange(forcats::fct_rev(RS)) %>%
        # Assign increasing values starting from 1.5
        dplyr::mutate(mid_y = seq(1.5, by = 1, length.out = dplyr::n()))



     # print(res)

    }else{

      #scale halfway
      res <- res %>%
        # Reverse the factor levels of RS and assign increasing values starting from 0.5
        dplyr::mutate(mid_y = seq(1, by = 1, length.out = dplyr::n()))

      res <- res %>%
        # Reverse the factor levels of RS and arrange the data frame by the reversed order
        dplyr::arrange(forcats::fct_rev(RS)) %>%
        # Assign increasing values starting from 1.5
        dplyr::mutate(mid_y = seq(1, by = 1, length.out = dplyr::n()))


    }

  #}



  res$Left_Plot_Value[res$RS == "-aaa-rs99999999"] <- ""
  res$Left_Plot_Value[res$RS == "---a"] <- ""
  res$Left_Plot_Value[res$RS == "-a-aaarModel"] <- Left_Title
  res$Plot_Value[res$RS == "-a-aaarModel"] <- res$Overall_Row_Number[res$RS == "-a-aaarModel"]


 # print(res)





  print("Getting Plot Ready")


  # Open a dummy PDF device to measure string widths
  grDevices::pdf(file = NULL)

  # Step 1: Calculate full string widths
  full_widths <- grid::convertWidth(grid::stringWidth(res$Left_Plot_Value), unitTo = "npc", valueOnly = TRUE)

  # Step 2: Extract parts before and after <br>
  parts <- strsplit(res$Left_Plot_Value, "<br>", fixed = TRUE)

  # Step 3: Initialize width vectors
  before_widths <- numeric(length(res$Left_Plot_Value))
  after_widths <- numeric(length(res$Left_Plot_Value))

  # Loop through each row to measure both parts
  for (i in seq_along(parts)) {
    if (length(parts[[i]]) == 2) {  # Ensure we have both parts
      before_widths[i] <- grid::convertWidth(grid::stringWidth(parts[[i]][1]), unitTo = "npc", valueOnly = TRUE)
      after_widths[i] <- grid::convertWidth(grid::stringWidth(parts[[i]][2]), unitTo = "npc", valueOnly = TRUE)
    }
  }

  # Step 4: Compute the difference and create a space string equivalent to half of it
  diff_widths <- abs(before_widths - after_widths) / 2  # Half of the difference
  res$DIF <- sapply(diff_widths, function(w) paste(rep(" ", round(w * 100)), collapse = ""))  # Convert to spaces

  # Close the dummy PDF device
  dev.off()

  # Print the new DIF column
  print(res$DIF)
  print(res$Left_Plot_Value)

  print("Getting Plot Ready again")



  # res$Left_Plot_Value <- ifelse(res$Left_Plot_Value == "SNP",
  #                               res$Left_Plot_Value,
  #                               paste0(res$Left_Plot_Value, "   "))
  #
  # # Measure the visual width of each string in res$Left_Plot_Value
  # string_widths <- grid::convertWidth(grid::stringWidth(res$Left_Plot_Value), unitTo = "npc", valueOnly = TRUE)
  #
  # # Find the maximum visual width
  # max_width <- max(string_widths, na.rm = TRUE)
  #
  # print(max_width)
  #
  #
  #
  # # Measure the width of a single underscore
  # underscore_width <- grid::convertWidth(grid::stringWidth("_"), unitTo = "npc", valueOnly = TRUE)
  #
  # # Calculate the number of underscores needed to match the maximum width
  # num_underscores <- ceiling(max_width / underscore_width)
  #
  # # Create the padded underscore string
  # padded_value <- paste(rep("_", num_underscores), collapse = "")
  #
  #
  #
  #
  # print(padded_value)

  # Replace "SNP" with the visually matched underscore string
#  res$Left_Plot_Value[res$Left_Plot_Value == "SNP"] <- padded_value


  Left_Spaces_Dub <- Left_Spaces * 3
 # print(Left_Spaces)
  Left_Spaces <- strrep("Z", Left_Spaces)  # Reassign as a repeated "Z" string
#  print(Left_Spaces)


  res$Left_Spaces_Dub <- paste0(res$DIF, Left_Spaces)


  if(Double_Label == F)
  {
  # Ensure strings except "SNP" have spaces added
  res$Left_Plot_Value <- ifelse(res$Left_Plot_Value == "fake",
                                res$Left_Plot_Value,
                                paste0(res$Left_Plot_Value, Left_Spaces))
}


  if(Double_Label == T)
  {
  # res$Left_Plot_Value <- gsub("(<br>)", paste0(Left_Spaces, "\\1"), res$Left_Plot_Value)  # Insert SPACES before <br>
  # res$Left_Plot_Value <- gsub("($)", Left_Spaces_Dub, res$Left_Plot_Value, perl = TRUE)  # Insert SPACES at the end of the string
  # res_plot$Left_Title[res_plot$RS == "-aaa-rs99999999"] <- Left_Title


    # Apply the modifications ONLY to rows where RS is NOT "-aaa-rs99999999"
    res$Left_Plot_Value[res$RS != "-a-aaarModel"] <- gsub(
      "(<br>)", paste0(Left_Spaces, "\\1"), res$Left_Plot_Value[res$RS != "-aaa-rs99999999"]
    )

    res$Left_Plot_Value[res$RS != "-a-aaarModel"] <- gsub(
      "($)", paste0(Left_Spaces, "\\1"), res$Left_Plot_Value[res$RS != "-aaa-rs99999999"], perl = TRUE
    )

    # Keep assigning Left_Title for the specific RS as per your original code
  #  res_plot$Left_Title[res_plot$RS == "-aaa-rs99999999"] <- Left_Title

    res$Left_Plot_Value[res$RS == "-a-aaarModel"] <- paste0(res$Left_Plot_Value[res$RS == "-a-aaarModel"], Left_Spaces)

    print(table(res$RS))

    print(res$Left_Plot_Value)

  print("Getting Plot Ready again yes" )
#print(res$Left_Plot_Value)

  #get first bit (longer for width)
  res$Left_Plot_Value_Mini <- sub("<br>.*", "", res$Left_Plot_Value)  # Extract only the part before <br>

  }
#print(res$Left_Plot_Value)

  print(res$Left_Plot_Value_Mini)



#no bother
  # Measure the visual width of each string in res$Left_Plot_Value
#string_widths <- grid::convertWidth(grid::stringWidth(res$Left_Plot_Value), unitTo = "npc", valueOnly = TRUE)
  if(Double_Label == T)
  {

  grDevices::pdf(file = NULL)  # Open a dummy PDF device
 string_widths <- grid::convertWidth(grid::stringWidth(res$Left_Plot_Value_Mini), unitTo = "npc", valueOnly = TRUE)
 dev.off()

  }else{
    grDevices::pdf(file = NULL)  # Open a dummy PDF device
    string_widths <- grid::convertWidth(grid::stringWidth(res$Left_Plot_Value), unitTo = "npc", valueOnly = TRUE)
    dev.off()


  }
 # string_widths <- 20

print(string_widths)



  print("Getting Plot Ready again yes ok" )


  # Find the maximum visual width
  max_width <- max(string_widths, na.rm = TRUE)



  print(max_width)



  print("Getting Plot Ready 2")
  #print(max_width)

  # Measure the width of a single underscore
 # underscore_width <- grid::convertWidth(grid::stringWidth("---"), unitTo = "npc", valueOnly = TRUE)
  grDevices::pdf(file = NULL)
  underscore_width <- grid::convertWidth(grid::stringWidth("---"), unitTo = "npc", valueOnly = TRUE)
  dev.off()
 # underscore_width <- 20


  print(underscore_width)





  # Calculate the exact number of underscores needed to match the maximum width
  exact_num_underscores <- max_width / underscore_width

  #/100 due to massive font size required with this render was 8

  numbar <- ceiling(exact_num_underscores * (SNP_Stat_Text_Size / 19.5))



#  numbar <- ceiling(exact_num_underscores)

  print(numbar)


#numbar <- 3



  # Create the visually matched underscore string by trimming the exact width
  underscores <- paste(rep("---", floor(exact_num_underscores)), collapse = "") # Start with more underscores
  final_padded_value <- substr(underscores, 1, round(exact_num_underscores))    # Trim to match the exact length

#  print(final_padded_value)


 # print(res)





#
#   # Replace "SNP" with the visually matched underscore string
#   res$Left_Plot_Value[res$RS == "---a"] <- final_padded_value
#   res$Plot_Value <- ifelse(
#     res$RS == "---a",              # Condition to check
#     res$Overall_Row_Number ,        # Value to assign when condition is TRUE
#     res$Plot_Value                # Retain the original value otherwise
#   )

#  res$Left_Plot_Value[res$RS == "---a"] <- final_padded_value
#  res$Overall_Row_Number[res$RS == "---a"] <- 43  # Set Overall_Row_Number to 45

  # Assign Plot_Value for "---a" as -1 of the corresponding neighboring value
  res$Plot_Value <- ifelse(
    res$RS == "---a",
    res$Plot_Value[res$RS == "-a-aaarModel"] - 1, # Subtract 1 from the Plot_Value of the neighboring RS
    res$Plot_Value # Retain the original value otherwise
  )

  # Assign Plot_Value for "-aaa-rs99999999" as +1 of the corresponding neighboring value
  res$Plot_Value <- ifelse(
    res$RS == "-aaa-rs99999999",
    res$Plot_Value[res$RS == "-a-aaarModel"] + 1, # Add 1 to the Plot_Value of the neighboring RS
    res$Plot_Value # Retain the original value otherwise
  )


  print("Getting Plot Ready 3")


  #numbar <- 100

 # bold_line_string <- paste(rep("\u2501", 2), collapse = "")
#  print(bold_line_string)



#  res$Left_Plot_Value[res$RS == "---a"] <- bold_line_string

#  res$Left_Plot_Value[res$RS == "-aaa-rs99999999"] <- bold_line_string

  labmatch <-  paste(rep("\u2501", numbar), collapse = "")
#  labmatch <-  paste("SNP", Left_Spaces)

  res$Left_Plot_Value <- ifelse(
    res$RS == "---a",
    paste(rep("\u2501", numbar), collapse = ""),  # Repeat \u2501 'num_repeats' times
    res$Left_Plot_Value
  )


  res$Left_Plot_Value <- ifelse(
    res$RS == "-aaa-rs99999999",
    paste(rep("\u2501", numbar), collapse = ""),  # Repeat \u2501 'num_repeats' times
    res$Left_Plot_Value
  )

  max_row_num <- max(res$Overall_Row_Number, na.rm = TRUE)
#  print(max_row_num)


  print("Getting Plot Ready 4")


  # Define the necessary columns and their values
  new_row <- data.frame(
    Overall_Row_Number = 0,            # y-axis value for the new row
    Left_Plot_Value = "Custom Label",  # Label for y = 0
    Plot_Value = 0,
    RS = "Custom Label" # Plot value for alignment
  )

  # Identify the columns missing in the new row and add them as NA
  missing_cols <- setdiff(names(res), names(new_row))
  for (col in missing_cols) {
    new_row[[col]] <- NA
  }

  # Ensure column order matches the original dataframe
  new_row <- new_row[names(res)]

  # Append the new row to the dataframe
  res <- rbind(new_row, res)


  res$Left_Plot_Value <- ifelse(
    res$Left_Plot_Value == "Custom Label",
    paste(rep("\u2501", numbar), collapse = ""),  # Repeat \u2501 'num_repeats' times
    res$Left_Plot_Value
  )



  # res$Left_Plot_Value[res$Left_Plot_Value == "SNP"] <- "____"
#
 #  p <-
 #    res |>
 #    ggplot2::ggplot(ggplot2::aes( y = Overall_Row_Number)) +  # Use Plot_Value for y axis
 #    ggplot2::theme_classic() +
 #
 # # Keep the y-axis as integers but customize the labels to float
 #    ggplot2::scale_y_continuous(
 #      breaks = seq(0, max_row_num, by = 0.5),  # Integer breaks from 0 to the number of rows
 #     labels =  function(x) ifelse(x %in% res$Plot_Value, res$Left_Plot_Value[match(x, res$Plot_Value)], ""),  # Only label at Plot_Value positions
 #      limits = c(1,( max_row_num)) ,  # Set y-axis limits from 0 to the number of rows
 #      expand = ggplot2::expansion(add = c(1, 0.04)) # No extra padding
 #    ) + ggplot2::theme(
 #    axis.text.y = ggplot2::element_text(family = "Courier", size = 8, color = "black") # Replace "Arial Unicode MS" with an available font
 #  )


  print("Beginning Plot Object")

  p <- res |>
    ggplot2::ggplot(ggplot2::aes(y = Overall_Row_Number)) +  # Use Plot_Value for y-axis
    ggplot2::theme_classic() +

    # Keep the y-axis as integers but customize the labels with HTML for U2501
    ggplot2::scale_y_continuous(
      breaks = seq(0, max_row_num, by = 0.5),  # Integer breaks from 0 to the number of rows
      labels = function(x) {
        labels <- ifelse(
          x %in% res$Plot_Value,
          res$Left_Plot_Value[match(x, res$Plot_Value)],
          ""
        )


       # formatted_labels <- gsub("Z", "<span style='color:white;'>Z</span>", labels)
        formatted_labels <- gsub("Z", "<span style='color:#ffffff00;'>Z</span>", labels)
      #  formatted_labels <- gsub("Z", "<span style='opacity:0;'>Z</span>", labels)
      #  formatted_labels <- gsub("Z", "<span style='visibility:hidden;'>Z</span>", labels)
      #  formatted_labels <- gsub("Z", "<span style='color: rgba(0,0,0,0);'>Z</span>", labels)


        # Apply HTML styling for U2501
        ifelse(
          grepl("\u2501", formatted_labels),  # Check if label contains U2501 character
          paste0("<span style='font-family: Arial;font-size:8pt; color:black'>", formatted_labels, "</span>"),
          paste0("<span style='font-family: Arial;font-size:", SNP_Stat_Text_Size, "pt; color:black'>", formatted_labels, "</span>")   # Default 12pt for all others
        )
      },
      limits = c(1, max_row_num),  # Set y-axis limits from 0 to the number of rows
      expand = ggplot2::expansion(add = c(1, 0.03)),  # No extra padding
     ) +
    ggplot2::theme(
      axis.text.y = ggtext::element_markdown(
        family = "Courier",
        margin = ggplot2::margin(r = 0),  # No space between labels and axis
        vjust  = 0.58  # Adjust vertical alignment to center labels on the tick
      ) ,
      ,
     axis.ticks.y = ggplot2::element_blank(),  # Remove y-axis tick marks
      axis.ticks.length.y = ggplot2::unit(0, "cm")  # Ensure tick length is 0
    )





  # longest_label <- max(nchar(res$Left_Plot_Value), na.rm = TRUE)
  #
  # # Longest label
  # #longest_label <- "Your longest label here"  # Replace with the longest label
  # longest_label_width <- grid::convertWidth(
  #   grid::stringWidth(longest_label, gp = grid::gpar(fontfamily = "Courier", fontsize = 16)),
  #   unitTo = "mm",
  #   valueOnly = TRUE
  # )
  #
  # # Single U2501 character rendered at 8pt
  # u2501_width <- grid::convertWidth(
  #   grid::stringWidth("\u2501", gp = grid::gpar(fontfamily = "Courier", fontsize = 8)),
  #   unitTo = "mm",
  #   valueOnly = TRUE
  # )
  #
  # # Calculate how many U2501 characters fit
  # u2501_fit <- floor(longest_label_width / u2501_width)
  #
  # u2501_fit



  # p <- res |>
  #   ggplot(aes(y = Overall_Row_Number)) +
  #   theme_classic() +
  #
  #   # Customize y-axis labels with HTML styling, avoiding NA values and handling U2501 character
  #   scale_y_continuous(
  #     breaks = res$Plot_Value,  # Only show labels at Plot_Value positions
  #     labels = function(x) {
  #       labels <- res$Left_Plot_Value[match(x, res$Plot_Value)]
  #
  #       # Remove NA labels (replace with empty string)
  #       labels[is.na(labels)] <- ""
  #
  #       # Apply different font sizes based on U2501 presence
  #       ifelse(
  #         grepl("\u2501", labels),  # Check if label contains U2501 character
  #         paste0("<span style='font-size:8pt; color:black'>", labels, "</span>"),
  #         paste0("<span style='font-size:12pt; color:black'>", labels, "</span>")  # Default 20pt for all others
  #       )
  #     },
  #     limits = c(1, max(res$Overall_Row_Number)),
  #     expand = expansion(add = c(1, 0.04))
  #   ) +
  #   theme(
  #     axis.text.y = element_markdown(
  #       margin = margin(r = 0)  # No space between labels and axis
  #     ),
  #     axis.ticks.y = element_blank(),           # Remove tick marks
  #     axis.ticks.length = unit(0, "cm")         # Ensure no tick length
  #   )
  # Ensure the ggtext package is installed

# pz <- p

# ggplot2::ggsave("TESTING.jpg", plot = pz, width = Width, height = Height, units = "in", dpi = Quality)

# print("hel")
# print(u2501_fit)




 # print(p)

 # print(res$Left_Plot_Value)


#  print(labmatch)



#  print(res$RS)


#  print("failed below here")

  #Modify the plot
  # p <- p +
  #   ggplot2::theme(
  #     axis.ticks.y = ggplot2::element_blank(),                          # Remove tick marks
  #     axis.ticks.length = grid::unit(0, "cm"),                         # Ensure no tick length
  #     axis.text.y = ggplot2::element_text(
  #       margin = ggplot2::margin(r = 0)#,  # No space between labels and axis
  #      # hjust = 1  #,
  #     #  vjust = 0.8 # Align text to the right to touch the axis
  #     )
  #   )

 # print(p)



  #print(p)





#
#     if(Test_Statistic == "BETA")
#     {
#       p <- p + ggplot2::scale_x_continuous(limits = c(mincalcLFull, mincalcRFull), breaks=seq(midmaxneg1dp, midmaxpos1dp, X_Axis_Separation), labels = scales::number_format(accuracy = 10^- X_Axis_Text_Resolution) )
#     }else{
#
#       p <- p +  ggplot2::scale_x_continuous(limits = c(mincalcLFull, mincalcRFull), breaks=seq(midmaxneg1dp, midmaxpos1dp, X_Axis_Separation), trans = "log10",  labels = scales::number_format(accuracy = 10^- X_Axis_Text_Resolution) )
#
#
#     }


  # if(Test_Statistic == "BETA") {
  #   # Ensure that 0 is always in the breaks for BETA
  #   breaks <- c(0, seq(midmaxneg1dp, midmaxpos1dp, X_Axis_Separation))
  #   p <- p + ggplot2::scale_x_continuous(
  #     limits = c(mincalcLFull, mincalcRFull),
  #     breaks = breaks,
  #     labels = scales::number_format(accuracy = 10^-X_Axis_Text_Resolution)
  #   )
  # } else {
  #   # Ensure that 1 is always in the breaks for OR (log10)
  #   breaks <- c(1, seq(midmaxneg1dp, midmaxpos1dp, X_Axis_Separation))
  #   p <- p + ggplot2::scale_x_continuous(
  #     limits = c(mincalcLFull, mincalcRFull),
  #     breaks = breaks,
  #     trans = "log10",
  #     labels = scales::number_format(accuracy = 10^-X_Axis_Text_Resolution)
  #   )
  # }





#
#   if(Test_Statistic == "BETA") {
#     # Ensure that 0 is always in the breaks for BETA, but avoid too many labels near 0
#     breaks <- seq(midmaxneg1dp, midmaxpos1dp, X_Axis_Separation)
#
#     # Exclude breaks very close to 0, only include 0
#     breaks <- breaks[abs(breaks) > Null_Buffer]  # Adjust threshold if necessary
#     breaks <- c(0, breaks)  # Ensure 0 is still included
#
#     p <- p + ggplot2::scale_x_continuous(
#       limits = c(mincalcLFull, mincalcRFull),
#       breaks = breaks,
#       labels = scales::number_format(accuracy = 10^-X_Axis_Text_Resolution)
#     )
#  # }



#print()



#  ggplot2::ggsave("TESTING.jpg", plot = p, width = Width, height = Height, units = "in", dpi = Quality)



#  print(res$Left_Plot_Value)


  # Dynamically adjust limits based on y-axis label size
  # Measure the maximum width of y-axis labels
  max_label_width <- max(grid::stringWidth(res$Left_Plot_Value))  # Replace `res$SNP` with your y-axis label column
 # print(max_label_width)
 # zzz
  # Convert the label width into axis scale units

#no bother
  #label_width_in_units <- as.numeric(grid::convertWidth(max_label_width, "npc", valueOnly = TRUE))
# grDevices::pdf(file = NULL)
# label_width_in_units <- as.numeric(grid::convertWidth(max_label_width, "npc", valueOnly = TRUE))
# dev.off()
  label_width_in_units <- 20


  # Adjust proportionally to the current axis range
  x_axis_range <- abs(mincalcR - mincalcL)
  adjustment_factor <- label_width_in_units * x_axis_range  # Scale the label width to the axis range

  # Adjust the left and right limits based on the scaled label width
  adjusted_min_limit <- mincalcL - adjustment_factor
  adjusted_max_limit <- mincalcR + adjustment_factor
#
#   print(adjusted_min_limit)
#   print(adjusted_max_limit)



  # # Update the plot with adjusted x-axis limits
  # p <- p + ggplot2::scale_x_continuous(
  #   limits = c(adjusted_min_limit, adjusted_max_limit),  # Dynamically adjusted x-axis limits
  #   breaks = displayed_breaks,                          # Reuse filtered breaks
  #   labels = displayed_labels                           # Reuse filtered labels
  # )
  #


###

 # print(res)

  # Load required package
 # library(dplyr)

  # Ensure data is sorted by Overall_Row_Number
  res <- dplyr::arrange(res, Overall_Row_Number)

  # Identify min and max values for Overall_Row_Number
  min_row <- min(res$Overall_Row_Number)  # Find the lowest value
  max_row <- max(res$Overall_Row_Number)  # Find the maximum value

  # Filter to exclude Overall_Row_Number == 0 and stop 2 before max
  filtered_res <- dplyr::filter(res, Overall_Row_Number > min_row & Overall_Row_Number < (max_row - 2))

  #print(filtered_res)


  # Initialize variables
  count_list <- list()  # Store counts per change

  if(Model_Reference == T)
  {
  current_value <- filtered_res$group[1]
  }else{
  current_value <- filtered_res$Left_Plot_Value[1]  # Track the current Left_Plot_Value
  }
  count <- 0  # Count occurrences
  group_index <- 1  # Track position in list

  # Iterate through rows in the filtered dataset
  for (i in seq_len(nrow(filtered_res))) {
    if(Model_Reference == F)
    {
    left_value <- filtered_res$Left_Plot_Value[i]
    }else{
      left_value <- filtered_res$group[i]
    }

    # If Left_Plot_Value changes, store the count and reset
    if (left_value != current_value) {
      count_list[[paste0("Group_", group_index)]] <- count
      group_index <- group_index + 1  # Move to next group
      current_value <- left_value  # Update tracked value
      count <- 1  # Reset count
    } else {
      count <- count + 1  # Continue counting
    }
  }

  # Store last counted group
  count_list[[paste0("Group_", group_index)]] <- count

  # Print results
 # print(count_list)



  blocksize_list <- count_list


#  print(blocksize_list)
#  print(blocksize_list[[1]])

  #basiaclly redundant?

  print("Beginning Blocks")


  blocksize <- length(unique(res$Study))
 # print(res$Study)
#  print(blocksize)
  #gets how many studies/size of block shaded bit
  blocks <- length(unique(res$ID))

  if(Model_Reference == T)
  {
    blocks <- length(unique(res$group))
  }


  blocks <- blocks -1
  #same for 1 NA study
  blocksize <- blocksize -1




  ymin <- 1



  if(Model_Reference == T)
  {

    #dont alter original plot data
  res2 <- res %>%
    dplyr::arrange(RS)



  res_reversed <- res2 %>% dplyr::arrange(desc(dplyr::row_number()))

  # Create an empty list to store the counts in order
  reference_count_list <- list()

  # Loop through the reversed dataframe and count occurrences of each Reference category
  for (ref in unique(res_reversed$group)) {
    if (!is.na(ref)) {
      count <- sum(res_reversed$group == ref, na.rm = TRUE)
      reference_count_list[[ref]] <- count
    }
  }

  # Print the resulting list of Reference counts in bottom-to-top order


  reference_count_vector <- unlist(reference_count_list)

  # Convert the counts to a numeric vector in the same order
  reference_count_numeric <- as.numeric(reference_count_vector)

  # Print the numeric vector

  blocksizes <- reference_count_numeric
  }





  res$BETA2 <- round(res$BETA, digits = 2)



  # wrangle results into pre-plotting table form
  res_plot <- res

  res_plot$Left_Title[res_plot$RS == "-aaa-rs99999999"] <- Left_Title
  res_plot$P_Value_Title[res_plot$RS == "-aaa-rs99999999"] <- P_Value_Title
  res_plot$Test_Stat_Se_Title[res_plot$RS == "-aaa-rs99999999"] <- Test_Stat_Se_Title

  unique_study <- unique(na.omit(res$Study))


  class(res_plot$BETA2)


  negative_sign_length <- nchar(sub("^-", "", as.character(min(res_plot$BETA2[res_plot$BETA2 < 0], na.rm = TRUE))))

  res_plot$BETA2
  res_plot$BETA2 <- as.numeric(res_plot$BETA2)
  res_plot$SE <- as.numeric(res_plot$SE)
  res_plot$BETA2 <- ifelse(res_plot$BETA2 >= 0  | res_plot$BETA2 == 0.00, paste(rep("", negative_sign_length), sprintf('%.2f', res_plot$BETA2)), sprintf('%.2f', res_plot$BETA2))



  res_plot$SE <- sprintf('%.2f', res_plot$SE)

  if(Display_Test_Stat_Se_Column == T){

    res_plot$BETA2 <- paste0(res_plot$BETA2, " (", res_plot$SE, ")")

  }

  if(Display_Test_Stat_CI_Column == T){

    res_plot$LL <-   sprintf('%.2f', res_plot$LL)

    res_plot$UL <-   sprintf('%.2f', res_plot$UL)
    #res_plot$LL <- round(  res_plot$LL, digits = 2)
    #res_plot$UL <- round(  res_plot$UL, digits = 2)

    res_plot$BETA2 <- paste0(res_plot$BETA2, " (", res_plot$LL, "-", res_plot$UL, ")")

  }


  res_plot$BETA2 <- gsub(" -0.00", "-0.00", res_plot$BETA2)






  res_plot$P[res_plot$RS == "-aaa-rs99999999"] <- "p-value"

  res_plot$BETA2[res_plot$RS == "-aaa-rs99999999"] <- " BETA (SE)"


  res_plot$BETA2[res_plot$BETA2 == "NA (NA)"] <- " BETA (SE)"



  #res <- res_plot
  res_plot$P[res_plot$RS == "zzz"] <- ""
  res_plot$P[res_plot$RS == "---a"] <- ""
  res_plot$BETA2[res_plot$RS == "zzz"] <- ""
  res_plot$BETA2[res_plot$RS == "---a"] <- ""
  res_plot$P[res_plot$RS == "-a-aaarModel"] <- paste0(P_Value_Title, "    ")
  res_plot$BETA2[res_plot$RS == "-a-aaarModel"] <- Test_Stat_Se_Title
  res_plot$P[res_plot$RS == "-aaa-rs99999999"] <- ""
  res_plot$BETA2[res_plot$RS == "-aaa-rs99999999"] <- ""
 # res_plot$BETA2[res_plot$P == "1.00e+00"] <- "" # has to go first due to below
#  res_plot$P[res_plot$P == "1.00e+00"] <- ""


  print(.Machine$double.xmin)

  #Make 000 as min floating point value
 # res_plot$P[res_plot$P == "0.00e+00"] <- .Machine$double.xmin



  # Dynamically construct the string representation for 1.00e+00 and 0.00e+00
  one_string <- sprintf(paste0("%.", P_Value_Resolution, "e"), 1)
  zero_string <- sprintf(paste0("%.", P_Value_Resolution, "e"), 0)

  #print(one_string)
  #print(zero_string)



  # Replace "1.00e+00" dynamically
  res_plot$BETA2[res_plot$P == one_string] <- "" # Ensure this runs first
  res_plot$P[res_plot$P == one_string] <- ""

  # Replace "0.00e+00" dynamically with the smallest floating-point value
  # Format .Machine$double.xmin to match P_Value_Resolution
  formatted_min_value <- sprintf(paste0("%.", P_Value_Resolution, "e"), .Machine$double.xmin)

  # Replace "0.00e+00" dynamically with the formatted smallest floating-point value
  res_plot$P[res_plot$P == zero_string] <- formatted_min_value
 # print(res_plot)



  #formatting




  # if(Display_P_Value_Column == T)
  # {
  #
  # res$P_BETA_SE <- paste0(res_plot$P, "         ", res_plot$BETA2)
  # }else
  # {
  #   res$P_BETA_SE <- paste0(res_plot$BETA2)
  # }

  #no first

  P_Calc <- P_Stat_Spaces -1

  P_Stat_Spaces <- strrep("Z", P_Stat_Spaces)



  P_Stat_Spaces_3_Dig  <- strrep("Z", P_Calc)

  #P_Stat_Spaces_Title <- strrep("Z", P_Stat_Spaces_Title)

  P_Stat_Spaces_Adj <- paste0(P_Stat_Spaces, "ZZZ") #add extra 3 for this


  P_Stat_Spaces_Title <- paste0(P_Stat_Spaces, "") #title spaces!


  RS_condition <- grepl("-a-aaarModel", res$RS)

  if(Display_P_Value_Column == T) {
    # Adjust the space padding based on the length of the exponent part
    res$P_BETA_SE <- ifelse(
      nchar(gsub(".*e[\\+\\-]([0-9]+)", "\\1", res_plot$P)) == 3,  # Check if the exponent has 3 digits
      paste0(res_plot$P, P_Stat_Spaces_3_Dig, res_plot$BETA2),   # paste0(res_plot$P, P_Stat_Spaces, "\u2009",  "\u200a", res_plot$BETA2),  # 9 spaces for three-digit exponents  #titles???? actually leave the same
      paste0(res_plot$P, P_Stat_Spaces, res_plot$BETA2)  # 10 spaces for shorter exponents
    )


    res$P_BETA_SE <- ifelse(
      RS_condition,
      paste0(res_plot$P, P_Stat_Spaces_Title, res_plot$BETA2),
      res$P_BETA_SE  # Keep unchanged if RS does not match
    )


     } else {
    res$P_BETA_SE <- res_plot$BETA2
  }


#  print(res)

#"      "
#"      "

  #****



  Right_Spaces <- strrep("Z", Right_Spaces)



  res$P_BETA_SE <- ifelse(res$P_BETA_SE == "fake",
                                res$P_BETA_SE,
                                paste0(Right_Spaces, res$P_BETA_SE))

  # Measure the visual width of each string in res$Left_Plot_Value
#  string_widths <- grid::convertWidth(grid::stringWidth(res$P_BETA_SE), unitTo = "npc", valueOnly = TRUE)
 grDevices::pdf(file = NULL)
  string_widths <- grid::convertWidth(grid::stringWidth(res$P_BETA_SE), unitTo = "npc", valueOnly = TRUE)
  dev.off()
 # string_widths <- 20


  # Find the maximum visual width
  max_width <- max(string_widths, na.rm = TRUE)

 # print(max_width)

  # Measure the width of a single underscore

#  underscore_width <- grid::convertWidth(grid::stringWidth("---"), unitTo = "npc", valueOnly = TRUE)
 grDevices::pdf(file = NULL)
  underscore_width <- grid::convertWidth(grid::stringWidth("---"), unitTo = "npc", valueOnly = TRUE)
  dev.off()
#  underscore_width <- 20

 # print(underscore_width)




  # Calculate the exact number of underscores needed to match the maximum width
  exact_num_underscores <- max_width / underscore_width

#  numbar <- ceiling(exact_num_underscores)

  numbar <- floor(exact_num_underscores * (SNP_Stat_Text_Size /19.5))

 # print(numbar)



  # Create the visually matched underscore string by trimming the exact width
  underscores <- paste(rep("---", floor(exact_num_underscores)), collapse = "") # Start with more underscores
  final_padded_value <- substr(underscores, 1, round(exact_num_underscores))    # Trim to match the exact length

 # print(final_padded_value)


 # print(res)





  #
  #   # Replace "SNP" with the visually matched underscore string
  #   res$Left_Plot_Value[res$RS == "---a"] <- final_padded_value
  #   res$Plot_Value <- ifelse(
  #     res$RS == "---a",              # Condition to check
  #     res$Overall_Row_Number ,        # Value to assign when condition is TRUE
  #     res$Plot_Value                # Retain the original value otherwise
  #   )

  #  res$Left_Plot_Value[res$RS == "---a"] <- final_padded_value
  #  res$Overall_Row_Number[res$RS == "---a"] <- 43  # Set Overall_Row_Number to 45

  # Assign Plot_Value for "---a" as -1 of the corresponding neighboring value

  #
  #  res$Plot_Value <- ifelse(
  #   res$RS == "---a",
  #   res$Plot_Value[res$RS == "-a-aaarModel"] - 1, # Subtract 1 from the Plot_Value of the neighboring RS
  #   res$Plot_Value # Retain the original value otherwise
  # )
  #
  # # Assign Plot_Value for "-aaa-rs99999999" as +1 of the corresponding neighboring value
  # res$Plot_Value <- ifelse(
  #   res$RS == "-aaa-rs99999999",
  #   res$Plot_Value[res$RS == "-a-aaarModel"] + 1, # Add 1 to the Plot_Value of the neighboring RS
  #   res$Plot_Value # Retain the original value otherwise
  # )
  #



  # bold_line_string <- paste(rep("\u2501", 2), collapse = "")
  #  print(bold_line_string)



  #  res$Left_Plot_Value[res$RS == "---a"] <- bold_line_string

  #  res$Left_Plot_Value[res$RS == "-aaa-rs99999999"] <- bold_line_string


  res$P_BETA_SE <- ifelse(
    res$RS == "---a",
    paste(rep("\u2501", numbar), collapse = ""),  # Repeat \u2501 'num_repeats' times
    res$P_BETA_SE
  )


  res$P_BETA_SE <- ifelse(
    res$RS == "-aaa-rs99999999",
    paste(rep("\u2501", numbar), collapse = ""),  # Repeat \u2501 'num_repeats' times
    res$P_BETA_SE
  )

  max_row_num <- max(res$Overall_Row_Number, na.rm = TRUE)
 # print(max_row_num)





  res$P_BETA_SE <- ifelse(
    res$RS == "Custom Label",
    paste(rep("\u2501", numbar), collapse = ""),  # Repeat \u2501 'num_repeats' times
    res$P_BETA_SE
  )




  # res$P_BETA_SE <- ifelse(
  #   res$RS == "Custom Label",
  #   paste0("<span style='font-size:50pt;'>TESTING</span>"),
  #   res$P_BETA_SE
  # )



 # print(res$P_BETA_SE)

#  ****










#
#
  # p <- p +   ggplot2::guides(y.sec = ggh4x::guide_axis_manual(
  #   breaks = res$Overall_Row_Number  , labels = res$P_BETA_SE))
  #
  #
  # print(p)
#
#
# #This is where size occurs
#
#
#   #need to mod before moving
#   p <- p+ ggplot2::theme(
#     # Remove original y-axis text (left side)
#     axis.ticks.y.left = ggplot2::element_blank(),
#     axis.ticks.y.right = ggplot2::element_blank(),
#     axis.text.y.left   = ggplot2::element_text(size = SNP_Stat_Text_Size),
#    axis.text.y.right   = ggplot2::element_text(size = SNP_Stat_Text_Size)
# # axis.text.y.right = ggtext::element_markdown(size = SNP_Stat_Text_Size)
#     # Remove original y-axis ticks (left side)
#   )   # Remove original y-axis line (left side)
#
#
#
#
#
#
#
#   print("adjusting right")
#
#   p <- p +
#     ggplot2::theme(
#       axis.ticks.y.right = ggplot2::element_blank(),                          # Remove tick marks
#       axis.ticks.length.right = grid::unit(0, "cm"),                         # Ensure no tick length
#       axis.text.y.right = ggplot2::element_text(
#         margin = ggplot2::margin(r = 0) #,  # No space between labels and axis
#       #  hjust = 1  #,
#         #  vjust = 0.8 # Align text to the right to touch the axis
#       )
#     )
#



#
  # p <- p +
  #   ggplot2::theme(
  #     # Remove right-side tick marks
  #     axis.ticks.y.right = ggplot2::element_blank(),
  #     axis.ticks.length.right = grid::unit(0, "cm"),
  #
  #     # Apply formatted y-axis text on the right
  #     axis.text.y.right = ggtext::element_markdown(
  #       margin = ggplot2::margin(r = 0),  # No space between labels and axis
  #       size = SNP_Stat_Text_Size,
  #       vjust = 0.58
  #     )
  #   ) +
  #
  #   # Apply formatted labels based on U2501 character for right-side axis
  #   ggplot2::guides(y.sec = ggh4x::guide_axis_manual(
  #     breaks = res$Overall_Row_Number,
  #     labels = function(x) {
  #       labels <- res$P_BETA_SE[match(x, res$Overall_Row_Number)]
  #
  #       # Remove NA labels
  #       labels[is.na(labels)] <- ""
  #
  #       formatted_labels <- gsub("Z", "<span style='color:#ffffff00;'>Z</span>", labels)
  #
  #       # Apply different font sizes based on U2501 presence
  #       ifelse(
  #         grepl("\u2501", formatted_labels),  # Check if label contains U2501 character
  #         paste0("<span style='font-size:8pt; color:black'>", formatted_labels, "</span>"),
  #         paste0("<span style='font-size:", SNP_Stat_Text_Size, "pt; color:black'>", formatted_labels, "</span>")   # Default 12pt for all others
  #       )
  #     }
  #   ))


#
# p <- p + ggplot2::scale_y_continuous(
#   sec.axis = ggplot2::sec_axis(
#     ~.,
#     breaks = res$Overall_Row_Number,
#     labels = function(x) {
#       labels <- res$P_BETA_SE[match(x, res$Overall_Row_Number)]
#       labels[is.na(labels)] <- ""
#       formatted_labels <- gsub("Z", "<span style='color:#ffffff00;'>Z</span>", labels)
#       formatted_labels <- gsub("\\.\\.", "<span style='color:#ffffff00;'>..</span>", formatted_labels)
#
#       ifelse(
#         grepl("\u2501", formatted_labels),
#         paste0("<span style='font-size:8pt; color:black'>", formatted_labels, "</span>"),
#         paste0("<span style='font-size:", SNP_Stat_Text_Size, "pt; color:black'>", formatted_labels, "</span>")
#       )
#     }
#   )
# )


# print( res$P_BETA_SE)


#print(res)

  print(res$Left_Plot_Value)

 # res$Left_Plot_Value <- 1

  print(res$Left_Plot_Value)

  #res$Left_Plot_Value <- trimws(res$Left_Plot_Value)

#THIS IS THE REAL AXIS.
p <- res |>
  ggplot2::ggplot(ggplot2::aes(y = Overall_Row_Number)) +  # Use Plot_Value for y-axis
  ggplot2::theme_classic() +

  # Single scale_y_continuous handling both left and right axes
  ggplot2::scale_y_continuous(
    breaks = seq(0, max_row_num, by = 0.5),  # Integer breaks from 0 to max rows
    labels = function(x) {
      labels <- ifelse(
        x %in% res$Plot_Value,
        res$Left_Plot_Value[match(x, res$Plot_Value)],
        ""
      )

      formatted_labels <- gsub("Z", "<span style='color:#ffffff00;'>Z</span>", labels)

      ifelse(#this bit also controls left
        grepl("\u2501", formatted_labels),
        paste0("<span style='font-family: Courier2; font-size:70pt; color:black'>", formatted_labels, "</span>"),
        paste0("<span style='font-family: Courier2; font-size:", SNP_Stat_Text_Size, "pt; color:black'>", formatted_labels, "</span>")

      )

    },
    limits = c(1, max_row_num),  # Y-axis limits from 1 to max rows
    expand = ggplot2::expansion(add = c(1, 0.03)),  # No extra padding

    sec.axis = ggplot2::sec_axis(
      trans = ~.,  # Keep transformation the same
      breaks = res$Overall_Row_Number,
      labels = function(x) {
        labels <- res$P_BETA_SE[match(x, res$Overall_Row_Number)]
        labels[is.na(labels)] <- ""
        formatted_labels <- gsub("Z", "<span style='color:#ffffff00;'>Z</span>", labels)
        formatted_labels <- gsub("\\.\\.", "<span style='color:#ffffff00;'>..</span>", formatted_labels)

        ifelse(
          grepl("\u2501", formatted_labels),
          paste0("<span style='font-family: Courier2; font-size:70pt; color:black'>", formatted_labels, "</span>"),
          paste0("<span style='font-family: Courier2; font-size:", SNP_Stat_Text_Size, "pt; color:black'>", formatted_labels, "</span>")
        #  paste0("<span style='font-family: Courier2; font-size:180pt; color:black'>", formatted_labels, "</span>")
        )
      }
    )
  ) +

  ggplot2::theme(
    axis.text.y = ggtext::element_markdown(
   family = "Courier2", #this bit controls left, above controls right!
      margin = ggplot2::margin(l = 0, r = 0),  # No space between labels and axis
      vjust  = 0.642,  # Adjust vertical alignment to center labels on the tick - bigger = more down
      hjust = 1
    ),
     axis.text.y.right = ggplot2::element_text(
       margin = ggplot2::margin(l = 0, r = 0),  # Removes gap to the right of secondary y-axis
       hjust = 0  # Ensures right-alignment
     ),
   axis.text.y.left = ggplot2::element_text(
     margin = ggplot2::margin(l = 0, r = -1.2),  # Removes gap to the right of secondary y-axis r= -0.4
     hjust = 1  # Ensures right-alignment
   ),
   axis.ticks.y = ggplot2::element_blank(),  # Remove y-axis tick marks
   axis.ticks.length.y = ggplot2::unit(0, "cm")  # Ensure tick length is 0
  )


  # ggplot2::theme(
  #   axis.text.y = ggtext::element_markdown(  # Use Markdown for LEFT axis for consistency
  #     family = "Courier2",  # Use standard Courier (not Courier2, which may not exist)
  #     margin = ggplot2::margin(l = 0, r = 0),  # Remove unnecessary spacing
  #     vjust  = 0.635,  # Align labels correctly
  #     hjust = 1  # Right-align text
  #   ),
  #    axis.text.y.right = ggtext::element_markdown(  # Use Markdown for RIGHT axis as well
  #      family = "Courier2",
  #      margin = ggplot2::margin(l = 0, r = 0),  # Match left axis spacing
  #      vjust = 0.635,  # Align correctly
  #      hjust = 0  # Left-align text for right axis
  #    )
  #   ,
  #    axis.text.y.left = ggtext::element_markdown(  # Ensure consistency with left axis
  #      family = "Courier2",
  #      margin = ggplot2::margin(l = 0, r = -1.2),  # Match right axis spacing
  #      vjust = 0.635,
  #      hjust = 1 #,
  #     # vjust = 0.5# Right-align text
  #
  # ),     #   ),
  #    axis.ticks.y = ggplot2::element_blank(),  # Remove y-axis tick marks
  #    axis.ticks.length.y.left  = ggplot2::unit(0, "cm"),  # Ensure tick length is 0
  # axis.ticks.length.y.right  = ggplot2::unit(0, "cm")
  # )


p <- p +

  ggplot2::geom_point(ggplot2::aes(x=BETA, color = STUDY), shape=res$Shape, size=3) +
  ggplot2::geom_linerange(ggplot2::aes(xmin=LL, xmax=UL))




#print(p)



#  pz <- p

if(Legend_On == FALSE)
{
  p <- p +  ggplot2::theme(legend.position = "none")
}


#return(p)

#
# p <- p +  ggplot2::guides(y.sec = ggh4x::guide_axis_manual(
#     breaks = res$Overall_Row_Number,
#     labels = function(x) {
#       labels <- res$P_BETA_SE[match(x, res$Overall_Row_Number)]
#
#       # Remove NA labels
#       labels[is.na(labels)] <- ""
#
#       formatted_labels <- gsub("Z", "<span style='color:#ffffff00;'>Z</span>", labels) #invis
#       formatted_labels <- gsub("\\.\\.", "<span style='color:#ffffff00;'>Z</span>", formatted_labels)
#
#       # Apply different font sizes based on U2501 presence
#       ifelse(
#         grepl("\u2501", formatted_labels),  # Check if label contains U2501 character
#         paste0("<span style='font-size:8pt; color:black'>", formatted_labels, "</span>"),
#         paste0("<span style='font-size:", SNP_Stat_Text_Size, "pt; color:black'>", formatted_labels, "</span>")   # Default 12pt for all others
#       )
#     }
#   )) +

    # ggplot2::theme(
    #   # Remove right-side tick marks
    #   axis.ticks.y.right = ggplot2::element_blank(),
    #   axis.ticks.length.right = grid::unit(0, "cm"),
    #
    #   # Apply formatted y-axis text on the right
    #   axis.text.y.right = ggtext::element_markdown(
    #     margin = ggplot2::margin(r = 0),  # No space between labels and axis
    #     size = SNP_Stat_Text_Size,
    #     vjust = 0.62,
    #     lineheight = 0.8,
    #     family = "Courier"
    #   )
    # )
    #


print("Scaling Range")


#auto adjust
if (is.null(X_Axis_Separation)) {
  range_width <- max(abs(mincalc), abs(maxcalc))  # Get the max range
  X_Axis_Separation <- 10^floor(log10(range_width) - 1)  # Auto-calculate step size
}

print("Scaling Breaks")

if (Test_Statistic == "BETA") {
  # Calculate the maximum absolute range
  max_range <- max(abs(mincalc), abs(maxcalc))  # Ensure symmetry around 0
  max_range <- ceiling(max_range / X_Axis_Separation) * X_Axis_Separation  # Round to nearest separation

  # Generate symmetrical breaks around 0
  breaks <- seq(-max_range, max_range, X_Axis_Separation)

  # Apply Null_Buffer: Remove breaks close to 0, but ensure 0 is included
  breaks <- breaks[abs(breaks) > Null_Buffer]
  breaks <- sort(c(0, breaks))  # Ensure 0 is included and sorted symmetrically

  # Filter out breaks and labels beyond midmaxneg1dp and midmaxpos1dp
  displayed_breaks <- breaks[breaks >= mincalc & breaks <= maxcalc]
  displayed_labels <- ifelse(
    breaks >= mincalc & breaks <= maxcalc,
    scales::number_format(accuracy = 10^-X_Axis_Text_Resolution)(breaks),
    ""  # Hide labels beyond the range
  )

  # Ensure `breaks` and `labels` are of the same length
  displayed_labels <- displayed_labels[breaks %in% displayed_breaks]


  buffer <- (maxcalc - mincalc) * Line_Space

  print(buffer)
  print("Max")
  print(maxcalc)
  print("Min")
  print(mincalc)


  # Set axis with symmetrical breaks and filtered labels and ticks
  p <- p + ggplot2::scale_x_continuous(
    limits = c(mincalc - buffer, maxcalc + buffer), # Symmetrical axis limits - enough for strips to fit just under
    breaks = displayed_breaks,          # Filtered breaks for ticks
    labels = displayed_labels,
    expand = c(0, 0)
  #  expand = ggplot2::expansion(mult = c(Line_Space, Line_Space)) # Filtered labels
  ) +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_line(size = 0.5, colour = "black")  # Regular tick size
    )
}


else {
  # Ensure that 1 is always in the breaks for OR (log10), but avoid too many labels close to 1
  # breaks <- seq(midmaxneg1dp, midmaxpos1dp, X_Axis_Separation)
  #
  # # Exclude breaks very close to 1, only include 1
  # breaks <- breaks[abs(breaks - 1) > Null_Buffer]  # Adjust threshold if necessary
  # breaks <- c(1, breaks)  # Ensure 1 is still included
  #
  # p <- p + ggplot2::scale_x_continuous(
  #   limits = c(mincalcLFull, mincalcRFull),
  #   breaks = breaks,
  #   trans = "log10",
  #   labels = scales::number_format(accuracy = 10^-X_Axis_Text_Resolution)
  # )

  max_range <- max(abs(mincalc), abs(maxcalc))  # Ensure symmetry around 0
  max_range <- ceiling(max_range / X_Axis_Separation) * X_Axis_Separation  # Round to nearest separation

  # Generate symmetrical breaks around 0
  breaks <- seq(-max_range, max_range, X_Axis_Separation)

  # Apply Null_Buffer: Remove breaks close to 0, but ensure 0 is included
  breaks <- breaks[abs(breaks) > Null_Buffer]
  breaks <- sort(c(0, breaks))  # Ensure 0 is included and sorted symmetrically

  # Filter out breaks and labels beyond midmaxneg1dp and midmaxpos1dp
  displayed_breaks <- breaks[breaks >= mincalc & breaks <= maxcalc]
  displayed_labels <- ifelse(
    breaks >= mincalc & breaks <= maxcalc,
    scales::number_format(accuracy = 10^-X_Axis_Text_Resolution)(breaks),
    ""  # Hide labels beyond the range
  )

  # Ensure `breaks` and `labels` are of the same length
  displayed_labels <- displayed_labels[breaks %in% displayed_breaks]

  buffer <- (log10(maxcalc) - log10(mincalc)) * Line_Space  # Apply scaling in log space
  buffer <- 10^buffer  # Convert back from log scale

#  print(Line_Space)
#  print(maxcalc)
#  print(mincalc)
#  print(buffer)
  #z

  # Set axis with symmetrical breaks and filtered labels and ticks
  p <- p + ggplot2::scale_x_continuous(
    limits = c(mincalc / buffer, maxcalc * buffer),   # Symmetrical axis limits - enough for strips to fit just under
    breaks = displayed_breaks,          # Filtered breaks for ticks
    labels = displayed_labels,
    trans = "log10"#,
#    expand = ggplot2::expansion(mult = c(Line_Space, Line_Space)) # Filtered labels
  ) +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_line(size = 0.5, colour = "black")  # Regular tick size
    )
}



print("Scaling Breaks Done")

#ggplot2::ggsave("TESTING.jpg", plot = p, width = Width, height = Height, units = "in", dpi = Quality)





  if(Test_Statistic == "BETAx")
  {

    shift_axis_x <- function(p, x = 0) {
      g <- ggplot2::ggplotGrob(p)
      dummy <- data.frame(x = x)

      # Extract the original y-axis grob (axis-l)
      ax <- g[["grobs"]][g$layout$name == "axis-l"][[1]]



      # Add the new axis at x = 2 and remove the old y-axis components
      p + ggplot2::annotation_custom(grid::grobTree(ax, vp = grid::viewport(x = 1, width = sum(ax$width))),
                                     xmax = x, xmin = x) +
        ggplot2::geom_vline(ggplot2::aes(xintercept = x), data = dummy) +

        # Remove only the original y-axis (left side)
        ggplot2::theme(
          axis.text.y.left = ggplot2::element_blank(),  # Remove original y-axis text (left side)
          axis.ticks.y.left = ggplot2::element_blank(), # Remove original y-axis ticks (left side)
          axis.line.y.left = ggplot2::element_blank()   # Remove original y-axis line (left side)
        )
    }


  }





  if(Test_Statistic == "ORx")
  {


    shift_axis_x <- function(p, x = 1) {  # Use 1 as a default for log10 scale (since log10(1) = 0)
      g <- ggplot2::ggplotGrob(p)
      dummy <- data.frame(x = log10(x))  # Transform x position if using log10

      # Extract the original y-axis grob (axis-l)
      ax <- g[["grobs"]][g$layout$name == "axis-l"][[1]]

      axis_width <- sum(ax$width)
     # axis_width_npc <- grid::convertWidth(sum(ax$width), unitTo = "npc", valueOnly = TRUE)
   #   grDevices::pdf(file = NULL)
    #  axis_width_npc <- grid::convertWidth(sum(ax$width), unitTo = "npc", valueOnly = TRUE)
    #  dev.off()

     axis_width_npc <- 20

      # Get the original x-limits of the plot (data coordinates)
      xlims <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x.range

      # Calculate the data range of the x-axis
      x_range <- diff(xlims)

      # Convert the npc width to data coordinates (using the proportion of npc width relative to the x range)
      axis_width_data_units <- axis_width_npc * x_range




      newl <- mincalcL - axis_width_data_units


      # Adjust axis grob to the correct transformed x position
      p + ggplot2::annotation_custom(grid::grobTree(ax, vp = grid::viewport(x = 1, width = sum(ax$width))),
                                     xmax = log10(x), xmin = log10(x)) +  # Apply log10 to the x-axis values
        ggplot2::geom_vline(ggplot2::aes(xintercept = x), data = dummy) +  # Adjust vline based on transformed x

        # Remove the original y-axis
        ggplot2::theme(
          axis.text.y.left = ggplot2::element_blank(),  # Remove original y-axis text
          axis.ticks.y.left = ggplot2::element_blank(), # Remove original y-axis ticks
          axis.line.y.left = ggplot2::element_blank()   # Remove original y-axis line
        )


    }


  }



 # ggplot2::ggsave("Test.jpg", plot = p, width = Width, height = Height, units = "in", dpi = Quality)







  # Example usage:
#  p <- shift_axis_x(p, mincalcL)


#  p <- p +
#    ggplot2::geom_hline(
#      yintercept = 2,
#      color = "red",     # Change the color of the line as needed
#      linetype = "dashed" # Optional: Use a dashed line
#    )+
#    ggplot2::coord_cartesian(clip = "off") # Ensure the line fills the full plot area

#  print(p)


  # Print the updated plot
#  print(p)


#  ggplot2::ggsave("Test.jpg", plot = p, width = Width, height = Height, units = "in", dpi = Quality)





#print(mincalcL)
#print(p)

#zzz


print("Still going")


  if(Test_Statistic == "BETA")
  {

    shift_axis_y_right <- function(p, x = 0) {
      g <- ggplot2::ggplotGrob(p)
      dummy <- data.frame(x = x)

      # Extract the original right y-axis grob (axis-r)
      ax <- g[["grobs"]][g$layout$name == "axis-r"][[1]]

      # Add the new right axis at the specified x position and remove the old right y-axis components
      p + ggplot2::annotation_custom(grid::grobTree(ax, vp = grid::viewport(x = 1, width = sum(ax$width))),
                                     xmax = x, xmin = x) +
        ggplot2::geom_vline(ggplot2::aes(xintercept = x), data = dummy) +

        # Remove the original right y-axis elements (text, ticks, and line)
        ggplot2::theme(
          axis.text.y.right = ggplot2::element_blank(),  # Remove original right y-axis text
          axis.ticks.y.right = ggplot2::element_blank(), # Remove original right y-axis ticks
          axis.line.y.right = ggplot2::element_blank()   # Remove original right y-axis line
        )
    }

  }



  if(Test_Statistic == "OR")
  {


    shift_axis_y_right <- function(p, x = 1) {  # Default x = 1 for log10 scale (log10(1) = 0)
      g <- ggplot2::ggplotGrob(p)
      dummy <- data.frame(x = log10(x))  # Transform x position if using log10

      # Extract the original right y-axis grob (axis-r)
      ax <- g[["grobs"]][g$layout$name == "axis-r"][[1]]

      # Adjust axis grob to the correct transformed x position
      p + ggplot2::annotation_custom(grid::grobTree(ax, vp = grid::viewport(x = 1, width = sum(ax$width))),
                                     xmax = log10(x), xmin = log10(x)) +  # Apply log10 to the x-axis values
        ggplot2::geom_vline(ggplot2::aes(xintercept = x), data = dummy) +  # Adjust vline based on transformed x

        # Remove the original right y-axis
        ggplot2::theme(
          axis.text.y.right = ggplot2::element_blank(),  # Remove original right y-axis text
          axis.ticks.y.right = ggplot2::element_blank(), # Remove original right y-axis ticks
          axis.line.y.right = ggplot2::element_blank()   # Remove original right y-axis line
        )
    }

  }
  # Example usage:
#  p <- shift_axis_y_right(p, mincalcR)



#dont use

if(2==1)
{

  for (block in 1:blocks) {

    # Calculate ymax

    if(Model_Reference == T)
    {

    blocksize <- blocksizes[block]

    }

    ymax <- ymin + (blocksize)



    if(block == 1)
    {
      ymax <- ymax -0.5
    }



    #if you dont want beta(se) or ps on the right then limit to rect limit
    if(Display_Test_Stat_Se_Column == F & Display_P_Value_Column == F & Display_Test_Stat_CI_Column == F){


      mincalcR <- midmaxpos
      mincalcRFull <- midmaxpos

    }

    if(Display_Test_Stat_Se_Column == T & Display_P_Value_Column == F){


      mincalcR <- mincalcR
      mincalcRFull <- mincalcRFull

    }

    if(Display_Test_Stat_Se_Column == F & Display_Test_Stat_CI_Column == F  & Display_P_Value_Column == T){


      mincalcR <- mincalcR
      mincalcRFull <- mincalcRFull

    }


    print("Still going2")




    if(Strips == T)
    {


  #    print("Adding Strips")




    # Add annotation for the rectangle
   if(block %% 2 ==0)
  #    if(1 == 1)
    {


  #    print("running stirps")



      p <- p + ggplot2::annotate("rect", xmin = midmaxneg, xmax = midmaxpos, ymin = ymin, ymax = ymax,
                        alpha = .1, fill = Strip_Colour)



      # print(ymin)
      # print(ymax)

    }

    }

    ymin <- ymax
  }




}
#  pz <- p


print("Still going3")


#print(blocksize_list)

for (i in seq_along(blocksize_list)) {
  #print(i)
}


#print("Hi")



for (block in seq_along(blocksize_list)) {

#  print(block)


  # Calculate ymax

  # if(Model_Reference == T)
  # {
  #
  #   blocksize <- blocksizes[block]
  #
  # }

#  print((blocksize_list[[block]]))
#  print("Hi")


  ymax <- ymin + (blocksize_list[[block]])



  if(block == 1)
  {
    ymax <- ymax -0.5
  }



  #if you dont want beta(se) or ps on the right then limit to rect limit
  # if(Display_Test_Stat_Se_Column == F & Display_P_Value_Column == F & Display_Test_Stat_CI_Column == F){
  #
  #
  #   mincalcR <- midmaxpos
  #   mincalcRFull <- midmaxpos
  #
  # }
  #
  # if(Display_Test_Stat_Se_Column == T & Display_P_Value_Column == F){
  #
  #
  #   mincalcR <- mincalcR
  #   mincalcRFull <- mincalcRFull
  #
  # }
  #
  # if(Display_Test_Stat_Se_Column == F & Display_Test_Stat_CI_Column == F  & Display_P_Value_Column == T){
  #
  #
  #   mincalcR <- mincalcR
  #   mincalcRFull <- mincalcRFull
  #
  # }






  if(Strips == T)
  {


    #    print("Adding Strips")




    # Add annotation for the rectangle
    if(block %% 2 ==0)     #just for alternating block
      #    if(1 == 1)
    {


      #    print("running stirps")



      p <- p + ggplot2::annotate("rect", xmin = midmaxneg, xmax = midmaxpos, ymin = ymin, ymax = ymax,
                                 alpha = .1, fill = Strip_Colour)



      # print(ymin)
      # print(ymax)

  }

  }

  ymin <- ymax
}



print("Still going4")

if (num_pattern_count > 2) {


Names <- sub("^[0-9]+_", "",   Names)


}

#Names <- sub("^[0-9]+_", "",   Names)

#Names <- paste0(Names, "1")
#Data_Set_Colours <- unlist(stringr::str_extract_all(Data_Set_Colours, "#[0-9A-Fa-f]+"))

#print(Data_Set_Colours)
#z
#print(Names)




#Names <- sort(Names)
#Names <- rev(sort(Names))

#STUDY col is the clean one!
values <- setNames(Data_Set_Colours, Names)



labels <- setNames(Names, Names)

#HERO



#print(values)



#print(labels)

#print(Names)


#print(Data_Set_Colours)
#print("Hi")


#print(res$STUDY)


print("Still going5")





  p <- p + ggplot2::scale_color_manual(
     values = values,
    labels = labels,
    breaks = Names

  )





  if(Legend_On == TRUE)

  {

    #names is in order

  #  p <- p + ggplot2::guides(color = ggplot2::guide_legend(title = Legend_Title, override.aes = list(shape = 15)))
#
#     p <- p + ggplot2::guides(
#       color = ggplot2::guide_legend(
#         title = Legend_Title,
#         override.aes = list(shape = c(15, 15, 18, rep(15, length(Names) - 3)))
#
#       )
#     )

    p <- p + ggplot2::guides(
      color = ggplot2::guide_legend(
        title = Legend_Title,
        override.aes = list(shape = Shapes)
      )
    )

}





  print("Still going6")



  blocksize <- blocksize
  #gets how many studies/size of block shaded bit
  blocks <- blocks


  end <- (blocksize)*(blocks) + 1

 # print(end)


  end <- max(res$Overall_Row_Number) - 2


 # print(res)
  #z

#  print(end)



  if(Model_Reference == T)
  {
    end <- sum(blocksizes) + 1
  }



  if(Test_Statistic == "BETA")
  {
    p <- p +


      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 0,  yend = end, y =  -Inf), color = Null_Line_Colour, linetype = Null_Line_Type)+


      ggplot2::labs(x=X_Axis_Title, y="YL")

  }else{
    p <- p +

      ggplot2::geom_segment(ggplot2::aes(x = 1, xend = 1,  yend = end, y = -Inf), color = Null_Line_Colour, linetype = Null_Line_Type)+


      ggplot2::labs(x=X_Axis_Title, y="YL")
    p
}











  if(Legend_On == FALSE)
  {
    p_mid <- p +
      ggplot2::theme(#axis.line.y = ggplot2::element_blank(),
                     axis.ticks.y= ggplot2::element_blank(),
                     legend.position = "none",
                     #    axis.text.y.left= ggplot2::element_text(size = 5),
                     axis.title.y = ggplot2::element_blank())
    #      axis.title.y= element_blank())

  }else{

    p_mid <- p +
      ggplot2::theme(#axis.line.y = ggplot2::element_blank(),
                     axis.ticks.y= ggplot2::element_blank(),
                     legend.position = "bottom",
                     legend.box = "horizontal",
                     #    axis.text.y.left= ggplot2::element_text(size = 5),
                     axis.title.y = ggplot2::element_blank())
    #      axis.title.y= element_blank())


  }




  p_mid <- p_mid +   ggplot2::geom_hline(yintercept = (end), linetype = "solid", color = "black")+# + geom_hline(yintercept = 0, linetype = "solid", color = "black")
 #   ggplot2::geom_hline(yintercept = 1, linetype = "solid", color = "black")+
    ggplot2::geom_hline(yintercept = end+2, linetype = "solid", color = "black")+
    #+ geom_segment(aes(x = -Inf, xend = Inf, y = 36, yend = 36), color = "black", linetype = "solid")
    ggplot2::theme(#axis.line.x = ggplot2::element_blank(),
         legend.text = ggplot2::element_text(size = Legend_Text_Size, family = "Courier2", color ="black"),
         legend.title = ggplot2::element_text(size = Legend_Title_Size, family = "Courier2", color ="black"),
         axis.text.x = ggplot2::element_text(size = X_Axis_Text_Size, family = "Courier2", color ="black"),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 15, b = -X_Axis_Title_Size), family = "Courier2", vjust = -3, color = "black", size = X_Axis_Title_Size ))#,
  #      axis.ticks.y= element_blank(),
  #      axis.text.y= element_blank(),
  #      axis.title.y= element_blank())




 # print(p_mid)


 p_mid <- p_mid + ggplot2::theme(
  axis.ticks.length.x  = ggplot2::unit(0.4,"cm"),
   axis.text.x = ggplot2::element_text(vjust = -1, size = X_Axis_Text_Size, family = "Courier2"),
   legend.margin=ggplot2::margin(10,10,10,10)
 ) # adds more


#  print(p_mid)




#print(end)
#print(end + 2)









 # print(p_mid)

  if(X_Axis_Label == T)
  {

    p_mid <- p_mid + ggplot2::xlab("BETA\n")
  }
  if(X_Axis_Label == F)
  {

   # p_mid <- p_mid + ggplot2::theme(axis.title.x = ggplot2::element_blank())
    p_mid <- p_mid + ggplot2::xlab(" ")
  }


  if(Test_Statistic == "OR")
  {

   # p_mid <- p_mid + ggplot2::scale_x_continuous(trans = "log10")
  }

 # print(p_mid)


  p_mid2 <- p_mid


   #scale halfway
   res_plot <- res_plot %>%
     # Reverse the factor levels of RS and assign increasing values starting from 0.5
     dplyr::mutate(mid_y = seq(1.5, by = 1, length.out = dplyr::n()))

   res_plot <- res_plot %>%
     # Reverse the factor levels of RS and arrange the data frame by the reversed order
     dplyr::arrange(forcats::fct_rev(RS)) %>%
     # Assign increasing values starting from 1.5
     dplyr::mutate(mid_y = seq(1.5, by = 1, length.out = dplyr::n()))

  # res_plot[res_plot$RS == "zzz" ] <- 0

  # print(res_plot)



   p_left <-
     res_plot  |>
     ggplot2::ggplot(ggplot2::aes(y = forcats::fct_rev(res_plot$RS)))
   p_left


#   print(p_left)




   # Subset the index list to get labels to make white

   #Find index
   middle_index <- ceiling(blocksize / 2)
   index_list <- levels(forcats::fct_rev(res_plot$RS))
   index_list <- setdiff(index_list, c("rMZA", "aModel", "-a-aaarModelNA", "rs99999999NA", "-aaab"))
   index_list <- index_list[middle_index:length(index_list)]
   index_list <- c(rep("NA", (blocksize-1)), index_list)

   labels_to_make_white <- index_list[seq_along(index_list) %% (blocksize) == 0]
   labels_to_make_white <- setdiff(labels_to_make_white, c("rMZA", "aModel", "-a-aaarModelNA", "rs99999999NA", "-a-abModel",  "-aaab"))
   # Assuming labels_to_make_white already exists
   labels_to_make_white <- c(labels_to_make_white, "Location (HG38)", "-aaab")

   labels_to_exclude <- c("-aaaLocation (HG38)")

   # Remove the specified labels
   labels_to_make_white <- labels_to_make_white[!labels_to_make_white %in% labels_to_exclude]


  # print(labels_to_make_white)

   res_plot$Left_Title[res_plot$RS == "-aaa-rs99999999"] <- Left_Title
   res_plot$P_Value_Title[res_plot$RS == "-aaa-rs99999999"] <- P_Value_Title
   res_plot$Test_Stat_Se_Title[res_plot$RS == "-aaa-rs99999999"] <- Test_Stat_Se_Title




   unique_study <- unique(na.omit(res$Study))


  name <- paste0(File_Name, ".", File_Type)



#  print(p_mid2)



#  png("test_plot.png", type = "cairo")
#  print(p_mid2)
#  dev.off()

  #setwd("C:/Users/callumon/Downloads/Simple/Plots")

  print("SAVING...")

  return(p_mid2)

#  ggplot2::ggsave(name, plot = p_mid2, width = Width, height = Height, units = "in", limitsize = F, dpi = Quality)

#here
#
#   img <- jpeg::readJPEG(paste0(name))
#
#   # Create the PDF and insert the image
#   pdf_file <- paste0(name, "_Scrollable.pdf")
#   grDevices::pdf(file = pdf_file, width = Width, height = Height)
#   grid::grid.raster(img)
#   grDevices::dev.off()


  #to here


  # ggplot2::ggsave(
  #   name,  # Ensure the file is saved as a PDF
  #   plot = p_mid2,
  #   width = Width,
  #   height = Height,
  #   units = "in",
  #   dpi = Quality,
  #   device = "pdf"  # Specify PDF as the output format
  # )

 # ggplot2::ggsave(paste0(name, ".pdf"), plot = p_mid2, width = Width, height = Height, units = "in", dpi = Quality, device = "pdf")


#  return(Combined_Processed_Data)
}


