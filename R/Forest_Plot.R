#' Title Create forest plots of test statistics of regression model covariates and GWAS summary statistics for SNPs of interest
#'
#' @param Data_Sets These are a list of GWAS summary statistics or munged model test statistics to be plotted; defaults to c("ModelSum", "ModelSum")
#' @param Names These are a list of names for the data sets to be labelled with in the forest plot - this also specifies the order down the page; defaults to c("Vals 1", "Vals 2")
#' @param Data_Set_Colours These are a list of colours for the data sets to be coloured with in the forest plot - this also specifies the order down the page; defaults to c("blue", "darkgreen")
#' @param Chromosome_Columns These are a list of manual chromosome column names for the data sets used in the order specified; defaults to c("chromosome", "chromosome")
#' @param Model_Reference Do you want to use statistics from an R regression model (T/F); defaults to T
#' @param Line_Space A continuous value specifying the boundary between the plot margin and the left and right lines; defaults to 1
#' @param Border_Space_Left A continuous value specifying the boundary between the left hand side plot labels and the left margin; defaults to 1
#' @param Border_Space_Right A continuous value specifying the boundary between the right hand side plot labels and the right margin; defaults to 1
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
#' @param Standard_Error_Columns These are a list of manual SE column names for the data sets used in the order specified; defaults to c("SE","SE"),
#' @param PValue_Columns These are a list of manual P column names for the data sets used in the order specified; defaults to c("P","P")
#' @param Match_Allele_Direction Do you want to match the allele directions to avoid flips (T/F); defaults to T
#' @param Match_Allele_Study If the above is T, which data sets should be the reference point for effect allele direction; defaults to "LbDementia_Sum_Stats"
#' @param Selected_SNPs A list of SNPs to be shown in the forest plot, automatically selected for from the summary statistics data set; defaults to c("rs2616526", "rs7974838", "rs59867714", "rs9571588", "rs79007041")
#' @param Selected_Covariates A list of covariates to be shown in the forest plot, automatically selected for from the munged model data set; defaults to c()
#' @param Reference_Alleles These are a list of manual reference allele column names for the data sets used in the order specified; defaults to c("other_allele","other_allel")
#' @param Effect_Alleles These are a list of manual effect allele column names for the data sets used in the order specified; defaults to c("effect_allele","effect_allele")
#' @param Upper_CI_Columns These are a list of manual CI UL column names for the data sets used in the order specified; defaults to c("UL","UL")
#' @param Lower_CI_Columns These are a list of manual CI LL column names for the data sets used in the order specified; defaults to c("LL","LL")
#' @param File_Name File name to save plot as; defaults to "Forest_Plot"
#' @param Width Width of saved plot; defaults to 10
#' @param Height Height of saved plot; defaults to 6
#' @param Quality Quality of saved plot (dpi); defaults to 900
#' @param File_Type File type of saved plot; defaults to "jpg"
#' @param X_Axis_Text_Resolution Number of decimal places to display on X axis text; defaults to 1
#'
#' @return Image of Single Forest Plot is saved to the current directory and ggplot object is saved
#' @export
#'
#' @examples Forest_Plot_Model_LM <- Forest_Plot(Data_Sets = c("ModelSumLM", "ModelSumLM"),
#'              Names = c("BMI LM (1)", "BMI LM (2)"),
#'              Model_Reference = TRUE,
#'              Test_Statistic = "BETA", #could be OR
#'              Display_P_Value_Column = TRUE,
#'              Display_Test_Stat_Se_Column = TRUE,
#'              X_Axis_Separation = 0.8,
#'              Border_Space_Left = 30,
#'              Border_Space_Right = 75,
#'              Strips = TRUE,
#'              Pre_Calculated_CIs = FALSE,
#'              Legend_Title = "Model",
#'              Left_Title = "Covariate",
#'              Test_Stat_Se_Title = "BETA (SE)",
#'              File_Name = "Forest_Plot", Width =10, Height = 6, Quality = 900,
#'              File_Type = "jpg"
#'              )
#'
#'
#'
#'
#' Forest_Plot_Model_GLM <- Forest_Plot(Data_Sets = c("ModelSumGLM", "ModelSumGLM"),
#'                                      Names = c("Dementia GLM (1)", "Dementia GLM (2)"),
#'                                      Model_Reference = TRUE,
#'                                      Test_Statistic = "OR",
#'                                      X_Axis_Text_Resolution = 2,
#'                                      Border_Space_Left = 2.8,
#'                                      Border_Space_Right = 4.7,
#'                                      Display_Test_Stat_CI_Column = TRUE,
#'                                      Display_P_Value_Column = TRUE,
#'                                      X_Axis_Separation = 0.02,
#'                                      Pre_Calculated_CIs = FALSE,
#'                                      Legend_Title = "Model",
#'                                      Left_Title = "Covariate",
#'                                      P_Value_Title = "p-value",
#'                                      Test_Stat_Se_Title = "OR (CI)",
#'                                      File_Name = "Forest_Plot", Width =10, Height = 6, Quality = 900,
#'                                      File_Type = "jpg"
#' )
#'
#'
#'
#' Forest_Plot_SNPs_BETA <- Forest_Plot(Data_Sets = c("Household_Income_Sum_Stats", "Intelligence_Sum_Stats"),
#'                                      Names = c("Income", "IQ"),
#'                                      Model_Reference = FALSE,
#'                                      Line_Space = 1.5,
#'                                      Border_Space_Right = 40,
#'                                      Border_Space_Left = 20,
#'                                      Test_Statistic = "BETA", #could be OR
#'                                      Display_Test_Stat_Se_Column = TRUE,
#'                                      Display_P_Value_Column = TRUE,
#'                                      X_Axis_Separation = 0.01,
#'                                      Pre_Calculated_CIs = FALSE,
#'                                      X_Axis_Text_Resolution = 2,
#'                                      Legend_Title = "Study",
#'                                      Left_Title = "SNP",
#'                                      P_Value_Title = "p-value",
#'                                      Test_Stat_Se_Title = "BETA (SE)",
#'                                      Match_Allele_Direction = TRUE,
#'                                      Match_Allele_Study = "Household_Income_Sum_Stats",
#'                                      Selected_SNPs = c("rs74832835",  "rs1157671",   "rs1790177",
#'                                                        "rs9508063",   "rs225682",  "rs56201315" ),
#'                                      File_Name = "Forest_Plot", Width =10, Height = 6, Quality = 900,
#'                                      File_Type = "jpg"
#' )
#'
#'
#'
#' Forest_Plot_SNPs_OR <- Forest_Plot(Data_Sets = c("LbDementia_Sum_Stats", "LbDementia_Sum_Stats"),
#'                                    Names = c("Dementia (1)", "Dementia (2)"),
#'                                    Model_Reference = FALSE,
#'                                    Line_Space = 1.1,
#'                                    Border_Space_Left = 7,
#'                                    Border_Space_Right = 30,
#'                                    Test_Statistic = "OR",
#'                                    Display_P_Value_Column = TRUE,
#'                                    Display_Test_Stat_CI_Column = TRUE,
#'                                    X_Axis_Separation = 0.2,
#'                                    Pre_Calculated_CIs = FALSE,
#'                                    Legend_Title = "Study",
#'                                    Left_Title = "SNP",
#'                                    P_Value_Title = "p-value",
#'                                    Test_Stat_Se_Title = "OR (CI)",
#'                                    Match_Allele_Direction = TRUE,
#'                                    Match_Allele_Study = "LbDementia_Sum_Stats",
#'                                    Selected_SNPs = c("rs59867714","rs7913723","rs79007041",
#'                                                      "rs34624328", "rs492457",
#'                                                      "rs7974838"  ),
#'                                    File_Name = "Forest_Plot", Width =10, Height = 6, Quality = 900,
#'                                    File_Type = "jpg"
#' )
#'
#'
#'
#'




Forest_Plot <- function(Data_Sets = c("ModelSum", "ModelSum"),
                        Names = c("Model 1", "Model 2"),
                        Data_Set_Colours = c("blue", "darkgreen"),
                        Chromosome_Columns = c(),
                        Model_Reference = TRUE,
                        Line_Space = 1,
                        Border_Space_Left = 2.75,
                        Border_Space_Right = 4.5,
                        Test_Statistic = "OR",
                        Display_Test_Stat_Se_Column = FALSE,
                        Display_Test_Stat_CI_Column = TRUE,
                        Display_P_Value_Column = TRUE,
                        Shapes = c("square", "diamond"),
                        Null_Line_Colour = "red",
                        Null_Line_Type = "dashed",
                        X_Axis_Title = "BETA",
                        X_Axis_Title_Size = 20,
                        X_Axis_Label = FALSE,
                        X_Axis_Separation = 0.05,
                        Strip_Colour = "lightblue",
                        Strips = TRUE,
                        X_Axis_Text_Resolution = 2,
                        Pre_Calculated_CIs = FALSE,
                        Legend_Title_Size = 15,
                        Legend_Text_Size = 12.5,
                        Legend_Title = "Study",
                        Left_Title = "SNP",
                        P_Value_Title = "p-value",
                        Test_Stat_Se_Title = "OR (CI)",
                        OR_Columns = c(),
                        Position_Columns = c() , SNP_ID_Columns = c(),
                        Beta_Columns = c(), Standard_Error_Columns = c(),
                        PValue_Columns = c(),
                        Match_Allele_Direction = FALSE,
                        Match_Allele_Study = "",
                        Selected_SNPs = c(),
                        Selected_Covariates = c(),
                        Reference_Alleles = c(),  Effect_Alleles = c(),
                        Upper_CI_Columns = c(), Lower_CI_Columns = c(),
                        File_Name = "Forest_Plot", Width =10, Height = 6, Quality = 600,
                        File_Type = "jpg"
                        )

  {


  Combined_Processed_Data <- data.frame()

 # print(Data_Sets)

  for (i in seq_along(Data_Sets)) {

    print(Data_Sets)

    # dataset_name <- Data_Sets[i]
    # Data <- dataset_name

    dataset_name <- Data_Sets[i]  # Get the dataset name

    print(dataset_name)

    Data <- get(dataset_name)  # Load the dataset using get()

    print(Data)

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


  print(Data)

  Data <- Data[Data$ID %in% Selected_SNPs,]


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

  print(Data)

  if(Model_Reference == F)
  {
  Data <- Data %>% dplyr::select(ID, ALLELE0, ALLELE1, CHROM, GENPOS, BETA, SE, P, STUDY, Shape)
  }

  if(Model_Reference == T)
  {
    Data <- Data %>% dplyr::select(ID,  BETA, SE, P, STUDY, Shape, group, Reference)

  }


  Combined_Processed_Data <- rbind(Combined_Processed_Data, Data)


  }






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



  if(Model_Reference == F)
  {

if(Match_Allele_Direction == T)

{

  if(Test_Statistic == "BETA")
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
            Combined_Processed_Data$BETA[j] <- Combined_Processed_Data$BETA[j] * -1
          }
        }
      }
    }
  }

  }


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

}

  Combined_Processed_Data$RS <- Combined_Processed_Data$ID

  res <- Combined_Processed_Data



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





  res$P <- sprintf("%.2e", res$P)
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

  if(Model_Reference == F)
  {
  res <- res %>%
    dplyr::mutate(RS = dplyr::case_when(
      RS %in% Selected_SNPs ~ paste(postfixes[match(RS, Selected_SNPs)], RS, sep = "-"),
      TRUE ~ RS  # Default if none match
    ))

  }else{
    res <- res %>%
      dplyr::mutate(RS = dplyr::case_when(
        group %in% Selected_SNPs ~ paste(postfixes[match(group, Selected_SNPs)], RS, sep = "-"),
        TRUE ~ RS  # Default if none match
      ))

}



  if(Test_Statistic == "BETA")
  {



    maxcalc <-  max(res$UL, na.rm = T) * 1
    mincalc <-  min(res$LL, na.rm = T) * 1
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
    mincalcL <- midmaxneg + (0.01 * midmaxneg) * (1 + (Line_Space/100))
    mincalcLFull <- mincalcL * 1 * (1 + (Border_Space_Left/100))
    mincalcR <- midmaxpos + (0.01 * midmaxneg * -1) *  (1 + (Line_Space/100))
    mincalcRFull <- mincalcR * 1 * (1 + (Border_Space_Right/100))

  }else{

    #maybe add exact =/- here adjusted for log10

    mincalcL <- midmaxneg - (0.01 * midmaxneg) * (1 + (Line_Space/100))
    mincalcLFull <- mincalcL / 1 / (1 + (Border_Space_Left/100))
    mincalcR <- midmaxpos + (0.01 * midmaxneg)  *  (1 + (Line_Space/100))
    mincalcRFull <- mincalcR * 1 * (1 + (Border_Space_Right/100))



  }


  midmaxneg1dp <- round(midmaxneg, 2)
  midmaxneg1dp <- min(res$LL, na.rm = T)

 # midmaxneg1dp <- round(midmaxneg1dp, 2)
  midmaxneg1dp <- ceiling(midmaxneg1dp * 100) / 100
  #round up/right to always give left space

  midmaxpos1dp <-  round(midmaxpos, 2)
  midmaxpos1dp <- max(res$UL, na.rm = T)
#  midmaxpos1dp <-  round(midmaxpos1dp, 2)
  #then down
  midmaxpos1dp <- floor(midmaxpos1dp * 100) / 100



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

        # Check if the current and next row have the same ID and both have non-NA row_number
        if (i < nrow(res2) && res2$ID[i] == res2$ID[i + 1] && !is.na(res2$row_number[i + 1])) {
          # If two consecutive rows have the same ID, take the median of their row_number
          median_value <- median(c(res2$row_number[i], res2$row_number[i + 1]))
          res2$Plot_Value[i] <- median_value
          res2$Plot_Value[i + 1] <- median_value
          # Skip the next row since it's already handled
          i <- i + 2
        } else {
          # If the row has no consecutive match or next row is NA, set Plot_Value equal to row_number
          res2$Plot_Value[i] <- res2$row_number[i]
          i <- i + 1
        }
      }

      # View the updated dataframe
     # print(res2)


      res2$Plot_Value <- res2$Plot_Value

      res <- res2


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



  p <-
    res |>
    ggplot2::ggplot(ggplot2::aes(y = Overall_Row_Number)) +  # Use Plot_Value for y axis
    ggplot2::theme_classic() +

    # Keep the y-axis as integers but customize the labels to float
    ggplot2::scale_y_continuous(
      breaks = seq(0, nrow(res), by = 0.5),  # Integer breaks from 0 to the number of rows
      labels = function(x) ifelse(x %in% res$Plot_Value, res$Left_Plot_Value[match(x, res$Plot_Value)], ""),  # Only label at Plot_Value positions
      limits = c(1,( nrow(res))),  # Set y-axis limits from 0 to the number of rows
      expand = ggplot2::expansion(add = c(1, 0.03))  # No extra padding
    )





    if(Test_Statistic == "BETA")
    {
      p <- p + ggplot2::scale_x_continuous(limits = c(mincalcLFull, mincalcRFull), breaks=seq(midmaxneg1dp, midmaxpos1dp, X_Axis_Separation), labels = scales::number_format(accuracy = 10^- X_Axis_Text_Resolution) )
    }else{

      p <- p +  ggplot2::scale_x_continuous(limits = c(mincalcLFull, mincalcRFull), breaks=seq(midmaxneg1dp, midmaxpos1dp, X_Axis_Separation), trans = "log10",  labels = scales::number_format(accuracy = 10^- X_Axis_Text_Resolution) )


    }






  p <- p +

    ggplot2::geom_point(ggplot2::aes(x=BETA, color = STUDY), shape=res$Shape, size=3) +
    ggplot2::geom_linerange(ggplot2::aes(xmin=LL, xmax=UL))







  blocksize <- length(unique(res$Study))
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
  res_plot$BETA2[res_plot$P == "1.00e+00"] <- "" # has to go first due to below
  res_plot$P[res_plot$P == "1.00e+00"] <- ""



  if(Display_P_Value_Column == T)
  {

  res$P_BETA_SE <- paste0(res_plot$P, "         ", res_plot$BETA2)
  }else
  {
    res$P_BETA_SE <- paste0(res_plot$BETA2)
  }

  p <- p +   ggplot2::guides(y.sec = ggh4x::guide_axis_manual(
    breaks = res$Overall_Row_Number  , labels = res$P_BETA_SE))



  #need to mod before moving
  p <- p+ ggplot2::theme(
    # Remove original y-axis text (left side)
    axis.ticks.y.left = ggplot2::element_blank(),
    axis.ticks.y.right = ggplot2::element_blank(),
    axis.text.y.left   = ggplot2::element_text(size = 12),
    axis.text.y.right   = ggplot2::element_text(size = 12)
    # Remove original y-axis ticks (left side)
  )   # Remove original y-axis line (left side)





  if(Test_Statistic == "BETA")
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

  if(Test_Statistic == "OR")
  {


    shift_axis_x <- function(p, x = 1) {  # Use 1 as a default for log10 scale (since log10(1) = 0)
      g <- ggplot2::ggplotGrob(p)
      dummy <- data.frame(x = log10(x))  # Transform x position if using log10

      # Extract the original y-axis grob (axis-l)
      ax <- g[["grobs"]][g$layout$name == "axis-l"][[1]]

      axis_width <- sum(ax$width)
      axis_width_npc <- grid::convertWidth(sum(ax$width), unitTo = "npc", valueOnly = TRUE)


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
  # Example usage:
  p <- shift_axis_x(p, mincalcL)






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
  p <- shift_axis_y_right(p, mincalcR)




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






    if(Strips == T)
    {



    # Add annotation for the rectangle
    if(block %% 2 ==0)
    {


      p <- p + ggplot2::annotate("rect", xmin = midmaxneg, xmax = midmaxpos, ymin = ymin, ymax = ymax,
                        alpha = .1, fill = Strip_Colour)


    }

    }

    ymin <- ymax
  }












#STUDY col is the clean one!
values <- setNames(Data_Set_Colours, Names)
labels <- setNames(Names, Names)





  p <- p + ggplot2::scale_color_manual(
     values = values,
    labels = labels,
    breaks = Names

  )



    p <- p + ggplot2::guides(color = ggplot2::guide_legend(title = Legend_Title, override.aes = list(shape = 15)))












  blocksize <- blocksize
  #gets how many studies/size of block shaded bit
  blocks <- blocks


  end <- (blocksize)*(blocks) + 1

  if(Model_Reference == T)
  {
    end <- sum(blocksizes) + 1
  }



  if(Test_Statistic == "BETA")
  {
    p <- p +


      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 0,  yend = end, y =  -Inf), color = Null_Line_Colour, linetype = Null_Line_Type)+


      ggplot2::labs(x=X_Axis_Title, y="YL")
    p
  }else{
    p <- p +

      ggplot2::geom_segment(ggplot2::aes(x = 1, xend = 1,  yend = end, y = -Inf), color = Null_Line_Colour, linetype = Null_Line_Type)+


      ggplot2::labs(x=X_Axis_Title, y="YL")
    p
}














  p_mid <- p +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(),
          axis.ticks.y= ggplot2::element_blank(),
          legend.position = "bottom",
          legend.box = "horizontal",
      #    axis.text.y.left= ggplot2::element_text(size = 5),
  axis.title.y = ggplot2::element_blank())
  #      axis.title.y= element_blank())





  p_mid <- p_mid + ggplot2::geom_hline(yintercept = (end), linetype = "solid", color = "black")+# + geom_hline(yintercept = 0, linetype = "solid", color = "black")
 #   ggplot2::geom_hline(yintercept = 1, linetype = "solid", color = "black")+
    ggplot2::geom_hline(yintercept = end+2, linetype = "solid", color = "black")+
    #+ geom_segment(aes(x = -Inf, xend = Inf, y = 36, yend = 36), color = "black", linetype = "solid")
    ggplot2::theme(#axis.line.x = ggplot2::element_blank(),
          legend.text = ggplot2::element_text(size = Legend_Text_Size),
          legend.title = ggplot2::element_text(size = Legend_Title_Size),
          axis.text.x = ggplot2::element_text(size = 15),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 15), vjust = -2, size = X_Axis_Title_Size ))#,
  #      axis.ticks.y= element_blank(),
  #      axis.text.y= element_blank(),
  #      axis.title.y= element_blank())





 p_mid <- p_mid + ggplot2::theme(
  axis.ticks.length.x  = ggplot2::unit(0.4,"cm"),
   axis.text.x = ggplot2::element_text(vjust = -1),
   legend.margin=ggplot2::margin(0,0,0,0)
 ) # adds more


 # print(p_mid)





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





#ggplot2::ggsave(name, plot = p_mid2, width = Width, height = Height, units = "in", dpi = Quality)


  return(p_mid2)
}


