#Need these to gen and mod data coming in package
library("MiamiR")
library("vroom")
library("dplyr")
library("devtools")

?vroom

setwd("C:/Users/callumon/Downloads")
Male <- vroom("Male_Stats_Plot.txt")
Female <- vroom("Female_Stats_Plot.txt")

Labelled_Data <- Annotate_Data(Data = Intelligence_Sum_Stats,
                        Chromosome_Column = "CHROM",
                        Position_Column = "GENPOS", SNP_ID_Column = "ID",
                        PValue_Column = "P", Build = "HG19")


print(head(Labelled_Data, n = 10))

Manhattan_Plot <- Single_Plot(Data = Male,
                        Chromosome_Column = "CHROM",
                        Position_Column = "GENPOS", SNP_ID_Column = "ID",
                        PValue_Column = "P",
                        Title = "Intelligence",
                        Chromosome_Colours = c("blue", "turquoise"),
                        File_Name = "Manhattan_Plot")


print(Manhattan_Plot)

Miami_Plot <- Miami_Plot(Top_Data = Male, Bottom_Data = Female,
                       Top_Chromosome_Column = "CHROM",
                       Top_Position_Column = "GENPOS", Top_SNP_ID_Column = "ID",
                       Top_PValue_Column = "P",
                       Top_Colour_One = "green", Top_Colour_Two = "purple",
                       Bottom_Colour_One = "green", Bottom_Colour_Two = "purple",
                       Bottom_Chromosome_Column = "CHROM",
                       Bottom_Position_Column = "GENPOS", Bottom_SNP_ID_Column = "ID",
                       Bottom_PValue_Column = "P",
                       Top_Title = "Intelligence", Bottom_Title = "Household Income",
                       File_Name = "Miami_Plot_Error", Width = 30, Height = 15, Quality = 900,
                       File_Type = "jpg")





Forest_Plot_Model_LM <- Forest_Plot(Data_Sets = c("ModelSumLM", "ModelSumLM"),
                                     Names = c("BMI LM (1)", "BMI LM (2)"),
                                     Model_Reference = TRUE,
                                     Test_Statistic = "BETA", #could be OR
                                     Display_P_Value_Column = TRUE,
                                     Display_Test_Stat_Se_Column = TRUE,
                                     X_Axis_Separation = 0.8,
                                     Border_Space_Left = 30,
                                     Border_Space_Right = 75,
                                     Strips = TRUE,
                                     Pre_Calculated_CIs = FALSE,
                                     Legend_Title = "Model",
                                     Left_Title = "Covariate",
                                     Test_Stat_Se_Title = "BETA (SE)",
                                     File_Name = "Forest_Plot", Width =10, Height = 6, Quality = 900,
                                     File_Type = "jpg"
)




Forest_Plot_Model_GLM <- Forest_Plot(Data_Sets = c("ModelSumGLM", "ModelSumGLM"),
                        Names = c("Dementia GLM (1)", "Dementia GLM (2)"),
                        Model_Reference = TRUE,
                        Test_Statistic = "OR",
                        X_Axis_Text_Resolution = 2,
                        Border_Space_Left = 2.8,
                        Border_Space_Right = 4.7,
                        Line_Space = 1,
                        Display_Test_Stat_CI_Column = TRUE,
                        Display_P_Value_Column = TRUE,
                        X_Axis_Separation = 0.02,
                        Pre_Calculated_CIs = FALSE,
                        Legend_Title = "Model",
                        Left_Title = "Covariate",
                        P_Value_Title = "p-value",
                        Test_Stat_Se_Title = "OR (CI)",
                        File_Name = "Forest_Plot", Width =10, Height = 6, Quality = 900,
                        File_Type = "jpg"
)



Forest_Plot_SNPs_BETA <- Forest_Plot(Data_Sets = c("Household_Income_Sum_Stats", "Intelligence_Sum_Stats"),
                        Names = c("Income", "IQ"),
                        Model_Reference = FALSE,
                        Line_Space = 1.5,
                        Border_Space_Right = 40,
                        Border_Space_Left = 20,
                        Test_Statistic = "BETA", #could be OR
                        Display_Test_Stat_Se_Column = TRUE,
                        Display_P_Value_Column = TRUE,
                        X_Axis_Separation = 0.01,
                        Pre_Calculated_CIs = FALSE,
                        X_Axis_Text_Resolution = 2,
                        Legend_Title = "Study",
                        Left_Title = "SNP",
                        P_Value_Title = "p-value",
                        Test_Stat_Se_Title = "BETA (SE)",
                        Match_Allele_Direction = TRUE,
                        Match_Allele_Study = "Household_Income_Sum_Stats",
                        Selected_SNPs = c("rs74832835",  "rs1157671",   "rs1790177",
                                          "rs9508063",   "rs225682",  "rs56201315" ),
                        File_Name = "Forest_Plot", Width =10, Height = 6, Quality = 900,
                        File_Type = "jpg"
)



Forest_Plot_SNPs_OR <- Forest_Plot(Data_Sets = c("LbDementia_Sum_Stats", "LbDementia_Sum_Stats"),
                                     Names = c("Dementia (1)", "Dementia (2)"),
                                     Model_Reference = FALSE,
                                     Line_Space = 1.1,
                                     Border_Space_Left = 7,
                                     Border_Space_Right = 30,
                                     Test_Statistic = "OR",
                                     Display_P_Value_Column = TRUE,
                                     Display_Test_Stat_CI_Column = TRUE,
                                     X_Axis_Separation = 0.2,
                                     Pre_Calculated_CIs = FALSE,
                                     Legend_Title = "Study",
                                     Left_Title = "SNP",
                                     P_Value_Title = "p-value",
                                     Test_Stat_Se_Title = "OR (CI)",
                                     Match_Allele_Direction = TRUE,
                                     Match_Allele_Study = "LbDementia_Sum_Stats",
                                     Selected_SNPs = c("rs59867714","rs7913723","rs79007041",
                                                       "rs34624328", "rs492457",
                                                       "rs7974838"  ),
                                     File_Name = "Forest_Plot", Width =10, Height = 6, Quality = 900,
                                     File_Type = "jpg"
)


ModelSumGLM <- Model_Munge(Model_Object = "Model", Model_Type = "glm") #Dementia
ModelSumLM <- Model_Munge(Model_Object = "Model2",  Model_Type = "lm") # BMI




setwd("C:/Users/callumon/Downloads")

Fake_Demo_Data <- vroom("fake_demographic_data.txt")

#use_data(Fake_Demo_Data)
#use_data(Model)

Model <- glm(Dementia ~ Ethnicity + Age + Sex + BMI + Location, data = Fake_Demo_Data)


Model2 <- lm(BMI ~ Ethnicity + Age + Sex + Dementia + Location, data = Fake_Demo_Data)

use_mit_license()

covariates <- attr(terms(Model), "term.labels")
ModelSum <- summary(Model)
ModelSum <- as.data.frame(ModelSum$coefficients)
ModelSum$Covariate <- rownames(ModelSum)
rownames(ModelSum) <- NULL
ModelSum <- ModelSum[ModelSum$Covariate != "(Intercept)",]
ModelSum$OR <- exp(ModelSum$Estimate)
ModelSum$Estimate <- NULL
ModelSum$SE <- ModelSum$`Std. Error`
ModelSum$`Std. Error` <- NULL
ModelSum$P <- ModelSum$`Pr(>|t|)`
ModelSum$`t value` <- NULL
ModelSum$`Pr(>|t|)` <- NULL

ModelSum$group <- sapply(ModelSum$Covariate, function(cov) {
  # Check for which covariate name is present in the Covariate string
  matched_group <- covariates[sapply(covariates, function(group) grepl(group, cov))]

  # If a match is found, return the group, otherwise return NA
  if (length(matched_group) > 0) {
    return(matched_group)
  } else {
    return(NA)
  }
})


# Now remove the group part from the Covariate string, but only if they are different
ModelSum$Covariate <- mapply(function(cov, group) {
  if (!is.na(group) && cov != group) {
    return(gsub(group, "", cov))  # Remove the group from Covariate string
  } else {
    return(cov)  # Don't change if covariate is the same as group
  }
}, ModelSum$Covariate, ModelSum$group)


# Add the 'Reference' column by comparing the full levels with those already in the Covariate column
ModelSum$Reference <- mapply(function(covariate_group, covariate) {
  # Get all levels for the current group from Model$xlevels
  all_levels <- Model$xlevels[[covariate_group]]

  # If it's not a categorical variable, return NA
  if (is.null(all_levels)) {
    return("None")
  }

  # Find the level that is not present in the Covariate column for this group
  non_reference_levels <- ModelSum$Covariate[ModelSum$group == covariate_group]
  reference_level <- setdiff(all_levels, non_reference_levels)

  # Return the reference level (if any)
  if (length(reference_level) > 0) {
    return(reference_level)
  } else {
    return("NA")
  }
}, ModelSum$group, ModelSum$Covariate)


use_data(ModelSum)
# View the updated dataframe



LbDementia_Sum_Stats <- vroom("GCST90001390_buildGRCh38.tsv.gz")
Intelligence_Sum_Stats <- vroom("SavageJansen_2018_intelligence_metaanalysis.txt")
Household_Income_Sum_Stats <- vroom("HillWD_31844048_household_Income.txt.gz")
#Create mini versions to save on space and load time
Intelligence_Sum_Stats <- sample_n(Intelligence_Sum_Stats, 100000)
Household_Income_Sum_Stats <- sample_n(Household_Income_Sum_Stats, 100000)
LbDementia_Sum_Stats <- sample_n(LbDementia_Sum_Stats, 100000)

setwd("C:/Users/callumon/Miami_Package_R/MiamiR")
use_data(Intelligence_Sum_Stats, overwrite = T)
use_data(Household_Income_Sum_Stats, overwrite = T)
use_data(LbDementia_Sum_Stats, overwrite = T)

#Packages

use_package("biomaRt")
use_package("dplyr")
use_package("ggmanh")
use_package("ggplot2")
use_package("ggtext")
use_package("tidyr")
usethis::use_pipe() # for %>%


#vignettes - must build and install from source first
use_vignette("Introduction_To_MiamiR")

