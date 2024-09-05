#Need these to gen and mod data coming in package
library("MiamiR")
library("vroom")
library("dplyr")
library("devtools")

Labelled_Data <- Annotate_Data(Data = Intelligence_Sum_Stats,
                        Chromosome_Column = "CHROM",
                        Position_Column = "GENPOS", SNP_ID_Column = "ID",
                        PValue_Column = "P", Build = "HG19")


print(head(Labelled_Data, n = 10))

Manhattan_Plot <- Single_Plot(Data = Labelled_Data,
                        Chromosome_Column = "CHROM",
                        Position_Column = "GENPOS", SNP_ID_Column = "ID",
                        PValue_Column = "P",
                        Title = "Intelligence",
                        Chromosome_Colours = c("blue", "turquoise"),
                        File_Name = "Manhattan_Plot")


print(Manhattan_Plot)

Miami_Plot <- Miami_Plot(Top_Data = Intelligence_Sum_Stats, Bottom_Data = Household_Income_Sum_Stats,
                       Top_Chromosome_Column = "CHROM",
                       Top_Position_Column = "GENPOS", Top_SNP_ID_Column = "ID",
                       Top_PValue_Column = "P",
                       Top_Colour_One = "green", Top_Colour_Two = "purple",
                       Bottom_Colour_One = "green", Bottom_Colour_Two = "purple",
                       Bottom_Chromosome_Column = "CHROM",
                       Bottom_Position_Column = "GENPOS", Bottom_SNP_ID_Column = "ID",
                       Bottom_PValue_Column = "P",
                       Top_Title = "Intelligence", Bottom_Title = "Household Income",
                       File_Name = "Miami_Plot", Width = 30, Height = 15, Quality = 900,
                       File_Type = "jpg")




print(Miami_Plot)




setwd("C:/Users/callumon/Downloads")
Intelligence_Sum_Stats <- vroom("SavageJansen_2018_intelligence_metaanalysis.txt")
Household_Income_Sum_Stats <- vroom("HillWD_31844048_household_Income.txt.gz")
#Create mini versions to save on space and load time
Intelligence_Sum_Stats <- sample_n(Intelligence_Sum_Stats, 100000)
Household_Income_Sum_Stats <- sample_n(Household_Income_Sum_Stats, 100000)

setwd("C:/Users/callumon/Miami_Package_R/MiamiR")
use_data(Intelligence_Sum_Stats, overwrite = T)
use_data(Household_Income_Sum_Stats, overwrite = T)


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

