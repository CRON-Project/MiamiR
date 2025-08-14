
#Establishing example datasets .rda
setwd("C:/Users/callumon/Downloads")

LbDementia_Sum_Stats <- vroom::vroom("GCST90001390_buildGRCh38.tsv.gz")
Intelligence_Sum_Stats <- vroom::vroom("SavageJansen_2018_intelligence_metaanalysis.txt")
Household_Income_Sum_Stats <- vroom::vroom("HillWD_31844048_household_Income.txt.gz")
Fake_Demo_Data <- vroom("fake_demographic_data.txt")

#Create mini versions to save on space and load time
Intelligence_Sum_Stats <- dplyr::sample_n(Intelligence_Sum_Stats, 100000)
Household_Income_Sum_Stats <- dplyr::sample_n(Household_Income_Sum_Stats, 100000)
LbDementia_Sum_Stats <- dplyr::sample_n(LbDementia_Sum_Stats, 100000)

#Any DF
usethis::use_data()

#Any Package
usethis::use_package()

#Pipe Package
usethis::use_pipe() # for %>%

#vignettes - must build and install from source first
usethis::use_vignette()

#License - open source, nice default
usethis::use_mit_license()
usethis::use_mit_license("Callum O'Neill")


Household_Income_Sum_Stats$ID <- seq_len(nrow(Household_Income_Sum_Stats))
