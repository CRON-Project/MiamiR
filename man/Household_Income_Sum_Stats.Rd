% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Household_Income_Sum_Stats.R
\docType{data}
\name{Household_Income_Sum_Stats}
\alias{Household_Income_Sum_Stats}
\title{100,000 SNP Sample Subset of GWAS Summary Statistics of Household Income}
\format{
A data frame with 100,000 rows (randomly selected) and 8 variables from: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009523/
}
\usage{
data(Household_Income_Sum_Stats)
}
\description{
A dataset containing a subset of the summary statistics from the Hill, W.D. 2019 Household Income GWAS
The variables are as follows:
}
\details{
\itemize{
\item Chr = chromosome number
\item SNP = SNP rs number
\item BPos = basepair position (GRCh37)
\item Non_effect_Allele = The non-effect allele
\item Effect_Allele = The effect allele
\item Beta = The regression weight
\item Standard_Error_of_Beta = The standard error of the Beta
\item P = The P value of the association test.

}
}
\keyword{datasets}
