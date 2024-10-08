% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LbDementia_Sum_Stats.R
\docType{data}
\name{LbDementia_Sum_Stats}
\alias{LbDementia_Sum_Stats}
\title{100,000 SNP Sample Subset of GWAS Summary Statistics of Lewy Body Dementia}
\format{
A data frame with 100,000 rows (randomly selected) and 10 variables from: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90001001-GCST90002000/GCST90001390/
}
\usage{
data(LbDementia_Sum_Stats)
}
\description{
A dataset containing a subset of the summary statistics from the Chia, R. 2021 Lewy Body Dementia GWAS
The variables are as follows:
}
\details{
\itemize{
\item chromosome = chromosome number
\item variant_id = SNP rs number
\item effect_allele_frequency = observed frequency of effect allele
\item base_pair_location = basepair position (GRCh38)
\item other_allele = The non-effect allele
\item effect_allele = The effect allele
\item beta = The regression weight
\item odds_ratio = The effect allele test statistic from a logistic regression model
\item standard_error = The standard error of the Beta
\item p_value = The P value of the association test.

}
}
\keyword{datasets}
