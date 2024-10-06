#' Data frame of covariate effects from lm object "Model2", created using the Model_Munge() function.
#'
#' A dataset containing processed test statistics from the BMI lm model
#' The variables are as follows:
#'
#' \itemize{
#'   \item Covariate = Name of covariate which test statistics represent
#'   \item BETA = BETA estimate/effect size
#'   \item SE = Standard Error of the effect estimate/OR
#'   \item P = p-value
#'   \item group = The overall group which covariates, particularly factor covarites fall into
#'   \item Reference = The reference level for factor covariates, derived from meta-data
#'
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ModelSumLM
#' @usage data(ModelSumLM)
#' @format A data frame with 9 rows and 6 variables (simulated)
NULL
