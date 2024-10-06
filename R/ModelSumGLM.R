#' Data frame of covariate effects from glm object "Model", created using the Model_Munge() function.
#'
#' A dataset containing processed test statistics from the dementia glm model
#' The variables are as follows:
#'
#' \itemize{
#'   \item Covariate = Name of covariate which test statistics represent
#'   \item OR = Odds Ratio
#'   \item SE = Standard Error of the effect estimate/OR
#'   \item P = p-value
#'   \item group = The overall group which covariates, particularly factor covarites fall into
#'   \item Reference = The reference level for factor covariates, derived from meta-data
#'
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ModelSumGLM
#' @usage data(ModelSumGLM)
#' @format A data frame with 9 rows and 6 variables (simulated)
NULL
