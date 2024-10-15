#' Title Process model output to data frame which can be used for plotting test statistics for different covariates
#'
#' @param Model_Object This is the name of the base R model generated from running an lm() or glm(); defaults to Model
#' @param Model_Type This specifies the type of model run; defaults to "glm"
#'
#' @return Outputs the munged model data, ready for parsing to the Forest_Plot function
#' @export
#'
#' @examples ModelSumGLM <- Model_Munge(Model_Object = "Model", Model_Type = "glm")
#'           ModelSumLM <- Model_Munge(Model_Object = "Model2",  Model_Type = "lm")
#'
#'

Model_Munge <- function(Model_Object = Model, Model_Type = "glm")

{

print("Processing model")


Model <- get(Model_Object)

print("Object Found")


print(Model)

covariates <- attr(terms(Model), "term.labels")
ModelSum <- summary(Model)
ModelSum <- as.data.frame(ModelSum$coefficients)
ModelSum$Covariate <- rownames(ModelSum)
rownames(ModelSum) <- NULL
ModelSum <- ModelSum[ModelSum$Covariate != "(Intercept)",]
if(Model_Type == "glm")
{
ModelSum$OR <- exp(ModelSum$Estimate)
}
if(Model_Type == "lm")
{
  ModelSum$BETA <- ModelSum$Estimate #straight input for lm
}
ModelSum$Estimate <- NULL
ModelSum$SE <- ModelSum$`Std. Error`
ModelSum$`Std. Error` <- NULL
ModelSum$P <- ModelSum$`Pr(>|t|)`
ModelSum$`t value` <- NULL
ModelSum$`Pr(>|t|)` <- NULL


print("Data Extracted")

ModelSum$group <- sapply(ModelSum$Covariate, function(cov) {
  # Check for which covariate name is present in the Covariate string
  matched_group <- covariates[sapply(covariates, function(group) grepl(group, cov))]
 #bug fix for similar names like PC1 and PC10
#   matched_group <- covariates[sapply(covariates, function(group) grepl(paste0("^", group, "$"), cov))]

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


print("Processed Data:")

print(ModelSum)

return(ModelSum)

}





