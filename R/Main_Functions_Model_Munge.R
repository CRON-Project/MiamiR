
#' Munged, Covariate & Modelling Data
#'
#' @param Model_Object Base R model generated from running an lm() or glm(); defaults to NULL
#' @param ... Shadow argument to detect arguments passed; defaults to list()
#' @param Verbose Prevent display of progress bar as function is running and instead show key milestone outputs/messages (mainly for debugging purposes); defaults to FALSE.
#'
#' @return Munged dataframe amenable to plotting test statistics which can then be passed to other MiamiR functions
#' @export
#'
#'
#' @examples  Model_Munged <- Model_Munge(Model_Object = "Model_One")

    Model_Munge <- function(

      Model_Object = NULL,
      Verbose = FALSE,
      ...)

    {

    Model_Type <- "Unknown"

    message("Processing model")

    Model <- get(Model_Object)

    if(inherits(Model, "glm"))
    {
      Model_Type <- "glm"
    }

    if(Model_Type == "glm")

    {
    if(Model$family$family == "binomial")
    {
      Special <- "binomial"

    }else{

     Special <- "none"

    }

    }else{

      Special <- "none"

    }

    if(inherits(Model, "lm") & Special == "none")
    {
      Model_Type <- "lm"
    }


    message("Object Found")

    covariates <- attr(terms(Model), "term.labels")
    ModelSum <- summary(Model)

    message("Extracting test statistics")

    ModelSum <- as.data.frame(ModelSum$coefficients)

    message("Delineating covariates")

    ModelSum$Covariate <- rownames(ModelSum)

    rownames(ModelSum) <- NULL

    message("Ignoring intercept")

    ModelSum <- ModelSum[ModelSum$Covariate != "(Intercept)",]

    message("Determining model type")

    if(Model_Type == "glm")
    {
    ModelSum$OR <- exp(ModelSum$Estimate)
    }
    if(Model_Type == "lm")
    {
      ModelSum$BETA <- ModelSum$Estimate #straight input for lm
    }

    message("Assigning test statistic values")

    ModelSum$Estimate <- NULL
    ModelSum$SE <- ModelSum$`Std. Error`
    ModelSum$`Std. Error` <- NULL

    if(Special == "binomial")
    {
    ModelSum$P <- ModelSum$`Pr(>|z|)`
    }else{
    ModelSum$P <- ModelSum$`Pr(>|t|)`
    }
    ModelSum$`t value` <- NULL
    if(Special == "binomial")
    {
    ModelSum$`Pr(>|z|)` <- NULL
    ModelSum$`z value` <- NULL
    }else
    {
    ModelSum$`Pr(>|t|)` <- NULL
    }



    message("Data Extracted")

    ModelSum$group <- sapply(ModelSum$Covariate, function(cov) {
      parts <- unlist(strsplit(cov, ":"))

      matched_groups <- covariates[sapply(covariates, function(group)
        any(grepl(group, parts))
      )]

      if (length(matched_groups) > 0) {

        return(paste(matched_groups, collapse = ":")) # Keep full interaction mapping

      } else {

        return(NA)

      }
    })


    # Now remove the group part from the Covariate string, but only if they are different

    message("Mapping covariates")

      ModelSum$Covariate <- mapply(function(cov, grp) {

      if (!is.na(grp)) {

        # Split both covariate and group by ':'
        cov_parts <- strsplit(cov, ":")[[1]]
        grp_parts <- strsplit(grp, ":")[[1]]

        # Keep numeric (continuous) covariates untouched, clean factors
        cleaned_parts <- mapply(function(c, g) {

          all_levels <- Model$xlevels[[g]]

          if (!is.null(all_levels)) {


            return(sub(paste0("^", g), "", c))

          } else {

            # Continuous variable - keep as is

            return(c)
          }
        }, cov_parts, grp_parts, SIMPLIFY = TRUE)

        paste(cleaned_parts, collapse = ":")
      } else {

        cov

      }
    }, ModelSum$Covariate, ModelSum$group)

    message("Determining reference levels")

    ModelSum$Reference <- mapply(function(covariate_group, covariate) {

      # Handle interaction terms (e.g., "Ethnicity:Sex")
      if (grepl(":", covariate_group)) {

        group_parts <- unlist(strsplit(covariate_group, ":"))

        # Get reference level for each part
        ref_levels <- sapply(group_parts, function(g) {

          all_levels <- Model$xlevels[[g]]
          if (is.null(all_levels)) {

            return("None")

          }

          # Determine reference level (level not present as dummy var)
          non_reference_levels <- ModelSum$Covariate[ModelSum$group == g]
          reference_level <- setdiff(all_levels, non_reference_levels)

          if (length(reference_level) > 0) reference_level else "NA"

        }, USE.NAMES = FALSE)

        # Combine into one reference string
        return(paste(ref_levels, collapse = ":"))
      }

      # Existing logic for single terms
      all_levels <- Model$xlevels[[covariate_group]]

      # If it's not a categorical variable, return NA
      if (is.null(all_levels)) {

        return("None")

      }

      # Find the level that is not present in the Covariate column for this group
      non_reference_levels <- ModelSum$Covariate[ModelSum$group == covariate_group]
      reference_level <- setdiff(all_levels, non_reference_levels)

      if (length(reference_level) > 0) {

        return(reference_level)

      } else {

        return("NA")
      }

    }, ModelSum$group, ModelSum$Covariate)

    message("Processed Data")

    return(ModelSum)

    }


    .Model_Munge_original <- Model_Munge

    Model_Munge <- function(..., session = NULL) {

      args <- list(...)
      args$session <- session

      # Find the name of the Data argument if present
      data_arg_name <- if ("Data" %in% names(args)) {

        deparse(substitute(Data), backtick = TRUE)

      } else if ("Model_Object" %in% names(args)) {

        args$Model_Object

      } else {

        NULL
      }

      if (!is.null(data_arg_name)) {

        message(sprintf("Processing dataset: %s", data_arg_name))

      }

      verbose_mode <- isTRUE(args$Verbose)

      if (verbose_mode) {

        return(do.call(.Model_Munge_original, args))

      } else {

        return(
          suppressMessages(
            suppressWarnings(
              run_with_counter(
                func    = .Model_Munge_original,
                args    = args,
                session = session
              )
            )
          )
        )
      }
    }
