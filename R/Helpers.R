#' @noRd
detect_column <- function(Data, Candidate_Column = NULL, Allowed_Names, column_type = "column") {
  if (is.null(Candidate_Column)) {
    for (name in Allowed_Names) {
      if (name %in% colnames(Data)) {
        message(paste0("Using ", name, " as ", column_type, " column (automatic)"))
        return(name)
      }
    }
    stop(paste("No valid", column_type, "column found."), call. = FALSE)
  } else {
    message(paste0("Using ", Candidate_Column, " as ", column_type, " column (manual)"))
    return(Candidate_Column)
  }
}

#' @noRd
detect_chromosome_column <- function(Data, Chromosome_Column = NULL) {
  allowed <- c("chromosome", "chrom", "chr", "CHROM", "Chromosome", "CHR", "Chr", "Chrom")
  detect_column(Data, Chromosome_Column, allowed, column_type = "chromosome")
}

#' @noRd
detect_position_column <- function(Data, Position_Column = NULL) {
  allowed <- c("POS", "pos", "Pos", "BP", "BPos", "bpos", "BPOS", "bPos",
               "Position", "position", "POSITION", "genpos", "GENPOS", "Genpos",
               "base_pair_location")
  detect_column(Data, Position_Column, allowed, column_type = "position")
}

#' @noRd
detect_snp_column <- function(Data, SNP_ID_Column = NULL) {
  allowed <- c("ID", "Id", "RsID", "RsId", "RSID", "snp", "SNP", "Snp", "snv",
               "SNV", "Snv", "RS", "rs", "variant_id")
  detect_column(Data, SNP_ID_Column, allowed, column_type = "SNP ID")
}

#' @noRd
detect_pvalue_column <- function(Data, PValue_Column = NULL) {
  allowed <- c("P", "p", "Pvalue", "pvalue", "P-Value", "p-value", "p-Value", "P-VALUE",
               "logp", "LogP", "LOGP", "Logp", "log10p", "Log10P", "LOG10P", "-LOG10P",
               "p_value")

  P_col <- detect_column(Data, PValue_Column, allowed, column_type = "p-value")

  # If 'P' not in columns, but log10P column is detected, compute and assign
  if (!("P" %in% colnames(Data)) &&
      P_col %in% c("logp", "LogP", "LOGP", "Logp", "log10p", "Log10P", "LOG10P", "-LOG10P")) {

    message("No P column detected, computing from log10 P-value column: ", P_col)
    Data$P <- 10^-(as.numeric(Data[[P_col]]))
    assign(deparse(substitute(Data)), Data, envir = parent.frame())
    return("P")
  }

  return(P_col)
}



#' @noRd
detect_reference_allele_column <- function(Data, Reference_Allele_Column = NULL) {
  allowed <- c(
    "A0", "A2", "REF", "other_allele", "a0", "a2", "ref",
    "Non_effect_Allele", "Reference", "reference",
    "allele0", "allele2", "ALLELE0", "ALLELE2"
  )
  detect_column(Data, Reference_Allele_Column, allowed, column_type = "reference allele")
}

#' @noRd
detect_effect_allele_column <- function(Data, Effect_Allele_Column = NULL) {
  allowed <- c(
    "A1", "ALT", "a1", "alt", "Alternate", "Effect_Allele",
    "effect_allele", "alternate", "allele1", "ALLELE1"
  )
  detect_column(Data, Effect_Allele_Column, allowed, column_type = "effect allele")
}


#' @noRd
get_name_fast <- function(x) {
  name <- substitute(x)
  if (is.symbol(name)) {
    # Very fast
    as.character(name)
  } else {
    # Fallback for complex expressions
    deparse(name)
  }
}


#' @noRd
run_with_counter <- function(func, args = list(), session = NULL, default_val = NULL) {
  exprs <- as.list(body(func))[-1]
  total <- length(exprs)
  env <- new.env(parent = environment(func))

  # Inject arguments and defaults
  # Inject arguments and defaults
  defaults <- formals(func)
  for (name in names(defaults)) {
    if (!is.null(args[[name]])) {
      assign(name, args[[name]], envir = env)
    } else if (!is.symbol(defaults[[name]]) && !is.language(defaults[[name]])) {
      assign(name, eval(defaults[[name]], envir = env), envir = env)
    }
  }

  result <- NULL

  tryCatch({
    for (i in seq_along(exprs)) {
      pct <- round(100 * i / total)
      bar_length <- 30
      filled_len <- round(bar_length * pct / 100)
      bar <- paste0("[", strrep("=", filled_len), strrep(" ", bar_length - filled_len), "] ", sprintf("%3d%%", pct))

      expr_text <- paste(deparse(exprs[[i]]), collapse = " ")
      is_vroom_line <- grepl("vroom", expr_text)

      if (!is_vroom_line) cat("\r", bar)

      # Evaluate each line and only assign to result if it's the last line
      if (i == length(exprs)) {
        result <- eval(exprs[[i]], envir = env)
      } else {
        suppressMessages(suppressWarnings(
          capture.output(eval(exprs[[i]], envir = env), type = "output")
        ))
      }

      if (!is.null(session)) {
        session$sendCustomMessage("plot_progress", list(
          pct = pct,
          msg = paste("Line", i, "of", total)
        ))
      }

      flush.console()
    }

    cat("\r[", strrep("=", 30), "] 100%\n", sep = "")
    invisible(result)
  }, error = function(e) {
    message("run_with_counter(): error in line ", i, " - ", e$message)
    if (!is.null(session)) {
      session$sendCustomMessage("plot_progress", list(
        pct = NA,
        msg = paste("Error on line", i, ":", e$message)
      ))
    }
    return(default_val)
  })
}

