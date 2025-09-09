

#' @noRd
detect_column <- function(Data, Candidate_Column = NULL, Allowed_Names, column_type = "column") {
  if (is.null(Candidate_Column)) {
    for (name in Allowed_Names) {
      if (name %in% colnames(Data)) {
        message(paste0("Using ", name, " as ", column_type, " column (automatic)"))
        return(name)
      }
    }
    message(paste("No valid", column_type, "column found."))
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
detect_upper_ci_column <- function(Data, Upper_CI_Column = NULL) {
  allowed <- c("UL", "Upper_CI", "CIU")
  detect_column(Data, Upper_CI_Column, allowed, column_type = "upper CI")
}

#' @noRd
detect_lower_ci_column <- function(Data, Lower_CI_Column = NULL) {
  allowed <- c("LL", "Lower_CI", "CIL")
  detect_column(Data, Lower_CI_Column, allowed, column_type = "lower CI")
}

#' @noRd
detect_beta_column <- function(Data, Beta_Column = NULL) {
  allowed <- c("BETA", "Beta", "beta", "B", "stdBeta")
  detect_column(Data, Beta_Column, allowed, column_type = "beta")
}

#' @noRd
detect_or_column <- function(Data, OR_Column = NULL) {
  allowed <- c("OR", "or", "Or", "Odds", "Odd Ratio", "odds_ratio")
  detect_column(Data, OR_Column, allowed, column_type = "OR")
}

#' @noRd
detect_se_column <- function(Data, SE_Column = NULL) {
  allowed <- c("SE", "se", "Se", "Standard_Error", "standard_error", "Standard_Error_of_Beta")
  detect_column(Data, SE_Column, allowed, column_type = "standard error")
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
ensure_P <- function(Data) {
  p_names   <- c("P","p","Pvalue","pvalue","P-Value","p-value","p-Value","P-VALUE","p_value")
  log_names <- c("logp","LogP","LOGP","Logp","log10p","Log10P","LOG10P","-LOG10P")

  # If any P-style column exists, keep as-is (optionally you could standardise to P here)
  hit_p <- intersect(p_names, names(Data))
  if (length(hit_p) > 0) return(Data)

  # Otherwise, if any log10P-style column exists, create P = 10^-(log10P)
  hit_log <- intersect(log_names, names(Data))
  if (length(hit_log) > 0) {
    src <- hit_log[1]                # take the first match; reorder log_names to set priority
    Data$P <- 10^-(as.numeric(Data[[src]]))
    return(Data)
  }

  stop("No p-value or log10(p) columns found. Looked for: ",
       paste(c(p_names, log_names), collapse = ", "))
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
run_with_counter <- function(func, args = list(), session = NULL, default_val = NULL) {
  exprs <- as.list(body(func))[-1]
  total <- length(exprs)
  env <- new.env(parent = environment(func))

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

#' @noRd
save_all_plots <- function(

    obj,
    base_dir = ".",
    dpi = 300,
    ext = "jpg",
    single_width = 30,
    regional_width = 30,
    forest_width = 10
) {

  #Helpers
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

  sanitize <- function(x) {

    x <- as.character(x)
    x <- gsub("[/\\:*?\"<>|]+", "_", x, perl = TRUE)
    x <- gsub("\\s+", "_", trimws(x))
    x <- gsub("_+", "_", x)
    x <- sub("^_+", "", x)
    x <- sub("_+$", "", x)
    x[!nzchar(x)] <- "unnamed"
    x

  }

  get_dyn_h <- function(x, default = 30) {

    h <- attr(x, "dynamic_height")

    if (is.numeric(h) && is.finite(h)) h else default

  }

  # choose a graphics device when can't ggsave
  open_device <- function(path, w, h, dpi) {

    ext <- tolower(tools::file_ext(path))

    if (ext %in% c("jpg","jpeg")) {

      grDevices::jpeg(path, width = w, height = h, units = "in", res = dpi, quality = 95)

    } else if (ext == "png") {

      grDevices::png(path, width = w, height = h, units = "in", res = dpi)

    } else if (ext %in% c("tif","tiff")) {

      grDevices::tiff(path, width = w, height = h, units = "in", res = dpi, compression = "lzw")

    } else if (ext == "pdf") {

      grDevices::pdf(path, width = w, height = h, onefile = FALSE)

    } else {

      grDevices::png(path, width = w, height = h, units = "in", res = dpi)

    }
  }

  # robust save that tolerates ggplot, patchwork, grob/gtable...
  save_any <- function(plot_obj, path, width, height, dpi) {

    ok <- TRUE

    tryCatch({

      ggplot2::ggsave(path, plot = plot_obj, width = width, height = height,
                      units = "in", dpi = dpi, limitsize = FALSE)
    }, error = function(e) ok <<- FALSE)

    if (ok) return(invisible(path))

    ok2 <- TRUE

    tryCatch({
      ggplot2::ggsave(path, plot = ggplotify::as.ggplot(plot_obj), width = width, height = height,
                      units = "in", dpi = dpi, limitsize = FALSE)
    }, error = function(e) ok2 <<- FALSE)

    if (ok2) return(invisible(path))

    open_device(path, width, height, dpi)
    grid::grid.newpage(); grid::grid.draw(plot_obj); grDevices::dev.off()
    invisible(path)
  }

  dataset_from_title <- function(title) {

    if (is.null(title)) return(NA_character_)

    title <- as.character(title)

    if (!nzchar(title)) return(NA_character_)

    if (grepl(":", title, fixed = TRUE))

      return(trimws(sub("^\\s*([^:]+):.*$", "\\1", title)))

    NA_character_

  }

  # Longest common prefix by '_' tokens across a vector of names

  common_stem <- function(nms) {
    if (length(nms) == 0) return("dataset")

    s <- sanitize(nms)
    if (length(s) == 1) return(s)

    split <- strsplit(s, "_", fixed = TRUE)
    min_len <- min(lengths(split))
    stem <- character(0)
    for (k in seq_len(min_len)) {

      tok <- unique(vapply(split, `[`, character(1), k))
      if (length(tok) == 1) stem <- c(stem, tok) else break

    }

    out <- paste(stem, collapse = "_")
    if (!nzchar(out)) "dataset" else out

  }

  # file extension normalize
  ext <- tolower(ext)
  if (!ext %in% c("jpg","jpeg","png","tif","tiff","pdf")) ext <- "jpg"

  # collect outputs
  out <- list(
    SinglePlot = character(0),
    RegionalPlot = character(0),
    ForestPlot = character(0),
    GWASModelPlot = character(0)
  )

  #SinglePlot
  if (!is.null(obj$SinglePlot)) {

    single_root <- file.path(base_dir, "SinglePlot")

    if (!dir.exists(single_root)) dir.create(single_root, recursive = TRUE, showWarnings = FALSE)

    sp <- obj$SinglePlot

    if (!is.list(sp) || inherits(sp, c("ggplot","grob","gTable","gtable","patchwork"))) {

      ds_dir <- file.path(single_root, "_misc")
      if (!dir.exists(ds_dir)) dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

      fn <- file.path(ds_dir, paste0("SinglePlot.", ext))
      h <- get_dyn_h(sp, default = 15)
      save_any(sp, fn, width = single_width, height = h, dpi = dpi)
      out$SinglePlot <- c(out$SinglePlot, fn)

    } else {

      nms <- names(sp)
      for (i in seq_along(sp)) {

        item <- sp[[i]]
        ds_name <- if (!is.null(nms) && nzchar(nms[i])) sanitize(nms[i]) else sprintf("dataset_%02d", i)

        ds_dir <- file.path(single_root, ds_name)

        if (!dir.exists(ds_dir)) dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

        if (inherits(item, c("ggplot","grob","gTable","gtable","patchwork"))) {

          fn <- file.path(ds_dir, paste0(ds_name, ".", ext))
          h <- get_dyn_h(item, default = 15)
          save_any(item, fn, width = single_width, height = h, dpi = dpi)
          out$SinglePlot <- c(out$SinglePlot, fn)

        } else if (is.list(item)) {

          subn <- names(item)
          for (j in seq_along(item)) {

            p <- item[[j]]
            p_name <- if (!is.null(subn) && nzchar(subn[j])) sanitize(subn[j]) else sprintf("%s_%02d", ds_name, j)
            fn <- file.path(ds_dir, paste0(p_name, ".", ext))
            h <- get_dyn_h(p, default = 15)
            save_any(p, fn, width = single_width, height = h, dpi = dpi)
            out$SinglePlot <- c(out$SinglePlot, fn)

          }
        }
      }
    }
  }

  #RegionalPlot (ONE folder per dataset)
  if (is.list(obj$RegionalPlot) && length(obj$RegionalPlot)) {

    reg_root <- file.path(base_dir, "RegionalPlot")

    if (!dir.exists(reg_root)) dir.create(reg_root, recursive = TRUE, showWarnings = FALSE)

    rp <- obj$RegionalPlot
    top_nms <- names(rp)

    #flat list of plots (each element is a plot)

    if (all(vapply(rp, function(x) inherits(x, c("ggplot","grob","gTable","gtable","patchwork")), logical(1)))) {

      ds_name <- common_stem(if (is.null(top_nms)) character(0) else top_nms)
      ds_dir  <- file.path(reg_root, ds_name)

      if (!dir.exists(ds_dir)) dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

      for (j in seq_along(rp)) {

        p <- rp[[j]]
        p_title <- if (!is.null(top_nms) && nzchar(top_nms[j])) top_nms[j] else sprintf("%s_Regional_%02d", ds_name, j)
        fn <- file.path(ds_dir, paste0(sanitize(p_title), ".", ext))
        h <- get_dyn_h(p, default = 30)
        save_any(p, fn, width = regional_width, height = h, dpi = dpi)
        out$RegionalPlot <- c(out$RegionalPlot, fn)

      }

    } else {

      #nested by dataset
      for (i in seq_along(rp)) {

        bucket <- rp[[i]]

        if (inherits(bucket, c("ggplot","grob","gTable","gtable","patchwork"))) {


          nm <- if (!is.null(top_nms) && nzchar(top_nms[i])) top_nms[i] else sprintf("dataset_%02d", i)

          ds_name <- common_stem(nm)
          ds_dir  <- file.path(reg_root, sanitize(ds_name))

          if (!dir.exists(ds_dir)) dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

          fn <- file.path(ds_dir, paste0("Regional_", sanitize(nm), ".", ext))
          h <- get_dyn_h(bucket, default = 30)
          save_any(bucket, fn, width = regional_width, height = h, dpi = dpi)
          out$RegionalPlot <- c(out$RegionalPlot, fn)

          next
        }

        if (!is.list(bucket) || !length(bucket)) next

        inner_nms <- names(bucket)
        ds_from_title <- if (!is.null(inner_nms) && nzchar(inner_nms[1])) dataset_from_title(inner_nms[1]) else NA_character_

        ds_guess <- if (!is.na(ds_from_title)) ds_from_title

        else if (!is.null(top_nms) && nzchar(top_nms[i])) top_nms[i]

        else sprintf("dataset_%02d", i)

        ds_dir <- file.path(reg_root, sanitize(ds_guess))

        if (!dir.exists(ds_dir)) dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

        for (j in seq_along(bucket)) {

          p <- bucket[[j]]
          p_title <- if (!is.null(inner_nms) && nzchar(inner_nms[j])) inner_nms[j] else sprintf("%s_Regional_%02d", ds_guess, j)

          fn <- file.path(ds_dir, paste0(sanitize(p_title), ".", ext))
          h <- get_dyn_h(p, default = 30)
          save_any(p, fn, width = regional_width, height = h, dpi = dpi)
          out$RegionalPlot <- c(out$RegionalPlot, fn)

        }
      }
    }
  }

  #ForestPlot
  if (!is.null(obj$ForestPlot)) {

    forest_root <- file.path(base_dir, "ForestPlot")

    if (!dir.exists(forest_root)) dir.create(forest_root, recursive = TRUE, showWarnings = FALSE)

    fp <- obj$ForestPlot
    fn <- file.path(forest_root, paste0("ForestPlot.", ext))
    h <- get_dyn_h(fp, default = 30)
    save_any(fp, fn, width = forest_width, height = h, dpi = dpi)
    out$ForestPlot <- c(out$ForestPlot, fn)

  }

  #GWASModelPlot
  if (!is.null(obj$GWASModelPlot)) {

    gwas_root <- file.path(base_dir, "GWASModelPlot")

    if (!dir.exists(gwas_root)) dir.create(gwas_root, recursive = TRUE, showWarnings = FALSE)

    gp <- obj$GWASModelPlot
    fn <- file.path(gwas_root, paste0("GWASModelPlot.", ext))
    h <- get_dyn_h(gp, default = 30)
    save_any(gp, fn, width = forest_width, height = h, dpi = dpi)
    out$GWASModelPlot <- c(out$GWASModelPlot, fn)

  }

  return(out)
}
