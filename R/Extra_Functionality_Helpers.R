
  # General helpers for detecting column types by common name

  # Don't need R documentation (Rd) for these helpers

  #' @noRd
  detect_column <- function(Data, Candidate_Column = NULL, Allowed_Names, column_type = "column") {

  # Inherited in the below functions

    .msg <- get0(".detect_msg",
                 ifnotfound = base::message,
                 inherits = TRUE)

    if (is.null(Candidate_Column)) {

      for (name in Allowed_Names) {

        if (name %in% colnames(Data)) {

          .msg(paste0("Using ", name, " as ", column_type, " column (automatic)"))

          return(name)

        }

      }

      .msg(paste("No valid", column_type, "column found."))

      return(NULL)

    } else {

      .msg(paste0("Using ", Candidate_Column, " as ", column_type, " column (manual)"))

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

    allowed <- c("UL", "Upper_CI", "CIU", "CI_Upper")

    detect_column(Data, Upper_CI_Column, allowed, column_type = "upper CI")

  }

  #' @noRd
  detect_lower_ci_column <- function(Data, Lower_CI_Column = NULL) {

    allowed <- c("LL", "Lower_CI", "CIL", "CI_Lower")

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

    allowed <- c(

      # Raw P Variants

      "P","p","Pvalue","pvalue","P-Value","p-value","p-Value","P-VALUE","p_value",

      # log / -log10(P) Variants (detected but NOT transformed) - handled internally in say Single_Plot()

      "logp","LogP","LOGP","Logp","log10p","Log10P","LOG10P","-LOG10P"

    )

    detect_column(Data, PValue_Column, allowed, column_type = "p-value")

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
  run_with_counter <- function(func, args = list(), session = NULL,
                               default_val = NULL,
                               inline_progress = TRUE,
                               capture_output = FALSE,
                               print_run_log = FALSE) {

    exprs <- as.list(body(func))[-1]
    n <- length(exprs)

    if (n == 0L) return(do.call(func, args))

    env <- new.env(parent = environment(func))
    defaults <- formals(func)

    for (nm in names(defaults)) {

      if (!is.null(args[[nm]])) {

        assign(nm, args[[nm]], envir = env)

      } else if (!is.symbol(defaults[[nm]]) && !is.language(defaults[[nm]])) {

        assign(nm, eval(defaults[[nm]], envir = env), envir = env)

      }
    }

    #Progress state

    last_pct <- -1L
    last_msg <- ""
    bar_active <- TRUE
    run_log <- character()

    # Counter helper defined - hidden from global

    log_line <- function(pct_i, msg_i = "") {

      if (!nzchar(msg_i)) return()

      run_log <<- c(run_log, sprintf("[%3d%%] %s", pct_i, msg_i))

    }

    # Counter helper defined - hidden from global

    .draw <- function(pct, msg = "") {

      pct_i <- floor(max(0, min(100, pct)))
      last_pct <<- pct_i; last_msg <<- msg

      # Live Bar

      if (inline_progress && bar_active) {

        bar_len <- 30L
        filled  <- round(bar_len * pct_i / 100)
        cat("\r[", strrep("=", filled), strrep(" ", bar_len - filled), "] ",
            sprintf("%3d%%", pct_i), sep = "")
        flush.console()

      } else if (!inline_progress) {

        message(sprintf("[%3d%%] %s", pct_i, msg))

      }

      # Always buffer a clean transcript line

      log_line(pct_i, msg)

      if (!is.null(session)) {

        session$sendCustomMessage("plot_progress", list(pct = pct_i, msg = msg))

      }

    }

    result <- NULL
    i <- 0L
    on.exit({

      if (inline_progress) cat("\n")

      if (print_run_log && length(run_log)) {

        cat(paste(run_log, collapse = "\n"), "\n", sep = "")

      }
    }, add = TRUE)

    tryCatch({

      for (i in seq_len(n)) {

        assign(".progress", function(frac, msg = "") {

          frac <- max(0, min(1, frac))
          .draw(100 * ((i - 1) + frac) / n, msg)

        }, envir = env)

        assign(".progress_done", function(msg = "") .draw(100 * i / n, msg), envir = env)

        assign(".progress_pause", function() { if (inline_progress && bar_active) cat("\n"); bar_active <<- FALSE }, envir = env)

        assign(".progress_resume", function() { bar_active <<- TRUE; .draw(last_pct, last_msg) }, envir = env)

        .draw(100 * (i - 1) / n, paste("Line", i, "of", n))

        if (i == n) {

          result <- eval(exprs[[i]], envir = env)

        } else {

          if (capture_output) {

            suppressMessages(suppressWarnings(
              capture.output(eval(exprs[[i]], envir = env), type = "output")
            ))

          } else {

            eval(exprs[[i]], envir = env)

          }

        }

        .draw(100 * i / n, paste("Line", i, "done"))
      }

      .draw(100, "Complete")
      invisible(result)

    }, error = function(e) {

      message("run_with_counter(): error in line ", i, " - ", e$message)

      if (!is.null(session)) {

        session$sendCustomMessage("plot_progress", list(pct = NA, msg = paste("Error on line", i, ":", e$message)))

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

    # Choose a graphics device when can't ggsave

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

    # Robust save that tolerates ggplot, patchwork, grob/gtable etc.

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

    # Longest common prefix by '_' across a vector of names

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

    # File extension normalisation

    ext <- tolower(ext)

    if (!ext %in% c("jpg","jpeg","png","tif","tiff","pdf")) ext <- "jpg"

    # Collect outputs

    out <- list(
      SinglePlot = character(0),
      RegionalPlot = character(0),
      ForestPlot = character(0),
      GWASModelPlot = character(0)
    )

    #Single_Plot

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

    #Regional_Plot (ONE! folder per dataset)

    if (is.list(obj$RegionalPlot) && length(obj$RegionalPlot)) {

      reg_root <- file.path(base_dir, "RegionalPlot")

      if (!dir.exists(reg_root)) dir.create(reg_root, recursive = TRUE, showWarnings = FALSE)

      rp <- obj$RegionalPlot
      top_nms <- names(rp)

      #Flat list of plots (each element is a plot)

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

        # Nested by dataset

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

    # Forest_Plot

    if (!is.null(obj$ForestPlot)) {

      forest_root <- file.path(base_dir, "ForestPlot")

      if (!dir.exists(forest_root)) dir.create(forest_root, recursive = TRUE, showWarnings = FALSE)

      fp <- obj$ForestPlot
      fn <- file.path(forest_root, paste0("ForestPlot.", ext))
      h <- get_dyn_h(fp, default = 30)
      save_any(fp, fn, width = forest_width, height = h, dpi = dpi)
      out$ForestPlot <- c(out$ForestPlot, fn)

    }

    # GWASModelPlot of Forest_Plot

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

  #Helpers for inset function (Single_With_Regional_Plot())

  #Transparent saver: make ggplot fully transparent & return grob + aspect ratio

  #' @noRd
  .transparentify_plot_to_grob <- function(p, width_in, height_in, dpi = 100) {

    stopifnot(requireNamespace("patchwork"), requireNamespace("ragg"), requireNamespace("png"))

    p_trans <- p &
      ggplot2::theme(
        plot.background        = ggplot2::element_rect(fill = "transparent", colour = NA),
        panel.background       = ggplot2::element_rect(fill = "transparent", colour = NA),
        legend.background      = ggplot2::element_rect(fill = "transparent", colour = NA),
        legend.box.background  = ggplot2::element_rect(fill = "transparent", colour = NA),
        strip.background       = ggplot2::element_rect(fill = "transparent", colour = NA)
      )

    tmp_png <- tempfile(fileext = ".png")

    ragg::agg_png(
      filename  = tmp_png,
      width     = width_in,
      height    = height_in,
      units     = "in",
      res       = dpi,
      background= "transparent"
    )

    print(p_trans)
    grDevices::dev.off()

    img_arr <- png::readPNG(tmp_png)
    unlink(tmp_png)

    g_inset <- grid::rasterGrob(

      img_arr,
      x = grid::unit(0, "npc"), y = grid::unit(0, "npc"),
      width = grid::unit(1, "npc"), height = grid::unit(1, "npc"),
      just = c("left","bottom"),
      interpolate = TRUE

    )

    ar_img <- ncol(img_arr) / nrow(img_arr)

    list(grob = g_inset, ar = ar_img)
  }

  #' @noRd
  .is_gg <- function(x) inherits(x, "ggplot") || inherits(x, "gg")

  #' @noRd
  .as_single_plot <- function(man_in) {

    if (.is_gg(man_in)) return(man_in)

    if (is.list(man_in)) {

      for (nm in c("plot","p","man","gg","ggplot","g")) {

        if (!is.null(man_in[[nm]]) && .is_gg(man_in[[nm]])) return(man_in[[nm]])

      }

      flat <- unlist(man_in, recursive = FALSE)
      ggs  <- Filter(.is_gg, flat)

      if (length(ggs) >= 1) return(ggs[[1]])

    }

    stop("`Single_Plot` / `LD_Legend_Plot` must be a ggplot or an object/list that contains ggplot.")

  }

  #' @noRd
  .as_reg_list <- function(reg) {

    is_gg <- function(x) inherits(x, "ggplot") || inherits(x, "gg")

    if (is_gg(reg)) return(list(reg))

    if (is.list(reg) && length(reg) > 0 && all(vapply(reg, is_gg, logical(1)))) return(reg)

    if (is.list(reg)) {

      for (nm in c("plots","Plots")) {

        if (!is.null(reg[[nm]]) && is.list(reg[[nm]]) &&

            length(reg[[nm]]) > 0 && all(vapply(reg[[nm]], is_gg, logical(1)))) {

          return(reg[[nm]])

        }
      }

      flat <- unlist(reg, recursive = FALSE)
      pulled <- Filter(is_gg, flat)

      if (length(pulled) > 0) return(pulled)

    }

    stop("`Regional_Plot` must be a ggplot or a list of ggplots")

  }

  #' @noRd
  .get_panel_ranges <- function(p) {

    gb <- ggplot2::ggplot_build(p)
    pp <- gb$layout$panel_params[[1]]

    xr <- if (!is.null(pp$x.range)) pp$x.range else if (!is.null(pp$x) && !is.null(pp$x$range)) pp$x$range else NULL

    yr <- if (!is.null(pp$y.range)) pp$y.range else if (!is.null(pp$y) && !is.null(pp$y$range)) pp$y$range else NULL

    if (is.null(xr) || is.null(yr)) stop("Couldn't extract panel ranges.")

    list(xr = xr, yr = yr)

  }

  #' @noRd
  .get_y_labels_num <- function(p) {

    gb <- ggplot2::ggplot_build(p)
    pp <- gb$layout$panel_params[[1]]

    br <- if (!is.null(pp$y.breaks)) pp$y.breaks else gb$layout$panel_scales_y[[1]]$get_breaks()

    lb <- if (!is.null(pp$y.labels)) pp$y.labels else gb$layout$panel_scales_y[[1]]$get_labels(br)

    lb_num <- suppressWarnings(as.numeric(lb))

    if (length(lb) != length(lb_num) || any(!is.finite(lb_num))) {

      stop("Y-axis labels are not purely numeric; cannot use labels to define BIG mapping.")

    }

    lb_num

  }

  #' @noRd
  .pick_y <- function(df) {

    if ("log10p"    %in% names(df)) return(as.numeric(df$log10p))

    if ("log10pval" %in% names(df)) return(as.numeric(df$log10pval))

    if ("P"         %in% names(df)) return(-log10(as.numeric(df$P)))

    stop("Need y as 'log10p' or 'log10pval' (or 'P' to compute).")

  }

  #' @noRd
  .to_big_with_label <- function(y, max_label) {

    y  <- as.numeric(y)
    ml <- as.numeric(max_label)[1]

    if (!is.finite(ml) || ml <= 0) {

      return(list(big = rep(0, length(y)), segment = rep("invalid", length(y))))

    }

    out <- numeric(length(y)); seg <- character(length(y))

    if (ml > 10) {

      low  <- y <= 10
      high <- y > 10
      out[low]  <- 8 * pmax(0, pmin(10, y[low]))
      seg[low]  <- "low(0-10)"
      span_hi   <- (ml - 10)
      slope_hi  <- 20 / span_hi
      y_clamped <- pmax(10, pmin(ml, y[high]))
      out[high] <- 80 + (y_clamped - 10) * slope_hi
      seg[high] <- sprintf("high(10-%.6g), slope=%.6f", ml, slope_hi)

    } else {

      slope_all <- 100 / ml
      out <- pmax(0, pmin(100, y * slope_all))
      seg[] <- sprintf("all(0-%.6g), slope=%.6f", ml, slope_all)

    }

    list(big = pmax(0, pmin(100, out)), segment = seg)

  }

  #' @noRd
  .to_big_one <- function(val, max_label) {

    as.numeric(.to_big_with_label(val, max_label)$big)[1]

  }

  #' @noRd
  .big_to_raw_with_label <- function(big, max_label) {

    b  <- as.numeric(big)
    ml <- as.numeric(max_label)[1]

    if (!is.finite(ml) || ml <= 0) return(rep(0, length(b)))

    if (ml > 10) {

      ifelse(b <= 80, b/8, 10 + (b - 80) * (ml - 10) / 20)

    } else {

      b * (ml / 100)

    }

  }


  #' @noRd
  .pt_to_mm <- function(pt) pt / 2.83465

  #' @noRd
  .per_inset_num <- function(val, idx, chr = NULL, default = 0) {

    if (is.null(val)) return(default)

    if (is.character(val) && length(val) == 1L && grepl(",", val)) {

      val <- trimws(strsplit(val, ",")[[1]])

    }

    if (is.list(val)) val <- unlist(val, use.names = TRUE)

    nm <- names(val)

    if (!is.null(chr) && !is.null(nm) && any(nzchar(nm))) {

      chr_str <- as.character(chr)

      if (chr_str %in% nm) {

        v <- suppressWarnings(as.numeric(val[[chr_str]]))

        return(ifelse(is.finite(v), v, default))

      }

    }

    if (length(val) == 1L) {

      v <- suppressWarnings(as.numeric(val))

      return(ifelse(is.finite(v), v, default))

    }

    j <- max(1, min(idx, length(val)))
    v <- suppressWarnings(as.numeric(val[[j]]))

    ifelse(is.finite(v), v, default)

  }

  #' @noRd
  .per_inset_val <- function(val, idx, chr = NULL, default = NULL) {

    if (is.null(val)) return(default)

    if (is.list(val)) val <- unlist(val, use.names = TRUE)

    if (is.character(val) && length(val) == 1L && grepl(",", val)) {

      val <- trimws(strsplit(val, ",")[[1]])

    }

    if (length(val) == 1L && is.null(names(val))) {

      v <- suppressWarnings(as.numeric(val))

      return(ifelse(is.finite(v), v, default))

    }

    nm <- names(val)

    if (!is.null(chr) && !is.null(nm) && any(nzchar(nm))) {

      chr_str <- as.character(chr)

      if (chr_str %in% nm) {

        v <- suppressWarnings(as.numeric(val[[chr_str]]))

        return(ifelse(is.finite(v), v, default))

      }

    }

    j <- max(1, min(idx, length(val)))
    v <- suppressWarnings(as.numeric(val[[j]]))
    ifelse(is.finite(v), v, default)

  }

  #' @noRd
  .parse_preset <- function(x) {

    if (is.null(x)) return(NULL)

    if (is.character(x) && length(x) == 1L && grepl(",", x)) {

      return(trimws(strsplit(x, ",")[[1]]))

    }

    x

  }

  #' @noRd
  .per_inset_preset <- function(preset, idx, chr = NULL, default = "Right") {

    if (is.null(preset)) return(default)

    p <- .parse_preset(preset)

    if (is.character(p) && !is.null(names(p)) && any(nzchar(names(p)))) {

      chr_str <- as.character(chr)

      if (chr_str %in% names(p)) return(p[[chr_str]])

    }

    if (length(p) == 1L) return(p[[1]])

    j <- ((idx - 1L) %% length(p)) + 1L

    p[[j]]

  }

  # Null-coalescing helper

  #' @noRd
  `%||%` <- function(a, b) if (!is.null(a)) a else b


  #' @noRd
  make_ld_legend_plot_from_spec <- function(ld_info) {

    if (is.null(ld_info)) return(NULL)

    stopifnot(requireNamespace("ggplot2"), requireNamespace("cowplot"))

    title_size <- ld_info$title_size %||% 30
    text_size  <- ld_info$text_size  %||% 20

    legend_margin <- ld_info$legend_margin %||% ggplot2::margin(0, 0, 0, 0)
    box_margin    <- ld_info$box_margin    %||% ggplot2::margin(0, 0, 0, 0)

    message("Building LD legend from LD_legend_spec:")
    message("  title:      ", ld_info$title)
    message("  levels:     ", paste(ld_info$levels, collapse = ", "))
    message("  title_size: ", title_size,

            if (!is.null(ld_info$title_size)) " (manual)" else " (default)")

    message("  text_size:  ", text_size,

            if (!is.null(ld_info$text_size))  " (manual)" else " (default)")

    df_legend <- data.frame(
      LD_Bin = factor(ld_info$levels, levels = ld_info$levels),
      x = 1,
      y = seq_along(ld_info$levels)
    )

    p_legend <- ggplot2::ggplot(df_legend,
                                ggplot2::aes(x = x, y = y, colour = LD_Bin)) +

      # no visible points in panel; only used to create legend entries

      ggplot2::geom_point(size = 0, alpha = 0, show.legend = TRUE) +
      ggplot2::scale_colour_manual(
        values = ld_info$legend_colors,
        drop   = FALSE
      ) +
      ggplot2::labs(colour = ld_info$title, x = NULL, y = NULL) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position      = "right",
        legend.justification = "top",
        legend.key.height    = grid::unit(1, "cm"),
        legend.title         = ggplot2::element_text(size = title_size, face = "bold"),
        legend.text          = ggplot2::element_text(size = text_size),
        legend.background    = ggplot2::element_rect(
          fill      = scales::alpha("white", 0.7),
          colour    = "black",
          linewidth = 0.5
        ),
        legend.margin        = legend_margin,
        legend.box.margin    = box_margin,
        legend.spacing.x     = grid::unit(0, "pt"),
        legend.spacing.y     = grid::unit(0, "pt"),
        plot.margin          = ggplot2::margin(0, 0, 0, 0)
      ) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          override.aes = list(
            shape  = unname(ld_info$shapes[ld_info$levels]),
            size   = 7,
            alpha  = 1,
            colour = unname(ld_info$legend_colors[ld_info$levels])
          )

        )

      )

    p_legend

  }

#' @noRd
make_gene_legend_plot_from_spec <- function(gene_info) {

  if (is.null(gene_info)) return(NULL)

  stopifnot(requireNamespace("ggplot2"),
            requireNamespace("cowplot"),
            requireNamespace("scales"))

  title_size <- gene_info$title_size %||% 30
  text_size  <- gene_info$text_size  %||% 20

  legend_margin <- gene_info$legend_margin %||% ggplot2::margin(0, 0, 0, 0)
  box_margin    <- gene_info$box_margin    %||% ggplot2::margin(0, 0, 0, 0)

  message("Building Gene legend from Gene_legend_spec:")
  message("  title:      ", gene_info$title)
  message("  levels:     ", paste(gene_info$levels, collapse = ", "))
  message("  title_size: ", title_size,

          if (!is.null(gene_info$title_size)) " (manual)" else " (default)")

  message("  text_size:  ", text_size,

          if (!is.null(gene_info$text_size))  " (manual)" else " (default)")

  n_levels <- length(gene_info$levels %||% character(0))

  col_vec <- gene_info$colors %||%
    gene_info$colours %||%
    gene_info$legend_colors

  if (!is.null(col_vec)) {

    col_vec <- unname(as.vector(col_vec))

  }

  if (is.null(col_vec) || length(col_vec) < n_levels) {

    message("Gene legend colours missing or too short; generating default palette via scales::hue_pal().")

    col_vec <- scales::hue_pal()(n_levels)

  }

  if (length(col_vec) > n_levels) {

    col_vec <- col_vec[seq_len(n_levels)]

  }

  df_legend <- data.frame(

    Gene_Type = factor(gene_info$levels, levels = gene_info$levels),
    x = 1,
    y = seq_along(gene_info$levels)
  )

  p_legend <- ggplot2::ggplot(df_legend,
                              ggplot2::aes(x = x, y = y, colour = Gene_Type)) +

    # invisible in panel; just defines legend entries

    ggplot2::geom_point(size = 0, alpha = 0, show.legend = TRUE) +
    ggplot2::scale_colour_manual(
      values = col_vec,
      drop   = FALSE
    ) +
    ggplot2::labs(colour = gene_info$title, x = NULL, y = NULL) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position      = "right",
      legend.justification = "top",
      legend.key.height    = grid::unit(1, "cm"),
      legend.title         = ggplot2::element_text(size = title_size, face = "bold"),
      legend.text          = ggplot2::element_text(size = text_size),
      legend.background    = ggplot2::element_rect(
        fill      = scales::alpha("white", 0.7),
        colour    = "black",
        linewidth = 0.5
      ),
      legend.margin        = legend_margin,
      legend.box.margin    = box_margin,
      legend.spacing.x     = grid::unit(0, "pt"),
      legend.spacing.y     = grid::unit(0, "pt"),
      plot.margin          = ggplot2::margin(0, 0, 0, 0)
    ) +
    ggplot2::guides(

      colour = ggplot2::guide_legend(

        override.aes = list(

          # square colour blocks in the key

          shape  = 15,
          size   = 7,
          alpha  = 1,
          colour = col_vec

        )

      )

    )

  p_legend

}


  #' @noRd
  if (!exists("run_with_counter", mode = "function")) {

    run_with_counter <- function(fn, args = list(), session = NULL) {

      do.call(fn, args)

    }

  }

  #' @noRd
  .silence_messages <- function(expr) {

    withCallingHandlers(expr, message = function(m) invokeRestart("muffleMessage"))

  }

  #' @noRd
  .silence_warnings <- function(expr) suppressWarnings(expr)


  #' @noRd
  .collect_gene_specs <- function(regs_list, Gene_Legend_Plot = NULL, include_explicit = TRUE) {

    specs <- list()

    # Optionally treat explicit Gene_Legend_Plot as just another source

    if (include_explicit && !is.null(Gene_Legend_Plot)) {

      gp <- .as_single_plot(Gene_Legend_Plot)
      gs <- .get_attr(gp, "Gene_legend_spec")

      if (!is.null(gs)) {

        specs[[length(specs)+1]] <- list(name = "Gene_Legend_Plot", spec = gs)

      }

    }

    for (i in seq_along(regs_list)) {

      g  <- .unwrap_plot(regs_list[[i]])
      gs <- .get_attr(g, "Gene_legend_spec")

      if (!is.null(gs)) {

        nm <- attr(g, "plot_title", exact = TRUE) %||%
          paste0("Regional_", i)
        specs[[length(specs)+1]] <- list(name = nm, spec = gs)

      }

    }

    specs

  }

  #' @noRd
  .merge_gene_specs <- function(spec_entries) {

    if (length(spec_entries) == 0) return(list(spec = NULL, provenance = NULL))

    # Union of all levels in appearance order

    levels_all <- character(0)
    push_levels <- function(x) {

      for (lv in x) if (!(lv %in% levels_all)) levels_all <<- c(levels_all, lv)

    }

    for (e in spec_entries) {

      lv <- e$spec$levels

      if (!is.null(lv)) push_levels(lv)

    }

    col_map <- list()
    shp_map <- list()
    prov    <- list()

    record <- function(level, field, from_name, from_idx, value) {

      key <- paste(level, field, sep = "::")
      prov[[key]] <<- data.frame(
        Level           = level,
        Field           = field,
        From_Spec_Name  = from_name,
        From_Spec_Index = from_idx,
        Value           = as.character(value),
        stringsAsFactors = FALSE
      )

    }

    for (i in seq_along(spec_entries)) {

      nm <- spec_entries[[i]]$name
      gs <- spec_entries[[i]]$spec

      cols <- gs$legend_colors %||% gs$colors %||% gs$colours
      shps <- gs$shapes
      levs <- gs$levels %||% names(cols) %||% names(shps)

      if (!is.null(cols)) {

        for (lv in intersect(levs, names(cols))) {

          if (is.null(col_map[[lv]]) && !is.null(cols[[lv]])) {

            col_map[[lv]] <- cols[[lv]]
            record(lv, "colour", nm, i, cols[[lv]])

          }

        }

      }

      if (!is.null(shps)) {

        for (lv in intersect(levs, names(shps))) {

          if (is.null(shp_map[[lv]]) && !is.null(shps[[lv]])) {

            shp_map[[lv]] <- shps[[lv]]
            record(lv, "shape", nm, i, shps[[lv]])

          }

        }

      }

    }

    first_spec <- spec_entries[[1]]$spec

    legend_colors <- if (length(col_map)) unlist(col_map[levels_all]) else NULL

    if (!is.null(legend_colors)) names(legend_colors) <- levels_all

    shapes <- if (length(shp_map)) unlist(shp_map[levels_all]) else NULL

    if (!is.null(shapes)) names(shapes) <- levels_all

    merged <- list(
      levels        = levels_all,
      legend_colors = legend_colors,
      shapes        = shapes,
      title         = first_spec$title      %||% "Gene Biotype",
      title_size    = first_spec$title_size %||% 30,
      text_size     = first_spec$text_size  %||% 20,
      legend_margin = first_spec$legend_margin %||% ggplot2::margin(0,0,0,0),
      box_margin    = first_spec$box_margin    %||% ggplot2::margin(0,0,0,0)
    )

    prov_tbl <- if (length(prov)) do.call(rbind, prov) else NULL

    list(spec = merged, provenance = prov_tbl)

  }

  #' @noRd
  .unwrap_plot <- function(x) {

    # If it's already a ggplot, return as-is

    if (inherits(x, "ggplot")) return(x)

    # If it's a one-item list that contains a ggplot, unwrap it

    if (is.list(x) && length(x) == 1L && inherits(x[[1]], "ggplot")) return(x[[1]])

    x

  }

  #' @noRd
  .get_attr <- function(obj, nm) {

    # Try attribute on the object

    a <- attr(obj, nm, exact = TRUE)

    if (!is.null(a)) return(a)

    # If object is a one-item list, try its child

    if (is.list(obj) && length(obj) == 1L) return(attr(obj[[1]], nm, exact = TRUE))

    NULL

  }

  #' @noRd
  .legend_to_grob_fixed <- function(make_plot_fn,
                                    spec,
                                    scale_val,
                                    xr, yr,
                                    Target_Width,
                                    Target_Height,
                                    base_height_frac = 0.22,
                                    dpi = 100,
                                    side = c("right", "left"),
                                    vert = c("top", "bottom"),
                                    Transparent_Legends = FALSE,
                                    trim_fuzz = "6%",
                                    inner_pad_pt = 6,
                                    outer_pad_pt = 0) {

    side <- match.arg(side)
    vert <- match.arg(vert)

    if (is.null(spec)) return(NULL)

    p_leg <- make_plot_fn(spec)

    if (is.null(p_leg)) return(NULL)

    title_sz <- (spec$title_size %||% 30)
    text_sz  <- (spec$text_size  %||% 20)

    p_leg <- p_leg +
      ggplot2::theme(
        legend.title = ggplot2::element_text(size = title_sz),
        legend.text  = ggplot2::element_text(size = text_sz)
      )

    if (isTRUE(Transparent_Legends)) {

      p_leg <- p_leg +
        ggplot2::theme(
          legend.background     = ggplot2::element_rect(fill = NA, colour = NA),
          legend.box.background = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.5)
        )

    } else {

      p_leg <- p_leg +
        ggplot2::theme(
          legend.background     = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.5),
          legend.box.background = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.5)
        )

    }

    message("Padding and spacing the legends")

    p_leg <- p_leg +
      ggplot2::theme(

        legend.margin        = ggplot2::margin(inner_pad_pt, inner_pad_pt, inner_pad_pt, inner_pad_pt),
        legend.box.margin    = ggplot2::margin(outer_pad_pt, outer_pad_pt, outer_pad_pt, outer_pad_pt),
        legend.spacing.y     = grid::unit(2, "pt"),
        legend.spacing.x     = grid::unit(2, "pt"),
        plot.margin          = ggplot2::margin(0, 0, 0, 0),

        legend.position      = if (side == "left") "left" else "right",

        legend.justification = c(if (side == "left") 0 else 1, if (vert == "top") 1 else 0)

      )

    sc <- suppressWarnings(as.numeric(scale_val[1])); if (!is.finite(sc)) sc <- 1
    sc <- max(0.1, min(sc, 3))

    dx <- diff(xr); dy <- diff(yr)
    ar_panel <- Target_Width / Target_Height

    # Render at natural base size then trim

    base_h_in <- Target_Height * base_height_frac
    base_w_in <- base_h_in * 0.8

    tf <- tempfile(fileext = ".png")
    ggplot2::ggsave(tf, p_leg, width = base_w_in, height = base_h_in, dpi = dpi, limitsize = FALSE)
    img <- magick::image_read(tf); unlink(tf)

    # Trim away padding outside the border

    img <- magick::image_trim(img, fuzz = trim_fuzz)

    # Post-trim pixel info

    info_trim <- magick::image_info(img)
    ar_img    <- info_trim$width / info_trim$height

    # Scale data-height by post-trim pixel height

    base_h_px <- base_h_in * dpi
    base_hdat <- base_height_frac * dy
    hdat0     <- base_hdat * sc
    hdat      <- hdat0 * (info_trim$height / base_h_px)

    # Width via aspect conversions

    wdat <- hdat * ar_img * (1 / ar_panel) * (dx / dy)

    # Clamp, if needed...

    if (hdat > dy) {

      hdat <- dy * 0.9
      wdat <- hdat * ar_img * (1 / ar_panel) * (dx / dy)

    }

    # Grob for annotation_custom

    grob <- grid::rasterGrob(
      as.raster(img),
      x           = grid::unit(0, "npc"),
      y           = grid::unit(0, "npc"),
      width       = grid::unit(1, "npc"),
      height      = grid::unit(1, "npc"),
      just        = c("left","bottom"),
      interpolate = TRUE
    )

    list(grob = grob, wdat = wdat, hdat = hdat)

  }

  #End of inset helper functions - could slim down in future.

  # SAVE helper - HLA

  Save_Segregate_HLA_Plots <- function(res, dir = ".", dpi = 150, units = "in") {

    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

    save_one <- function(p) {

      fn <- attr(p, "suggested_filename")

      w  <- attr(p, "recommended_width")

      h  <- attr(p, "recommended_height")

      if (is.null(fn) || is.null(w) || is.null(h)) return(invisible(NULL))

      ggplot2::ggsave(

        filename  = file.path(dir, fn),
        plot      = p,
        width     = w,
        height    = h,
        units     = units,
        dpi       = dpi,
        limitsize = FALSE

      )

      invisible(fn)

    }

    if (!is.null(res$HLA_plot))  save_one(res$HLA_plot)

    if (!is.null(res$AA_plot))   save_one(res$AA_plot)

    if (!is.null(res$SNPS_plot)) save_one(res$SNPS_plot)

    if (!is.null(res$RSID_plots) && length(res$RSID_plots)) {

      for (nm in names(res$RSID_plots)) save_one(res$RSID_plots[[nm]])

    }

    invisible(TRUE)
  }

