
  # options(shiny.useragg = FALSE)   # force Shiny NOT to use ragg
  # options(bitmapType = "cairo")    # more reliable headless rendering
  # source("renv/activate.R")


  # rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)

   use_wrapper <- FALSE

  # devtools::load_all(".")

  suppressPackageStartupMessages({

    library(MiamiR)
    library(shiny)
    library(grid)
    library(ggplot2)
    library(ggplotify)
    library(vroom)
    library(slickR)
    library(ggiraph)
    library(rlang)
    library(tools)
    library(dplyr)
    library(stringr)
    library(stringi)
    library(ggtext)
    library(ggh4x)
    library(patchwork)
    library(viridis)
    library(gtable)
    library(plotly)
    library(htmlwidgets)
    library(ggrastr)
    library(cowplot)

  })


  options(shiny.maxRequestSize = -1)

#Load functions from package

if (requireNamespace("MiamiR", quietly = TRUE)) {

  Miami_Plot        <- MiamiR::Miami_Plot
  Single_Plot       <- MiamiR::Single_Plot
  Regional_Plot     <- MiamiR::Regional_Plot
  Annotate_Data     <- MiamiR::Annotate_Data
  METASOFT_File_Gen <- MiamiR::METASOFT_File_Gen
  Forest_Plot       <- MiamiR::Forest_Plot
  Model_Munge       <- MiamiR::Model_Munge

  .Single_Plot_original       <- tryCatch(getFromNamespace(".Single_Plot_original", "MiamiR"),

                                          error = function(e) MiamiR::Single_Plot)

  .Regional_Plot_original     <- tryCatch(getFromNamespace(".Regional_Plot_original", "MiamiR"),

                                          error = function(e) MiamiR::Regional_Plot)

  .METASOFT_File_Gen_original <- tryCatch(getFromNamespace("METASOFT_File_Gen", "MiamiR"),

                                          error = function(e) MiamiR::METASOFT_File_Gen)

  .Forest_Plot_original       <- tryCatch(getFromNamespace(".Forest_Plot_original", "MiamiR"),

                                          error = function(e) MiamiR::Forest_Plot)

  .Model_Munge_original       <- tryCatch(getFromNamespace("Model_Munge", "MiamiR"),

                                          error = function(e) MiamiR::Model_Munge)

  Segregate_HLA_Plot <- MiamiR::Segregate_HLA_Plot

  .Segregate_HLA_Plot_original <- tryCatch(

    getFromNamespace(".Segregate_HLA_Plot_original", "MiamiR"),

    error = function(e) MiamiR::Segregate_HLA_Plot

  )

} else {

  if (file.exists("R/Segregate_HLA_Plot.R")) source("R/Segregate_HLA_Plot.R")

  stopifnot(exists("Segregate_HLA_Plot", mode="function"))

  .Segregate_HLA_Plot_original <- if (exists(".Segregate_HLA_Plot_original")) .Segregate_HLA_Plot_original else Segregate_HLA_Plot

  if (file.exists("R/Miami_Plot.R"))        source("R/Miami_Plot.R")
  if (file.exists("R/Single_Plot.R"))       source("R/Single_Plot.R")
  if (file.exists("R/Regional_Plot.R"))     source("R/Regional_Plot.R")
  if (file.exists("R/Annotate_Data.R"))     source("R/Annotate_Data.R")
  if (file.exists("R/METASOFT_File_Gen.R")) source("R/METASOFT_File_Gen.R")
  if (file.exists("R/Forest_Plot.R"))       source("R/Forest_Plot.R")
  if (file.exists("R/Model_Munge.R"))       source("R/Model_Munge.R")

  stopifnot(

    exists("Single_Plot", mode="function"),
    exists("Regional_Plot", mode="function"),
    exists("Annotate_Data", mode="function"),
    exists("METASOFT_File_Gen", mode="function"),
    exists("Forest_Plot", mode="function"),
    exists("Model_Munge", mode="function")

  )

  .Single_Plot_original       <- if (exists(".Single_Plot_original")) .Single_Plot_original else Single_Plot

  .Regional_Plot_original     <- if (exists(".Regional_Plot_original")) .Regional_Plot_original else Regional_Plot

  .METASOFT_File_Gen_original <- METASOFT_File_Gen

  .Forest_Plot_original       <- if (exists(".Forest_Plot_original")) .Forest_Plot_original else Forest_Plot

  .Model_Munge_original       <- Model_Munge

}

  # Shared helpers

  # Regional interactive helpers (prefixed)

  `%||%` <- function(a,b) if (!is.null(a)) a else b

  reg_BASE_SVG_WIDTH_IN <- 30
  reg_DISPLAY_WIDTH_PX  <- 1200
  reg_LEFT_PAD_PT       <- 80
  reg_CAPTURE_SIZE_TOP    <- 18
  reg_CAPTURE_SIZE_BOTTOM <- 14

 reg_split_regional <- function(x) {

  if (is.list(x) && all(c("object","object2") %in% names(x))) return(list(top = x$object, bottom = x$object2))

  if (inherits(x, "patchwork")) {

    pls <- NULL

    if (!is.null(x$plots)) pls <- x$plots

    if (is.null(pls) && !is.null(x$patches) && !is.null(x$patches$plots)) pls <- x$patches$plots

    if (is.null(pls)) pls <- tryCatch(patchwork:::plots(x), error = function(e) NULL)

    if (!is.null(pls) && length(pls) >= 2) return(list(top = pls[[1]], bottom = pls[[2]]))

  }

  if (is.list(x) && length(x) == 2) return(list(top = x[[1]], bottom = x[[2]]))

  stop("Couldn't split: unknown structure")

 }

 reg_find_xy_for_overlay <- function(p) {

  pick <- function(nm, map, dat) {

    v <- if (!is.null(map[[nm]])) rlang::as_label(map[[nm]]) else NULL

    if (!is.null(v) && v %in% names(dat)) v else NULL

  }

  for (lyr in p$layers) {

    dat <- if (!is.null(lyr$data) && NROW(lyr$data)) lyr$data else p$data

    if (is.null(dat) || !NROW(dat)) next

    map <- modifyList(as.list(p$mapping), as.list(lyr$mapping))
    x <- pick("x", map, dat); y <- pick("y", map, dat)

    if (!is.null(x) && !is.null(y)) return(list(df = dat, x = x, y = y, xlab = x, ylab = y))

    x1 <- pick("x", map, dat); x2 <- pick("xend", map, dat)
    y1 <- pick("y", map, dat); y2 <- pick("yend", map, dat)

    if (!is.null(x1) && !is.null(x2)) {

      dat$.__xmid <- (dat[[x1]] + dat[[x2]])/2

      if (!is.null(y1) && !is.null(y2)) {

        dat$.__ymid <- (dat[[y1]] + dat[[y2]])/2

        return(

          list(

            df = dat, x = ".__xmid", y = ".__ymid",
                    xlab = paste0(x1,"..",x2), ylab = paste0(y1,"..",y2)

            )

          )

      }

      if (!is.null(y1)) {

        dat$.__ymid <- dat[[y1]]

        return(

          list(

            df = dat, x = ".__xmid", y = ".__ymid",
                    xlab = paste0(x1,"..",x2), ylab = y1)

               )

        }

    }

    xmin <- pick("xmin", map, dat); xmax <- pick("xmax", map, dat)
    ymin <- pick("ymin", map, dat); ymax <- pick("ymax", map, dat)

    if (!is.null(xmin) && !is.null(xmax)) {

      dat$.__xmid <- (dat[[xmin]] + dat[[xmax]])/2

      if (!is.null(ymin) && !is.null(ymax)) {

        dat$.__ymid <- (dat[[ymin]] + dat[[ymax]])/2

        return(list(df = dat, x = ".__xmid", y = ".__ymid",
                    xlab = paste0(xmin,"..",xmax), ylab = paste0(ymin,"..",ymax)))

      }

    }

  }

  NULL

}

 reg_build_tooltips <- function(df, xlab, ylab, xcol = NULL, ycol = NULL) {

  if ("Hover_Info" %in% names(df)) return(as.character(df$Hover_Info))

  pick <- function(nm) if (nm %in% names(df)) nm else NULL

  idcol  <- pick("SNP") %||% pick("ID") %||% pick("RS") %||% pick("rsid")

  xnm <- if (!is.null(xcol) && xcol %in% names(df)) xcol else if (!is.null(xlab) && xlab %in% names(df)) xlab else NULL

  ynm <- if (!is.null(ycol) && ycol %in% names(df)) ycol else if (!is.null(ylab) && ylab %in% names(df)) ylab else NULL

  format_field <- function(cl, v) {

    if (length(v) == 0 || is.na(v)) return("")

    if (!is.null(cl) && cl %in% c("P", "LOG10P")) {

      num <- suppressWarnings(as.numeric(v))

      if (is.finite(num)) return(format(signif(num, 3), trim = TRUE, scientific = TRUE))

      return(as.character(v))

    }

    if (is.numeric(v)) {

      if (is.finite(v) && v == floor(v)) return(formatC(v, format = "d", big.mark = ","))

      return(format(signif(v, 6), trim = TRUE, scientific = FALSE, big.mark = ","))

    }

    as.character(v)

  }

  parts <- list(

    if (!is.null(idcol)) "ID: %s" else NULL,

    if (!is.null(pick("CHR"))) "CHR: %s" else NULL,

    if (!is.null(pick("POS"))) "POS: %s" else NULL,

    if (!is.null(pick("P")))   "P: %s" else NULL,

    if (!is.null(pick("LOG10P"))) "LOG10P: %s" else NULL,

    if (!is.null(pick("REF"))) "REF: %s" else NULL,

    if (!is.null(pick("ALT"))) "ALT: %s" else NULL,

    if (!is.null(xnm)) paste0((xlab %||% xnm), ": %s") else NULL,

    if (!is.null(ynm)) paste0((ylab %||% ynm), ": %s") else NULL

  )

  parts <- parts[!vapply(parts, is.null, logical(1))]

  if (!length(parts)) return(rep("point", NROW(df)))

  fmt  <- paste(parts, collapse = "<br/>")

  cols <- c(

    if (!is.null(idcol)) idcol,

    intersect(c("CHR","POS","P","LOG10P","REF","ALT"), names(df)),

    c(xnm, ynm)

  )

  cols <- cols[!duplicated(cols) & cols %in% names(df)]

  vapply(seq_len(NROW(df)), function(i) {

    vals <- vapply(cols, function(cl) format_field(cl, df[[cl]][i]), character(1))

    do.call(sprintf, c(fmt, as.list(vals)))

  }, character(1))

}

 reg_add_big_capture_layer <- function(p, size_pts = 60, tag = "panel") {

  info <- reg_find_xy_for_overlay(p)

  if (is.null(info)) return(p)

  df <- info$df
  xlab <- p$labels$x %||% info$xlab
  ylab <- p$labels$y %||% info$ylab
  df$.tooltip <- reg_build_tooltips(df, xlab, ylab, xcol = info$x, ycol = info$y)

  idcol <- if ("ID" %in% names(df)) "ID" else if ("SNP" %in% names(df)) "SNP" else if ("RS" %in% names(df)) "RS" else if ("rsid" %in% names(df)) "rsid" else NULL

  if (!is.null(idcol)) df$.id <- as.character(df[[idcol]]) else df$.id <- paste0(tag, "-", seq_len(NROW(df)))

  p + ggiraph::geom_point_interactive(
    data = df,
    mapping = ggplot2::aes(
      x = .data[[info$x]], y = .data[[info$y]],
      tooltip = .data$.tooltip, data_id = .data$.id
    ),
    size  = size_pts, alpha = 0.002, stroke = 0,
    inherit.aes = FALSE, show.legend = FALSE
  )

}

reg_get_xlim_from_plot <- function(p) {

  if (!inherits(p, "ggplot")) return(NULL)

  b <- ggplot2::ggplot_build(p)
  rng <- NULL

  pp <- try(b$layout$panel_params[[1]], silent = TRUE)

  if (!inherits(pp, "try-error")) {

    for (cand in list(try(pp$x.range, silent=TRUE),
                      try(pp$x$range$range, silent=TRUE),
                      try(pp$x$range$rng, silent=TRUE))) {

      if (!inherits(cand, "try-error") && length(cand)==2 && all(is.finite(cand))) {

        rng <- as.numeric(cand); break
      }

    }

  }

  if (is.null(rng)) {

    xs <- numeric()

    for (lyr in p$layers) {

      dat <- if (!is.null(lyr$data) && NROW(lyr$data)) lyr$data else p$data

      if (is.null(dat) || !NROW(dat)) next

      map <- modifyList(as.list(p$mapping), as.list(lyr$mapping))

      pick <- function(nm) {

        v <- if (!is.null(map[[nm]])) rlang::as_label(map[[nm]]) else NULL

        if (!is.null(v) && v %in% names(dat)) dat[[v]] else NULL

      }

      for (nm in c("x","xmin","xmax","xend")) {

        v <- pick(nm); if (!is.null(v) && is.numeric(v)) xs <- c(xs, v)

      }

    }

    if (length(xs)) rng <- range(xs, finite = TRUE)

  }

  rng

}

reg_overlay_fullpanel <- function(obj, xlim = NULL, size_pts = 60, tag = "panel") {

  if (is.null(xlim) || length(xlim) != 2 || any(!is.finite(xlim))) xlim <- c(0, 1)

  overlay_df <- data.frame(x = seq(xlim[1], xlim[2], length.out = 400), y = 0.5)
  overlay_df$.tooltip <- sprintf("x: %s", signif(overlay_df$x, 6))
  overlay_df$.id <- paste0(tag, "-grob-", seq_len(nrow(overlay_df)))
  overlay_plot <- ggplot2::ggplot() +
    ggiraph::geom_point_interactive(
      data = overlay_df,
      ggplot2::aes(x = x, y = y, tooltip = .tooltip, data_id = .id),
      size = size_pts, alpha = 0.002, stroke = 0, show.legend = FALSE
    ) +
    ggplot2::scale_x_continuous(limits = xlim, expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_y_continuous(limits = c(0,1), expand = ggplot2::expansion(mult = 0)) +
    ggplot2::theme_void()
  base <- if (inherits(obj, "ggplot")) obj else ggplotify::as.ggplot(obj)
  base + patchwork::inset_element(overlay_plot, left = 0, bottom = 0, right = 1, top = 1, align_to = "panel")

}

reg_ensure_hover <- function(p, size_pts = 60, tag = "panel") {

  if (inherits(p, "ggplot")) {

    cand <- reg_add_big_capture_layer(p, size_pts = size_pts, tag = tag)

    if (identical(length(cand$layers), length(p$layers))) {

      return(reg_overlay_fullpanel(p, xlim = reg_get_xlim_from_plot(p), size_pts = size_pts, tag = tag))

    }

    return(cand)

  }

  reg_overlay_fullpanel(p, xlim = c(0,1), size_pts = size_pts, tag = tag)

}

reg_extract_rel_heights <- function(pw) {

  if (!inherits(pw, "patchwork")) return(NULL)

  h <- NULL
  h <- h %||% pw$patches$layout$heights
  h <- h %||% pw$heights

  if (!is.null(h)) {

    tryCatch({
      v <- as.numeric(h)
      if (all(is.finite(v)) && length(v) >= 2) return(v[1:2])
      NULL
    }, error = function(e) NULL)

  } else NULL

}

reg_inject_gene_hover <- function(p, capture_size = reg_CAPTURE_SIZE_BOTTOM) {

  if (!inherits(p, "ggplot")) return(p)

  layer_dfs <- list()

  for (idx in seq_along(p$layers)) {

    lyr <- p$layers[[idx]]

    dat <- if (!is.null(lyr$data) && NROW(lyr$data)) lyr$data else p$data

    if (is.data.frame(dat) && all(c("start","end","y") %in% names(dat))) {

      feat <- {

        v <- rep(NA_character_, nrow(dat))
        txt_cols <- intersect(c("feature","type","class","label","region","annotation"), names(dat))

        if (length(txt_cols)) {

          txt <- tolower(trimws(apply(dat[, txt_cols, drop=FALSE], 1, paste, collapse=" ")))
          v[grepl("\\bintron\\b", txt)] <- "INTRON"
          v[is.na(v) & grepl("\\bexon\\b|\\butr\\b", txt)] <- "EXON"

        }

        v[is.na(v)] <- if (isTRUE((lyr$aes_params$alpha %||% 1) < 0.5)) "INTRON" else "EXON"
        v

      }

      dat$.__feature <- feat
      layer_dfs[[length(layer_dfs) + 1L]] <- dat

    }

  }

  if (!length(layer_dfs)) return(p)

  all_cols <- Reduce(union, lapply(layer_dfs, names))
  layer_dfs <- lapply(layer_dfs, function(d) {

    miss <- setdiff(all_cols, names(d))

    if (length(miss)) d[miss] <- NA

    d[, all_cols, drop = FALSE]

  })

  df <- do.call(rbind, layer_dfs)

  #Concise Hover

  common_labels <- c(
    gene_id       = "Gene ID",
    label         = "Gene",
    gene_biotype  = "Biotype",
    strand        = "Strand",
    tx_start      = "Tx start",
    tx_end        = "Tx end",
    tx_length     = "Tx length"
  )

  hover_keys_all <- c("gene_id","label","gene_biotype","strand")
  present_hover_all <- intersect(hover_keys_all, names(df))

  format_any <- function(v) {

    if (length(v) == 0 || is.na(v)) return("")

    if (is.numeric(v)) {

      if (is.finite(v) && v == floor(v)) return(formatC(v, format = "d", big.mark = ","))

      return(format(signif(v, 6), trim = TRUE, scientific = FALSE, big.mark = ","))

    }

    as.character(v)

  }

  df$Hover_Info <- vapply(seq_len(nrow(df)), function(i) {

    feat_txt <- ifelse(toupper(df$.__feature[i]) == "INTRON", "INTRON", "EXON")
    header  <- sprintf("<b>%s</b>", feat_txt)
    keys_i <- present_hover_all

    if (identical(feat_txt, "INTRON")) keys_i <- setdiff(keys_i, "transcript_id")

    lines <- character(0)

    if (length(keys_i)) {

      vals <- vapply(keys_i, function(cl) format_any(df[[cl]][i]), character(1))
      labs <- unname(common_labels[keys_i])
      lines <- sprintf("<b>%s:</b> %s", labs, vals)

    }

    paste0(header, if (length(lines)) paste0("<br/>", paste(lines, collapse = "<br/>")) else "")

  }, character(1))

  ok <- is.finite(df$start) & is.finite(df$end) & is.finite(df$y)
  cap <- df[ok, , drop = FALSE]

  if (!nrow(cap)) return(p)

  cap$.mid <- (cap$start + cap$end) / 2
  cap$.id  <- paste0("gene-", seq_len(nrow(cap)))

  p2 <- p + ggiraph::geom_point_interactive(
    data = cap,
    ggplot2::aes(x = .mid, y = y, tooltip = Hover_Info, data_id = .id),
    size = capture_size, alpha = 0.002, stroke = 0,
    inherit.aes = FALSE, show.legend = FALSE
  )
  attr(p2, ".__injected_hover") <- TRUE
  attr(p2, ".__capture_df") <- cap
  p2

}

# HLA ones

hla_extract_slides <- function(out) {
  slides <- list()

  add_slide <- function(plot_obj, label) {
    if (is.null(plot_obj)) return()
    if (inherits(plot_obj, "ggplot") || tube_is_plot_like(plot_obj) || inherits(plot_obj, "patchwork")) {
      slides[[length(slides) + 1L]] <<- list(plot = plot_obj, label = label)
    }
  }

  # Case 1: single result object: class "SegregateHLAPlots"
  if (is.list(out) && any(c("HLA_plot","AA_plot","SNPS_plot","RSID_plots") %in% names(out))) {
    add_slide(out$HLA_plot,  "Classical HLA Alleles")
    add_slide(out$AA_plot,   "AA Sets")
    add_slide(out$SNPS_plot, "SNP Sets")

    # RSID_plots can be named list of patchworks
    if (is.list(out$RSID_plots) && length(out$RSID_plots)) {
      for (nm in names(out$RSID_plots)) {
        add_slide(out$RSID_plots[[nm]], paste0("RSID Regional: ", nm))
      }
    }

    return(slides)
  }

  # Case 2: wrapper returns a list of SegregateHLAPlots (multi datasets)
  # We'll treat each element as "its own output"; caller should handle dataset splitting.
  slides
}


hla_slide_height_in <- function(p, fallback = 10) {
  h <- attr(p, "recommended_height")
  if (is.numeric(h) && is.finite(h) && h > 0) return(as.numeric(h))

  h <- attr(p, "dynamic_height")
  if (is.numeric(h) && is.finite(h) && h > 0) return(as.numeric(h))

  fallback
}

hla_save_slide_png <- function(p, file, width_in, dpi) {
  h <- hla_slide_height_in(p, fallback = 10)
  tube_save_plot_png(p, file, width_in = width_in, height_in = h, dpi = dpi)
}


#Build an invisible ggiraph overlay from a ggplot

hover_overlay_from_plot <- function(p, prefix = "") {

  df <- as.data.frame(p$data)

  if (!nrow(df)) return(p)

  pick <- function(cands) { for (nm in cands) if (nm %in% names(df)) return(nm); NULL }

  xcol <- pick(c("new_pos","GENPOS","POS","BP","position"))

  if (is.null(xcol)) stop("Could not find genomic position column in plot$data.")

  y_vals <- NULL
  y_alt  <- pick(c("LOG10P","MLP","log10P","Log10P"))

  if ("P" %in% names(df) && is.numeric(df$P)) y_vals <- -log10(df$P)

  if (is.null(y_vals) && !is.null(y_alt))      y_vals <- df[[y_alt]]

  if (is.null(y_vals)) stop("No P or log10P-like column in plot$data to build hover y.")

  #significance filter (P < 5e-8)

  p_cut      <- 5e-8
  logp_thr   <- -log10(p_cut)  # ~7.30103
  keep_mask  <- rep(TRUE, nrow(df))

  if ("P" %in% names(df) && is.numeric(df$P)) {

    keep_mask <- is.finite(df$P) & (df$P < p_cut)

  } else if (!is.null(y_alt) && is.numeric(df[[y_alt]])) {


    keep_mask <- is.finite(df[[y_alt]]) & (df[[y_alt]] >= logp_thr)

  }

  df <- df[keep_mask, , drop = FALSE]
  y_vals <- y_vals[keep_mask]

  if (!nrow(df)) return(p)

  idcol  <- pick(c("ID","SNP","RSID","RS","RS_NUMBER"))
  chrcol <- pick(c("CHROM","CHR","chrom","Chromosome"))
  poscol <- pick(c("GENPOS","POS","BP","position"))
  refcol <- pick(c("REF","REF_ALLELE","ALLELE0","A1","EA","EFFECT_ALLELE"))
  altcol <- pick(c("ALT","ALT_ALLELE","ALLELE1","A2","NEA","NON_EFFECT_ALLELE"))

  #tooltip text

  if ("Hover_Info" %in% names(df)) {

    tooltip_vec <- as.character(df$Hover_Info)
    tooltip_vec[is.na(tooltip_vec)] <- ""

  } else {

    ref_lab <- "REF"; alt_lab <- "ALT"

    mkline <- function(lbl, col) if (!is.null(col)) paste0(lbl, df[[col]]) else NULL

    tt_parts_list <- list(
      mkline("SNP: ", idcol),
      mkline("CHR: ", chrcol),
      mkline("POS: ", poscol),

      if ("P" %in% names(df)) paste0("P: ", signif(df$P, 4)) else NULL,

      if (!is.null(y_alt))    paste0("-log10P: ", signif(df[[y_alt]], 4)) else NULL,

      if (!is.null(refcol))   mkline(paste0(ref_lab, ": "), refcol) else NULL,

      if (!is.null(altcol))   mkline(paste0(alt_lab, ": "), altcol) else NULL

    )
    tooltip_vec <- vapply(seq_len(nrow(df)), function(i) {

      parts_i <- lapply(tt_parts_list, function(p) if (length(p)) p[i] else NULL)
      paste0(Filter(Negate(is.null), parts_i), collapse = "\n")

    }, "", USE.NAMES = FALSE)

  }

  src_idx <- seq_len(nrow(df))
  hover_df <- data.frame(
    .x  = df[[xcol]],
    .y  = y_vals,
    .id_raw  = if (!is.null(idcol)) as.character(df[[idcol]]) else as.character(src_idx),
    .tt      = tooltip_vec,
    .src_idx = src_idx,
    stringsAsFactors = FALSE
  )

  ok <- is.finite(hover_df$.x) & is.finite(hover_df$.y)
  hover_df <- hover_df[ok, , drop = FALSE]

  if (!nrow(hover_df)) return(p)

  hover_df$.id <- paste0(prefix, hover_df$.id_raw)

  rows_map <- setNames(
    lapply(hover_df$.src_idx, function(i) {

      vals <- vapply(df[i, , drop = FALSE], function(v) {

        v <- as.character(v); v[is.na(v)] <- ""; paste(v, collapse = ", ")

      }, "")
      as.list(vals)

    }),
    hover_df$.id
  )

  p <- p + ggiraph::geom_point_interactive(
    data = hover_df,
    ggplot2::aes(x = .x, y = .y, tooltip = .tt, data_id = .id),
    inherit.aes = FALSE,
    size = 8, shape = 16, stroke = 0, alpha = 0
  )
  attr(p, "rows_map") <- rows_map
  attr(p, "cols")     <- names(df)
  p

}


sanitize <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
squash_progress_lines <- function(lines) {

  if (is.null(lines) || !length(lines)) return(lines)

  n <- gsub("\r", "", lines, fixed = TRUE)

  n2 <- vapply(n, function(s) {

    m <- gregexpr("\\[[= >-]*\\][[:space:]]*\\d+%|\\b\\d+%\\b", s)
    idx <- m[[1]]

    if (length(idx) > 0 && idx[1] != -1) {

      lens <- attr(m[[1]], "match.length")
      k <- length(idx)
      substr(s, idx[k], idx[k] + lens[k] - 1)

    } else {

      s

    }
  }, character(1), USE.NAMES = FALSE)

  is_bar <- grepl("^[[:space:]]*(\\[[= >-]*\\][[:space:]]*\\d+%|\\d+%)[[:space:]]*$", n2)

  if (!any(is_bar)) return(n2)

  c(n2[!is_bar], tail(n2[is_bar], 1))

}


make_input <- function(name, default) {

  def_str <- paste(deparse(default, width.cutoff = 500L), collapse = "")

  if (isTRUE(default) || isFALSE(default)) checkboxInput(name, name, value = default)

  else if (is.numeric(default) && length(default) == 1) numericInput(name, name, value = default)

  else textInput(name, name, value = def_str)

}

safe_eval_default_reg <- function(expr, fun_env) tryCatch(eval(expr, envir = fun_env), error = function(e) NULL)

make_input_like_meta <- function(prefix, name, default) {

  inputId <- paste0(prefix, name)

  if (is.logical(default) && length(default) == 1)

    checkboxInput(inputId, label = name, value = default)

  else if (is.numeric(default) && length(default) == 1)

    numericInput(inputId, label = name, value = default)

  else if (is.character(default) && length(default) == 1)

    textInput(inputId, label = name, value = default)

  else
    textInput(inputId, label = paste0(name, ""),
              value = deparse(default, width.cutoff = 60))

}

parse_input_val_meta <- function(x) {

  if (is.null(x)) return(NULL)

  if (isTRUE(x) || identical(x, "TRUE")) return(TRUE)

  if (identical(x, FALSE) || identical(x, "FALSE")) return(FALSE)

  tryCatch(eval(parse(text = x)), error = function(e) x)

}

format_arg_val_meta <- function(v) {

  if (is.null(v)) return(NULL)

  if (is.logical(v) && length(v) == 1) return(ifelse(v, "TRUE", "FALSE"))

  if (is.numeric(v)) return(paste0("c(", paste(v, collapse = ","), ")"))

  if (is.character(v)) {

    if (length(v) == 1) return(paste0("\"", v, "\""))

    return(paste0("c(", paste(paste0("\"", v, "\""), collapse = ","), ")"))

  }

  deparse(v, width.cutoff = 60)

}

zip_it <- function(files, zipfile) {

  if (!length(files)) stop("No files to zip")

  if (requireNamespace("zip", quietly = TRUE)) {

    zip::zipr(zipfile, files, include_directories = FALSE)

  } else {

    old <- setwd(dirname(files[1])); on.exit(setwd(old), add = TRUE)

    utils::zip(zipfile, basename(files))

  }

}

reload_miamir <- function() {


  # if (requireNamespace("devtools", quietly = TRUE)) {
  #
  #   try(devtools::load_all(".", quiet = TRUE), silent = TRUE)
  #
  # } else if (requireNamespace("pkgload", quietly = TRUE)) {
  #
  #   try(pkgload::load_all(".", quiet = TRUE, reset = TRUE), silent = TRUE)
  #
  # }

}

#Tube_Alloys helpers

get_tube_fun <- function() {

  assign("use_wrapper", TRUE, envir = .GlobalEnv)

  if (exists("Tube_Alloys", mode = "function", envir = .GlobalEnv)) {

    get("Tube_Alloys", envir = .GlobalEnv)

  } else if (requireNamespace("MiamiR", quietly = TRUE)) {

    MiamiR::Tube_Alloys

  } else {

    stop("Tube_Alloys() not found. Define it in .GlobalEnv or install/load MiamiR.")

  }

}


tube_is_plot_like <- function(x) {

  inherits(x, "ggplot") || inherits(x, c("grob", "gTree", "gtable", "gTable"))

}

tube_get_dyn_height <- function(p, fallback = 7) {

  h <- attr(p, "dynamic_height")

  if (is.numeric(h) && is.finite(h)) as.numeric(h) else fallback

}

tube_save_plot_png <- function(p, file, width_in, height_in, dpi = 100) {

  if (inherits(p, "ggplot")) {

    ggsave(file, plot = p, width = width_in, height = height_in, dpi = dpi,
           units = "in", limitsize = FALSE)

  } else {

    gg <- ggplotify::as.ggplot(p)
    ggsave(file, plot = gg, width = width_in, height = height_in, dpi = dpi,
           units = "in", limitsize = FALSE)

  }

}

tube_save_plot_pdf <- function(p, file, width_in, height_in) {

  if (inherits(p, "ggplot")) {

    ggsave(file, plot = p, width = width_in, height = height_in, device = cairo_pdf)
  } else {

    gg <- ggplotify::as.ggplot(p)
    ggsave(file, plot = gg, width = width_in, height = height_in, device = cairo_pdf)

  }

}

#Forest helpers

get_forest_fun <- function() {

  if (exists("Forest_Plot", mode = "function", envir = .GlobalEnv)) {

    get("Forest_Plot", envir = .GlobalEnv)

  } else if (requireNamespace("MiamiR", quietly = TRUE)) {

    MiamiR::Forest_Plot

  } else stop("Forest_Plot() not found.")


}

safe_eval_default_any <- function(expr, env) tryCatch(eval(expr, envir = env), error = function(e) NULL)

make_arg_input <- function(id_prefix, name, default) {

  inputId <- paste0(id_prefix, name); lab <- name

  if (is.logical(default) && length(default) == 1) checkboxInput(inputId, lab, value = isTRUE(default))

  else if (is.numeric(default) && length(default) == 1) numericInput(inputId, lab, value = default %||% NA_real_)

  else if (is.character(default) && length(default) == 1) textInput(inputId, lab, value = default %||% "")

  else textInput(inputId, paste0(lab, ""),

                 value = if (is.null(default)) "" else paste(deparse(default, width.cutoff = 500L), collapse = ""))

}


parse_arg_val <- function(raw) {

  if (is.null(raw)) return(NULL)

  if (isTRUE(raw) || identical(raw, "TRUE")) return(TRUE)

  if (identical(raw, FALSE) || identical(raw, "FALSE")) return(FALSE)

  if (is.character(raw)) {

    if (!nzchar(raw)) return(NULL)

    out <- tryCatch(eval(parse(text = raw), envir = .GlobalEnv), error = function(e) raw)

    return(out)

  }

  raw

}

format_arg_preview <- function(v) {

  if (is.null(v)) return(NULL)

  if (is.logical(v) && length(v) == 1) return(ifelse(v,"TRUE","FALSE"))

  if (is.numeric(v)) return(paste0("c(", paste(v, collapse=","), ")"))

  if (is.character(v)) {

    if (length(v) == 1) return(paste0("\"", v, "\""))

    return(paste0("c(", paste(paste0("\"", v, "\""), collapse=","), ")"))

  }

  paste(deparse(v, width.cutoff = 60), collapse = "")

}

#UI
ui <- fluidPage(

  titlePanel("MiamiR Shiny App"),
  tags$head(

    tags$style(HTML("

                /* Light background panel behind tiles */

                .home-section {
                  max-width: 1400px;
                  margin: 24px auto;
                  padding: 32px 28px 40px;
                  background: #f7fafc;
                  border: 1px solid #e5e7eb;
                  border-radius: 20px;
                }

                /* Grid and larger tiles */

                .tile-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(360px, 1fr)); gap: 28px; }
                .tile-grid .shiny-input-container, .tile-grid .form-group { margin: 0; }

                .tile-grid .shiny-input-container > .btn.tile,
                .tile-grid .form-group > .btn.tile {
                  display: block; width: 100%; text-align: left;
                  padding: 34px 32px; border-radius: 22px;
                  border: 1px solid #e5e7eb; background: #fff;
                  box-shadow: 0 4px 14px rgba(0,0,0,.08);
                  font-size: 24px; font-weight: 800; white-space: normal;
                  min-height: 200px; line-height: 1.22; color: #1f2937;
                  transition: transform .08s ease, box-shadow .12s ease, filter .12s ease;
                }
                .tile-caption { display:block; margin-top:12px; font-size:16px; font-weight:500; color:rgba(31,41,55,.75); }
                .tile-grid .shiny-input-container > .btn.tile:hover,
                .tile-grid .form-group > .btn.tile:hover { filter: brightness(0.985); box-shadow: 0 12px 26px rgba(0,0,0,.12); transform: translateY(-1px); }

                /* Per-tile pastel colors (stable regardless of order) */

                .btn.tile.tile--miami    { background:#FDE2E4; border-color:#F8C9CE; }
                .btn.tile.tile--single   { background:#E8F7F2; border-color:#CFECE2; }
                .btn.tile.tile--regional { background:#E7F0FF; border-color:#CCDBFF; }
                .btn.tile.tile--annotate { background:#F5E6FF; border-color:#E2CCFF; }
                .btn.tile.tile--meta     { background:#FFF4D6; border-color:#FBE6A9; }
                .btn.tile.tile--forest   { background:#EAF3EA; border-color:#D4E7D4; }
                .btn.tile.tile--model    { background:#FDE7F3; border-color:#F7CDE6; }
                .btn.tile.tile--tube     { background:#EAF0FF; border-color:#D3E0FF; }
                .btn.tile.tile--hla { background:#EFFFF6; border-color:#CDEFE0; }


                /* Wider screens: even larger */

                @media (min-width: 1500px) {
                  .tile-grid { grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 32px; }
                  .tile-grid .shiny-input-container > .btn.tile,
                  .tile-grid .form-group > .btn.tile { min-height: 220px; font-size: 26px; padding: 38px 34px; }
                }

              ")),

    tags$style(HTML("
            /* Slightly lower the Miami save-type dropdown so it lines up with buttons */

            #miami_ext + .selectize-control { margin-top: 2px; }

            .btnbar .shiny-input-container { margin: 0; } /* remove extra vertical margin in toolbars */

          ")),

    tags$style(HTML("

            /* a smidge more space after the Run button */

            .btn--gap-right { margin-right: 14px !important; }

          ")),

    tags$style(HTML("

            .btnbar--spaced { margin-bottom: 22px; }  /* extra breathing room below the buttons */

          ")),


    tags$style(HTML("

          /* METASOFT preview: don't overflow the panel */

          #meta_preview_wrap { max-width:100%; overflow-x:auto; -webkit-overflow-scrolling:touch; }
          #meta_preview_wrap table { width:auto !important; table-layout:auto; white-space:nowrap; }

        ")),

    tags$style(HTML("

              /* FORCE the hover/press animation on the tiles */

              .tile-grid .btn.tile {
                cursor: pointer;
                will-change: transform, box-shadow, filter;
                transition:
                  transform .18s cubic-bezier(.2,.8,.2,1),
                  box-shadow .18s cubic-bezier(.2,.8,.2,1),
                  filter .18s cubic-bezier(.2,.8,.2,1) !important;
                transform: translateZ(0); /* GPU nudge */
              }

              .tile-grid .btn.tile:hover,
              .tile-grid .btn.tile:focus-visible {
                transform: translateY(-10px) scale(1.03) !important;
                box-shadow: 0 18px 36px rgba(0,0,0,.22) !important;
                filter: brightness(0.985) !important;
              }

              .tile-grid .btn.tile:active {
                transform: translateY(-2px) scale(0.995) !important;
                box-shadow: 0 10px 20px rgba(0,0,0,.18) !important;
              }

              /* Motion-safe fallback */

              @media (prefers-reduced-motion: reduce) {
                .tile-grid .btn.tile { transition: none !important; }
                .tile-grid .btn.tile:hover,
                .tile-grid .btn.tile:focus-visible {
                  transform: none !important;
                  box-shadow: 0 8px 16px rgba(0,0,0,.12) !important;

                }

              }

            ")),

    #THEME - per-tab color variables

    tags$style(HTML("
                :root {
                  --accent-bg: #ffffff;
                  --accent-border: #e5e7eb;
                  --accent-ink: #1f2937;
                }
                :root[data-theme='miami']    { --accent-bg:#FDE2E4; --accent-border:#F8C9CE; }
                :root[data-theme='single']   { --accent-bg:#E8F7F2; --accent-border:#CFECE2; }
                :root[data-theme='regional'] { --accent-bg:#E7F0FF; --accent-border:#CCDBFF; }
                :root[data-theme='annotate'] { --accent-bg:#F5E6FF; --accent-border:#E2CCFF; }
                :root[data-theme='meta']     { --accent-bg:#FFF4D6; --accent-border:#FBE6A9; }
                :root[data-theme='forest']   { --accent-bg:#EAF3EA; --accent-border:#D4E7D4; }
                :root[data-theme='model']    { --accent-bg:#FDE7F3; --accent-border:#F7CDE6; }
                :root[data-theme='tube']     { --accent-bg:#EAF0FF; --accent-border:#D3E0FF; }
                :root[data-theme='hla'] { --accent-bg:#EFFFF6; --accent-border:#CDEFE0; }


                /* Themed wrapper that tints each tool tab to match its tile */

                .themed-tab {
                  background: var(--accent-bg);
                  border: 1px solid var(--accent-border);
                  border-radius: 20px;
                  padding: 24px;
                  margin: 24px auto;
                  max-width: 1400px;
                  box-shadow: 0 4px 14px rgba(0,0,0,.04);
                }

                /* make little bits feel cohesive */

                .themed-tab .link-back { color:#111827; opacity:.75; text-decoration:none; }
                .themed-tab .link-back:hover { opacity:1; text-decoration:underline; }
                .themed-tab .form-control, .themed-tab .selectize-input {
                  background: #fff; border-color: var(--accent-border);
                }
                .themed-tab hr { border-top-color: var(--accent-border); }

              ")),

    # JS hook: switch theme when the active tab changes
    tags$script(HTML("

                Shiny.addCustomMessageHandler('set-theme', function(theme) {
                  const root = document.documentElement;
                  if (!theme || theme === 'default') root.removeAttribute('data-theme');
                  else root.setAttribute('data-theme', theme);

                });

              ")),


    tags$style(HTML('

            /* Split dropdown (button + caret) */

            .btnbar { display:flex; gap:10px; align-items:center; flex-wrap:wrap; }
            .split { position:relative; display:inline-flex; border-radius:10px; overflow:hidden; }
            .split > button { border:1px solid #e5e7eb; background:#111827; color:#fff; padding:8px 14px; font-weight:600; }
            .split > button:first-child { border-right:none; border-top-left-radius:10px; border-bottom-left-radius:10px; }
            .split > button:last-child { width:40px; border-left:1px solid rgba(255,255,255,.25); border-top-right-radius:10px; border-bottom-right-radius:10px; }
            .split.open > .menu { display:block; }
            .split .menu { display:none; position:absolute; top:100%; right:0; min-width:220px; background:#fff; border:1px solid #e5e7eb; border-radius:12px; box-shadow:0 10px 30px rgba(0,0,0,.12); z-index:1000; }
            .split .menu a { display:block; padding:10px 14px; color:#111827; text-decoration:none; }
            .split .menu a:hover { background:#f3f4f6; }
            .details-wrap details { border:1px solid var(--accent-border); background:#fff; border-radius:12px; padding:10px 14px; }
            .details-wrap summary { font-weight:700; cursor:pointer; }
            .details-wrap details + details { margin-top:10px; }

            ')),

    tags$style(HTML("

          /* --- Toolbar consistency: inline, fixed-size controls --- */

          .btnbar { display:flex; align-items:center; gap:10px; flex-wrap:wrap; }
          .btnbar--nowrap { flex-wrap:nowrap; }                 /* add to toolbars that must not wrap */
          .btnbar--spaced { margin-bottom:22px; }               /* optional breathing room below a bar */

          /* make each control size-to-content and sit inline */

          .btnbar .shiny-input-container,
          .btnbar .form-group { width:auto !important; margin:0; flex:0 0 auto; }

          /* keep selectize compact & vertically aligned in any toolbar */

          .btnbar .selectize-control {
            margin-top:2px;            /* align with buttons */
            width:110px !important;    /* match selectInput width */
            min-width:110px;           /* prevent flex grow */
          }

          /* robustness: whether selectize wrapper is after <select> or <label> */

          #miami_ext + .selectize-control,
          label[for='miami_ext'] + .selectize-control,
          #single_ext + .selectize-control,
          label[for='single_ext'] + .selectize-control,
          #regional_ext + .selectize-control,
          label[for='regional_ext'] + .selectize-control,
          #annotate_ext + .selectize-control,
          label[for='annotate_ext'] + .selectize-control {
            margin-top:2px;
            width:110px !important;
            min-width:110px;
          }


          /* a touch more space after the Run button */

          .btn--gap-right { margin-right:14px !important; }

          "))

    ,

    #Add this tiny JS once

    tags$script(HTML('

            (function(){
            document.addEventListener("click", function(e){
            const btn = e.target.closest(".split > button:last-child");
            document.querySelectorAll(".split.open").forEach(el=>{ if(!el.contains(e.target)) el.classList.remove("open"); });
            if(btn){ const root = btn.parentElement; root.classList.toggle("open"); }
            });
            })();
            ')),


    tags$style(HTML("

          /* Regional carousels: scale giant images to fit like other tabs */

          .themed-tab .slick-slide { display:flex; align-items:center; justify-content:center; }
          .themed-tab .slick-slide img {
            max-width: 100%;
            max-height: 660px;  /* ~680px container minus a little padding */
            width: auto;
            height: auto;
            display: block;
          }

          ")),

    tags$style(HTML("

          /* DROP-IN: Fix Model_Munge toolbar crowding/overlap */

          [data-value='model'] .btnbar--nowrap { flex-wrap: wrap !important; }

          /* Keep controls compact and aligned on one row when space allows */

          [data-value='model'] .btnbar .shiny-input-container,
          [data-value='model'] .btnbar .form-group {
            margin: 0 !important;
            width: auto !important;
            flex: 0 0 auto !important;
          }

          /* Give the format dropdown a little extra width */

          [data-value='model'] #mm_ext + .selectize-control {
            width: 140px !important;
            min-width: 140px !important;
          }

        ")),

    tags$script(HTML("

          // Map split-menu anchors -> hidden downloadButton IDs

          (function(){
            const map = {
              'tube_dl_single':   { png: 'tube_dl_single_zip_png',   pdf: 'tube_dl_single_zip_pdf' },
              'tube_dl_regional': { png: 'tube_dl_regional_zip_png', pdf: 'tube_dl_regional_zip_pdf' },
              'tube_dl_forest':   { png: 'tube_dl_forest_png',       pdf: 'tube_dl_forest_pdf' },
              'tube_dl_gwas':     { png: 'tube_dl_gwas_png',         pdf: 'tube_dl_gwas_pdf' }
            };

            document.addEventListener('click', function(e){
              const link = e.target.closest('.split .menu a');
              if(!link) return;

              // IDs are like 'tube_dl_single-png' or 'tube_dl_gwas-pdf'
              const m = link.id && link.id.match(/^(.*)-(png|pdf)$/);
              if(!m) return;

              e.preventDefault();
              const base = m[1], type = m[2];
              const targetId = map[base] && map[base][type];
              const target = targetId && document.getElementById(targetId);
              if(target){
                // close any open menu
                const root = link.closest('.split');
                if(root) root.classList.remove('open');
                target.click();

              }

            });

          })();

        ")),


    # Slightly lower the save-type dropdown
    tags$style(HTML("

          #miami_ext + .selectize-control,
          #single_ext + .selectize-control,
          #regional_ext + .selectize-control,
          #annotate_ext + .selectize-control,
          #tube_ext + .selectize-control {
            margin-top: 2px;
          }

        ")),

    # Keep selectise compact

    tags$style(HTML("

          .btnbar .selectize-control {
            margin-top:2px;
            width:110px !important;
            min-width:110px;
          }
          #miami_ext + .selectize-control,
          label[for='miami_ext'] + .selectize-control,
          #single_ext + .selectize-control,
          label[for='single_ext'] + .selectize-control,
          #regional_ext + .selectize-control,
          label[for='regional_ext'] + .selectize-control,
          #annotate_ext + .selectize-control,
          label[for='annotate_ext'] + .selectize-control,
          #tube_ext + .selectize-control,                 /* <-- added */
          label[for='tube_ext'] + .selectize-control {    /* <-- added */
            margin-top:2px;
            width:110px !important;
            min-width:110px;
          }

        ")),

    tags$script(HTML("

          Shiny.addCustomMessageHandler('open-window', function(url){
            window.open(url, '_blank');
          });

        ")),

    tags$script(HTML("

  Shiny.addCustomMessageHandler('click-el', function(x){
    var el = document.getElementById(x.id);
    if (el) el.click();

  });

")),


    tags$style(HTML("

      .modal-content {
        border-radius: 16px; border: none; box-shadow: 0 6px 20px rgba(0,0,0,0.2);
      }
      .modal-title { font-weight: 700; font-size: 20px; color: #2c3e50; }
      .modal-footer { justify-content: space-between; }
      .btn-google {
        background-color:#4285F4; color:#fff; border-radius:8px; font-weight:700;
      }
      .btn-google:hover { background-color:#3367D6; color:#fff; }
      .btn-dbsnp {
        background-color:#28a745; color:#fff; border-radius:8px; font-weight:700;
      }
      .btn-dbsnp:hover { background-color:#218838; color:#fff; }

    "))

  )

  ,

  # Hidden tabset; start on Home (tiles)

  tabsetPanel(

    id = "main_tabs",
    type = "hidden",
    selected = "home",

    #Home (tiles)

    tabPanel(
      "Home", value = "home",
      div(class = "home-section",
          h2(""),

          # Main tiles (UNCHANGED)
          div(class = "tile-grid",
              actionButton("tile_miami",    HTML("Miami_Plot<span class='tile-caption'>Mirrored Manhattan</span>"),      class = "btn tile tile--miami"),
              actionButton("tile_single",   HTML("Single_Plot<span class='tile-caption'>Manhattan(s)</span>"),          class = "btn tile tile--single"),
              actionButton("tile_regional", HTML("Regional_Plot<span class='tile-caption'>Regional Plot(s)</span>"),    class = "btn tile tile--regional"),
              actionButton("tile_annotate", HTML("Annotate_Data<span class='tile-caption'>Query RSIDs</span>"),         class = "btn tile tile--annotate"),
              actionButton("tile_meta",     HTML("METASOFT_File_Gen<span class='tile-caption'>Meta-analysis Input Curation</span>"), class = "btn tile tile--meta"),
              actionButton("tile_forest",   HTML("Forest_Plot<span class='tile-caption'>Evaluate Select SNPs</span>"),  class = "btn tile tile--forest"),
              actionButton("tile_model",    HTML("Model_Munge<span class='tile-caption'>Fit & Munge Models</span>"),    class = "btn tile tile--model"),
              actionButton("tile_tube",     HTML("Tube_Alloys<span class='tile-caption'>One-click Overview Pipeline</span>"), class = "btn tile tile--tube"),
              actionButton("tile_hla",
                           HTML("Segregate_HLA_Plot<span class='tile-caption'>HLA segregation plots</span>"),
                           class = "btn tile tile--hla"
              )

          ),

          # Slight drop / separate section
          div(style = "margin-top: 28px;"),

          tags$h4(style="margin: 0 0 14px 4px; font-weight: 800; color: rgba(31,41,55,.85);",
                  "Resources"),

          div(class = "tile-grid",
              actionButton(
                "help_btn",
                HTML("Help<span class='tile-caption'>Browse MiamiR Documentation</span>"),
                class = "btn tile tile--help"
              ),
              downloadButton(
                "download_demo_data",
                label = HTML("Download demo data<span class='tile-caption'>MiamiR_Demo_Data.tar</span>"),
                class = "btn tile tile--download"
              )
          )
      )
    )
    ,

    #Seg HLA

    tabPanel("Segregate_HLA_Plot", value = "hla",
             div(class = "themed-tab",
                 sidebarLayout(
                   sidebarPanel(
                     list(
                       actionLink("back_home_hla", "\u2190 All tools", class = "link-back"),
                       h4("Upload data (multiple allowed)"),
                       fileInput("hla_files", "Data file(s)", multiple = TRUE,
                                 accept = c(".csv", ".tsv", ".txt", ".gz", ".bgz")),

                       div(class = "btnbar btnbar--spaced",
                           actionButton("run_hla", "Run Segregate_HLA_Plot", class = "btn-primary btn--gap-right"),
                           downloadButton("hla_download_simple", "Download"),
                           selectInput("hla_ext", NULL,
                                       choices = c("PNG"="png","PDF"="pdf","JPG"="jpg","SVG"="svg"),
                                       selected = "png", width = "110px")
                       ),

                       div(class = "details-wrap",
                           tags$details(
                             tags$summary("Segregate_HLA_Plot settings"),
                             uiOutput("hla_arg_inputs")
                           ),
                           tags$details(
                             tags$summary("Render options"),
                             numericInput("hla_save_width",  "Save Width (inches)",  value = 30, min = 4,  max = 60, step = 1),
                             numericInput("hla_save_height", "Save Height (inches)", value = 10, min = 3,  max = 40, step = 0.5),
                             numericInput("hla_save_dpi",    "Save DPI",             value = 100, min = 72, max = 1200, step = 10)
                           ),
                           tags$details(
                             tags$summary("Show call preview"),
                             div(class="call-preview", verbatimTextOutput("hla_call_preview", placeholder = TRUE))
                           )
                       )
                     )
                   ),
                   mainPanel(
                     uiOutput("hla_dynamic_tabs")
                   )
                 )
             )
    ),


    #Miami_Plot

    tabPanel("Miami_Plot", value = "miami",
             div(class = "themed-tab",
                 sidebarLayout(
                   sidebarPanel(
                     list(
                       actionLink("back_home_miami", "\u2190 All tools", class = "link-back"),
                       h4("Upload data"),
                       fileInput("top_file", "Top data file", accept = c(".csv", ".tsv", ".txt", ".gz", ".bgz")),
                       fileInput("bottom_file", "Bottom data file", accept = c(".csv", ".tsv", ".txt", ".gz", ".bgz")),
                       div(class = "btnbar btnbar--spaced",
                           actionButton("run", "Run Miami_Plot", class = "btn-primary btn--gap-right"),
                           actionButton("miami_interactive_view", "Interactive View", class = "btn btn-info btn--gap-right"),  # <-- add this
                           downloadButton("miami_download_png_simple", "Download"),
                           selectInput("miami_ext", NULL,
                                       choices  = c("PNG"="png","PDF"="pdf","JPG"="jpg","SVG"="svg"),
                                       selected = "png", width = "110px")
                       )
                       ,

                       div(class = "details-wrap",

                           tags$details(
                             tags$summary("Miami_Plot settings"),
                             h4("Miami_Plot specific arguments"),
                             textInput("top_title", "Top_Title",    placeholder = "Leave blank to infer"),
                             textInput("bottom_title", "Bottom_Title", placeholder = "Leave blank to infer"),

                             tags$hr(),
                             h4("Argument mode"),
                             radioButtons(
                               "miami_mode", NULL,
                               choices = c("Simple (Apply to Top & Bottom)" = "simple",
                                           "Advanced (Separate to Top & Bottom)"      = "advanced"),
                               selected = "advanced", inline = TRUE
                             ),

                             conditionalPanel("input.miami_mode == 'simple'",
                                              uiOutput("miami_args_simple")
                             ),

                             conditionalPanel("input.miami_mode == 'advanced'",
                                              div(id = "miami-adv-switch", class = "btnbar",
                                                  tags$label("Edit panel"),
                                                  radioButtons(
                                                    "miami_adv_panel", label = NULL,
                                                    choices  = c("Top" = "top", "Bottom" = "bottom"),
                                                    selected = "top", inline = TRUE
                                                  )
                                              ),
                                              conditionalPanel("input.miami_adv_panel == 'top'",
                                                               uiOutput("miami_args_ui_top")),
                                              conditionalPanel("input.miami_adv_panel == 'bottom'",
                                                               uiOutput("miami_args_ui_bottom"))
                             )

                           ),

                           tags$details(
                             tags$summary("Render options"),
                             numericInput("save_width",  "Save Width (inches)",  value = 30, min = 4,  max = 60, step = 1),
                             numericInput("save_height", "Save Height (inches)", value = 15, min = 3,  max = 40, step = 0.5),
                             numericInput("save_dpi",    "Save DPI",             value = 100, min = 72, max = 1200, step = 10)
                           ),

                           tags$details(
                             tags$summary("Show call preview"),
                             div(class = "call-preview", verbatimTextOutput("miami_call_preview", placeholder = TRUE))
                           )
                       )
                       ,

                       hidden = list(
                         downloadButton("download_png", "PNG", style = "display:none"),
                         downloadButton("download_pdf", "PDF", style = "display:none"),
                         downloadButton("download_jpg", "JPG", style = "display:none")
                       )
                     )
                   ),
                   mainPanel(
                     imageOutput("miami_plot")
                   )

                 )
             )
    ),

    #Single_Plot

    tabPanel("Single_Plot", value = "single",
             div(class = "themed-tab",
                 sidebarLayout(
                   sidebarPanel(
                     list(
                       actionLink("back_home_single", "\u2190 All tools", class = "link-back"),
                       h4("Upload data (multiple allowed)"),
                       fileInput("single_file", "Data file(s)", multiple = TRUE,
                                 accept = c(".csv", ".tsv", ".txt", ".gz", ".bgz")),
                       div(class = "btnbar btnbar--spaced",
                           actionButton("run_single", "Run Single_Plot", class = "btn-primary btn--gap-right"),
                           actionButton("single_interactive_view", "Interactive View", class = "btn btn-info btn--gap-right"),
                           downloadButton("single_download_simple", "Download"),
                           selectInput(
                             "single_ext", NULL,
                             choices  = c("PNG" = "png", "PDF" = "pdf", "JPG" = "jpg", "SVG" = "svg"),
                             selected = "png", width = "110px"
                           )
                       ),
                       div(class = "details-wrap",
                           tags$details(
                             tags$summary("Single_Plot settings"),
                             uiOutput("single_args_ui2")
                           ),
                           tags$details(
                             tags$summary("Render options"),
                             numericInput("save_width_single",  "Save Width (inches)",  value = 30, min = 4,  max = 60, step = 1),
                             numericInput("save_height_single", "Save Height (inches)", value = 10, min = 3,  max = 40, step = 0.5),
                             numericInput("save_dpi_single",    "Save DPI",             value = 100, min = 72, max = 1200, step = 10)
                           ),
                           tags$details(tags$summary("Show call preview"),
                                        div(class = "call-preview", verbatimTextOutput("single_call_preview", placeholder = TRUE))
                           )
                       )
                       ,
                       hidden = list(
                         downloadButton("download_png_single", "PNG", style = "display:none"),
                         downloadButton("download_pdf_single", "PDF", style = "display:none"),
                         downloadButton("download_jpg_single", "JPG", style = "display:none")
                       )

                     )
                   ),
                   mainPanel(uiOutput("single_dynamic_tabs"))

                 )
             )
    ),

    #Regional_Plot

    tabPanel("Regional_Plot", value = "regional",
             div(class = "themed-tab",
                 sidebarLayout(
                   sidebarPanel(
                     list(
                       actionLink("back_home_regional", "\u2190 All tools", class = "link-back"),
                       h4("Upload data (multiple allowed)"),
                       fileInput("regional_datafile", "Data file(s)", multiple = TRUE,
                                 accept = c(".csv", ".tsv", ".txt", ".gz", ".bgz")),
                       div(
                         class = "btnbar btnbar--spaced",
                         actionButton("run_regional", "Run Regional_Plot", class = "btn-primary btn--gap-right"),
                         actionButton("regional_interactive_view", "Interactive View", class = "btn btn-info btn--gap-right"),
                         tags$div(style = "flex-basis:100%; height:0;"),  # ← hard line break in the flex row
                         downloadButton("regional_download_simple", "Download"),
                         selectInput(
                           "regional_ext", NULL,
                           choices  = c("PNG" = "png", "PDF" = "pdf", "JPG" = "jpg", "SVG" = "svg"),
                           selected = "png", width = "110px"
                         )
                       )
                       ,

                       div(class = "details-wrap",
                           tags$details(
                             tags$summary("Regional_Plot settings"),
                             numericInput("regional_chrom_input", "Chromosomes", value = 4, min = 1, step = 1),
                             uiOutput("regional_arg_inputs")
                           ),
                           tags$details(
                             tags$summary("Render options"),
                             numericInput("regional_save_width", "Save Width (inches)", value = 30, min = 4, max = 60, step = 1),
                             numericInput("regional_save_dpi",   "Save DPI",            value = 100, min = 72, max = 1200, step = 10)
                           ),
                           tags$details(tags$summary("Show call preview"),
                                        div(class = "call-preview", verbatimTextOutput("regional_call_preview", placeholder = TRUE))
                           )
                       )

                     )
                   ),
                   mainPanel(uiOutput("regional_dynamic_tabs"))
                 )
             )
    ),

    #Annotate_Data

    tabPanel("Annotate_Data", value = "annotate",
             div(class = "themed-tab",
                 sidebarLayout(
                   sidebarPanel(
                     list(
                       actionLink("back_home_annotate", "\u2190 All tools", class = "link-back"),
                       h4("Upload data (multiple allowed)"),
                       fileInput("annotate_files", "Data file(s)", multiple = TRUE,
                                 accept = c(".csv", ".tsv", ".txt", ".gz", ".bgz")),
                       div(class = "btnbar btnbar--spaced btnbar--nowrap",
                           actionButton("run_annotate", "Run Annotate_Data", class = "btn-primary btn--gap-right"),
                           downloadButton("annotate_download_simple", "Download"),
                           selectInput(
                             "annotate_ext", NULL,
                             choices  = c("CSV" = "csv", "TSV" = "tsv", "TXT" = "txt"),
                             selected = "csv", width = "110px"
                           )
                       ),
                       div(class = "details-wrap",
                           tags$details(tags$summary("Annotate_Data settings"), uiOutput("annotate_arg_inputs")),
                           tags$details(tags$summary("Show call preview"),
                                        div(class = "call-preview", verbatimTextOutput("annotate_call_preview", placeholder = TRUE))
                           )
                       )
                     )
                   ),
                   mainPanel(uiOutput("annotate_dynamic_tabs"))
                 )
             )
    ),

    #METASOFT_File_Gen

    tabPanel("METASOFT_File_Gen", value = "meta",
             div(class = "themed-tab",
                 sidebarLayout(
                   sidebarPanel(
                     width = 5,
                     list(
                       actionLink("back_home_meta", "\u2190 All tools", class = "link-back"),
                       h4("Upload data (min. 2 allowed)"),
                       fileInput("meta_files", "Data file(s)", multiple = TRUE,
                                 accept = c(".csv", ".tsv", ".txt", ".gz", ".bgz")),
                       div(class = "btnbar btnbar--spaced btnbar--nowrap",
                           actionButton("run_meta", "Run METASOFT_File_Gen", class = "btn-primary btn--gap-right"),
                           downloadButton("meta_download_simple", "Download"),
                           selectInput("meta_ext", NULL,
                                       choices = c("CSV" = "csv", "TSV" = "tsv", "TXT" = "txt"),
                                       selected = "csv", width = "110px")
                       ),
                       div(class = "details-wrap",
                           tags$details(tags$summary("METASOFT_File_Gen settings"), uiOutput("meta_arg_inputs")),
                           tags$details(tags$summary("Show call preview"),
                                        div(class = "call-preview", verbatimTextOutput("meta_call_preview", placeholder = TRUE))
                           )
                       )
                     )
                   ),
                   mainPanel(
                     width = 7,
                     div(id = "meta_preview_wrap",
                         tableOutput("meta_output_preview")
                     )
                   )

                 )
             )
    )
    ,

    #Forest_Plot

    tabPanel("Forest_Plot", value = "forest",
             div(class = "themed-tab",
                 sidebarLayout(
                   sidebarPanel(
                     list(
                       actionLink("back_home_forest", "\u2190 All tools", class = "link-back"),
                       h4("Upload data (multiple allowed)"),
                       fileInput("forest_files", "Data file(s)", multiple = TRUE,
                                 accept = c(".csv", ".tsv", ".txt", ".gz", ".bgz")),
                       div(class = "btnbar btnbar--spaced",
                           actionButton("run_forest", "Run Forest_Plot", class = "btn-primary btn--gap-right"),
                           downloadButton("forest_download_simple", "Download"),
                           selectInput("forest_ext", NULL,
                                       choices  = c("PNG"="png","PDF"="pdf","JPG"="jpg","SVG"="svg"),
                                       selected = "png", width = "110px")
                       ),

                       div(class = "details-wrap",
                           tags$details(
                             tags$summary("Forest_Plot settings"),
                             uiOutput("forest_arg_inputs")
                           ),
                           tags$details(
                             tags$summary("Render options"),
                             numericInput("forest_save_width", "Save Width (inches)", value = 15, min = 3, step = 0.5),
                             numericInput("forest_save_dpi",   "Save DPI",            value = 100, min = 72, max = 1200, step = 10)
                           ),
                           tags$details(
                             tags$summary("Show call preview"),
                             div(class = "call-preview", verbatimTextOutput("forest_call_preview", placeholder = TRUE))
                           )
                       )
                     )
                   ),
                   mainPanel(
                     uiOutput("forest_dynamic_tabs")
                   )
                 )
             )
    )
    ,

    #Model_Munge

    tabPanel("Model_Munge", value = "model",
             div(class = "themed-tab",
                 sidebarLayout(
                   sidebarPanel(
                     list(
                       actionLink("back_home_model", "\u2190 All tools", class = "link-back"),
                       h4("Upload data"),
                       fileInput("mm_file", "Data file", multiple = FALSE,
                                 accept = c(".csv", ".tsv", ".txt", ".gz", ".bgz")),
                       div(class = "btnbar btnbar--spaced btnbar--nowrap",
                           actionButton("mm_fit", "Fit model & save", class = "btn-primary btn--gap-right"),
                           actionButton("mm_run_munge", "Run Model_Munge", class = "btn-success btn--gap-right"),
                           downloadButton("mm_download_simple", "Download munged"),
                           selectInput("mm_ext", NULL,
                                       choices = c("CSV" = "csv", "TSV" = "tsv", "TXT" = "txt"),
                                       selected = "csv", width = "110px")
                       )
                       ,
                       div(class = "details-wrap",
                           tags$details(tags$summary("Model builder"),
                                        textInput("mm_dataset_name", "Dataset object name", value = ""),
                                        radioButtons("mm_model_type", "Model type", choices = c("lm","glm"), inline = TRUE),
                                        uiOutput("mm_var_ui"),
                                        checkboxInput("mm_intercept", "Include intercept", TRUE),
                                        selectInput("mm_glm_family", "GLM family",
                                                    choices = c("gaussian","binomial","poisson","Gamma","inverse.gaussian","quasipoisson","quasibinomial"),
                                                    selected = "gaussian"),
                                        textAreaInput("mm_formula_manual", "Or custom formula (optional, overrides builder)", rows = 2, placeholder = "e.g. y ~ x1 + x2"),
                                        textInput("mm_model_name", "Model object name", value = "")
                           ),
                           tags$details(tags$summary("Model_Munge settings"), uiOutput("mm_arg_inputs")),
                           tags$details(
                             tags$summary("Show call preview"),
                             div(class = "call-preview",
                                 verbatimTextOutput("mm_call_preview", placeholder = TRUE))
                           )
                       )

                     )
                   ),
                   mainPanel(
                     tabsetPanel(
                       tabPanel("Data preview", tableOutput("mm_data_preview")),
                       tabPanel("Model summary", verbatimTextOutput("mm_model_summary")),
                       tabPanel("Munge output (first 10 rows)", tableOutput("mm_munge_preview"))
                     )
                   )
                 )
             )
    ),

    #Tube_Alloys

    tabPanel("Tube_Alloys", value = "tube",

             div(class = "themed-tab",
                 sidebarLayout(
                   sidebarPanel(
                     list(
                       tags$div(style = "display:none;",
                                actionButton("tube_single_interactive_view",   "proxy"),
                                actionButton("tube_regional_interactive_view", "proxy")
                       ),
                       actionLink("back_home_tube", "\u2190 All tools", class = "link-back"),
                       h4("Upload data"),
                       fileInput("tube_data", "Data (required)", multiple = TRUE,
                                 accept = c(".csv",".tsv",".txt",".gz",".bgz")),
                       fileInput("tube_pheno", "Phenos (optional)", accept = c(".csv",".tsv",".txt",".gz",".bgz")),
                       fileInput("tube_covar", "Covars (optional)", accept = c(".csv",".tsv",".txt",".gz",".bgz")),
                       fileInput("tube_pheno_covar", "Phenos_Covars (optional)", accept = c(".csv",".tsv",".txt",".gz",".bgz")),
                       div(
                         class = "btnbar btnbar--spaced",
                         actionButton("tube_run", "Run Tube_Alloys", class = "btn btn-primary btn--gap-right"),
                         actionButton("tube_interactive_view", "Interactive View", class = "btn btn-info btn--gap-right"),

                         downloadButton("tube_download", "Download"),
                         selectInput("tube_ext", NULL, choices = c("PNG" = "png", "PDF" = "pdf"), selected = "png", width = "110px")
                       )

                       ,

                       div(class = "details-wrap",

                           tags$details(
                             tags$summary("Tube_Alloys settings"),
                             uiOutput("tube_arg_inputs")
                           ),

                           tags$details(
                             tags$summary("Render options"),
                             tabsetPanel(
                               id = "tube_render_tabs",
                               tabPanel("Single",
                                        numericInput("tube_single_width",  "Width (in)",  30, min = 4,  max = 60, step = 1),
                                        numericInput("tube_single_height", "Height (in)", 10, min = 3,  max = 40, step = 0.5),
                                        numericInput("tube_single_dpi",    "DPI",        100, min = 72, max = 1200, step = 10)
                               ),
                               tabPanel("Regional",
                                        numericInput("tube_regional_width","Width (in)",  30, min = 4,  max = 60, step = 1),
                                        numericInput("tube_regional_dpi",  "DPI",        100, min = 72, max = 1200, step = 10),
                                        helpText("Height is dynamic per plot.")
                               ),
                               tabPanel("Forest",
                                        numericInput("tube_forest_width",  "Width (in)",  15, min = 3,  max = 40, step = 0.5),
                                        numericInput("tube_forest_dpi",    "DPI",        100, min = 72, max = 1200, step = 10),
                                        helpText("Height is dynamic.")
                               ),
                               tabPanel("GWAS Model",
                                        numericInput("tube_gwas_width",    "Width (in)",  15, min = 3,  max = 40, step = 0.5),
                                        numericInput("tube_gwas_dpi",      "DPI",        100, min = 72, max = 1200, step = 10),
                                        helpText("Height is dynamic.")
                               )
                             )
                           )
                           ,

                           tags$details(
                             tags$summary("Show call preview"),
                             div(class = "call-preview", verbatimTextOutput("tube_call_preview", placeholder = TRUE))
                           )
                       )

                     )
                   ),

                   mainPanel(
                     tabsetPanel(
                       id = "tube_tabs",
                       tabPanel("Single_Plot",   uiOutput("tube_single_ui")),
                       tabPanel("Regional_Plot", uiOutput("tube_regional_ui")),
                       tabPanel("Forest_Plot",   imageOutput("tube_forest_img")),
                       tabPanel("GWAS_Model",    imageOutput("tube_gwas_img"))

                     )
                   )

                 )
             )
    )
  )
)


#Server

server <- function(input, output, session) {

  # server.R (or inside server function)
  help_dir <- system.file("app_help", package = "MiamiR")
  shiny::addResourcePath("miamir_help", help_dir)

  topics <- sub("\\.html$", "", list.files(help_dir, pattern = "\\.html$"))

  observeEvent(input$help_btn, {
    showModal(modalDialog(
      title = "MiamiR Help",
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close"),

      selectInput("help_topic", "Choose a function or dataset for info", choices = sort(topics), selected = "Single_Plot"),
      uiOutput("help_frame")
    ))
  })

  output$help_frame <- renderUI({
    req(input$help_topic)
    tags$iframe(
      src = paste0("miamir_help/", input$help_topic, ".html"),
      style = "width:100%; height:80vh; border:none;"
    )
  })


  # Server
  output$download_demo_data <- downloadHandler(
    filename = function() "MiamiR Demo Data.tar.gz",
    content = function(file) {

      pkg_data_dir <- system.file("data", package = "MiamiR")
      pkg_ext_dir  <- system.file("extdata", package = "MiamiR")

      stopifnot(dir.exists(pkg_data_dir), dir.exists(pkg_ext_dir))

      stage_root <- tempfile("miamir_demo_stage_")
      dir.create(stage_root, recursive = TRUE)

      bundle_root <- file.path(stage_root, "MiamiR Demo Data")
      dir.create(bundle_root, recursive = TRUE)

      copy_tree <- function(src_dir, dest_dir) {
        if (!dir.exists(src_dir)) return(invisible(FALSE))
        dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

        files <- list.files(src_dir, full.names = TRUE, recursive = TRUE, all.files = TRUE, no.. = TRUE)
        files <- files[file.exists(files) & !dir.exists(files)]
        if (!length(files)) return(invisible(TRUE))

        src_norm <- normalizePath(src_dir, winslash = "/")
        for (src in files) {
          rel  <- sub(paste0("^", src_norm, "/?"), "", normalizePath(src, winslash = "/"))
          dest <- file.path(dest_dir, rel)
          dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
          file.copy(src, dest, overwrite = TRUE)
        }
        invisible(TRUE)
      }

      # 1) package data -> MiamiR Demo Data/data
      copy_tree(pkg_data_dir, file.path(bundle_root, "data"))

      # 2) inst/extdata folders -> MiamiR Demo Data/extdata/<folder>
      copy_tree(file.path(pkg_ext_dir, "Example_Data_Raw"),
                file.path(bundle_root, "extdata", "Example_Data_Raw"))

      copy_tree(file.path(pkg_ext_dir, "Fake_GWAS_Files_Raw"),
                file.path(bundle_root, "extdata", "Fake_GWAS_Files_Raw"))

      # Create archive containing the single folder "MiamiR Demo Data"
      tmp_tar <- tempfile(fileext = ".tar.gz")
      oldwd <- getwd()
      on.exit(setwd(oldwd), add = TRUE)
      setwd(stage_root)

      utils::tar(
        tarfile = tmp_tar,
        files = "MiamiR Demo Data",
        compression = "gzip"
      )

      file.copy(tmp_tar, file, overwrite = TRUE)
    },
    contentType = "application/gzip"
  )

  output$miami_args_simple <- renderUI({
    sa <- formals(.Single_Plot_original)
    keep <- setdiff(names(sa), "Data")
    lapply(
      keep,
      function(arg) make_input_like_meta(
        prefix = "miami_s_",
        name   = arg,
        default = safe_eval_default_reg(sa[[arg]], environment(.Single_Plot_original))
      )
    )
  })

  output$miami_args_advanced <- renderUI({
    sa <- formals(.Single_Plot_original)
    keep <- setdiff(names(sa), "Data")
    tagList(
      h5("Top panel"),
      lapply(
        keep,
        function(arg) make_input_like_meta(
          prefix = "top_",
          name   = arg,
          default = safe_eval_default_reg(sa[[arg]], environment(.Single_Plot_original))
        )
      ),
      tags$hr(),
      h5("Bottom panel"),
      lapply(
        keep,
        function(arg) make_input_like_meta(
          prefix = "bottom_",
          name   = arg,
          default = safe_eval_default_reg(sa[[arg]], environment(.Single_Plot_original))
        )
      )
    )
  })

  output$miami_args_ui_top <- renderUI({
    single_args <- formals(.Single_Plot_original)
    keep <- setdiff(names(single_args), "Data")
    lapply(keep, function(arg){
      make_input_like_meta(
        prefix = "top_",
        name   = arg,
        default = safe_eval_default_reg(single_args[[arg]], environment(.Single_Plot_original))
      )
    })
  })

  output$miami_args_ui_bottom <- renderUI({
    single_args <- formals(.Single_Plot_original)
    keep <- setdiff(names(single_args), "Data")
    lapply(keep, function(arg){
      make_input_like_meta(
        prefix = "bottom_",
        name   = arg,
        default = safe_eval_default_reg(single_args[[arg]], environment(.Single_Plot_original))
      )
    })
  })


  #Map tab

  tab_to_theme <- c(
    home="default", miami="miami", single="single", regional="regional",
    annotate="annotate", meta="meta", forest="forest", model="model",
    tube="tube",
    hla="hla"
  )




  observeEvent(input$tube_interactive_view, {

    cur <- isolate(input$tube_tabs %||% "Single_Plot")

    target_id <- if (identical(cur, "Regional_Plot"))

      "tube_regional_interactive_view"
    else
      "tube_single_interactive_view"

    session$sendCustomMessage("click-el", list(id = target_id))

  })


  observeEvent(input$regional_interactive_view, {

    res <- regional_results()

    if (!length(res)) {

      showNotification("Run Regional_Plot first.", type = "warning", duration = 6)

      return(invisible(NULL))

    }

    curTab <- input$regional_output_tabs

    if (is.null(curTab) || !(curTab %in% names(res))) {

      showNotification("Pick a Regional dataset tab first.", type = "warning", duration = 6)

      return(invisible(NULL))

    }


    idx_input <- paste0("regional_current_slide_", curTab)
    idx <- isolate(input[[idx_input]])

    if (is.null(idx) || !is.finite(idx)) idx <- 1L

    idx <- max(1L, min(idx, length(res[[curTab]]$plots)))

    rp <- res[[curTab]]$plots[[idx]]
    req(!is.null(rp))

    rel_h <- reg_extract_rel_heights(rp)
    dyn_h_in <- as.numeric(attr(rp, "dynamic_height"))

    if (!is.finite(dyn_h_in)) dyn_h_in <- 20

    dyn_h_in <- max(8, min(60, dyn_h_in))

    sp <- reg_split_regional(rp)

    top <- if (inherits(sp$top, "ggplot")) sp$top else ggplotify::as.ggplot(sp$top)

    bot <- if (inherits(sp$bottom, "ggplot")) sp$bottom else ggplotify::as.ggplot(sp$bottom)

    top_info <- reg_find_xy_for_overlay(top)

    top_df   <- if (!is.null(top_info)) top_info$df else NULL

    top_i    <- reg_ensure_hover(top, size_pts = reg_CAPTURE_SIZE_TOP, tag = "top")

    bot_i <- reg_inject_gene_hover(bot, capture_size = reg_CAPTURE_SIZE_BOTTOM)
    gene_cap_df <- attr(bot_i, ".__capture_df")

    top_i <- top_i +
      ggplot2::theme(
        plot.margin  = ggplot2::margin(t = 5, r = 10, b = 0, l = reg_LEFT_PAD_PT, unit = "pt"),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8))
      ) +
      ggplot2::coord_cartesian(clip = "off")

    bot_i <- bot_i +
      ggplot2::theme(
        plot.margin  = ggplot2::margin(t = 0, r = 10, b = 5, l = reg_LEFT_PAD_PT, unit = "pt"),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8))
      ) +
      ggplot2::coord_cartesian(clip = "off")

    combined <- (top_i / bot_i) + patchwork::plot_layout(guides = "keep", heights = rel_h %||% NULL)

    #Build widget

    g <- ggiraph::girafe(
      ggobj      = combined,
      width_svg  = reg_BASE_SVG_WIDTH_IN,
      height_svg = dyn_h_in,
      options = list(
        ggiraph::opts_sizing(rescale = TRUE, width = 1),
        ggiraph::opts_hover(css = "opacity:1; stroke:black; stroke-width:3px; r:6;"),
        ggiraph::opts_tooltip(css = "background:white;color:black;padding:8px;border-radius:6px;font-size:14px;box-shadow:0 2px 6px rgba(0,0,0,0.2);white-space:pre-line;"),
        ggiraph::opts_toolbar(saveaspng = TRUE, position = "topright"),
        ggiraph::opts_zoom(min = 0.5, max = 2),
        ggiraph::opts_selection(type = "single")
      )
    )


    #Build click-detail maps

    TOP_ROWS <- list(); TOP_COLS <- character(0)

    if (!is.null(top_df) && NROW(top_df)) {

      idcol <- intersect(c("ID","SNP","RS","rsid","rsID"), names(top_df))[1]

      ids <- if (!is.na(idcol) && nzchar(idcol)) as.character(top_df[[idcol]]) else sprintf("top-%d", seq_len(NROW(top_df)))

      TOP_COLS <- names(top_df)

      for (i in seq_len(NROW(top_df))) {

        idv <- ids[i]; if (!nzchar(idv)) next

        row <- top_df[i, , drop = FALSE]

        TOP_ROWS[[idv]] <- lapply(row, function(x){ x <- as.character(x); x[is.na(x)] <- ""; paste(x, collapse = ", ") })

      }

    }

    BOT_ROWS <- list(); BOT_COLS <- character(0)

    if (!is.null(gene_cap_df) && NROW(gene_cap_df)) {

      BOT_COLS <- names(gene_cap_df)

      for (i in seq_len(NROW(gene_cap_df))) {

        idv <- as.character(gene_cap_df$.id[i]); if (!nzchar(idv)) next

        row <- gene_cap_df[i, , drop = FALSE]

        BOT_ROWS[[idv]] <- lapply(row, function(x){ x <- as.character(x); x[is.na(x)] <- ""; paste(x, collapse = ", ") })

      }

    }

    #Styles

    g <- htmlwidgets::prependContent(
      g,
      htmltools::tags$style(htmltools::HTML("

      html, body { margin:0; background:#fff; }
      .girafe_container_std { margin:0 !important; padding:0 !important; text-align:left !important; }
      .girafe_container_std svg { display:block; overflow:visible !important; }

      #snp-modal { position: fixed; inset:0; display:none; align-items:center; justify-content:center;
                   background:rgba(0,0,0,.4); z-index:99999; }
      #snp-modal .content {
        position:relative; background:#fff; width:min(720px,92vw);
        border-radius:16px; padding:22px 24px 18px;
        box-shadow:0 18px 40px rgba(0,0,0,.25); border:1px solid #e5e7eb;
      }
      #snp-title { margin:0 0 12px; font-weight:800; font-size:20px; color:#2c3e50; }
      #snp-close { position:absolute; top:10px; right:10px; width:28px; height:28px;
                   border-radius:9999px; border:1px solid #e5e7eb; background:#f3f4f6; color:#374151;
                   font-weight:700; cursor:pointer; }
      #snp-close:hover { background:#e5e7eb; }

      #snp-row { max-height:50vh; overflow:auto; margin:6px 0 12px; }
      table.snp-table { width:100%; border-collapse:collapse; table-layout:fixed; }
      table.snp-table th, table.snp-table td {
        border:1px solid #e5e7eb; padding:6px 8px; vertical-align:top; word-break:break-word;
      }
      table.snp-table th { background:#f9fafb; width:30%; text-align:left; color:#374151; }
      table.snp-table td { background:#ffffff; color:#111827; }

      .snp-actions { display:flex; gap:10px; margin-top:8px; flex-wrap:wrap; }
      .btn-google, .btn-dbsnp, .btn-ncbi, .btn-ncbi-id {
        display:inline-block; text-decoration:none; padding:10px 14px; font-weight:700; color:#fff; border-radius:8px;
      }
      .btn-google { background:#4285F4; }  .btn-google:hover { background:#3367D6; }
      .btn-dbsnp  { background:#28a745; }  .btn-dbsnp:hover  { background:#218838; }
      .btn-ncbi   { background:#0c7bdc; }  .btn-ncbi:hover   { background:#0a67b6; }
      .btn-ncbi-id{ background:#0b5fa5; }  .btn-ncbi-id:hover{ background:#094e88; }

    "))
    )

    g <- htmlwidgets::prependContent(
      g,
      htmltools::tags$script(HTML(sprintf(
        "window.__TOP_ROWS__=%s; window.__TOP_COLS__=%s; window.__BOT_ROWS__=%s; window.__BOT_COLS__=%s;",
        jsonlite::toJSON(TOP_ROWS, auto_unbox=TRUE, null="null", digits=NA),
        jsonlite::toJSON(TOP_COLS, auto_unbox=TRUE, null="null", digits=NA),
        jsonlite::toJSON(BOT_ROWS, auto_unbox=TRUE, null="null", digits=NA),
        jsonlite::toJSON(BOT_COLS, auto_unbox=TRUE, null="null", digits=NA)
      )))
    )

    g <- htmlwidgets::onRender(

      g,
      '
  function(el){
    // Build modal once (same style/structure as standalone Regional)
    var modal = document.getElementById("snp-modal");
    if(!modal){
      modal = document.createElement("div");
      modal.id = "snp-modal";
      modal.innerHTML = `
        <div class="content">
          <button id="snp-close" aria-label="Close">×</button>
          <h2 id="snp-title"></h2>
          <div id="snp-row"></div>
          <div class="snp-actions">
            <a id="btn-google"  class="btn-google" target="_blank" rel="noopener">🔍 Google</a>
            <a id="btn-dbsnp"   class="btn-dbsnp"  target="_blank" rel="noopener">🧬 dbSNP</a>
            <a id="btn-ncbi"    class="btn-ncbi"   target="_blank" rel="noopener" style="display:none;">🔎 NCBI Gene (search)</a>
            <a id="btn-ncbi-id" class="btn-ncbi-id"target="_blank" rel="noopener" style="display:none;">🆔 NCBI Gene (ID)</a>
          </div>
        </div>`;
      document.body.appendChild(modal);
      modal.querySelector("#snp-close").addEventListener("click", function(){ modal.style.display="none"; });
      modal.addEventListener("click", function(e){ if(e.target===modal) modal.style.display="none"; });
      document.addEventListener("keydown", function(e){ if(e.key==="Escape") modal.style.display="none"; });
    }

    var titleEl = modal.querySelector("#snp-title");
    var rowEl   = modal.querySelector("#snp-row");
    var aGoogle = modal.querySelector("#btn-google");
    var aDbSnp  = modal.querySelector("#btn-dbsnp");
    var aNcbi   = modal.querySelector("#btn-ncbi");
    var aNcbiId = modal.querySelector("#btn-ncbi-id");

    var TOP_ROWS = window.__TOP_ROWS__ || {};
    var TOP_COLS = window.__TOP_COLS__ || [];
    var BOT_ROWS = window.__BOT_ROWS__ || {};
    var BOT_COLS = window.__BOT_COLS__ || [];

    function renderTable(obj, order){
      var tbl = document.createElement("table"); tbl.className = "snp-table";
      var tb  = document.createElement("tbody");
      var keys = (Array.isArray(order) && order.length) ? order.slice() : Object.keys(obj || {});

      var chromIdx = keys.findIndex(function(k){ return String(k).toUpperCase() === "CHROM"; });
      if (chromIdx >= 0) {
        var kept = keys.slice(0, chromIdx); // exclude CHROM itself
        if (kept.length === 0) {
          // CHROM was first → keep all non-CHROM keys so table isn’t empty
          kept = keys.filter(function(k){ return String(k).toUpperCase() !== "CHROM"; });
        }
        keys = kept;
      }

      keys.forEach(function(k){
        if (!(k in obj)) return;
        var v = obj[k]; if (v == null) return;
        var vs = String(v); if (!vs) return;

        var tr = document.createElement("tr");
        var th = document.createElement("th"); th.textContent = k;
        var td = document.createElement("td"); td.textContent = vs;
        tr.appendChild(th); tr.appendChild(td);
        tb.appendChild(tr);
      });

      tbl.appendChild(tb);
      return tbl;
    }

    function getFirst(obj, keys){
      for (var i=0;i<keys.length;i++){
        var k = keys[i];
        if (obj && obj[k] != null && String(obj[k]).length) return String(obj[k]);
      }
      return "";
    }
    function onlyDigits(s){ var m = String(s||"").match(/\\d+/); return m?m[0]:""; }

    function openGene(id){
      var row = BOT_ROWS[id]; if(!row){ return; }
      var geneSym = getFirst(row, ["label","gene_symbol","symbol","hgnc_symbol","GeneSymbol","gene"]);
      var geneId  = getFirst(row, ["gene_id","ensembl_gene_id","ENSEMBL","Ensembl","ENSG","ensg","gene_id_version","gene_stable_id","GeneID","geneId"]);
      var entrez  = getFirst(row, ["ENTREZID","entrez","entrez_id","ncbi_gene_id"]);

      titleEl.textContent = geneSym || geneId || "Gene";

      var GENE_KEEP = ["gene_id","tx_start","tx_end","tx_length","strand","gene_biotype","label","start","end"];
      var order = GENE_KEEP.filter(function(k){ return Object.prototype.hasOwnProperty.call(row, k); });

      rowEl.innerHTML = "";
      rowEl.appendChild(renderTable(row, order)); // table builder is shared

      // Buttons (identical behavior)
      var qGoogle = encodeURIComponent(geneSym || geneId || "gene");
      aGoogle.href = "https://www.google.com/search?q=" + qGoogle;
      aDbSnp.style.display = "none"; // hide dbSNP for gene popup

      var qNcbi = encodeURIComponent(geneId || "");
      aNcbi.style.display = "inline-block";
      aNcbi.textContent = "🔎 NCBI Gene (search)";
      aNcbi.href = "https://www.ncbi.nlm.nih.gov/gene/?term=" + (qNcbi || encodeURIComponent(geneSym || "gene"));

      var ncbiId = onlyDigits(entrez);
      if(ncbiId){
        aNcbiId.style.display = "inline-block";
        aNcbiId.href = "https://www.ncbi.nlm.nih.gov/gene/" + ncbiId;
      } else {
        aNcbiId.style.display = "none";
      }

      modal.style.display = "flex";
    }

    function openSnp(id){
      var row = TOP_ROWS[id]; if(!row){ return; }
      titleEl.textContent = "SNP: " + id;

      var chr = getFirst(row, ["CHR","CHROM","Chromosome","chromosome"]);
      var pos = getFirst(row, ["POS","GENPOS","Position","BP","position"]);
      var term = (chr && pos) ? (chr + ":" + pos) : id;

      rowEl.innerHTML = "";
      rowEl.appendChild(renderTable(row, TOP_COLS));

      function looksLikeRS(x){
        if(!x) return false;
        var s=String(x); if (s.length<3) return false;
        if (s.slice(0,2).toLowerCase()!=="rs") return false;
        for (var i=2;i<s.length;i++){ var c=s.charCodeAt(i); if (c<48 || c>57) return false; }
        return true;
      }
      var rs = looksLikeRS(id) ? id : "";
      var googleQ = encodeURIComponent(rs || term || id);
      aGoogle.href = "https://www.google.com/search?q=" + googleQ;

      aDbSnp.style.display = "inline-block";
      aDbSnp.href = rs
        ? ("https://www.ncbi.nlm.nih.gov/snp/" + rs)
        : ("https://www.ncbi.nlm.nih.gov/snp/?term=" + googleQ);

      aNcbi.style.display = "none";
      aNcbiId.style.display = "none";

      modal.style.display = "flex";
    }

    el.addEventListener("click", function(e){
      var node = e.target.closest("[data-id]");
      if(!node) return;
      var id = node.getAttribute("data-id") || "";
      if (id.indexOf("gene-") === 0) openGene(id);
      else openSnp(id);
    });
  }'
    )

    tmphtml <- tempfile(pattern = "regional_girafe_", fileext = ".html")
    htmlwidgets::saveWidget(g, tmphtml, selfcontained = TRUE)
    prefix <- paste0("rview_", as.integer(runif(1, 1e6, 1e9)))
    shiny::addResourcePath(prefix, dirname(tmphtml))
    rel_url <- paste(prefix, basename(tmphtml), sep = "/")
    session$sendCustomMessage("open-window", rel_url)

  })


  observeEvent(input$main_tabs, {

    theme <- unname(tab_to_theme[[input$main_tabs]] %||% "default")
    session$sendCustomMessage("set-theme", theme)

  }, ignoreInit = FALSE)


  #Tile - tab navigation

  nav_map <- c(
    tile_miami="miami", tile_single="single", tile_regional="regional",
    tile_annotate="annotate", tile_meta="meta", tile_forest="forest",
    tile_model="model", tile_tube="tube",
    tile_hla="hla"
  )
  lapply(names(nav_map), function(id) {

    observeEvent(input[[id]], {

      updateTabsetPanel(session, "main_tabs", selected = nav_map[[id]])

    }, ignoreInit = TRUE)

  })

  #Back links - Home
  back_ids <- c(
    "back_home_miami","back_home_single","back_home_regional",
    "back_home_annotate","back_home_meta","back_home_forest",
    "back_home_model","back_home_tube",
    "back_home_hla"
  )


  lapply(back_ids, function(id) {
    observeEvent(input[[id]], {
      updateTabsetPanel(session, "main_tabs", selected = "home")
    }, ignoreInit = TRUE)
  })

  #Tube_Alloys (merged tab)

  tube_data_labels <- reactiveVal(character(0))

  tube_log_txt <- reactiveVal("")
  tube_append_log <- function(...) isolate(tube_log_txt(paste0(tube_log_txt(), format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = ""), "\n")))
  output$tube_log <- renderText(tube_log_txt())

  tube_fun <- reactiveVal(NULL)
  tube_defaults <- reactiveVal(NULL)
  observe({

    f <- get_tube_fun()
    tube_fun(f)
    fmls <- formals(f)
    drops <- c("Data","Phenos","Covars","Phenos_Covars","session","...")
    keep  <- setdiff(names(fmls), drops)
    defs  <- lapply(fmls[keep], function(x) safe_eval_default_any(x, environment(f)))
    names(defs) <- keep
    tube_defaults(defs)

  })

  output$tube_arg_inputs <- renderUI({

    defs <- tube_defaults(); if (is.null(defs)) return(NULL)

    lapply(names(defs), function(nm) make_arg_input("tube_", nm, defs[[nm]]))

  })

  #per-dataset sets for Tube outputs

  tube_single_sets   <- reactiveVal(list())
  tube_regional_sets <- reactiveVal(list())

  tube_key_from_name <- function(nm, kind) {

    nm <- tolower(nm %||% "")
    k  <- sub(sprintf("^(%s)[^A-Za-z0-9]*", kind), "", nm)
    k  <- gsub("[^A-Za-z0-9]+", "_", k)
    k  <- sub("^_+", "", k)

    if (nzchar(k)) k else nm

  }

  tube_label_from_key <- function(k) {

    if (!nzchar(k)) return("Set")

    tools::toTitleCase(gsub("_", " ", k))

  }


  output$tube_download <- downloadHandler(

    filename = function() {

      ext <- tolower(input$tube_ext %||% "png")

      n_single <- !is.null(tube_single_plot())
      n_gwas   <- !is.null(tube_gwas_plot())
      n_forest <- !is.null(tube_forest_plot())
      n_reg    <- length(tube_regionals())

      total <- as.integer(n_single) + as.integer(n_gwas) + as.integer(n_forest) + n_reg

      if (total > 1) return(paste0("Tube_All_", ext, ".zip"))

      sub <- input$tube_tabs %||% "Single_Plot"
      base <- switch(sub,
                     "Single_Plot"   = "Single_Plot",
                     "Regional_Plot" = "Regional_all",
                     "Forest_Plot"   = "Forest",
                     "GWAS_Model"    = "GWAS_Model",
                     "Tube"
      )

      if (sub == "Regional_Plot") paste0(base, "_", ext, ".zip") else paste0(base, ".", ext)

    },

    content = function(file) {

      ext <- tolower(input$tube_ext %||% "png")

      save_one <- function(p, path, w_in, h_in, dpi = 100, ext = "png") {

        if (ext == "pdf") tube_save_plot_pdf(p, path, width_in = w_in, height_in = h_in)

        else              tube_save_plot_png(p, path, width_in = w_in, height_in = h_in, dpi = dpi)

      }

      has_single <- !is.null(tube_single_plot())
      has_gwas   <- !is.null(tube_gwas_plot())
      has_forest <- !is.null(tube_forest_plot())
      regs       <- tube_regionals()
      n_reg      <- length(regs)

      total <- as.integer(has_single) + as.integer(has_gwas) + as.integer(has_forest) + n_reg

      #ZIP if multiple artifacts

      if (total > 1) {

        tmpdir <- tempfile("tube_all_"); dir.create(tmpdir)
        files <- character(0)

        if (has_single) {

          p  <- tube_single_plot()
          h  <- tube_single_hin()
          fn <- file.path(tmpdir, paste0("Single_Plot.", ext))
          save_one(p, fn, input$tube_single_width, h, input$tube_single_dpi, ext); files <- c(files, fn)

        }

        if (n_reg) {

          for (i in seq_along(regs)) {

            p  <- regs[[i]]
            h  <- tube_get_dyn_height(p, fallback = 7)
            fn <- file.path(tmpdir, sprintf("Regional_%02d.%s", i, ext))
            save_one(p, fn, input$tube_regional_width, h, input$tube_regional_dpi, ext); files <- c(files, fn)

          }

        }

        if (has_forest) {

          p  <- tube_forest_plot()
          h  <- tube_forest_hin()
          fn <- file.path(tmpdir, paste0("Forest.", ext))
          save_one(p, fn, input$tube_forest_width, h, input$tube_forest_dpi, ext); files <- c(files, fn)

        }

        if (has_gwas) {

          p  <- tube_gwas_plot()
          h  <- tube_gwas_hin()
          fn <- file.path(tmpdir, paste0("GWAS_Model.", ext))
          save_one(p, fn, input$tube_gwas_width, h, input$tube_gwas_dpi, ext); files <- c(files, fn)

        }

        zip_it(files, file)
        return(invisible())

      }

      #Otherwise - export just the active tab

      sub <- input$tube_tabs %||% "Single_Plot"

      if (sub == "Single_Plot") {

        p <- tube_single_plot(); req(tube_is_plot_like(p))
        h <- tube_single_hin()
        save_one(p, file, input$tube_single_width, h, input$tube_single_dpi, ext)

        return(invisible())

      }

      if (sub == "GWAS_Model") {

        p <- tube_gwas_plot(); req(tube_is_plot_like(p))
        h <- tube_gwas_hin()
        save_one(p, file, input$tube_gwas_width, h, input$tube_gwas_dpi, ext)

        return(invisible())

      }

      if (sub == "Forest_Plot") {

        p <- tube_forest_plot(); req(tube_is_plot_like(p))
        h <- tube_forest_hin()
        save_one(p, file, input$tube_forest_width, h, input$tube_forest_dpi, ext)

        return(invisible())

      }

      if (sub == "Regional_Plot") {

        lst <- tube_regionals(); req(length(lst))
        tmpdir <- tempfile("tube_regional_"); dir.create(tmpdir)

        files <- vapply(seq_along(lst), function(i) {

          p  <- lst[[i]]
          h  <- tube_get_dyn_height(p, fallback = 7)
          fn <- file.path(tmpdir, sprintf("Regional_%02d.%s", i, ext))
          save_one(p, fn, input$tube_regional_width, h, input$tube_regional_dpi, ext)
          fn

        }, character(1))

        zip_it(files, file)

        return(invisible())

      }

      stop("Nothing to download yet — run Tube_Alloys first.")

    }

  )


  tube_assign_file <- function(file) {

    if (is.null(file)) return(NULL)
    obj <- sanitize(tools::file_path_sans_ext(basename(file$name)))
    df  <- vroom::vroom(file$datapath, show_col_types = FALSE)
    assign(obj, df, envir = .GlobalEnv)
    tube_append_log("Assigned object: ", obj)
    obj

  }

  tube_assign_files <- function(files) {

    if (is.null(files)) return(character(0))

    vapply(seq_len(nrow(files)), function(i) {

      obj <- sanitize(tools::file_path_sans_ext(basename(files$name[i])))
      df  <- vroom::vroom(files$datapath[i], show_col_types = FALSE)
      assign(obj, df, envir = .GlobalEnv)
      tube_append_log("Assigned object: ", obj)
      obj
    }, character(1))

  }

  tube_build_call_preview <- reactive({

    nm_data_vec <- if (!is.null(input$tube_data$name))

      sanitize(tools::file_path_sans_ext(input$tube_data$name)) else character(0)

    nm_pheno <- if (!is.null(input$tube_pheno$name))

      sanitize(tools::file_path_sans_ext(input$tube_pheno$name)) else NULL

    nm_covar <- if (!is.null(input$tube_covar$name))

      sanitize(tools::file_path_sans_ext(input$tube_covar$name)) else NULL

    nm_phc   <- if (!is.null(input$tube_pheno_covar$name))

      sanitize(tools::file_path_sans_ext(input$tube_pheno_covar$name)) else NULL

    defs <- tube_defaults() %||% list()
    arg_vals <- lapply(names(defs), function(nm) parse_arg_val(input[[paste0("tube_", nm)]]))
    names(arg_vals) <- names(defs)
    shown <- Filter(Negate(is.null), arg_vals)

    data_str <- if (!length(nm_data_vec)) "Data"

    else paste0("c(", paste(nm_data_vec, collapse = ", "), ")")

    args <- c(sprintf("Data = %s", data_str))

    if (!is.null(nm_pheno)) args <- c(args, sprintf("Phenos = %s", nm_pheno))

    if (!is.null(nm_covar)) args <- c(args, sprintf("Covars = %s", nm_covar))

    if (!is.null(nm_phc))   args <- c(args, sprintf("Phenos_Covars = %s", nm_phc))

    if (length(shown)) {

      for (nm in names(shown)) {

        args <- c(args, sprintf("%s = %s", nm, format_arg_preview(shown[[nm]])))

      }

    }

    paste0("Tube_Alloys(\n  ", paste(args, collapse = ",\n  "), "\n)")

  })

  output$tube_call_preview <- renderText(tube_build_call_preview())

  tube_single_plot  <- reactiveVal(NULL); tube_single_png <- reactiveVal(NULL); tube_single_hin <- reactiveVal(10)
  tube_gwas_plot    <- reactiveVal(NULL); tube_gwas_png   <- reactiveVal(NULL); tube_gwas_hin   <- reactiveVal(8)
  tube_regionals    <- reactiveVal(list()); tube_regionals_pngs <- reactiveVal(character(0))
  tube_forest_plot  <- reactiveVal(NULL); tube_forest_png <- reactiveVal(NULL); tube_forest_hin <- reactiveVal(8)
  tube_data_frames  <- reactiveVal(list())

  tube_classify_kind <- function(nm, val) {

    nm_low <- tolower(nm %||% "")

    if (grepl("gwas",    nm_low)) return("gwas")

    if (grepl("forest",  nm_low)) return("forest")

    if (grepl("regional",nm_low)) return("regional")

    if (grepl("single",  nm_low)) return("single")

    if (tube_is_plot_like(val)) return("single")

    if (is.list(val) && any(vapply(val, tube_is_plot_like, logical(1)))) return("regional")

    "other"

  }

  # Run Tube_Alloys

  observeEvent(input$tube_run, {

    assign("use_wrapper", TRUE, envir = .GlobalEnv)
    reload_miamir()

    f <- get_tube_fun()
    tube_fun(f)

    fmls  <- formals(f)
    drops <- c("Data","Phenos","Covars","Phenos_Covars","session","...")
    keep  <- setdiff(names(fmls), drops)
    defs  <- lapply(fmls[keep], function(x) safe_eval_default_any(x, environment(f)))
    names(defs) <- keep
    tube_defaults(defs)

    req(input$tube_data)

    tube_append_log("Reading and assigning inputs to .GlobalEnv...")

    nm_data_vec <- tube_assign_files(input$tube_data)
    nm_pheno    <- tube_assign_file(input$tube_pheno)
    nm_covar    <- tube_assign_file(input$tube_covar)
    nm_phc      <- tube_assign_file(input$tube_pheno_covar)

    tube_data_labels(nm_data_vec)


    arg_vals <- lapply(names(defs), function(nm) parse_arg_val(input[[paste0("tube_", nm)]]))

    names(arg_vals) <- names(defs)

    assign(".NF_Orig_Names", nm_data_vec, envir = .GlobalEnv)

    on.exit(try(rm(".NF_Orig_Names", envir = .GlobalEnv), silent = TRUE), add = TRUE)

    data_expr <- as.call(c(list(as.name("c")), lapply(nm_data_vec, as.name)))

    call_list <- c(
      list(as.name("Tube_Alloys"), Data = data_expr),

      if (!is.null(nm_pheno)) list(Phenos = as.name(nm_pheno)) else NULL,

      if (!is.null(nm_covar)) list(Covars = as.name(nm_covar)) else NULL,

      if (!is.null(nm_phc))   list(Phenos_Covars = as.name(nm_phc)) else NULL,

      arg_vals
    )

    ca <- as.call(call_list)

    tube_append_log("Calling:\n", tube_build_call_preview())

    out <- tryCatch(

      eval(ca, envir = .GlobalEnv),

      error = function(e) {

        tube_append_log("Tube_Alloys error: ", e$message)
        showNotification(paste("Tube_Alloys error:", e$message), type = "error", duration = 8)
        NULL
      },

      finally = { try(rm(".NF_Orig_Names", envir = .GlobalEnv), silent = TRUE) }

    )

    req(!is.null(out))

    tube_data_frames(Filter(is.data.frame, out))

    tube_single_sets(list())
    tube_regional_sets(list())
    tube_regionals(list())

    tube_single_plot(NULL); tube_single_png(NULL); tube_single_hin(input$tube_single_height %||% 10)
    tube_gwas_plot(NULL);   tube_gwas_png(NULL);   tube_gwas_hin(8)
    tube_forest_plot(NULL); tube_forest_png(NULL); tube_forest_hin(8)

    tube_extract_plots <- function(x) {


      if (tube_is_plot_like(x)) return(list(x))

      if (!is.list(x)) return(list())

      out <- list()

      for (e in x) {

        if (tube_is_plot_like(e)) out <- c(out, list(e))

        else if (is.list(e)) {

          kids <- Filter(tube_is_plot_like, e)

          if (length(kids)) out <- c(out, kids)

        }

      }

      out

    }

    tube_group_plots <- function(x) {


      if (!is.list(x)) return(list(default = tube_extract_plots(x)))


      if (all(vapply(x, tube_is_plot_like, logical(1)))) {

        return(list(default = x))

      }

      groups <- list()
      nms <- names(x)

      for (i in seq_along(x)) {

        grp <- tube_extract_plots(x[[i]])

        if (length(grp)) {

          nm <- nms[i]
          if (is.null(nm) || !nzchar(nm)) nm <- paste0("set_", i)
          groups[[nm]] <- grp

        }

      }

      if (!length(groups)) list(default = tube_extract_plots(x)) else groups

    }

    tube_make_unique_key <- function(existing, key) {

      if (key %in% existing) {

        make.unique(c(existing, key))[length(existing) + 1]

      } else key

    }


    for (nm in names(out)) {

      val  <- out[[nm]]
      nmlo <- tolower(nm %||% "")

      kind <- if (grepl("gwas",    nmlo)) "gwas"

      else if (grepl("forest",   nmlo)) "forest"

      else if (grepl("regional", nmlo)) "regional"

      else if (grepl("single",   nmlo) || tube_is_plot_like(val)) "single"

      else "other"

      if (kind == "single") {

        ds_labels <- tube_data_labels()

        sets      <- tube_single_sets()

        make_unique <- function(base, existing) {

          k <- base; i <- 2
          while (k %in% existing) { k <- paste0(base, "_", i); i <- i + 1 }
          k

        }

        add_one <- function(obj, j_idx = 1) {

          if (!tube_is_plot_like(obj)) return(invisible(NULL))

          h   <- tube_get_dyn_height(obj, fallback = 10)
          png <- tempfile(fileext = ".png")
          tube_save_plot_png(
            obj, png,
            width_in  = input$tube_single_width,
            height_in = input$tube_single_height %||% h,
            dpi       = input$tube_single_dpi
          )

          base_lab <- if (j_idx <= length(ds_labels)) ds_labels[[j_idx]] else paste0("Set ", j_idx)

          pretty   <- tools::toTitleCase(gsub("_", " ", base_lab))
          base_key <- paste0("single_", sanitize(base_lab))
          key      <- make_unique(base_key, names(sets))

          sets[[key]] <<- list(plot = obj, png = png, hin = input$tube_single_height %||% h, label = pretty)
        }

        if (tube_is_plot_like(val)) {

          add_one(val, 1)


        } else if (is.list(val) && length(val)) {


          is_plot <- vapply(val, tube_is_plot_like, logical(1))

          if (any(is_plot)) {

            for (j in which(is_plot)) add_one(val[[j]], j)

          } else {

            for (j in seq_along(val)) {

              grp <- val[[j]]

              if (is.list(grp) && length(grp)) {

                k <- which(vapply(grp, tube_is_plot_like, logical(1)))

                if (length(k)) add_one(grp[[k[1]]], j)

              }

            }

          }

        }

        tube_single_sets(sets)


        if (is.null(tube_single_plot()) && length(sets)) {

          first_key <- names(sets)[1]
          tube_single_plot(sets[[first_key]]$plot)
          tube_single_png(sets[[first_key]]$png)
          tube_single_hin(sets[[first_key]]$hin)

        }

        next

      }


      #REGIONAL

      if (kind == "regional") {


        groups <- if (is.list(val) && length(val)) val else list(val)

        rsets <- tube_regional_sets()
        flat  <- tube_regionals()
        ds_labels <- tube_data_labels()

        make_unique <- function(base, existing) {

          k <- base; i <- 2

          while (k %in% existing) { k <- paste0(base, "_", i); i <- i + 1 }

          k

        }

        for (j in seq_along(groups)) {

          grp <- groups[[j]]

          keep <- if (is.list(grp)) Filter(tube_is_plot_like, grp) else if (tube_is_plot_like(grp)) list(grp) else list()

          if (!length(keep)) next

          base_label <- if (j <= length(ds_labels)) ds_labels[[j]] else paste0("Set ", j)

          pretty     <- tools::toTitleCase(gsub("_", " ", base_label))

          base_key <- paste0("regional_", sanitize(base_label))
          key      <- make_unique(base_key, names(rsets))

          pngs <- vapply(seq_along(keep), function(i) {

            p  <- keep[[i]]
            h  <- tube_get_dyn_height(p, fallback = 7)
            fn <- tempfile(fileext = sprintf("_regional_%02d.png", i))
            tube_save_plot_png(
              p, fn,
              width_in  = input$tube_regional_width,
              height_in = h,
              dpi       = input$tube_regional_dpi
            )
            fn
          }, character(1))

          rsets[[key]] <- list(plots = keep, pngs = pngs, label = pretty)
          flat <- c(flat, keep)

        }

        tube_regional_sets(rsets)
        tube_regionals(flat)

        next

      }

      #FOREST

      if (kind == "forest") {

        candidate <- NULL

        if (tube_is_plot_like(val)) candidate <- val

        if (is.null(candidate) && is.list(val)) {

          idx <- which(vapply(val, tube_is_plot_like, logical(1)))

          if (length(idx)) candidate <- val[[idx[1]]]

        }

        if (!is.null(candidate)) {

          h <- tube_get_dyn_height(candidate, fallback = 8)
          tmp <- tempfile(fileext = ".png")
          tube_save_plot_png(
            candidate, tmp,
            width_in = input$tube_forest_width, height_in = h, dpi = input$tube_forest_dpi
          )
          tube_forest_plot(candidate); tube_forest_png(tmp); tube_forest_hin(h)

        }

        next

      }

      #GWAS MODEL

      if (kind == "gwas") {

        candidate <- NULL

        if (tube_is_plot_like(val)) candidate <- val

        if (is.null(candidate) && is.list(val)) {

          idx <- which(vapply(val, tube_is_plot_like, logical(1)))

          if (length(idx)) candidate <- val[[idx[1]]]

        }

        if (!is.null(candidate)) {

          h <- tube_get_dyn_height(candidate, fallback = 8)
          tmp <- tempfile(fileext = ".png")
          tube_save_plot_png(
            candidate, tmp,
            width_in = input$tube_gwas_width, height_in = h, dpi = input$tube_gwas_dpi
          )
          tube_gwas_plot(candidate); tube_gwas_png(tmp); tube_gwas_hin(h)
        }

        next

      }

    }

    tube_append_log("Done.")
  })

  #Tube_Alloys - Single tab UI

  output$tube_single_ui <- renderUI({

    sets <- tube_single_sets()

    tabs <- lapply(names(sets), function(k) {

      imgId <- paste0("tube_single_img_", k)

      local({

        kk <- k
        output[[imgId]] <- renderImage({
          list(src = tube_single_sets()[[kk]]$png, contentType = "image/png", width = "100%")
        }, deleteFile = FALSE)

      })

      tabPanel(title = sets[[k]]$label, value = k, imageOutput(imgId, height = "auto"))

    })

    do.call(tabsetPanel, c(list(id = "tube_single_tabs"), tabs))

  })

  output$tube_regional_ui <- renderUI({

    rsets <- tube_regional_sets()

    tabs <- lapply(names(rsets), function(k) {

      carId <- paste0("tube_regional_carousel_", k)

      local({

        kk <- k

        output[[carId]] <- renderSlickR({

          imgs <- tube_regional_sets()[[kk]]$pngs

          if (!length(imgs)) return(NULL)

          w <- slickR(imgs, slideType = "img") +
            settings(dots = TRUE, infinite = TRUE, arrows = TRUE, adaptiveHeight = FALSE)

          htmlwidgets::onRender(

            w,
            sprintf(

              "function(el,x){
               var $s = $(el).find('.slick-slider');
               // initialize
               Shiny.setInputValue('%1$s', 1, {priority:'event'});
               $s.on('afterChange', function(e, slick, current){
                 Shiny.setInputValue('%1$s', current+1, {priority:'event'});
               });
             }",

              paste0('tube_regional_current_slide_', kk)

            )

          )

        })

      })

      tabPanel(title = rsets[[k]]$label, value = k,
               slickROutput(carId, width = "100%", height = "680px"))

    })

    do.call(tabsetPanel, c(list(id = "tube_regional_tabs"), tabs))

  })


  observeEvent(input$tube_single_interactive_view, {

    sets <- tube_single_sets()
    validate(need(length(sets) > 0, "Run Tube_Alloys first (Single_Plot output required)."))

    cur <- isolate(input$tube_single_tabs)

    if (is.null(cur) || !(cur %in% names(sets))) cur <- names(sets)[1]

    p <- sets[[cur]]$plot
    validate(need(!is.null(p), "No plot found for the selected Tube_Alloys Single tab."))


    p_i <- tryCatch(

      hover_overlay_from_plot(p, prefix = "tube:"),
      error = function(e) {

        showNotification(paste("Interactive overlay failed:", e$message), type = "error", duration = 8)
        NULL

      }

    )

    validate(need(!is.null(p_i), "Unable to build interactive overlay."))

    rows_map  <- attr(p_i, "rows_map")
    col_order <- attr(p_i, "cols")

    g <- ggiraph::girafe(
      ggobj      = p_i,
      width_svg  = 30,
      height_svg = tube_single_hin(),
      options = list(
        ggiraph::opts_sizing(rescale = TRUE, width = 1),
        ggiraph::opts_hover(css = "opacity:1; stroke:black; stroke-width:3px; r:6;"),
        ggiraph::opts_selection(type = "single"),
        ggiraph::opts_tooltip(css = "background:white;color:black;padding:8px;border-radius:6px;font-size:14px;box-shadow:0 2px 6px rgba(0,0,0,0.2);white-space:pre-line;"),
        ggiraph::opts_toolbar(saveaspng = TRUE, position = "topright"),
        ggiraph::opts_zoom(min = 0.5, max = 2)
      )
    )

    g <- htmlwidgets::prependContent(
      g,
      htmltools::tags$style(htmltools::HTML("

      html, body { margin:0; background:#fff; }
      .girafe_container_std { margin:0 !important; padding:0 !important; text-align:left !important; }
      .girafe_container_std svg { display:block; overflow:visible !important; }

      #snp-modal { position: fixed; inset:0; display:none; align-items:center; justify-content:center;
                   background:rgba(0,0,0,.4); z-index:99999; }
      #snp-modal .content {
        position:relative; background:#fff; width:min(720px,92vw);
        border-radius:16px; padding:22px 24px 18px;
        box-shadow:0 18px 40px rgba(0,0,0,.25); border:1px solid #e5e7eb;
      }
      #snp-title { margin:0 0 12px; font-weight:800; font-size:20px; color:#2c3e50; }
      #snp-close { position:absolute; top:10px; right:10px; width:28px; height:28px;
                   border-radius:9999px; border:1px solid #e5e7eb; background:#f3f4f6; color:#374151;
                   font-weight:700; cursor:pointer; }
      #snp-close:hover { background:#e5e7eb; }

      #snp-row { max-height:50vh; overflow:auto; margin:6px 0 12px; }
      table.snp-table { width:100%; border-collapse:collapse; table-layout:fixed; }
      table.snp-table th, table.snp-table td {
        border:1px solid #e5e7eb; padding:6px 8px; vertical-align:top; word-break:break-word;
      }
      table.snp-table th { background:#f9fafb; width:30%; text-align:left; color:#374151; }
      table.snp-table td { background:#ffffff; color:#111827; }

      .snp-actions { display:flex; gap:10px; margin-top:8px; flex-wrap:wrap; }
      .btn-google, .btn-dbsnp {
        display:inline-block; text-decoration:none; padding:10px 14px; font-weight:700; color:#fff; border-radius:8px;
      }
      .btn-google { background:#4285F4; } .btn-google:hover { background:#3367D6; }
      .btn-dbsnp  { background:#28a745; } .btn-dbsnp:hover  { background:#218838; }
    "))
    )

    g <- htmlwidgets::onRender(
      g,
      '
    function(el, x, data){

      // Build modal once

      var modal = document.getElementById("snp-modal");
      if(!modal){
        modal = document.createElement("div");
        modal.id = "snp-modal";
        modal.innerHTML = `
          <div class="content">
            <button id="snp-close" aria-label="Close">×</button>
            <h2 id="snp-title"></h2>
            <div id="snp-row"></div>
            <div class="snp-actions">
              <a id="snp-google" class="btn-google" target="_blank" rel="noopener">🔍 Google</a>
              <a id="snp-dbsnp"  class="btn-dbsnp"  target="_blank" rel="noopener">🧬 dbSNP</a>
            </div>
          </div>`;
        document.body.appendChild(modal);
        modal.querySelector("#snp-close").addEventListener("click", function(){ modal.style.display="none"; });
        modal.addEventListener("click", function(e){ if(e.target===modal) modal.style.display="none"; });
        document.addEventListener("keydown", function(e){ if(e.key==="Escape") modal.style.display="none"; });
      }

      var titleEl = modal.querySelector("#snp-title");
      var rowEl   = modal.querySelector("#snp-row");
      var aGo     = modal.querySelector("#snp-google");
      var aDb     = modal.querySelector("#snp-dbsnp");

      var MAP     = (data && data.rows_map)  || {};
      var CORDER  = (data && data.col_order) || [];

      // EXACT same helpers as standalone Single tab
      function looksLikeRS(x){
        if(!x) return false;
        var s = String(x);
        if (s.length < 3) return false;
        var p = s.slice(0,2).toLowerCase();
        if (p !== "rs") return false;
        for (var i = 2; i < s.length; i++){
          var c = s.charCodeAt(i);
          if (c < 48 || c > 57) return false; // 0-9
        }
        return true;

      }

function renderTable(obj, order){

  var tbl = document.createElement("table"); tbl.className = "snp-table";
  var tb  = document.createElement("tbody");

  // Preserve provided order; otherwise discover keys

  var keys = (Array.isArray(order) && order.length) ? order.slice() : Object.keys(obj || {});


  var chromIdx = keys.findIndex(function(k){
    return String(k).toUpperCase() === "CHROM";

  });

  if (chromIdx >= 0) {

    var kept = keys.slice(0, chromIdx); // exclude CHROM itself

    if (kept.length === 0) {

      kept = keys.filter(function(k){ return String(k).toUpperCase() !== "CHROM"; });
    }

  keys = kept;

}

keys.forEach(function(k){

  if (!(k in obj)) return;
  var v = obj[k]; if (v == null) return;
  var vs = String(v); if (!vs) return;

  var tr = document.createElement("tr");
  var th = document.createElement("th"); th.textContent = k;
  var td = document.createElement("td"); td.textContent = vs;
  tr.appendChild(th); tr.appendChild(td);
  tb.appendChild(tr);

});

tbl.appendChild(tb);
return tbl; // return to caller
}

      function openSnp(id){

        var row = MAP[id]; if(!row) return;
        titleEl.textContent = "SNP: " + id;

        // CHR:POS search term if present

        var chr = row.CHR || row.CHROM || row.Chromosome || row.chromosome || "";
        var pos = row.POS || row.GENPOS || row.Position || row.BP || row.position || "";
        var term = (chr && pos) ? (String(chr) + ":" + String(pos)) : id;

        // Table
        rowEl.innerHTML = "";
        rowEl.appendChild(renderTable(row, CORDER));

        // Links

        var rs = looksLikeRS(id) ? id : "";
        var q  = encodeURIComponent(rs || term || id);
        aGo.href = "https://www.google.com/search?q=" + q;
        aDb.href = rs
          ? ("https://www.ncbi.nlm.nih.gov/snp/" + rs)
          : ("https://www.ncbi.nlm.nih.gov/snp/?term=" + q);

        modal.style.display = "flex";

      }

      // Click → open

      el.addEventListener("click", function(e){
        var node = e.target.closest("[data-id]");
        if(!node) return;
        var id = node.getAttribute("data-id") || "";
        openSnp(id);

      });

    }',

data = list(rows_map = rows_map, col_order = col_order)

    )

    tmphtml <- tempfile(pattern = "tube_single_girafe_", fileext = ".html")
    htmlwidgets::saveWidget(g, tmphtml, selfcontained = TRUE)
    prefix <- paste0("rview_", as.integer(runif(1, 1e6, 1e9)))
    shiny::addResourcePath(prefix, dirname(tmphtml))
    rel_url <- paste(prefix, basename(tmphtml), sep = "/")
    session$sendCustomMessage("open-window", rel_url)


  })


  observeEvent(input$tube_regional_interactive_view, {

    rsets <- tube_regional_sets(); req(length(rsets))

    curTab <- input$tube_regional_tabs

    if (is.null(curTab) || !(curTab %in% names(rsets))) curTab <- names(rsets)[1]

    idx_input <- paste0("tube_regional_current_slide_", curTab)
    idx <- isolate(input[[idx_input]]); if (is.null(idx) || !is.finite(idx)) idx <- 1L

    rp <- rsets[[curTab]]$plots[[idx]]; req(!is.null(rp))

    rel_h <- reg_extract_rel_heights(rp)
    dyn_h_in <- as.numeric(attr(rp, "dynamic_height")); if (!is.finite(dyn_h_in)) dyn_h_in <- 20

    dyn_h_in <- max(8, min(60, dyn_h_in))

    sp  <- reg_split_regional(rp)

    top <- if (inherits(sp$top, "ggplot")) sp$top else ggplotify::as.ggplot(sp$top)

    bot <- if (inherits(sp$bottom, "ggplot")) sp$bottom else ggplotify::as.ggplot(sp$bottom)

    top_info <- reg_find_xy_for_overlay(top)

    top_df   <- if (!is.null(top_info)) top_info$df else NULL

    top_i    <- reg_ensure_hover(top, size_pts = reg_CAPTURE_SIZE_TOP, tag = "top")

    bot_i <- reg_inject_gene_hover(bot, capture_size = reg_CAPTURE_SIZE_BOTTOM)

    gene_cap_df <- attr(bot_i, ".__capture_df")

    top_i <- top_i +
      ggplot2::theme(
        plot.margin  = ggplot2::margin(t = 5, r = 10, b = 0, l = reg_LEFT_PAD_PT, unit = "pt"),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8))
      ) + ggplot2::coord_cartesian(clip = "off")
    bot_i <- bot_i +
      ggplot2::theme(
        plot.margin  = ggplot2::margin(t = 0, r = 10, b = 5, l = reg_LEFT_PAD_PT, unit = "pt"),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8))
      ) + ggplot2::coord_cartesian(clip = "off")

    combined <- (top_i / bot_i) + patchwork::plot_layout(guides = "keep", heights = rel_h %||% NULL)

    g <- ggiraph::girafe(
      ggobj      = combined,
      width_svg  = reg_BASE_SVG_WIDTH_IN,
      height_svg = dyn_h_in,
      options = list(
        ggiraph::opts_sizing(rescale = TRUE, width = 1),
        ggiraph::opts_hover(css = "opacity:1; stroke:black; stroke-width:3px; r:6;"),
        ggiraph::opts_tooltip(css = "background:white;color:black;padding:8px;border-radius:6px;font-size:14px;box-shadow:0 2px 6px rgba(0,0,0,0.2);white-space:pre-line;"),
        ggiraph::opts_toolbar(saveaspng = TRUE, position = "topright"),
        ggiraph::opts_zoom(min = 0.5, max = 2),
        ggiraph::opts_selection(type = "single")
      )
    )

    TOP_ROWS <- list(); TOP_COLS <- character(0)

    if (!is.null(top_df) && NROW(top_df)) {

      idcol <- intersect(c("ID","SNP","RS","rsid","rsID"), names(top_df))[1]

      ids <- if (!is.na(idcol) && nzchar(idcol)) as.character(top_df[[idcol]]) else sprintf("top-%d", seq_len(NROW(top_df)))

      TOP_COLS <- {

        cols <- names(top_df)

        idx  <- match("CHROM", cols, nomatch = 0L)

        if (idx > 0L) cols[seq_len(idx - 1L)] else cols

      }

      for (i in seq_len(NROW(top_df))) {

        idv <- ids[i]; if (!nzchar(idv)) next

        row <- top_df[i, , drop = FALSE]

        TOP_ROWS[[idv]] <- lapply(row, function(x){ x <- as.character(x); x[is.na(x)] <- ""; paste(x, collapse = ", ") })

      }

    }

    BOT_ROWS <- list(); BOT_COLS <- character(0)

    if (!is.null(gene_cap_df) && NROW(gene_cap_df)) {

      cols <- names(gene_cap_df)

      idx <- match("end", cols, nomatch = 0L)

      BOT_COLS <- if (idx > 0L) cols[seq_len(idx)] else cols

      BOT_COLS <- setdiff(BOT_COLS, c("y", "transcript_id"))

      for (i in seq_len(nrow(gene_cap_df))) {

        idv <- as.character(gene_cap_df$.id[i]); if (!nzchar(idv)) next


        row_vis <- gene_cap_df[i, BOT_COLS, drop = FALSE]

        BOT_ROWS[[idv]] <- lapply(row_vis, function(x){

          x <- as.character(x); x[is.na(x)] <- ""; paste(x, collapse = ", ")

        })

      }

    }


    g <- htmlwidgets::prependContent(

      g,
      htmltools::tags$style(htmltools::HTML("
      html, body { margin:0; background:#fff; }
      .girafe_container_std { margin:0 !important; padding:0 !important; text-align:left !important; }
      .girafe_container_std svg { display:block; overflow:visible !important; }
      #snp-modal { position: fixed; inset:0; display:none; align-items:center; justify-content:center; background:rgba(0,0,0,.4); z-index:99999; }
      #snp-modal .content { position:relative; background:#fff; width:min(720px,92vw); border-radius:16px; padding:22px 24px 18px; box-shadow:0 18px 40px rgba(0,0,0,.25); border:1px solid #e5e7eb; }
      #snp-title { margin:0 0 12px; font-weight:800; font-size:20px; color:#2c3e50; }
      #snp-close { position:absolute; top:10px; right:10px; width:28px; height:28px; border-radius:9999px; border:1px solid #e5e7eb; background:#f3f4f6; color:#374151; font-weight:700; cursor:pointer; }
      #snp-close:hover { background:#e5e7eb; }
      #snp-row { max-height:50vh; overflow:auto; margin:6px 0 12px; }
      table.snp-table { width:100%; border-collapse:collapse; table-layout:fixed; }
      table.snp-table th, table.snp-table td { border:1px solid #e5e7eb; padding:6px 8px; vertical-align:top; word-break:break-word; }
      table.snp-table th { background:#f9fafb; width:30%; text-align:left; color:#374151; }
      table.snp-table td { background:#ffffff; color:#111827; }
      .snp-actions { display:flex; gap:10px; margin-top:8px; flex-wrap:wrap; }
      .btn-google, .btn-dbsnp, .btn-ncbi, .btn-ncbi-id {
        display:inline-block; text-decoration:none; padding:10px 14px; font-weight:700; color:#fff; border-radius:8px;
      }
      .btn-google { background:#4285F4; }  .btn-google:hover { background:#3367D6; }
      .btn-dbsnp  { background:#28a745; }  .btn-dbsnp:hover  { background:#218838; }
      .btn-ncbi   { background:#0c7bdc; }  .btn-ncbi:hover   { background:#0a67b6; }
      .btn-ncbi-id{ background:#0b5fa5; }  .btn-ncbi-id:hover{ background:#094e88; }
    "))

    )

    g <- htmlwidgets::prependContent(

      g,
      htmltools::tags$script(HTML(sprintf(
        "window.__TOP_ROWS__=%s; window.__TOP_COLS__=%s; window.__BOT_ROWS__=%s; window.__BOT_COLS__=%s;",

        jsonlite::toJSON(TOP_ROWS, auto_unbox=TRUE, null='null', digits=NA),
        jsonlite::toJSON(TOP_COLS, auto_unbox=TRUE, null='null', digits=NA),
        jsonlite::toJSON(BOT_ROWS, auto_unbox=TRUE, null='null', digits=NA),
        jsonlite::toJSON(BOT_COLS, auto_unbox=TRUE, null='null', digits=NA)

      )))
    )

    g <- htmlwidgets::onRender(

      g,
      "
      function(el){

        // Build modal once

        var modal = document.getElementById('snp-modal');
        if(!modal){
          modal = document.createElement('div');
          modal.id = 'snp-modal';
          modal.innerHTML = `
            <div class='content'>
              <button id='snp-close' aria-label='Close'>×</button>
              <h2 id='snp-title'></h2>
              <div id='snp-row'></div>
              <div class='snp-actions'>
                <a id='btn-google' class='btn-google' target='_blank' rel='noopener'>🔍 Google</a>
                <a id='btn-dbsnp'  class='btn-dbsnp'  target='_blank' rel='noopener'>🧬 dbSNP</a>
                <a id='btn-ncbi'   class='btn-ncbi'   target='_blank' rel='noopener' style='display:none;'>🔎 NCBI Gene (search)</a>
                <a id='btn-ncbi-id' class='btn-ncbi-id' target='_blank' rel='noopener' style='display:none;'>🆔 NCBI Gene (ID)</a>
              </div>
            </div>`;
          document.body.appendChild(modal);
          modal.querySelector('#snp-close').addEventListener('click', function(){ modal.style.display='none'; });
          modal.addEventListener('click', function(e){ if(e.target===modal) modal.style.display='none'; });
          document.addEventListener('keydown', function(e){ if(e.key==='Escape') modal.style.display='none'; });
        }

        var titleEl = modal.querySelector('#snp-title');
        var rowEl   = modal.querySelector('#snp-row');
        var aGoogle = modal.querySelector('#btn-google');
        var aDbSnp  = modal.querySelector('#btn-dbsnp');
        var aNcbi   = modal.querySelector('#btn-ncbi');
        var aNcbiId = modal.querySelector('#btn-ncbi-id');

        var TOP_ROWS = window.__TOP_ROWS__ || {};
        var TOP_COLS = window.__TOP_COLS__ || [];
        var BOT_ROWS = window.__BOT_ROWS__ || {};
        var BOT_COLS = window.__BOT_COLS__ || [];

        function renderTable(obj, order){

          rowEl.innerHTML = '';
          if(!obj) return;
          var tbl=document.createElement('table'); tbl.className='snp-table';
          var tb =document.createElement('tbody');
          var keys=(order && order.length)?order.slice():Object.keys(obj||{}).sort();
          keys.forEach(function(k){
            if(!(k in obj)) return;
            var v=obj[k]; if(v==null) return;
            var vs=String(v); if(!vs) return;
            var tr=document.createElement('tr');
            var th=document.createElement('th'); th.textContent=k;
            var td=document.createElement('td'); td.textContent=vs;
            tr.appendChild(th); tr.appendChild(td); tb.appendChild(tr);
          });
          tbl.appendChild(tb); rowEl.appendChild(tbl);

        }

        function looksLikeRS(x){
          if(!x) return false; var s=String(x); if(s.length<3) return false;
          if(s.slice(0,2).toLowerCase()!=='rs') return false;
          for(var i=2;i<s.length;i++){ var c=s.charCodeAt(i); if(c<48||c>57) return false; }
          return true;
        }
        function getFirst(obj, keys){ for(var i=0;i<keys.length;i++){ var k=keys[i]; if(obj && obj[k]) return String(obj[k]); } return ''; }
        function onlyDigits(s){ return (String(s||'').match(/\\d+/)||[''])[0]; }

        function openGene(id){
          var row = BOT_ROWS[id]; if(!row){ return; }
          var geneSym = getFirst(row, ['label','gene_symbol','symbol','hgnc_symbol','GeneSymbol','gene']);
          var geneId  = getFirst(row, ['gene_id','ensembl_gene_id','ENSEMBL','Ensembl','ENSG','ensg','gene_id_version','gene_stable_id','GeneID','geneId']);
          var entrez  = getFirst(row, ['ENTREZID','entrez','entrez_id','ncbi_gene_id']);
          titleEl.textContent = geneSym || geneId || 'Gene';

          // Keep important cols first
          var order = BOT_COLS.slice();
          var keep = ['gene_id','tx_start','tx_end','tx_length','strand','gene_biotype','label','start','end'];
          order = keep.filter(function(k){ return order.indexOf(k) >= 0; }).concat(order.filter(function(k){ return keep.indexOf(k) < 0; }));
          renderTable(row, order);

          aDbSnp.style.display = 'none';
          var qGoogle = encodeURIComponent(geneSym || geneId || 'gene');
          aGoogle.href = 'https://www.google.com/search?q=' + qGoogle;

          var qNcbi = encodeURIComponent(geneId || '');
          aNcbi.style.display = 'inline-block';
          aNcbi.textContent   = '🔎 NCBI Gene (search)';
          aNcbi.href          = 'https://www.ncbi.nlm.nih.gov/gene/?term=' + (qNcbi || encodeURIComponent(geneSym || 'gene'));

          var ncbiId = onlyDigits(entrez);
          if(ncbiId){
            aNcbiId.style.display = 'inline-block';
            aNcbiId.href          = 'https://www.ncbi.nlm.nih.gov/gene/' + ncbiId;
          } else {
            aNcbiId.style.display = 'none';
          }
          modal.style.display = 'flex';
        }

        function openSnp(id){
          var row = TOP_ROWS[id]; if(!row){ return; }
          titleEl.textContent = looksLikeRS(id) ? ('SNP: '+id) : (row.SNP||row.ID||id||'Variant');
          renderTable(row, TOP_COLS);

          var chr = row.CHR || row.CHROM || row.Chromosome || '';
          var pos = row.POS || row.GENPOS || row.BP || row.position || '';
          var term = (chr && pos) ? (chr+':'+pos) : id;
          var q    = encodeURIComponent(looksLikeRS(id)?id:term);
          aGoogle.href = 'https://www.google.com/search?q=' + q;
          aDbSnp.style.display = 'inline-block';
          aDbSnp.href = looksLikeRS(id) ? ('https://www.ncbi.nlm.nih.gov/snp/' + id)
                                        : ('https://www.ncbi.nlm.nih.gov/snp/?term=' + q);
          aNcbi.style.display = 'none';
          aNcbiId.style.display = 'none';
          modal.style.display = 'flex';
        }

        el.addEventListener('click', function(e){
          var node = e.target.closest('[data-id]'); if(!node) return;
          var id = node.getAttribute('data-id') || '';
          if(id.indexOf('gene-') === 0) openGene(id);
          else openSnp(id);
        });
      }
    "
    )


    tmphtml <- tempfile(pattern = "tube_regional_girafe_", fileext = ".html")
    htmlwidgets::saveWidget(g, tmphtml, selfcontained = TRUE)
    prefix <- paste0("rview_", as.integer(runif(1, 1e6, 1e9)))
    shiny::addResourcePath(prefix, dirname(tmphtml))
    session$sendCustomMessage("open-window", paste0("/", prefix, "/", basename(tmphtml)))

  })

  #Forest_Plot

  output$tube_forest_img <- renderImage({

    req(tube_forest_png())
    list(src = tube_forest_png(), contentType = "image/png", width = "100%")

  }, deleteFile = FALSE)

  #GWAS_Model

  output$tube_gwas_img <- renderImage({

    req(tube_gwas_png())
    list(src = tube_gwas_png(), contentType = "image/png", width = "100%")

  }, deleteFile = FALSE)

  output$meta_download_simple <- downloadHandler(

    filename = function() {

      ext <- tolower(input$meta_ext %||% "csv")
      paste0("METASOFT_File_Gen_output.", ext)

    },
    content = function(file) {

      df  <- meta_results(); req(!is.null(df))
      ext <- tolower(input$meta_ext %||% "csv")

      if (ext == "csv") {

        utils::write.csv(df, file, row.names = FALSE, na = "")

      } else if (ext %in% c("tsv", "txt")) {

        utils::write.table(df, file, sep = "\t", row.names = FALSE,
                           col.names = TRUE, quote = FALSE, na = "")

      } else {


        utils::write.csv(df, file, row.names = FALSE, na = "")

      }

    }

  )

  #Data preview

  output$tube_data_out <- renderTable({

    dfs <- tube_data_frames()

    if (!length(dfs)) return(NULL)

    head(dfs[[1]], 10)
  }, striped = TRUE, bordered = TRUE, hover = TRUE, rownames = FALSE)

  #Miami_Plot state

  log_text <- reactiveVal(""); append_log <- function(...) isolate(log_text(paste0(log_text(), format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse=""), "\n")))

  current_grob <- reactiveVal(NULL); preview_file <- reactiveVal(NULL)

  #Single_Plot state

  single_sources <- reactiveVal(list())
  sp_girafe_df     <- reactiveVal(NULL)
  sp_selected_row  <- reactiveVal(NULL)

  # Track the active Single tab key
  single_selected <- reactiveVal(NULL)

  observeEvent(input$single_output_tabs, {

    if (!is.null(input$single_output_tabs)) single_selected(input$single_output_tabs)

  }, ignoreInit = FALSE)

  log_text_single <- reactiveVal(""); append_log_single <- function(...) isolate(log_text_single(paste0(log_text_single(), format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse=""), "\n")))
  single_grobs <- reactiveVal(list()); single_previews <- reactiveVal(list()); single_labels <- reactiveVal(list())
  current_grob_single <- reactiveVal(NULL)

  observeEvent({

    list(input$save_width_single, input$save_height_single, input$save_dpi_single)
  }, {

    grobs <- single_grobs()

    if (!length(grobs)) return(invisible(NULL))

    new_previews <- list()

    for (k in names(grobs)) {

      tmp <- tempfile(fileext = ".png")
      gg  <- ggplotify::as.ggplot(grobs[[k]])
      ggsave(
        filename = tmp, plot = gg,
        width = input$save_width_single,
        height = input$save_height_single,
        dpi = input$save_dpi_single,
        units = "in", limitsize = FALSE
      )

      new_previews[[k]] <- tmp
    }

    single_previews(new_previews)
  }, ignoreInit = TRUE)


  #Regional state

  regional_results  <- reactiveVal(list())
  regional_previews <- reactiveVal(list())

  reg_last_w   <- reactiveVal(30)
  reg_last_dpi <- reactiveVal(100)

  observeEvent(input$regional_save_width, {

    w <- suppressWarnings(as.numeric(input$regional_save_width))

    if (is.finite(w) && w > 0) reg_last_w(w)

  }, ignoreInit = TRUE)

  observeEvent(input$regional_save_dpi, {

    d <- suppressWarnings(as.numeric(input$regional_save_dpi))

    if (is.finite(d) && d > 0) reg_last_dpi(d)

  }, ignoreInit = TRUE)


  #Annotate state

  annotate_results <- reactiveVal(list())

  #METASOFT state

  meta_results <- reactiveVal(NULL)

  #Forest state

  forest_log <- reactiveVal("")

  forest_append_log <- function(...) isolate(forest_log(paste0(

    forest_log(), format(Sys.time(), "%H:%M:%S"), " | ",
    paste0(..., collapse = ""), "\n"

  )))

  forest_fun <- reactiveVal(NULL)

  forest_defaults <- reactiveVal(NULL)

  # multiple outputs
  forest_grobs     <- reactiveVal(list())
  forest_previews  <- reactiveVal(list())
  forest_labels    <- reactiveVal(list())
  forest_heights   <- reactiveVal(numeric(0))

  forest_last_w   <- reactiveVal(15)
  forest_last_dpi <- reactiveVal(100)
  observeEvent(input$forest_save_width, {

    w <- suppressWarnings(as.numeric(input$forest_save_width))

    if (is.finite(w) && w > 0) forest_last_w(w)

  }, ignoreInit = TRUE)

  observeEvent(input$forest_save_dpi, {

    d <- suppressWarnings(as.numeric(input$forest_save_dpi))

    if (is.finite(d) && d > 0) forest_last_dpi(d)

  }, ignoreInit = TRUE)

  #Model_Munge state

  mm_df <- reactiveVal(NULL)

  mm_model <- reactiveVal(NULL)

  mm_munge_df <- reactiveVal(NULL)

  #Build UI for inherited Single_Plot args

  output$single_args_ui <- renderUI({

    single_args <- formals(.Single_Plot_original); keep <- setdiff(names(single_args), "Data")

    lapply(keep, function(arg) make_input(arg, single_args[[arg]]))

  })

  output$single_args_ui2 <- renderUI({

    single_args <- formals(.Single_Plot_original)

    lapply(names(single_args), function(arg) if (arg != "Data") make_input(arg, single_args[[arg]]))

  })


  output$hla_arg_inputs <- renderUI({
    fmls <- formals(.Segregate_HLA_Plot_original)
    drops <- c("Data","data","session","...")  # keep flexible
    keep <- setdiff(names(fmls), drops)

    lapply(keep, function(arg){
      make_input_like_meta(
        prefix  = "hla_",
        name    = arg,
        default = safe_eval_default_reg(fmls[[arg]], environment(.Segregate_HLA_Plot_original))
      )
    })
  })

  hla_build_call_preview <- reactive({
    fmls <- formals(.Segregate_HLA_Plot_original)
    drops <- c("Data","data","session","...")
    keep <- setdiff(names(fmls), drops)

    nm_vec <- if (!is.null(input$hla_files$name))
      sanitize(tools::file_path_sans_ext(input$hla_files$name))
    else character(0)

    data_str <- if (!length(nm_vec)) "Data" else paste0("c(", paste(nm_vec, collapse=", "), ")")

    arg_vals <- lapply(keep, function(nm) parse_input_val_meta(input[[paste0("hla_", nm)]]))
    names(arg_vals) <- keep
    shown <- Filter(Negate(is.null), arg_vals)

    args <- c(sprintf("Data = %s", data_str))
    if (length(shown)) {
      for (nm in names(shown)) {
        args <- c(args, sprintf("%s = %s", nm, format_arg_val_meta(shown[[nm]])))
      }
    }

    paste0("Segregate_HLA_Plot(\n  ", paste(args, collapse=",\n  "), "\n)")
  })

  output$hla_call_preview <- renderText(hla_build_call_preview())


  hla_results <- reactiveVal(list())  # named by dataset/tab

  # ---- Shared saver for downloads (png/pdf/jpg/svg) ----
  # ---- Shared saver for downloads (png/pdf/jpg/svg) ----
  hla_save_any <- function(p, file, ext, width_in, height_in, dpi = 100) {
    gg <- if (inherits(p, "ggplot")) p else ggplotify::as.ggplot(p)

    ext <- tolower(ext %||% "png")

    if (ext == "png") {
      ggplot2::ggsave(file, plot = gg, width = width_in, height = height_in,
                      units = "in", dpi = dpi, limitsize = FALSE)
      return(invisible(TRUE))
    }

    if (ext %in% c("jpg","jpeg")) {
      ggplot2::ggsave(file, plot = gg, width = width_in, height = height_in,
                      units = "in", dpi = dpi, device = "jpeg", limitsize = FALSE)
      return(invisible(TRUE))
    }

    if (ext == "pdf") {
      ggplot2::ggsave(file, plot = gg, width = width_in, height = height_in,
                      units = "in", device = cairo_pdf, limitsize = FALSE)
      return(invisible(TRUE))
    }

    if (ext == "svg") {
      if (!requireNamespace("svglite", quietly = TRUE)) {
        stop("SVG export requires the 'svglite' package. Install it with install.packages('svglite').")
      }
      ggplot2::ggsave(file, plot = gg, width = width_in, height = height_in,
                      units = "in", device = svglite::svglite, limitsize = FALSE)
      return(invisible(TRUE))
    }

    stop("Unknown export type: ", ext)
  }


  # ---- Segregate_HLA_Plot tab: RUN handler (full observeEvent block) ----
  # Assumes you already created:
  #   - hla_results <- reactiveVal(list())
  #   - UI inputs: hla_files, run_hla, hla_save_width, hla_save_dpi
  #   - Helpers: sanitize(), tube_save_plot_png(), `%||%`,
  #              hla_extract_slides(), hla_save_slide_png() (from earlier)

  observeEvent(input$run_hla, {

    req(input$hla_files)

    out_res <- list()

    assign_one_file <- function(file_row) {
      obj <- sanitize(tools::file_path_sans_ext(basename(file_row$name)))
      df  <- vroom::vroom(file_row$datapath, show_col_types = FALSE)
      assign(obj, df, envir = .GlobalEnv)
      obj
    }

    # Build arg list from UI inputs
    fun_run <- Segregate_HLA_Plot
    fmls_run <- formals(fun_run)
    drops <- c("Data","data","session","...")

    keep <- setdiff(names(fmls_run), drops)

    ui_args <- lapply(keep, function(nm) parse_input_val_meta(input[[paste0("hla_", nm)]]))
    names(ui_args) <- keep
    ui_args <- Filter(Negate(is.null), ui_args)

    # Optional: if function accepts session, pass it; otherwise omit
    accepts_session <- "session" %in% names(fmls_run)

    # Loop each uploaded dataset
    for (i in seq_len(nrow(input$hla_files))) {

      file_row <- input$hla_files[i, , drop = FALSE]
      obj_name <- assign_one_file(file_row)

      # Make unique dataset key
      ds_key <- obj_name
      if (ds_key %in% names(out_res)) ds_key <- make.unique(c(names(out_res), ds_key))[length(out_res) + 1L]

      args <- c(
        list(Data = get(obj_name, envir = .GlobalEnv)),
        ui_args
      )
      if (accepts_session) args$session <- session

      out <- tryCatch(
        do.call(fun_run, args),
        error = function(e) {
          showNotification(
            paste0("Segregate_HLA_Plot error for ", ds_key, ": ", e$message),
            type = "error", duration = 8
          )
          NULL
        }
      )
      if (is.null(out)) next

      slides <- hla_extract_slides(out)

      pngs   <- character(0)
      labels <- character(0)
      plots  <- list()

      if (length(slides)) {

        for (j in seq_along(slides)) {

          p   <- slides[[j]]$plot
          lab <- slides[[j]]$label %||% paste0("Plot ", j)

          # Preview images for carousel are always PNG
          png <- tempfile(fileext = ".png")

          # Use dynamic/recommended height if available, otherwise fall back to input box
          h_dyn <- hla_slide_height_in(p, fallback = input$hla_save_height %||% 10)
          w_in  <- input$hla_save_width %||% 30
          dpi   <- input$hla_save_dpi   %||% 100

          tryCatch(
            {
              tube_save_plot_png(p, png, width_in = w_in, height_in = h_dyn, dpi = dpi)
            },
            error = function(e) {
              showNotification(paste0("Failed to render preview slide ", j, " (", ds_key, "): ", e$message),
                               type = "error", duration = 8)
            }
          )

          pngs   <- c(pngs, png)
          labels <- c(labels, lab)
          plots[[j]] <- p
        }
      }

      out_res[[ds_key]] <- list(
        plots  = plots,
        pngs   = pngs,
        labels = labels
      )
    }

    if (!length(out_res)) {
      showNotification("No outputs produced. Check input files / arguments.", type = "warning", duration = 6)
      return(invisible(NULL))
    }

    hla_results(out_res)

    # default to first dataset tab
    isolate({
      updateTabsetPanel(session, "hla_output_tabs", selected = names(out_res)[1])
    })

  })




  # ---- Segregate_HLA_Plot: dynamic UI (tabs + slickR carousels) ----
  output$hla_dynamic_tabs <- renderUI({

    res <- hla_results()
    if (!length(res)) return(NULL)

    tabs <- lapply(names(res), function(k) {

      carId <- paste0("hla_carousel_", k)

      local({
        kk <- k
        output[[carId]] <- renderSlickR({

          imgs <- hla_results()[[kk]]$pngs
          if (!length(imgs)) return(NULL)

          w <- slickR(imgs, slideType = "img") +
            settings(dots = TRUE, infinite = TRUE, arrows = TRUE, adaptiveHeight = FALSE)

          # Track current slide index per dataset (optional; useful if you later add "Interactive View")
          htmlwidgets::onRender(
            w,
            sprintf(
              "function(el,x){
               var $s = $(el).find('.slick-slider');
               Shiny.setInputValue('%1$s', 1, {priority:'event'});
               $s.on('afterChange', function(e, slick, current){
                 Shiny.setInputValue('%1$s', current+1, {priority:'event'});
               });
             }",
              paste0('hla_current_slide_', kk)
            )
          )
        })
      })

      tabPanel(
        title = tools::toTitleCase(gsub("_", " ", k)),
        value = k,
        slickROutput(carId, width = "100%", height = "680px")
      )
    })

    do.call(tabsetPanel, c(list(id = "hla_output_tabs"), tabs))
  })


  observe({
    res <- hla_results()
    if (!length(res)) return(invisible(NULL))

    lapply(names(res), function(nm){
      local({
        id <- paste0("hla_slick_", nm)

        output[[id]] <- slickR::renderSlickR({
          pngs   <- res[[nm]]$pngs
          labels <- res[[nm]]$labels

          req(length(pngs))

          # Data URI avoids resourcePath headaches
          mk_src <- function(path) {
            if (requireNamespace("base64enc", quietly = TRUE)) {
              base64enc::dataURI(file = path, mime = "image/png")
            } else {
              # fallback: plain file path (may not load in browser without resourcePath)
              path
            }
          }

          slides <- lapply(seq_along(pngs), function(i){
            htmltools::tags$div(
              style = "text-align:center;",
              htmltools::tags$img(src = mk_src(pngs[i]), style = "max-width:100%; max-height:660px; width:auto; height:auto;"),
              htmltools::tags$div(style="margin-top:10px; font-weight:700; color:#111827;",
                                  labels[i] %||% "")
            )
          })

          slickR::slickR(
            obj = slides,
            slideId = paste0("hla_carousel_", nm)
          ) + slickR::settings(
            dots = TRUE,
            arrows = TRUE,
            adaptiveHeight = TRUE
          )
        })
      })
    })
  })



  output$hla_download_simple <- downloadHandler(
    filename = function() {
      ext <- tolower(input$hla_ext %||% "png")
      res <- hla_results()
      total_plots <- sum(vapply(res, function(x) length(x$plots), integer(1)))

      if (total_plots > 1) return(paste0("Segregate_HLA_Plot_", ext, ".zip"))
      paste0("Segregate_HLA_Plot.", ext)
    },
    content = function(file) {
      ext <- tolower(input$hla_ext %||% "png")
      res <- hla_results()
      req(length(res))

      plots_all <- list()
      for (nm in names(res)) {
        if (length(res[[nm]]$plots)) {
          for (j in seq_along(res[[nm]]$plots)) {
            plots_all[[paste0(nm, "_", j)]] <- res[[nm]]$plots[[j]]
          }
        }
      }
      req(length(plots_all))

      save_one <- function(p, path) {
        if (ext == "pdf") tube_save_plot_pdf(p, path, input$hla_save_width, input$hla_save_height)
        else tube_save_plot_png(p, path, input$hla_save_width, input$hla_save_height, input$hla_save_dpi, dpi = input$hla_save_dpi)
      }

      if (length(plots_all) == 1) {
        save_one(plots_all[[1]], file)
        return(invisible())
      }

      tmpdir <- tempfile("hla_dl_"); dir.create(tmpdir)
      files <- character(0)

      for (nm in names(plots_all)) {
        out <- file.path(tmpdir, paste0("HLA_", sanitize(nm), ".", ext))
        save_one(plots_all[[nm]], out)
        files <- c(files, out)
      }

      zip_it(files, file)
    }
  )



  #Miami_Plot - call preview & run

  miami_call_txt <- reactive({

    top_label    <- if (!is.null(input$top_file$name))    tools::file_path_sans_ext(basename(input$top_file$name))    else "Top_Data"

    bottom_label <- if (!is.null(input$bottom_file$name)) tools::file_path_sans_ext(basename(input$bottom_file$name)) else "Bottom_Data"

    args <- list(

      Top_Title    = if (shiny::isTruthy(input$top_title))    input$top_title

      else if (!is.null(input$top_file$name))    tools::file_path_sans_ext(input$top_file$name)    else NULL,

      Bottom_Title = if (shiny::isTruthy(input$bottom_title)) input$bottom_title

      else if (!is.null(input$bottom_file$name)) tools::file_path_sans_ext(input$bottom_file$name) else NULL,

      Verbose      = isTRUE(input$verbose)

    )

    if (isTRUE(input$use_adjust)) args$Adjust <- input$adjust

    single_args <- formals(.Single_Plot_original)
    keep <- setdiff(names(single_args), "Data")
    mode <- input$miami_mode

    if (is.null(mode)) mode <- "advanced"

    if (identical(mode, "simple")) {

      for (arg in keep) {

        v <- parse_input_val_meta(input[[paste0("miami_s_", arg)]])

        if (!is.null(v)) args[[arg]] <- v

      }

    }
    else {


      for (arg in keep) {

        tv <- parse_input_val_meta(input[[paste0("top_", arg)]])
        bv <- parse_input_val_meta(input[[paste0("bottom_", arg)]])

        if (!is.null(tv)) args[[paste0("Top_", arg)]] <- tv

        if (!is.null(bv)) args[[paste0("Bottom_", arg)]] <- bv

      }

    }

    shown <- Filter(Negate(is.null), args)

    parts <- vapply(names(shown), function(nm) paste0(nm, " = ", format_arg_val_meta(shown[[nm]])), character(1))

    paste0(
      "Miami_Plot(\n  Top_Data = ", top_label,
      ",\n  Bottom_Data = ", bottom_label,

      if (length(parts)) paste0(",\n  ", paste(parts, collapse = ",\n  ")) else "",

      "\n)"
    )

  })

  output$miami_call_preview <- renderText(miami_call_txt())

  observeEvent(input$run, {

    req(input$top_file, input$bottom_file)

    args <- list(
      Top_Data    = input$top_file$datapath,
      Bottom_Data = input$bottom_file$datapath,
      Verbose     = isTRUE(input$verbose)
    )

    # Titles (safe)
    if (shiny::isTruthy(input$top_title))    args$Top_Title    <- input$top_title

    else if (!is.null(input$top_file$name))  args$Top_Title    <- tools::file_path_sans_ext(input$top_file$name)

    if (shiny::isTruthy(input$bottom_title)) args$Bottom_Title <- input$bottom_title

    else if (!is.null(input$bottom_file$name)) args$Bottom_Title <- tools::file_path_sans_ext(input$bottom_file$name)


    if (isTRUE(input$use_adjust)) args$Adjust <- input$adjust

    single_args <- formals(.Single_Plot_original)
    keep <- setdiff(names(single_args), "Data")
    mode <- input$miami_mode

    if (is.null(mode)) mode <- "advanced"

    if (identical(mode, "simple")) {


      for (arg in keep) {

        v <- parse_input_val_meta(input[[paste0("miami_s_", arg)]])

        if (!is.null(v)) args[[arg]] <- v

      }

    } else {

      for (arg in keep) {

        tv <- parse_input_val_meta(input[[paste0("top_", arg)]])
        bv <- parse_input_val_meta(input[[paste0("bottom_", arg)]])

        if (!is.null(tv)) args[[paste0("Top_", arg)]] <- tv

        if (!is.null(bv)) args[[paste0("Bottom_", arg)]] <- bv

      }

    }

    append_log("Running Miami_Plot() with per-panel args ...")

    gp <- NULL
    out_txt <- capture.output({

      gp <- tryCatch(

        withCallingHandlers(

          do.call(Miami_Plot, args),
          message = function(m) { append_log(conditionMessage(m)); invokeRestart("muffleMessage") },
          warning = function(w) { append_log(paste("WARNING:", conditionMessage(w))); invokeRestart("muffleWarning") }

        ),
        error = function(e) {

          append_log("ERROR: ", conditionMessage(e))
          showNotification(paste("Error:", conditionMessage(e)), type = "error")
          NULL

        }

      )

    }, type = "output")

    if (length(out_txt)) append_log(paste(squash_progress_lines(out_txt), collapse = "\n"))

    if (!is.null(gp)) {

      current_grob(gp)
      tmp <- tempfile(fileext = ".png")
      gg  <- ggplotify::as.ggplot(gp)
      ggsave(
        filename = tmp, plot = gg,
        width = input$save_width, height = input$save_height,
        dpi = input$save_dpi, units = "in", limitsize = FALSE
      )
      preview_file(tmp)
      append_log("Done.")

    }
  })

  # ---- Segregate_HLA_Plot: Download handler ----
  output$hla_download_simple <- downloadHandler(

    filename = function() {
      ext <- tolower(input$hla_ext %||% "png")
      res <- hla_results()
      n_plots <- sum(vapply(res, function(x) length(x$plots), integer(1)))

      if (n_plots > 1) {
        paste0("Segregate_HLA_Plot_", ext, ".zip")
      } else {
        paste0("Segregate_HLA_Plot.", ext)
      }
    },

    content = function(file) {

      ext <- tolower(input$hla_ext %||% "png")
      res <- hla_results()
      req(length(res))

      n_plots <- sum(vapply(res, function(x) length(x$plots), integer(1)))
      w_in  <- input$hla_save_width  %||% 30
      dpi   <- input$hla_save_dpi    %||% 100

      # If single output, write directly
      if (n_plots <= 1) {
        ds <- names(res)[1]
        p  <- res[[ds]]$plots[[1]]
        req(!is.null(p))

        h_in <- input$hla_save_height %||% hla_slide_height_in(p, fallback = 10)
        hla_save_any(p, file, ext = ext, width_in = w_in, height_in = h_in, dpi = dpi)
        return(invisible())
      }

      # Otherwise zip
      tmpdir <- tempfile("hla_dl_"); dir.create(tmpdir, recursive = TRUE)
      files <- character(0)

      for (ds in names(res)) {
        for (j in seq_along(res[[ds]]$plots)) {
          p <- res[[ds]]$plots[[j]]
          if (is.null(p)) next

          lab <- res[[ds]]$labels[j] %||% paste0("Plot_", j)
          fn  <- file.path(
            tmpdir,
            paste0(sanitize(ds), "__", sprintf("%02d", j), "__", sanitize(lab), ".", ext)
          )

          h_in <- input$hla_save_height %||% hla_slide_height_in(p, fallback = 10)
          hla_save_any(p, fn, ext = ext, width_in = w_in, height_in = h_in, dpi = dpi)

          files <- c(files, fn)
        }
      }

      zip_it(files, file)
    }
  )

  output$miami_plot <- renderImage({ req(preview_file()); list(src = preview_file(), contentType = "image/png", width = "100%") }, deleteFile = FALSE)
  output$log <- renderText({ log_text() })

  output$download_png <- downloadHandler(
    filename = function() paste0("MiamiPlot_", sanitize(file_path_sans_ext(input$top_file$name)), "_vs_", sanitize(file_path_sans_ext(input$bottom_file$name)), ".png"),
    content  = function(file) { req(current_grob()); gg <- ggplotify::as.ggplot(current_grob()); ggsave(file, gg, width=input$save_width, height=input$save_height, dpi=input$save_dpi) },
    contentType = "image/png"
  )

  output$download_pdf <- downloadHandler(
    filename = function() paste0("MiamiPlot_", sanitize(file_path_sans_ext(input$top_file$name)), "_vs_", sanitize(file_path_sans_ext(input$bottom_file$name)), ".pdf"),
    content  = function(file) { req(current_grob()); gg <- ggplotify::as.ggplot(current_grob()); ggsave(file, gg, width=input$save_width, height=input$save_height, dpi=input$save_dpi, device=cairo_pdf) },
    contentType = "application/pdf"
  )

  output$regional_download_simple <- downloadHandler(
    filename = function() {

      ext <- tolower(input$regional_ext %||% "png")
      res <- regional_results()
      total <- if (!length(res)) 0 else sum(vapply(res, function(x) length(x$plots), integer(1)))

      if (total <= 1) {

        curTab <- input$regional_output_tabs

        base   <- if (!is.null(curTab) && curTab %in% names(regional_results())) {

          sanitize(regional_results()[[curTab]]$label)

        } else "RegionalPlot"

        paste0("RegionalPlot_", base, ".", ext)

      } else {

        paste0("Regional_ALL_", ext, ".zip")

      }

    },

    content = function(file) {

      w   <- reg_last_w()
      dpi <- reg_last_dpi()

      if (!is.finite(w) || w <= 0)  stop("Invalid saved width.")

      if (!is.finite(dpi) || dpi <= 0) stop("Invalid saved DPI.")

      res <- regional_results()
      ext <- tolower(input$regional_ext %||% "png")
      req(length(res) > 0)

      save_one <- function(p, path, w, h, dpi) {

        gg <- if (inherits(p, c("gg","ggplot"))) p else ggplotify::as.ggplot(p)

        if (ext == "png") {

          ggsave(path, gg, width = w, height = h, dpi = dpi, units = "in", limitsize = FALSE)

        } else if (ext == "jpg") {

          ggsave(path, gg, device = "jpeg", width = w, height = h, dpi = dpi, units = "in", limitsize = FALSE)

        } else if (ext == "pdf") {

          ggsave(path, gg, device = cairo_pdf, width = w, height = h)

        } else if (ext == "svg") {

          if (!requireNamespace("svglite", quietly = TRUE)) stop("SVG export requires 'svglite'.")

          ggsave(path, gg, device = svglite::svglite, width = w, height = h)

        } else {

          ggsave(path, gg, width = w, height = h, dpi = dpi, units = "in", limitsize = FALSE)

        }

      }

      totals <- vapply(res, function(x) length(x$plots), integer(1))
      total  <- sum(totals)

      if (total == 1) {

        obj <- NULL
        lbl <- "RegionalPlot"

        if (length(res) == 1) {

          lbl <- res[[1]]$label
          obj <- res[[1]]$plots[[1]]

        }

        if (is.null(obj)) {

          curTab <- input$regional_output_tabs

          if (!is.null(curTab) && curTab %in% names(res)) {

            lbl <- res[[curTab]]$label
            idx_input <- paste0("regional_current_slide_", curTab)
            idx <- isolate(input[[idx_input]])

            if (is.null(idx) || !is.finite(idx)) idx <- 1

            obj <- res[[curTab]]$plots[[idx]]

          }

        }

        req(!is.null(obj))

        h <- attr(obj, "dynamic_height"); if (is.null(h) || !is.finite(h)) h <- 7

        save_one(obj, file, w, h, dpi)

        return(invisible())

      }

      tmpdir <- tempfile("regional_all_"); dir.create(tmpdir)

      files <- character()

      for (k in names(res)) {

        lbl <- sanitize(res[[k]]$label)
        plots <- res[[k]]$plots
        nd <- nchar(as.character(length(plots)))

        for (i in seq_along(plots)) {

          p <- plots[[i]]
          h <- attr(p, "dynamic_height"); if (is.null(h) || !is.finite(h)) h <- 7

          fn <- file.path(tmpdir, sprintf("Regional_%s_%0*d.%s", lbl, nd, i, ext))
          save_one(p, fn, w, h, dpi)

          files <- c(files, fn)

        }

      }

      zip_it(files, file)

    }

  )

  output$annotate_download_simple <- downloadHandler(

    filename = function() {

      ext <- tolower(input$annotate_ext %||% "csv")
      res <- annotate_results()
      n   <- length(res)

      if (n <= 1) {


        key <- if (n == 1) {

          names(res)[1]

        } else {


          cur <- input$annotate_output_tabs

          if (!is.null(cur) && cur %in% names(res)) cur else names(res)[1]

        }

        nm <- if (!is.null(res[[key]]$label)) res[[key]]$label else "Annotated"

        paste0("Annotated_", sanitize(nm), ".", ext)

      } else {

        paste0("Annotated_ALL_", ext, ".zip")

      }
    },

    content = function(file) {

      res <- annotate_results()
      ext <- tolower(input$annotate_ext %||% "csv")
      n   <- length(res)
      req(n > 0)

      write_one <- function(df, path) {

        if (ext == "csv") {

          utils::write.csv(df, path, row.names = FALSE, na = "")

        } else if (ext == "tsv") {

          utils::write.table(df, path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")

        } else {

          utils::write.table(df, path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")

        }
      }

      if (n == 1) {

        key <- names(res)[1]
        write_one(res[[key]]$df, file)

      } else {

        tmpdir <- tempfile("annotated_all_"); dir.create(tmpdir)

        files <- vapply(names(res), function(k) {

          lbl <- if (!is.null(res[[k]]$label)) res[[k]]$label else k

          fn  <- file.path(tmpdir, paste0("Annotated_", sanitize(lbl), ".", ext))

          write_one(res[[k]]$df, fn)
          fn
        }, character(1))

        zip_it(files, file)

      }

    }

  )



  output$single_download_simple <- downloadHandler(

    filename = function() {

      ext <- tolower(input$single_ext %||% "png")
      n   <- length(single_grobs())

      if (n <= 1) {

        cur <- input$single_output_tabs

        nm  <- if (!is.null(cur) && cur %in% names(single_labels()))

          single_labels()[[cur]] else "SinglePlot"

        paste0("SinglePlot_", sanitize(nm), ".", ext)

      } else {

        paste0("Single_ALL_", ext, ".zip")

      }

    },

    content = function(file) {

      grobs    <- single_grobs()
      previews <- single_previews()
      labels   <- single_labels()
      ext      <- tolower(input$single_ext %||% "png")
      n        <- length(grobs)
      req(n > 0)

      if (ext == "png") {

        if (n == 1) {

          cur <- input$single_output_tabs

          key <- if (!is.null(cur) && cur %in% names(previews)) cur else names(previews)[1]

          file.copy(previews[[key]], file, overwrite = TRUE)

          return(invisible())

        } else {

          tmpdir <- tempfile("single_all_png_"); dir.create(tmpdir)

          files <- vapply(names(previews), function(k) {

            dst <- file.path(tmpdir, paste0("SinglePlot_", sanitize(labels[[k]]), ".png"))
            file.copy(previews[[k]], dst, overwrite = TRUE)

            dst

          }, character(1))

          zip_it(files, file)
          return(invisible())

        }

      }

      write_one <- function(g, path) {

        gg <- ggplotify::as.ggplot(g)

        if (ext == "jpg") {

          ggsave(path, gg, device = "jpeg",
                 width = input$save_width_single, height = input$save_height_single,
                 dpi = input$save_dpi_single, units = "in", limitsize = FALSE)

        } else if (ext == "pdf") {

          ggsave(path, gg, device = cairo_pdf,
                 width = input$save_width_single, height = input$save_height_single)

        } else if (ext == "svg") {

          if (!requireNamespace("svglite", quietly = TRUE)) stop("SVG export requires 'svglite'.")
          ggsave(path, gg, device = svglite::svglite,
                 width = input$save_width_single, height = input$save_height_single)

        } else {

          ggsave(path, gg,
                 width = input$save_width_single, height = input$save_height_single,
                 dpi = input$save_dpi_single, units = "in", limitsize = FALSE)

        }

      }

      if (n == 1) {

        cur <- input$single_output_tabs

        key <- if (!is.null(cur) && cur %in% names(grobs)) cur else names(grobs)[1]

        write_one(grobs[[key]], file)

      } else {

        tmpdir <- tempfile("single_all_"); dir.create(tmpdir)

        files <- vapply(names(grobs), function(k) {

          fn <- file.path(tmpdir, paste0("SinglePlot_", sanitize(labels[[k]]), ".", ext))
          write_one(grobs[[k]], fn); fn
        }, character(1))
        zip_it(files, file)

      }

    }

  )


  output$download_png_single <- downloadHandler(

    filename = function() {

      current <- input$single_output_tabs

      nm <- if (!is.null(current) && current %in% names(single_labels()))

        single_labels()[[current]] else "SinglePlot"

      paste0("SinglePlot_", sanitize(nm), ".png")
    },
    content  = function(file) {

      previews <- single_previews()
      cur <- input$single_output_tabs

      key <- if (!is.null(cur) && cur %in% names(previews))

        cur else names(previews)[1]

      file.copy(previews[[key]], file, overwrite = TRUE)

    },

    contentType = "image/png"

  )

  output$miami_download_png_simple <- downloadHandler(

    filename = function() {

      ext <- tolower(input$miami_ext %||% "png")

      top    <- if (!is.null(input$top_file$name))  sanitize(tools::file_path_sans_ext(input$top_file$name))  else "Top_Data"

      bottom <- if (!is.null(input$bottom_file$name)) sanitize(tools::file_path_sans_ext(input$bottom_file$name)) else "Bottom_Data"

      paste0("MiamiPlot_", top, "_vs_", bottom, ".", ext)
    },

    content = function(file) {

      req(current_grob())
      gg <- ggplotify::as.ggplot(current_grob())
      ext <- tolower(input$miami_ext %||% "png")

      if (ext == "png") {

        ggsave(file, gg,
               width  = input$save_width, height = input$save_height,
               dpi    = input$save_dpi, units = "in", limitsize = FALSE)

        return(invisible())

      }

      if (ext == "jpg") {

        ggsave(file, gg, device = "jpeg",
               width  = input$save_width, height = input$save_height,
               dpi    = input$save_dpi, units = "in", limitsize = FALSE)

        return(invisible())

      }

      if (ext == "pdf") {

        ggsave(file, gg, device = cairo_pdf,
               width  = input$save_width, height = input$save_height)

        return(invisible())

      }

      if (ext == "svg") {

        if (!requireNamespace("svglite", quietly = TRUE)) {

          stop("SVG export requires the 'svglite' package. Install it or choose another format.")

        }

        ggsave(file, gg, device = svglite::svglite,
               width  = input$save_width, height = input$save_height)

        return(invisible())

      }

      ggsave(file, gg,
             width  = input$save_width, height = input$save_height,
             dpi    = input$save_dpi, units = "in", limitsize = FALSE)
    }
  )


  #Single_Plot

  single_call_txt <- reactive({

    lbls <- if (!is.null(input$single_file$name)) file_path_sans_ext(basename(input$single_file$name)) else character(0)

    file_label <- if (length(lbls) == 1) lbls else if (length(lbls) > 1) "[multiple files]" else "Data"
    args <- list()
    single_args <- formals(.Single_Plot_original)

    for (arg in setdiff(names(single_args), "Data")) {

      val <- input[[arg]]

      if (!is.null(val) && !(is.character(val) && val == "")) args[[arg]] <- tryCatch(eval(parse(text = val)), error = function(e) val)

    }

    if (!("Title" %in% names(args)) || is.null(args$Title) || (is.character(args$Title) && !nzchar(args$Title))) args$Title <- file_label

    shown <- Filter(Negate(is.null), args)
    parts <- vapply(names(shown), function(nm) paste0(nm, " = ", format_arg_val_meta(shown[[nm]])), "")

    paste0("Single_Plot(\n  Data = ", file_label, if (length(parts)) paste0(",\n  ", paste(parts, collapse = ",\n  ")) else "", "\n)")

  })

  output$single_call_preview <- renderText(single_call_txt())

  observeEvent(input$run_single, {

    req(input$single_file)
    grobs <- list(); previews <- list(); labels <- list()
    src_map <- list()

    single_args_formals <- formals(.Single_Plot_original)

    for (i in seq_len(nrow(input$single_file))) {

      datapath <- input$single_file$datapath[i]
      dispname <- file_path_sans_ext(input$single_file$name[i])
      key <- sanitize(dispname)

      args <- list(Data = datapath)

      for (arg in setdiff(names(single_args_formals), "Data")) {

        val <- input[[arg]]

        if (!is.null(val) && !(is.character(val) && val == "")) args[[arg]] <- tryCatch(eval(parse(text = val)), error = function(e) val)

      }

      if (!("Title" %in% names(args)) || is.null(args$Title) || (is.character(args$Title) && !nzchar(args$Title))) args$Title <- dispname

      append_log_single("Running Single_Plot for '", dispname, "' ...")
      gp <- tryCatch({ do.call(Single_Plot, args) }, error = function(e) { append_log_single("ERROR (", dispname, "): ", conditionMessage(e)); NULL })

      if (!is.null(gp)) {

        grobs[[key]] <- gp; labels[[key]] <- dispname
        tmp <- tempfile(fileext = ".png")
        gg <- ggplotify::as.ggplot(gp)
        ggsave(tmp, gg, width=input$save_width_single, height=input$save_height_single, dpi=input$save_dpi_single)
        previews[[key]] <- tmp
        append_log_single("Done: ", dispname)

        src_map[[key]] <- datapath

      }

    }

    single_grobs(grobs); single_labels(labels); single_previews(previews)
    single_sources(src_map)

    if (length(grobs)) current_grob_single(grobs[[1]]) else current_grob_single(NULL)

  })


  output$single_dynamic_tabs <- renderUI({

    grobs    <- single_grobs()
    previews <- single_previews()
    labels   <- single_labels()

    current <- isolate(single_selected())

    if (is.null(current) || !(current %in% names(grobs))) {

      current <- if (length(grobs)) names(grobs)[1] else NULL

    }

    tabs <- lapply(names(grobs), function(k) {

      tabPanel(
        title = labels[[k]], value = k,
        div(imageOutput(paste0("single_plot_", k)))
      )

    })

    for (k in names(previews)) local({

      kk <- k
      output[[paste0("single_plot_", kk)]] <- renderImage({

        list(src = single_previews()[[kk]], contentType = "image/png", width = "100%")

      }, deleteFile = FALSE)

    })

    do.call(tabsetPanel, c(list(id = "single_output_tabs", selected = current), tabs))

  })

  output$log_single <- renderText({ log_text_single() })

  output$download_png_single <- downloadHandler(

    filename = function() {

      current <- input$single_output_tabs

      nm <- if (!is.null(current) && current %in% names(single_labels())) single_labels()[[current]] else "SinglePlot"

      paste0("SinglePlot_", sanitize(nm), ".png")

    },

    content  = function(file) {

      current <- input$single_output_tabs; req(!is.null(current), current %in% names(single_grobs()))
      gg <- ggplotify::as.ggplot(single_grobs()[[current]])
      ggsave(file, gg, width=input$save_width_single, height=input$save_height_single, dpi=input$save_dpi_single)

    },
    contentType = "image/png"
  )
  output$download_pdf_single <- downloadHandler(

    filename = function() {

      current <- input$single_output_tabs
      nm <- if (!is.null(current) && current %in% names(single_labels())) single_labels()[[current]] else "SinglePlot"
      paste0("SinglePlot_", sanitize(nm), ".pdf")

    },

    content  = function(file) {

      current <- input$single_output_tabs; req(!is.null(current), current %in% names(single_grobs()))
      gg <- ggplotify::as.ggplot(single_grobs()[[current]])
      ggsave(file, gg, width=input$save_width_single, height=input$save_height_single, dpi=input$save_dpi_single, device=cairo_pdf)

    },
    contentType = "application/pdf"

  )

  observeEvent(input$miami_interactive_view, {
    req(input$top_file, input$bottom_file)

    mi_args <- list(
      Top_Data    = input$top_file$datapath,
      Bottom_Data = input$bottom_file$datapath,
      Verbose     = isTRUE(input$verbose)
    )

    if (shiny::isTruthy(input$top_title))      mi_args$Top_Title    <- input$top_title

    else if (!is.null(input$top_file$name))    mi_args$Top_Title    <- tools::file_path_sans_ext(input$top_file$name)

    if (shiny::isTruthy(input$bottom_title))   mi_args$Bottom_Title <- input$bottom_title

    else if (!is.null(input$bottom_file$name)) mi_args$Bottom_Title <- tools::file_path_sans_ext(input$bottom_file$name)

    if (isTRUE(input$use_adjust)) mi_args$Adjust <- input$adjust


    sp_fmls <- formals(.Single_Plot_original)
    keep <- setdiff(names(sp_fmls), "Data")
    mode <- input$miami_mode %||% "advanced"

    if (identical(mode, "simple")) {

      for (arg in keep) {

        v <- parse_input_val_meta(input[[paste0("miami_s_", arg)]])

        if (!is.null(v)) mi_args[[arg]] <- v

      }

    } else {

      for (arg in keep) {

        tv <- parse_input_val_meta(input[[paste0("top_", arg)]])
        bv <- parse_input_val_meta(input[[paste0("bottom_", arg)]])

        if (!is.null(tv)) mi_args[[paste0("Top_", arg)]]    <- tv

        if (!is.null(bv)) mi_args[[paste0("Bottom_", arg)]] <- bv

      }

    }

    #Miami_Plot

    miami <- tryCatch(
      do.call(Miami_Plot, mi_args),
      error = function(e) {

        showNotification(paste("Interactive Miami_Plot error:", e$message),
                         type = "error", duration = 8)
        NULL

      }
    )
    req(!is.null(miami), is.list(miami), length(miami) >= 2)

    #Build hover overlays

    safe_hover <- function(p, prefix) {

      tryCatch(

        hover_overlay_from_plot(p, prefix = prefix),

        error = function(e) {

          showNotification(paste("Interactive hover fallback:", e$message),
                           type = "warning", duration = 6)
          reg_ensure_hover(p, size_pts = 50, tag = gsub(":", "", prefix))

        }

      )
    }

    top_hover    <- safe_hover(miami[[1]], "top:")
    bottom_hover <- safe_hover(miami[[2]], "bot:")

    #Overlay maps

    ov_top_map <- attr(top_hover, "rows_map")
    ov_bot_map <- attr(bottom_hover, "rows_map")

    top_raw <- tryCatch(

      vroom::vroom(input$top_file$datapath, show_col_types = FALSE, progress = FALSE),
      error = function(e) NULL

    )

    bot_raw <- tryCatch(

      vroom::vroom(input$bottom_file$datapath, show_col_types = FALSE, progress = FALSE),
      error = function(e) NULL

    )

    top_cols <- if (!is.null(top_raw)) names(top_raw) else attr(top_hover, "cols")

    bot_cols <- if (!is.null(bot_raw)) names(bot_raw) else attr(bottom_hover, "cols")

    pick_id_col <- function(df) {

      intersect(c("ID","SNP","RSID","RS","RS_NUMBER"), names(df))[1]

    }

    as_rowlist <- function(row_df) {


      vapply(row_df, function(x) {
        x <- as.character(x); x[is.na(x)] <- ""; paste(x, collapse = ", ")
      }, character(1), USE.NAMES = TRUE)

    }
    build_map_from_original <- function(df, ids_raw, prefix) {

      if (is.null(df) || !length(ids_raw)) return(list())
      idcol <- pick_id_col(df)
      if (is.na(idcol) || !nzchar(idcol)) return(list())
      idx <- match(ids_raw, as.character(df[[idcol]]))
      ok  <- which(!is.na(idx))
      if (!length(ok)) return(list())
      setNames(lapply(ok, function(i) as.list(as_rowlist(df[idx[i], , drop = FALSE]))),
               paste0(prefix, ids_raw[ok]))

    }

    top_ids_raw <- if (length(ov_top_map)) sub("^top:", "", names(ov_top_map)) else character(0)

    bot_ids_raw <- if (length(ov_bot_map)) sub("^bot:", "", names(ov_bot_map)) else character(0)


    map_top_orig <- if (!is.null(top_raw)) build_map_from_original(top_raw, top_ids_raw, "top:") else list()

    map_bot_orig <- if (!is.null(bot_raw)) build_map_from_original(bot_raw, bot_ids_raw, "bot:") else list()

    rows_map <- c(ov_top_map, ov_bot_map)
    rows_map[names(map_top_orig)] <- map_top_orig
    rows_map[names(map_bot_orig)] <- map_bot_orig

    col_order <- list(top = top_cols, bot = bot_cols)

    combined_try <- try({

      (top_hover / bottom_hover) + patchwork::plot_layout(heights = c(1, 1))
    }, silent = TRUE)

    g <- ggiraph::girafe(
      ggobj      = if (!inherits(combined_try, "try-error")) combined_try else list(top_hover, bottom_hover),
      width_svg  = 30,
      height_svg = 16,
      options = list(
        ggiraph::opts_sizing(width = 1, rescale = TRUE),
        ggiraph::opts_hover(css = "opacity:1; stroke:black; stroke-width:3px; r:6;"),
        ggiraph::opts_selection(type = "single"),
        ggiraph::opts_tooltip(css = "background:white;color:black;padding:8px;border-radius:6px;font-size:14px;box-shadow:0 2px 6px rgba(0,0,0,0.2);white-space:pre-line;"),
        ggiraph::opts_toolbar(saveaspng = TRUE, position = "topright"),
        ggiraph::opts_zoom(min = 0.5, max = 2)
      )
    )

    g <- htmlwidgets::prependContent(

      g,
      htmltools::tags$style(htmltools::HTML("
      html, body { margin:0; background:#fff; }
      .girafe_container_std { margin:0 !important; padding:0 !important; text-align:left !important; }
      .girafe_container_std svg { display:block; overflow:visible !important; }
      #snp-modal { position: fixed; inset:0; display:none; align-items:center; justify-content:center; background:rgba(0,0,0,.4); z-index:99999; }
      #snp-modal .content { position:relative; background:#fff; width:min(720px,92vw); border-radius:16px; padding:22px 24px 18px; box-shadow:0 18px 40px rgba(0,0,0,.25); border:1px solid #e5e7eb; }
      #snp-title { margin:0 0 12px; font-weight:800; font-size:20px; color:#2c3e50; }
      #snp-close { position:absolute; top:10px; right:10px; width:28px; height:28px; border-radius:9999px; border:1px solid #e5e7eb; background:#f3f4f6; color:#374151; font-weight:700; cursor:pointer; }
      #snp-close:hover { background:#e5e7eb; }
      #snp-row { max-height:50vh; overflow:auto; margin:6px 0 12px; }
      table.snp-table { width:100%; border-collapse:collapse; table-layout:fixed; }
      table.snp-table th, table.snp-table td { border:1px solid #e5e7eb; padding:6px 8px; vertical-align:top; word-break:break-word; }
      table.snp-table th { background:#f9fafb; width:30%; text-align:left; color:#374151; }
      table.snp-table td { background:#ffffff; color:#111827; }
      .snp-actions { display:flex; gap:10px; margin-top:8px; flex-wrap:wrap; }
      .btn-google, .btn-dbsnp { display:inline-block; text-decoration:none; padding:10px 14px; font-weight:700; color:#fff; border-radius:8px; }
      .btn-google { background:#4285F4; } .btn-google:hover { background:#3367D6; }
      .btn-dbsnp  { background:#28a745; } .btn-dbsnp:hover  { background:#218838; }
    "))
    )

    g <- htmlwidgets::onRender(
      g,
      "function(el, x, data){
      var modal = document.getElementById('snp-modal');
      if(!modal){
        modal = document.createElement('div'); modal.id = 'snp-modal';
        modal.innerHTML = `
          <div class='content'>
            <button id='snp-close' aria-label='Close'>×</button>
            <h2 id='snp-title'></h2>
            <div id='snp-row'></div>
            <div class='snp-actions'>
              <a id='snp-google' class='btn-google' target='_blank' rel='noopener'>🔍 Google</a>
              <a id='snp-dbsnp'  class='btn-dbsnp'  target='_blank' rel='noopener'>🧬 dbSNP</a>
            </div>
          </div>`;
        document.body.appendChild(modal);
        modal.querySelector('#snp-close').addEventListener('click', function(){ modal.style.display='none'; });
        modal.addEventListener('click', function(e){ if(e.target===modal) modal.style.display='none'; });
        document.addEventListener('keydown', function(e){ if(e.key==='Escape') modal.style.display='none'; });
      }
      var titleEl = modal.querySelector('#snp-title');
      var rowEl   = modal.querySelector('#snp-row');
      var aGo     = modal.querySelector('#snp-google');
      var aDb     = modal.querySelector('#snp-dbsnp');

      var MAP  = (data && data.rows_map)  || {};
      var CORD = (data && data.col_order) || {};

      function renderTable(container, rowMap, order){
        container.innerHTML = '';
        if(!rowMap) return;
        // keep only keys that exist in the row and keep EXACT order from 'order'
        var keys = Array.isArray(order) ? order.filter(function(k){ return Object.prototype.hasOwnProperty.call(rowMap, k); })
                                        : Object.keys(rowMap || {});
        var tbl = document.createElement('table'); tbl.className='snp-table';
        var tb  = document.createElement('tbody');
        keys.forEach(function(k){
          var v = rowMap[k]; if (v == null) return;
          var vs = String(v); if (!vs) return;
          var tr = document.createElement('tr');
          var th = document.createElement('th'); th.textContent = k;
          var td = document.createElement('td'); td.textContent = vs;
          tr.appendChild(th); tr.appendChild(td); tb.appendChild(tr);
        });
        tbl.appendChild(tb); container.appendChild(tbl);
      }

      function findRS(s){
        var m = String(s||'').match(/\\brs\\d+\\b/i);
        return m ? m[0] : '';
      }

      function clickOpen(id, tip){
        var row   = MAP[id] || {};
        var which = id.indexOf('top:')===0 ? 'top' : (id.indexOf('bot:')===0 ? 'bot' : 'top');
        var order = (CORD && CORD[which]) || Object.keys(row);
        titleEl.textContent = id ? ('ID: ' + id) : 'Variant';
        renderTable(rowEl, row, order);

        var rs = findRS(id) || findRS(tip);
        var q  = encodeURIComponent(rs || id);
        aGo.href = 'https://www.google.com/search?q=' + q;
        aDb.href = rs ? ('https://www.ncbi.nlm.nih.gov/snp/' + rs)
                      : ('https://www.ncbi.nlm.nih.gov/snp/?term=' + q);

        modal.style.display = 'flex';
      }

      el.addEventListener('click', function(e){
        var node = e.target.closest('[data-id]'); if(!node) return;
        var id   = node.getAttribute('data-id') || '';
        var tip  = node.getAttribute('data-tooltip') || node.getAttribute('title') ||
                   (node.querySelector && (node.querySelector('title,desc') && node.querySelector('title,desc').textContent)) || '';
        clickOpen(id, tip);
      });
    }",
      data = list(rows_map = rows_map, col_order = col_order)
    )

    # Save temp widget and open in a new tab
    tmphtml <- tempfile(pattern = "miami_girafe_", fileext = ".html")
    htmlwidgets::saveWidget(g, tmphtml, selfcontained = TRUE)
    prefix <- paste0("rview_", as.integer(runif(1, 1e6, 1e9)))
    shiny::addResourcePath(prefix, dirname(tmphtml))
    rel_url <- paste(prefix, basename(tmphtml), sep = "/")
    session$sendCustomMessage("open-window", rel_url)

  })

  observeEvent(input$single_interactive_view, {

    grobs <- single_grobs()

    if (!length(grobs)) {

      showNotification("Run Single_Plot first.", type = "warning", duration = 6)

      return(invisible(NULL))

    }

    cur <- input$single_output_tabs

    if (is.null(cur) || !(cur %in% names(grobs))) cur <- names(grobs)[1]

    p <- grobs[[cur]]

    if (!inherits(p, "ggplot")) p <- ggplotify::as.ggplot(p)

    local_overlay_sig <- function(p, prefix = "sp:") {

      df <- as.data.frame(p$data)

      if (!nrow(df)) stop("Single plot has no data.")

      pick <- function(cands) { for (nm in cands) if (nm %in% names(df)) return(nm); NULL }

      xcol <- pick(c("new_pos","GENPOS","POS","BP","position"))

      if (is.null(xcol)) stop("Could not find genomic position column (POS/BP/GENPOS).")


      y_alt <- pick(c("LOG10P","MLP","log10P","Log10P"))
      y_vals <- NULL

      if ("P" %in% names(df)) {

        suppressWarnings({ y_vals <- -log10(as.numeric(df$P)) })

      }

      if (is.null(y_vals) && !is.null(y_alt)) {

        suppressWarnings({ y_vals <- as.numeric(df[[y_alt]]) })

      }

      if (is.null(y_vals)) stop("No P or log10P-like column found for interactive overlay.")

      #significance mask (P < 5e-8 OR log10P > 7.30103)

      sig_thr <- -log10(5e-8)  # ~7.30103
      pass <- rep(TRUE, nrow(df))
      if ("P" %in% names(df) && is.numeric(df$P)) {

        suppressWarnings({ pass <- is.finite(df$P) & (as.numeric(df$P) < 5e-8) })

      } else if (!is.null(y_alt) && is.numeric(df[[y_alt]])) {

        pass <- is.finite(df[[y_alt]]) & (as.numeric(df[[y_alt]]) > sig_thr)

      } else {

        pass <- is.finite(y_vals) & (y_vals > sig_thr)

      }

      pass <- pass & is.finite(df[[xcol]]) & is.finite(y_vals)

      if (!any(pass)) {

        attr(p, "rows_map") <- list()
        attr(p, "cols") <- names(df)

        return(p)

      }

      df2 <- df[pass, , drop = FALSE]
      y2  <- y_vals[pass]

      idcol  <- pick(c("ID","SNP","RSID","RS","RS_NUMBER"))
      chrcol <- pick(c("CHR","CHROM","Chromosome","chrom"))
      poscol <- pick(c("GENPOS","POS","BP","position"))
      refcol <- pick(c("REF","REF_ALLELE","ALLELE0","A1","EA","EFFECT_ALLELE"))
      altcol <- pick(c("ALT","ALT_ALLELE","ALLELE1","A2","NEA","NON_EFFECT_ALLELE"))

      if ("Hover_Info" %in% names(df2)) {

        tooltip_vec <- as.character(df2$Hover_Info); tooltip_vec[is.na(tooltip_vec)] <- ""

      } else {

        mk <- function(lbl, col) if (!is.null(col)) paste0(lbl, df2[[col]]) else NULL

        tt_parts <- list(
          mk("SNP: ", idcol),
          mk("CHR: ", chrcol),
          mk("POS: ", poscol),

          if ("P" %in% names(df2)) paste0("P: ", signif(as.numeric(df2$P), 4)) else NULL,

          if (!is.null(y_alt))     paste0("log10P: ", signif(as.numeric(df2[[y_alt]]), 4)) else NULL,

          mk("REF: ", refcol),
          mk("ALT: ", altcol)
        )

        tooltip_vec <- vapply(seq_len(nrow(df2)), function(i) {

          paste0(Filter(Negate(is.null), lapply(tt_parts, `[`, i)), collapse = "\n")

        }, character(1))
      }

      src_idx <- which(pass)
      hover_df <- data.frame(
        .x  = df2[[xcol]],
        .y  = y2,
        .id_raw  = if (!is.null(idcol)) as.character(df2[[idcol]]) else as.character(src_idx),
        .tt      = tooltip_vec,
        .src_idx = src_idx,
        stringsAsFactors = FALSE
      )
      hover_df$.id <- paste0(prefix, hover_df$.id_raw)

      rows_map <- setNames(
        lapply(hover_df$.src_idx, function(i) {
          vals <- vapply(df[i, , drop = FALSE], function(v) {
            v <- as.character(v); v[is.na(v)] <- ""; paste(v, collapse = ", ")
          }, character(1))
          as.list(vals)
        }),
        hover_df$.id
      )

      p2 <- p + ggiraph::geom_point_interactive(
        data = hover_df,
        ggplot2::aes(x = .x, y = .y, tooltip = .tt, data_id = .id),
        inherit.aes = FALSE,
        size = 8, stroke = 0, alpha = 0  # big invisible capture
      )
      attr(p2, "rows_map") <- rows_map
      attr(p2, "cols")     <- names(df)
      p2
    }

    p_i <- try(local_overlay_sig(p, prefix = "sp:"), silent = TRUE)

    if (inherits(p_i, "try-error")) {

      showNotification(paste("Interactive overlay error:", as.character(p_i)), type = "error", duration = 8)

      return(invisible(NULL))
    }

    rows_map <- attr(p_i, "rows_map")
    if (is.null(rows_map) || !length(rows_map)) {

      showNotification("No points pass P < 5e-8 for this plot.", type = "message", duration = 6)

    }

    #Build widget

    h_in <- input$save_height_single %||% 10
    g <- ggiraph::girafe(
      ggobj      = p_i,
      width_svg  = 30,
      height_svg = h_in,
      options = list(
        ggiraph::opts_sizing(width = 1, rescale = TRUE),
        ggiraph::opts_hover(css = "opacity:1; stroke:black; stroke-width:3px; r:6;"),
        ggiraph::opts_selection(type = "single"),
        ggiraph::opts_tooltip(css = "background:white;color:black;padding:8px;border-radius:6px;font-size:14px;box-shadow:0 2px 6px rgba(0,0,0,0.2);white-space:pre-line;"),
        ggiraph::opts_toolbar(saveaspng = TRUE, position = "topright"),
        ggiraph::opts_zoom(min = 0.5, max = 2)
      )
    )

    g <- htmlwidgets::prependContent(

      g,

      htmltools::tags$style(htmltools::HTML("

      html, body { margin:0; background:#fff; }
      #sp-modal { position: fixed; inset:0; display:none; align-items:center; justify-content:center; background:rgba(0,0,0,.4); z-index:99999; }
      #sp-modal .content { position:relative; background:#fff; width:min(720px,92vw); border-radius:16px; padding:22px 24px 18px; box-shadow:0 18px 40px rgba(0,0,0,.25); border:1px solid #e5e7eb; }
      #sp-title { margin:0 0 12px; font-weight:800; font-size:20px; color:#2c3e50; }
      #sp-close { position:absolute; top:10px; right:10px; width:28px; height:28px; border-radius:9999px; border:1px solid #e5e7eb; background:#f3f4f6; color:#374151; font-weight:700; cursor:pointer; }
      #sp-close:hover { background:#e5e7eb; }
      #sp-row { max-height:50vh; overflow:auto; margin:6px 0 12px; }
      #sp-row table { width:100%; border-collapse:collapse; table-layout:fixed; }
      #sp-row th, #sp-row td { border:1px solid #e5e7eb; padding:6px 8px; vertical-align:top; word-break:break-word; }
      #sp-row th { background:#f9fafb; width:30%; text-align:left; color:#374151; }
      #sp-row td { background:#ffffff; color:#111827; }
      #sp-actions { display:flex; gap:10px; margin-top:8px; }
      #sp-actions a { display:inline-block; text-decoration:none; padding:10px 14px; font-weight:700; color:#fff; border-radius:8px; }
      #sp-google { background:#4285F4; } #sp-google:hover { background:#3367D6; }
      #sp-dbsnp  { background:#28a745; } #sp-dbsnp:hover  { background:#218838; }
    "))
    )

    data_payload <- list(rows_map = rows_map, cols = attr(p_i, "cols"))

    g <- htmlwidgets::onRender(

      g,
      "
    function(el, x, data){
      var modal = document.getElementById('sp-modal');
      if(!modal){
        modal = document.createElement('div'); modal.id='sp-modal';
        modal.innerHTML = `
          <div class='content'>
            <button id='sp-close' aria-label='Close'>×</button>
            <h2 id='sp-title'></h2>
            <div id='sp-row'></div>
            <div id='sp-actions'>
              <a id='sp-google' target='_blank' rel='noopener'>🔍 Google</a>
              <a id='sp-dbsnp'  target='_blank' rel='noopener'>🧬 dbSNP</a>
            </div>
          </div>`;
        document.body.appendChild(modal);
        modal.querySelector('#sp-close').addEventListener('click', function(){ modal.style.display='none'; });
        modal.addEventListener('click', function(e){ if(e.target===modal) modal.style.display='none'; });
        document.addEventListener('keydown', function(e){ if(e.key==='Escape') modal.style.display='none'; });
      }

      var titleEl = modal.querySelector('#sp-title');
      var rowEl   = modal.querySelector('#sp-row');
      var aGo     = modal.querySelector('#sp-google');
      var aDb     = modal.querySelector('#sp-dbsnp');

      var ROWS = (data && data.rows_map) || {};
      var COLS = (data && data.cols) || [];

    function renderRowTable(container, row, order){
  container.innerHTML = '';
  if(!row){ return; }
  var tbl=document.createElement('table');
  var tb=document.createElement('tbody');


  // start from provided column order (if any), else all keys

  var keys = (order && order.length) ? order.slice() : Object.keys(row).sort();

  // cut at CHROM (exclude CHROM and everything after it)

  var cut = keys.findIndex(function(k){ return String(k) === 'CHROM'; });
  if (cut >= 0) keys = keys.slice(0, cut);

  keys.forEach(function(k){
    var v=row[k]; if(v===null || v===undefined) return;
    var vs=String(v); if(!vs) return;
    var tr=document.createElement('tr');
    var th=document.createElement('th'); th.textContent=k;
    var td=document.createElement('td'); td.textContent=vs;
    tr.appendChild(th); tr.appendChild(td); tb.appendChild(tr);
  });

  tbl.appendChild(tb);
  container.appendChild(tbl);

}
      function findRS(s){ var m=String(s||'').match(/rs\\d+/i); return m?m[0]:''; }

      el.addEventListener('click', function(e){
        var node = e.target.closest('[data-id]'); if(!node) return;
        var id   = node.getAttribute('data-id') || '';
        var bare = id.replace(/^sp:/,'');
        var row  = ROWS[id];

        titleEl.textContent = bare ? ('SNP: ' + bare) : 'Variant';
        renderRowTable(rowEl, row, COLS);

        // Build links (prefer rsID; else CHR:POS)
        var rs = findRS(bare) || (row && (findRS(row.ID) || findRS(row.SNP) || findRS(row.RS)));
        var term = rs;
        if(!term && row){
          var chr = row.CHR || row.CHROM || row.chrom || row.Chromosome;
          var pos = row.POS || row.GENPOS || row.Position || row.BP || row.position;
          if(chr && pos) term = String(chr)+':'+String(pos);
        }

        var q = encodeURIComponent(term || bare || 'variant');
        aGo.href = 'https://www.google.com/search?q=' + q;
        aDb.href = rs ? ('https://www.ncbi.nlm.nih.gov/snp/' + rs)
                      : ('https://www.ncbi.nlm.nih.gov/snp/?term=' + q);

        modal.style.display = 'flex';
      });

    }

    ",

    data = data_payload

    )

    tmphtml <- tempfile(pattern = "single_girafe_", fileext = ".html")
    htmlwidgets::saveWidget(g, tmphtml, selfcontained = TRUE)
    prefix <- paste0("sview_", as.integer(runif(1, 1e6, 1e9)))
    shiny::addResourcePath(prefix, dirname(tmphtml))
    rel_url <- paste(prefix, basename(tmphtml), sep = "/")
    session$sendCustomMessage("open-window", rel_url)

  })

  #Regional_Plot


  get_target_fun <- function() {

    if (exists(".Regional_Plot_original", mode = "function")) get(".Regional_Plot_original")

    else if (exists("Regional_Plot", mode = "function")) get("Regional_Plot")

    else stop("Regional_Plot not found")

  }

  get_single_fun <- function() {

    if (exists(".Single_Plot_original", mode = "function")) get(".Single_Plot_original")

    else if (exists("Single_Plot", mode = "function")) get("Single_Plot")

    else stop("Single_Plot not found")

  }

  make_input_regional <- function(name, val) {

    inputId <- paste0("arg_", name)

    if (is.logical(val) && length(val) == 1) checkboxInput(inputId, label = name, value = val)

    else if (is.numeric(val) && length(val) == 1) numericInput(inputId, label = name, value = val)

    else if (is.character(val) && length(val) == 1) textInput(inputId, label = name, value = val)

    else textInput(inputId, label = paste0(name, ""), value = deparse(val, width.cutoff = 60))

  }

  parse_input_val_reg <- function(x) {

    if (is.null(x)) return(NULL)

    if (isTRUE(x) || identical(x, "TRUE")) return(TRUE)

    if (identical(x, FALSE) || identical(x, "FALSE")) return(FALSE)

    tryCatch(eval(parse(text = x)), error = function(e) x)

  }

  format_arg_val_reg <- function(v) {

    if (is.null(v)) return(NULL)

    if (is.logical(v) && length(v) == 1) return(ifelse(v, "TRUE", "FALSE"))

    if (is.numeric(v)) return(paste0("c(", paste(v, collapse = ","), ")"))

    if (is.character(v)) {

      if (length(v) == 1) return(paste0("\"", v, "\""))

      return(paste0("c(", paste(paste0("\"", v, "\""), collapse = ","), ")"))

    }

    deparse(v, width.cutoff = 60)

  }

  regional_defaults <- reactive({

    f_reg <- get_target_fun(); f_sin <- get_single_fun()
    reg_fmls <- formals(f_reg); sin_fmls <- formals(f_sin)
    drops <- c("Data","Top_Data","Bottom_Data","...",".dots")

    keep_reg <- setdiff(names(reg_fmls), drops); keep_sin <- setdiff(names(sin_fmls), drops)

    reg_def <- lapply(reg_fmls[keep_reg], function(x) safe_eval_default_reg(x, environment(f_reg))); names(reg_def) <- keep_reg

    sin_def <- lapply(sin_fmls[keep_sin], function(x) safe_eval_default_reg(x, environment(f_sin))); names(sin_def) <- keep_sin
    combined <- modifyList(sin_def, reg_def)

    if (!("Chromosomes" %in% names(combined)) && ("Chromosomes" %in% names(reg_fmls))) combined[["Chromosomes"]] <- safe_eval_default_reg(reg_fmls[["Chromosomes"]], environment(f_reg))

    if (!("Chromosome"  %in% names(combined)) && ("Chromosome"  %in% names(reg_fmls))) combined[["Chromosome"]]  <- safe_eval_default_reg(reg_fmls[["Chromosome"]],  environment(f_reg))

    combined

  })

  output$regional_arg_inputs <- renderUI({

    lapply(names(regional_defaults()), function(nm) if (!(nm %in% c("Chromosomes","Chromosome"))) make_input_regional(nm, regional_defaults()[[nm]]))

  })

  regional_arg_vals <- reactive({

    vals <- lapply(names(regional_defaults()), function(nm) {

      if (nm %in% c("Chromosomes","Chromosome")) c(input$regional_chrom_input)

      else {

        val <- input[[paste0("arg_", nm)]]

        if (is.null(val)) regional_defaults()[[nm]] else parse_input_val_reg(val)

      }

    })

    names(vals) <- names(regional_defaults()); vals

  })

  build_call_preview_reg <- function(arg_list, data_label = "user_data") {

    shown <- Filter(Negate(is.null), arg_list)

    parts <- vapply(names(shown), function(nm) paste0(nm, " = ", format_arg_val_reg(shown[[nm]])), "")

    paste0("Regional_Plot(\n  Data = ", data_label,

           if (length(parts)) paste0(",\n  ", paste(parts, collapse = ",\n  ")) else "",

           "\n)")

  }

  regional_call_txt <- reactive({

    labels_vec <- if (!is.null(input$regional_datafile$name)) file_path_sans_ext(basename(input$regional_datafile$name)) else character(0)

    file_label <- if (length(labels_vec) == 1) labels_vec else if (length(labels_vec) > 1) "[multiple files]" else "user_data"

    vals <- regional_arg_vals()

    if (!any(names(vals) %in% c("Chromosomes","Chromosome"))) vals[["Chromosomes"]] <- c(input$regional_chrom_input)

    build_call_preview_reg(vals, data_label = file_label)

  })

  output$regional_call_preview <- renderText(regional_call_txt())

  observeEvent(input$run_regional, {

    w_run   <- isolate(reg_last_w())
    dpi_run <- isolate(reg_last_dpi())

    if (!is.finite(w_run) || w_run <= 0) {

      showNotification("Save Width must be a positive number (inches).", type = "error", duration = 6)

      return(invisible(NULL))

    }

    if (!is.finite(dpi_run) || dpi_run <= 0) {

      showNotification("DPI must be a positive number.", type = "error", duration = 6)

      return(invisible(NULL))

    }

    req(input$regional_datafile)

    results <- list()
    previews <- list()

    vals <- regional_arg_vals()

    if (!any(names(vals) %in% c("Chromosomes","Chromosome"))) {

      vals[["Chromosomes"]] <- c(input$regional_chrom_input)

    }

    for (i in seq_len(nrow(input$regional_datafile))) {

      datapath <- input$regional_datafile$datapath[i]
      dispname <- file_path_sans_ext(input$regional_datafile$name[i])
      key <- sanitize(dispname)

      df <- vroom::vroom(datapath, show_col_types = FALSE)
      assign(dispname, df, envir = .GlobalEnv)
      assign("user_data", get(dispname, envir = .GlobalEnv), envir = .GlobalEnv)

      call_txt <- build_call_preview_reg(vals, data_label = dispname)

      plots <- tryCatch(

        { eval(parse(text = call_txt), envir = .GlobalEnv) },

        error = function(e) {

          showNotification(paste("Error (", dispname, "):", e$message), type = "error", duration = 8)
          NULL

        }

      )

      if (is.null(plots) || !length(plots)) next

      results[[key]] <- list(label = dispname, plots = plots)

      img_files <- lapply(seq_along(plots), function(j) {

        plot_obj <- plots[[j]]
        dynamic_height <- attr(plot_obj, "dynamic_height")

        if (is.null(dynamic_height) || !is.finite(dynamic_height)) dynamic_height <- 7

        tmpfile <- tempfile(fileext = ".png")

        gg <- if (inherits(plot_obj, c("gg","ggplot"))) plot_obj else ggplotify::as.ggplot(plot_obj)

        ggsave(

          filename = tmpfile,
          plot     = gg,
          width    = w_run,
          height   = dynamic_height,
          dpi      = dpi_run,
          units    = "in",
          limitsize = FALSE
        )
        tmpfile

      })

      previews[[key]] <- unlist(img_files, use.names = FALSE)
    }

    regional_results(results)
    regional_previews(previews)

  })

  output$regional_dynamic_tabs <- renderUI({

    res <- regional_results(); prevs <- regional_previews()

    tabs <- lapply(names(res), function(k) {

      lbl <- res[[k]]$label

      tabPanel(

        title = lbl, value = k,

        div(

          slickROutput(paste0("regional_carousel_", k), width = "100%", height = "750px"),
          tags$br()

        )

      )

    })

    for (k in names(res)) local({

      kk <- k
      output[[paste0("regional_carousel_", kk)]] <- renderSlickR({

        imgs <- regional_previews()[[kk]]

        if (!length(imgs)) return(NULL)

        w <- slickR(imgs, slideType = "img") +
          settings(dots = TRUE, infinite = TRUE, arrows = TRUE, adaptiveHeight = FALSE)

        htmlwidgets::onRender(

          w,
          sprintf("
          function(el){
            var outId = 'regional_current_slide_%s';
            function push(idx){ if(window.Shiny){ Shiny.setInputValue(outId, idx+1, {priority:'event'}); } }
            var $el = window.jQuery ? window.jQuery(el) : null;
            if(!$el) return;
            // after the slick is initialized
            $el.on('init', function(e, slick){ push(slick.currentSlide || 0); });
            $el.on('afterChange', function(e, slick, current){ push(current); });
          }

        ", kk)

        )

      })

      output[[paste0("download_zip_png_regional_", kk)]] <- downloadHandler(

        filename = function() {

          paste0("Regional_", sanitize(regional_results()[[kk]]$label), "_plots_png.zip")

        },

        content  = function(file) {

          plots <- regional_results()[[kk]]$plots
          req(length(plots) > 0)

          w   <- reg_last_w()
          dpi <- reg_last_dpi()

          if (!is.finite(w) || w <= 0 || !is.finite(dpi) || dpi <= 0) {

            stop("Invalid saved width/DPI.")

          }

          tmpdir <- tempfile(paste0("regional_png_", kk, "_")); dir.create(tmpdir)
          nd <- nchar(as.character(length(plots)))

          files <- vapply(seq_along(plots), function(i) {

            p <- plots[[i]]
            h <- attr(p, "dynamic_height"); if (is.null(h) || !is.finite(h)) h <- 7
            f <- file.path(
              tmpdir,
              sprintf("Regional_%s_%0*d.png", sanitize(regional_results()[[kk]]$label), nd, i)
            )

            gg <- if (inherits(p, c("gg","ggplot"))) p else ggplotify::as.ggplot(p)

            ggsave(

              filename = f, plot = gg,
              width = w, height = h, dpi = dpi,
              units = "in", limitsize = FALSE

            )

            f

          }, character(1))

          zip_it(files, file)

        },

        contentType = "application/zip"#

      )

      # PDF ZIP per dataset (full replacement)

      output[[paste0("download_zip_pdf_regional_", kk)]] <- downloadHandler(

        filename = function() {

          paste0("Regional_", sanitize(regional_results()[[kk]]$label), "_plots_pdf.zip")

        },

        content  = function(file) {

          plots <- regional_results()[[kk]]$plots

          req(length(plots) > 0)

          w <- reg_last_w()

          if (!is.finite(w) || w <= 0) {

            stop("Invalid saved width.")

          }

          tmpdir <- tempfile(paste0("regional_pdf_", kk, "_")); dir.create(tmpdir)
          nd <- nchar(as.character(length(plots)))

          files <- vapply(seq_along(plots), function(i) {

            p <- plots[[i]]
            h <- attr(p, "dynamic_height"); if (is.null(h) || !is.finite(h)) h <- 7
            f <- file.path(
              tmpdir,
              sprintf("Regional_%s_%0*d.pdf", sanitize(regional_results()[[kk]]$label), nd, i)
            )
            gg <- if (inherits(p, c("gg","ggplot"))) p else ggplotify::as.ggplot(p)

            ggsave(
              filename = f, plot = gg,
              width = w, height = h,
              device = cairo_pdf
            )
            f
          }, character(1))

          zip_it(files, file)
        },
        contentType = "application/zip"
      )

    })
    do.call(tabsetPanel, c(list(id = "regional_output_tabs"), tabs))
  })

  # Annotate_Data

  annotate_defaults <- reactive({

    f_ann <- get("Annotate_Data"); fmls <- tryCatch(formals(f_ann), error = function(e) NULL)

    drops <- c("Data", "...", ".dots", "id", "input", "output", "session", "Debug")

    keep  <- if (!is.null(fmls)) setdiff(names(fmls), drops) else character(0)

    if (length(keep) == 0) {

      return(list(

        Chromosome_Column       = NULL,
        Position_Column         = NULL,
        SNP_ID_Column           = NULL,
        PValue_Column           = NULL,
        Reference_Allele_Column = NULL,
        Effect_Allele_Column    = NULL,
        Genome_Build            = "grch38",
        Verbose                 = FALSE

      ))

    }

    defs <- lapply(fmls[keep], function(x) safe_eval_default_reg(x, environment(f_ann))); names(defs) <- keep; defs

  })

  output$annotate_arg_inputs <- renderUI({

    lapply(names(annotate_defaults()), function(nm) make_input_like_meta("ann_", nm, annotate_defaults()[[nm]]))

  })

  annotate_arg_vals <- reactive({

    vals <- lapply(names(annotate_defaults()), function(nm) {

      val <- input[[paste0("ann_", nm)]]; if (is.null(val)) annotate_defaults()[[nm]] else parse_input_val_meta(val)

    })

    names(vals) <- names(annotate_defaults()); vals

  })

  build_call_preview_annotate <- function(arg_list, data_label = "user_data") {

    shown <- Filter(Negate(is.null), arg_list)

    parts <- vapply(names(shown), function(nm) paste0(nm, " = ", format_arg_val_meta(shown[[nm]])), "")

    paste0("Annotate_Data(\n  Data = ", data_label,

           if (length(parts)) paste0(",\n  ", paste(parts, collapse = ",\n  ")) else "",

           "\n)")

  }

  annotate_call_txt <- reactive({

    labels_vec <- if (!is.null(input$annotate_files$name)) file_path_sans_ext(basename(input$annotate_files$name)) else character(0)

    file_label <- if (length(labels_vec) == 1) labels_vec else if (length(labels_vec) > 1) "[multiple files]" else "user_data"

    build_call_preview_annotate(annotate_arg_vals(), data_label = file_label)

  })

  output$annotate_call_preview <- renderText(annotate_call_txt())

  observeEvent(input$run_annotate, {

    req(input$annotate_files)
    res <- list(); vals <- annotate_arg_vals()

    for (i in seq_len(nrow(input$annotate_files))) {

      datapath <- input$annotate_files$datapath[i]

      dispname <- file_path_sans_ext(input$annotate_files$name[i])

      key <- sanitize(dispname)

      df <- vroom::vroom(datapath, show_col_types = FALSE)

      assign(dispname, df, envir = .GlobalEnv)
      assign("user_data", get(dispname, envir = .GlobalEnv), envir = .GlobalEnv)
      call_txt <- build_call_preview_annotate(vals, data_label = dispname)

      out_df <- tryCatch({ eval(parse(text = call_txt), envir = .GlobalEnv) }, error = function(e) { showNotification(paste("Annotate_Data error (", dispname, "):", e$message), type="error", duration = 8); NULL })

      if (!is.null(out_df)) res[[key]] <- list(label = dispname, df = as.data.frame(out_df))

    }

    annotate_results(res)

  })

  output$annotate_dynamic_tabs <- renderUI({

    res <- annotate_results()

    tabs <- lapply(names(res), function(k) {

      lbl <- res[[k]]$label

      tabPanel(

        title = lbl, value = k,
        div(

          h4("Preview (first 10 rows)"),
          tableOutput(paste0("annotate_table_", k))

        )

      )

    })

    for (k in names(res)) local({

      kk <- k
      output[[paste0("annotate_table_", kk)]] <- renderTable({
        df <- annotate_results()[[kk]]$df

        if ("Lab" %in% names(df)) {

          lab_chr <- as.character(df$Lab)
          keep <- !is.na(lab_chr) & nzchar(trimws(lab_chr))
          df <- df[keep, , drop = FALSE]

        }

        df_preview <- head(df, 10)

        if (ncol(df_preview) > 7) {

          last_idx <- (ncol(df_preview) - 6):ncol(df_preview)
          df_preview <- df_preview[, last_idx, drop = FALSE]

        }

        df_preview
      }, striped = TRUE, bordered = TRUE, hover = TRUE, rownames = FALSE)

    })

    do.call(tabsetPanel, c(list(id = "annotate_output_tabs"), tabs))

  })

  .METASOFT_File_Gen_original <- tryCatch(

    getFromNamespace(".METASOFT_File_Gen_original", "MiamiR"),

    error = function(e) {

      tryCatch(getFromNamespace("METASOFT_File_Gen", "MiamiR"),

               error = function(e2) MiamiR::METASOFT_File_Gen)

    }

  )

  .METASOFT_File_Gen_original <-

    if (exists(".METASOFT_File_Gen_original")) get(".METASOFT_File_Gen_original")

  else METASOFT_File_Gen


  # METASOFT

  meta_defaults <- reactive({

    fmls  <- tryCatch(formals(.METASOFT_File_Gen_original), error = function(e) NULL)

    drops <- c("Data", "...", ".dots", "session", "input", "output")

    keep  <- if (!is.null(fmls)) setdiff(names(fmls), drops) else character(0)

    defs  <- lapply(fmls[keep], function(x)

      tryCatch(eval(x, envir = environment(.METASOFT_File_Gen_original)), error = function(e) NULL))

    names(defs) <- keep
    defs

  })

  output$meta_arg_inputs <- renderUI({

    defs <- meta_defaults()

    if (is.null(defs) || !length(defs))

      return(helpText("No additional METASOFT settings found in this build."))

    lapply(names(defs), function(nm) make_input_like_meta("meta_", nm, defs[[nm]]))

  })


  meta_arg_vals <- reactive({

    vals <- lapply(names(meta_defaults()), function(nm) {

      val <- input[[paste0("meta_", nm)]]; if (is.null(val)) meta_defaults()[[nm]] else parse_input_val_meta(val)

    })

    names(vals) <- names(meta_defaults()); vals

  })

  build_call_preview_meta <- function(arg_list, file_labels) {

    shown <- Filter(Negate(is.null), arg_list)

    parts <- vapply(names(shown), function(nm) paste0(nm, " = ", format_arg_val_meta(shown[[nm]])), "")
    paste0("METASOFT_File_Gen(\n  Data = c(", paste(paste0("\"", file_labels, "\""), collapse = ", "), ")",


           if (length(parts)) paste0(",\n  ", paste(parts, collapse = ",\n  ")) else "", "\n)")

  }

  output$meta_call_preview <- renderText({

    lbls <- if (!is.null(input$meta_files$name)) file_path_sans_ext(basename(input$meta_files$name)) else character(0)

    if (length(lbls) < 2) return("Need at least two files")

    build_call_preview_meta(meta_arg_vals(), lbls)

  })

  observeEvent(input$run_meta, {

    req(input$meta_files, nrow(input$meta_files) >= 2)
    lbls <- file_path_sans_ext(basename(input$meta_files$name))

    for (i in seq_along(input$meta_files$datapath)) {

      df <- vroom::vroom(input$meta_files$datapath[i], show_col_types = FALSE)
      objname <- lbls[i]; assign(objname, df, envir = .GlobalEnv)

    }

    args <- meta_arg_vals(); args$Data <- lbls

    out_df <- tryCatch({ do.call(.METASOFT_File_Gen_original, args) }, error = function(e) { showNotification(paste("Error:", e$message), type = "error"); NULL })

    if (!is.null(out_df)) meta_results(as.data.frame(out_df))

  })

  output$meta_output_preview <- renderTable({ head(meta_results(), 10) }, striped = TRUE, bordered = TRUE, hover = TRUE, rownames = FALSE)

  output$download_meta_csv <- downloadHandler(

    filename = function() "METASOFT_File_Gen_output.csv",

    content = function(file) { utils::write.csv(meta_results(), file, row.names = FALSE, na = "") },

    contentType = "text/csv"

  )

  #Forest_Plot

  observe({

    f <- get_forest_fun(); forest_fun(f)
    fmls <- formals(f); drops <- c("Data", "...", "session", "input", "output", "id", "Debug")
    keep  <- setdiff(names(fmls), drops)

    defs  <- lapply(fmls[keep], function(x) safe_eval_default_any(x, environment(f)))

    names(defs) <- keep; forest_defaults(defs)

  })

  # Forest_Plot: collect args from the settings

  collect_forest_args <- function() {

    fun <- tryCatch(.Forest_Plot_original, error = function(e) NULL)

    if (!is.function(fun)) fun <- get_forest_fun()

    sa   <- formals(fun)
    keep <- setdiff(names(sa), c("Data", "session", "..."))

    args <- list()

    for (arg in keep) {

      v <- parse_input_val_meta(input[[paste0("forest_", arg)]])
      if (!is.null(v)) args[[arg]] <- v

    }

    list(fun = fun, args = args)

  }

  # Forest_Plot: build settings panel inputs

  output$forest_arg_inputs <- renderUI({

    fun <- tryCatch(.Forest_Plot_original, error = function(e) NULL)
    if (!is.function(fun)) fun <- get_forest_fun()

    sa <- formals(fun)

    # Don't generate controls for these arguments

    drop <- c("Data", "session", "...")
    keep <- setdiff(names(sa), drop)

    lapply(

      keep,

      function(arg) make_input_like_meta(

        prefix  = "forest_",
        name    = arg,
        default = safe_eval_default_any(sa[[arg]], environment(fun))

      )

    )

  })


  observeEvent(input$forest_files, {

    req(input$forest_files$name)
    lbls <- tools::file_path_sans_ext(basename(input$forest_files$name))

    update_if_exists <- function(id, value) { if (!is.null(input[[id]])) updateTextInput(session, id, value = value) }

    as_vec_code <- function(chr) paste0("c(", paste(paste0("\"", chr, "\""), collapse = ", "), ")")

    n <- length(lbls)

    update_if_exists("forest_Names", as_vec_code(lbls))
    update_if_exists("forest_Shapes", as_vec_code(rep("square", n)))
    update_if_exists("forest_Styles", as_vec_code(rep("normal", n)))
    update_if_exists("forest_Data_Set_Colours", sprintf("viridis::viridis(%d)", n))

  }, ignoreInit = TRUE)

  forest_current_args <- reactive({

    defs <- forest_defaults(); if (is.null(defs)) return(list())

    vals <- lapply(names(defs), function(nm) parse_arg_val(input[[paste0("forest_", nm)]]))

    names(vals) <- names(defs); Filter(Negate(is.null), vals)

  })

  output$forest_call_preview <- renderText({

    fun <- tryCatch(.Forest_Plot_original, error = function(e) NULL)

    if (!is.function(fun)) fun <- get_forest_fun()

    sa <- formals(fun)

    drop <- c("Data", "session", "...")
    keep <- setdiff(names(sa), drop)

    args <- list()

    for (arg in keep) {

      v <- parse_input_val_meta(input[[paste0("forest_", arg)]])
      if (!is.null(v)) args[[arg]] <- v

    }

    shown <- Filter(Negate(is.null), args)
    parts <- vapply(names(shown),

                    function(nm) sprintf("%s = %s", nm, format_arg_val_meta(shown[[nm]])),

                    character(1))

    paste0("Forest_Plot(\n  Data = <file(s)>",
           if (length(parts)) paste0(",\n  ", paste(parts, collapse = ",\n  ")) else "",
           "\n)")
  })

  # Forest_Plot: run with user-specified args

  observeEvent(input$run_forest, {
    req(input$forest_files)

    fa          <- collect_forest_args()
    fun         <- fa$fun
    extra_args  <- fa$args

    paths       <- input$forest_files$datapath
    labels_vec  <- tools::file_path_sans_ext(input$forest_files$name)
    key_vec     <- vapply(labels_vec, sanitize, character(1))

    # Single call: pass files in one go

    args <- c(list(Data = paths), extra_args)

    forest_append_log("Running Forest_Plot with ", length(paths), " file(s) ...")

    out <- tryCatch(

      do.call(fun, args),

      error = function(e) { forest_append_log("ERROR: ", conditionMessage(e)); NULL }

    )
    req(!is.null(out))

    # Normalize output

    plots <- list()

    if (tube_is_plot_like(out)) {

      plots <- list(out)
      names(plots) <- key_vec[1]

    } else if (is.list(out)) {

      # keep only plot-like entries; if nested, flatten shallowly

      cand <- Filter(tube_is_plot_like, out)
      if (!length(cand) && any(vapply(out, is.list, logical(1)))) {
        cand <- Filter(tube_is_plot_like, unlist(out, recursive = FALSE))

      }

      plots <- cand

    }
    validate(need(length(plots) > 0, "Forest_Plot returned no plot objects."))

    # Build previews/tabs

    grobs <- list(); previews <- list(); labels <- list(); heights <- numeric(0); i <- 0

    for (nm in names(plots)) {

      i <- i + 1
      gp  <- plots[[nm]]
      key <- nm %||% key_vec[pmin(i, length(key_vec))]
      lab <-        labels_vec[pmin(i, length(labels_vec))]

      grobs[[key]]  <- gp
      labels[[key]] <- lab
      heights[key]  <- tube_get_dyn_height(gp, fallback = 8)

      tmp <- tempfile(fileext = ".png")

      ggsave(

        filename = tmp, plot = ggplotify::as.ggplot(gp),
        width = forest_last_w(), height = heights[key],
        dpi = forest_last_dpi(), units = "in", limitsize = FALSE

      )
      previews[[key]] <- tmp
    }

    forest_grobs(grobs)
    forest_previews(previews)
    forest_labels(labels)
    forest_heights(heights)
    forest_append_log("Done.")
  }, ignoreInit = TRUE)



  output$forest_dynamic_tabs <- renderUI({

    grobs <- forest_grobs(); previews <- forest_previews(); labels <- forest_labels()


    tabs <- lapply(names(grobs), function(k) {


      tabPanel(
        title = labels[[k]], value = k,
        div(imageOutput(paste0("forest_plot_", k), height = "auto"))

      )

    })

    for (k in names(previews)) local({

      kk <- k

      output[[paste0("forest_plot_", kk)]] <- renderImage({

        list(src = forest_previews()[[kk]], contentType = "image/png", width = "100%")

      }, deleteFile = FALSE)

    })

    do.call(tabsetPanel, c(list(id = "forest_output_tabs"), tabs))

  })

  output$forest_log <- renderText(forest_log())

  output$forest_download_simple <- downloadHandler(

    filename = function() {

      ext <- tolower(input$forest_ext %||% "png")
      n   <- length(forest_grobs())

      if (n <= 1) {

        cur <- input$forest_output_tabs

        nm  <- if (!is.null(cur) && cur %in% names(forest_labels()))

          forest_labels()[[cur]] else "Forest"

        paste0("Forest_", sanitize(nm), ".", ext)

      } else {

        paste0("Forest_ALL_", ext, ".zip")

      }

    },

    content = function(file) {

      grobs    <- forest_grobs()
      previews <- forest_previews()
      labels   <- forest_labels()
      heights  <- forest_heights()
      ext      <- tolower(input$forest_ext %||% "png")
      n        <- length(grobs)
      req(n > 0)


      if (ext == "png") {

        if (n == 1) {

          cur <- input$forest_output_tabs

          key <- if (!is.null(cur) && cur %in% names(previews)) cur else names(previews)[1]

          file.copy(previews[[key]], file, overwrite = TRUE)

        } else {

          tmpdir <- tempfile("forest_all_png_"); dir.create(tmpdir)
          files <- vapply(names(previews), function(k) {
            dst <- file.path(tmpdir, paste0("Forest_", sanitize(labels[[k]]), ".png"))
            file.copy(previews[[k]], dst, overwrite = TRUE); dst
          }, character(1))

          zip_it(files, file)
        }

        return(invisible())

      }

      write_one <- function(g, path, h) {

        gg <- if (inherits(g, c("gg","ggplot"))) g else ggplotify::as.ggplot(g)

        if (ext == "jpg") {

          ggsave(path, gg, device = "jpeg",
                 width = forest_last_w(), height = h,
                 dpi = forest_last_dpi(), units = "in", limitsize = FALSE)

        } else if (ext == "pdf") {

          ggsave(path, gg, device = cairo_pdf,
                 width = forest_last_w(), height = h)

        } else if (ext == "svg") {

          if (!requireNamespace("svglite", quietly = TRUE)) stop("SVG export requires 'svglite'.")

          ggsave(path, gg, device = svglite::svglite,
                 width = forest_last_w(), height = h)

        } else {

          ggsave(path, gg,
                 width = forest_last_w(), height = h,
                 dpi = forest_last_dpi(), units = "in", limitsize = FALSE)

        }

      }

      if (n == 1) {

        cur <- input$forest_output_tabs

        key <- if (!is.null(cur) && cur %in% names(grobs)) cur else names(grobs)[1]

        h   <- as.numeric(heights[[key]]); if (!is.finite(h)) h <- 8
        write_one(grobs[[key]], file, h)

      } else {

        tmpdir <- tempfile("forest_all_"); dir.create(tmpdir)

        files <- vapply(names(grobs), function(k) {

          h  <- as.numeric(heights[[k]]); if (!is.finite(h)) h <- 8
          fn <- file.path(tmpdir, paste0("Fo
                                         rest_", sanitize(labels[[k]]), ".", ext))
          write_one(grobs[[k]], fn, h); fn
        }, character(1))

        zip_it(files, file)

      }

    }

  )


  output$forest_plot_img <- renderImage({ req(forest_preview_png()); list(src = forest_preview_png(), contentType = "image/png", width = "100%") }, deleteFile = FALSE)

  output$forest_log <- renderText(forest_log())

  make_forest_base <- function() if (!is.null(input$forest_files$name)) paste0("Forest_", paste0(tools::file_path_sans_ext(basename(input$forest_files$name)), collapse = "_")) else "Forest"

  output$download_forest_png <- downloadHandler(

    filename = function() paste0(make_forest_base(), ".png"),

    content = function(file) {

      gp <- forest_result_grob(); req(gp); h <- forest_dynamic_height()

      if (inherits(gp, c("gg","ggplot"))) ggsave(file, plot = gp, width = input$forest_save_width, height = h, dpi = input$forest_save_dpi, units = "in", limitsize = FALSE)

      else { gg <- ggplotify::as.ggplot(gp); ggsave(file, plot = gg, width = input$forest_save_width, height = h, dpi = input$forest_save_dpi, units = "in", limitsize = FALSE) }

    }, contentType = "image/png"
  )

  output$download_forest_pdf <- downloadHandler(

    filename = function() paste0(make_forest_base(), ".pdf"),

    content = function(file) {

      gp <- forest_result_grob(); req(gp)

      if (inherits(gp, c("gg","ggplot"))) ggsave(file, plot = gp, width = input$forest_save_width, height = forest_dynamic_height(), device = cairo_pdf)

      else { gg <- ggplotify::as.ggplot(gp); ggsave(file, plot = gg, width = input$forest_save_width, height = forest_dynamic_height(), device = cairo_pdf) }

    }, contentType = "application/pdf"
  )
  output$download_forest_rds <- downloadHandler(

    filename = function() paste0(make_forest_base(), ".rds"),

    content = function(file) { gp <- forest_result_grob(); req(gp); saveRDS(gp, file) }

  )

  #Model_Munge

  mm_defaults <- reactive({

    fmls <- tryCatch(formals(.Model_Munge_original), error = function(e) NULL)

    drops <- c("Model_Object", "...", "session", "input", "output", "id", "Debug")

    keep  <- if (!is.null(fmls)) setdiff(names(fmls), drops) else character(0)

    defs  <- lapply(fmls[keep], function(x) { tryCatch(eval(x, envir = environment(.Model_Munge_original)), error = function(e) NULL) })

    names(defs) <- keep; defs
  })

  # Model_Munge: build settings panel inputs

  output$mm_arg_inputs <- renderUI({

    # Prefer the original if present; otherwise use the public function

    fun <- tryCatch(.Model_Munge_original, error = function(e) NULL)

    if (!is.function(fun)) fun <- Model_Munge

    if (!is.function(fun)) {

      return(div(helpText("Model_Munge() not found.")))

    }

    sa   <- formals(fun)

    drop <- c("Data", "Model", "session", "...")

    keep <- setdiff(names(sa), drop)

    # If there are no configurable args, show a helpful note

    if (length(keep) == 0) {

      return(tagList(

        div(class = "text-muted",
            "Data/Model already considered")

      ))

    }

    lapply(

      keep,

      function(arg) make_input_like_meta(

        prefix  = "mm_",
        name    = arg,
        default = safe_eval_default_any(sa[[arg]], environment(fun))

      )

    )

  })

  mm_arg_vals <- reactive({

    defs <- mm_defaults(); if (is.null(defs)) return(list())

    vals <- lapply(names(defs), function(nm) {

      val <- input[[paste0("mm_arg_", nm)]]; if (is.null(val)) defs[[nm]] else parse_input_val_meta(val)

    })

    names(vals) <- names(defs); vals

  })


  mm_call_txt <- reactive({

    model_name <- input$mm_model_name

    if (!nzchar(model_name)) model_name <- "<model_object>"

    args <- mm_arg_vals()
    shown <- Filter(Negate(is.null), args)

    parts <- vapply(
      names(shown),
      function(nm) paste0(nm, " = ", format_arg_val_meta(shown[[nm]])),

      character(1)

    )

    paste0(

      "Model_Munge(\n  Model_Object = ", model_name,

      if (length(parts)) paste0(",\n  ", paste(parts, collapse = ",\n  ")) else "",

      "\n)"
    )
  })

  output$mm_call_preview <- renderText(mm_call_txt())

  # Keep the uploaded dataset around for the builder/preview
  observeEvent(input$mm_file, {
    req(input$mm_file)
    df <- vroom::vroom(input$mm_file$datapath, show_col_types = FALSE)
    mm_df(df)
  }, ignoreInit = TRUE)

  # --- Model builder (response/predictors) ---
  output$mm_var_ui <- renderUI({
    df <- mm_df()
    if (is.null(df)) return(div(helpText("Upload a data file to pick variables.")))

    cols <- names(df)
    tagList(
      selectInput("mm_response",   "Response (y)",  choices = cols),
      selectInput("mm_predictors", "Predictors (x)", choices = cols, multiple = TRUE)
    )
  })


  output$mm_data_preview <- renderTable({

    df <- mm_df(); if (is.null(df)) return(NULL)

    head(df, 10)

  }, striped = TRUE, bordered = TRUE, hover = TRUE, rownames = FALSE)

  output$mm_var_ui <- renderUI({

    df <- mm_df(); if (is.null(df)) return(NULL)

    cols <- names(df)
    tagList(
      selectInput("mm_response", "Response (y)", choices = cols, selected = cols[1]),
      selectizeInput("mm_predictors", "Predictors (x)", choices = cols, selected = setdiff(cols, cols[1])[1], multiple = TRUE)

    )

  })

  mm_formula_built <- reactive({

    df <- mm_df(); req(df); y <- input$mm_response; xs <- input$mm_predictors

    if (is.null(y) || !length(xs)) return("")

    rhs <- paste(xs, collapse = " + ")

    if (!isTRUE(input$mm_intercept)) rhs <- paste("0 +", rhs)

    paste(y, "~", rhs)

  })


  observeEvent(input$mm_fit, {

    df <- mm_df(); req(df)
    ds_name <- input$mm_dataset_name; if (!nzchar(ds_name)) ds_name <- "mm_data"

    assign(ds_name, df, envir = .GlobalEnv)

    formula_str <- trimws(input$mm_formula_manual)

    if (!nzchar(formula_str)) formula_str <- mm_formula_built()

    req(nzchar(formula_str))

    form <- tryCatch(as.formula(formula_str), error = function(e) { showNotification(paste("Invalid formula:", e$message), type="error"); NULL })

    req(!is.null(form))

    model_name <- input$mm_model_name; if (!nzchar(model_name)) model_name <- paste0(ds_name, "_model")

    fit <- tryCatch({

      if (identical(input$mm_model_type, "glm")) {

        fam_fun <- switch(input$mm_glm_family,
                          gaussian = gaussian, binomial = binomial, poisson = poisson,
                          Gamma = Gamma, `inverse.gaussian` = inverse.gaussian,
                          quasipoisson = quasipoisson, quasibinomial = quasibinomial, gaussian)
        glm(form, data = get(ds_name, envir = .GlobalEnv), family = fam_fun())

      } else {


        lm(form, data = get(ds_name, envir = .GlobalEnv))

      }

    }, error = function(e) { showNotification(paste("Model fit error:", e$message), type="error"); NULL })

    req(!is.null(fit))

    assign(model_name, fit, envir = .GlobalEnv)
    mm_model(fit)
    showNotification(paste0("Model '", model_name, "' saved in .GlobalEnv"), type = "message")

  })


  output$mm_model_summary <- renderPrint({

    fit <- mm_model(); if (is.null(fit)) return(cat("Fit a model first."))
    summary(fit)

  })

  observeEvent(input$mm_run_munge, {

    model_name <- input$mm_model_name; req(nzchar(model_name))
    if (!exists(model_name, envir = .GlobalEnv)) {
      showNotification(paste0("Object '", model_name, "' not found in .GlobalEnv. Fit model first."), type = "error"); return(NULL)
    }

    args <- mm_arg_vals()
    args$Model_Object <- model_name
    out_df <- tryCatch({ do.call(.Model_Munge_original, args) },
                       error = function(e) { showNotification(paste("Model_Munge error:", e$message), type="error"); NULL })

    if (!is.null(out_df)) mm_munge_df(as.data.frame(out_df))

  })

  output$mm_munge_preview <- renderTable({

    df <- mm_munge_df(); if (is.null(df)) return(NULL)

    head(df, 10)

  }, striped = TRUE, bordered = TRUE, hover = TRUE, rownames = FALSE)

  output$mm_download_csv <- downloadHandler(

    filename = function() {

      nm <- if (nzchar(input$mm_model_name)) sanitize(input$mm_model_name) else "Model_Munge"

      paste0(nm, "_munged.csv")

    },
    content = function(file) {

      df <- mm_munge_df(); req(df)
      utils::write.csv(df, file, row.names = FALSE, na = "")

    },
    contentType = "text/csv"
  )

  output$mm_download_simple <- downloadHandler(

    filename = function() {

      ext <- tolower(input$mm_ext %||% "csv")

      nm  <- if (nzchar(input$mm_model_name)) sanitize(input$mm_model_name) else "Model_Munge"

      paste0(nm, "_munged.", ext)
    },

    content = function(file) {

      df <- mm_munge_df(); req(df)
      ext <- tolower(input$mm_ext %||% "csv")

      if (ext == "csv") {

        utils::write.csv(df, file, row.names = FALSE, na = "")

      } else if (ext %in% c("tsv", "txt")) {

        utils::write.table(df, file, sep = "\t",
                           row.names = FALSE, col.names = TRUE,
                           quote = FALSE, na = "")

      } else {


        utils::write.csv(df, file, row.names = FALSE, na = "")

      }

    }

  )


}

shinyApp(ui, server)




