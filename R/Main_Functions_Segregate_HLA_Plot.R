
#' Annotated, Customised HLA-WAS Plots
#'
#' @param Data These are the summary statistics (either file path(s), environmental object(s) or a variable(s) containing such contents) to be plotted; defaults to NULL
#' @param Verbose Prevent display of progress bar as function is running and instead show key milestone outputs/messages (mainly for debugging purposes); defaults to FALSE
#' @param Title Manually specify title to be displayed; defaults to NULL, leading to auto-detection of dataframe object/file name
#' @param Title_On Toggle on/off whether Title (either manually or automatically determined) displays on final plot; defaults to TRUE
#' @param Title_Size Size of title; defaults to 28
#' @param X_Axis_Title Manually specify X axis title; defaults to "HLA Genes"
#' @param Y_Axis_Title Manually specify Y axis title; defaults to "-log\u2081\u2080(P)" - shows as "-log₁₀(P)"
#' @param X_Axis_Title_Size Size of X axis title; defaults to 26
#' @param Y_Axis_Title_Size Size of Y axis title; defaults to 26
#' @param X_Axis_Title_Vjust Number of units of blank space between X_Axis_Title element area and above axis; defaults to -3
#' @param Y_Axis_Title_Vjust Number of units of blank space between Y_Axis_Title element area and adjacent axis; defaults to 6
#' @param X_Axis_Text_Vjust Number of units of blank space between X_Axis_Text element area and above axis; defaults to -0.2
#' @param Y_Axis_Text_Vjust Number of units of blank space between Y_Axis_Text element area and adjacent axis; defaults to 0.2
#' @param X_Axis_Text_Size Size of X axis text; defaults to 20
#' @param Y_Axis_Text_Size Size of Y axis text; defaults to 20
#' @param Label_Index Annotate the index SNPs with ID provided; defaults to TRUE
#' @param Label_Size Size of index labels if Label_Index is TRUE; defaults to 3
#' @param Label_Height Spacing between the point and the label position; defaults to 2
#' @param Segment_Width Half-width of each per-point horizontal segment (total length = 2 * Segment_Width); defaults to 0.18
#' @param Segment_Thickness Thickness/line width for the per-point horizontal segments; defaults to 0.9
#' @param Sig_Line_Colour Colour of significance threshold line; defaults to "red"
#' @param Sig_Line_Type Type of significance threshold line; defaults to "dashed"
#' @param Sig_Threshold Value of significance threshold; defaults to standard Bonferroni correction of 5e-8
#' @param Sig_Line_Width Width of significance threshold line; defaults to 0.6
#'
#' @returns Image of HLA type-specific Line Plot(s) is allocated to specified object and the resulting ggplots object can then be saved to images
#' @export
#'
#' @examples HLA_Plots <- Segregate_HLA_Plot(Data = "Fake_HLA_100K")
#'

        Segregate_HLA_Plot <- function(

        Data = NULL,
        Verbose = FALSE,
        Title = NULL,
        Title_On = TRUE,
        Title_Size = 28,
        X_Axis_Title = "HLA Genes",
        Y_Axis_Title = "-log\u2081\u2080(P)",
        X_Axis_Title_Size = 26,
        Y_Axis_Title_Size = 26,
        X_Axis_Title_Vjust = -3,
        Y_Axis_Title_Vjust = 6,
        X_Axis_Text_Vjust = -0.2,
        Y_Axis_Text_Vjust = 0.2,
        X_Axis_Text_Size = 20,
        Y_Axis_Text_Size = 20,
        Label_Index  = TRUE,
        Label_Size   = 3,
        Label_Height = 2,
        Segment_Width     = 0.18,
        Segment_Thickness = 0.9,
        Sig_Line_Colour = "red",
        Sig_Line_Type   = "dashed",
        Sig_Threshold   = 5e-8,
        Sig_Line_Width  = 0.6

    ) {

      # Progress hooks

      if (!exists(".progress", inherits = TRUE))       .progress        <- function(...) invisible(NULL)

      if (!exists(".progress_done", inherits = TRUE))  .progress_done   <- function(...) invisible(NULL)

      if (!exists(".progress_pause", inherits = TRUE)) .progress_pause  <- function(...) invisible(NULL)

      if (!exists(".progress_resume", inherits = TRUE)).progress_resume <- function(...) invisible(NULL)

      .allow_msgs <- isTRUE(Verbose)

      message <- function(..., domain = NULL, appendLF = TRUE) {

        if (isTRUE(get0(".allow_msgs", ifnotfound = FALSE, inherits = TRUE))) {

          base::message(..., domain = domain, appendLF = appendLF)

        }

        invisible(NULL)

      }

      # Silence message() calls

      .quiet_call <- function(fun, ..., verbose = isTRUE(Verbose)) {

        if (verbose) return(fun(...))

        withCallingHandlers(fun(...), message = function(m) invokeRestart("muffleMessage"))

      }

      # Pause progress bar

      .progress_break <- function(expr, msg, pad = 1L) {

        .progress_pause()
        cat("\r\n", sep = "")
        cat(strrep("\n", pad), msg, strrep("\n\n", pad), "\n", sep = "")
        utils::flush.console()
        on.exit(.progress_resume(), add = TRUE)
        force(expr)

      }

      message("Arguments received, beginning Segregate_HLA_Plot()")

      if (is.null(Data)) stop("Please provide Data.", call. = FALSE)

      pheno <- "PHENO"

      message("Resolving input (data.frame / file path / object name)")

      if (is.character(Data) && length(Data) == 1L) {

        if (file.exists(Data)) {

          file_path <- Data
          pheno <- tools::file_path_sans_ext(basename(file_path))

          message("Dataset absent from environment")

          message("Reading data from file: ", file_path)

          message("Loading data using vroom...")

          Data <- suppressMessages(

            suppressWarnings(

              vroom::vroom(file_path, show_col_types = FALSE, progress = FALSE)

            )

          )

          message("Finished reading")

        } else {

          if (exists(Data, inherits = TRUE)) {

            pheno <- Data
            Data  <- get(Data, inherits = TRUE)

            message("Input resolved to object: ", pheno)

          } else {

            stop("The provided character string does not point to an existing file or object: ", Data, call. = FALSE)

          }

        }

      } else {

        stopifnot(is.data.frame(Data))

        # Prefer wrapper

        pheno <- attr(Data, ".pheno", exact = TRUE)

        if (is.null(pheno) || !nzchar(pheno)) {

          pheno <- tryCatch(deparse(substitute(Data)), error = function(...) "PHENO")

        }

        message("Input is a data.frame: ", pheno)

      }

      # detect key columns

      message("Detecting key columns")

      ID_Column <- .quiet_call(detect_snp_column, Data, NULL)
      P_Column  <- .quiet_call(detect_pvalue_column, Data, NULL)

      if (is.null(ID_Column) || is.null(P_Column)) {

        stop(

          "Could not detect required columns (ID + P/LOG10P). Provide standard columns or pass a harmonised dataset.",
          call. = FALSE

        )

      }

      # harmonise to ID + P + LOG10P - don't need all usuals

      message("Harmonising columns")

      df <- Data

      if (!identical(ID_Column, "ID")) {

        df[["ID"]] <- df[[ID_Column]]

      }

      p_vec <- df[[P_Column]]

      p_vec <- if (is.numeric(p_vec)) p_vec else as.numeric(p_vec)

      if (grepl("log", P_Column, ignore.case = TRUE)) {

        message("Converting log value to plain P (Single_Plot-style)")

        P_raw <- exp(-p_vec * 2.302585092994046)

      } else {

        P_raw <- p_vec

      }

      df[["P"]]      <- P_raw
      df[["LOG10P"]] <- -log10(P_raw)

      # Checks

      message("Checking required columns after harmonisation")

      req_cols <- c("ID", "LOG10P")
      miss <- setdiff(req_cols, names(df))

      if (length(miss)) stop("Data is missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)

      message("Required columns present")

      # helpers for hla tpyes

      get_hla_type <- function(x) stringr::str_match(x, "^HLA_([^*]+)")[, 2]

      get_aa_cat   <- function(x) stringr::str_match(x, "^(?:AA_|INDEL_AA_)([^_]+)_")[, 2]

      get_snps_cat <- function(x) stringr::str_match(x, "^(?:SNPS_|INDEL_SNPS_)([^_]+)_")[, 2]

      is_rsid <- function(x) stringr::str_detect(x, "^rs\\d+")

      is_hla  <- function(x) stringr::str_starts(x, "HLA_")

      is_aa   <- function(x) stringr::str_starts(x, "AA_")   | stringr::str_starts(x, "INDEL_AA_")

      is_snps <- function(x) stringr::str_starts(x, "SNPS_") | stringr::str_starts(x, "INDEL_SNPS_")

      make_pal <- function(n) scales::hue_pal(l = 55, c = 100)(n)

      # sig threshold (from args)

      thr <- -log10(Sig_Threshold)

      # classify

      message("Classifying variants into HLA / AA / SNPS / RSID")

      classed <- dplyr::mutate(

        df,

        CLASS = dplyr::case_when(

          is_hla(ID)  ~ "HLA",
          is_aa(ID)   ~ "AA",
          is_snps(ID) ~ "SNPS",
          is_rsid(ID) ~ "RSID",
          TRUE        ~ "OTHER"

        )

      )

      # classification counts out of total rows + non-classified

      n_total <- nrow(classed)
      cls_tab <- as.data.frame(table(classed$CLASS), stringsAsFactors = FALSE)
      names(cls_tab) <- c("CLASS", "N")
      ord <- match(cls_tab$CLASS, c("HLA", "AA", "SNPS", "RSID", "OTHER"))
      cls_tab <- cls_tab[order(ord, na.last = TRUE), , drop = FALSE]
      counts_str <- paste0(cls_tab$CLASS, "=", cls_tab$N, collapse = " | ")

      message("Classification: ", counts_str)

      message("Classified rows: ", sum(cls_tab$N), " / ", n_total)

      other_df <- dplyr::filter(classed, CLASS == "OTHER")

      if (nrow(other_df) > 0) {

        message("Unmatched rows (CLASS == 'OTHER') head():")

        show_cols <- intersect(c("ID", "LOG10P", "P", "CLASS"), names(other_df))

        if (length(show_cols) == 0) show_cols <- names(other_df)

        print(utils::head(other_df[, show_cols, drop = FALSE], 10))

      } else {

        message("No unmatched rows (CLASS == 'OTHER').")

      }

      # split for plots

      df_hla  <- dplyr::filter(classed, CLASS == "HLA") |>

        dplyr::mutate(

          HLA_TYPE = get_hla_type(ID),
          Allele   = sub("^HLA_[^*]*\\*", "", ID)

        )

      df_aa   <- dplyr::filter(classed, CLASS == "AA")   |> dplyr::mutate(AA_CAT = get_aa_cat(ID))
      df_snps <- dplyr::filter(classed, CLASS == "SNPS") |> dplyr::mutate(SNPS_CAT = get_snps_cat(ID))
      df_rsid <- dplyr::filter(classed, CLASS == "RSID")

      unique_hla_types <- sort(unique(df_hla$HLA_TYPE))
      unique_aa_cats   <- sort(unique(df_aa$AA_CAT))
      unique_snps_cats <- sort(unique(df_snps$SNPS_CAT))

      message(

        paste0(

          "Counts: HLA loci=", length(unique_hla_types),
          " | AA sets=", length(unique_aa_cats),
          " | SNP sets=", length(unique_snps_cats),
          " | RSIDs=", nrow(df_rsid)

        )

      )

      title_prefix <- if (is.null(Title) || !nzchar(Title)) pheno else Title

      .base_theme <- function() {

        ggplot2::theme_classic() +

          ggplot2::theme(

            plot.title   = ggplot2::element_text(hjust = 0.5, size = Title_Size, colour = "black"),
            axis.title.x = ggplot2::element_text(size = X_Axis_Title_Size, vjust = X_Axis_Title_Vjust, colour = "black"),
            axis.title.y = ggplot2::element_text(size = Y_Axis_Title_Size, vjust = Y_Axis_Title_Vjust, colour = "black"),
            axis.text.x  = ggplot2::element_text(
              size = X_Axis_Text_Size, vjust = X_Axis_Text_Vjust, colour = "black",
              margin = ggplot2::margin(t = 6)

            ),

            axis.text.y  = ggplot2::element_text(size = Y_Axis_Text_Size, vjust = Y_Axis_Text_Vjust, colour = "black"),
            axis.line    = ggplot2::element_line(colour = "black"),
            plot.margin  = ggplot2::margin(t = 14, r = 14, b = 12, l = 12)

          )
      }


      # HLA

      p_hla <- NULL

      if (length(unique_hla_types)) {

        message("Building HLA locus plot")

        pal_hla <- make_pal(length(unique_hla_types))

        hla_plot_df <- dplyr::mutate(

          df_hla,
          Locus   = factor(HLA_TYPE, levels = unique_hla_types),
          x_num   = as.numeric(Locus),
          x_left  = x_num - Segment_Width,
          x_right = x_num + Segment_Width,
          LABEL   = sub("^HLA_[^*]*\\*", "", ID)

        )

        hla_labels <- hla_plot_df |>
          dplyr::group_by(Locus) |>
          dplyr::slice_max(LOG10P, n = 1, with_ties = FALSE) |>
          dplyr::ungroup() |>
          dplyr::mutate(label_x = x_num, label_y = LOG10P)

        # Keep forced label behaviour for future

        # hla_labels_force <- dplyr::filter(hla_plot_df, ID %in% c("HLA_B*08", "HLA_B*44")) |>
        #   dplyr::mutate(label_x = x_num, label_y = LOG10P)
        #
        # if (nrow(hla_labels_force)) {
        #
        #   hla_labels <- dplyr::bind_rows(
        #
        #     hla_labels,
        #     dplyr::anti_join(hla_labels_force, hla_labels, by = c("ID"))
        #
        #   )
        #
        #   message("Added forced labels for: ", paste(hla_labels_force$ID, collapse = ", "))
        #
        # }

        hla_title <- paste0(title_prefix, " — Classical HLA Alleles")

        if (!isTRUE(Title_On)) hla_title <- NULL

        p_hla <- ggplot2::ggplot(hla_plot_df, ggplot2::aes(color = Locus)) +
          ggplot2::geom_segment(
            ggplot2::aes(x = x_left, xend = x_right, y = LOG10P, yend = LOG10P),
            linewidth = Segment_Thickness,
            alpha = 0.55
          ) +

          ggplot2::geom_hline(

            yintercept = thr,
            linetype  = Sig_Line_Type,
            colour    = Sig_Line_Colour,
            linewidth = Sig_Line_Width

          ) +

          ggplot2::geom_text(

            data = if (isTRUE(Label_Index)) hla_labels else NULL,
            ggplot2::aes(x = label_x, y = label_y, label = LABEL),
            hjust  = 0.5,
            vjust  = -Label_Height,
            size   = Label_Size,
            colour = "black",
            alpha  = 0.9,
            show.legend = FALSE

          ) +

          ggplot2::scale_color_manual(values = pal_hla, guide = "none") +

          ggplot2::scale_x_continuous(

            breaks = seq_along(unique_hla_types),
            labels = unique_hla_types,
            expand = ggplot2::expansion(mult = c(0.05, 0.10))

          ) +

          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.18))) +
          ggplot2::labs(x = X_Axis_Title, y = Y_Axis_Title, title = hla_title) +
          ggplot2::coord_cartesian(clip = "off") +

          .base_theme()

        attr(p_hla, "recommended_width")  <- 12
        attr(p_hla, "recommended_height") <- max(6, length(unique_hla_types) * 0.6)
        attr(p_hla, "suggested_filename") <- paste0(pheno, "_HLA_Locus.jpg")

        message("HLA plot complete")

      }

      else {

        message("No HLA alleles detected; skipping HLA plot")

      }


      # AA plot

      p_aa <- NULL

      if (length(unique_aa_cats)) {

        message("Building AA sets plot")

        pal_aa <- make_pal(length(unique_aa_cats))

        aa_plot_df <- dplyr::mutate(

          df_aa,
          AA_CAT  = factor(AA_CAT, levels = unique_aa_cats),
          x_num   = as.numeric(AA_CAT),
          x_left  = x_num - Segment_Width,
          x_right = x_num + Segment_Width

        )

        aa_labels <- aa_plot_df |>
          dplyr::group_by(AA_CAT) |>
          dplyr::slice_max(LOG10P, n = 1, with_ties = FALSE) |>
          dplyr::ungroup() |>
          dplyr::mutate(

            LABEL = dplyr::case_when(

              stringr::str_starts(ID, "INDEL_AA_") ~ sub("^INDEL_AA_[^_]+_(.*)$", "\\1", ID),
              TRUE                                 ~ sub("^AA_[^_]+_(.*)$", "\\1", ID)

            ),

            label_x = x_num,
            label_y = LOG10P

          )

        aa_title <- paste0(title_prefix, " — AA Sets")

        if (!isTRUE(Title_On)) aa_title <- NULL

        p_aa <- ggplot2::ggplot(aa_plot_df, ggplot2::aes(color = AA_CAT)) +

          ggplot2::geom_segment(

            ggplot2::aes(x = x_left, xend = x_right, y = LOG10P, yend = LOG10P),
            linewidth = Segment_Thickness,
            alpha = 0.55

          ) +
          ggplot2::geom_hline(

            yintercept = thr,
            linetype  = Sig_Line_Type,
            colour    = Sig_Line_Colour,
            linewidth = Sig_Line_Width

          ) +
          ggplot2::geom_text(

            data = if (isTRUE(Label_Index)) aa_labels else NULL,
            ggplot2::aes(x = label_x, y = label_y, label = LABEL),
            hjust  = 0.5,
            vjust  = -Label_Height,
            size   = Label_Size,
            colour = "black",
            alpha  = 0.9,
            show.legend = FALSE

          ) +

          ggplot2::scale_color_manual(values = pal_aa, guide = "none") +
          ggplot2::scale_x_continuous(

            breaks = seq_along(unique_aa_cats),
            labels = unique_aa_cats,
            expand = ggplot2::expansion(mult = c(0.05, 0.10))

          ) +

          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.18))) +
          ggplot2::labs(x = X_Axis_Title, y = Y_Axis_Title, title = aa_title) +
          ggplot2::coord_cartesian(clip = "off") +
          .base_theme()

        attr(p_aa, "recommended_width")  <- 12
        attr(p_aa, "recommended_height") <- max(6, length(unique_aa_cats) * 0.6)
        attr(p_aa, "suggested_filename") <- paste0(pheno, "_AA_Categories.jpg")

        message("AA plot complete")

      } else {

        message("No AA sets detected; skipping AA plot")

      }

      # SNPS plot

      p_snps <- NULL

      if (length(unique_snps_cats)) {

        message("Building SNP sets plot")

        pal_snps <- make_pal(length(unique_snps_cats))

        snps_plot_df <- dplyr::mutate(
          df_snps,
          SNPS_CAT = factor(SNPS_CAT, levels = unique_snps_cats),
          x_num    = as.numeric(SNPS_CAT),
          x_left   = x_num - Segment_Width,
          x_right  = x_num + Segment_Width
        )

        snps_labels <- snps_plot_df |>
          dplyr::group_by(SNPS_CAT) |>
          dplyr::slice_max(LOG10P, n = 1, with_ties = FALSE) |>
          dplyr::ungroup() |>
          dplyr::mutate(

            LABEL = dplyr::case_when(
              stringr::str_starts(ID, "INDEL_SNPS_") ~ sub("^INDEL_SNPS_[^_]+_(.*)$", "\\1", ID),
              TRUE                                    ~ sub("^SNPS_[^_]+_(.*)$", "\\1", ID)

            ),

            label_x = x_num,
            label_y = LOG10P

          )

        snps_title <- paste0(title_prefix, " — SNP Sets")

        if (!isTRUE(Title_On)) snps_title <- NULL

        p_snps <- ggplot2::ggplot(snps_plot_df, ggplot2::aes(color = SNPS_CAT)) +
          ggplot2::geom_segment(

            ggplot2::aes(x = x_left, xend = x_right, y = LOG10P, yend = LOG10P),
            linewidth = Segment_Thickness,
            alpha = 0.55

          ) +
          ggplot2::geom_hline(

            yintercept = thr,
            linetype  = Sig_Line_Type,
            colour    = Sig_Line_Colour,
            linewidth = Sig_Line_Width

          ) +
          ggplot2::geom_text(

            data = if (isTRUE(Label_Index)) snps_labels else NULL,
            ggplot2::aes(x = label_x, y = label_y, label = LABEL),
            hjust  = 0.5,
            vjust  = -Label_Height,
            size   = Label_Size,
            colour = "black",
            alpha  = 0.9,
            show.legend = FALSE

          ) +

          ggplot2::scale_color_manual(values = pal_snps, guide = "none") +

          ggplot2::scale_x_continuous(

            breaks = seq_along(unique_snps_cats),
            labels = unique_snps_cats,
            expand = ggplot2::expansion(mult = c(0.05, 0.10))

          ) +

          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.18))) +
          ggplot2::labs(x = X_Axis_Title, y = Y_Axis_Title, title = snps_title) +
          ggplot2::coord_cartesian(clip = "off") +
          .base_theme()

        attr(p_snps, "recommended_width")  <- 12
        attr(p_snps, "recommended_height") <- max(6, length(unique_snps_cats) * 0.6)
        attr(p_snps, "suggested_filename") <- paste0(pheno, "_SNPS_Categories.jpg")

        message("SNP plot complete")

      }

      else {

        message("No SNP sets detected; skipping SNP plot")

      }

      # RSID Regional_Plot

      RSID_Plots <- NULL

      if (nrow(df_rsid) > 0 && exists("Regional_Plot")) {

        message("Generating RSID Regional_Plot() outputs")

        utils::capture.output({

          RSID_Plots <- suppressMessages(suppressWarnings(

            Regional_Plot(

              Data = df_rsid,
              Gene_Tracks = TRUE,
              Separation_Distance = 1e6,
              Region_Window = 1e5,
              Verbose = TRUE

            )

          ))

        }, type = "output")

        if (is.list(RSID_Plots) && length(RSID_Plots)) {

          for (nm in names(RSID_Plots)) {

            dh <- attr(RSID_Plots[[nm]], "dynamic_height")

            if (!is.null(dh) && is.finite(dh) && dh > 0) {

              attr(RSID_Plots[[nm]], "recommended_height") <- dh

            }

            attr(RSID_Plots[[nm]], "recommended_width")  <- 30

            attr(RSID_Plots[[nm]], "suggested_filename") <- paste0(

              pheno, "_", gsub("[^A-Za-z0-9]+", "_", nm), ".jpg"

            )

          }

          message("RSID Regional_Plot outputs complete (", length(RSID_Plots), " plot(s))")

        } else {

          message("RSID Regional_Plot returned no plots")

        }

      } else {

        if (nrow(df_rsid) == 0) message("No RSIDs detected; skipping Regional_Plot()")

        if (nrow(df_rsid) > 0 && !exists("Regional_Plot")) message("Regional_Plot() not found; skipping RSID plots")

      }

      message("Packaging outputs")

      out <- list(

        pheno       = pheno,
        HLA_plot    = p_hla,
        AA_plot     = p_aa,
        SNPS_plot   = p_snps,
        RSID_plots  = RSID_Plots
      )

      class(out) <- c("SegregateHLAPlots", class(out))

      message("Segregate_HLA_Plot() complete")

      out

    }

    # Wrapper: multi-input + per-dataset args + progress behaviour

    .Segregate_HLA_Plot_original <- Segregate_HLA_Plot

    if (!exists("use_wrapper")) use_wrapper <- TRUE

    if (isTRUE(use_wrapper)) {

      Segregate_HLA_Plot <- local({

        orig <- .Segregate_HLA_Plot_original

        resolve_one <- function(x, default_name = "DATA") {

          if (is.data.frame(x)) return(list(item = x, name = default_name))

          if (is.character(x) && length(x) == 1L) {

            if (file.exists(x)) {

              nm <- tools::file_path_sans_ext(basename(x))

              return(list(item = x, name = nm))

            }

            if (exists(x, inherits = TRUE)) {

              return(list(item = get(x, inherits = TRUE), name = x))

            }

            stop("Character input is neither an existing file nor an object name: ", x, call. = FALSE)

          }

          stop("Unsupported input type. Use a data.frame, a file path, or an object name string.", call. = FALSE)

        }

        # Per-dataset arg selection

        pick_for_dataset <- function(val, i, n, nm) {

          if (inherits(val, "AsIs")) return(val)

          if (is.list(val)) {

            nms <- names(val)

            if (!is.null(nms) && nzchar(nm) && nm %in% nms) return(val[[nm]])

            if (length(val) == n) return(val[[i]])

            return(val)

          }

          if (is.atomic(val) && length(val) > 1L && length(val) == n) {

            return(val[i])

          }

          val
        }

        # Fully evaluated args list + per-dataset slicing

        build_full_args_for_dataset <- function(user_call, data_value, caller_env, i, n_items, ds_name) {

          fmls <- formals(orig)
          fml_names <- names(fmls)

          user_list <- as.list(user_call)[-1]

          user_list$session <- NULL

          eval_env <- new.env(parent = environment(orig))

          out <- vector("list", length(fml_names))
          names(out) <- fml_names

          for (nm in fml_names) {

            if (nm == "Data") {

              val <- data_value
              out[[nm]] <- val
              assign(nm, val, envir = eval_env)

              next

            }

            if (!is.null(user_list[[nm]])) {

              raw_val <- eval(user_list[[nm]], envir = caller_env)
              val <- pick_for_dataset(raw_val, i = i, n = n_items, nm = ds_name)

            } else {

              def <- fmls[[nm]]

              if (is.symbol(def) && identical(def, quote(expr = ))) next

              val <- eval(def, envir = eval_env)

            }

            out[[nm]] <- val
            assign(nm, val, envir = eval_env)

          }

          out

        }

        f <- function(Data, session = NULL) {

          cl <- match.call(expand.dots = TRUE)
          caller_env <- parent.frame()
          data_expr <- cl$Data

          build_items <- function() {

            # Data as bare symbol

            if (!is.null(data_expr) && is.name(data_expr)) {

              nm <- deparse(data_expr)

              if (exists(nm, inherits = TRUE)) {

                return(list(list(item = get(nm, inherits = TRUE), name = nm)))

              }

            }

            x <- eval(data_expr, envir = caller_env)

            if (is.list(x) && !is.data.frame(x)) {

              if (!length(x)) stop("Data list is empty.", call. = FALSE)

              out <- vector("list", length(x))

              for (j in seq_along(x)) {

                out[[j]] <- resolve_one(x[[j]], default_name = paste0("DATA_", j))

              }

              return(out)

            }

            if (is.character(x) && length(x) > 1L) {

              out <- vector("list", length(x))

              for (j in seq_along(x)) {

                out[[j]] <- resolve_one(x[[j]], default_name = paste0("DATA_", j))

              }

              return(out)

            }

            list(resolve_one(x, default_name = "DATA_1"))
          }

          items <- build_items()
          n_items <- length(items)

          run_one <- function(resolved, i) {

            item <- resolved$item
            nm   <- resolved$name

            # (Verbose or not)

            base::message(sprintf("Processing dataset: %s", nm))

            # Stamp name so core uses it for pheno when Data is a df

            if (is.data.frame(item)) attr(item, ".pheno") <- nm

            args <- build_full_args_for_dataset(
              user_call  = cl,
              data_value = item,
              caller_env = caller_env,
              i          = i,
              n_items    = n_items,
              ds_name    = nm
            )

            verbose_val <- isTRUE(args$Verbose)

            if (isTRUE(verbose_val)) {

              return(do.call(orig, args))

            } else {

              args$Verbose <- FALSE

              suppressWarnings(

                run_with_counter(

                  func    = orig,
                  args    = args,
                  session = session

                )

              )

            }

          }

          if (n_items == 1L) return(run_one(items[[1]], 1L))

          out <- vector("list", n_items)
          nms <- vapply(items, `[[`, character(1), "name")

          for (i in seq_along(items)) out[[i]] <- run_one(items[[i]], i)

          names(out) <- make.unique(nms)
          class(out) <- c("SegregateHLAPlotsList", class(out))
          out

        }

        # Sync wrapper signature to original + session

        formals(f) <- c(formals(orig), alist(session = NULL))

        f

      })

    }
