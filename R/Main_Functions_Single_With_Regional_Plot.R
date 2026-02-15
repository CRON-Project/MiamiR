
#' Customised Manhattan Plot with Inset Regional Plots
#'
#' @param Single_Plot Object returned from Single_Plot(); defaults to NULL
#' @param Regional_Plot Object(s) returned from Regional_Plot(); defaults to NULL
#' @param LD_Legend_Plot Specific Regional_Plot() object to derive LD_Legend from; defaults to NULL, leading to first object being utilised
#' @param Gene_Legend_Plot Specific Regional_Plot() object to derive Gene_Legend from; defaults to NULL, leading to merged legend being created across and from any objects supplied
#' @param LD_Legend_Scale Proportion of default sizing to use when displaying LD key on Single_Plot() background derived from inset Regional_Plot(); defaults to 1
#' @param Gene_Legend_Scale Proportion of default sizing to use when displaying Gene Biotype key on Single_Plot() background derived from inset Regional_Plot(); defaults to 1
#' @param LD_Legend_Height Proportion of default sizing height to use when anchoring LD key on Single_Plot() background; defaults to NULL
#' @param Gene_Legend_Height Proportion of default sizing height to use when anchoring Gene Biotype key on Single_Plot() background; defaults to NULL, leading to max anchoring relative to default or specified Gene_Legend_Preset
#' @param LD_Legend_Preset Automatically set anchoring of LD key on Single_Plot() background relative to common positioning; defaults to NULL, leading to "TopRight" being applied
#' @param Gene_Legend_Preset Automatically set anchoring of Gene Biotype key on Single_Plot() background relative to common positioning; defaults to NULL, leading to "TopRight" being applied
#' @param Target_Width Target save width in inches which the user intends on applying when using ggsave on the combined object; defaults to 30
#' @param Target_Height Target save height in inches which the user intends on applying when using ggsave on the combined object; defaults to 15
#' @param Inset_Scale Proportion of the full default size of the Regional_Plot() panels to use when overlaying as insets; defaults to 0.17
#' @param Inset_Chromosomes Allocated chromosomes where inset based on lead signal/set of SNPs is to be inserted; defaults to NULL leading to all Regional_Plot() objects provided to Regional_Plot argument being inserted in ascending order of significant signal on combined plot and by current order of regional objects
#' @param Preset_Tweak_X Add further spacing along the X-axis to the predefined Inset_Preset direction of Regional_Plot() anchoring around signal of interest; defaults to 0
#' @param Preset_Tweak_Y Add further spacing along the Y-axis to the predefined Inset_Preset direction of Regional_Plot() anchoring around signal of interest; defaults to 0
#' @param Preset_Regional_Height_Tweak Proportion of the full default Inset_Scale height of the Regional_Plot() panels to use when overlaying as insets allowing for further control in 2-Dimensions; defaults to 0.3
#' @param Inset_Preset Relative predefined position of Regional_Plot() inset relative to peak of interest; defaults to "Right"
#' @param Chromosome_Column Manually specify chromosome column; defaults to NULL, leading to auto-detection
#' @param Peak_Highlight_On Draw box around significant signals of interest from which inset expands from; defaults to TRUE
#' @param Peak_Highlight_Colour Colour of box around significant signals of interest from which inset expands from where Peak_Highlight_On is TRUE; defaults to "grey40"
#' @param Peak_Highlight_Linetype Line type of box around significant signals of interest from which inset expands from where Peak_Highlight_On is TRUE; defaults to "solid"
#' @param Peak_Highlight_Thickness Thickness of box around Regional_Plot() inset; defaults to 0.45
#' @param Regional_Inset_Highlight_On Draw box around Regional_Plot() inset; defaults to TRUE
#' @param Regional_Inset_Highlight_Colour Colour of box around Regional_Plot() inset; defaults to "grey40"
#' @param Regional_Inset_Highlight_Linetype Line type of box around Regional_Plot() inset; defaults to "solid"
#' @param Regional_Inset_Highlight_Thickness Line thickness of box around Regional_Plot() inset; defaults to 0.45
#' @param Inset_Lines_On Draw lines from box around significant signals of interest from which inset expands from to the box around the Regional_Plot() inset; defaults to TRUE
#' @param Inset_Lines_Colour Colour of lines from box around significant signals of interest from which inset expands from to the box around the Regional_Plot() inset; defaults to "black"
#' @param Inset_Lines_Linetype Type of lines from box around significant signals of interest from which inset expands from to the box around the Regional_Plot() inset; defaults to "dotted"
#' @param Inset_Lines_Thickness Thickness of lines from box around significant signals of interest from which inset expands from to the box around the Regional_Plot() inset; defaults to 0.45
#' @param Inset_Lines_Gap_Shading_On Shade space between lines from box around significant signals of interest from which inset expands from to the box around the Regional_Plot() inset; defaults to TRUE
#' @param Inset_Lines_Gap_Shading_Fill Colour of shaded space between lines from box around significant signals of interest from which inset expands from to the box around the Regional_Plot() inset; defaults to "grey70"
#' @param Inset_Lines_Gap_Shading_Alpha Opacity of shaded space between lines from box around significant signals of interest from which inset expands from to the box around the Regional_Plot() inset; defaults to 0.45
#' @param Peak_Highlight_To_Zero When Peak_Highlight_On is TRUE, ensure that the box drawn around significant signals of interest from which inset expands from reach to the baseline of background Single_Plot(); defaults to TRUE
#' @param Transparent_Regionals Ensure that Regional_Plot() insets are transparent, therefore allowing any overlapped elements on Single_Plot() background to be visible; defaults to FALSE
#' @param Preset_Regional_Width_Tweak Proportion of the full default Inset_Scale width of the Regional_Plot() panels to use when overlaying as insets allowing for further control in 2-Dimensions; defaults to 1
#' @param Verbose Prevent display of progress bar as function is running and instead show key milestone outputs/messages (mainly for debugging purposes); defaults to FALSE
#' @param Peak_Highlight_Tweak_Y_Up Add/restrict spacing (matching to scale) along the maximum position on the Y-axis to the predefined Peak Highlight box around signal(s) of interest; defaults to 0.15
#' @param Peak_Highlight_Tweak_Y_Down Add/restrict spacing (matching to scale) along the minimum position on the Y-axis to the predefined Peak Highlight box around signal(s) of interest; defaults to 2
#' @param Peak_Highlight_Tweak_X Add/restrict spacing (matching to scale) along the midpoint on the X-axis to the predefined Peak Highlight box around signal(s) of interest; defaults to 0.0035
#' @param Inset_Lines_To_Zero When Inset_Lines_On is TRUE, ensure that lines from box around significant signals of interest from which inset expands from are also anchored to the baseline of the Single_Plot() background when Peak_Highlight_To_Zero is TRUE; defaults to TRUE
#' @param Transparent_Legends Ensure that Regional_Plot() legend insets are transparent, therefore allowing any overlapped elements on Single_Plot() background to be visible; defaults to FALSE
#' @param LD_Legend_On Show Auto_LD Legend derived from Regional_Plot object(s), when Auto_LD was TRUE; defaults to TRUE
#' @param Gene_Legend_On Show Gene_Legend, derived from Regional_Plot object(s); defaults to TRUE
#' @param Legend_Order When both LD_Legend_On & Gene_Legend_On are TRUE, decide the vertical order of the legends when using default spacing options; defaults to c("LD","Gene")
#' @param Legend_Gap_Frac Fixed vertical gap between the two legend boxes (as fraction of panel height); defaults to 0.04
#' @param LD_Legend_Horizontal Proportion of default sizing horizontal positioning to use when anchoring LD  key on Single_Plot() background; defaults to NULL, leading to max anchoring relative to default or specified LD_Legend_Preset
#' @param Gene_Legend_Horizontal Proportion of default sizing horizontal positioning to use when anchoring Gene Biotype key on Single_Plot() background; defaults to NULL, leading to max anchoring relative to default or specified Gene_Legend_Preset
#' @param Inset_Lines_Lower_Frac Proportion of the default Peak_Highlight_On box to draw the lower Inset_Line from; defaults to 0, allowing for full inset spacing
#'
#' @returns Image of Manhattan Plot with inset Regional_Plot(s) is allocated to specified object and the resulting ggplot object can then be saved to an image
#'
#' @export
#'
#' @examples
#'
#'
#' Manhattan_Plot <- Single_Plot(Data = Household_Income_Sum_Stats, Label_Index = FALSE)
#'
#' Regional_Plot <- Regional_Plot(Data = c(Household_Income_Sum_Stats, Household_Income_Sum_Stats),
#'                                Chromosomes = c(9,3),
#'                                Separation_Distance = 1e9,
#'                                Region_Window = c(3e5, 8e4),
#'                                Auto_LD = TRUE,
#'                                Genome_Build = "grch37",
#'                                Title_On = FALSE,
#'                                X_Axis_Title_Size = 100,
#'                                X_Axis_Title_Vjust = 0,
#'                                Y_Axis_Title_Vjust = 5.5,
#'                                Y_Axis_Title_Size = 100,
#'                                Y_Axis_Text_Size = 90,
#'                                Y_Axis_Title = "\n-log₁₀(P)",
#'                                Chromosome_Label_Size = 90,
#'                                Sig_Line_Width = 5,
#'                                Recombination_Line = TRUE,
#'                                Recombination_Line_Thickness = 7.5,
#'                                Recombination_Axis_Text_Size = 90,
#'                                Recombination_Axis_Title_Size= 100,
#'                                Recombination_Axis_Title_Vjust = 5.5,
#'                                Recombination_Axis_Title     = "\nRR (cM/Mb)",
#'                                Point_Size = 20,
#'                                Label_Size = 25,
#'                                Label_Height = 40,
#'                                Label_Angle  = 2.5,
#'                                Panel_Ratio  = "1:1",
#'                                Gene_Label_Vertical_Spacing = 0.25,
#'                                Gene_Track_Vertical_Spacing =  0.55,
#'                                LD_Legend_Title_Size = 20,
#'                                LD_Legend_Text_Size = 15,
#'                                Gene_Legend_Title_Size = 20,
#'                                Gene_Legend_Text_Size = 15,
#'                                LD_Legend_On = FALSE,
#'                                Gene_Legend_On = FALSE,
#'                                Manual_Decimal_Places  = c(3,2),
#'                                Target_Position_Breaks = 3,
#'                                Gene_Label_Buffer = 6.5,
#'                                Gene_Label_Size = 20,
#'                                Gene_Structure_Size_Exons = 24,
#'                                Gene_Structure_Size_Introns = 5,
#'                                Sense_Arrow_Body_Length = 0.02,
#'                                Sense_Arrow_Gene_Gap = 0.015,
#'                                Sense_Arrow_Head_Size = 7)
#'
#'
#' Annotated_Man <- Single_With_Regional_Plot(
#'
#'   Single_Plot   = Manhattan_Plot,
#'   Regional_Plot = Regional_Plot,
#'   Inset_Chromosomes = c(3,9),
#'   Inset_Preset = c("left", "right"),
#'   Preset_Tweak_Y    = c(0.52, 0.4))
#'
#'
#'
#'

   Single_With_Regional_Plot <- function(   Single_Plot = NULL,
                                            Regional_Plot = NULL,
                                            LD_Legend_Plot    = NULL,
                                            Gene_Legend_Plot  = NULL,
                                            LD_Legend_Scale   = 1,
                                            Gene_Legend_Scale = 1,
                                            LD_Legend_Height  = NULL,
                                            Gene_Legend_Height= NULL,
                                            LD_Legend_Horizontal  = NULL,
                                            Gene_Legend_Horizontal =  NULL,
                                            LD_Legend_Preset  = NULL,
                                            Gene_Legend_Preset = NULL,
                                            Target_Width  = 30,
                                            Target_Height = 15,
                                            Inset_Scale = 0.17,
                                            Inset_Chromosomes = NULL,
                                            Preset_Tweak_X = 0,
                                            Preset_Tweak_Y = 0,
                                            Preset_Regional_Height_Tweak = 0.3,
                                            Preset_Regional_Width_Tweak  = 1,
                                            Inset_Preset = "Right",
                                            Chromosome_Column = NULL,
                                            Peak_Highlight_On = TRUE,
                                            Peak_Highlight_Colour = "grey40",
                                            Peak_Highlight_Linetype = "solid",
                                            Peak_Highlight_Thickness = 0.45,
                                            Regional_Inset_Highlight_On = TRUE,
                                            Regional_Inset_Highlight_Colour = "grey40",
                                            Regional_Inset_Highlight_Linetype = "solid",
                                            Regional_Inset_Highlight_Thickness = 0.45,
                                            Inset_Lines_On = TRUE,
                                            Inset_Lines_Colour = "black",
                                            Inset_Lines_Linetype = "dotted",
                                            Inset_Lines_Thickness = 0.45,
                                            Inset_Lines_Gap_Shading_On = TRUE,
                                            Inset_Lines_Gap_Shading_Fill = "grey70",
                                            Inset_Lines_Gap_Shading_Alpha = 0.45,
                                            Peak_Highlight_To_Zero = TRUE,
                                            Peak_Highlight_Tweak_Y_Up =  0.15,
                                            Peak_Highlight_Tweak_Y_Down = 2,
                                            Peak_Highlight_Tweak_X = 0.0035,
                                            Inset_Lines_To_Zero = TRUE,
                                            Transparent_Regionals = FALSE,
                                            Transparent_Legends = FALSE,
                                            LD_Legend_On  = TRUE,
                                            Gene_Legend_On = TRUE,
                                            Legend_Order = c("LD","Gene"),
                                            Legend_Gap_Frac = 0.04,
                                            Inset_Lines_Lower_Frac = 0,
                                            Verbose = FALSE )


    {

    man <- .as_single_plot(Single_Plot)

    message("Extracted Single_Plot as ggplot object successfully")

    df <- if (!is.null(man$data) && nrow(man$data)) man$data else {

      grabbed <- NULL

      for (ly in man$layers) if (is.data.frame(ly$data) && nrow(ly$data)) { grabbed <- ly$data; break }

      if (is.null(grabbed)) grabbed <- ggplot2::ggplot_build(man)$plot$data

      grabbed

    }

    chr_col <- detect_chromosome_column(df, Chromosome_Column)

    x_col <- "new_pos"

    y <- .pick_y(df)
    invisible(ggplot2::ggplot_build(man))
    ranges <- .get_panel_ranges(man); xr <- ranges$xr; yr <- ranges$yr
    dx <- diff(xr); dy <- diff(yr)

    y_labels_num <- .get_y_labels_num(man)

    attached <- attr(man, "expanded_labels", exact = TRUE)

    if (!is.null(attached)) {

      attached_num <- suppressWarnings(as.numeric(gsub("[^0-9\\.\\-]+", "", as.character(attached))))
      attached_num <- attached_num[is.finite(attached_num)]

      max_label <- if (length(attached_num)) max(attached_num, na.rm = TRUE) else max(y_labels_num, na.rm = TRUE)

    } else {

      max_label <- max(y_labels_num, na.rm = TRUE)

    }

    message("Gathering Regional_Plot(s)")

    regs_list <- .as_reg_list(Regional_Plot)

    ar_panel  <- Target_Width / Target_Height
    thr       <- -log10(5e-8)

    chr <- df[[chr_col]]
    x   <- as.numeric(df[[x_col]])

    if (!is.null(Inset_Chromosomes)) Inset_Chromosomes <- as.character(Inset_Chromosomes)

    p_out   <- man
    reg_idx <- 0L

    boxes_df        <- list()
    inset_boxes_df  <- list()
    connectors_df   <- list()
    shade_poly_df   <- list()
    inset_log       <- list()

    for (this_chr in unique(chr)) {

      idx_chr <- which(chr == this_chr & is.finite(x) & is.finite(y))

      if (!length(idx_chr)) next

      ok <- idx_chr[y[idx_chr] >= thr]

      if (!length(ok)) { ; next }

      if (!is.null(Inset_Chromosomes) && !(as.character(this_chr) %in% Inset_Chromosomes)) {

        next

      }

      xmin_raw <- min(x[ok], na.rm = TRUE)
      xmax_raw <- max(x[ok], na.rm = TRUE)
      ymin_raw <- thr
      ymax_raw <- max(y[ok], na.rm = TRUE)

      ymin_buf <- ymin_raw
      ymax_buf <- ymax_raw

      # Highlight box center & original half-width

      mid_x   <- (xmin_raw + xmax_raw) / 2
      half_w0 <- (xmax_raw - xmin_raw) / 2

      val <- as.numeric(Peak_Highlight_Tweak_X %||% 0)

      half_w <- half_w0 + val * dx

      half_w <- max(half_w, 0.001 * dx)
      half_w <- min(half_w, 0.50  * dx)

      xmin_buf <- mid_x - half_w
      xmax_buf <- mid_x + half_w

      xmin_buf <- max(xr[1], xmin_buf)
      xmax_buf <- min(xr[2], xmax_buf)

      if (xmin_buf > xmax_buf) { tmp <- xmin_buf; xmin_buf <- xmax_buf; xmax_buf <- tmp }

      next_ins_idx <- reg_idx + 1L
      ins_scale <- .per_inset_val(Inset_Scale,     next_ins_idx, this_chr, default = 0.10)
      ins_gap   <- .per_inset_val(0.1,  next_ins_idx, this_chr, default = 0.40)
      ins_yoff  <- .per_inset_val(0.00, next_ins_idx, this_chr, default = 0.00)
      preset_i  <- .per_inset_preset(Inset_Preset, next_ins_idx, this_chr, default = "Right")

      tweak_x  <- .per_inset_num(Preset_Tweak_X, next_ins_idx, this_chr, default = 0)
      tweak_y  <- .per_inset_num(Preset_Tweak_Y, next_ins_idx, this_chr, default = 0)

      inset_lines_lower_frac_i <- .per_inset_num(

        Inset_Lines_Lower_Frac,
        next_ins_idx,
        this_chr,
        default = 0

      )

      # clamp to [0,1] so it always behaves as a fraction

      if (!is.finite(inset_lines_lower_frac_i)) inset_lines_lower_frac_i <- 0
      inset_lines_lower_frac_i <- max(0, min(1, inset_lines_lower_frac_i))

      # per-peak highlight tweaks (match inset index / chromosome)

      ph_y_up <- .per_inset_num(

        Peak_Highlight_Tweak_Y_Up,
        next_ins_idx,
        this_chr,
        default = 0.15

      )

      reg_idx <- next_ins_idx
      reg_i   <- regs_list[[reg_idx]]
      reg_t <- reg_i

      dyn_h <- suppressWarnings(as.numeric(attr(reg_i, "dynamic_height", exact = TRUE)))[1]

      if (!is.finite(dyn_h) || dyn_h <= 0) dyn_h <- 20.1

      height_tweak <- .per_inset_val(
        Preset_Regional_Height_Tweak,
        next_ins_idx,
        this_chr,
        default = 1
      )

      width_tweak <- .per_inset_val(
        Preset_Regional_Width_Tweak,
        next_ins_idx,
        this_chr,
        default = 1
      )

      ins_idx  <- as.integer(next_ins_idx)
      chr_lbl  <- as.character(this_chr)

      message(sprintf("[Building Inset %d, CHR=%s] ", ins_idx, chr_lbl))

      if(height_tweak > 0)

      {

      scale_h <- (dyn_h / Target_Height) * height_tweak

      }

      else{

      scale_h <- (dyn_h / Target_Height)

      }

      saved_w <- Target_Width
      saved_h <- dyn_h * scale_h

      if (isTRUE(Transparent_Regionals)) {

        message("Making regional transparent")

        tr <- .transparentify_plot_to_grob(
          p         = reg_t,
          width_in  = saved_w,
          height_in = saved_h,
          dpi       = 100
        )

        g_inset <- tr$grob
        ar_img  <- tr$ar

      } else {

        tmp_png <- tempfile(fileext = ".png")
        ggplot2::ggsave(tmp_png, reg_t, width = saved_w, height = saved_h, limitsize = FALSE, dpi = 100)
        img <- magick::image_read(tmp_png); unlink(tmp_png)

        g_inset <- grid::rasterGrob(

          as.raster(img),

          x      = grid::unit(0,"npc"),
          y      = grid::unit(0,"npc"),
          width  = grid::unit(1,"npc"),
          height = grid::unit(1,"npc"),
          just   = c("left","bottom"),
          interpolate=TRUE

        )

        info   <- magick::image_info(img)
        ar_img <- info$width / info$height

      }

      height_npc  <- ins_scale * ar_panel / ar_img
      width_data  <- (ins_scale * dx) * width_tweak
      height_data <- height_npc  * dy

      top_big  <- .to_big_one(ymax_buf, max_label)
      bot_big  <- .to_big_one(ymin_buf, max_label)
      ymid_big <- 0.5 * (top_big + bot_big)
      xmid_main <- 0.5 * (xmin_buf + xmax_buf)

      pset <- tolower(preset_i)

      message("Positioning inset relative to peak")

      if (pset == "right") {

        xmin_ins <- xmax_buf + (ins_gap + tweak_x) * width_data
        xmax_ins <- xmin_ins + width_data
        ymax_ins_big <- ymid_big + (0.5 + ins_yoff + tweak_y) * height_data
        ymin_ins_big <- ymax_ins_big - height_data

      } else if (pset == "left") {

        xmax_ins <- xmin_buf - (ins_gap + tweak_x) * width_data
        xmin_ins <- xmax_ins - width_data
        ymax_ins_big <- ymid_big + (0.5 + ins_yoff + tweak_y) * height_data
        ymin_ins_big <- ymax_ins_big - height_data

      } else if (pset == "top") {

        xmid_ins   <- xmid_main + tweak_x * width_data
        xmin_ins   <- xmid_ins - 0.5 * width_data
        xmax_ins   <- xmid_ins + 0.5 * width_data
        ymin_ins_big <- top_big + (ins_gap + ins_yoff + tweak_y) * height_data
        ymax_ins_big <- ymin_ins_big + height_data

      } else if (pset == "bottom") {

        xmid_ins   <- xmid_main + tweak_x * width_data
        xmin_ins   <- xmid_ins - 0.5 * width_data
        xmax_ins   <- xmid_ins + 0.5 * width_data
        ymax_ins_big <- bot_big - (ins_gap - ins_yoff - tweak_y) * height_data
        ymin_ins_big <- ymax_ins_big - height_data

      }

      p_out <- p_out + ggplot2::annotation_custom(

        grob = g_inset,
        xmin = xmin_ins, xmax = xmax_ins,
        ymin = ymin_ins_big, ymax = ymax_ins_big

      )

      top_y_inset_raw <- .big_to_raw_with_label(ymax_ins_big, max_label)
      bot_y_inset_raw <- .big_to_raw_with_label(ymin_ins_big, max_label)

      if (Peak_Highlight_On) {

        message("Highlighting peak of interest")

        ymin_peak <- if (isTRUE(Peak_Highlight_To_Zero)) max(yr[1], 0) else ymin_buf

        ymin_peak <- ymin_peak + Peak_Highlight_Tweak_Y_Down
        ymax_buf  <- ymax_buf  + ph_y_up

        boxes_df[[length(boxes_df)+1]] <- data.frame(

          CHROM = this_chr, reg_index = reg_idx,
          xmin = xmin_buf, xmax = xmax_buf,
          ymin = ymin_peak, ymax = ymax_buf

        )

      }

      if (Regional_Inset_Highlight_On) {

        message("Highlighting Inset")

        inset_boxes_df[[length(inset_boxes_df)+1]] <- data.frame(

          CHROM=this_chr, reg_index=reg_idx,
          xmin=xmin_ins, xmax=xmax_ins, ymin=bot_y_inset_raw, ymax=top_y_inset_raw

        )

      }

      if (Inset_Lines_On || Inset_Lines_Gap_Shading_On) {

        message("Drawing inset connectors")

        if (pset %in% c("right","left")) {

          # decide lower anchor for connectors inside the peak box

          # start at bottom of peak box

          anchor_low <- ymin_peak

          if (isTRUE(Inset_Lines_To_Zero)) {

            # force to baseline of plot when requested

            anchor_low <- max(yr[1], 0)

          } else {

            # per-peak fraction up from the bottom of the peak box

            f <- inset_lines_lower_frac_i

            anchor_low <- ymin_peak + f * (ymax_buf - ymin_peak)

          }

          if (pset == "right") {

            # RIGHT-SIDE

            connectors_df[[length(connectors_df)+1]] <- data.frame(

              CHROM = this_chr, reg_index = reg_idx,
              x    = c(xmax_buf, xmax_buf),
              y    = c(ymax_buf, anchor_low),
              xend = c(xmin_ins, xmin_ins),
              yend = c(top_y_inset_raw, bot_y_inset_raw)

            )

            if (Inset_Lines_Gap_Shading_On) {

              message("Shading inset connection geom (right)")

              shade_poly_df[[length(shade_poly_df)+1]] <- data.frame(

                CHROM = this_chr, reg_index = reg_idx,
                x = c(xmax_buf, xmin_ins, xmin_ins, xmax_buf),
                y = c(ymax_buf, top_y_inset_raw, bot_y_inset_raw, anchor_low),
                gid = paste0("shade_", this_chr, "_", reg_idx)

              )

            }

          } else if (pset == "left") {

            # LEFT-SIDE BEHAVIOUR

            connectors_df[[length(connectors_df)+1]] <- data.frame(

              CHROM = this_chr, reg_index = reg_idx,
              x    = c(xmin_buf, xmin_buf),
              y    = c(ymax_buf, anchor_low),
              xend = c(xmax_ins, xmax_ins),
              yend = c(top_y_inset_raw, bot_y_inset_raw)

            )

            if (Inset_Lines_Gap_Shading_On) {

              message("Shading inset connection geom (left)")

              shade_poly_df[[length(shade_poly_df)+1]] <- data.frame(

                CHROM = this_chr, reg_index = reg_idx,
                x = c(xmin_buf, xmax_ins, xmax_ins, xmin_buf),
                y = c(ymax_buf, top_y_inset_raw, bot_y_inset_raw, anchor_low),
                gid = paste0("shade_", this_chr, "_", reg_idx)

              )

            }

          }

        }

          else if (pset == "top") {

          connectors_df[[length(connectors_df)+1]] <- data.frame(

            CHROM=this_chr, reg_index=reg_idx,
            x    = c(xmin_buf, xmax_buf),
            y    = c(ymax_buf, ymax_buf),
            xend = c(xmin_ins, xmax_ins),
            yend = c(bot_y_inset_raw, bot_y_inset_raw)

          )

          if (Inset_Lines_Gap_Shading_On) {

            shade_poly_df[[length(shade_poly_df)+1]] <- data.frame(

              CHROM=this_chr, reg_index=reg_idx,
              x = c(xmin_buf, xmin_ins, xmax_ins, xmax_buf),
              y = c(ymax_buf, bot_y_inset_raw, bot_y_inset_raw, ymax_buf),
              gid = paste0("shade_", this_chr, "_", reg_idx)

            )

          }

        } else if (pset == "bottom") {

          connectors_df[[length(connectors_df)+1]] <- data.frame(

            CHROM=this_chr, reg_index=reg_idx,
            x    = c(xmin_buf, xmax_buf),
            y    = c(ymin_buf, ymin_buf),
            xend = c(xmin_ins, xmax_ins),
            yend = c(top_y_inset_raw, top_y_inset_raw)

          )

          if (Inset_Lines_Gap_Shading_On) {

            shade_poly_df[[length(shade_poly_df)+1]] <- data.frame(

              CHROM=this_chr, reg_index=reg_idx,
              x = c(xmin_buf, xmin_ins, xmax_ins, xmax_buf),
              y = c(ymin_buf, top_y_inset_raw, top_y_inset_raw, ymin_buf),
              gid = paste0("shade_", this_chr, "_", reg_idx)

            )

          }

        }

      }

      inset_log[[length(inset_log)+1]] <- data.frame(

        CHROM=this_chr, reg_index=reg_idx,
        inset_xmin=xmin_ins, inset_xmax=xmax_ins,
        inset_ymin_raw=bot_y_inset_raw, inset_ymax_raw=top_y_inset_raw,
        Inset_Scale_Used=ins_scale, Inset_X_Spacing_Used=ins_gap, Inset_Y_Spacing_Used=ins_yoff,
        Inset_Preset=preset_i, Preset_Tweak_X=tweak_x, Preset_Tweak_Y=tweak_y,
        stringsAsFactors=FALSE

      )

    }

    boxes_df        <- if (length(boxes_df)) do.call(rbind, boxes_df) else data.frame()

    inset_boxes_df  <- if (length(inset_boxes_df)) do.call(rbind, inset_boxes_df) else data.frame()

    connectors_df   <- if (length(connectors_df)) do.call(rbind, connectors_df) else data.frame()

    inset_log       <- if (length(inset_log)) do.call(rbind, inset_log) else data.frame()

    shade_poly_df   <- if (length(shade_poly_df)) do.call(rbind, shade_poly_df) else data.frame()

    if (Inset_Lines_Gap_Shading_On && nrow(shade_poly_df)) {

      p_out <- p_out + ggplot2::geom_polygon(

        data = shade_poly_df,
        ggplot2::aes(x = x, y = y, group = gid),
        fill = Inset_Lines_Gap_Shading_Fill,
        alpha = Inset_Lines_Gap_Shading_Alpha,
        colour = NA, inherit.aes = FALSE

      )

    }

    if (Peak_Highlight_On && nrow(boxes_df)) {

      p_out <- p_out + ggplot2::geom_rect(

        data = boxes_df,
        ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
        fill = NA,
        colour = Peak_Highlight_Colour,
        linetype = Peak_Highlight_Linetype,
        linewidth = Peak_Highlight_Thickness,
        inherit.aes = FALSE

      )

    }

    if (Regional_Inset_Highlight_On && nrow(inset_boxes_df)) {

      p_out <- p_out + ggplot2::geom_rect(

        data = inset_boxes_df,
        ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
        fill = NA,
        colour = Regional_Inset_Highlight_Colour,
        linetype = Regional_Inset_Highlight_Linetype,
        linewidth = Regional_Inset_Highlight_Thickness,
        inherit.aes = FALSE

      )

    }

    if (Inset_Lines_On && nrow(connectors_df)) {

      p_out <- p_out + ggplot2::geom_segment(

        data = connectors_df,
        ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
        colour = Inset_Lines_Colour,
        linetype = Inset_Lines_Linetype,
        linewidth = Inset_Lines_Thickness,
        lineend = "round",
        inherit.aes = FALSE

      )

    }

    gene_info_combined <- NULL

    {

      gene_specs <- lapply(regs_list, function(g) attr(g, "Gene_legend_spec", exact = TRUE))

      gene_specs <- Filter(Negate(is.null), gene_specs)

      if (length(gene_specs) > 0) {

        all_levels <- unique(unlist(lapply(gene_specs, function(gs) gs$levels)))

        col_map   <- list()
        shape_map <- list()

        for (gs in gene_specs) {

          levs <- gs$levels
          cols <- gs$legend_colors %||% gs$colors %||% gs$colours
          shps <- gs$shapes

          if (!is.null(cols)) {

            for (lv in levs) {

              if (!is.null(cols[[lv]]) && is.null(col_map[[lv]])) {

                col_map[[lv]] <- cols[[lv]]

              }

            }

          }

          if (!is.null(shps)) {

            for (lv in levs) {

              if (!is.null(shps[[lv]]) && is.null(shape_map[[lv]])) {

                shape_map[[lv]] <- shps[[lv]]

              }

            }

          }

        }

        legend_colors <- unlist(col_map[all_levels])

        if (length(legend_colors) == 0) legend_colors <- NULL

        names(legend_colors) <- all_levels

        shapes <- unlist(shape_map[all_levels])

        if (length(shapes) == 0) shapes <- NULL

        names(shapes) <- all_levels

        first <- gene_specs[[1]]

        gene_info_combined <- list(

          levels        = all_levels,
          legend_colors = legend_colors,
          shapes        = shapes,
          title         = first$title %||% "Gene Biotype",
          title_size    = first$title_size %||% 30,
          text_size     = first$text_size %||% 20,
          legend_margin = first$legend_margin %||% ggplot2::margin(0, 0, 0, 0),
          box_margin    = first$box_margin    %||% ggplot2::margin(0, 0, 0, 0)

        )

      }

    }

    regs_list <- lapply(.as_reg_list(Regional_Plot), .unwrap_plot)

    if (!is.null(LD_Legend_Plot) || !is.null(Gene_Legend_Plot) || length(regs_list) > 0) {

    message("Merging and inserting legend...")

    # Default LD legend source (first regional if not given)

    ld_plot <- if (!is.null(LD_Legend_Plot)) .as_single_plot(LD_Legend_Plot)

    else if (length(regs_list) > 0) .unwrap_plot(regs_list[[1]]) else NULL

    ld_info <- if (!is.null(ld_plot)) .get_attr(ld_plot, "LD_legend_spec") else NULL

    # Merge Gene legend across all regionals (and explicit Gene_Legend_Plot if supplied)

    spec_entries <- .collect_gene_specs(regs_list, Gene_Legend_Plot = Gene_Legend_Plot, include_explicit = TRUE)
    gene_merge   <- .merge_gene_specs(spec_entries)
    gene_info    <- gene_merge$spec

    # Presets

    side_from_preset <- function(preset_raw) {

      p <- tolower(preset_raw)

      if (p %in% c("topright","tr","right","bottomright","br")) "right"

      else if (p %in% c("topleft","tl","left","bottomleft","bl")) "left"

      else "right"

    }

    vert_from_preset <- function(preset_raw) {

      p <- tolower(preset_raw)

      if (p %in% c("bottomright","br","bottom","bottomleft","bl")) "bottom" else "top"

    }

    ld_preset_raw   <- LD_Legend_Preset   %||% "TopRight"

    gene_preset_raw <- Gene_Legend_Preset %||% ld_preset_raw

    side_ld   <- side_from_preset(ld_preset_raw)
    side_gene <- side_from_preset(gene_preset_raw)
    vert_ld   <- vert_from_preset(ld_preset_raw)
    vert_gene <- vert_from_preset(gene_preset_raw)

    # Use one shared base height for both legends

    legend_base_frac <- 1

    ld_g <- .legend_to_grob_fixed(

      make_plot_fn     = make_ld_legend_plot_from_spec,
      spec             = ld_info,
      scale_val        = LD_Legend_Scale,
      xr               = xr, yr = yr,
      Target_Width     = Target_Width,
      Target_Height    = Target_Height,
      base_height_frac = legend_base_frac,
      side             = side_ld,
      vert             = vert_ld,
      Transparent_Legends = Transparent_Legends

    )

    gene_g <- .legend_to_grob_fixed(

      make_plot_fn     = make_gene_legend_plot_from_spec,
      spec             = gene_info,
      scale_val        = Gene_Legend_Scale %||% LD_Legend_Scale,
      xr               = xr, yr = yr,
      Target_Width     = Target_Width,
      Target_Height    = Target_Height,
      base_height_frac = legend_base_frac,
      side             = side_gene,
      vert             = vert_gene,
      Transparent_Legends = Transparent_Legends

    )

      have_ld   <- !is.null(ld_g)   && isTRUE(LD_Legend_On)
      have_gene <- !is.null(gene_g) && isTRUE(Gene_Legend_On)

      if (!have_ld && !have_gene) {

        message("Legend specs present but failed to render \u2013 skipping.")

      } else {

        left_x  <- xr[1]
        right_x <- xr[2]
        dy_pan  <- diff(yr)

        # Optional absolute tops (fractions)

        top_frac_ld   <- if (!is.null(LD_Legend_Height))   as.numeric(LD_Legend_Height[1])   else NA_real_

        top_frac_gene <- if (!is.null(Gene_Legend_Height)) as.numeric(Gene_Legend_Height[1]) else NA_real_


        xfrac_ld   <- if (!is.null(LD_Legend_Horizontal))   as.numeric(LD_Legend_Horizontal[1])   else NA_real_

        xfrac_gene <- if (!is.null(Gene_Legend_Horizontal)) as.numeric(Gene_Legend_Horizontal[1]) else NA_real_

        if (is.finite(xfrac_ld))   xfrac_ld   <- max(0, min(1, xfrac_ld))

        if (is.finite(xfrac_gene)) xfrac_gene <- max(0, min(1, xfrac_gene))

        if (is.finite(top_frac_ld))   top_frac_ld   <- max(0, min(1, top_frac_ld))

        if (is.finite(top_frac_gene)) top_frac_gene <- max(0, min(1, top_frac_gene))

        # Utilities

        x_range_for <- function(side, wdat) if (side == "right") c(right_x - wdat, right_x) else c(left_x, left_x + wdat)

        clamp_y <- function(ymin, ymax) {

          if ((ymax - ymin) > dy_pan) { ymin <- yr[1]; ymax <- yr[2] }

          if (ymin < yr[1]) { ymax <- min(yr[2], ymax + (yr[1] - ymin)); ymin <- yr[1] }

          if (ymax > yr[2]) { ymin <- max(yr[1], yr[2] - (ymax - ymin)); ymax <- yr[2] }

          c(ymin, ymax)

        }

        clamp_x <- function(xmin, xmax) {

          if ((xmax - xmin) > dx) { xmin <- xr[1]; xmax <- xr[2] }

          if (xmin < xr[1]) { xmax <- min(xr[2], xmax + (xr[1] - xmin)); xmin <- xr[1] }

          if (xmax > xr[2]) { xmin <- max(xr[1], xr[2] - (xmax - xmin)); xmax <- xr[2] }

          c(xmin, xmax)
        }


        add_box <- function(gr, side, vert, top_frac_opt, wdat, hdat, xfrac_opt = NA_real_) {

          # Y placement

          if (is.finite(top_frac_opt)) {

            y_max <- yr[1] + top_frac_opt * dy_pan
            y_min <- y_max - hdat

          } else if (vert == "top") {

            y_max <- yr[2]; y_min <- y_max - hdat

          } else {

            y_min <- yr[1]; y_max <- y_min + hdat

          }

          ybox <- clamp_y(y_min, y_max)

          # X placement

          if (is.finite(xfrac_opt)) {

            x_center <- xr[1] + xfrac_opt * dx
            x_min    <- x_center - wdat / 2
            x_max    <- x_center + wdat / 2

          } else if (side == "right") {

            x_min <- xr[2] - wdat; x_max <- xr[2]

          } else {

            x_min <- xr[1];        x_max <- xr[1] + wdat

          }

          xbox <- clamp_x(x_min, x_max)

          list(xmin = xbox[1], xmax = xbox[2], ymin = ybox[1], ymax = ybox[2], grob = gr)
        }

        # Stack legends with fixed gap

        # If forced horizontal placement for either legend, skip stacking and place independently

        force_independent_x <- is.finite(xfrac_ld) || is.finite(xfrac_gene)

        if (!force_independent_x &&
            have_ld && have_gene &&
            side_ld == side_gene &&
            vert_ld == vert_gene &&
            !is.finite(top_frac_ld) && !is.finite(top_frac_gene)) {

          message("Stacking and spacing legends")

          # Parse Legend_Order

          order_vec <- Legend_Order

          if (length(order_vec) == 1L && grepl(",", order_vec)) order_vec <- trimws(strsplit(order_vec, ",")[[1]])

          order_vec <- toupper(order_vec); order_vec <- order_vec[order_vec %in% c("LD","GENE")]

          if (!length(order_vec)) order_vec <- c("LD","GENE")

          if (length(order_vec) == 1L) order_vec <- c(order_vec, setdiff(c("LD","GENE"), order_vec))

          top_name    <- order_vec[1]

          second_name <- if (top_name == "LD") "GENE" else "LD"

          # Raw aligned boxes for X

          gene_box_raw <- add_box(gene_g$grob, side_gene, vert_gene, NA_real_, gene_g$wdat, gene_g$hdat)
          ld_box_raw   <- add_box(ld_g$grob,   side_ld,   vert_ld,   NA_real_, ld_g$wdat,   ld_g$hdat)

          if (side_ld == "left") {

            gene_box_raw$xmin <- left_x;  gene_box_raw$xmax <- left_x  + gene_g$wdat
            ld_box_raw$xmin   <- left_x;  ld_box_raw$xmax   <- left_x  + ld_g$wdat

          } else {

            gene_box_raw$xmax <- right_x; gene_box_raw$xmin <- right_x - gene_g$wdat
            ld_box_raw$xmax   <- right_x; ld_box_raw$xmin   <- right_x - ld_g$wdat

          }

          # Fixed gap in data units (order-invaria...)

          gap <- as.numeric(Legend_Gap_Frac %||% 0.04) * dy_pan

          # Build TOP anchored to panel edge

          if (top_name == "LD") {

            top_xmin <- ld_box_raw$xmin; top_xmax <- ld_box_raw$xmax; top_h <- ld_g$hdat;  top_grob <- ld_g$grob

          } else {

            top_xmin <- gene_box_raw$xmin; top_xmax <- gene_box_raw$xmax; top_h <- gene_g$hdat; top_grob <- gene_g$grob

          }

          if (vert_ld == "top") {

            top_box <- list(xmin = top_xmin, xmax = top_xmax, ymax = yr[2], ymin = yr[2] - top_h)

          } else {

            top_box <- list(xmin = top_xmin, xmax = top_xmax, ymin = yr[1], ymax = yr[1] + top_h)

          }

          ct <- clamp_y(top_box$ymin, top_box$ymax); top_box$ymin <- ct[1]; top_box$ymax <- ct[2]

          # SECOND strictly relative to TOP

          if (second_name == "LD") {

            second_xmin <- ld_box_raw$xmin; second_xmax <- ld_box_raw$xmax; second_h <- ld_g$hdat;  second_grob <- ld_g$grob

          } else {

            second_xmin <- gene_box_raw$xmin; second_xmax <- gene_box_raw$xmax; second_h <- gene_g$hdat; second_grob <- gene_g$grob

          }

          if (vert_ld == "top") {

            second_box <- list(xmin = second_xmin, xmax = second_xmax,
                               ymax = top_box$ymin - gap,
                               ymin = top_box$ymin - gap - second_h)

          } else {

            second_box <- list(xmin = second_xmin, xmax = second_xmax,
                               ymin = top_box$ymax + gap,
                               ymax = top_box$ymax + gap + second_h)

          }

          cs <- clamp_y(second_box$ymin, second_box$ymax); second_box$ymin <- cs[1]; second_box$ymax <- cs[2]

          # Draw (top then second)

          p_out <- p_out +
            ggplot2::annotation_custom(top_grob,    xmin = top_box$xmin,    xmax = top_box$xmax,
                                       ymin = top_box$ymin,    ymax = top_box$ymax) +
            ggplot2::annotation_custom(second_grob, xmin = second_box$xmin, xmax = second_box$xmax,
                                       ymin = second_box$ymin, ymax = second_box$ymax)

          message(sprintf("[Legends] Stacked on %s-%s; fixed gap: %.3f (data units).", side_ld, vert_ld, gap))

        } else {

          # Different sides or explicit heights = independent placement

          if (have_gene) {

            gb <- add_box(

              gene_g$grob, side_gene, vert_gene, top_frac_gene,
              gene_g$wdat, gene_g$hdat, xfrac_opt = xfrac_gene

            )

            p_out <- p_out + ggplot2::annotation_custom(gb$grob, xmin = gb$xmin,
                                                        xmax = gb$xmax, ymin = gb$ymin, ymax = gb$ymax)

          }

          if (have_ld) {

            lb <- add_box(

              ld_g$grob, side_ld, vert_ld, top_frac_ld,
              ld_g$wdat, ld_g$hdat, xfrac_opt = xfrac_ld

            )

            p_out <- p_out + ggplot2::annotation_custom(lb$grob, xmin = lb$xmin, xmax = lb$xmax,
                                                        ymin = lb$ymin, ymax = lb$ymax)

          }

        }

      }

    } else {

      message("No LD/Gene legend info \u2013 no legend insets added.")

    }

    message("Finalising and returning plot object")

          plot  <-  p_out

      }

      .Single_With_Regional_Plot_original <- Single_With_Regional_Plot

      Single_With_Regional_Plot <- function(..., Verbose = FALSE, session = NULL) {

        dots <- list(...)

        fn_formals <- formals(.Single_With_Regional_Plot_original)
        valid_args <- names(fn_formals)

        if ("verbose" %in% names(dots)) dots$verbose <- NULL

        clean_args <- dots[names(dots) %in% valid_args]

        for (arg in valid_args) {

          if (!arg %in% names(clean_args)) {

            default_expr <- fn_formals[[arg]]

            clean_args[[arg]] <- tryCatch(

              eval(default_expr, envir = environment(.Single_With_Regional_Plot_original)),

              error = function(e) eval(default_expr)

            )

          }

        }

        if (isTRUE(Verbose)) {

          return(invisible(.silence_warnings(do.call(.Single_With_Regional_Plot_original, clean_args))))

        }

        .silence_messages(

          .silence_warnings(

            run_with_counter(.Single_With_Regional_Plot_original, args = clean_args, session = session)

          )

        )

      }
