
#' Annotated, Customised Miami Plots
#'
#' @param Top_Data These are the summary statistics (either file path(s),environmental object(s) or a variable(s) containing such contents) to be plotted on the top; defaults to NULL
#' @param Bottom_Data These are the summary statistics (either file path(s),environmental object(s) or a variable(s) containing such contents) to be plotted on the bottom; defaults to NULL
#' @param Chromosome_Column Manually specify chromosome column; defaults to NULL, leading to auto-detection
#' @param Position_Column Manually specify chromosome genomic position column; defaults to NULL, leading to auto-detection
#' @param SNP_ID_Column Manually specify SNP ID column; defaults to NULL, leading to auto-detection
#' @param PValue_Column Manually specify P Value column; defaults to NULL, leading to auto-detection
#' @param Reference_Allele_Column Manually specify reference allele column; defaults to NULL, leading to auto-detection
#' @param Effect_Allele_Column Manually specify effect allele column; defaults to NULL, leading to auto-detection
#' @param Verbose Prevent display of progress bar as function is running and instead show key milestone outputs/messages (mainly for debugging purposes); defaults to FALSE.
#' @param Adjust Shift shared X_Axis_Text down relative to the top plot, the larger the value the greater the shift; defaults to NULL and if not passed then calculated from expression mapped to Chromosome_Label_Size inherted from Single_Plot()
#' @param ... Another shadow argument to detect inherited arguments passed and modified from the Single_Plot() function; defaults to list()
#' @param .dots Shadow argument to detect inherited arguments passed and modified from the Single_Plot() function; defaults to list()
#' @param session Blank session definitions; defaults to NULL
#'
#' @returns Image of Miami Plot(s) is allocated to specified object and the resulting combined grob object can then be saved to an image
#' @export
#'
#' @examples Miami_Plot <- Miami_Plot(Top_Data = Household_Income_Sum_Stats, Bottom_Data = Intelligence_Sum_Stats)

  Miami_Plot <- function(Top_Data = NULL,
                       Bottom_Data = NULL,
                       Chromosome_Column = NULL,
                       Position_Column = NULL,
                       SNP_ID_Column = NULL,
                       PValue_Column = NULL,
                       Reference_Allele_Column = NULL,
                       Effect_Allele_Column = NULL,
                       Verbose = FALSE,
                       Adjust = NULL,
                       ...,
                       session = NULL,
                       .dots = list()) {

  user_args <- .dots

  message("Separating Top and Bottom inputs")

  # Separate Top_ and Bottom_ arguments
  top_args <- user_args[grepl("^Top_", names(user_args))]
  bottom_args <- user_args[grepl("^Bottom_", names(user_args))]

  # Remove prefix
  names(top_args) <- sub("^Top_", "", names(top_args))
  names(bottom_args) <- sub("^Bottom_", "", names(bottom_args))

  # Shared args apply to both
  shared_args <- user_args[!grepl("^(Top_|Bottom_)", names(user_args))]

  # Merge shared into each
  top_args <- modifyList(shared_args, top_args)
  bottom_args <- modifyList(shared_args, bottom_args)

  # Defaults from Single_Plot
  defaults <- formals(.Single_Plot_original)
  eval_defaults <- lapply(defaults, function(x) {
    if (is.symbol(x) || is.language(x)) {
      tryCatch(eval(x), error = function(e) NULL)
    } else {
      x
    }
  })

  final_top_args <- modifyList(eval_defaults, top_args)
  final_bottom_args <- modifyList(eval_defaults, bottom_args)

  if (is.null(Adjust)) {

    message("Calculating adjust parameter for later")

    base_size <- 35
    base_adjust <- 0.66
    current_size <- final_top_args$Chromosome_Label_Size

    if (is.null(current_size)) {
      current_size <- base_size  # fallback
    }

    # Calculate adjustment: as label size decreases, Adjust increases
    Adjust <- base_adjust + (base_size - current_size) * 0.0045
  }




  if(is.null(bottom_args$Bottom_Anchor_Label))
  {
  final_bottom_args$Anchor_Label <- "left mirror"
  }

  if (!("Bottom_Anchor_Label" %in% names(user_args)) &&
      !("Anchor_Label" %in% names(user_args))) {

    message("Automatically reversing label view settings for bottom plot")

    final_bottom_args$Anchor_Label <- "left mirror"

  }

  if (is.character(Top_Data) && length(Top_Data) == 1 && file.exists(Top_Data)) {

    message("Loading Top Data from file path")

    Top_Data <- vroom::vroom(Top_Data, show_col_types = FALSE)

  }

  if (is.character(Bottom_Data) && length(Bottom_Data) == 1 && file.exists(Bottom_Data)) {

    message("Loading Bottom Data from file path")

    Bottom_Data <- vroom::vroom(Bottom_Data, show_col_types = FALSE)

  }

  message("Deducing key column names for Top Data")

  Chromosome_Column <- detect_chromosome_column(Top_Data, Chromosome_Column)

  Top_Data <- Top_Data %>%
    dplyr::mutate(!!Chromosome_Column := ifelse(.data[[Chromosome_Column]] == "23", "X", .data[[Chromosome_Column]]))

  PValue_Column     <- detect_pvalue_column(Top_Data, PValue_Column)
  Position_Column   <- detect_position_column(Top_Data, Position_Column)
  SNP_ID_Column     <- detect_snp_column(Top_Data, SNP_ID_Column)
  Ref_Allele_Column     <- detect_reference_allele_column(Top_Data, Reference_Allele_Column)
  Alt_Allele_Column     <- detect_effect_allele_column(Top_Data, Effect_Allele_Column)

  message("Storing before examining bottom dataset")

  Top_Chromosome_Column <- Chromosome_Column
  Top_PValue_Column     <- PValue_Column
  Top_Position_Column   <- Position_Column
  Top_SNP_ID_Column     <- SNP_ID_Column
  Top_Ref_Allele_Column <- Ref_Allele_Column
  Top_Alt_Allele_Column <- Alt_Allele_Column

  message("Assigning key columns for Top Data to standardised nomenclature")

  Top_Data$CHROM <- Top_Data[[Chromosome_Column]]
  Top_Data$GENPOS <- Top_Data[[Position_Column]]
  Top_Data$ID <- Top_Data[[SNP_ID_Column]]
  Top_Data$P <- Top_Data[[PValue_Column]]
  Top_Data$ALLELE0 <- Top_Data[[Ref_Allele_Column]]
  Top_Data$ALLELE1 <- Top_Data[[Alt_Allele_Column]]

  message("Resetting auto detect for bottom dataset")

  Chromosome_Column <- NULL
  PValue_Column <- NULL
  Position_Column <- NULL
  SNP_ID_Column <- NULL
  Ref_Allele_Column <- NULL
  Alt_Allele_Column<- NULL

  message("Deducing key column names for Bottom Data")

  Chromosome_Column <- detect_chromosome_column(Bottom_Data, Chromosome_Column)

  Bottom_Data <- Bottom_Data %>%
    dplyr::mutate(!!Chromosome_Column := ifelse(.data[[Chromosome_Column]] == "23", "X", .data[[Chromosome_Column]]))

  PValue_Column     <- detect_pvalue_column(Bottom_Data, PValue_Column)
  Position_Column   <- detect_position_column(Bottom_Data, Position_Column)
  SNP_ID_Column     <- detect_snp_column(Bottom_Data, SNP_ID_Column)
  Ref_Allele_Column     <- detect_reference_allele_column(Bottom_Data, Reference_Allele_Column)
  Alt_Allele_Column     <- detect_effect_allele_column(Bottom_Data, Effect_Allele_Column)

  Bottom_Chromosome_Column <- Chromosome_Column
  Bottom_PValue_Column     <- PValue_Column
  Bottom_Position_Column   <- Position_Column
  Bottom_SNP_ID_Column     <- SNP_ID_Column
  Bottom_Ref_Allele_Column <- Ref_Allele_Column
  Bottom_Alt_Allele_Column <- Alt_Allele_Column

  message("Assigning key columns for Top Data to standardised nomenclature")

  Bottom_Data$CHROM <- Bottom_Data[[Chromosome_Column]]
  Bottom_Data$GENPOS <- Bottom_Data[[Position_Column]]
  Bottom_Data$ID <- Bottom_Data[[SNP_ID_Column]]
  Bottom_Data$P <- Bottom_Data[[PValue_Column]]
  Bottom_Data$ALLELE0 <- Bottom_Data[[Ref_Allele_Column]]
  Bottom_Data$ALLELE1 <- Bottom_Data[[Alt_Allele_Column]]

  message("Create CHROM, GENPOS ID agnostic location keys for both data sets")

  Top_Data <- Top_Data %>%
    dplyr::mutate(combined_key = paste(GENPOS, CHROM, sep = "_"))

  Bottom_Data <- Bottom_Data %>%
    dplyr::mutate(combined_key = paste(GENPOS, CHROM, sep = "_"))

  message("Find rows in Bottom_Data that are not in Top_Data")

  unique_in_bottom <- Bottom_Data %>%
    dplyr::filter(!(combined_key %in% Top_Data$combined_key))

  message("Find rows in Top_Data that are not in Bottom_Data")

  unique_in_top <- Top_Data %>%
    dplyr::filter(!(combined_key %in% Bottom_Data$combined_key))

  if(nrow(unique_in_bottom) != 0)

  {

    message("Determining and adding buffer to Top_Data for missing Bottom_Data locations")

    new_rows <- unique_in_bottom %>%
      dplyr::select(GENPOS, CHROM)

    empty_cols <- Top_Data %>%
      dplyr::select(-GENPOS, -CHROM) %>%
      dplyr::summarise_all(~NA) %>%
      dplyr::slice(rep(1, nrow(new_rows)))

    new_data <-  dplyr::bind_cols(new_rows, empty_cols)
    new_data$FAKE <- 1
    new_data$P <- 1
    new_data[[Top_Position_Column]] <- new_data$GENPOS
    new_data[[Top_Chromosome_Column]] <- new_data$CHROM
    new_data[[Top_PValue_Column]] <- new_data$P

    Top_Data[Top_Chromosome_Column] <- as.character(Top_Data[Top_Chromosome_Column])
    new_data[Top_Chromosome_Column] <- as.character(new_data[Top_Chromosome_Column])
    Top_Data$CHROM <- as.character(Top_Data$CHROM)
    new_data$CHROM <- as.character(new_data$CHROM)

    Top_Data <-  dplyr::bind_rows(Top_Data, new_data)

  }

  Top_Data[Top_Chromosome_Column] <- NULL #weird character too much space

  if(nrow(unique_in_top) != 0)

  {

    message("Determining and adding buffer to Bottom_Data for missing Top_Data locations")

    new_rows <- unique_in_top %>%
      dplyr::select(GENPOS, CHROM)

    empty_cols <- Bottom_Data %>%
      dplyr::select(-GENPOS, -CHROM) %>%
      dplyr::summarise_all(~NA) %>%
      dplyr::slice(rep(1, nrow(new_rows)))

    new_data <-  dplyr::bind_cols(new_rows, empty_cols)
    new_data$FAKE <- 1
    new_data$P <- 1
    new_data[[Bottom_Position_Column]] <- new_data$GENPOS
    new_data[[Bottom_Chromosome_Column]] <-  new_data$CHROM
    new_data[[Bottom_PValue_Column]] <- new_data$P

    Bottom_Data[Bottom_Chromosome_Column] <- as.character(Bottom_Data[Bottom_Chromosome_Column])
    new_data[Bottom_Chromosome_Column] <- as.character(new_data[Bottom_Chromosome_Column])
    Bottom_Data$CHROM <- as.character(Bottom_Data$CHROM)
    new_data$CHROM <- as.character(new_data$CHROM)

    Bottom_Data <-  dplyr::bind_rows(Bottom_Data, new_data)

  }

  Bottom_Data[Bottom_Chromosome_Column] <- NULL #weird character too much space

  #Don't need this later anyway

  message("Injecting key arguments")

  final_top_args$Data <- Top_Data

  final_bottom_args$Data <- Bottom_Data

  final_top_args$Title <- user_args$Top_Title

  final_bottom_args$Title <- user_args$Bottom_Title

  message("Calling Single_Plot() on Top Data")

  Top_Plot <- do.call(.Single_Plot_original, final_top_args) +
    ggplot2::labs(x = NULL)

  message("Calling Single_Plot() on Bottom Data")

  Bottom_Plot <- do.call(.Single_Plot_original, final_bottom_args) +
    ggplot2::labs(x = NULL)

  message("Formatting top plot for Miami stacking")

  Top_Plot <-  Top_Plot + ggplot2::theme(panel.border =  ggplot2::element_blank(),
                 panel.grid.major =  ggplot2::element_blank(),
                 panel.grid.minor =  ggplot2::element_blank(),
                 panel.background =  ggplot2::element_blank(),
                 axis.line =  ggplot2::element_line(),
                 axis.title.x =  ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(size = final_top_args$Chromosome_Label_Size, vjust = Adjust, colour = "black"),
                 axis.ticks.length.x  =  ggplot2::unit(0.8,"cm"),
                 legend.position = "none",
                 plot.margin =  ggplot2::margin(t = 20,  # Top margin
                                                r = 20,  # Right margin
                                                b = 0,  # Bottom margin
                                                l = 75)) #Left margin


  message("Calculating unified X Axis format...")

  df <- as.data.frame(Top_Plot$data)

  middle_new_pos_values <- numeric()

  for (chrom in c(1:22, "X")) {

      df_filtered <- df %>%
      dplyr::filter(CHROM == chrom)

    # Calculate the maximum and minimum of 'new_pos' if the Top_Dataframe is not empty
    if (nrow(df_filtered) > 0) {

      max_new_pos <- max(df_filtered$new_pos, na.rm = TRUE)
      min_new_pos <- min(df_filtered$new_pos, na.rm = TRUE)
      # Compute the average of max and min new_pos
      middle_new_pos <- (max_new_pos + min_new_pos) / 2

    } else {

      middle_new_pos <- NA  # Assign NA if no rows are found for the chromosome

    }

    # Append the calculated middle 'new_pos' value to the vector
    middle_new_pos_values <- c(middle_new_pos_values, middle_new_pos)

  }

  if (is.null(final_top_args$Chromosome_Label_Drops)){

    message("Remembering default chromosome label drops")

    if (is.null(final_top_args$Chromosome_Label_Drops) && all(!(unique(Top_Data$CHROM) %in% c("X", "x", 23))))
    {

      message("No X")

      final_top_args$Chromosome_Label_Drops <- c(21)

    }else{

      message("X Present")

      final_top_args$Chromosome_Label_Drops <- c(21,22)

    }

  }

  # Ensure labels_vec is character - labels are inherited and specified in Single_Plot()
  labels_vec <- as.character(final_top_args$Chromosome_Labels)

  # Ensure breaks_vec is character for matching
  breaks_vec <- as.character(final_top_args$Chromosome_Label_Drops)

  # Make Miami_Labels with condition
  Miami_Labels <- ifelse(labels_vec %in% breaks_vec, "", labels_vec)

  message("Finalising unified scale")

  suppressMessages(
    suppressWarnings(

  Top_Plot <- Top_Plot +  ggplot2::scale_x_continuous(

                                      breaks = middle_new_pos_values,
                                      labels = factor(Miami_Labels),
                                      position = "bottom",
                                      expand =  ggplot2::expansion(mult = c(0.01, 0.01)) )

    ))

  message("Formatting bottom plot for Miami stacking")

  Bottom_Plot <- Bottom_Plot + ggplot2::theme(panel.border =  ggplot2::element_blank(),
                                              panel.grid.major =  ggplot2::element_blank(),
                                              panel.grid.minor =  ggplot2::element_blank(),
                                              panel.background =  ggplot2::element_blank(),
                                              plot.title =  ggplot2::element_blank(),
                                              axis.text.x =  ggplot2::element_blank(),
                                              axis.title.x.top =  ggplot2::element_blank(),
                                              axis.title.x.bottom =  ggplot2::element_text(size = final_bottom_args$Title_Size, hjust = 0.5),
                                              axis.ticks.x.bottom  =  ggplot2::element_blank() ,
                                              axis.line.x.bottom =  ggplot2::element_blank(),
                                              axis.text.x.top  =  ggplot2::element_blank(),
                                              axis.text.x.bottom =  ggplot2::element_blank(),
                                              axis.line =  ggplot2::element_line(),
                                              axis.ticks.length.x.top  =  ggplot2::unit(0.8,"cm"),
                                              legend.position = "none",
                                              plot.margin =  ggplot2::margin(t = 0,  # Top margin - ticks gap - bigger more
                                                                             r = 20,  # Right margin
                                                                             b = 55,  # Bottom margin - actual title area
                                                                             l = 35)) #Left margin

  message("Inverting and aligning bottom plot")

  suppressMessages(
    suppressWarnings(

  Bottom_Plot <- Bottom_Plot +  ggplot2::coord_flip() +   ggplot2::coord_trans(y = "reverse",  clip = "off")

    ))

  message("Matching unified scale to bottom plot")

  df <- as.data.frame(Bottom_Plot$data)

  middle_new_pos_values <- numeric()

  for (chrom in c(1:22, "X")) {

    # Filter for the current chromosome
    df_filtered <- df %>%
      dplyr::filter(CHROM == chrom)

    # Calculate the maximum and minimum of 'new_pos' if the dataframe is not empty
    if (nrow(df_filtered) > 0) {

      max_new_pos <- max(df_filtered$new_pos, na.rm = TRUE)
      min_new_pos <- min(df_filtered$new_pos, na.rm = TRUE)

      # Compute the average of max and min new_pos
      middle_new_pos <- (max_new_pos + min_new_pos) / 2

    } else {

      middle_new_pos <- NA  # Assign NA if no rows are found for the chromosome

    }

    # Append the calculated middle 'new_pos' value to the vector
    middle_new_pos_values <- c(middle_new_pos_values, middle_new_pos)
  }

  message("Finalising unified scale matching for bottom plot")


  suppressMessages(
    suppressWarnings(

  Bottom_Plot <- Bottom_Plot +  ggplot2::scale_x_continuous(

                                          breaks = middle_new_pos_values,
                                          labels = Miami_Labels,
                                          position = "top",
                                          expand =  ggplot2::expansion(mult = c(0.01, 0.01)),
                                          sec.axis =  ggplot2::dup_axis(name = final_bottom_args$Title))


    ))

  message("Grobbing Top")
  #gA <-  ggplot2::ggplotGrob(Top_Plot)

  message("Grobbing Bottom")
  #gB <-  ggplot2::ggplotGrob(Bottom_Plot)

  message("Combining top and bottom grobs")
  #combined_grob <- rbind(gA, gB, size = "first")

  combined_grob <- Top_Plot / Bottom_Plot +
  patchwork::plot_layout(heights = c(1, 1))  # same sizes

  return(invisible(combined_grob))

  #return(invisible(combined_grob))

}

.Miami_Plot_original <- Miami_Plot

Miami_Plot <- function(..., session = NULL) {
  call <- match.call()
  args <- list(...)
  args$session <- session

  # Compute inferred titles BEFORE assigning .dots:
  if (is.null(args$Top_Title) && !is.null(args$Top_Data)) {

    if (is.character(args$Top_Data) && length(args$Top_Data) == 1 && file.exists(args$Top_Data)) {

      args$Top_Title <- basename(args$Top_Data)
    } else if (!is.null(call$Top_Data)) {

      args$Top_Title <- deparse(call$Top_Data)

    } else {

      args$Top_Title <- "Top Data"

    }

  }

  if (is.null(args$Bottom_Title) && !is.null(args$Bottom_Data)) {

    if (is.character(args$Bottom_Data) && length(args$Bottom_Data) == 1 && file.exists(args$Bottom_Data)) {

      args$Bottom_Title <- basename(args$Bottom_Data)

    } else if (!is.null(call$Bottom_Data)) {

      args$Bottom_Title <- deparse(call$Bottom_Data)

    } else {

      args$Bottom_Title <- "Bottom Data"

    }
  }

  args$.dots <- args

  message("Processing top dataset: ", args$Top_Title)
  message("Processing bottom dataset: ", args$Bottom_Title)

  verbose_mode <- if ("Verbose" %in% names(args)) isTRUE(args$Verbose) else FALSE

  if (verbose_mode) {

    return(do.call(.Miami_Plot_original, args))

  } else {

    return(
      suppressMessages(
        suppressWarnings(
          run_with_counter(.Miami_Plot_original, args = args, session = session)
        )
      )
    )

  }
}

