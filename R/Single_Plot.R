
#' Single, Annotated Manhattan Plot
#'
#' @param Data These are the summary statistics (either file path or environmental variable) to be plotted; defaults to NULL
#' @param Chromosome_Column Manually specify chromosome column; defaults to NULL, leading to auto-detection
#' @param Position_Column Manually specify chromosome position column; defaults to NULL, leading to auto-detection
#' @param SNP_ID_Column Manually specify SNP ID column; defaults to NULL, leading to auto-detection
#' @param PValue_Column Manually specify P Value column; defaults to NULL, leading to auto-detection
#' @param Reference_Allele_Column Manually specify reference allele column; defaults to NULL, leading to auto-detection
#' @param Effect_Allele_Column Manually specify effect allele column; defaults to NULL, leading to auto-detection
#' @param Title Title to be displayed; defaults to NULL, leading to auto-detection of data set/file name
#' @param Y_Axis_Title Y axis title; defaults to "-log₁₀P"
#' @param X_Axis_Title X axis title; defaults to "Chromosome"
#' @param Title_Size Size of title; defaults to 35
#' @param Y_Axis_Text_Size Size of Y axis text; defaults to 30
#' @param Chromosome_Label_Size Size of chromosome (X axis text) numbers; defaults to 30
#' @param Y_Axis_Title_Size Size of Y axis title; defaults to 35
#' @param X_Axis_Title_Size Size of X axis title; defaults to 35
#' @param Chromosome_Colours Alternating list of two colours to allocate for data points (odd/even); defaults to c("blue", "turquoise")
#' @param Sig_Line_Colour Colour of significance threshold; defaults to "red"
#' @param Sig_Line_Type Type of significance threshold line; defaults to "dashed"
#' @param Sig_Threshold Value of significance threshold; defaults to standard Bonferroni correction of 5e-8
#' @param Sig_Line_Width Width of significance threshold line; defaults to 0.5
#' @param Point_Size Size of data points; defaults to 2.5
#' @param Label_Index Label the index SNPs with ID provided; defaults to TRUE
#' @param Label_Size Size of index labels if used; defaults to 5
#' @param Label_Angle Angle of index labels if used; defaults to 45
#' @param Label_Colour Colour of index SNP label if used; defaults to "black"
#' @param Colour_Index Highlight the index SNP; defaults to TRUE
#' @param Colour_Of_Index Colour of index SNP highlight if used; defaults to "red"
#' @param Chromosome_Labels List of chromosomes to show index SNP label on, given that Label_Index is TRUE; defaults to c(1:22, "X")
#' @param Chromosome_Index  List of chromosomes to show index SNP highlight/colour on, given that Colour_Index is TRUE; defaults to c(1:22, "X")
#' @param Condense_Scale Condense scale; defaults to TRUE
#' @param Break_Point Value where y axis condenses; defaults to 1e-10, frequent max value observed in a small sample GWAS
#' @param Draft_Plot Specify whether a random subset of SNPs, the number of which is specified by Random_Selection below, is plotted; defaults to FALSE
#' @param Random_Selection Manually select a random subset of SNPs to plot (more for internal debugging and testing), given that Draft_Plot is TRUE; defaults to 10000
#' @param Chromosome_Label_Drops Manually specify which chromosome/X axis labels (as a list) should not be shown to prevent crowding; defaults to c(21,22)
#' @param Interactive View the output as an interactive plotly object, which can be used to inspect data, with key overlays (more for shiny App usage); defaults to FALSE
#'
#' @return Image of Single Manhattan Plot is saved to specified object and the resulting ggplot object can then be saved to an image
#' @export
#'
#' @examples  Manhattan_Plot <- Single_Plot(Data = Intelligence_Sum_Stats)
#' @examples  Manhattan_Plot <- Single_Plot(Data = "C:/Users/callumon/Miami_Package_R/MiamiR/Intelligence_Sum_Stats_Mini.txt")

Single_Plot<- function(Data = NULL,
                       Random_Selection = 10000,
                       Draft_Plot = FALSE,
                       Title = NULL,
                       Title_Size = 35,
                       Y_Axis_Text_Size = 30,
                       Y_Axis_Title_Size = 35,
                       X_Axis_Title_Size = 35,
                       Chromosome_Label_Size = 30,
                       Chromosome_Label_Drops = NULL,  # can be numeric or character
                       X_Axis_Title = NULL,
                       Y_Axis_Title = "-log₁₀(P)",
                       Chromosome_Colours = c("blue", "turquoise"),
                       Sig_Line_Colour = "red",
                       Sig_Line_Type = "dashed",
                       Sig_Threshold = 5e-8,
                       Sig_Line_Width = 0.5,
                       Point_Size = 2.5,
                       Index_Size = NULL,
                       Index_Highlight_Shape_Size = 50,
                       Index_Thickness = 1,
                       Chromosome_Labels = c(1:22, "X"),
                       Chromosome_Index = c(1:22, "X"),
                       Chromosome_Diamond = c(1:22, "X"),
                       Anchor_Label = "left",
                       Label_Index = TRUE,
                       Label_Size = 5, Label_Angle = 45,
                       Label_Colour = "black",
                       Label_Height = 6,
                       Colour_Of_Index = "darkred",
                       Colour_Index = TRUE,
                       Diamond_Index = NULL,
                       Colour_Of_Diamond = "purple",
                       Diamond_Index_Size = NULL,
                       Condense_Scale = TRUE,
                       Break_Point = 1e-10,
                       Chromosome_Column = NULL,
                       Position_Column = NULL,
                       SNP_ID_Column = NULL,
                       PValue_Column = NULL,
                       Reference_Allele_Column = NULL,
                       Effect_Allele_Column = NULL,
                       Interactive = FALSE,
                       Lab = NULL,
                       Verbose = FALSE)

{



  is_called_by_regional_plot <- any(vapply(sys.calls(), function(x) {
    "Regional_Plot" %in% deparse(x[[1]])
  }, logical(1)))


  #Points need to be bigger if called as part of Regional_Plot()
  #If left as default this will occur:

  if (is_called_by_regional_plot && missing(Point_Size)) {
    Point_Size <- 10
  }
  if (is_called_by_regional_plot && missing(Diamond_Index)) {
    Diamond_Index <- TRUE
  }
  if (is_called_by_regional_plot && missing(Label_Size)) {
    Label_Size <- 9
  }
  if (is_called_by_regional_plot && missing(Label_Height)) {
    Label_Height <- 25
  }
  if (is_called_by_regional_plot && missing(Index_Thickness)) {
    Index_Thickness <- 4
  }
  if (is_called_by_regional_plot && missing(Diamond_Index_Size)) {
    Diamond_Index_Size <- Point_Size * 2
  }



  #Standard X2 looks good, but can be modified/
  if(is.null(Index_Size))
  {
    Index_Size <- Point_Size * 2
  }


  if(is.null(Diamond_Index_Size))
  {
    Diamond_Index_Size <- Point_Size * 2
  }

  if(is.null(Diamond_Index))
  {
    Diamond_Index <- FALSE
  }

  if (is.null(Data)) {
    stop("Please provide data.", call. = FALSE)
  }




#
#   if (is.null(Title)) {
#     if (is.character(Data) && length(Data) == 1 && file.exists(Data)) {
#       # Extract file name without path or extension
#
#
#       Title <-   tools::file_path_sans_ext(basename(Data))
#     } else {
#
#
#       # Fallback to variable name if Data is not a path
#      Title <- (substitute(Data)) #removed desparse
#
#
#
# print(Title)
#
#      Title <- if (is.symbol(Title)) as.character(Title) else NULL
#
#
#   #    Title <-  paste(deparse(substitute(Data)), collapse = "") #get_name_fast(Data)
#     }
#   }
#
#
#   print(Title)
# z



if(is.null(X_Axis_Title))
{
  X_Axis_Title <- "Chromosome"
}

    if (is.character(Data) && length(Data) == 1) {
    file_path <- Data



     message("Dataset absent from environment")

    if (file.exists(file_path)) {

      message("Reading data from file: ", file_path)


      message("Loading data using vroom...")

      Data <- vroom::vroom(file_path, show_col_types = FALSE)

      message("✅ Finished reading")




    } else {
      stop("The provided character string does not point to an existing file: ", file_path, call. = FALSE)
    }
    }


  # Handle multiple file paths (vector of strings)
  # if (is.character(Data) && length(Data) > 1) {
  #
  #
  #
  #   message("Multiple file paths detected. Converting to named list of data frames.")
  #
  #   # Set default names from filenames if not named already
  #   names(Data) <- basename(tools::file_path_sans_ext(Data))
  #
  #   # Read each file into a data frame
  #   Data <- lapply(Data, function(path) {
  #     if (!file.exists(path)) stop("File does not exist: ", path)
  #     vroom::vroom(path, show_col_types = FALSE)
  #   })
  #
  #   # Recursively call outer Single_Plot with the list of loaded data
  #   return(Single_Plot(Data = Data))
  # }






  if(Draft_Plot == T)
  {

    message("Generating random subset of points specified")

    Data <- dplyr::sample_n(Data, Random_Selection)

  }






  Chromosome_Column <- detect_chromosome_column(Data, Chromosome_Column)


  Data <- Data %>%
    dplyr::mutate(!!Chromosome_Column := ifelse(.data[[Chromosome_Column]] == "23", "X", .data[[Chromosome_Column]]))


  PValue_Column     <- detect_pvalue_column(Data, PValue_Column)
  Position_Column   <- detect_position_column(Data, Position_Column)
  SNP_ID_Column     <- detect_snp_column(Data, SNP_ID_Column)
  Ref_Allele_Column     <- detect_reference_allele_column(Data, Reference_Allele_Column)
  Alt_Allele_Column     <- detect_effect_allele_column(Data, Effect_Allele_Column)




  #Manually assign columns as certain names for ease of use
  Data$CHROM <- Data[[Chromosome_Column]]
  Data$GENPOS <- Data[[Position_Column]]
  Data$ID <- Data[[SNP_ID_Column]]
  Data$P <- Data[[PValue_Column]]
  Data$ALLELE0 <- Data[[Ref_Allele_Column]]
  Data$ALLELE1 <- Data[[Alt_Allele_Column]]



  if (is.null(Chromosome_Label_Drops)){


  if (is.null(Chromosome_Label_Drops) && all(!(unique(Data$CHROM) %in% c("X", "x", 23))))
  {

    Chromosome_Label_Drops <- c(21)

  }else{


    Chromosome_Label_Drops <- c(21,22)
  }


  }

  max_logp <- max(-log10(Data[[PValue_Column]]), na.rm = TRUE)



  #fake anchoring around single point to prevent scale missing

  if (nrow(Data) == 1) {

    message("Only one point provided - anchoring")

    # Original row
    original_row <- Data[1, ]

    # Dynamically get the column names
    pos_col <- Position_Column
    pval_col <- PValue_Column

    # Create a copy with position - 1
    fake_minus <- original_row
    fake_minus[[pos_col]] <- fake_minus[[pos_col]] - 1
    fake_minus[[pval_col]] <- 0.9
    fake_minus$ID <- "FAKEMINUS1"

    # Create a copy with position + 1
    fake_plus <- original_row
    fake_plus[[pos_col]] <- fake_plus[[pos_col]] + 1
    fake_plus[[pval_col]] <- 0.9
    fake_plus$ID <- "FAKEPLUS1"

    # Combine all
    Data <- dplyr::bind_rows(original_row, fake_minus, fake_plus)
  }


  if(!is.null(Lab))
  {

  }

  if(("Lab" %in% colnames(Data))) # if annotate used before then skip
  {
    message("It looks like you've applied Annotate_Data() or provided specific lead IDs")
  }


  if(!("Lab" %in% colnames(Data))) # if annotate used before then skip
  {

  message("Discerning Lead SNP Per Chromosome to annotate")

  Data <- Data %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(min_P_GENPOS = GENPOS[which.min(P)]) %>% #P is formed in function
    dplyr::ungroup()




  # Label the minimum P value for each CHROM
  Data <- Data %>%
    dplyr::mutate(Lab = ifelse(GENPOS == min_P_GENPOS & P < 5e-8, ID,
                        ifelse(abs(GENPOS - min_P_GENPOS) > 2000000000 & P == min(P[abs(GENPOS - min_P_GENPOS) > 2000000000]),
                               ID, "")))



  }


  #return(Data)



  if(Diamond_Index == TRUE)
  {



    #top already from regional sometimes
    if (!("top" %in% colnames(Data))) {
      Data <- Data %>%
        dplyr::group_by(CHROM) %>%
        dplyr::mutate(
          min_P_GENPOS = GENPOS[which.min(P)],
          top = (GENPOS == min_P_GENPOS)
        ) %>%
        dplyr::ungroup()
    }



  #  print(table(Data$top))

  Data2 <- Data
  # Store rows where top == TRUE
  TOPS <- dplyr::filter(Data, top == TRUE)

  # Remove rows where top == TRUE from Data
 # Data <- dplyr::filter(Data, top == FALSE)

  }

  Data$Lab[Data$Lab == ""] <- NA
  Data$COLOUR[is.na(Data$Lab)] <- NA
  Data$COLOUR[!is.na(Data$Lab)] <- 4 #slightly smaller ring as bigger shape
  Data$COLOURTRY <- NA
  Data$COLOURTRY[is.na(Data$Lab)] <- NA
  Data$COLOURTRY[!is.na(Data$Lab)] <- "pink"



  Data$Lab[!(Data$CHROM %in% Chromosome_Labels)] <- NA



  Data$COLOUR[!(Data$CHROM %in% Chromosome_Index)] <- NA



  Title <- paste0(Title, "\n \n ") # plotly only interprets for y if actual content on line in app



 # Title <- "1"
#  system.time(Title <- paste0(Title, "\n\n"))




  X_Axis_Title <- paste0("\n ", X_Axis_Title) # needs to go above it to push down opp to above.
  Y_Axis_Title <- paste0(Y_Axis_Title, "\n ") # needs to go left above it to push down opp to above.





  #Basic Plot




  message("Plotting Top Framework")




  if(Interactive == TRUE)

  {

    Data$Hover_Info <- paste0(
      "SNP: ", Data$ID, "\n",
      "CHR: ", Data$CHROM, "\n",
      "POS: ", Data$GENPOS, "\n",
      "P: ", signif(Data$P, 4), "\n",
      "REF: ", Data$ALLELE0, " ALT: ", Data$ALLELE1
    )



}




# Compute max -log10(P)
max_logp <- max(-log10(Data[[PValue_Column]]), na.rm = TRUE)


Break_Point_LOG <- -log10(Break_Point)

# Conditionally turn off Condense_Scale
if (max_logp < Break_Point_LOG) {
  Condense_Scale <- FALSE
}

# if (nrow(Data) <= 1)
# {
#
#   message("ggmanh kinda engaged:")
#      a <- ggmanh::manhattan_plot(x = Data, preserve.position = T, plot.title = Title,
#                          chr.colname = Chromosome_Column, pos.colname = Position_Column, label.colname = NULL,
#                          pval.colname = PValue_Column, annotateTop = FALSE,
#                          chr.order = c(10),
#                          chr.col = Chromosome_Colours, chrlabs = c(1:22, "X"), rescale = Condense_Scale,
#                          signif = Break_Point, rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
#                          signif.col = c("transparent"),  point.size = Point_Size, x.label = "", y.label = "")
#
# }else{


message("ggmanh engaged:")



if(length(unique(Data$CHROM)) == 1)

{

  chroms <- unique(Data$CHROM)

#
#   set.seed(123)  # Optional reproducibility
#   n <- nrow(Data)
#   half_n <- floor(n / 2)
#   labels <- c(rep("yes", half_n), rep("no", n - half_n))
#   Data$is_highlight <- sample(labels)  # Shuffle
#
#   # Check distribution:
#   table(Data$is_highlight)
#

  #pos preserve needs false as normalised coords bad here and slighlty off!

  a <- ggmanh::manhattan_plot(x = Data, preserve.position = T, plot.title = Title,
                      chr.colname = Chromosome_Column, pos.colname = Position_Column, label.colname = NULL,
                      pval.colname = PValue_Column, annotateTop = FALSE,
                #      chr.order = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"),
               chr.order = c(chroms),
                      chr.col = Chromosome_Colours, chrlabs = c(1:22, "X"), rescale = Condense_Scale,
                      signif = Break_Point, rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
                      signif.col = c("transparent"),  point.size = Point_Size, x.label = "", y.label = "")






  # Define alternating chromosome colours:
  unique_chr <- sort(unique(a$data$CHROM))
  chr_colour_map <- setNames(
    rep(Chromosome_Colours, length.out = length(unique_chr)),
    unique_chr
  )



  # Assign alternating colours:
  a$data$colour_group <- chr_colour_map[as.character(a$data$CHROM)]

  if(Diamond_Index == TRUE)
  {
  # Override for top==TRUE:
  a$data$colour_group[a$data$top == TRUE] <- "transparent"
}
  # Update mapping:
  a$mapping$colour <- rlang::quo(colour_group)

  # Auto-generate scale based on actual values present:
  colour_vals <- unique(a$data$colour_group)
  colour_map <- setNames(colour_vals, colour_vals)


  suppressMessages({

  a <- a + ggplot2::scale_colour_manual(values = colour_map)

})


  if(Diamond_Index == TRUE)
  {

  a2 <- ggmanh::manhattan_plot(x = Data2, preserve.position = T, plot.title = Title,
                              chr.colname = Chromosome_Column, pos.colname = Position_Column, label.colname = NULL,
                              pval.colname = PValue_Column, annotateTop = FALSE,
                              #      chr.order = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"),
                              chr.order = c(chroms),
                              chr.col = Chromosome_Colours, chrlabs = c(1:22, "X"), rescale = Condense_Scale,
                              signif = Break_Point, rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
                              signif.col = c("transparent"),  point.size = Point_Size, x.label = "", y.label = "")

  }

}else{



  set.seed(123)  # Optional reproducibility
  n <- nrow(Data)
  half_n <- floor(n / 2)
  labels <- c(rep("yes", half_n), rep("no", n - half_n))
  Data$is_highlight <- sample(labels)  # Shuffle

  # Check distribution:
  table(Data$is_highlight)




  a <- ggmanh::manhattan_plot(x = Data, preserve.position = T, plot.title = Title,
                              chr.colname = Chromosome_Column, pos.colname = Position_Column, label.colname = NULL,
                              pval.colname = PValue_Column, annotateTop = FALSE,
                                    chr.order = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"),
                            #  chr.order = c(chroms),
                              chr.col = Chromosome_Colours, chrlabs = c(1:22, "X"), rescale = Condense_Scale,
                              signif = Break_Point, rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
                              signif.col = c("transparent"),  point.size = Point_Size, x.label = "", y.label = "")





  # # 1️⃣ Ensure CHROM/CHR consistency:
  # if ("CHR" %in% colnames(TOPS)) {
  #   TOPS$CHR <- "DMD"
  # }
  #
  # if ("CHROM" %in% colnames(TOPS)) {
  #   TOPS$CHROM <- "DMD"
  # }
  #
  #
  # print(TOPS)
  # print(a$data)
  #
  # # 2️⃣ Append to `a$data`:
  # a$data <- dplyr::bind_rows(a$data, TOPS)
  #
  # print(a$data)

#
#   unique_chr <- unique(a$data$CHR)
#   # Example palette for all CHRs:
#   # (This assumes alternating blue/turquoise for CHRs and pink for CHR==6)
#
#   base_colours <- c("blue", "turquoise")
#
#   chr_colours <- setNames(
#     rep(base_colours, length.out = length(unique_chr)),
#     unique_chr
#   )
#
#   # Force CHR==6 to pink:
#   chr_colours["DMD"] <- "pink"
#   a <- a + ggplot2::scale_colour_manual(values = chr_colours)
#
#

  # a$data$colour_group <- ifelse(a$data$log10pval > 5, "high", "low")
  # a$mapping$colour <- rlang::quo(colour_group)
  # a <- a +
  #   ggplot2::scale_colour_manual(values = c("high" = "yellow", "low" = "red"))


  # Define alternating chromosome colours:
  unique_chr <- sort(unique(a$data$CHROM))
  chr_colour_map <- setNames(
    rep(Chromosome_Colours, length.out = length(unique_chr)),
    unique_chr
  )

  # Assign alternating colours:
  a$data$colour_group <- chr_colour_map[as.character(a$data$CHROM)]

  if(Diamond_Index == TRUE)
  {
  # Override for top==TRUE:
  a$data$colour_group[a$data$top == TRUE & a$data$CHROM %in% Chromosome_Diamond] <- "transparent"
 # return(a)
  #Z
}

  # Update mapping:
  a$mapping$colour <- rlang::quo(colour_group)

  # Auto-generate scale based on actual values present:
  colour_vals <- unique(a$data$colour_group)
  colour_map <- setNames(colour_vals, colour_vals)


  suppressMessages({
    a <- a + ggplot2::scale_colour_manual(values = colour_map)
  })


  # ,
  # highlight.colname = "is_highlight",
  # highlight.col = c("yes" = "blue", "no" = "grey"),
  # color.by.highlight = TRUE

 # return(a)

  if(Diamond_Index == TRUE)
  {

  a2 <- ggmanh::manhattan_plot(x = Data2, preserve.position = T, plot.title = Title,
                              chr.colname = Chromosome_Column, pos.colname = Position_Column, label.colname = NULL,
                              pval.colname = PValue_Column, annotateTop = FALSE,
                              chr.order = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"),
                              #  chr.order = c(chroms),
                              chr.col = Chromosome_Colours, chrlabs = c(1:22, "X"), rescale = Condense_Scale,
                              signif = Break_Point, rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
                              signif.col = c("transparent"),  point.size = Point_Size, x.label = "", y.label = "")

  }

}



  # a <- ggmanh::manhattan_plot(x = Data,
  # # , preserve.position = T, plot.title = Title,
  #                              chr.colname = Chromosome_Column, pos.colname = Position_Column, label.colname = NULL,
  #                              pval.colname = PValue_Column, annotateTop = FALSE, #)#,
  #                              chr.order = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"))#,
  # #                             chr.col = Chromosome_Colours, chrlabs = c(1:22, "X"), rescale = Condense_Scale,
  # #                             signif = Break_Point, rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
  # #                             signif.col = c("transparent"),  point.size = Point_Size, x.label = "", y.label = "")

  #Significance Line


#  Data <- head(Data, 1000)

    # a <- ggplot(Data, aes(x = GENPOS, y = -log10(P), text = Hover_Info)) +
    # geom_point(size = 1.5, alpha = 0.7) +
    # theme_minimal() +
    # labs(
    #   title = "GENPOS vs -log10(P)",
    #   x = "Genomic Position",
    #   y = "-log10(P-value)"
    # )
    #


#Plot_Outcome <- a

#return(Plot_Outcome)





  message("Adding Significance Line")

  b <- a + ggplot2::geom_hline(yintercept = -(log10(Sig_Threshold)), linetype = Sig_Line_Type,
                      color = Sig_Line_Colour, size = Sig_Line_Width)

  if(Diamond_Index == TRUE)
  {

  b2 <- a2 + ggplot2::geom_hline(yintercept = -(log10(Sig_Threshold)), linetype = Sig_Line_Type,
                               color = Sig_Line_Colour, size = Sig_Line_Width)

  }
  #Theme/formatting


  message("Formatting")

  p <- b + ggplot2::theme(panel.border =  ggplot2::element_blank(),
                 panel.grid.major =  ggplot2::element_blank(),
                 panel.grid.minor =  ggplot2::element_blank(),
                 panel.background =  ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(hjust = 0.5, size = Title_Size),
                 axis.text.x = ggplot2::element_text(size = Chromosome_Label_Size, vjust = -0.5, colour = "black"),
                 axis.line =  ggplot2::element_line(),
                 axis.title.x =  ggplot2::element_text(size = X_Axis_Title_Size, vjust = -7 ),
                 axis.title.y =  ggplot2::element_text(size = Y_Axis_Title_Size, vjust = 7),
                 axis.text.y =  ggplot2::element_text(size = Y_Axis_Text_Size, colour = "black"),
                 axis.ticks.length.x  =  ggplot2::unit(0.8,"cm"),
                 legend.position = "none",
                 plot.margin =  ggplot2::margin(t = 20,  # Top margin
                                      r = 40,  # Right margin
                                      b = 60,  # Bottom margin
                                      l = 60)) #Left margin


  if(Diamond_Index == TRUE)
  {
  p2 <- b2 + ggplot2::theme(panel.border =  ggplot2::element_blank(),
                          panel.grid.major =  ggplot2::element_blank(),
                          panel.grid.minor =  ggplot2::element_blank(),
                          panel.background =  ggplot2::element_blank(),
                          plot.title = ggplot2::element_text(hjust = 0.5, size = Title_Size),
                          axis.text.x = ggplot2::element_text(size = Chromosome_Label_Size, vjust = -0.5, colour = "black"),
                          axis.line =  ggplot2::element_line(),
                          axis.title.x =  ggplot2::element_text(size = X_Axis_Title_Size, vjust = -7 ),
                          axis.title.y =  ggplot2::element_text(size = Y_Axis_Title_Size, vjust = 7),
                          axis.text.y =  ggplot2::element_text(size = Y_Axis_Text_Size, colour = "black"),
                          axis.ticks.length.x  =  ggplot2::unit(0.8,"cm"),
                          legend.position = "none",
                          plot.margin =  ggplot2::margin(t = 20,  # Top margin
                                                         r = 40,  # Right margin
                                                         b = 60,  # Bottom margin
                                                         l = 60)) #Left margin

  }

  p <- p +  ggplot2::xlab(X_Axis_Title) + ggplot2::theme(axis.title.x =  ggplot2::element_text(vjust=-5, size = X_Axis_Title_Size, hjust = 0.5))
  p <- p +  ggplot2::ylab(Y_Axis_Title) + ggplot2::theme(axis.title.y =  ggplot2::element_text(vjust=5, size = Y_Axis_Title_Size))


  if(Diamond_Index == TRUE)
  {

  p2 <- p2 +  ggplot2::xlab(X_Axis_Title) + ggplot2::theme(axis.title.x =  ggplot2::element_text(vjust=-5, size = X_Axis_Title_Size, hjust = 0.5))
  p2 <- p2 +  ggplot2::ylab(Y_Axis_Title) + ggplot2::theme(axis.title.y =  ggplot2::element_text(vjust=5, size = Y_Axis_Title_Size))


  }
 # p <- p +  ggplot2::ylab(Y_Axis_Title) + ggplot2::theme(axis.title.y =  ggplot2::element_text(vjust=5, size = Y_Axis_Title_Size))

#
#   p <- p +  ggplot2::ylab("WAAA-log₁₀P\u2002\u2002\u2002\u2002\u2002\u2002\u2002\u2002\u2002\u2002\u2002\u2002\u2002\u2002\u2002\u2002\u2002\u2002\u2002")+  # or paste0(Y_Axis_Title, "~~~~")
#   ggplot2::theme(axis.title.y =  ggplot2::element_text(vjust=5, size = Y_Axis_Title_Size))
#
#
 #  p <- p +  labs(y = "-log₁₀P\n \n ")  # Add a "W" just to see effect

  if(Interactive == TRUE)
  {
  dot_padding <- "\n \n "  # You can add more lines if needed
 Y_Axis_Title <- paste0(Y_Axis_Title, dot_padding)
 X_Axis_Title <- paste0(dot_padding, X_Axis_Title)


  p <- p +  labs(y = Y_Axis_Title, x = X_Axis_Title)

  if(Diamond_Index == TRUE)
  {
  p2 <- p2 +  labs(y = Y_Axis_Title, x = X_Axis_Title)

  }
  }
  #La belling interact

  if(Interactive == TRUE)

  {
  p <- p + ggplot2::aes(text = Hover_Info)
  if(Diamond_Index == TRUE)
  {
  p2 <- p2 + ggplot2::aes(text = Hover_Info)
  }
}


  if(Anchor_Label == "left")
  {
    HJUST = 0
    VJUST = 0
  }

  if(Anchor_Label == "centre")
  {
    HJUST = 0
    VJUST = 0.5
  }
  if(Anchor_Label == "right")
  {
    HJUST = 0
    VJUST = 1
  }


  if(Label_Index == TRUE)
{
    message("Adding Labels")
    p <- p +  ggtext::geom_richtext( ggplot2::aes(label = Lab), hjust = HJUST, vjust = VJUST, label.size = 0,
                                    label.padding = grid::unit(rep(Label_Height, 4), "pt"), nudge_y = 0, fill = NA,
                                    nudge_x = 0, color = Label_Colour, label.color = NA, size = Label_Size, angle = Label_Angle)

    if(Diamond_Index == TRUE)
    {

    p2 <- p2 +  ggtext::geom_richtext( ggplot2::aes(label = Lab), hjust = HJUST, vjust = VJUST, label.size = 0,
                                     label.padding = grid::unit(rep(Label_Height, 4), "pt"), nudge_y = 0, fill = NA,
                                     nudge_x = 0, color = Label_Colour, label.color = NA, size = Label_Size, angle = Label_Angle)

    }

  }




  suppressMessages(suppressWarnings({


  if(Colour_Index == TRUE & Interactive == FALSE & Diamond_Index == FALSE)
  {


    df <- as.data.frame(p$data)

    special_points <- df[df$COLOUR == 4, ]

    if(is_called_by_regional_plot == TRUE)

    {



    p <- p + ggplot2::geom_point(
      data = special_points,
      ggplot2::aes(x = GENPOS, y = log10pval),  #need GENPOS for inherited in regional
      color = Colour_Of_Index,
      size = Index_Size,
      stroke = Index_Thickness,
      pch = 21
    )

    }else{

      p <- p + ggplot2::geom_point(
        data = special_points,
        ggplot2::aes(x = new_pos, y = log10pval),
        color = Colour_Of_Index,
        size = Index_Size,
        stroke = Index_Thickness,
        pch = 21
      )



    }

}


    if(Diamond_Index == TRUE)
    {
    # Add them back (at the end, if row order is not important)
#    Data <- dplyr::bind_rows(Data, TOPS)


   # print(TOPS)
    message("Diamond Index SNPs")



    df <- as.data.frame(p2$data)




    special_points <- df[df$top == TRUE, ]


    if (!("log10pval" %in% colnames(special_points))) {
      special_points <- special_points %>%
        dplyr::mutate(log10pval = 5)
    }



#    print(special_points)




    #filter for only ones to be added




    special_points <- special_points[special_points$CHROM %in% Chromosome_Diamond, ]


  #  special_points$log10pval <- 24



    if(is_called_by_regional_plot == TRUE)

    {



      suppressMessages({

      p <- p + ggplot2::geom_point(
        data = special_points,
        ggplot2::aes(x = GENPOS, y = log10pval),  #need GENPOS for inherited in regional
        color = Colour_Of_Diamond,
        size = Diamond_Index_Size,
        alpha = 0.5,
        pch = 18
      )

      })

    }else{


      suppressMessages({

      p <- p + ggplot2::geom_point(
        data = special_points,
        ggplot2::aes(x = new_pos, y = log10pval),
        color = Colour_Of_Diamond,
        size = Diamond_Index_Size,
        alpha = 0.5,
        pch = 18
      )
})


    }


  }



  p <- p +  ggplot2::coord_cartesian(clip = "off")


  }))

  message("Scaling Axes")



  df <- as.data.frame(p$data)



  middle_new_pos_values <- numeric()


  #currently need to adjust to X manually update in future.

  for (chrom in c(1:22, "X")) {
    # Filter for the current chromosome
    df_filtered <- df %>%
      dplyr::filter(CHROM == chrom)

    # Calculate the maximum and minimum of 'new_pos' if the Dataframe is not empty
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






  # Get unique chromosomes from the data
  chroms_in_data <- unique(as.character(Data$CHROM))

  # Define full label set (assuming you want 1-22 and X)
  all_chromosomes <- as.character(c(1:22, "X"))

  # Filter out dropped and missing chromosomes
  chromosome_labels <- setdiff(all_chromosomes, as.character(Chromosome_Label_Drops))
  chromosome_labels <- intersect(chromosome_labels, chroms_in_data)

  # Ensure label order matches middle_new_pos_values
  # We'll assume your middle_new_pos_values vector is already in chromosomal order:
  final_labels <- ifelse(all_chromosomes %in% chromosome_labels, all_chromosomes, "")


  suppressMessages(suppressWarnings({

  d <- p +  ggplot2::scale_x_continuous(breaks = middle_new_pos_values,
                              # Midpoints or specific points for chromosome labels
                            #  labels = factor(c(1:20, "", "", "X")),
                            labels = final_labels,

                              position = "bottom",
                              expand =  ggplot2::expansion(mult = c(0.01, 0.01)) )




  }))

  # d <- d + ggplot2::coord_cartesian(
  #   ylim = c(0, max_logp + 100),
  #   clip = "off"
  # )

  Plot_Outcome <- d



  if (is_called_by_regional_plot) {


  idx <- which(sapply(Plot_Outcome$scales$scales, function(s) "y" %in% s$aesthetics))
  Plot_Outcome$scales$scales[[idx]]$expand <- ggplot2::expansion(mult = c(0.05, 0))

  }


  message("Plot Complete")


#  Plot_Outcome <- Plot_Outcome + ggplot2::scale_y_continuous(
 #   expand = ggplot2::expansion(mult = c(0.05, 0))  # 5% padding at bottom, none at top # 0.5 units at bottom only
#  )


  return(invisible(Plot_Outcome)) # invisible prevents slow print auto at end.

}



.Single_Plot_original <- Single_Plot


Single_Plot <- function(..., session = NULL) {
  call_expr <- match.call()
  fn_formals <- formals(.Single_Plot_original)
  valid_args <- names(fn_formals)
  dots <- list(...)
  clean_args <- dots[names(dots) %in% valid_args]

  # Capture raw unevaluated expression for Data
  raw_data_expr <- call_expr[["Data"]]

  # Handle Data = c(...) specifically and convert to named list
  # Inside Single_Plot wrapper (after call_expr and raw_data_expr)

  # Handle c("file1", "file2") input
  # Inside Single_Plot wrapper
  if (!is.null(raw_data_expr) &&
      is.call(raw_data_expr) &&
      identical(raw_data_expr[[1]], as.name("c"))) {

    message("Detected c(...) input")

    arg_exprs <- as.list(raw_data_expr[-1])
    arg_vals <- lapply(arg_exprs, function(e) eval(e, parent.frame()))

    # Case 1: All are data frames
    if (all(vapply(arg_vals, is.data.frame, logical(1)))) {
      message("Handling multiple data frames from c(...)")
      arg_names <- vapply(arg_exprs, deparse, character(1))
      names(arg_vals) <- arg_names
      clean_args$Data <- arg_vals

      # Case 2: All are character paths to existing files
    } else if (all(vapply(arg_vals, is.character, logical(1)))) {

      file_vec <- unlist(arg_vals)

      if (all(file.exists(file_vec))) {
        message("Handling multiple file paths from c(...)")

        file_names <- basename(tools::file_path_sans_ext(file_vec))
        data_list <- lapply(file_vec, function(f) {
          vroom::vroom(f, show_col_types = FALSE)
        })
        names(data_list) <- file_names
        clean_args$Data <- data_list

      } else {
        stop("❌ One or more file paths in `c(...)` do not exist.", call. = FALSE)
      }

    } else {
      stop("❌ Mixed or unrecognized types in `c(...)` — expected all data.frames or all file paths.", call. = FALSE)
    }
  }



  Data <- clean_args$Data

  is_df <- function(x) is.data.frame(x)

  # Handle character vector of multiple file paths
  if (is.character(Data) && length(Data) > 1 && all(file.exists(Data))) {
    message("Reading multiple file paths into list of data frames")
    Data_list <- lapply(Data, function(path) {
      vroom::vroom(path, show_col_types = FALSE)
    })
    names(Data_list) <- basename(tools::file_path_sans_ext(Data))
    clean_args$Data <- Data_list
    Data <- Data_list
  }


  # Multiple data frames case
  if (is.list(Data) && all(vapply(Data, is_df, logical(1)))) {
    message("Multiple data frames detected. Running Single_Plot on each.")

    data_names <- names(Data)
    if (is.null(data_names) || any(data_names == "")) {
      data_names <- vapply(Data, function(d) {
        deparse(substitute(d))
      }, character(1))
      names(Data) <- data_names
    }
    plots <- lapply(seq_along(Data), function(i) {
      df <- Data[[i]]
      name <- names(Data)[i]

      args_i <- clean_args
      args_i$Data <- df

      # Set inferred title if not already set
      if (is.null(args_i$Title) || !nzchar(args_i$Title)) {
        args_i$Title <- name
   #     message("📢 [Silent] Title being passed: ", args_i$Title)
      }

      message("Calling plot for: ", args_i$Title)

      # Ensure all defaults are set (including Chromosome_Labels, etc.)
      for (arg in setdiff(valid_args, names(args_i))) {
        if (arg == "Title" && !is.null(args_i$Title) && nzchar(args_i$Title)) next

        default_val <- fn_formals[[arg]]
        if (!identical(default_val, quote(expr = ))) {
          val <- tryCatch(
            eval(default_val, envir = environment(.Single_Plot_original)),
            error = function(e) NULL
          )
          args_i[[arg]] <- val
        }
      }

      if (isTRUE(args_i$Verbose)) {
        do.call(.Single_Plot_original, args_i)
      } else {
        run_with_counter(.Single_Plot_original, args = args_i, session = session)
      }
    })


    names(plots) <- names(Data)
    return(invisible(plots))
  }

  # Handle single dataframe title inference
  if (is.null(clean_args$Title) && !is.null(raw_data_expr)) {
    clean_args$Title <- deparse(raw_data_expr)
    if (isTRUE(clean_args$Verbose)) {
  #    message("Title auto-set to: ", clean_args$Title)
    } else {
  #    print(paste("📢 [Silent] Title being passed:", clean_args$Title))
    }
  }

  # Evaluate missing defaults (skip overriding inferred Title)
  for (arg in setdiff(valid_args, names(clean_args))) {
    if (arg == "Title" && !is.null(clean_args$Title) && nzchar(clean_args$Title)) next

    default_val <- fn_formals[[arg]]
    if (!identical(default_val, quote(expr = ))) {
      val <- tryCatch(
        eval(default_val, envir = environment(.Single_Plot_original)),
        error = function(e) NULL
      )
      clean_args[[arg]] <- val
    }
  }


  clean_args$session <- session

  if (isTRUE(clean_args$Verbose)) {
 #   message("🚀 Calling .Single_Plot_original directly")
    return(do.call(.Single_Plot_original, clean_args))
  } else {
    message(paste("Calling plot for:", clean_args$Title))
    return(run_with_counter(.Single_Plot_original, args = clean_args, session = session))
  }
}


#
#
#
#
# Single_Plot <- function(..., session = NULL) {
#   dots <- list(...)
#   fn_formals <- formals(.Single_Plot_original)
#   valid_args <- names(fn_formals)
#
#   clean_args <- dots[names(dots) %in% valid_args]
#
#   # Evaluate missing defaults
#   for (arg in setdiff(valid_args, names(clean_args))) {
#     default_val <- fn_formals[[arg]]
#     if (!identical(default_val, quote(expr = ))) {
#       val <- tryCatch(
#         eval(default_val, envir = environment(.Single_Plot_original)),
#         error = function(e) NULL
#       )
#       clean_args[[arg]] <- val
#     }
#   }
#
#   Data <- clean_args$Data
#
#   # Check if multiple data frames were passed
#   is_df <- function(x) is.data.frame(x)
#   if (is.list(Data) && all(vapply(Data, is_df, logical(1)))) {
#     message("Multiple data frames detected. Running Single_Plot on each.")
#
#     plots <- lapply(seq_along(Data), function(i) {
#       df <- Data[[i]]
#       name <- names(Data)[i]
#
#       args_i <- clean_args
#       args_i$Data <- df
#
#       # Auto-title if not specified
#       if (is.null(args_i$Title) && !is.null(name) && nzchar(name)) {
#         args_i$Title <- name
#       }
#
#       if (isTRUE(args_i$Verbose)) {
#         do.call(.Single_Plot_original, args_i)
#       } else {
#         run_with_counter(.Single_Plot_original, args = args_i, session = session)
#       }
#     })
#
#     # Attach names if they exist
#     if (!is.null(names(Data))) {
#       names(plots) <- names(Data)
#     }
#
#     return(invisible(plots))
#   }
#
#   # Single data frame fallback
#   clean_args$session <- session
#   if (isTRUE(clean_args$Verbose)) {
#     return(do.call(.Single_Plot_original, clean_args))
#   } else {
#     return(run_with_counter(.Single_Plot_original, args = clean_args, session = session))
#   }
# }
