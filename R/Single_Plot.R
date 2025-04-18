
#' Single, Annotated Manhattan Plot
#'
#' @param Data These are the summary statistics to be plotted; defaults to NULL
#' @param Chromosome_Column Manually specify chromosome column; defaults to NULL
#' @param Position_Column Manually specify chromosome position column; defaults to NULL
#' @param SNP_ID_Column Manually specify SNP ID column; defaults to NULL
#' @param PValue_Column Manually specify P Value column; defaults to NULL
#' @param Reference_Allele_Column These are  manual reference allele column names for the data set
#' @param Effect_Allele_Column These are  manual effect allele column names for the data set
#' @param Title Title to be displayed; defaults to NULL
#' @param Title_Size Size of title; defaults to 35
#' @param Y_Axis_Text_Size Size of Y axis text; defaults to 30
#' @param Y_Axis_Title_Size Size of Y axis title; defaults to 35
#' @param X_Axis_Title_Size Size of X axis title; defaults to 35
#' @param X_Axis_Title X axis title; defaults to "Chromosome"
#' @param X_Axis_Title Y axis title; defaults to "-log₁₀P"
#' @param Chromosome_Label_Size Size of chromosome numbers; defaults to 30
#' @param Chromosome_Colours Alternating list of two colours for chromosome points (odd/even); defaults to c("blue", "turquoise")
#' @param Sig_Line_Colour Colour of significance threshold; defaults to "red"
#' @param Sig_Line_Type Type of significance threshold line; defaults to "dashed"
#' @param Sig_Threshold Value of significance threshold; defaults to 5e-8
#' @param Sig_Line_Width Width of significance threshold line; defaults to 0.5
#' @param Point_Size Size of points; defaults to 2.5
#' @param Label_Index Label the index SNPs with ID provided; defaults to T
#' @param Label_Size Size of labels used; defaults to 5
#' @param Label_Angle Angle of labels used; defaults to 45
#' @param Label_Colour Colour of index SNP label; defaults to "black"
#' @param Colour_Of_Index Colour of index SNP highlight; defaults to "red"
#' @param Colour_Index Highlight the index SNP; defaults to T
#' @param Condense_Scale Condense scale; defaults to T
#' @param Break_Point Value where y axis condenses; defaults to 1e-10
#' @param Random_Selection Manually select a random subset of SNPs to plot, given that Draft_Plot is TRUE; defaults to 10000
#' @param Draft_Plot Specify whether a subset of SNPs, specified by Random_Selection is plotted; defaults to F
#' @param Chromosome_Labels Manually specify which chromosome index SNPs should be labelled on the plot; defaults to  c(1:22, "X")
#' @param Chromosome_Index Manually specify which chromosome index SNPs should be circled/highlighted on the plot; defaults to c(1:22, "X")
#' @param Chromosome_Label_Drops Manually specify which chromosome/X axis labels should not be shown; defaults to c(21,22)
#'
#' @return Image of Single Manhattan Plot is saved to specified object and ggplot object can then be saved
#' @export
#'
#' @examples  Manhattan_Plot <- Single_Plot(Data = Intelligence_Sum_Stats)
#'

Single_Plot<- function(Data = NULL,
                       Random_Selection = 10000,
                       Draft_Plot = FALSE,
                       Chromosome_Labels = c(1:22, "X"),
                       Chromosome_Index = c(1:22, "X"),
                       Title = NULL,
                       Title_Size = 35,
                       Y_Axis_Text_Size = 30,
                       Y_Axis_Title_Size = 35,
                       X_Axis_Title_Size = 35,
                       Chromosome_Label_Size = 30,
                       Chromosome_Label_Drops = c(21,22),  # can be numeric or character
                       X_Axis_Title = "Chromosome",
                       Y_Axis_Title = "-log₁₀(P)",
                       Chromosome_Colours = c("blue", "turquoise"),
                       Sig_Line_Colour = "red",
                       Sig_Line_Type = "dashed",
                       Sig_Threshold = 5e-8,
                       Sig_Line_Width = 0.5,
                       Point_Size = 2.5,
                       Label_Index = TRUE,
                       Label_Size = 5, Label_Angle = 45,
                       Label_Colour = "black",
                       Colour_Of_Index = "red",
                       Colour_Index = TRUE,
                       Condense_Scale = TRUE,
                       Break_Point = 1e-10,
                       Chromosome_Column = NULL,
                       Position_Column = NULL, SNP_ID_Column = NULL,
                       PValue_Column = NULL,
                       Reference_Allele_Column = NULL,  Effect_Allele_Column = NULL,
                       Interactive = TRUE)

{





  if (is.null(Data)) {
    stop("Please provide data.", call. = FALSE)
  }


  if (is.null(Title)) {
    if (is.character(Data) && length(Data) == 1 && file.exists(Data)) {
      # Extract file name without path or extension
      Title <- tools::file_path_sans_ext(basename(Data))
    } else {
      # Fallback to variable name if Data is not a path
      Title <- deparse(substitute(Data))
    }
  }


  #print(Title)




  if (is.character(Data) && length(Data) == 1) {
    file_path <- Data

   # print("Reading")

    if (file.exists(file_path)) {
      message("Reading data from file: ", file_path)

      Data <- tryCatch(
        vroom::vroom(file_path, show_col_types = FALSE),
        error = function(e) stop("Failed to read file: ", conditionMessage(e), call. = FALSE)
      )
    } else {
      stop("The provided character string does not point to an existing file: ", file_path, call. = FALSE)
    }
  }

#  print(Data)


  if(Draft_Plot == T)
  {

    Data <- dplyr::sample_n(Data, Random_Selection)

  }



  Chromosome_Column <- detect_chromosome_column(Data, Chromosome_Column)


  #Need this for ggmanh conformation
#  Data$CHR[Data$CHR == 23] <- "X"
  Data[[Chromosome_Column]][as.character(Data[[Chromosome_Column]]) == "23"] <- "X"


  PValue_Column     <- detect_pvalue_column(Data, PValue_Column)
  Position_Column   <- detect_position_column(Data, Position_Column)
  SNP_ID_Column     <- detect_snp_column(Data, SNP_ID_Column)
  Ref_Allele_Column     <- detect_reference_allele_column(Data, Reference_Allele_Column)
  Alt_Allele_Column     <- detect_effect_allele_column(Data, Effect_Allele_Column)

#  print(Chromosome_Column)


  #Manually assign columns for ease of use
  Data$CHROM <- Data[[Chromosome_Column]]
  Data$GENPOS <- Data[[Position_Column]]
  Data$ID <- Data[[SNP_ID_Column]]
  Data$P <- Data[[PValue_Column]]
  Data$ALLELE0 <- Data[[Ref_Allele_Column]]
  Data$ALLELE1 <- Data[[Alt_Allele_Column]]

  Data$ALLELE0 <- toupper(Data$ALLELE0 )
  Data$ALLELE1 <- toupper(Data$ALLELE1 )


  if(("Lab" %in% colnames(Data))) # if annotate used before then skip
  {
    print("It looks like you've applied Annotate_Data()")
  }

  if(!("Lab" %in% colnames(Data))) # if annotate used before then skip
  {

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

  Data$Lab[Data$Lab == ""] <- NA
  Data$COLOUR[is.na(Data$Lab)] <- NA
  Data$COLOUR[!is.na(Data$Lab)] <- 4



  Data$Lab[!(Data$CHROM %in% Chromosome_Labels)] <- NA
  Data$COLOUR[!(Data$CHROM %in% Chromosome_Index)] <- NA

  Title <- paste0(Title, "\n \n ") # plotly only interprets for y if actual content on line in app

  X_Axis_Title <- paste0("\n ", X_Axis_Title) # needs to go above it to push down opp to above.
  Y_Axis_Title <- paste0(Y_Axis_Title, "\n ") # needs to go left above it to push down opp to above.


  #Basic Plot

#  print(Chromosome_Colours)
#  print(Chromosome_Index)
#  print(Chromosome_Labels)

  message("Plotting Top Framework")

 # print(is.numeric(Break_Point))


  if(Interactive == TRUE)

  {
#
#   Data$Hover_Info <- paste0(
#     "SNP: ", Data$ID, "\n",
#     "CHR: ", Data$CHROM, "\n",
#     "POS: ", Data$GENPOS, "\n",
#     "P: ", signif(Data$P, 4)
#   )
#
    Data$Hover_Info <- paste0(
      "SNP: ", Data$ID, "\n",
      "CHR: ", Data$CHROM, "\n",
      "POS: ", Data$GENPOS, "\n",
      "P: ", signif(Data$P, 4), "\n",
      "REF: ", Data$ALLELE0, " ALT: ", Data$ALLELE1
    )



}

  a <- ggmanh::manhattan_plot(x = Data, preserve.position = T, plot.title = Title,
                      chr.colname = Chromosome_Column, pos.colname = Position_Column, label.colname = NULL,
                      pval.colname = PValue_Column, annotateTop = FALSE,
                      chr.order = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"),
                      chr.col = Chromosome_Colours, chrlabs = c(1:22, "X"), rescale = Condense_Scale,
                      signif = Break_Point, rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
                      signif.col = c("transparent"),  point.size = Point_Size, x.label = "", y.label = "")

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

 # return(a)




  message("Adding Significance Line")

  b <- a + ggplot2::geom_hline(yintercept = -(log10(Sig_Threshold)), linetype = Sig_Line_Type,
                      color = Sig_Line_Colour, size = Sig_Line_Width)

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


  p <- p +  ggplot2::xlab(X_Axis_Title) + ggplot2::theme(axis.title.x =  ggplot2::element_text(vjust=-5, size = X_Axis_Title_Size, hjust = 0.5))
  p <- p +  ggplot2::ylab(Y_Axis_Title) + ggplot2::theme(axis.title.y =  ggplot2::element_text(vjust=5, size = Y_Axis_Title_Size))
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

  }
  #La belling interact

  if(Interactive == TRUE)

  {
  p <- p + ggplot2::aes(text = Hover_Info)
}



  if(Label_Index == TRUE)
{
    message("Adding Labels")
    p <- p +  ggtext::geom_richtext( ggplot2::aes(label = Lab), hjust = 0, vjust = 0, label.size = 0,
                                    label.padding = grid::unit(rep(6, 4), "pt"), nudge_y = 0, fill = NA,
                                    nudge_x = 0, color = Label_Colour, label.color = NA, size = Label_Size, angle = Label_Angle)

  }



  suppressMessages(suppressWarnings({


  if(Colour_Index == TRUE & Interactive == FALSE)
  {


  message("Highlighting Index SNPs")
   p <- p +   ggplot2::geom_point( ggplot2::aes(size=COLOUR),
                        pch=21, fill=NA,  colour=Colour_Of_Index, stroke=1)

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

  Plot_Outcome <- d


  # Plot_Outcome <- Plot_Outcome +  theme(
  #   plot.margin = margin(t = 60, r = 30, b = 30, l = 30, unit = "pt")
  # )

  message("Plot Complete")

#  Overall_Name <- paste0(File_Name, ".", File_Type)

#  print("Saving Single Manhattan Plot")


#
#    ggplot2::ggsave(Overall_Name, plot = Plot_Outcome, width = Width,
#           height = Height, units = "in", dpi = Quality)
#



  return(Plot_Outcome)
}
