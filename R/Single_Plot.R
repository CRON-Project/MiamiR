
#' Single, Annotated Manhattan Plot
#'
#' @param Data These are the summary statistics to be plotted; defaults to Intelligence_Sum_Stats
#' @param Chromosome_Column Manually specify chromosome column; defaults to "CHROM"
#' @param Position_Column Manually specify chromosome position column; defaults to "GENPOS"
#' @param SNP_ID_Column Manually specify SNP ID column; defaults to "ID"
#' @param PValue_Column Manually specify P Value column; defaults to "P"
#' @param Title Title to be displayed; defaults to "Plot1"
#' @param Title_Size Size of title; defaults to 35
#' @param Y_Axis_Text_Size Size of Y axis text; defaults to 30
#' @param Y_Axis_Title_Size Size of Y axis title; defaults to 35
#' @param X_Axis_Title_Size Size of X axis title; defaults to 35
#' @param X_Axis_Title X axis title; defaults to "Chromosome"
#' @param Chromosome_Label_Size Size of chromosome numbers; defaults to 30
#' @param Chromosome_Colours Alternating list of two colours for chromosome points (odd/even); defaults to c("blue", "turquoise")
#' @param Sig_Line_Colour Colour of significance threshold; defaults to "red"
#' @param Sig_Line_Type Type of significance threshold line; defaults to "dashed"
#' @param Sig_Threshold Value of significance threshold; defaults to 5e-8
#' @param Sig_Line_Width Width of significance threshold line; defaults to 0.5
#' @param Point_Size Size of points; defaults to 2.5
#' @param Label_Index Label the index SNPs with ID provided; defaults to T
#' @param Label_Size Size of labels used; defaults to 6
#' @param Label_Angle Angle of labels used; defaults to 30
#' @param Label_Colour Colour of index SNP label; defaults to "black"
#' @param Colour_Of_Index Colour of index SNP highlight; defaults to "red"
#' @param Colour_Index Highlight the index SNP; defaults to T
#' @param Condense_Scale Condense scale; defaults to T
#' @param Break_Point Value where y axis condenses; defaults to 1e-10
#' @param File_Name File name to save plot as; defaults to "Manhattan_Plot"
#' @param Width Width of saved plot; defaults to 30
#' @param Height Height of saved plot; defaults to 15
#' @param Quality Quality of saved plot (dpi); defaults to 900
#' @param File_Type File type of saved plot; defaults to "jpg"
#' @param Random_Selection Manually select a random subset of SNPs to plot, given that Draft_Plot is TRUE; defaults to 10000
#' @param Draft_Plot Specify whether a subset of SNPs, specified by Random_Selection is plotted; defaults to F
#' @param Chromosome_Labels Manually specify which chromosome index SNPs should be labelled on the plot; defaults to  c(2,6,12,14,16,20)
#' @param Chromosome_Index Manually specify which chromosome index SNPs should be circled/highlighted on the plot; defaults to  c(2,6,12,14,16,20)
#'
#' @return Image of Single Manhattan Plot is saved to the current directory and ggplot object is saved
#' @export
#'
#' @examples  Manhattan_Plot <- Single_Plot(Data = Intelligence_Sum_Stats,
#'                              Chromosome_Column = "CHROM",
#'                              Position_Column = "GENPOS", SNP_ID_Column = "ID",
#'                              PValue_Column = "P",
#'                              Title = "Intelligence",
#'                              Chromosome_Colours = c("blue", "turquoise"),
#'                              File_Name = "Manhattan_Plot")
#'

Single_Plot <- function(Data = Intelligence_Sum_Stats,
                       Chromosome_Column = "CHROM",
                       Random_Selection = 10000, Draft_Plot = FALSE,
                       Position_Column = "GENPOS", SNP_ID_Column = "ID",
                       PValue_Column = "P", Chromosome_Labels = c(2,6,12,14,16,20),
                       Chromosome_Index = c(2,6,12,14,16,20),
                       Title = "Plot1",
                       Title_Size = 35, Y_Axis_Text_Size = 30,  Y_Axis_Title_Size = 35,
                       X_Axis_Title_Size = 35,
                       X_Axis_Title = "Chromosome",
                       Chromosome_Label_Size = 30,
                       Chromosome_Colours = c("blue", "turquoise"),
                       Sig_Line_Colour = "red",
                       Sig_Line_Type = "dashed",
                       Sig_Threshold = 5e-8,
                       Sig_Line_Width = 0.5,
                       Point_Size = 2.5,
                       Label_Index = TRUE,
                       Label_Size = 6, Label_Angle = 30,
                       Label_Colour = "black",
                       Colour_Of_Index = "red",
                       Colour_Index = TRUE,
                       Condense_Scale = TRUE,
                       Break_Point = 1e-10,
                       File_Name = "Manhattan_Plot", Width = 30, Height = 15, Quality = 600,
                       File_Type = "jpg"
)

{

  if(Draft_Plot == T)
  {

    Data <- dplyr::sample_n(Data, Random_Selection)

  }


  allowed_names_chromosomes <- c("chromosome", "chrom", "chr", "CHROM", "Chromosome", "CHR", "Chr", "Chrom")

  for (allowed_name in allowed_names_chromosomes) {
    if (allowed_name %in% colnames(Data)) {
      usable_chrom_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as chromosome column"))
  Chromosome_Column <- allowed_name

  allowed_names_pvals <- c("P", "p", "Pvalue", "pvalue", "P-Value", "p-value", "p-Value",
                           "P-VALUE", "logp","LogP", "LOGP", "Logp", "log10p","Log10P",
                           "LOG10P", "Log10p",
                           "log10p", "LOG10P", "-LOG10P", "")

  for (allowed_name in allowed_names_pvals) {
    if (allowed_name %in% colnames(Data)) {
      usable_p_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as P-Val column"))
  PValue_Column <- allowed_name


  if (!("P" %in% colnames(Data)) & any(colnames(Data) %in% c("logp", "LogP", "LOGP", "Logp",
                                                                     "log10p", "Log10P", "LOG10P",
                                                                     "Log10p", "-LOG10P"))) {

    print("No P Value in first dataset, Calculating from LOG10P column detected")
    Data$P <- 10^-(as.numeric(Data[[PValue_Column]]))
    PValue_Column <- "P"

  }



  allowed_names_pos <- c("POS", "pos", "Pos", "BP", "BPos", "bpos", "BPOS", "bPos",
                         "Position", "position", "POSITION", "genpos", "GENPOS",
                         "Genpos")


  for (allowed_name in allowed_names_pos) {
    if (allowed_name %in% colnames(Data)) {
      usable_pos_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as position column"))
  Position_Column <- allowed_name


  allowed_names_SNP <- c("ID", "Id", "ID", "RsID", "RsId","RSID", "snp", "SNP", "Snp",
                         "snv" ,"SNV" , "Snv", "RS", "rs")

  for (allowed_name in allowed_names_SNP) {
    if (allowed_name %in% colnames(Data)) {
      usable_snp_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as SNP column"))
  SNP_ID_Column <- allowed_name



  #Manually assign columns for ease of use
  Data$CHROM <- Data[[Chromosome_Column]]
  Data$GENPOS <- Data[[Position_Column]]
  Data$ID <- Data[[SNP_ID_Column]]
  Data$P <- Data[[PValue_Column]]


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

  Title <- paste0(Title, "\n\n")

  X_Axis_Title <- paste0(X_Axis_Title, "\n\n")

  #Basic Plot

  print(Chromosome_Colours)
  print(Chromosome_Index)
  print(Chromosome_Labels)

  print("Plotting Top Framework")

  print(is.numeric(Break_Point))


  a <- ggmanh::manhattan_plot(x = Data, preserve.position = T, plot.title = Title,
                      chr.colname = Chromosome_Column, pos.colname = Position_Column, label.colname = NULL,
                      pval.colname = PValue_Column, annotateTop = FALSE,
                      chr.order = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"),
                      chr.col = Chromosome_Colours, chrlabs = c(1:22, "X"), rescale = Condense_Scale,
                      signif = Break_Point, rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
                      signif.col = c("transparent"),  point.size = Point_Size, x.label = "")

  #Significance Line

  print("Adding Significance Line")

  b <- a + ggplot2::geom_hline(yintercept = -(log10(Sig_Threshold)), linetype = Sig_Line_Type,
                      color = Sig_Line_Colour, size = Sig_Line_Width)

  #Theme/formatting


  print("Formatting")

  p <- b + ggplot2::theme(panel.border =  ggplot2::element_blank(),
                 panel.grid.major =  ggplot2::element_blank(),
                 panel.grid.minor =  ggplot2::element_blank(),
                 panel.background =  ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(hjust = 0.5, size = Title_Size),
                 axis.text.x = ggplot2::element_text(size = Chromosome_Label_Size, vjust = -0.5, colour = "black"),
                 axis.line =  ggplot2::element_line(),
                 axis.title.x =  ggplot2::element_text(size = X_Axis_Title_Size, vjust = -7 ),
                 axis.title.y =  ggplot2::element_text(size = Y_Axis_Title_Size, vjust = 40),
                 axis.text.y =  ggplot2::element_text(size = Y_Axis_Text_Size, colour = "black"),
                 axis.ticks.length.x  =  ggplot2::unit(0.8,"cm"),
                 legend.position = "none",
                 plot.margin =  ggplot2::margin(t = 20,  # Top margin
                                      r = 40,  # Right margin
                                      b = 60,  # Bottom margin
                                      l = 35)) #Left margin


  p <- p +  ggplot2::xlab(X_Axis_Title) + ggplot2::theme(axis.title.x =  ggplot2::element_text(vjust=-5, size = 35, hjust = 0.5))


  #Labelling




  if(Label_Index == TRUE)
  {
    print("Adding Labels")
    p <- p +  ggtext::geom_richtext( ggplot2::aes(label = Lab), hjust = 0, vjust = 0, label.size = 0,
                                    label.padding = grid::unit(rep(6, 4), "pt"), nudge_y = 0, fill = NA,
                                    nudge_x = 0, color = Label_Colour, label.color = NA, size = Label_Size, angle = Label_Angle)

  }



  if(Colour_Index == TRUE)
  {

    print("Highlighting Index SNPs")
    p <- p +   ggplot2::geom_point( ggplot2::aes(size=COLOUR),
                          pch=21, fill=NA,  colour=Colour_Of_Index, stroke=1)

  }

  p <- p +  ggplot2::coord_cartesian(clip = "off")


  print("Scaling Axes")



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





  d <- p +  ggplot2::scale_x_continuous(breaks = middle_new_pos_values,
                              # Midpoints or specific points for chromosome labels
                              labels = factor(c(1:20, "", "22", "X")),

                              position = "bottom",
                              expand =  ggplot2::expansion(mult = c(0.01, 0.01)) )




  Plot_Outcome <- d


  print("Plot Complete")

  Overall_Name <- paste0(File_Name, ".", File_Type)

  print("Saving Single Manhattan Plot")



   ggplot2::ggsave(Overall_Name, plot = Plot_Outcome, width = Width,
          height = Height, units = "in", dpi = Quality)



  return(Plot_Outcome)
}
