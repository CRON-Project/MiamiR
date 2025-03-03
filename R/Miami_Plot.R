#' Annotated Miami Plot
#'
#' @param Top_Data These are the summary statistics to be plotted on the top panel; defaults to Intelligence_Sum_Stats
#' @param Bottom_Data These are the summary statistics to be plotted on the bottom panel; defaults to Household_Income_Sum_Stats
#' @param Top_Title Title to be displayed on the top panel; defaults to "Plot1"
#' @param Bottom_Title Title to be displayed on the bottom panel; defaults to "Plot2"
#' @param Bottom_Title_Size Size of bottom panel title; defaults to 35
#' @param Top_Title_Size Size of top panel title; defaults to 35
#' @param Y_Axis_Text_Size Size of Y axes text; defaults to 30
#' @param Y_Axis_Title_Size Size of Y axes titles; defaults to 35
#' @param Chromosome_Label_Size Size of chromosome numbers; defaults to 30
#' @param Top_Sig_Line_Colour Colour of significance threshold on top panel; defaults to "red"
#' @param Bottom_Sig_Line_Colour Colour of significance threshold on bottom panel; defaults to "red"
#' @param Top_Sig_Line_Type Type of significance threshold line on top panel; defaults to "dashed"
#' @param Bottom_Sig_Line_Type Type of significance threshold line on bottom panel; defaults to "dashed"
#' @param Top_Sig_Threshold Value of significance threshold on top panel; defaults to 5e-8
#' @param Bottom_Sig_Threshold Value of significance threshold on bottom panel; defaults to 5e-8
#' @param Top_Sig_Line_Width Width of significance threshold line on top panel; defaults to 0.5
#' @param Bottom_Sig_Line_Width Width of significance threshold line on bottom panel; defaults to 0.5
#' @param Point_Size Size of points; defaults to 2.5
#' @param Top_Label_Index Label the index SNPs with ID provided on top panel; defaults to T
#' @param Bottom_Label_Index Label the index SNPs with ID provided on bottom panel; defaults to T
#' @param Label_Size Size of labels used; defaults to 6
#' @param Label_Angle Angle of labels used; defaults to 30
#' @param Top_Label_Colour Colour of index SNP label on top panel; defaults to "black"
#' @param Bottom_Label_Colour Colour of index SNP label on bottom panel; defaults to "black"
#' @param Top_Colour_Of_Index Colour of index SNP highlight on top panel; defaults to "red"
#' @param Bottom_Colour_Of_Index Colour of index SNP highlight on bottom panel; defaults to "red"
#' @param Top_Colour_Index Highlight the index SNP on top panel; defaults to T
#' @param Bottom_Colour_Index Highlight the index SNP on bottom panel; defaults to T
#' @param Top_Condense_Scale Condense scale on top panel; defaults to T
#' @param Bottom_Condense_Scale Condense scale on bottom panel; defaults to T
#' @param Top_Break_Point Value where y axis condenses on top panel; defaults to "red"
#' @param Bottom_Break_Point Value where y axis condenses on bottom panel; defaults to "red"
#' @param File_Name File name to save plot as; defaults to "Miami_Plot"
#' @param Width Width of saved plot; defaults to 30
#' @param Height Height of saved plot; defaults to 15
#' @param Quality Quality of saved plot (dpi); defaults to 900
#' @param File_Type File type of saved plot; defaults to "jpg"
#' @param Top_Chromosome_Column Manually specify chromosome column for top panel; defaults to "CHROM"
#' @param Top_Position_Column Manually specify chromosome position column for top panel; defaults to "GENPOS"
#' @param Top_SNP_ID_Column Manually specify SNP ID column for top panel; defaults to "ID"
#' @param Top_PValue_Column Manually specify P Value column for top panel; defaults to "P"
#' @param Bottom_Chromosome_Column Manually specify chromosome column for bottom panel; defaults to "CHROM"
#' @param Bottom_Position_Column Manually specify chromosome position column for bottom panel; defaults to "GENPOS"
#' @param Bottom_SNP_ID_Column Manually specify SNP ID column for bottom panel; defaults to "ID"
#' @param Bottom_PValue_Column Manually specify P Value column for bottom panel; defaults to "P"
#' @param Random_Selection Manually select a random subset of SNPs to plot, given that Draft_Plot is TRUE; defaults to 10000
#' @param Draft_Plot Specify whether a subset of SNPs, specified by Random_Selection is plotted; defaults to F
#' @param Top_Colour_One The first colour (odd chromosome numbers) to use on the top plot; defaults to "green"
#' @param Top_Colour_Two The second colour (even chromosome numbers) to use on the top plot; defaults to "purple"
#' @param Bottom_Colour_One  The first colour (odd chromosome numbers) to use on the bottom plot; defaults to "green"
#' @param Bottom_Colour_Two The second colour (even chromosome numbers) to use on the bottom plot; defaults to "purple"
#' @param Top_Chromosome_Labels Manually specify which chromosome index SNPs should be labelled on the top plot; defaults to  c(2,6,12,14,16,20)
#' @param Top_Chromosome_Index Manually specify which chromosome index SNPs should be circled/highlighted on the top plot; defaults to  c(2,6,12,14,16,20)
#' @param Bottom_Chromosome_Labels Manually specify which chromosome index SNPs should be labelled on the bottom plot; defaults to  c(2,6,12,14,16,20)
#' @param Bottom_Chromosome_Index Manually specify which chromosome index SNPs should be circled/highlighted on the bottom plot; defaults to  c(2,6,12,14,16,20)
#'
#' @return Image of Miami Plot is saved to the current directory and ggplot object is saved for top and bottom plots in addition to a combined grob
#' @export
#'
#' @examples Miami_Plot <- Miami_Plot(Top_Data = Intelligence_Sum_Stats, Bottom_Data = Household_Income_Sum_Stats,
#'                         Top_Chromosome_Column = "CHROM",
#'                         Top_Position_Column = "GENPOS", Top_SNP_ID_Column = "ID",
#'                         Top_PValue_Column = "P",
#'                         Top_Colour_One = "green", Top_Colour_Two = "purple",
#'                         Bottom_Colour_One = "green", Bottom_Colour_Two = "purple",
#'                         Bottom_Chromosome_Column = "CHROM",
#'                         Bottom_Position_Column = "GENPOS", Bottom_SNP_ID_Column = "ID",
#'                         Bottom_PValue_Column = "P",
#'                         Top_Title = "Intelligence", Bottom_Title = "Household Income",
#'                         File_Name = "Miami_Plot")
#'
#'

Miami_Plot <- function(Top_Data = Intelligence_Sum_Stats, Bottom_Data = Household_Income_Sum_Stats,
                       Random_Selection = 10000, Draft_Plot = FALSE,
                       Top_Chromosome_Column = "CHROM",
                       Top_Position_Column = "GENPOS", Top_SNP_ID_Column = "ID",
                       Top_PValue_Column = "P",
                       Top_Colour_One = "green", Top_Colour_Two = "purple",
                       Bottom_Colour_One = "green", Bottom_Colour_Two = "purple",
                       Bottom_Chromosome_Column = "CHROM",
                       Bottom_Position_Column = "GENPOS", Bottom_SNP_ID_Column = "ID",
                       Bottom_PValue_Column = "P",
                       Top_Title = "Plot1", Bottom_Title = "Plot2",
                       Bottom_Title_Size = 35,
                       Top_Title_Size = 35, Y_Axis_Text_Size = 30,  Y_Axis_Title_Size = 35,
                       Chromosome_Label_Size = 30,
                       Top_Sig_Line_Colour = "red", Bottom_Sig_Line_Colour = "red",
                       Top_Sig_Line_Type = "dashed", Bottom_Sig_Line_Type = "dashed",
                       Top_Sig_Threshold = 5e-8, Bottom_Sig_Threshold = 5e-8,
                       Top_Sig_Line_Width = 0.5, Bottom_Sig_Line_Width = 0.5,
                       Point_Size = 2.5,
                       Top_Label_Index = TRUE, Bottom_Label_Index = TRUE,
                       Label_Size = 6, Label_Angle = 30,
                       Top_Label_Colour = "black",   Bottom_Label_Colour = "black",
                       Top_Colour_Of_Index = "red",
                       Top_Chromosome_Labels = c(2,6,12,14,16,20),
                       Top_Chromosome_Index = c(2,6,12,14,16,20),
                       Bottom_Chromosome_Labels = c(2,6,12,14,16,20),
                       Bottom_Chromosome_Index = c(2,6,12,14,16,20),
                       Bottom_Colour_Of_Index = "red",
                       Top_Colour_Index = TRUE,  Bottom_Colour_Index = TRUE,
                       Top_Condense_Scale = TRUE, Bottom_Condense_Scale = TRUE,
                       Top_Break_Point = 1e-10, Bottom_Break_Point = 1e-10,
                       File_Name = "Miami_Plot", Width = 30, Height = 15, Quality = 600,
                       File_Type = "jpg"
)

{

  #add file loading later like others.

  if(Draft_Plot == T)
  {
  Top_Data <- dplyr::sample_n(Top_Data, Random_Selection)
  Bottom_Data <- dplyr::sample_n(Bottom_Data, Random_Selection)
  }

  allowed_names_chromosomes <- c("chromosome", "chrom", "chr", "CHROM", "Chromosome", "CHR", "Chr", "Chrom")

  for (allowed_name in allowed_names_chromosomes) {
    if (allowed_name %in% colnames(Top_Data)) {
      usable_chrom_top <- colnames(Top_Data)[which(colnames(Top_Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as chromosome column"))
  Top_Chromosome_Column <- allowed_name

  allowed_names_pvals <- c("P", "p", "Pvalue", "pvalue", "P-Value", "p-value", "p-Value",
                           "P-VALUE", "logp","LogP", "LOGP", "Logp", "log10p","Log10P",
                           "LOG10P", "Log10p",
                           "log10p", "LOG10P", "-LOG10P", "")

  for (allowed_name in allowed_names_pvals) {
    if (allowed_name %in% colnames(Top_Data)) {
      usable_p_top <- colnames(Top_Data)[which(colnames(Top_Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as P-Val column"))
  Top_PValue_Column <- allowed_name


  if (!("P" %in% colnames(Top_Data)) & any(colnames(Top_Data) %in% c("logp", "LogP", "LOGP", "Logp",
                                                                     "log10p", "Log10P", "LOG10P",
                                                                     "Log10p", "-LOG10P"))) {

    print("No P Value in first dataset, Calculating from LOG10P column detected")
    Top_Data$P <- 10^-(as.numeric(Top_Data[[Top_PValue_Column]]))
    Top_PValue_Column <- "P"

  }


  allowed_names_pos <- c("POS", "pos", "Pos", "BP", "BPos", "bpos", "BPOS", "bPos",
                         "Position", "position", "POSITION", "genpos", "GENPOS",
                         "Genpos")

  for (allowed_name in allowed_names_pos) {
    if (allowed_name %in% colnames(Top_Data)) {
      usable_pos_top <- colnames(Top_Data)[which(colnames(Top_Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as position column"))
  Top_Position_Column <- allowed_name


  allowed_names_SNP <- c("ID", "Id", "ID", "RsID", "RsId","RSID", "snp", "SNP", "Snp",
                         "snv" ,"SNV" , "Snv", "RS", "rs")

  for (allowed_name in allowed_names_SNP) {
    if (allowed_name %in% colnames(Top_Data)) {
      usable_snp_top <- colnames(Top_Data)[which(colnames(Top_Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as SNP column"))
  Top_SNP_ID_Column <- allowed_name



  #Manually assign columns for ease of use
  Top_Data$CHROM <- Top_Data[[Top_Chromosome_Column]]
  Top_Data$GENPOS <- Top_Data[[Top_Position_Column]]
  Top_Data$ID <- Top_Data[[Top_SNP_ID_Column]]
  Top_Data$P <- Top_Data[[Top_PValue_Column ]]


  if(("Lab" %in% colnames(Top_Data))) # if annotate used before then skip
  {
    print("It looks like you've applied Annotate_Data()")
  }

  if(!("Lab" %in% colnames(Top_Data))) # if annotate used before then skip
  {


  print("It looks like you haven't applied Annotate_Data()")

  Top_Data <- Top_Data %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(min_P_GENPOS = GENPOS[which.min(P)]) %>% #P is formed in function
    dplyr::ungroup()



  # Label the minimum P value for each CHROM
  Top_Data <- Top_Data %>%
    dplyr::mutate(Lab = ifelse(GENPOS == min_P_GENPOS & P < 5e-8, ID, # only sigs
                        ifelse(abs(GENPOS - min_P_GENPOS) > 2000000000 & P == min(P[abs(GENPOS - min_P_GENPOS) > 2000000000]),
                               ID, "")))

  }

  #Have to do this to scale searate builds and SNPs run.

  Top_Data$Lab[Top_Data$Lab == ""] <- NA
  Top_Data$COLOUR[is.na(Top_Data$Lab)] <- NA
  Top_Data$COLOUR[!is.na(Top_Data$Lab)] <- 5.5

  Top_Data$Lab[!(Top_Data$CHROM %in% Top_Chromosome_Labels)] <- NA
  Top_Data$COLOUR[!(Top_Data$CHROM %in% Top_Chromosome_Index)] <- NA

  Top_Data$COLOUR2[Top_Data$P == 1] <- 0 # random seed everytime if using that function
  Top_Data$COLOUR2[Top_Data$P != 1 & !(Top_Data$CHROM %in% c(1,3,5,7,9,11,13,15,17,19,21,23))] <- Point_Size
  Top_Data$COLOUR3[Top_Data$P == 1] <- 0 # random seed everytime if using that function
  Top_Data$COLOUR3[Top_Data$P != 1 & (Top_Data$CHROM %in% c(1,3,5,7,9,11,13,15,17,19,21,23))] <- Point_Size


 # return(Top_Data)


  Top_Title <- paste0(Top_Title, "\n")
  Bottom_Title <- paste0( "\n", Bottom_Title)

  #Basic Plot


  allowed_names_chromosomes <- c("chromosome", "chrom", "chr", "CHROM", "Chromosome", "CHR", "Chr", "Chrom")

  for (allowed_name in allowed_names_chromosomes) {
    if (allowed_name %in% colnames(Bottom_Data)) {
      usable_chrom_top <- colnames(Bottom_Data)[which(colnames(Bottom_Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as chromosome column"))
  Bottom_Chromosome_Column <- allowed_name

  allowed_names_pvals <- c("P", "p", "Pvalue", "pvalue", "P-Value", "p-value", "p-Value",
                           "P-VALUE", "logp","LogP", "LOGP", "Logp", "log10p","Log10P",
                           "LOG10P", "Log10p",
                           "log10p", "LOG10P", "-LOG10P", "")

  for (allowed_name in allowed_names_pvals) {
    if (allowed_name %in% colnames(Bottom_Data)) {
      usable_p_top <- colnames(Bottom_Data)[which(colnames(Bottom_Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as P-Val column"))
  Bottom_PValue_Column <- allowed_name


  if (!("P" %in% colnames(Bottom_Data)) & any(colnames(Bottom_Data) %in% c("logp", "LogP", "LOGP", "Logp",
                                                                     "log10p", "Log10P", "LOG10P",
                                                                     "Log10p", "-LOG10P"))) {

      print("No P Value in second dataset, Calculating from LOG10P column detected")
      Bottom_Data$P <- 10^-(as.numeric(Bottom_Data[[allowed_name]]))
      Bottom_PValue_Column <- "P"


  }


  allowed_names_pos <- c("POS", "pos", "Pos", "BP", "BPos", "bpos", "BPOS", "bPos",
                         "Position", "position", "POSITION", "genpos", "GENPOS",
                         "Genpos")

  for (allowed_name in allowed_names_pos) {
    if (allowed_name %in% colnames(Bottom_Data)) {
      usable_pos_top <- colnames(Bottom_Data)[which(colnames(Bottom_Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as position column"))
  Bottom_Position_Column <- allowed_name


  allowed_names_SNP <- c("ID", "Id", "ID", "RsID", "RsId","RSID", "snp", "SNP", "Snp",
                         "snv" ,"SNV" , "Snv")

  for (allowed_name in allowed_names_SNP) {
    if (allowed_name %in% colnames(Bottom_Data)) {
      usable_snp_top <- colnames(Bottom_Data)[which(colnames(Bottom_Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as position column"))
  Bottom_SNP_ID_Column <- allowed_name



  Bottom_Data$CHROM <- Bottom_Data[[Bottom_Chromosome_Column]]
  Bottom_Data$GENPOS <- Bottom_Data[[Bottom_Position_Column]]
  Bottom_Data$ID <- Bottom_Data[[Bottom_SNP_ID_Column]]
  Bottom_Data$P <- Bottom_Data[[Bottom_PValue_Column]]

  if(("Lab" %in% colnames(Bottom_Data))) # if annotate used before then skip
  {
     print("It looks like you've applied Annotate_Data()")
  }

  if(!("Lab" %in% colnames(Bottom_Data))) # if annotate used before then skip
  {

    print("It looks like you haven't applied Annotate_Data()")

  Bottom_Data <- Bottom_Data %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(min_P_GENPOS = GENPOS[which.min(P)]) %>% #P is formed in function
    dplyr::ungroup()



  # Label the minimum P value for each CHROM


  Bottom_Data <- Bottom_Data %>%
    dplyr::mutate(Lab = ifelse(GENPOS == min_P_GENPOS & P < 5e-8, ID,
                        ifelse(abs(GENPOS - min_P_GENPOS) > 2000000000 & P == min(P[abs(GENPOS - min_P_GENPOS) > 2000000000]),
                               ID, "")))

}

  Bottom_Data$Lab[Bottom_Data$Lab == ""] <- NA
  Bottom_Data$COLOUR[is.na(Bottom_Data$Lab)] <- NA
  Bottom_Data$COLOUR[!is.na(Bottom_Data$Lab)] <- 5.5

  Bottom_Data$Lab[!(Bottom_Data$CHROM %in% Bottom_Chromosome_Labels)] <- NA
  Bottom_Data$COLOUR[!(Bottom_Data$CHROM %in% Bottom_Chromosome_Index)] <- NA

  Bottom_Data$COLOUR2[Bottom_Data$P == 1] <- 0 # random seed everytime if using that function
  Bottom_Data$COLOUR2[Bottom_Data$P != 1 & !(Bottom_Data$CHROM %in% c(1,3,5,7,9,11,13,15,17,19,21,23))] <- Point_Size
  Bottom_Data$COLOUR3[Bottom_Data$P == 1] <- 0 # random seed everytime if using that function
  Bottom_Data$COLOUR3[Bottom_Data$P != 1 & (Bottom_Data$CHROM %in% c(1,3,5,7,9,11,13,15,17,19,21,23))] <- Point_Size




  Bottom_Title <- paste0("\n", Bottom_Title)

  #Basic Plots




  Top_Data <- Top_Data %>%
    dplyr::mutate(GENPOS = as.numeric(GENPOS), CHROM = as.numeric(CHROM))

#  return(Top_Data)

  Bottom_Data <- Bottom_Data %>%
    dplyr::mutate(GENPOS = as.numeric(GENPOS), CHROM = as.numeric(CHROM))

  # Create combined keys in both datasets
  Top_Data <- Top_Data %>%
    dplyr::mutate(combined_key = paste(GENPOS, CHROM, sep = "_"))

  Bottom_Data <- Bottom_Data %>%
    dplyr::mutate(combined_key = paste(GENPOS, CHROM, sep = "_"))

  # Find rows in Bottom_Data that are not in Top_Data
  unique_in_bottom <- Bottom_Data %>%
    dplyr::filter(!(combined_key %in% Top_Data$combined_key))

  # Find rows in Top_Data that are not in Bottom_Data
  unique_in_top <- Top_Data %>%
    dplyr::filter(!(combined_key %in% Bottom_Data$combined_key))

  print(nrow(unique_in_bottom))
  print(nrow(unique_in_top))


#  return(Top_Data)

  # Display results

  cat("Rows in Bottom_Data but not in Top_Data:\n")
  print("Adjusting Scales")

  cat("Rows in Top_Data but not in Bottom_Data:\n")
  print("Adjusting Scales")

  if(nrow(unique_in_bottom) != 0)

  {

  new_rows <- unique_in_bottom %>%
    dplyr::select(GENPOS, CHROM)

  empty_cols <- Top_Data %>%
    dplyr::select(-GENPOS, -CHROM) %>%
    dplyr::summarise_all(~NA) %>%
    dplyr::slice(rep(1, nrow(new_rows)))


  new_data <-  dplyr::bind_cols(new_rows, empty_cols)
  new_data$P <- 1
  new_data[[Top_Position_Column]] <- new_data$GENPOS
  new_data[[Top_Chromosome_Column]] <- new_data$CHROM
  new_data[[Top_PValue_Column]] <- new_data$P

  print("Error after this")

  print(Top_Data)
  print(new_data)

  #In case dfs identical skip this step as nothing to adjust for.

  Top_Data <-  dplyr::bind_rows(Top_Data, new_data)

  }

 # return(Top_Data)

  print("Top Config")

  if(nrow(unique_in_top) != 0)

  {

  new_rows <- unique_in_top %>%
    dplyr::select(GENPOS, CHROM)

  empty_cols <- Bottom_Data %>%
    dplyr::select(-GENPOS, -CHROM) %>%
    dplyr::summarise_all(~NA) %>%
    dplyr::slice(rep(1, nrow(new_rows)))


  new_data <-  dplyr::bind_cols(new_rows, empty_cols)
  new_data$P <- 1
  new_data[[Bottom_Position_Column]] <- new_data$GENPOS
  new_data[[Bottom_Chromosome_Column]] <-  new_data$CHROM
  new_data[[Bottom_PValue_Column]] <- new_data$P



  Bottom_Data <-  dplyr::bind_rows(Bottom_Data, new_data)

  }


  print("Bottom Config")



  print("Plotting Top Framework")

  print(Top_Data)

  print(Top_Title)
  print(Top_Chromosome_Column)
  print(Top_Position_Column)
  print(Top_PValue_Column)
  print(Top_Condense_Scale)
  print(Top_Break_Point)
  print(Point_Size)



  print("CHROMS")
  print(table(Top_Data$CHROM))
  Top_Data$CHROM[Top_Data$CHROM == 23] <- "X"
  print(table(Top_Data$CHROM))

  a <- ggmanh::manhattan_plot(x = Top_Data, preserve.position = T, plot.title = Top_Title,
                              chr.colname = Top_Chromosome_Column, pos.colname = Top_Position_Column, label.colname = NULL,
                              pval.colname = Top_PValue_Column, annotateTop = FALSE,
                              chr.order = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"),
                              chr.col = c("transparent"), chrlabs = c(1:22, "X"), rescale = Top_Condense_Scale,
                              signif = Top_Break_Point, rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
                              signif.col = c("transparent"),  point.size = Point_Size, x.label = "" )

  #Significance Line

  print("Adding Significance Line")

  b <- a + ggplot2::geom_hline(yintercept = -(log10(Top_Sig_Threshold)), linetype = Top_Sig_Line_Type,
                               color = Top_Sig_Line_Colour, size = Top_Sig_Line_Width)

  #Theme/formatting


  print("Formatting")

  p <- b + ggplot2::theme(panel.border =  ggplot2::element_blank(),
                          panel.grid.major =  ggplot2::element_blank(),
                          panel.grid.minor =  ggplot2::element_blank(),
                          panel.background =  ggplot2::element_blank(),
                          plot.title =  ggplot2::element_text(hjust = 0.5, size = Top_Title_Size),
                          axis.text.x =  ggplot2::element_text(size = Chromosome_Label_Size, vjust = -1.1, colour = "black"),
                          axis.line =  ggplot2::element_line(),
                          axis.title.x =  ggplot2::element_blank(),
                          axis.title.y =  ggplot2::element_text(size = Y_Axis_Title_Size, vjust = 2),
                          axis.text.y =  ggplot2::element_text(size = Y_Axis_Text_Size, colour = "black"),
                          axis.ticks.length.x  =  ggplot2::unit(0.8,"cm"),
                          legend.position = "none",
                          plot.margin =  ggplot2::margin(t = 20,  # Top margin
                                               r = 20,  # Right margin
                                               b = 0,  # Bottom margin
                                               l = 35)) #Left margin


  #Labelling



  if(Top_Label_Index == TRUE)
  {
    print("Adding Labels")
    p <- p +  ggtext::geom_richtext( ggplot2::aes(label = Lab), hjust = 0, vjust = 0, label.size = 0,
                                    label.padding = grid::unit(rep(6, 4), "pt"), nudge_y = 0, fill = NA,
                                    nudge_x = 0, color = Top_Label_Colour, label.color = NA, size = Label_Size, angle = Label_Angle)

  }



  if(Top_Colour_Index == TRUE)
  {

    print("Highlighting Index SNPs")
    p <- p +   ggplot2::geom_point( ggplot2::aes(size=COLOUR),
                                   pch=21, fill=NA,  colour=Top_Colour_Of_Index, stroke=1)


    p <- p + ggplot2::scale_size_identity()
  }

  print("Filling Background Noise")


  p <- p +   ggplot2::geom_point( ggplot2::aes(size=COLOUR2),
                                  fill=Top_Colour_One,  colour=Top_Colour_One)
  p <- p +   ggplot2::geom_point( ggplot2::aes(size=COLOUR3),
                                  fill=Top_Colour_Two,  colour=Top_Colour_Two)

  p <- p + ggplot2::scale_size_identity()

  p <- p +  ggplot2::coord_cartesian(clip = "off")


  print("Scaling Axes")




  df <- as.data.frame(p$data)

  middle_new_pos_values <- numeric()





  for (chrom in c(1:22, "X")) {
    # Filter for the current chromosome

    print(chrom)

     df_filtered <- df %>%
      dplyr::filter(CHROM == chrom)

    print(df_filtered)

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


  print(middle_new_pos_values)

#  return(p)


  d <- p +  ggplot2::scale_x_continuous(breaks = middle_new_pos_values,
                              # Midpoints or specific points for chromosome labels
                              labels = factor(c(1:20, "", "22", "X")),

                              position = "bottom",
                              expand =  ggplot2::expansion(mult = c(0.01, 0.01)) )



  #return(d)

  Top_Plot_Outcome <- d

 # return(Top_Plot_Outcome)
  #  Top_Data$CHROM[Top_Data$CHROM == 23] <- "X" #needs to go X later before calc due to ggmanh func
  #return(Top_Plot_Outcome) #also update for 23 files input
#  return(Top_Plot_Outcome)

  print("Top Plot Made...")






  print("Plotting Bottom Framework")

  Bottom_Data$CHROM[Bottom_Data$CHROM == 23] <- "X"

  a1 <- ggmanh::manhattan_plot(x = Bottom_Data, preserve.position = T, plot.title =  ggplot2::ggtitle(""),
                       chr.colname = Bottom_Chromosome_Column, pos.colname = Bottom_Position_Column, label.colname = NULL,
                       pval.colname = Bottom_PValue_Column, annotateTop = FALSE,
                       chr.order = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"),
                       chr.col = c("transparent"), chrlabs = c(1:22, "X"), rescale = Bottom_Condense_Scale,
                       signif = Bottom_Break_Point, rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
                       signif.col = c("transparent"),  point.size = Point_Size, x.label = "" )

  #Significance Line

  print("Adding Significance Line")

  b1 <- a1 + ggplot2::geom_hline(yintercept = -(log10(Bottom_Sig_Threshold)), linetype = Bottom_Sig_Line_Type,
                        color = Bottom_Sig_Line_Colour, size = Bottom_Sig_Line_Width)

  #Theme/formatting

  print("Formatting")

  p1 <- b1 + ggplot2::theme(panel.border =  ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_blank(),
                   panel.grid.minor =  ggplot2::element_blank(),
                   panel.background =  ggplot2::element_blank(),
                   #  plot.title = element_text(hjust = 0.5, size = Bottom_Title_Size),
                   plot.title =  ggplot2::element_blank(),
                   # axis.text.x = element_text(size = X_Axis_Text_Size, vjust = -0.5),
                   axis.text.x =  ggplot2::element_blank(),
                   #   axis.title.x = element_blank(),
                   axis.title.x.top =  ggplot2::element_blank(),
                   axis.title.x.bottom =  ggplot2::element_text(size = Bottom_Title_Size, hjust = 0.5),
                   #  axis.title.x.bottom = element_blank(),
                   axis.ticks.x.bottom  =  ggplot2::element_blank() ,
                   axis.line.x.bottom =  ggplot2::element_blank(),
                   axis.text.x.top  =  ggplot2::element_blank(),
                   axis.text.x.bottom =  ggplot2::element_blank(),
                   axis.line =  ggplot2::element_line(),
                   axis.title.y =  ggplot2::element_text(size = Y_Axis_Title_Size, vjust = 2),
                   axis.text.y =  ggplot2::element_text(size = Y_Axis_Text_Size, colour = "black"),
                   axis.ticks.length.x.top  =  ggplot2::unit(0.8,"cm"),
                   legend.position = "none",
                   plot.margin =  ggplot2::margin(t = 20,  # Top margin
                                        r = 20,  # Right margin
                                        b = 0,  # Bottom margin
                                        l = 35)) #Left margin


  #Labelling



  if(Bottom_Label_Index == TRUE)
  {
    print("Adding Labels")
    p1 <- p1 +  ggtext::geom_richtext( ggplot2::aes(label = Lab), hjust = 1, vjust = 1, label.size = 0,
                                      label.padding = grid::unit(rep(6, 4), "pt"), nudge_y = 0, fill = NA,
                                      nudge_x = 0, color = Bottom_Label_Colour, label.color = NA, size = Label_Size, angle = Label_Angle) +  ggplot2::theme(axis.title.x =  ggplot2::element_blank())

  }


  if(Bottom_Colour_Index == TRUE)
  {
    print("Highlighting Index SNPs")
    p1 <- p1 +   ggplot2::geom_point( ggplot2::aes(size=COLOUR),
                            pch=21, fill=NA,  colour=Bottom_Colour_Of_Index, stroke=1)

   p1 <- p1 + ggplot2::scale_size_identity()
  }

  print("Filling Background Noise")


  p1 <- p1 +   ggplot2::geom_point( ggplot2::aes(size=COLOUR2),
                                  fill=Bottom_Colour_One,  colour=Bottom_Colour_One)
  p1 <- p1 +   ggplot2::geom_point( ggplot2::aes(size=COLOUR3),
                                  fill=Bottom_Colour_Two,  colour=Bottom_Colour_Two)

  p1 <- p1 + ggplot2::scale_size_identity()

  p1 <- p1 +  ggplot2::coord_flip() +   ggplot2::coord_trans(y = "reverse",  clip = "off")



  print("Scaling Axes")

  df <- as.data.frame(p1$data)

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






  d1 <- p1 +  ggplot2::scale_x_continuous(breaks = middle_new_pos_values,
                                # Midpoints or specific points for chromosome labels
                                labels = factor(c(1:22, "X")),
                                position = "top",
                                expand =  ggplot2::expansion(mult = c(0.01, 0.01)),
                                sec.axis =  ggplot2::dup_axis(name = Bottom_Title))






  Bottom_Plot_Outcome <- d1

  print("Bottom Plot Made...")




  gA <-  ggplot2::ggplotGrob(Top_Plot_Outcome)
  gB <-  ggplot2::ggplotGrob(Bottom_Plot_Outcome)

  # Combine the two grobs using rbind
  combined_grob <- rbind(gA, gB, size = "first")





  print("Combined Plots")

  Overall_Name <- paste0(File_Name, ".", File_Type)

  print("Saving Miami Plot")

#programme in bitmap and other adjustments which are annoying in future.

  ggplot2::ggsave(Overall_Name, plot = combined_grob, width = Width,
         height = Height, units = "in", dpi = Quality)


  return(combined_grob) # has to be at end
}



