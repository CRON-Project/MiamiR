
#' Annotated, Customised Manhattan Plots
#'
#' @param Data These are the summary statistics (either file path(s), environmental object(s) or a variable(s) containing such contents) to be plotted; defaults to NULL
#' @param Chromosome_Column Manually specify chromosome column; defaults to NULL, leading to auto-detection
#' @param Position_Column Manually specify chromosome genomic position column; defaults to NULL, leading to auto-detection
#' @param SNP_ID_Column Manually specify SNP ID column; defaults to NULL, leading to auto-detection
#' @param PValue_Column Manually specify P Value column; defaults to NULL, leading to auto-detection
#' @param Reference_Allele_Column Manually specify reference allele column; defaults to NULL, leading to auto-detection
#' @param Effect_Allele_Column Manually specify effect allele column; defaults to NULL, leading to auto-detection
#' @param Title Manually specify title to be displayed; defaults to NULL, leading to auto-detection of dataframe object/file name
#' @param Y_Axis_Title Manually specify Y axis title; defaults to "-log₁₀P"
#' @param X_Axis_Title Manually specify X axis title; defaults to NULL, leading to assignment as Chromosome when called solely in Single_Plot() but Genomic Position (Mb) within Regional_Plot() by default.
#' @param Title_Size Size of title; defaults to 35
#' @param Y_Axis_Text_Size Size of Y axis text; defaults to 30
#' @param Chromosome_Label_Size Size of chromosome numbers (X axis text); defaults to 30
#' @param Y_Axis_Title_Size Size of Y axis title; defaults to 35
#' @param X_Axis_Title_Size Size of X axis title; defaults to 35
#' @param Chromosome_Colours Alternating list of two colours allocated to data points (odd/even chromosome numbers); defaults to c("blue", "turquoise") - specifying the same colour twice will lead to uniform colouring.
#' @param Sig_Line_Colour Colour of significance threshold; defaults to "red"
#' @param Sig_Line_Type Type of significance threshold line; defaults to "dashed"
#' @param Sig_Threshold Value of significance threshold; defaults to standard Bonferroni correction of 5e-8
#' @param Sig_Line_Width Width of significance threshold line; defaults to 0.5
#' @param Point_Size Size of data points; defaults to 2.5 but defaults to 10 when called within Regional_Plot()
#' @param Label_Index Annotate the index SNPs with ID provided; defaults to TRUE
#' @param Label_Size Size of index labels if Label_Index is TRUE; defaults to 5 but defaults to 9 when called within Regional_Plot()
#' @param Label_Angle Angle of index labels if Label_Index is TRUE; defaults to 60 (degrees)
#' @param Label_Colour Colour of index labels if Label_Index is TRUE; defaults to "black"
#' @param Colour_Index Highlight the index SNP with a circular ring around the original point; defaults to TRUE
#' @param Colour_Of_Index Colour of index SNP circular highlight if Colour_Index is TRUE; defaults to "darkred"
#' @param Chromosome_Labels List of chromosomes to show index SNP labels on, given that Label_Index is TRUE; defaults to c(1:22, "X")
#' @param Chromosome_Index  List of chromosomes to show index SNP circular highlight/colour on, given that Colour_Index is TRUE; defaults to c(1:22, "X")
#' @param Condense_Scale Log transform and condense scale above allotted break point to enhance visualisation of 'hits' with vastly different test statistics; defaults to TRUE
#' @param Break_Point Value of P where y axis begins to condense; defaults to 1e-10, guided by frequent max value observed in a small sample GWAS
#' @param Draft_Plot Specify whether a random subset of SNPs, the number of which is specified by Random_Selection, is plotted; defaults to FALSE
#' @param Random_Selection Manually select a random subset of SNPs to plot (more for internal debugging and testing), given that Draft_Plot is TRUE; defaults to 10,000
#' @param Chromosome_Label_Drops Manually specify which chromosome/X axis labels (as a list) should not be shown to prevent text crowding; defaults to NULL leading to c(21,22) being passed when X chromosome data is present and c(21) when only autosomes are present.
#' @param Interactive View the output as amenable to interactive plotly object, which can be used to inspect data, with key overlays (more for shiny App usage); defaults to FALSE
#' @param Index_Size Size of the outer ring diameter used to highlight the index SNPs, when Colour_Index is TRUE; defaults to NULL and if not manually passed is always scaled to 2X the passed or default Point_Size
#' @param Index_Thickness Thickness of the outer ring circumference used to highlight the index SNPs, when Colour_Index is TRUE; defaults to 1 but defaults to 4 when called within Regional_Plot()
#' @param Chromosome_Diamond List of chromosomes to show index SNP highlight/swapped diamond shape on, given that Diamond_Index is TRUE; defaults to c(1:22, "X")
#' @param Anchor_Label Which side of the point to anchor/position the index label to, when Label_Index is TRUE; defaults to "left"
#' @param Label_Height Spacing between the point and anchoring position of the index SNP label and the beginning of the label print area; defaults to 6 but defaults to 25 when called within Regional_Plot()
#' @param Diamond_Index Highlight the index SNP by swapping the standard point shape to a (different by default) coloured diamond; defaults to FALSE, unless called within Regional_Plot, then defaults to TRUE.
#' @param Colour_Of_Diamond Colour of diamond shape swapped for standard point shape ar index SNPs, when Diamond_Index is TRUE; defaults to "purple"
#' @param Diamond_Index_Size Size of the index SNP swapped diamond shape used to highlight the index SNPs, when Diamond_Index is TRUE; defaults to NULL and if not passed is scaled to 2X the according default or passed Point_Size
#' @param Lab Specify an additional info column with information to be labelled regarding Index SNPs, when Label_Index is TRUE; defaults to NULL leading to the values in the SNP ID column being used (Coords/RS codes) by default
#' @param Verbose Prevent display of progress bar as function is running and instead show key milestone outputs/messages (mainly for debugging purposes); defaults to FALSE.
#' @param Auto_Lab Allow Annotate_Data() to run on the data frame(s) provided so Lab can have RSIDs assigned, if only coordinates exist; defaults to FALSE
#' @param Genome_Build Reference genome to base required RSID annotations around; defaults to grch38
#'
#' @return Image of Manhattan Plot(s) is allocated to specified object and the resulting ggplot object can then be saved to an image
#' @export
#'
#' @examples  Manhattan_Plot <- Single_Plot(Data = Intelligence_Sum_Stats)

Single_Plot<- function(Data = NULL,
                       Random_Selection = 10000,
                       Draft_Plot = FALSE,
                       Title = NULL,
                       Title_Size = 35,
                       Y_Axis_Text_Size = 30,
                       X_Axis_Title_Size = 35,
                       Y_Axis_Title_Size = 35,
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
                       Index_Thickness = 1,
                       Chromosome_Labels = c(1:22, "X"),
                       Chromosome_Index = c(1:22, "X"),
                       Chromosome_Diamond = c(1:22, "X"),
                       Anchor_Label = "left",
                       Label_Index = TRUE,
                       Label_Size = 5, Label_Angle = 60,
                       Label_Colour = "black",
                       Label_Height = 6,
                       Colour_Of_Index = "darkred",
                       Colour_Index = TRUE,
                       Diamond_Index = FALSE,
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
                       Auto_Lab = FALSE,
                       Genome_Build = "grch38",
                       Verbose = TRUE)

{


  message("Arguments received, beginning Single_Plot() Script")

  if (is.null(Data)) {
    stop("Please provide data.", call. = FALSE)
  }
  is_called_by_regional_plot <- any(vapply(sys.calls(), function(x) {
    "Regional_Plot" %in% deparse(x[[1]])
  }, logical(1)))

  if (is.null(X_Axis_Title) &  is_called_by_regional_plot == FALSE) {
    X_Axis_Title <- "Chromosome"
  }

  is_called_by_miami_plot <- any(vapply(sys.calls(), function(x) {
    "Miami_Plot" %in% deparse(x[[1]])
  }, logical(1)))

  message(paste0("Called by Regional_Plot()?"), " ", is_called_by_regional_plot)
  message(paste0("Called by Miami_Plot()?"), " ", is_called_by_miami_plot)

  if(is_called_by_regional_plot)
  {
    message("Adjusting unpassed select defaults for Regional_Plot() call of Single_Plot()")
  }

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

  if(is.null(Index_Size))
  {
    message("Scaling Default Index Size")
    Index_Size <- Point_Size * 2
  }

  if(is.null(Diamond_Index))
  {
    message("Checking Diamond Index")
    Diamond_Index <- FALSE
  }

  if(Diamond_Index == TRUE & is.null(Diamond_Index_Size))
  {
    message("Scaling Default Diamond Index Size")
  }

  if(is.null(Diamond_Index_Size))
  {
    Diamond_Index_Size <- Point_Size * 2
  }

  if (is.character(Data) && length(Data) == 1)
  {
    file_path <- Data

    message("Dataset absent from environment")

    if (file.exists(file_path))
    {
     message("Reading data from file: ", file_path)

     message("Loading data using vroom...")

     Data <- vroom::vroom(file_path, show_col_types = FALSE)

     message("Finished reading")

    } else {

      stop("The provided character string does not point to an existing file: ", file_path, call. = FALSE)

         }
    }

  if(Draft_Plot == T)
  {
    message(paste0("Generating random subset of ", Random_Selection, " ", "points, specified as Draft_Plot activated"))
    Data <- dplyr::sample_n(Data, Random_Selection)
  }

  message("Deducing key column names")

  Chromosome_Column <- detect_chromosome_column(Data, Chromosome_Column)

  Data <- Data %>%
    dplyr::mutate(!!Chromosome_Column := ifelse(.data[[Chromosome_Column]] == "23", "X", .data[[Chromosome_Column]]))

  PValue_Column      <- detect_pvalue_column(Data, PValue_Column)
  Position_Column    <- detect_position_column(Data, Position_Column)
  SNP_ID_Column      <- detect_snp_column(Data, SNP_ID_Column)
  Ref_Allele_Column  <- detect_reference_allele_column(Data, Reference_Allele_Column)
  Alt_Allele_Column  <- detect_effect_allele_column(Data, Effect_Allele_Column)

  message("Assigning key columns to standardised nomenclature")

  Data$CHROM <- Data[[Chromosome_Column]]
  Data$GENPOS <- Data[[Position_Column]]
  Data$ID <- Data[[SNP_ID_Column]]
  Data$P <- Data[[PValue_Column]]
  Data$ALLELE0 <- Data[[Ref_Allele_Column]]
  Data$ALLELE1 <- Data[[Alt_Allele_Column]]


  if (is.null(Chromosome_Label_Drops)){

    message("Assigning auto Chromosome_Label_Drops")

  if (is.null(Chromosome_Label_Drops) && all(!(unique(Data$CHROM) %in% c("X", "x", 23, "chr23", "CHRX", "CHR23", "CHRx"))))
  {
    Chromosome_Label_Drops <- c(21)

  }else{

    Chromosome_Label_Drops <- c(21,22)

  }

  }

  max_logp <- max(-log10(Data[[PValue_Column]]), na.rm = TRUE)

  #fake anchoring around single point to prevent scale missing for single point plot/regional via inheritance

  if (nrow(Data) == 1) {

    message("Only one point provided - anchoring")

    # Original row
    original_row <- Data[1, ]

    pos_col <- Position_Column
    pval_col <- PValue_Column

    # Create a copy with position - 1
    fake_minus <- original_row
    fake_minus[[pos_col]] <- fake_minus[[pos_col]] - 1
    fake_minus[[pval_col]] <- 0.0000000000001
    fake_minus$ID <- "FAKEMINUS1"

    # Create a copy with position + 1
    fake_plus <- original_row
    fake_plus[[pos_col]] <- fake_plus[[pos_col]] + 1
    fake_plus[[pval_col]] <- 0.0000000000001
    fake_plus$ID <- "FAKEPLUS1"

    # Combine all
    Data <- dplyr::bind_rows(original_row, fake_minus, fake_plus)

  }


  if(!is.null(Lab))
  {
    message("No Lab/INFO column provided")
  }

  if(("Lab" %in% colnames(Data))) # if annotate used before then skip
  {
    message("It looks like you've applied Annotate_Data() or provided specific lead IDs/extra INFO")
  }

  if(Auto_Lab == TRUE) # if annotate used before then skip
  {
    message("Clearing existing Lab values")

    Data$Lab <- NA

    Data <- do.call(.Annotate_Data_original, list(Data = Data, Genome_Build = Genome_Build ))

    message("Running Auto_Lab from Annotate_Data()")

  }


  if(!("Lab" %in% colnames(Data))) # if annotate used before then skip
  {

  message("Discerning Lead SNP Per Chromosome to annotate")

  Data <- Data %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(min_P_GENPOS = GENPOS[which.min(P)]) %>% #P is formed in function
    dplyr::ungroup()

  # Label the minimum P value for each CHROM

  suppressMessages(suppressWarnings({

  Data <- Data %>%
    dplyr::mutate(Lab = ifelse(GENPOS == min_P_GENPOS & P < 5e-8, ID,
                        ifelse(abs(GENPOS - min_P_GENPOS) > 2000000000 & P == min(P[abs(GENPOS - min_P_GENPOS) > 2000000000]),
                               ID, "")))

  }))

  }


  if(Diamond_Index == TRUE)
  {

    message("Assigning Diamond Indices")



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


  #Backup df to draw from later
  Data2 <- Data

  TOPS <- dplyr::filter(Data, top == TRUE)

  }

  message("Setting Index Labels")

  Data$Lab[Data$Lab == ""] <- NA

  message("Setting Index Circle Highlights")

  suppressMessages(suppressWarnings({

  Data$COLOUR[is.na(Data$Lab)] <- NA
  Data$COLOUR[!is.na(Data$Lab)] <- 4

  }))

  message("Filtering Index Labels")

  Data$Lab[!(Data$CHROM %in% Chromosome_Labels)] <- NA

  message("Filtering Index Circle Highlights")

  Data$COLOUR[!(Data$CHROM %in% Chromosome_Index)] <- NA

  message("Formatting title spacing beneath")

  Title <- paste0(Title, "\n \n ") # plotly only interprets for y if actual content on line in app]

  message("Formatting X Axis Title spacing above")

  X_Axis_Title <- paste0("\n", X_Axis_Title) # needs to go above it to push down opp to above.

  if(Interactive == TRUE)

  {

    message("Setting Interactive hover INFO")

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

message("Plotting Top Framework")

if(length(unique(Data$CHROM)) == 1)

{

  message("Adapting for Regional_Plot()")

  suppressMessages(suppressWarnings({
  chroms <- unique(Data$CHROM)

  a <- ggmanh::manhattan_plot(x = Data, preserve.position = T, plot.title = Title,
                      chr.colname = Chromosome_Column, pos.colname = Position_Column, label.colname = NULL,
                      pval.colname = PValue_Column, annotateTop = FALSE, chr.order = c(chroms), chr.col = Chromosome_Colours,
                      chrlabs = c(1:22, "X"), rescale = Condense_Scale, signif = Break_Point, rescale.ratio.threshold = 0.0,
                      signif.rel.pos = 0.8, signif.col = c("transparent"),  point.size = Point_Size, x.label = "", y.label = "")

  # Define alternating chromosome colours:
  unique_chr <- sort(unique(a$data$CHROM))
  chr_colour_map <- setNames(
    rep(Chromosome_Colours, length.out = length(unique_chr)),
    unique_chr
  )


  # Assign alternating colours:
  a$data$colour_group <- chr_colour_map[as.character(a$data$CHROM)]

  }))
  if(Diamond_Index == TRUE)
  {

  message("Swapping lead SNPs for diamond index point shape")

    suppressMessages(suppressWarnings({

  a$data$colour_group[a$data$top == TRUE] <- "transparent"

    }))

  if(is_called_by_miami_plot){

    message("Shielding buffer range")

  }

    suppressMessages(suppressWarnings({

  a$data$colour_group[a$data$FAKE == 1] <- "transparent"

    }))

  }

  message("Verifying colour scale")

  # Update mapping:
  suppressMessages(suppressWarnings({

  a$mapping$colour <- rlang::quo(colour_group)

  # Auto-generate scale based on actual values present:
  colour_vals <- unique(a$data$colour_group)
  colour_map <- setNames(colour_vals, colour_vals)

  }))

  suppressMessages({

  a <- a + ggplot2::scale_colour_manual(values = colour_map)

})


  if(Diamond_Index == TRUE)
  {

  message("Creating backup reference plot for diamond index")

  suppressMessages(suppressWarnings({

  a2 <- ggmanh::manhattan_plot(x = Data2, preserve.position = T, plot.title = Title,
                              chr.colname = Chromosome_Column, pos.colname = Position_Column, label.colname = NULL,
                              pval.colname = PValue_Column, annotateTop = FALSE,
                              #      chr.order = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"),
                              chr.order = c(chroms),
                              chr.col = Chromosome_Colours, chrlabs = c(1:22, "X"), rescale = Condense_Scale,
                              signif = Break_Point, rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
                              signif.col = c("transparent"),  point.size = Point_Size, x.label = "", y.label = "")

  }))

  }

}else{

  message("Using standard for Single_Plot()")

  suppressMessages(suppressWarnings({

  a <- ggmanh::manhattan_plot(x = Data, preserve.position = T, plot.title = Title,
                              chr.colname = Chromosome_Column, pos.colname = Position_Column, label.colname = NULL,
                              pval.colname = PValue_Column, annotateTop = FALSE,
                              chr.order = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"),
                              chr.col = Chromosome_Colours, chrlabs = c(1:22, "X"), rescale = Condense_Scale,
                              signif = Break_Point, rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
                              signif.col = c("transparent"),  point.size = Point_Size, x.label = "", y.label = "")

  # Define alternating chromosome colours:
  unique_chr <- sort(unique(a$data$CHROM))
  chr_colour_map <- setNames(
    rep(Chromosome_Colours, length.out = length(unique_chr)),
    unique_chr
  )

  }))

  message("Verifying colour scale")

  # Assign alternating colours:
  a$data$colour_group <- chr_colour_map[as.character(a$data$CHROM)]

  message("Swapping lead SNPs for diamond index point shape")

  suppressMessages(suppressWarnings({

  # Override for top==TRUE:
  a$data$colour_group[a$data$top == TRUE & a$data$CHROM %in% Chromosome_Diamond] <- "transparent"

  }))

  if(is_called_by_miami_plot){

    message("Shielding buffer range")

  }

  suppressMessages(suppressWarnings({

  a$data$colour_group[a$data$FAKE == 1 & a$data$P == 1] <- "transparent"

  }))

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



  message("Creating backup reference plot for diamond index")

    suppressMessages(suppressWarnings({

  a2 <- ggmanh::manhattan_plot(x = Data2, preserve.position = T, plot.title = Title,
                              chr.colname = Chromosome_Column, pos.colname = Position_Column, label.colname = NULL,
                              pval.colname = PValue_Column, annotateTop = FALSE,
                              chr.order = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"),
                              #  chr.order = c(chroms),
                              chr.col = Chromosome_Colours, chrlabs = c(1:22, "X"), rescale = Condense_Scale,
                              signif = Break_Point, rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
                              signif.col = c("transparent"),  point.size = Point_Size, x.label = "", y.label = "")

    }))

  }

}

  message("Adding Significance Line")

  suppressMessages(suppressWarnings({

  b <- a + ggplot2::geom_hline(yintercept = -(log10(Sig_Threshold)), linetype = Sig_Line_Type,
                      color = Sig_Line_Colour, size = Sig_Line_Width)

  }))


  if(Diamond_Index == TRUE)
  {

  suppressMessages(suppressWarnings({

  b2 <- a2 + ggplot2::geom_hline(yintercept = -(log10(Sig_Threshold)), linetype = Sig_Line_Type,
                               color = Sig_Line_Colour, size = Sig_Line_Width)

  }))

  }

  message("Formatting Plot")

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

  message("Formatting axes and labels")

  p <- p +  ggplot2::xlab(X_Axis_Title) + ggplot2::theme(axis.title.x =  ggplot2::element_text(vjust=-5, size = X_Axis_Title_Size, hjust = 0.5))
  p <- p +  ggplot2::ylab(Y_Axis_Title) + ggplot2::theme(axis.title.y =  ggplot2::element_text(vjust=7, size = Y_Axis_Title_Size))


  if(Diamond_Index == TRUE)
  {

  p2 <- p2 +  ggplot2::xlab(X_Axis_Title) + ggplot2::theme(axis.title.x =  ggplot2::element_text(vjust=-5, size = X_Axis_Title_Size, hjust = 0.5))
  p2 <- p2 +  ggplot2::ylab(Y_Axis_Title) + ggplot2::theme(axis.title.y =  ggplot2::element_text(vjust=5, size = Y_Axis_Title_Size))


  }

  if(Interactive == TRUE)
  {

  message("Adjusting axes and labels for interactive mode")

  dot_padding <- "\n \n "  # You can add more lines if needed
  Y_Axis_Title <- paste0(Y_Axis_Title, dot_padding)
  X_Axis_Title <- paste0(dot_padding, X_Axis_Title)

  p <- p +  labs(y = Y_Axis_Title, x = X_Axis_Title)

  if(Diamond_Index == TRUE)
  {

  p2 <- p2 +  labs(y = Y_Axis_Title, x = X_Axis_Title)

  }
  }


  if(Interactive == TRUE)

  {

  message("Etching information for interactive mode")

  p <- p + ggplot2::aes(text = Hover_Info)

  if(Diamond_Index == TRUE)
  {

  p2 <- p2 + ggplot2::aes(text = Hover_Info)

  }

  }

  if(Label_Index == TRUE)

  {
   message("Anchoring Index Label text")
  }

  message("Anchoring ")

  if(Anchor_Label == "left")
  {
    HJUST = 0
    VJUST = 0
  }
  if (Anchor_Label == "left mirror") {
    HJUST <- 1   # mirror horizontally
    VJUST <- 1   # keep vertical alignment unchanged, or adjust if needed
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

    suppressMessages(suppressWarnings({

    p <- p +  ggtext::geom_richtext( ggplot2::aes(label = Lab), hjust = HJUST, vjust = VJUST, label.size = 0,
                                    label.padding = grid::unit(rep(Label_Height, 4), "pt"), nudge_y = 0, fill = NA,
                                    nudge_x = 0, color = Label_Colour, label.color = NA, size = Label_Size, angle = Label_Angle)

    }))

    if(Diamond_Index == TRUE)
    {

    message("Adding Diamonds")

    suppressMessages(suppressWarnings({

    p2 <- p2 +  ggtext::geom_richtext( ggplot2::aes(label = Lab), hjust = HJUST, vjust = VJUST, label.size = 0,
                                     label.padding = grid::unit(rep(Label_Height, 4), "pt"), nudge_y = 0, fill = NA,
                                     nudge_x = 0, color = Label_Colour, label.color = NA, size = Label_Size, angle = Label_Angle)

    }))

    }

  }


  suppressMessages(suppressWarnings({


  if(Colour_Index == TRUE & Interactive == FALSE & Diamond_Index == FALSE)
  {

    message("Colouring Index circle")

    df <- as.data.frame(p$data)

    #special_points <- df[df$COLOUR == 4, ]
    special_points <- df %>% dplyr::filter(COLOUR == 4)

    if(is_called_by_regional_plot == TRUE)

    {

    p <- p + ggrastr::geom_point_rast(
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

    message("Colouring Diamond Index SNPs")

    df <- as.data.frame(p2$data)

    suppressMessages(suppressWarnings({

    special_points <- df[df$top == TRUE, ]

    }))

    if (all(c("P", "FAKE") %in% colnames(df))) {

    chroms_all_fake <- df %>%
       dplyr::group_by(CHROM) %>%
       dplyr::filter(all(P == 1 & FAKE == 1)) %>%
       dplyr::pull(CHROM) %>%
       unique()

     # Filter out such chromosomes from special_points
     special_points <- special_points[!(special_points$CHROM %in% chroms_all_fake), ]

}

    message("Filtering Diamond Index SNPs")

    special_points <- special_points[special_points$CHROM %in% Chromosome_Diamond, ]


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

  message("Resolving Axes Area")

  p <- p +  ggplot2::coord_cartesian(clip = "off")


  }))

  message("Scaling Axes")

  df <- as.data.frame(p$data)

  middle_new_pos_values <- numeric()

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

  # Define maximum full label set
  all_chromosomes <- as.character(c(1:22, "X"))

  # Filter out dropped and missing chromosomes
  chromosome_labels <- setdiff(all_chromosomes, as.character(Chromosome_Label_Drops))
  chromosome_labels <- intersect(chromosome_labels, chroms_in_data)

  # Ensure label order matches middle_new_pos_values
  final_labels <- ifelse(all_chromosomes %in% chromosome_labels, all_chromosomes, "")

  suppressMessages(suppressWarnings({

  d <- p +  ggplot2::scale_x_continuous(breaks = middle_new_pos_values,
                                        labels = final_labels,
                                        position = "bottom",
                                       expand =  ggplot2::expansion(mult = c(0.01, 0.01)) )
  }))


  Plot_Outcome <- d

  if (is_called_by_regional_plot) {

  idx <- which(sapply(Plot_Outcome$scales$scales, function(s) "y" %in% s$aesthetics))
  Plot_Outcome$scales$scales[[idx]]$expand <- ggplot2::expansion(mult = c(0.05, 0))

  }

  message("Plot Complete")

  return(invisible(Plot_Outcome)) # invisible prevents slow print auto at end.

}

.Single_Plot_original <- Single_Plot

Single_Plot <- function(..., session = NULL) {
  call_expr <- match.call()
  fn_formals <- formals(.Single_Plot_original)
  valid_args <- names(fn_formals)
  dots <- list(...)
  clean_args <- dots[names(dots) %in% valid_args]

  raw_data_expr <- call_expr[["Data"]]

  if (!is.null(raw_data_expr) && is.call(raw_data_expr) && identical(raw_data_expr[[1]], as.name("c"))) {


    arg_exprs <- as.list(raw_data_expr[-1])
    arg_vals <- lapply(arg_exprs, function(e) eval(e, parent.frame()))

    if (all(vapply(arg_vals, is.data.frame, logical(1)))) {

      names(arg_vals) <- vapply(arg_exprs, function(e) deparse(e), character(1))
      clean_args$Data <- arg_vals

    } else if (all(vapply(arg_vals, is.character, logical(1)))) {

      if (all(file.exists(unlist(arg_vals)))) {

        data_list <- lapply(unlist(arg_vals), function(f) vroom::vroom(f, show_col_types = FALSE))
        names(data_list) <- basename(tools::file_path_sans_ext(unlist(arg_vals)))
        clean_args$Data <- data_list
      } else {

        stop("One or more file paths in c(...) do not exist.", call. = FALSE)

      }
    } else {

      stop("Mixed or invalid types in c(...) for Data.", call. = FALSE)

    }
  }

  Data <- clean_args$Data
  is_df <- function(x) is.data.frame(x)

  if (is.character(Data) && length(Data) > 1 && all(file.exists(Data))) {

    Data_list <- lapply(Data, vroom::vroom, show_col_types = FALSE)
    names(Data_list) <- basename(tools::file_path_sans_ext(Data))
    clean_args$Data <- Data_list
    Data <- Data_list

  }

  if (is.list(Data) && all(vapply(Data, is_df, logical(1)))) {

    n_data <- length(Data)
    data_names <- names(Data)

    plots <- lapply(seq_along(Data), function(i) {

      df <- Data[[i]]

      name <- if (!is.null(data_names) && nzchar(data_names[[i]])) data_names[[i]] else paste0("Data_", i)

      message(sprintf("Processing dataset: %s", name))

      args_i <- clean_args
      args_i$Data <- df

      user_supplied_title <- "Title" %in% names(dots) && !is.null(dots$Title) && nzchar(dots$Title)

      if (!user_supplied_title) {
        args_i$Title <- name
      } else if (length(dots$Title) == n_data) {
        args_i$Title <- dots$Title[[i]]
      } else {
        args_i$Title <- dots$Title
      }

      for (arg in valid_args) {

        if (arg %in% c("Data", "Title")) next
        supplied_val <- clean_args[[arg]]
        default_val <- fn_formals[[arg]]

        if (!is.null(supplied_val)) {
          if (length(supplied_val) == n_data) {
            args_i[[arg]] <- supplied_val[[i]]
          } else {
            args_i[[arg]] <- supplied_val
          }
        } else {
          val <- tryCatch(
            eval(default_val, envir = environment(.Single_Plot_original)),
            error = function(e) NULL
          )
          args_i[[arg]] <- val
        }

      }

      if (isTRUE(clean_args$Verbose)) {

        do.call(.Single_Plot_original, args_i)

      } else {

        run_with_counter(.Single_Plot_original, args = args_i, session = session)

      }
    })

    names(plots) <- data_names
    return(invisible(plots))
  }

  if (is.character(Data) && length(Data) == 1 && file.exists(Data)) {

    name <- basename(tools::file_path_sans_ext(Data))
    message(sprintf("Processing file: %s", name))
    Data_df <- vroom::vroom(Data, show_col_types = FALSE)
    clean_args$Data <- Data_df

    if (is.null(clean_args$Title) || !nzchar(clean_args$Title)) {

      clean_args$Title <- name

    }

  } else if (is.data.frame(Data)) {

    expr <- raw_data_expr
    name <- if (!is.null(expr)) deparse(expr) else "Data_1"
    message(sprintf("Processing dataset: %s", name))

    if (is.null(clean_args$Title) || !nzchar(clean_args$Title)) {

      clean_args$Title <- name

    }

  }

  for (arg in valid_args) {

    if (arg == "Data") next

    if (is.null(clean_args[[arg]])) {

      default_val <- fn_formals[[arg]]
      val <- tryCatch(
        eval(default_val, envir = environment(.Single_Plot_original)),
        error = function(e) NULL
      )
      clean_args[[arg]] <- val

    }
  }

  if (isTRUE(clean_args$Verbose)) {

    do.call(.Single_Plot_original, clean_args)
  } else {

    run_with_counter(.Single_Plot_original, args = clean_args, session = session)
  }
}

