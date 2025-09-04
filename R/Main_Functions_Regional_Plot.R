
#' Annotated, Customised Regional Plots
#'
#' @param Data These are the summary statistics (either file path(s), environmental object(s) or a variable(s) containing such contents) to be plotted; defaults to NULL
#' @param Chromosome_Column Manually specify chromosome column; defaults to NULL, leading to auto-detection
#' @param Position_Column Manually specify chromosome genomic position column; defaults to NULL, leading to auto-detection
#' @param SNP_ID_Column Manually specify SNP ID column; defaults to NULL, leading to auto-detection
#' @param PValue_Column Manually specify P Value column; defaults to NULL, leading to auto-detection
#' @param Chromosomes Manually specify which chromosomes from the summary statistics to run the function on; defaults to NULL, leading to all being included/absent filtering
#' @param MIN_GENE_GAP Minimum genomic distance between adjacent annotation anchoring/positioning for Gene_Tracks body names, when Gene_Tracks = TRUE; defaults to NULL and if not passed is calculated from 0.3 X the Region_Window size specified.
#' @param Separation_Distance The minimum genomic distance between significant regions of interest/lead SNPs to be defined as distinct peaks, defined by P_threshold value; defaults to 1e6
#' @param Region_Window Distance to be included up and downstream of the lead SNP in a region of interest; defaults to 5e5
#' @param Gene_Label_Size Size of the text allocated to the bodies displayed in the gene panel, when Gene_Tracks is TRUE; defaults to 10
#' @param Gene_Structure_Size_Exons Size of the exon object allocated to the bodies displayed in the gene panel, when Gene_Tracks is TRUE; defaults to 12
#' @param Gene_Structure_Size_Introns Size of the intron object allocated to the bodies displayed in the gene panel, when Gene_Tracks is TRUE; defaults to 2
#' @param Gene_Structure_Colour_Exons Colour of the exon object allocated to the bodies displayed in the gene panel, when Gene_Tracks is TRUE; defaults to black
#' @param Gene_Structure_Colour_Introns Colour of the intron object allocated to the bodies displayed in the gene panel, when Gene_Tracks is TRUE; defaults to black
#' @param Gene_Panel_Background_Colour Colour of the background panel surrounding the bodies displayed in the gene panel, when Gene_Tracks is TRUE; defaults to azure
#' @param Sense_Arrow Show direction of transcription arrow in the gene panel, when Gene_Tracks is TRUE; defaults to TRUE
#' @param Sense_Arrow_Colour Colour of direction of transcription arrow in the gene panel, when Gene_Tracks is TRUE; defaults to black
#' @param Sense_Arrow_Head_Size Size of direction of transcription arrow head/tip in the gene panel, when Gene_Tracks is TRUE; defaults to 10
#' @param Sense_Arrow_Body_Thickness Width of direction of transcription arrow body in the gene panel, when Gene_Tracks is TRUE; defaults to NULL and when not passed is calculated based on Gene_Structure_Size_Introns X 0.6
#' @param Sense_Arrow_Body_Length Length of direction of transcription arrow body in the gene panel, when Gene_Tracks is TRUE; defaults to 0.02
#' @param Sense_Arrow_Gene_Gap Minimum genomic distance (spacing) between end of Gene_Tracks bodies and beginning of direction of transcription arrow rendering, when Gene_Tracks = TRUE; defaults to 0.005 (factor of view window, proxy for Region_Window)
#' @param Recombination_Line Display recombination rate via line plot overlay across region shown; defaults to TRUE
#' @param Recombination_Line_Colour Colour of recombination rate line across region shown when Recombination Line is TRUE; defaults to darkred
#' @param Recombination_Axis_Title Y Axis Title of recombination rate scale when Recombination_Line is TRUE; defaults to Recombination Rate (cM/Mb)
#' @param Recombination_Axis_Text_Size Y Axis text size of recombination rate scale when Recombination_Line is TRUE; defaults to NULL and if not passed is based on inherited Single_Plot() Y_Axis_Text_Size and its default value
#' @param Recombination_Axis_Title_Size Y Axis title size of recombination rate scale when Recombination_Line is TRUE; defaults to NULL and if not passed is based on inherited Single_Plot() Y_Axis_Title_Size and its default value
#' @param Recombination_Line_Thickness Width of recombination rate line across region shown when Recombination_Line is TRUE; defaults to 1
#' @param Recombination_Line_Type Line style of recombination rate line across region shown when Recombination_Line is TRUE; defaults to dotted
#' @param P_threshold Maximum value of P (minimum LOG10P) allowed for a SNP to be classified as initiating detection as the lead SNP of a peak in the function; defaults to 5e-8
#' @param Point_Colour Colour of data points shown; defaults to darkblue
#' @param LD_Legend_Title_Size Size of title in Auto_LD Legend when both LD_Legend_On and Auto_LD are TRUE; defaults to 30
#' @param LD_Legend_Text_Size Size of text in Auto_LD Legend key when both LD_Legend_On and Auto_LD are TRUE; defaults to 20
#' @param LD_Legend_Location Relative location of LD Legend in the top panel when both LD_Legend_On and Auto_LD are TRUE; defaults to "Top Right"
#' @param LD_Legend_On Show Auto_LD Legend, when Auto_LD is TRUE; defaults to TRUE
#' @param Gene_Legend_Location Relative location of Gene_Tracks Legend in the bottom panel when both Gene_Legend_On and Gene_Tracks are TRUE; defaults to "Right Outside"
#' @param Gene_Legend_On Show Gene_Legend, when Gene_Tracks is TRUE; defaults to TRUE
#' @param Gene_Legend_Title_Size Size of title in the Gene_Tracks legend when both Gene_Legend_On and Gene_Tracks are TRUE; defaults to 30
#' @param Gene_Legend_Text_Size Size of text in the Gene_Tracks Legend key when both Gene_Legend_On and Gene_Tracks are TRUE; defaults to 20
#' @param Gene_Legend_Title Gene_Tracks Legend heading when both Gene_Legend_On and Gene_Tracks are TRUE; defaults to Gene Biotype
#' @param LD_Legend_Title Auto_LD Legend heading when both LD_Legend_On and Auto_LD are TRUE; defaults to R² LD
#' @param Genome_Build Reference genome to base required plot panels and annotations around (e.g. Auto_LD and Gene_Tracks functionality); defaults to grch38
#' @param Population Reference population for Auto_LD query; defaults to "EUR"
#' @param Gene_Biotype Filter bodies shown when Gene_Tracks is TRUE by biotype, e.g. coding; defaults to NULL allowing for all available biotypes to be shown.
#' @param Gene_Biotype_Colour Colour of bodies shown when Gene_Tracks is TRUE, either overall or stratified by biotype, e.g. c("blue", "red); defaults to NULL allowing for default colouring panel to be allotted
#' @param Gene_Name Filter bodies shown when Gene_Tracks is TRUE by body name, e.g. "DMD"; defaults to NULL allowing for all available bodies to be shown.
#' @param Auto_LD_Region_Size Region to query centered around the ascertained lead SNP in the peak region, when Auto_LD is TRUE; defaults to 1e6
#' @param Auto_LD Query for LD via the included LDlinkR::LDproxy functionality, based around Genome_Build & Population arguments, centered around the ascertained lead SNP in the region shown; defaults to FALSE
#' @param Reference_Allele_Column Manually specify reference allele column; defaults to NULL, leading to auto-detection
#' @param Effect_Allele_Column Manually specify effect allele column; defaults to NULL, leading to auto-detection
#' @param Verbose Prevent display of progress bar as function is running and instead show key milestone outputs/messages (mainly for debugging purposes); defaults to FALSE.
#' @param Interactive View the output as amenable to interactive plotly object, which can be used to inspect data, with key overlays (more for shiny App usage); defaults to FALSE
#' @param Gene_Tracks Display gene panel below plot with specified bodies within the genomic region shown; defaults to TRUE
#' @param .dots Shadow argument to detect inherited arguments passed and modified from the Single_Plot() function; defaults to list()
#'
#' @returns Image of Regional Plot(s) is allocated to specified object and the resulting ggplot object can then be saved to an image
#' @export
#'
#' @examples  Regional_Plot <- Regional_Plot(Data = Household_Income_Sum_Stats, Chromosomes = 1)

Regional_Plot <- function(Data = NULL,
                          Chromosome_Column = NULL,
                          Position_Column = NULL, SNP_ID_Column = NULL,
                          PValue_Column = NULL,
                          Chromosomes =  NULL,
                          MIN_GENE_GAP = NULL,
                          Separation_Distance = 1e6,
                          Region_Window = 5e5,
                          Gene_Label_Size = 10,
                          Gene_Structure_Size_Exons = 12,
                          Gene_Structure_Size_Introns = 2,
                          Gene_Structure_Colour_Exons = "black",
                          Gene_Structure_Colour_Introns = "black",
                          Gene_Panel_Background_Colour = "azure",
                          Sense_Arrow = TRUE,
                          Sense_Arrow_Colour = "black",
                          Sense_Arrow_Head_Size = 10,
                          Sense_Arrow_Body_Thickness = NULL,
                          Sense_Arrow_Body_Length = 0.02,
                          Sense_Arrow_Gene_Gap = 0.005,
                          Recombination_Line = TRUE,
                          Recombination_Line_Colour = "darkred",
                          Recombination_Axis_Title = "Recombination Rate (cM/Mb)",
                          Recombination_Axis_Text_Size = NULL,
                          Recombination_Axis_Title_Size = NULL,
                          Recombination_Line_Thickness = 1,
                          Recombination_Line_Type = "dotted",
                          P_threshold = 5e-8,
                          Point_Colour = "darkblue",
                          LD_Legend_Title_Size = 30,
                          LD_Legend_Text_Size = 20,
                          LD_Legend_Location = "Top Right",
                          LD_Legend_On = TRUE,
                          Gene_Legend_Location = "Right Outside",
                          Gene_Legend_On = TRUE,
                          Gene_Legend_Title_Size = 30,
                          Gene_Legend_Text_Size = 20,
                          Gene_Legend_Title = "Gene Biotype",
                          LD_Legend_Title = "R² LD",
                          Genome_Build = "grch38",
                          Population = "EUR",
                          Gene_Biotype = NULL,
                          Gene_Biotype_Colour = NULL,
                          Gene_Name = NULL,
                          Auto_LD_Region_Size = 1e6,
                          Auto_LD = FALSE,
                          Reference_Allele_Column = NULL,  Effect_Allele_Column = NULL,
                          Verbose = FALSE,
                          Interactive = FALSE, # full T/F to TRUE/FALSE critical
                          Gene_Tracks = TRUE,  .dots = list()) {


  if(Gene_Tracks == TRUE)
  {
   message("Assigning later required recommended MIN_GENE_GAP")
  }

  default_gap_fraction <- 0.3

  MIN_GENE_GAP <- if (!is.null(MIN_GENE_GAP)) {
    MIN_GENE_GAP
  } else {
    default_gap_fraction * Region_Window
  }

  message("Determining Plot Titles Fron Data Types")

  Title <- if (!is.null(attr(Data, "RegionalPlot_Title"))) {
    attr(Data, "RegionalPlot_Title")
  } else if (is.character(Data) && length(Data) == 1) {
    basename(Data)
  } else if (!is.null(attr(Data, "name"))) {
    attr(Data, "name")
  } else {
    "Untitled Regional Plot"
  }

  if (is.character(Data) && length(Data) == 1) {
    file_path <- Data

    message("Dataset absent from environment")

    if (file.exists(file_path)) {

    message("Reading data from file: ", file_path)

    message("Loading data using vroom...")

    Data <- vroom::vroom(file_path, show_col_types = FALSE)

    message("Finished reading")

    } else {
      stop("The provided character string does not point to an existing file: ", file_path, call. = FALSE)
    }
  }

  message("Deducing key column names")

  Chromosome_Column <- detect_chromosome_column(Data)

  Data[[Chromosome_Column]] <- as.character(Data[[Chromosome_Column]])

  PValue_Column     <- detect_pvalue_column(Data, PValue_Column)
  Position_Column   <- detect_position_column(Data, Position_Column)
  SNP_ID_Column     <- detect_snp_column(Data, SNP_ID_Column)
  Ref_Allele_Column     <- detect_reference_allele_column(Data, Reference_Allele_Column)
  Alt_Allele_Column     <- detect_effect_allele_column(Data, Effect_Allele_Column)

  message("Assigning key columns to standardised nomenclature")

  Data$CHROM <- Data[[Chromosome_Column]]
  Data$GENPOS <- Data[[Position_Column]]
  Data$ID <- Data[[SNP_ID_Column]]
  Data$P <- Data[[PValue_Column]]
  Data$ALLELE0 <- Data[[Ref_Allele_Column]]
  Data$ALLELE1 <- Data[[Alt_Allele_Column]]
  Data$ALLELE0 <- toupper(Data$ALLELE0 )
  Data$ALLELE1 <- toupper(Data$ALLELE1 )

  filtered_data <- Data  #Create Backup

  filtered_data$CHROM <- as.numeric(filtered_data$CHROM)

  if(!is.null(Chromosomes))

  {

  message(paste0("Filtering for key chromosomes: ", Chromosomes, "\n"))

  filtered_data <- filtered_data[filtered_data$CHROM %in% Chromosomes,]

  }

  message("Defining function to find peaks based on provided perameters")

  find_peaks <- function(data, Separation_Distance, P_threshold) {

    # Group by chromosome (or single group if only one CHROM)
    chrom_list <- if (length(unique(data$CHROM)) == 1) list(data) else base::split(data, data$CHROM)

    results <- base::lapply(chrom_list, function(chrom_data) {

      # Filter to only SNPs below genome-wide threshold - MAXIMUM P Value for significant hit/minimum LOG10P
      chrom_data <- chrom_data[chrom_data$P < P_threshold, , drop = FALSE]

      # If none pass threshold, return nothing
      if (nrow(chrom_data) == 0) return(chrom_data)

      # Sort by smallest P first
      chrom_data <- chrom_data[base::order(chrom_data$P), ]

      selected <- chrom_data[0, , drop = FALSE]

      for (i in base::seq_len(nrow(chrom_data))) {
        snp <- chrom_data[i, , drop = FALSE]

        if (nrow(selected) == 0) {
          selected <- snp
        } else {
          # Only keep if this SNP is at least Separation_Distance away from all previously selected - basically loop through peaks.
          too_close <- base::any(base::abs(snp$GENPOS - selected$GENPOS) < Separation_Distance)
          if (!too_close) {
            selected <- dplyr::bind_rows(selected, snp)
          }
        }
      }

      selected
    })

    dplyr::bind_rows(results)
  }

  message("Calling find_peaks algorithm based on user supplied significance and distinct peak distance threholds")

  peaks <- find_peaks(filtered_data, Separation_Distance = Separation_Distance, P_threshold = P_threshold)

  message("Lead SNPs of defined peaks ascertained")

  if(Auto_LD == TRUE)
  {

  message("Creating query coordinates for LD")

  peaks <- peaks %>%
    dplyr::mutate(
      LD_Query = paste0("chr", CHROM, ":", GENPOS)
    )

  ld_results_list <- list()

  }

  if(Auto_LD == TRUE)
  {

  for (i in seq_len(nrow(peaks))) {

    message("Calculating LD for lead SNP from peak ", i, " ...")

    query <- peaks$LD_Query[i]

    token <- "a3a5b2b4d4c5"

    result <- LDlinkR::LDproxy(
      snp = query,      # Query SNP
      pop = Population, # Population: EUR, AFR, EAS, AMR, SAS
      r2d = "r2",       # LD metric: "r2" or "d"
      token = token,    # Or use get_token() if you've set it
      win_size = Auto_LD_Region_Size,
      genome_build = Genome_Build
    )

    #INFO:
    #https://www.rdocumentation.org/packages/LDlinkR/versions/1.4.0/topics/LDproxy
    #Choose between one of the three options...`grch37` for genome build GRCh37 (hg19), `grch38` for GRCh38 (hg38),
    #or `grch38_high_coverage` for GRCh38 High Coverage (hg38) 1000 Genome Project data sets. Default is GRCh37 (hg19)

    # Store result in the list with identifier
    ld_results_list[[i]] <- result %>%
     dplyr::mutate(Peak_Index = i, Lead_SNP = peaks$SNP[i])

  }

  }

  message("Defining function to get region around now defined peak lead SNPs based on provided perameters")

  get_region_around_all_peaks <- function(filtered_data, peaks, Region_Window) {
    if (nrow(peaks) == 0) return(filtered_data[0, , drop = FALSE])  # return empty df if no peaks

    # Store list of regional data frames
    region_list <- lapply(seq_len(nrow(peaks)), function(i) {
      peak <- peaks[i, ]
      chrom <- peak$CHROM
      pos <- peak$GENPOS

      # Subset to same chromosome
      chr_data <- filtered_data[filtered_data$CHROM == chrom, , drop = FALSE]

      # Subset to window around peak (mirrored where available)
      region_data <- chr_data[chr_data$GENPOS >= (pos - Region_Window) & chr_data$GENPOS <= (pos + Region_Window), ]

      # Add a Peak_ID column
      region_data$Peak_ID <- paste0("peak", i)

      return(region_data)
    })

    # Combine all regions into one data frame
    combined_regions <- do.call(rbind, region_list)

    return(combined_regions)
  }

  message("Calling get_region_around_all_peaks algorithm based on lead SNPs uncovered and specified user Region_Window")

  region_df <- get_region_around_all_peaks(filtered_data, peaks, Region_Window = Region_Window)

  message("Adding peak IDs")
  region_df$pretty_peak_id_singleton <- gsub("peak(\\d+)", "Peak \\1", region_df$Peak_ID)
  region_df$Singleton_Name <- paste0("Chromosome ", region_df$CHROM, " ", region_df$pretty_peak_id_singleton, " Regional Plot")

  message("Finding Peak IDs and storing for Title matching later...")

  # Identify singleton Peak_IDs
  singleton_peaks <- region_df %>%
    dplyr::count(Peak_ID) %>%
    dplyr::filter(n == 1) %>%
    dplyr::pull(Peak_ID)

  # Get the original position column name
  pos_col <- Position_Column
  id_col <- SNP_ID_Column
  p_col <- PValue_Column

  # Map from Peak_ID → Position for singleton peaks
  singleton_pos_table <- region_df %>%
    dplyr::filter(Peak_ID %in% singleton_peaks) %>%
    dplyr::select(Singleton_Name, !!rlang::sym(pos_col)) %>%
    tibble::deframe() %>%
    as.list()

  # Build complete named list matching plot titles
  all_titles <- unique(region_df$Singleton_Name)

  singleton_positions <- setNames(
    lapply(all_titles, function(title) {
      if (title %in% names(singleton_pos_table)) singleton_pos_table[[title]] else NA
    }),
    all_titles
  )

  filtered_data <- region_df

  message("Dealing with inherited Single_Plot() defaults and passed arguments from within Regional_Plot")

  # Capture user arguments
  user_args <- .dots

  # Pull defaults from Single_Plot
  default_args <- formals(.Single_Plot_original) #changed to incorporate loading bar

  evaluated_defaults <- lapply(default_args, function(x) {
    if (missing(x)) return(NULL)  # just in case, though rarely hit
    if (is.symbol(x) || is.language(x)) {
      tryCatch(eval(x), error = function(e) NULL)
    } else {
      x
    }
  })

  # Merge user args over defaults (user overrides default)
  args <- modifyList(evaluated_defaults, user_args)

  message("Injecting and extracting required processed arguments and defaults for later")

  Colour_Of_Diamond <- args$Colour_Of_Diamond

  # Inject required args
  args$Data <- filtered_data

  if(Interactive == TRUE)
  {
    args$Interactive <- TRUE

  }else{

    args$Interactive <- FALSE
  }

  message("Defining function to plot defined peak lead SNPs and region windows in addition to any annotation, customisation and panels required")

  plot_all_peak_regions <- function(region_df, user_args) {
    if (!"Peak_ID" %in% names(region_df)) {
      stop("The dataframe must include a 'Peak_ID' column.")
    }

    region_df$pretty_peak_id <- gsub("peak(\\d+)", "Peak \\1",  region_df$Peak_ID)

    region_df$Title_Name <- paste0(Title, ": ", "Chromosome ", region_df$CHROM, " ", region_df$pretty_peak_id, " " ,"Regional Plot") # real title with spaces

    region_df$Save_Name <- paste0(Title, ": ",  "Chromosome ", region_df$CHROM, " ", region_df$pretty_peak_id, " " ,"Regional Plot") # real title with spaces

    Title_Names <- unique(region_df$Save_Name)

    message("Defining what peaks exist to loop apply logic to")

    peak_ids <- unique(region_df$Peak_ID)

    CHROM_Lab <- unique(region_df$CHROM)

    default_args <- formals(Single_Plot)

    # Remove missing arguments (i.e., no default)
    default_args <- formals(Single_Plot)
    default_args <- Filter(function(x) !identical(x, quote(expr = )), default_args)

    evaluated_defaults <- lapply(default_args, function(x) {
      if (is.symbol(x) || is.language(x)) {
        tryCatch(eval(x), error = function(e) NULL)
      } else {
        x
      }
    })

    args_base <- modifyList(evaluated_defaults, user_args)

    message("Within this apply function to list plots when multiple peaks... likely")

    # Create list of plots to execute logic on
    plots <- lapply(seq_along(peak_ids), function(i) {

      peak_id <- peak_ids[i]

      message("Segregating data to individual peak")

      filtered_data <- region_df[region_df$Peak_ID == peak_id, , drop = FALSE]

      message("Replicating LD query column creation for any LD information joining to regional plot")

      #need both due to data typing
      filtered_data <- filtered_data %>%
        dplyr::mutate(Coord = paste0("chr", CHROM, ":", as.integer(GENPOS))) %>%
        dplyr::mutate(
          LD_Query = paste0("chr", CHROM, ":", GENPOS)
        )

      if(Auto_LD == TRUE)

      {

      message("retrieving LD results")

      # Get corresponding LD result
      ld_df <- ld_results_list[[i]]

      }

      if(Auto_LD == TRUE)

      {
        if (!is.null(ld_df)) {

          message("Preparing LD results for join")

          ld_df_clean <- ld_df
          ld_df_clean$Proxy_SNP <- ld_df_clean$RS_Number
          ld_df_clean$R2_LD <- ld_df_clean$R2

          if (!"Coord" %in% colnames(ld_df_clean)) {

            message("Allotting NA as nothing found")

            ld_df_clean$Coord <- NA_character_
            ld_df_clean$R2_LD <- NA_character_ #need na values for later to grey
        }

        #sometimes lead snp missing from ld database and coord fails return in that case add NA coord
        # Join with filtered_data

          message("Preparing LD results for join")

           filtered_data <- filtered_data %>%
            dplyr::left_join(ld_df_clean, by = "Coord")

        } else {

          message("No LD results to join, perhaps lead SNP absent, or region too small for pairwise")

          filtered_data$R2_LD <- NA
          filtered_data$Proxy_SNP <- NA
        }

      }

      if (nrow(filtered_data) == 0) return(NULL)  # skip empty df

      #save real args for later injections
      args_real <- args

      #args_real used to control outputs

       if(is.null(args_real$X_Axis_Title))
      {

        args_real$X_Axis_Title <- "\nGenomic Position (Mb)"

      }else{

        args_real$X_Axis_Title <- paste0("\n", args_real$X_Axis_Title)

      }

      # Merge args and inject filtered data
      args <- args_base

      if(Auto_LD == TRUE)

      {

      filtered_data$R2_LD <- as.numeric(filtered_data$R2_LD)

      }

      message("Double check lead SNP for indexing in plot")

      lead_snp <- filtered_data %>%
        dplyr::filter(P == min(P, na.rm = TRUE)) %>%
        dplyr::pull(ID)  # or whatever your SNP ID column is called

      if(Auto_LD == TRUE)

      {

      message("Binning LD values for colouring")

      filtered_data$LD_Bin <- dplyr::case_when(
        filtered_data$SNP == lead_snp ~ "Lead",
        is.na(filtered_data$R2_LD) ~ "NA",
        filtered_data$R2_LD < 0.2 ~ "<0.2",
        filtered_data$R2_LD >= 0.2 & filtered_data$R2_LD < 0.4 ~ "0.2-0.4",
        filtered_data$R2_LD >= 0.4 & filtered_data$R2_LD < 0.6 ~ "0.4-0.6",
        filtered_data$R2_LD >= 0.6 & filtered_data$R2_LD < 0.8 ~ "0.6-0.8",
        filtered_data$R2_LD >= 0.8 & filtered_data$R2_LD < 1 ~ "0.8-<1",
        filtered_data$R2_LD == 1 ~ "1"
      )

      }

      if (1 == 1) {

        message("Creating hover text for interactive version")

        if(Auto_LD == TRUE)

        {

        filtered_data$Hover_Info <- paste0(
          "SNP: ", filtered_data$ID, "\n",
          "CHR: ", filtered_data$CHROM, "\n",
          "POS: ", filtered_data$GENPOS, "\n",
          "P: ", signif(filtered_data$P, 4), "\n",
          "REF: ", filtered_data$ALLELE0, "\n",
          "ALT: ", filtered_data$ALLELE1, "\n",
          "LD: ",  signif(filtered_data$R2_LD, 2)
        )

        }else{

          filtered_data$Hover_Info <- paste0(
            "SNP: ", filtered_data$ID, "\n",
            "CHR: ", filtered_data$CHROM, "\n",
            "POS: ", filtered_data$GENPOS, "\n",
            "P: ", signif(filtered_data$P, 4), "\n",
            "REF: ", filtered_data$ALLELE0, "\n",
            "ALT: ", filtered_data$ALLELE1)

        }

        #match single plot format
        #Double check if necessary when testing
     #   filtered_data$Hover_Info <- NA_character_
      }

      message("Indexing top variant")

      #need diamond for each as not always lead chrom peak
      filtered_data <- filtered_data %>%
        dplyr::mutate(top = (P == min(P, na.rm = TRUE)))

      pretty_peak_id <- gsub("peak(\\d+)", "Peak \\1", peak_id)

      Title_Name <- unique(filtered_data$Title_Name)

      message("Spacing Title")

      message(Title_Name)

  #    if (is.null(args$Title)) {

          args$Title <- paste0(Title_Name, "\n\n\n")

  #        }

      args <- modifyList(args, list(
        Data = filtered_data
      ))

      args$session <- NULL

      message("Assigning single colour")

      args$Chromosome_Colours <- Point_Colour

      if(Recombination_Line == TRUE)
      {

      message("Formatting recombination axes")

      }

      if(is.null(Recombination_Axis_Title_Size))
      {

        Recombination_Axis_Title_Size <- args$Y_Axis_Title_Size

      }
      if(is.null(Recombination_Axis_Text_Size))
      {

        Recombination_Axis_Text_Size <- args$Y_Axis_Text_Size

      }

      message("Calling inherited Single_Plot()")

      p <- do.call(.Single_Plot_original, args)

      message("Checking chromosome being processed")

      chrom <- unique(p$data$CHROM)

      message("Determining X axis optimal breaks and labels")

      pos_range <- range(filtered_data$GENPOS, na.rm = TRUE)
      breaks <- pretty(pos_range, n = 6)

      span_bp <- diff(pos_range)
      span_mb <- span_bp / 1e6

      decimal_places <- dplyr::case_when(
          span_mb < 0.1 ~ 3,
          span_mb < 1 ~ 2,
          span_mb < 10 ~ 1,
          TRUE ~ 0
        )

      labels <- paste0(round(breaks / 1e6, decimal_places), " Mb")

        # Apply pos_range as limits (optional)
        limits <- pos_range

        # Filter breaks and labels to be within limits
        valid_idx <- breaks >= limits[1] & breaks <= limits[2]
        breaks <- breaks[valid_idx]
        labels <- labels[valid_idx]

        if(Genome_Build == "grch38")
        {

           if(Gene_Tracks == TRUE)
          {

            message("Loading gene annotation data")

          }

          #setwd("C:/Users/callumon/Miami_Package_R/MiamiR/data/Processed_Gene_Data")
          #read_name <- paste0("Chromosome_", chrom, "_", "Gene_Data_HG38_Processed.json")
          #gene_data <- jsonlite::read_json(read_name, simplifyVector = TRUE)

          read_name <- paste0("Chromosome_", chrom, "_Gene_Data_HG38_Processed.json")
          json_path <- system.file("extdata", "Processed_Gene_Data", read_name, package = "MiamiR")
          gene_data <- jsonlite::read_json(json_path, simplifyVector = TRUE)


          if(Recombination_Line == TRUE)
          {

            message("Loading recombination annotation data")

          }

          #setwd("C:/Users/callumon/Miami_Package_R/MiamiR/data/Processed_Recombination_Data")
          #recomb_data <- read.csv("Processed_AVG_Recomb_HG38.csv")

          recomb_data <- read.csv(
            system.file("extdata", "Processed_Recombination_Data", "Processed_AVG_Recomb_HG38.csv", package = "MiamiR")
          )

        }

        if(Genome_Build == "grch37")
        {

          if(Gene_Tracks == TRUE)
          {

            message("Loading gene annotation data")

          }

          #setwd("C:/Users/callumon/Miami_Package_R/MiamiR/data/Processed_Gene_Data")
          #read_name <- paste0("Chromosome_", chrom, "_", "Gene_Data_HG19_Processed.json")
          #gene_data <- jsonlite::read_json(read_name, simplifyVector = TRUE)

          read_name <- paste0("Chromosome_", chrom, "_Gene_Data_HG19_Processed.json")
          json_path <- system.file("extdata", "Processed_Gene_Data", read_name, package = "MiamiR")
          gene_data <- jsonlite::read_json(json_path, simplifyVector = TRUE)

          if(Recombination_Line == TRUE)
          {

            message("Loading recombination annotation data")

          }

          #setwd("C:/Users/callumon/Miami_Package_R/MiamiR/data/Processed_Recombination_Data")
          #recomb_data <- read.csv("Processed_AVG_Recomb_HG19.csv")

          recomb_data <- read.csv(
            system.file("extdata", "Processed_Recombination_Data", "Processed_AVG_Recomb_HG19.csv", package = "MiamiR")
          )

        }

        message("Backing up existing P/Y axis scaling")

        idx <- which(sapply(p$scales$scales, function(s) "y" %in% s$aesthetics))

        if (length(idx) > 0) {

          yscale <- p$scales$scales[[idx]]

          trans_obj <- yscale$trans
          limits_obj <- yscale$limits
          breaks_obj <- yscale$breaks
          expand_obj <- yscale$expand

          actual_max <- max(p$data$log10pval, na.rm = TRUE)
          limits_obj[2] <- actual_max

          }
         else {

          trans_obj <- NULL
          limits_obj <- NULL
          breaks_obj <- NULL
          expand_obj <- NULL

        }

        suppressMessages(suppressWarnings({

        if(Recombination_Line == TRUE)

        {

          message("Adding recombination line based on midpoints")

          chrom <- unique(p$data$CHROM)
          genpos <- p$data$GENPOS
          xrange <- range(genpos, na.rm = TRUE)

          recomb_df <- dplyr::filter(
            recomb_data,
            CHROM == paste0("chr", chrom),
            end >= xrange[1],
            start <= xrange[2]
          ) %>%
            dplyr::mutate(
              midpoint = (start + end) / 2
            )

          message("Scaling to match P axis size")

          max_p <- max(p$data$log10pval, na.rm = TRUE)
          max_rr <- max(recomb_df$score, na.rm = TRUE)
          scale_factorz <- max_p / max_rr


          message("Lifting over recombination axes")

            p <- p +
              ggplot2::geom_line(
                data = recomb_df,
                mapping = ggplot2::aes(x = midpoint, y = score * scale_factorz),
                color = Recombination_Line_Colour,
                size = Recombination_Line_Thickness,
                linetype = Recombination_Line_Type
              ) +
              ggplot2::scale_y_continuous(
                name = args_real$Y_Axis_Title,   # Preserves original y-axis title
                trans = trans_obj,               # Preserve transformation (like ggmanh manhattan_scale)
                limits = limits_obj,             # Preserve limits
                breaks = breaks_obj,             # Preserve breaks
                expand = expand_obj,             # Preserve expansion
                sec.axis = ggplot2::sec_axis(
                  ~ . / scale_factorz,
                  name = Recombination_Axis_Title
                )
              )+
              ggplot2::theme(
                axis.title.y.right = ggplot2::element_text(color = "black", vjust = 7, size = Recombination_Axis_Title_Size ),
                axis.text.y.right =  ggplot2::element_text(size = Recombination_Axis_Text_Size, colour = "black"),
              )
        }


        }))

        #Might need to expand in future or have partial matching mechanism, especially between HG38 and HG19 files

        if (!is.null(Gene_Biotype)) {

          message("Pulling out key boptypes for gene annotation")

          # Alias table
          biotype_alias_table <- list(
            protein_coding = c("coding", "protein", "protein_coding", "proteincoding"),
            miRNA          = c("miRNA", "mirna", "mir", "microrna"),
            snoRNA         = c("snoRNA", "snorna", "sno"),
            snRNA          = c("snrna", "snr", "snRNA"),
            rRNA           = c("rRNA"),
            rRNA_pseudogene = c("rRNA_pseudogene"),
            scaRNA         = c("scarna", "scaRNA"),
            lncRNA         = c("lncrna", "longnoncoding","lncRNA"),
            processed_pseudogene = c("processedpseudo", "ppseudo","processed_pseudogene"),
            unprocessed_pseudogene = c("unprocessedpseudo", "upseudo","unprocessed_pseudogene"),
            transcribed_processed_pseudogene = c("transcribedprocessedpseudo", "transcribed_processed_pseudogene"),
            transcribed_unprocessed_pseudogene = c("transcribedunprocessedpseudo","transcribed_unprocessed_pseudogene"),
            transcribed_unitary_pseudogene = c("transcribedunitarypseudo","transcribed_unitary_pseudogene"),
            unitary_pseudogene = c("unitarypseudo","unitary_pseudogene"),
            misc_RNA = c("misc", "misc_rna","misc_RNA"),
            TEC = c("tec","TEC")
          )

          resolve_biotypes <- function(user_inputs, alias_table) {
            user_inputs <- tolower(trimws(user_inputs))
            matched_biotypes <- unique(unlist(
              lapply(user_inputs, function(term) {
                matches <- names(alias_table)[sapply(alias_table, function(aliases) term %in% aliases)]
                if (length(matches) > 0) return(matches)
                if (term %in% names(alias_table)) return(term)
                return(NULL)
              })
            ))
            return(matched_biotypes)
          }

          valid_biotypes <- resolve_biotypes(Gene_Biotype, biotype_alias_table)

          gene_data <- gene_data %>%
            dplyr::filter(gene_biotype %in% valid_biotypes)

          if (!is.null(Gene_Biotype_Colour)) {

            message("Custom colouring gene annotations")

            if (length(Gene_Biotype_Colour) == 1) {
              # If single color, apply to all
              gene_data$Gene_Biotype_Colour <- Gene_Biotype_Colour

            } else if (length(Gene_Biotype_Colour) >= length(valid_biotypes)) {
              # Multiple colors - assign based on order of valid_biotypes
              #Must match length of unique biotypes selected
              color_map <- setNames(Gene_Biotype_Colour[seq_along(valid_biotypes)], valid_biotypes)
              gene_data$Gene_Biotype_Colour <- color_map[gene_data$gene_biotype]

            } else {

              # Not enough colors provided - default single color
              warning("Not enough colors for Gene_Biotype; using first color for all.")

              gene_data$Gene_Biotype_Colour <- Gene_Biotype_Colour[1]
            }

          } else {

            # If Gene_Biotype_Colour not provided then assign default color
            gene_data$Gene_Biotype_Colour <- "black"

          }

        }

        if(!is.null(Gene_Name))

        {

         message("Pulling out key bodies by name for gene annotation")

         gene_data <- gene_data %>%
            dplyr::filter(gene_name %in% Gene_Name)

        }

        top <- max(pos_range, na.rm = TRUE)
        bottom <- min(pos_range, na.rm = TRUE)

        message("Filtering remaining gene data for range of regional plot")

        gene_data_top <- gene_data[gene_data$start <= top, , drop = FALSE]
        gene_data_filtered <- gene_data_top[gene_data_top$end >= bottom, , drop = FALSE]

        gene_data <- gene_data_filtered

        edge_buffer <- 0.2 * diff(pos_range)

        message("Configuring coordinates for neat alignment and presentation")

        gene_data <- gene_data %>%
          dplyr::mutate(
            visible_start = pmax(start, pos_range[1]),
            visible_end   = pmin(end,   pos_range[2]),
            visible_span  = visible_end - visible_start,
            view_width    = diff(pos_range),
            frac_visible         = visible_span / (end - start),
            frac_visible_of_view = visible_span / view_width,
            mid_visible = (visible_start + visible_end) / 2,

            label_nudge = pmin(0.0125 * view_width + 800 * nchar(label), 0.12 * view_width),

            label_x = dplyr::case_when(

              frac_visible_of_view > 0.05 ~ mid_visible,


              start < pos_range[1] + edge_buffer ~ visible_end + label_nudge,


              end > pos_range[2] - edge_buffer ~ visible_start - label_nudge,


              frac_visible_of_view > 0.5 ~ mid_visible,

              #Redundant, but keep for now

              # Clipped left
              # frac_visible < 0.5 & start < pos_range[1] ~ visible_end + label_nudge,

              # Clipped right
              # frac_visible < 0.5 & end > pos_range[2] ~ visible_start - label_nudge,

              TRUE ~ mid_visible
            ),

            label_hjust = dplyr::case_when(
              start < pos_range[1] + edge_buffer ~ 1,
              end > pos_range[2] - edge_buffer ~ 0,
              frac_visible_of_view > 0.5 ~ 0.5,
              frac_visible < 0.5 & start < pos_range[1] ~ 1,
              frac_visible < 0.5 & end > pos_range[2] ~ 0,
              TRUE ~ 0.5
            )
          )

        # Minimum gap (in base pairs) to allow genes on the same line
        min_gene_gap <- MIN_GENE_GAP

        # Sort by genomic start
        gene_data <- gene_data[order(gene_data$start), ]

        if(  nrow(gene_data) != 0)

        {

          message("Determining direction of transcription")

          #remove redundant debugging arrow directions
          gene_data$label <- gene_data$label |>
            gsub("^\\s*<-\\s*", "", x = _) |>
            gsub("\\s*(->|<-)\\s*$", "", x = _)

        # Initialize y-tier for stacking or annotations
        gene_data$y <- NA
        tiers <- list()

        message("algorithmically stacking gene annotation")

        for (i in seq_len(nrow(gene_data))) {

          gene_start <- gene_data$start[i]
          gene_end <- gene_data$end[i]

          placed <- FALSE

          for (tier_idx in seq_along(tiers)) {

                last_end <- tiers[[tier_idx]]

            if (gene_start > (last_end + min_gene_gap)) {

              gene_data$y[i] <- tier_idx * .75
              tiers[[tier_idx]] <- gene_end
              placed <- TRUE
              break

            }
          }

          if (!placed) {

            tier_number <- length(tiers) + 1
            gene_data$y[i] <- tier_number * .75
            tiers[[tier_number]] <- gene_end

          }
        }

        # Sort genes by position
        gene_data <- gene_data[order((gene_data$start + gene_data$end) / 2), ]

        }else{

          message("No bodies in region window, dispensing of Gene_Tracks panel")

          Gene_Tracks <- FALSE

        }

         p <- p + ggplot2::labs(x = NULL)



        suppressMessages(suppressWarnings({

        if(Auto_LD == TRUE & Interactive == TRUE)

        {

          message("Adjusting hover INFO for interactive mode")

        p <- p +  ggplot2::geom_point(
            data = filtered_data,
            ggplot2::aes(x = GENPOS, y = -log10(P), colour = LD_Bin,text = Hover_Info),
            size = args$Point_Size # gives mapping space under cursor same as point
          )


        }


        }))


      if(Auto_LD == TRUE)
      {

        message("Colouring Points by LD")

      }

        suppressMessages(suppressWarnings({

        if(Auto_LD == TRUE)
        {


          # Ensure LD_Bin is a factor with correct level ordering:
          ld_levels <- c("NA", "<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-<1", "Lead")
          p$data$LD_Bin <- factor(p$data$LD_Bin, levels = ld_levels)

          # Plot color scale
          plot_colors <- c(
            "NA" = "darkgrey",
            "<0.2" = "darkblue",
            "0.2-0.4" = "lightblue",
            "0.4-0.6" = "green",
            "0.6-0.8" = "yellow",
            "0.8-<1" = "red",
            "Lead" = "transparent"  # Plot Lead as invisible, will add later
          )

          legend_colors <- c(
            "NA" = "darkgrey",
            "<0.2" = "darkblue",
            "0.2-0.4" = "lightblue",
            "0.4-0.6" = "green",
            "0.6-0.8" = "yellow",
            "0.8-<1" = "red",
            "Lead" = Colour_Of_Diamond  # Legend Lead = black diamond
          )

          p <- p + ggplot2::aes(colour = LD_Bin)

          p <- p + ggplot2::scale_colour_manual(
            values = plot_colors,
            na.value = "darkgrey"
          ) +
            ggplot2::labs(colour = LD_Legend_Title)

          # Detect which categories are present:
          cats_present <- levels(droplevels(p$data$LD_Bin))

          # Shape override for legend:
          shape_map <- rep(16, length(cats_present))
          names(shape_map) <- cats_present

          if ("Lead" %in% cats_present) shape_map["Lead"] <- 18  # Diamond

          # Legend override:
          p <- p + ggplot2::guides(
            colour = ggplot2::guide_legend(
              override.aes = list(
                shape = unname(shape_map),
                size = 7,
                colour = legend_colors[cats_present]  # Explicitly map legend color, not plot color!
              )
            ),
            shape = "none",
            size = "none"
          )


        }


        }))

        if (LD_Legend_Location == "Top Right") {
          legend_pos <- c(1, 1)
          legend_just <- c(1, 1)
        } else if (LD_Legend_Location == "Top Left") {
          legend_pos <- c(0, 1)
          legend_just <- c(0, 1)
        } else if (LD_Legend_Location == "Bottom Right") {
          legend_pos <- c(1, 0)
          legend_just <- c(1, 0)
        } else if (LD_Legend_Location == "Bottom Left") {
          legend_pos <- c(0, 0)
          legend_just <- c(0, 0)
        } else if (LD_Legend_Location == "Centre Right") {
          legend_pos <- c(1, 0.5)
          legend_just <- c(1, 0.5)
        } else if (LD_Legend_Location == "Centre Left") {
          legend_pos <- c(0, 0.5)
          legend_just <- c(0, 0.5)
        } else if (LD_Legend_Location == "Top Centre") {
          legend_pos <- c(0.5, 1)
          legend_just <- c(0.5, 1)
        } else if (LD_Legend_Location == "Bottom Centre") {
          legend_pos <- c(0.5, 0)
          legend_just <- c(0.5, 0)
        } else if (LD_Legend_Location == "Right Outside") {
          legend_pos <- "right"
          legend_just <- NULL

        }

        #need separate logic or early trigger
        if (LD_Legend_On == FALSE) {
          legend_pos <- "none"  # restores ggplot2 default behavior outside plot
          legend_just <- NULL    # no justification needed for "right"

        }

        if (Auto_LD == TRUE) {

        message("Formatting LD Legend")

        p <- p +
          ggplot2::theme(
            legend.key.height = ggplot2::unit(1, "cm"),
            legend.position = legend_pos,
            legend.justification = legend_just,
            legend.box.just = "left",
            legend.spacing.y = ggplot2::unit(0.5, "cm"),  # Reduce vertical space
            legend.title = ggplot2::element_text(
              margin = ggplot2::margin(b = 30),
              size = LD_Legend_Title_Size,
              face = "bold",
              family = "sans"
            ),
            legend.text = ggplot2::element_text(
              size = LD_Legend_Text_Size,
              family = "sans",
              margin = ggplot2::margin(l = 8, b = 2)
            ),
            legend.key = ggplot2::element_blank(),
            legend.background = ggplot2::element_rect(
              colour = "black",   # Outline color
              linewidth = 0.5,    # Thickness of border
              fill = scales::alpha("white", 0.7)

            )
          )

}

        suppressMessages(suppressWarnings({


        if(Gene_Tracks == TRUE & Interactive == FALSE)
        {
          message("Formatting margins on main plot to prepare for gene panel join")

          p <- p + ggplot2::theme(axis.line.x = ggplot2::element_blank(), plot.margin =  ggplot2::margin(t = 75, b = 0, r = 75, l = 75),
                         axis.ticks.x = ggplot2::element_blank(),
                         axis.title.x =    ggplot2::element_blank(),
                         axis.text.x =    ggplot2::element_blank())
        }


        if(Gene_Tracks == TRUE & Interactive == TRUE)
        {

          message("Formatting margins on main plot, contextually for interactive mode to prepare for gene panel join")

          p <- p + ggplot2::theme(plot.margin =  ggplot2::margin(t = 75, b = 0, r = 75, l = 75),
                                  axis.ticks.x = ggplot2::element_blank(),
                                  axis.title.x =    ggplot2::element_blank(),
                                  axis.text.x =    ggplot2::element_blank())

        }



        }))

        suppressMessages(suppressWarnings({


        if(Gene_Tracks == FALSE)
        {

          message("Formatting margins on main plot, in absence of gene panel join")

        }

        if(Gene_Tracks == FALSE & (nrow(args$Data) > 1))
        {

            p <- p + ggplot2::scale_x_continuous(
              breaks = breaks,
              labels = labels,
              position = "bottom" ,
              expand = c(0,0)
            )

        }


        if(Gene_Tracks == FALSE )

        {

        p <- p +
          ggplot2::labs(x = args_real$X_Axis_Title)+ ggplot2::theme(
          plot.margin =    ggplot2::margin(t = 75, b = 75, r = 75, l = 75),
          axis.title.x = ggplot2::element_text(size = args_real$X_Axis_Title_Size, vjust = -7, colour = "black"),
          axis.text.x  = ggplot2::element_text(size = args_real$Chromosome_Label_Size, vjust = -2, colour = "black" )
          )


        }


      if(Auto_LD == FALSE)
      {

        message("No LD legend required")

        p <- p + ggplot2::theme(legend.position = "none")

      }

        }))

      p_main_plot <- p

      if( nrow(gene_data) != 0){

      message("Defining introns and exons")



      exons_df2 <-  suppressMessages( gene_data %>%
        dplyr::select(transcript_id, gene_id, tx_start, tx_end, tx_length, strand,gene_biotype,  label, y, exons) %>%
        tidyr::unnest(exons)
      )

      introns_df2 <-  suppressMessages( gene_data %>%
        dplyr::select(transcript_id, gene_id, tx_start, tx_end, tx_length, strand,gene_biotype,  label, y, introns) %>%
        tidyr::unnest(introns, names_repair = "unique"))

      introns_df2$start <-  introns_df2$intron_start
      introns_df2$end <-  introns_df2$intron_end

      gene_data$gene_biotype <- factor(gene_data$gene_biotype)

      levels(gene_data$gene_biotype) <- gsub("_", " ", levels(gene_data$gene_biotype))

      safely_strip_attributes <- function(x) {
        keep <- c("names", "row.names", "class")
        attr_list <- attributes(x)
        for (nm in setdiff(names(attr_list), keep)) {
          attr(x, nm) <- NULL
        }
        x
      }


      introns_df2 <- safely_strip_attributes(introns_df2)
      exons_df2 <- safely_strip_attributes(exons_df2)


      }







      suppressMessages(suppressWarnings({

        if( nrow(gene_data) != 0)
        {

          message("Detailing gene panel with exons")

          p_genes <- ggplot2::ggplot(gene_data) +
            ggplot2::geom_segment(data = exons_df2,
                                  ggplot2::aes(x = start, xend = end, y = y, yend = y),
                                  size = Gene_Structure_Size_Exons,
                                  color = Gene_Structure_Colour_Exons)


          n_gene_rows <- max(gene_data$y)

          # Dynamic scaling of total height

          scale_factor <- 0.6 * (Gene_Structure_Size_Exons / 7) + 0.4 * (Gene_Label_Size / 5)

          total_height <- (8 + (n_gene_rows * 2)) * scale_factor

          message("Considering strand and arrow")

          gene_data <- gene_data %>%
            dplyr::mutate(
              arrow_char = ifelse(strand == "+", "→", "←"),
              arrow_x = ifelse(strand == "+", end, start),
              arrow_y = (y + 0.013) * (scale_factor/1.8285) # 👈 Small upward nudge to center better
            )

          if(nrow(introns_df2) > 0)
          {

            message("Detailing gene panel with introns")

            p_genes <-  p_genes +  ggplot2::geom_segment(data = introns_df2,
                                                         ggplot2::aes(x = start, xend = end, y = y, yend = y),
                                                         color = Gene_Structure_Colour_Introns,
                                                         alpha = 0.1,
                                                         size = Gene_Structure_Size_Introns)



            if(Sense_Arrow == TRUE)
            {

              message("Finalising direction arrow characteristics")

              # Compute view width for arrow segment and gap
      view_width <- diff(pos_range)
      arrow_body_length <- Sense_Arrow_Body_Length * view_width     # Total arrow line length
      arrow_gap         <- Sense_Arrow_Gene_Gap * view_width    # Gap from gene edge

      threshold <- 10  # small buffer near plot edges

      plot_xmin <- pos_range[1]
      plot_xmax <- pos_range[2]

      threshold <- 0.01 * view_width

      strand_highlight_df <- gene_data %>%
        dplyr::mutate(
          near_right_edge = end + arrow_gap + arrow_body_length > (plot_xmax - threshold),
          near_left_edge  = start - arrow_gap - arrow_body_length < (plot_xmin + threshold),

          xstart = dplyr::case_when(
            strand == "+" & near_right_edge ~ start - arrow_gap - arrow_body_length,
            strand == "-" & near_left_edge  ~ end + arrow_gap + arrow_body_length,
            strand == "+" ~ end + arrow_gap,
            strand == "-" ~ start - arrow_gap
          ),
          xend = dplyr::case_when(
            strand == "+" & near_right_edge ~ start - arrow_gap,
            strand == "-" & near_left_edge  ~ end + arrow_gap,
            strand == "+" ~ end + arrow_gap + arrow_body_length,
            strand == "-" ~ start - arrow_gap - arrow_body_length
          )
        )


      if(is.null(Sense_Arrow_Body_Thickness))
      {

        message("Determining arrow thickness")
        Sense_Arrow_Body_Thickness <- Gene_Structure_Size_Introns * 0.6

      }


      # Draw arrowed segments with visible shaft and head
      message("Printing arrows")

      p_genes <- p_genes + ggarchery::geom_arrowsegment(
        data = strand_highlight_df,
        ggplot2::aes(x = xstart, xend = xend, y = y, yend = y),
        arrows = grid::arrow(
          type = "closed",
          length = grid::unit(Sense_Arrow_Head_Size, "pt")  # Adjust tip size if needed
        ),
        inherit.aes = FALSE,
        color = Sense_Arrow_Colour,
        fill = Sense_Arrow_Colour,
        size = Sense_Arrow_Body_Thickness
      )


      }

      }

      suppressMessages(suppressWarnings({


        if(!is.null(Gene_Biotype_Colour))
        {

          message("Adding custom gene body colours")

          valid_biotypes_pretty <- gsub("_", " ", valid_biotypes)

          gene_data$gene_biotype <- factor(
            gene_data$gene_biotype,
            levels = valid_biotypes_pretty
          )

          full_color_map <- setNames(Gene_Biotype_Colour, valid_biotypes_pretty)

          present_biotypes <- levels(droplevels(gene_data$gene_biotype))

          color_map_present <- full_color_map[present_biotypes]

          message("Adding bodies and formatting")

          p_genes <- p_genes +
            ggplot2::geom_text(
              ggplot2::aes(
                x = label_x, hjust = label_hjust,
                y = y + 0.375,
                label = label
              ),
              color = gene_data$Gene_Biotype_Colour,  # manual color assignment
              show.legend = FALSE,
              angle = 0,
              family = "sans",
              fontface = "plain",
              vjust = 0.5,
              hjust = 0.5,
              size = Gene_Label_Size
            ) +

            ggplot2::geom_point(
              data = unique(gene_data[c("gene_biotype", "start", "y")]),
              ggplot2::aes(x = start, y = y, color = gene_biotype),
              shape = 15, size = 0, alpha = 1, inherit.aes = FALSE, show.legend = TRUE
            ) +

            ggplot2::scale_color_manual(
              name = Gene_Legend_Title,
              values = color_map_present,
              limits = present_biotypes
            )+

            ggplot2::coord_cartesian(clip = "off") +

            ggplot2::guides(
              color = ggplot2::guide_legend(
                override.aes = list(
                  shape = 15,
                  size = 7,
                  stroke = 0,
                  vjust = 0.5
                )
              )
            )

        }else{

          message("Adding suitable gene body colours, bodies themselves and formatting")

          p_genes <-  p_genes +  ggplot2::geom_text(
            ggplot2::aes(
              x = label_x, hjust = label_hjust,
              y = y + 0.375,
              label = label ,
              color = gene_biotype
            ),
            show.legend = FALSE,
            angle = 0,
            family = "sans",
            fontface = "plain",
            vjust = 0.5,
            hjust = 0.5,
            size = Gene_Label_Size
          ) + ggplot2::scale_color_brewer(
            palette = "Paired", #wide enough to cover max selection
            name = Gene_Legend_Title,
            labels = levels(gene_data$gene_biotype)
            ) +
            ggplot2::coord_cartesian(clip = "off")+
            ggplot2::geom_point(
              data = unique(gene_data[c("gene_biotype", "start", "y")]),
              ggplot2::aes(x = start, y = y, color = gene_biotype),
              shape = 15, size = 0, alpha = 1, show.legend = TRUE, inherit.aes = FALSE
            )+
            ggplot2::guides(
              color = ggplot2::guide_legend(
                override.aes = list(
                  shape = 15,
                  size = 7, #size of squares/key
                  stroke = 0,
                  vjust = 0.5
                )
              )
            )

        }

}))


      message("Formatting gene panel")

      p_genes <- p_genes +
          ggplot2::scale_x_continuous(
          breaks = breaks,
          labels = labels,
          sec.axis = ggplot2::dup_axis( #create top ticks
            labels = labels,
            name = args_real$X_Axis_Title
          )
        )+ ggplot2::coord_cartesian(xlim = pos_range)

      p_genes <- p_genes +

        ggplot2::labs(x = args_real$X_Axis_Title, y = NULL) +
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = args_real$X_Axis_Title_Size, vjust = -7, colour = "black"),
          axis.text.x  = ggplot2::element_text(size = args_real$Chromosome_Label_Size, vjust = -2, colour = "black" ),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.title.x.top = ggplot2::element_blank(),
          axis.text.x.top = ggplot2::element_blank(),
          axis.ticks.length.x  = ggplot2::unit(0.3, "cm"), #match bottom and top ticks in opposite directions
          axis.ticks.length.x.top = grid::unit(-0.3, "cm"),
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(),
          plot.title = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(t = 0, r = 75, b = 75, l = 75),
          plot.title.position = "plot", #could allow for customisation here in future
          axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.5), #could also adjust
          legend.justification = "left",
          legend.title = ggplot2::element_text(
          margin = ggplot2::margin(b = 30), #to control space between legends.
          size = Gene_Legend_Title_Size,
          face = "bold",
          family = "sans"
        ),
         legend.text = ggplot2::element_text(
          size = Gene_Legend_Text_Size,
          vjust = 0,
          family = "sans",
          margin = ggplot2::margin(l = 8, b = 8)# match font for entries too
        ),
         legend.text.align = 0,
         legend.key = ggplot2::element_rect(fill = NA)  # ensures square color blocks
          )
         p_genes <- p_genes + ggplot2::theme(panel.background = ggplot2::element_rect(fill = Gene_Panel_Background_Colour, colour = NA))


      if(Interactive == FALSE)

      {

      message("Adding top and bottom buffering to gene panel area")

      p_genes <- p_genes + ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(add = c(0.3, 0.3))  # bottom, top 0.5, 0.5
        )

      }else{

        y_min <- 1.25 #smaller for bigger gap between
        y_max <- max(gene_data$y, na.rm = TRUE)

        pad <- 1

        y_lims <- c(y_min - pad, y_max + pad)

        p_genes <- p_genes + ggplot2::scale_y_continuous(
          limits = y_lims,
          expand = c(0, 0)  # turn off auto-padding
        )

      }


      if (Gene_Legend_Location == "Top Right") {
        legend_pos_gene <- c(1, 1)
        legend_just_gene <- c(1, 1)
      } else if (Gene_Legend_Location == "Top Left") {
        legend_pos_gene <- c(0, 1)
        legend_just_gene <- c(0, 1)
      } else if (Gene_Legend_Location == "Bottom Right") {
        legend_pos_gene <- c(1, 0)
        legend_just_gene <- c(1, 0)
      } else if (Gene_Legend_Location == "Bottom Left") {
        legend_pos_gene <- c(0, 0)
        legend_just_gene <- c(0, 0)
      } else if (Gene_Legend_Location == "Centre Right") {
        legend_pos_gene <- c(1, 0.5)
        legend_just_gene <- c(1, 0.5)
      } else if (Gene_Legend_Location == "Centre Left") {
        legend_pos_gene <- c(0, 0.5)
        legend_just_gene <- c(0, 0.5)
      } else if (Gene_Legend_Location == "Top Centre") {
        legend_pos_gene <- c(0.5, 1)
        legend_just_gene <- c(0.5, 1)
      } else if (Gene_Legend_Location == "Bottom Centre") {
        legend_pos_gene <- c(0.5, 0)
        legend_just_gene <- c(0.5, 0)
      } else if (Gene_Legend_Location == "Right Outside") {
        legend_pos_gene <- "right"
        legend_just_gene <- NULL

      }


      #need separate logic or early trigger
      if (Gene_Legend_On == FALSE) {

        message("Hiding gene legend")

        legend_pos_gene <- "none"
        legend_just_gene <- NULL  #now unnecessary

      }else{

        message("Formatting and positioning gene legend")

      }

      p_genes <- p_genes +
        ggplot2::theme(
          legend.box.margin = ggplot2::margin(t = 20, l = 20),
          legend.position = legend_pos_gene,
          legend.justification = legend_just_gene,
          legend.key = ggplot2::element_blank(),
          legend.background = ggplot2::element_rect(
            fill = scales::alpha("white", 0.7),
            colour = "black",   # Outline color
            linewidth = 0.5,    # Thickness of border
          )
        )

        }

        }))


        if (Auto_LD == TRUE) {

        if(LD_Legend_Location == "Top Right")

        {

          message("Positioning LD legend")

        p <- p +
          ggplot2::theme(
            legend.box.margin = ggplot2::margin(t = 40,  r = 40),
            legend.position = legend_pos,
            legend.justification = legend_just,
            legend.key = ggplot2::element_blank(),
            legend.background = ggplot2::element_rect(
              fill = scales::alpha("white", 0.7),
              colour = "black",   # Outline color
              linewidth = 0.5,    # Thickness of border   # Legend background fill (optional)
            )
          )

        }

        if(LD_Legend_Location == "Top Left")

        {
          message("Positioning LD legend")

          p <- p +
            ggplot2::theme(
              legend.box.margin = ggplot2::margin(t = 40,  l = 40),
              legend.position = legend_pos,
              legend.justification = legend_just,
              legend.key = ggplot2::element_blank(),
              legend.background = ggplot2::element_rect(
                fill = scales::alpha("white", 0.7),
                colour = "black",   # Outline color
                linewidth = 0.5,    # Thickness of border   # Legend background fill (optional)
              )
            )

        }

}

    if(Gene_Tracks == TRUE)
    {

      message("Adjusting plot proportions for gene panel size")

      # Base height for the top plot:
      base_top_height <- 3

      # Base height for the gene panel (minimum):
      base_bottom_height <- 1

      # Additional height per gene row:
      extra_per_gene_row <- 0.5

      # Compute dynamic height for gene panel:
      gene_panel_height <- base_bottom_height + extra_per_gene_row * max(n_gene_rows, 1)

      # Then build plot:
      plot <- cowplot::plot_grid(
        p, p_genes,
        ncol = 1,
        rel_heights = c(base_top_height, gene_panel_height),
        align = "v",
        axis = "lr"
      )

      plot <- (p / p_genes) +
        patchwork::plot_layout(heights = c(base_top_height, gene_panel_height))


    }


      if(Gene_Tracks == FALSE)
      {

        message("Adjusting unnecessary as gene panel absent here")

        plot <- p

      }


      if(Interactive == "L")

      {

        message("Adjusting plot objects for interactive mode before output")

        #need different ticks style - estimate looks good, adjust later

       tick_y_base <- max(gene_data$y, na.rm = TRUE)
       tick_df <- data.frame(
         x = breaks,
         y = tick_y_base + 1,
         yend = tick_y_base + 1 - 0.1
       )

       p_genes <- p_genes +
         ggplot2::geom_segment(
           data = tick_df,
           ggplot2::aes(x = x, xend = x, y = y, yend = yend),
           inherit.aes = FALSE,
           color = "black",
           linewidth = 0.2
         )

      }


    if(Gene_Tracks == TRUE)
    {

    n_gene_rows <- max(gene_data$y)

    # Dynamic scaling of total height
    total_height <- (5 + (n_gene_rows * 2)) * (Gene_Structure_Size_Exons/7) * (Gene_Label_Size/5)  # Base height + per-row bonus

    scale_factor <- 0.6 * (Gene_Structure_Size_Exons / 7) + 0.4 * (Gene_Label_Size / 5)
    total_height <- (6.5 + (n_gene_rows * 2)) * scale_factor

    attr(plot, "dynamic_height") <- total_height

    }else{

         attr(plot, "dynamic_height") <- 15

    }

      if (Interactive == 12) {

        message("Adding interactive save dimensions")

        # Always do this:
        attr(p, "main_ggplot") <- p
        attr(p, "gene_track_plot") <- p_genes
        attr(p, "dynamic_height") <- total_height
        plot <- p

      }

        return(invisible(plot))
    })


    message("Naming plot objects informatively")

    message(Title_Names)

    names(plots) <- Title_Names

    return(invisible(plots))

   }

  message("Calling plot_all_peak_regions algorithm on region(s) of interest")

  region_plots <- plot_all_peak_regions(region_df, user_args = user_args)

  message("Concatenating results")

  return(invisible(region_plots))

}

.Regional_Plot_original <- Regional_Plot

Regional_Plot <- function(..., session = NULL) {

  call_expr <- match.call()
  fn_formals <- formals(.Regional_Plot_original)
  valid_args <- names(fn_formals)
  dots <- list(...)

  single_plot_formals <- names(formals(.Single_Plot_original))

  is_list_of_c <- function(x) {

    is.list(x) && all(vapply(x, function(e) is.call(e) && identical(e[[1]], as.name("c")), logical(1)))

  }

  grouped_args_detected <- any(vapply(dots, is_list_of_c, logical(1)))

  if (grouped_args_detected) {

    first_grouped <- which(vapply(dots, is_list_of_c, logical(1)))[1]
    n_groups <- length(dots[[first_grouped]])

    plots <- lapply(seq_len(n_groups), function(i) {
      args_i <- list()

      for (argname in names(dots)) {

        argval <- dots[[argname]]

        if (is_list_of_c(argval)) {

          c_expr <- argval[[i]]
          evaluated <- lapply(as.list(c_expr[-1]), function(e) eval(e, parent.frame()))
          args_i[[argname]] <- evaluated
        } else if (is.list(argval) && length(argval) == n_groups) {

          args_i[[argname]] <- argval[[i]]

        } else {

          args_i[[argname]] <- argval

        }
      }

      single_plot_args <- args_i[setdiff(names(args_i), valid_args)]
      regional_args <- args_i[intersect(names(args_i), valid_args)]
      regional_args$.dots <- single_plot_args

      for (arg in setdiff(valid_args, names(regional_args))) {

        default_val <- fn_formals[[arg]]
        if (!identical(default_val, quote(expr = ))) {

          val <- tryCatch(eval(default_val, envir = environment(.Regional_Plot_original)), error = function(e) NULL)
          regional_args[[arg]] <- val

        }

      }

      regional_args$session <- session

      if (isTRUE(regional_args$Verbose)) {

        return(invisible(do.call(.Regional_Plot_original, regional_args)))

      }

      run_with_counter(.Regional_Plot_original, args = regional_args, session = session)

    })

    return(invisible(plots))
  }

  raw_data_expr <- call_expr[["Data"]]

  if (!is.null(raw_data_expr) && is.call(raw_data_expr) && identical(raw_data_expr[[1]], as.name("c"))) {

    arg_exprs <- as.list(raw_data_expr[-1])
    arg_vals <- lapply(arg_exprs, function(e) eval(e, parent.frame()))

    if (all(vapply(arg_vals, is.data.frame, logical(1)))) {

      names(arg_vals) <- vapply(arg_exprs, function(e) deparse(e), character(1))
      dots$Data <- arg_vals

    } else if (all(vapply(arg_vals, is.character, logical(1)))) {

      if (all(file.exists(unlist(arg_vals)))) {

        data_list <- lapply(unlist(arg_vals), function(f) vroom::vroom(f, show_col_types = FALSE))
        names(data_list) <- basename(tools::file_path_sans_ext(unlist(arg_vals)))
        dots$Data <- data_list

      } else {

        stop("One or more file paths in c(...) do not exist.", call. = FALSE)

      }
    } else {
      stop("Mixed or invalid types in c(...) for Data.", call. = FALSE)
    }
  }

  Data <- dots$Data
  is_df <- function(x) is.data.frame(x)

  if (is.character(Data) && length(Data) > 1 && all(file.exists(Data))) {

    Data_list <- lapply(Data, function(f) vroom::vroom(f, show_col_types = FALSE))
    names(Data_list) <- basename(tools::file_path_sans_ext(Data))
    dots$Data <- Data_list
    Data <- Data_list

  }

  if (is.list(Data) && all(vapply(Data, is_df, logical(1)))) {

    if (is.null(names(Data)) || any(names(Data) == "")) {

      inferred_names <- vapply(Data, function(df) {
        if (!is.null(attr(df, "name"))) attr(df, "name")
        else paste0("Dataset ", as.character(sample(1000:9999, 1)))
      }, character(1))
      names(Data) <- inferred_names

    }

    n_data <- length(Data)

    plots <- lapply(seq_along(Data), function(i) {

      df <- Data[[i]]
      dataset_name <- names(Data)[i]

      args_i <- list()
      for (arg in names(dots)) {
        if (arg == "Data") next
        val <- dots[[arg]]
        if (length(val) == 1) {
          args_i[[arg]] <- val
        } else if (length(val) == n_data) {
          args_i[[arg]] <- val[[i]]
        } else {
          args_i[[arg]] <- val[[1]]
        }
      }

      args_i$Data <- df

      attr(args_i$Data, "RegionalPlot_Title") <- dataset_name

      single_plot_args <- args_i[setdiff(names(args_i), valid_args)]
      regional_args <- args_i[intersect(names(args_i), valid_args)]
      regional_args$.dots <- single_plot_args

      for (arg in setdiff(valid_args, names(regional_args))) {

        default_val <- fn_formals[[arg]]
        if (!identical(default_val, quote(expr = ))) {
          val <- tryCatch(eval(default_val, envir = environment(.Regional_Plot_original)), error = function(e) NULL)
          regional_args[[arg]] <- val

        }
      }

      regional_args$session <- session

      if (isTRUE(regional_args$Verbose)) {

        message(sprintf("Processing dataset: \"%s\"", dataset_name))
        return(invisible(do.call(.Regional_Plot_original, regional_args)))

          }

      message(sprintf("Processing dataset: \"%s\"", dataset_name))
      run_with_counter(.Regional_Plot_original, args = regional_args, session = session)
    })

    return(invisible(plots))
  }

  # SINGLE dataset fallback:
  single_plot_args <- dots[setdiff(names(dots), valid_args)]
  regional_args <- dots[intersect(names(dots), valid_args)]
  regional_args$.dots <- single_plot_args

  for (arg in setdiff(valid_args, names(regional_args))) {

    default_val <- fn_formals[[arg]]
    if (!identical(default_val, quote(expr = ))) {

      val <- tryCatch(eval(default_val, envir = environment(.Regional_Plot_original)), error = function(e) NULL)
      regional_args[[arg]] <- val

    }

  }

  if (!is.null(raw_data_expr)) {

    if (!is.call(raw_data_expr)) {

      dataset_name <- deparse(raw_data_expr)
      if (is.data.frame(dots$Data)) {

        attr(regional_args$Data, "RegionalPlot_Title") <- dataset_name
        message(sprintf("Processing dataset: %s", dataset_name))
      } else if (is.character(dots$Data) && length(dots$Data) == 1) {

        file_path <- dots$Data
        short_name <- tools::file_path_sans_ext(basename(file_path))
        attr(regional_args$Data, "RegionalPlot_Title") <- short_name
        message(sprintf("Processing dataset: %s", short_name))

      }
    } else if (is.call(raw_data_expr) && identical(raw_data_expr[[1]], as.name("["))) {

      dataset_val <- dots$Data
      if (is.character(dataset_val) && length(dataset_val) == 1) {

        short_name <- tools::file_path_sans_ext(basename(dataset_val))
        attr(regional_args$Data, "RegionalPlot_Title") <- short_name
        message(sprintf("Processing dataset: %s", short_name))
      } else if (is.data.frame(dataset_val)) {

        message("Processing dataset: <dataframe>")
      }
    }
  }

  regional_args$session <- session

  if (isTRUE(regional_args$Verbose)) {

    return(invisible(do.call(.Regional_Plot_original, regional_args)))
  }

  run_with_counter(.Regional_Plot_original, args = regional_args, session = session)
}
