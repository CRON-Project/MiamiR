
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
#' @param Y_Axis_Title Manually specify Y axis title; defaults to "-log\u2081\u2080(P)" - shows as "-log₁₀P"
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
#' @param Point_Size Size of data points; defaults to NULL and is modified depending on inherited functionality, unless manually assigned
#' @param Label_Index Annotate the index SNPs with ID provided; defaults to TRUE
#' @param Label_Size Size of index labels if Label_Index is TRUE; defaults to 5 but is automatically scaled in inherited functions
#' @param Label_Angle Angle of index labels if Label_Index is TRUE; defaults to 45 (degrees)
#' @param Label_Colour Colour of index labels if Label_Index is TRUE; defaults to "black"
#' @param Colour_Index Highlight the index SNP with a circular ring around the original point; defaults to TRUE
#' @param Colour_Of_Index Colour of index SNP circular highlight if Colour_Index is TRUE; defaults to "darkred"
#' @param Chromosome_Labels List of chromosomes to show index SNP labels on, given that Label_Index is TRUE; defaults to c(1:22, "X", "Y", "M")
#' @param Chromosome_Index  List of chromosomes to show index SNP circular highlight/colour on, given that Colour_Index is TRUE; defaults to c(1:22, "X", "Y", "M")
#' @param Condense_Scale Log transform and condense scale above allotted break point to enhance visualisation of 'hits' with vastly different test statistics; defaults to TRUE
#' @param Break_Point Value of P where y axis begins to condense; defaults to 1e-10, guided by frequent max value observed in a small sample GWAS
#' @param Draft_Plot Specify whether a random subset of SNPs, the number of which is specified by Random_Selection, is plotted; defaults to FALSE
#' @param Random_Selection Manually select a random subset of SNPs to plot (more for internal debugging and testing), given that Draft_Plot is TRUE; defaults to 4,000
#' @param Chromosome_Label_Drops Manually specify which chromosome/X axis labels (as a list) should not be shown to prevent text crowding; defaults to NULL leading to c(21,22) being passed when X chromosome data is present and c(21) when only autosomes are present.
#' @param Index_Size Size of the outer ring diameter used to highlight the index SNPs, when Colour_Index is TRUE; defaults to NULL and if not manually passed is always scaled to 2X the passed or default Point_Size
#' @param Index_Thickness Thickness of the outer ring circumference used to highlight the index SNPs, when Colour_Index is TRUE; defaults to 1 but is automatically scaled in inherited functions
#' @param Chromosome_Diamond List of chromosomes to show index SNP highlight/swapped diamond shape on, given that Diamond_Index is TRUE; defaults to c(1:22, "X", "Y", "M")
#' @param Anchor_Label Which side of the point to anchor/position the index label to, when Label_Index is TRUE; defaults to "left"
#' @param Label_Height Spacing between the point and anchoring position of the index SNP label and the beginning of the label print area; defaults to NULL but is assigned and otherwise automatically scaled in inherited functions
#' @param Diamond_Index Highlight the index SNP by swapping the standard point shape to a (different by default) coloured diamond; defaults to FALSE, unless called within Regional_Plot, then defaults to TRUE
#' @param Colour_Of_Diamond Colour of diamond shape swapped for standard point shape ar index SNPs, when Diamond_Index is TRUE; defaults to "purple"
#' @param Diamond_Index_Size Size of the index SNP swapped diamond shape used to highlight the index SNPs, when Diamond_Index is TRUE; defaults to NULL and if not passed is scaled to 2X the according default or passed Point_Size
#' @param Lab Specify an additional info column with information to be labelled regarding Index SNPs, when Label_Index is TRUE; defaults to NULL leading to the values in the SNP ID column being used (Coords/RS codes) by default
#' @param Verbose Prevent display of progress bar as function is running and instead show key milestone outputs/messages (mainly for debugging purposes); defaults to FALSE
#' @param Auto_Lab Allow Annotate_Data() to run on the data frame(s) provided so Lab can have RSIDs assigned, if only coordinates exist; defaults to FALSE
#' @param Genome_Build Reference genome to base required RSID annotations around; defaults to grch38
#' @param Chromosomes Manually specify which chromosomes from the summary statistics to run the function on; defaults to NULL, leading to all being included/absent filtering
#' @param Lower_Mult Number of units of coordinate buffering space between bottom of plot area and X axis; defaults to NULL, leading to automatic assignment unless specified
#' @param X_Axis_Title_Vjust Number of units of blank space between X_Axis_Title element area and above axis; defaults to -5
#' @param Y_Axis_Title_Vjust Number of units of blank space between Y_Axis_Title element area and adjacent axis; defaults to 5
#' @param X_Axis_Text_Vjust Number of units of blank space between X_Axis_Text element area and above axis; defaults to -0.5
#' @param X_Axis_Title_Drop_Spaces Number of lines of blank space between X_Axis_Title element area and above axis; defaults to 0, allowing for a tight fit
#' @param Title_On Toggle on/off whether Title (either manually or automatically determined) displays on final plot; defaults to TRUE
#' @param Title_Drop_Spaces Number of lines of blank space between Title element area and plot below; defaults to 0, allowing for a tight fit
#' @param Diamond_Index_SNP Manually specify a lead/index SNP of interest to base diamond styling on, rather than automatically allocated local maximum; defaults to NULL
#' @param Top_Expand Number of units of coordinate blank spacing between end of top of Y axis line and entire plot area; defaults to NULL, and requires a Y value to engage - particularly useful for Single_With_Regional() for additional spacing of insets
#' @param Upper_Mult Number of units of coordinate buffering space between top of plot area and end of Y axis line; defaults to NULL, leading to automatic assignment unless specified
#' @param Right_Mult Number of units of coordinate buffering space between right of plot area and end of X axis line; defaults to NULL, leading to automatic assignment unless specified
#' @param Left_Mult Number of units of coordinate buffering space between left of plot area and Y axis line; defaults to NULL, leading to automatic assignment unless specified
#' @param Y_Axis_Text_Vjust Number of units of blank space between Y_Axis_Text element area and adjacent axis; defaults to 0.5
#' @param Y_Axis_Title_Drop_Spaces Number of lines of blank space between Y_Axis_Title element area and adjacent axis; defaults to 0, allowing for a tight fit
#' @param Interactive Create plot amenable to interactive viewing - leads to slight modifications in elements used for compatibility and primarily calling in Shiny App; defaults to FALSE
#'
#' @return Image of Manhattan Plot(s) is allocated to specified object and the resulting ggplot object can then be saved to an image
#' @export
#'
#' @examples  Manhattan_Plot <- Single_Plot(Data = Intelligence_Sum_Stats)
#'

Single_Plot<- function(Data = NULL,
                       Random_Selection = 4000,
                       Draft_Plot = FALSE,
                       X_Axis_Title_Vjust = -5,
                       Y_Axis_Title_Vjust = 5,
                       X_Axis_Text_Vjust = -0.5,
                       Y_Axis_Text_Vjust = 0.5,
                       X_Axis_Title_Drop_Spaces = 0,
                       Y_Axis_Title_Drop_Spaces = 0,
                       Title = NULL,  # Null generally is an adjustment for Region or inherited
                       Title_Drop_Spaces = 1,
                       Title_On = TRUE,
                       Title_Size = 35,
                       X_Axis_Title_Size = 35,
                       Y_Axis_Title_Size = 35,
                       X_Axis_Title = NULL,
                       Y_Axis_Title = "-log\u2081\u2080(P)",
                       Chromosome_Label_Size = 30, # X
                       Y_Axis_Text_Size = 30,
                       Chromosomes = NULL,
                       Chromosome_Label_Drops = NULL,
                       Chromosome_Colours = c("blue", "turquoise"),
                       Sig_Line_Colour = "red",
                       Sig_Line_Type = "dashed",
                       Sig_Threshold = 5e-8,
                       Sig_Line_Width = 0.5,
                       Point_Size = NULL,
                       Chromosome_Labels = c(1:22, "X", "Y", "M"),
                       Chromosome_Index = c(1:22, "X", "Y", "M"),
                       Chromosome_Diamond = c(1:22, "X", "Y", "M"),
                       Label_Index = TRUE,
                       Label_Size = 5,
                       Label_Angle = 45,
                       Label_Colour = "black",
                       Label_Height = NULL,
                       Anchor_Label = "left",
                       Colour_Index = TRUE,
                       Colour_Of_Index = "darkred",
                       Index_Size = NULL,
                       Index_Thickness = 1,
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
                       Lab = NULL,
                       Auto_Lab = FALSE,
                       Diamond_Index_SNP = NULL,
                       Genome_Build = "grch38",
                       Right_Mult = NULL,
                       Left_Mult = NULL,
                       Lower_Mult = NULL,
                       Upper_Mult = NULL,
                       Top_Expand = NULL,
                       Interactive = FALSE,
                       Verbose = FALSE)

{

  # Pause progress bar, print a transient status line, run expr, then resume on the next line

  .progress_break <- function(expr, msg, pad = 1L) {

    .progress_pause()

    # Safely end the current progress line

    cat("\r\n", sep = "")
    cat(strrep("\n", pad), msg, strrep("\n\n", pad), "\n", sep = "")
    utils::flush.console()

    on.exit(.progress_resume(), add = TRUE)

    force(expr)

  }

  # progress hooks stay safe if runner isn't used

  if (!exists(".progress", inherits = TRUE))       .progress        <- function(...) invisible(NULL)

  if (!exists(".progress_done", inherits = TRUE))  .progress_done   <- function(...) invisible(NULL)

  if (!exists(".progress_pause", inherits = TRUE)) .progress_pause  <- function(...) invisible(NULL)

  if (!exists(".progress_resume", inherits = TRUE)).progress_resume <- function(...) invisible(NULL)

  # hide messages by default

  .allow_msgs <- isTRUE(Verbose)

  message <- function(..., domain = NULL, appendLF = TRUE) {

    if (isTRUE(get0(".allow_msgs", ifnotfound = FALSE, inherits = TRUE))) {

      base::message(..., domain = domain, appendLF = appendLF)
    }

    invisible(NULL)

  }

  # Silence nested base::message() calls when Verbose = FALSE

  .quiet_call <- function(fun, ..., verbose = isTRUE(Verbose)) {

    if (verbose) return(fun(...))

    withCallingHandlers(
      fun(...),
      message = function(m) invokeRestart("muffleMessage")

    )

  }

  message("Arguments received, beginning Single_Plot() Script")

  if (is.null(Data)) {

    stop("Please provide data.", call. = FALSE)

  }

  is_called_by_regional_plot <- any(vapply(sys.calls(), function(x) {

    "Regional_Plot" %in% deparse(x[[1]])

  }, logical(1)))

  message("Determining how many chromosomes in the Data...")

  Chromosome_Column <- .quiet_call(detect_chromosome_column, Data, Chromosome_Column)

  if (!is.null(Chromosomes)) {

    message("Filtering for specific chromosomes supplied")

    Data <- Data[ Data[[Chromosome_Column]] %in% Chromosomes, , drop = FALSE ]

  }

  if (is.null(X_Axis_Title) & is_called_by_regional_plot == FALSE & length(unique(Data[[Chromosome_Column]])) > 1) {

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

  if (is_called_by_regional_plot && missing(Point_Size))

  {

    Point_Size <- 5

  }

  if (!is_called_by_regional_plot && missing(Point_Size))

  {

    Point_Size <- 2.5

  }

  if(is.null(Point_Size))

  {

    Point_Size <- 2.5

  }

  if (is_called_by_regional_plot && missing(Label_Height))

  {

    Label_Height <- 20

  }

  if (!is_called_by_regional_plot && missing(Label_Height))

  {

    Label_Height <- ( (7 * Point_Size/2.5) * 0.8 )

  }

  if (is.null(Label_Height)) {

    Label_Height <- ( (7 * Point_Size/2.5) * 0.8 )

  }

  if (is_called_by_regional_plot && missing(Diamond_Index) )

  {

    Diamond_Index <- TRUE

  }

  if (is_called_by_regional_plot && missing(Label_Size))

  {

    Label_Size <- 9

  }

  if (!is_called_by_regional_plot && missing(Label_Size))

  {

    Label_Size <- 25

  }

  if (is_called_by_regional_plot && missing(Index_Thickness))

  {

    Index_Thickness <- 4

  }

  if (is_called_by_regional_plot && missing(Diamond_Index_Size))

  {

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

      suppressMessages(

        suppressWarnings(

     Data <- vroom::vroom(file_path, show_col_types = FALSE, progress = F)

        ))

     message("Finished reading")

    } else {

      stop("The provided character string does not point to an existing file: ", file_path, call. = FALSE)

         }
    }

  if(Draft_Plot == T & (Random_Selection < nrow(Data)))

  {

    message(paste0("Generating random subset of ", Random_Selection, " ", "points, specified as Draft_Plot activated"))

    Data <- dplyr::sample_n(Data, Random_Selection)

  }

  message("Deducing other key column names:")

  PValue_Column     <- .quiet_call(detect_pvalue_column, Data, PValue_Column)
  Position_Column   <- .quiet_call(detect_position_column, Data, Position_Column)
  SNP_ID_Column     <- .quiet_call(detect_snp_column, Data, SNP_ID_Column)
  Ref_Allele_Column <- .quiet_call(detect_reference_allele_column, Data, Reference_Allele_Column)
  Alt_Allele_Column <- .quiet_call(detect_effect_allele_column, Data, Effect_Allele_Column)

  if (!is.null(PValue_Column) && grepl("log", PValue_Column, ignore.case = TRUE) & !is_called_by_regional_plot ) { # miami doesn't do

    message("Converting log value to plain P")

    x <- Data[[PValue_Column]]

    x <- if (is.numeric(x)) x else as.numeric(x)

    Data[[PValue_Column]] <- exp(-x * 2.302585092994046)  # 2.30258... = log(10)

  }

  before <- ncol(Data)

  message("Pruning down to key data")

  # Need to keep derivatives for regional

  if(!is_called_by_regional_plot & !is_called_by_miami_plot) # messes up buffer point

  {

  keep <- c(Chromosome_Column, PValue_Column, Position_Column, SNP_ID_Column,
            Ref_Allele_Column, Alt_Allele_Column, "Lab")

  Data[ setdiff(names(Data), keep) ] <- NULL

  after <- ncol(Data)

  message(sprintf("Columns before: %d", before))

  message(sprintf("Columns now:    %d  (%d removed)", after, before - after))

  }

  message("Addressing variable chromosome X nomenclature")

  # Use the detected column

  col <- Chromosome_Column

  # Pull once, mutate locally, assign back once

  v <- Data[[col]]

  if (is.factor(v)) {

    lv <- levels(v)
    hit <- (lv == "23") | (lv == 23L)

    if (any(hit)) levels(v)[hit] <- "X"

    Data[[col]] <- v

  } else {

    vc <- as.character(v)

    # map 23/23L -> X

    idx <- (vc == "23") | (vc == "23L")

    if (any(idx)) vc[idx] <- "X"

    Data[[col]] <- vc

  }

  if (is.null(X_Axis_Title) && (length(unique(Data[[Chromosome_Column]])) == 1)) {

    message("Adjusting default X axis labelling for single choromosome plot")

    X_Axis_Title <- "Genomic Position (Mb)"

  }

  message("Ensure uppercase alleles")

  cols <- c(Ref_Allele_Column, Alt_Allele_Column)

  for (j in cols) {

    v <- Data[[j]]

    if (is.factor(v)) {

      lv <- levels(v); lv2 <- chartr("atgc","ATGC", lv)

      if (!identical(lv, lv2)) levels(v) <- lv2

      Data[[j]] <- v

    } else {

      vc <- as.character(v)

      if (any(grepl("[atgc]", vc, perl = TRUE))) {

        Data[[j]] <- chartr("atgc","ATGC", vc)

      }

    }

  }

  message("Forging missing IDs")

  id_col    <- SNP_ID_Column
  chrom_col <- Chromosome_Column
  pos_col   <- Position_Column
  ref_col   <- Ref_Allele_Column
  alt_col   <- Alt_Allele_Column

  x <- Data[[id_col]]

  empty_id <- is.na(x) | x == "." | x == ""

  # Only really need top indexes, but keep extra for regionals I suppose

  empty_id <- ( (is.na(x) | x == "." | x == "") & (Data[[PValue_Column]] < 0.0001) ) #5e-8)

  if (!any(empty_id, na.rm = TRUE)) {

    message("No missing/placeholder IDs; nothing to forge")

  } else {

    idx <- which(empty_id); n <- length(idx)
    chunk <- max(1L, floor(n / 20))

    .progress_pause()
    .allow_msgs <- TRUE

    message("Updating ", n, " missing IDs...")

    for (i in seq(1, n, by = chunk)) {

      j  <- seq.int(i, min(i + chunk - 1L, n))
      ii <- idx[j]

      # coerce once to character to avoid factor/integer pitfalls

      chr <- as.character(Data[[chrom_col]][ii])
      pos <- as.character(Data[[pos_col]][ii])
      ref <- as.character(Data[[ref_col]][ii])
      alt <- as.character(Data[[alt_col]][ii])

      Data[[id_col]][ii] <- paste0("chr", chr, ":", pos, ":", ref, ":", alt)

      message(sprintf("  %.0f%% done", 100 * (i + length(j) - 1) / n))

    }

    message("Finished updating IDs.")

    .allow_msgs <- FALSE
    .progress_resume()

  }

  # FAST: set genomic order without sorting a big integer vector

  message("Ordering chromosome information 1..22,X,Y,M (fast)")

  ch <- Chromosome_Column
  vals <- as.character(Data[[ch]])

  desired <- c(as.character(1:22), "X", "Y", "M")       # canonical order
  present <- unique(vals)

  # keep only levels that actually appear, in canonical order

  lvl_order <- desired[desired %in% present]

  message("Applying chromosome factoring")

  Data[[ch]] <- factor(vals, levels = lvl_order, ordered = TRUE)

  message("Chromosomes ordered")

  if (is.null(Chromosome_Label_Drops)){

    message("Assigning auto Chromosome_Label_Drops")

    if (is.null(Chromosome_Label_Drops) &&

        all(!(unique(Data[[ch]]) %in% c("X", "x", 23, "chr23", "CHRX", "CHR23", "CHRx")))) {

      Chromosome_Label_Drops <- c(21)

    } else {

      Chromosome_Label_Drops <- c(21, 22, "Y")

    }

  }

  max_logp <- max(-log10(Data[[PValue_Column]]), na.rm = TRUE)

  # fake anchoring around single point to prevent scale missing for single point plot/regional via inheritance

  if(!is_called_by_regional_plot) # doesn't need to be available there as of now

  {

    if (!is.null(Top_Expand)) {

    pval_col <- PValue_Column

    # take the first row as a template

    fake_top <- Data[1, , drop = FALSE]

    # set P so that -log10(P) = 70

    fake_top[[pval_col]] <- 10^(-Top_Expand)

    # give it a distinct SNP id (optional but recommended)

    fake_top[[SNP_ID_Column]] <- "FAKETOP70" # doesn't matter initial name

    Data_Orig <- Data

    # append

    Data <- dplyr::bind_rows(Data, fake_top)

  }

  }

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
    fake_minus[[SNP_ID_Column]] <- "FAKEMINUS1"

    # Create a copy with position + 1

    fake_plus <- original_row
    fake_plus[[pos_col]] <- fake_plus[[pos_col]] + 1
    fake_plus[[pval_col]] <- 0.0000000000001
    fake_plus[[SNP_ID_Column]] <- "FAKEPLUS1"

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

  if (!("Lab" %in% names(Data))) {

    message("Discerning Lead SNP Per Chromosome to annotate (fast, direct P column)")

    g   <- Data[[Chromosome_Column]]
    P   <- Data[[PValue_Column]]
    POS <- Data[[Position_Column]]
    ID  <- Data[[SNP_ID_Column]]

    # Indexing chromosomes (first-seen order)

    gfac <- factor(g, levels = unique(g), ordered = TRUE)
    gi   <- as.integer(gfac)

    # Preparing P for min-search

    P2 <- P
    P2[is.na(P2)] <- Inf

    # Ordering once and selecting per-chromosome minimums

    ord <- order(gi, P2)
    lead_row <- ord[!duplicated(gi[ord])]   # one row index per chromosome

    # Assigning lead positions back to all rows

    K <- max(gi)
    min_pos_chr <- rep_len(NA_real_, K)
    min_pos_chr[ gi[lead_row] ] <- POS[lead_row]
    min_pos_vec <- min_pos_chr[gi]

    # Broadcast lead positions back to all rows (fast integer indexing)

    Data$min_P_GENPOS <- min_pos_chr[gi]

    # Label only the lead rows that pass threshold — no POS==min_pos_vec scan

    thr <- 5e-8
    lab <- rep_len("", nrow(Data))

    lead_ok <- P[lead_row] < thr

    if (any(lead_ok)) {

      # be robust to factor IDs

      id_chr <- if (is.factor(ID)) as.character(ID) else ID

      lab[ lead_row[lead_ok] ] <- id_chr[ lead_row[lead_ok] ]

    }

    Data$Lab <- lab

    message("Lead SNP annotation complete")

  }

  if(Diamond_Index == TRUE)

  {

    message("Assigning Diamond Indices")

    # top already from regional sometimes

    if (!("top" %in% colnames(Data))) {

      chr_col <- Chromosome_Column
      pos_col <- Position_Column
      p_col   <- if (!is.null(PValue_Column)) PValue_Column else "P"

      Data <- Data %>%
        dplyr::group_by(.data[[chr_col]]) %>%
        dplyr::mutate(
          min_P_GENPOS = .data[[pos_col]][ which.min(.data[[p_col]]) ],
          top = (.data[[pos_col]] == min_P_GENPOS)
        ) %>%
        dplyr::ungroup()

    }

    if(!is.null(Diamond_Index_SNP))

    {

    Data$top <- FALSE
    Data$top <- as.character(Data[[SNP_ID_Column]]) %in% as.character(Diamond_Index_SNP)

    }

  # Backup df to draw from later

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

  Data$Lab[!(Data[[Chromosome_Column]] %in% Chromosome_Labels)] <- NA

  message("Filtering Index Circle Highlights")

  Data$COLOUR[!(Data[[Chromosome_Column]] %in% Chromosome_Index)] <- NA

  message("Formatting title spacing beneath")

  if (length(unique(Data[[Chromosome_Column]])) == 1 & !is_called_by_regional_plot) {

    Chrom_Avail <- unique(Data[[Chromosome_Column]])
    Title <- paste0(Title, " - ", "Chromosome ", Chrom_Avail)

  }

  # e.g. Title_Drop_Spaces = 2  -> "\n \n "

  Title_Drop_Spaces <- as.integer(Title_Drop_Spaces)

  if (is.na(Title_Drop_Spaces) || Title_Drop_Spaces < 0) Title_Drop_Spaces <- 2

  Title <- paste0(Title, paste(rep("\n ", Title_Drop_Spaces), collapse = ""))

  message("Formatting X Axis Title spacing above")

  X_Axis_Title_Drop_Spaces <- as.integer(X_Axis_Title_Drop_Spaces)

  if (is.na(X_Axis_Title_Drop_Spaces) || X_Axis_Title_Drop_Spaces < 0) {

    X_Axis_Title_Drop_Spaces <- 1

  }

  X_Axis_Title <- paste0(
    paste(rep("\n ", X_Axis_Title_Drop_Spaces), collapse = ""),
    X_Axis_Title
  )

  message("Formatting Y Axis Title spacing beneath")

  Y_Axis_Title_Drop_Spaces <- as.integer(Y_Axis_Title_Drop_Spaces)

  if (is.na(Y_Axis_Title_Drop_Spaces) || Y_Axis_Title_Drop_Spaces < 0) {

    Y_Axis_Title_Drop_Spaces <- 0

  }

  Y_Axis_Title <- paste0(

    Y_Axis_Title,
    paste(rep("\n ", Y_Axis_Title_Drop_Spaces), collapse = "")

  )

  if(Interactive == TRUE)

  {

    message("Setting Interactive hover INFO")

    Data$Hover_Info <- paste0(
      "SNP: ", as.character(Data[[SNP_ID_Column]]), "\n",
      "CHR: ", as.character(Data[[Chromosome_Column]]), "\n",
      "POS: ", as.character(Data[[Position_Column]]), "\n",
      "P: ",  signif(as.numeric(Data[[PValue_Column]]), 2), "\n",
      "REF: ", as.character(Data[[Ref_Allele_Column]]), "\n",
      "ALT: ", as.character(Data[[Alt_Allele_Column]])

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

  if (length(unique(Data[[Chromosome_Column]])) == 1)

  {

  message("Adapting for Regional_Plot()")

  suppressMessages(suppressWarnings({

  chroms <- unique(Data[[Chromosome_Column]])

  a <- ggmanh::manhattan_plot(x = Data, preserve.position = T, plot.title = Title,
                      chr.colname = Chromosome_Column, pos.colname = Position_Column, label.colname = NULL,
                      pval.colname = PValue_Column, annotateTop = FALSE, chr.order = c(chroms), chr.col = Chromosome_Colours,
                      chrlabs = c(1:22, "X", "Y", "M"), rescale = Condense_Scale, signif = Break_Point, rescale.ratio.threshold = 0.0,
                      signif.rel.pos = 0.8, signif.col = c("transparent"),  point.size = Point_Size, x.label = "", y.label = "")


  # Use detected chromosome column instead of hard-coded CHROM

  chr_field <- Chromosome_Column

  # Define alternating chromosome colours:

  unique_chr <- sort(unique(a$data[[chr_field]]))

  chr_colour_map <- setNames(

    rep(Chromosome_Colours, length.out = length(unique_chr)),
    unique_chr

  )

  # Assign alternating colours:

  a$data$colour_group <- chr_colour_map[as.character(a$data[[chr_field]])]

  #   Was testing here for included HLA plot - more amenable to split for now as band plot more useful. Keep code here for future.

  #   if ("CLASS" %in% names(a$data)) {
  #
  #   # Fixed colours per CLASS
  #
  #   cls_colour_map <- c(
  #     AA   = "#1f77b4",  # blue
  #     HLA  = "#d62728",  # red
  #     RSID = "#2ca02c",  # green
  #     SNPS = "#9467bd",  # purple
  #     OTHER = "grey60"   # fallback
  #   )
  #
  #   # Ensure CLASS is a factor with the defined levels; map NAs to OTHER
  #
  #   a$data$CLASS <- as.character(a$data$CLASS)
  #   a$data$CLASS[is.na(a$data$CLASS)] <- "OTHER"
  #   a$data$CLASS <- factor(a$data$CLASS, levels = names(cls_colour_map))
  #
  #
  #
  #   a$data$colour_group <- unname(cls_colour_map[ as.character(a$data$CLASS) ])
  #
  # }

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
                              chr.order = c(chroms), # for single
                              chr.col = Chromosome_Colours, chrlabs = c(1:22, "X", "Y", "M"), rescale = Condense_Scale,
                              signif = Break_Point, rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
                              signif.col = c("transparent"),  point.size = Point_Size, x.label = "", y.label = "")

  }))

  }

  }else{

  normalize_chrom <- function(x) {

    y <- as.character(x)
    y <- toupper(y)
    y <- gsub("^CHR", "", y)  # drop "chr"/"CHR" prefix

    y[y == "MT" | y == "MITO" | y == "MTDNA"] <- "M"
    y[y == "23"] <- "X"
    y[y == "24"] <- "Y"
    y[y == "25"] <- "M"
    y

  }

  message("Determining available chromosomes")

  desired_order <- c(as.character(1:22), "X", "Y", "M")

  x <- as.character(Data[[Chromosome_Column]])   # already normalized

  idx <- match(x, desired_order, nomatch = 0L)     # map each row to desired_order

  tight_order <- desired_order[ tabulate(idx, nbins = length(desired_order)) > 0L ]

  message("Using standard for Single_Plot()")

  # Pretty count like 123,456

  .point_count <- formatC(nrow(Data), format = "d", big.mark = ",")


  a <- .progress_break(

    expr = suppressMessages(suppressWarnings(

      ggmanh::manhattan_plot(
        x = Data, preserve.position = TRUE, plot.title = Title,
        chr.colname = Chromosome_Column, pos.colname = Position_Column, label.colname = NULL,
        pval.colname = PValue_Column, annotateTop = FALSE,
        chr.order = tight_order, chr.col = Chromosome_Colours,
        chrlabs = c(1:22, "X", "Y", "M"),
        rescale = Condense_Scale, signif = Break_Point,
        rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
        signif.col = c("transparent"), point.size = Point_Size,
        x.label = "", y.label = ""

      )

    )),

    msg = sprintf("Rendering %s points\u2026", .point_count)

  )

  if (!is.null(Top_Expand)) {

  a_Orig <- .progress_break(
    expr = suppressMessages(suppressWarnings(
      ggmanh::manhattan_plot(
        x = Data_Orig, preserve.position = TRUE, plot.title = Title,
        chr.colname = Chromosome_Column, pos.colname = Position_Column, label.colname = NULL,
        pval.colname = PValue_Column, annotateTop = FALSE,
        chr.order = tight_order, chr.col = Chromosome_Colours,
        chrlabs = c(1:22, "X", "Y", "M"),
        rescale = Condense_Scale, signif = Break_Point,
        rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
        signif.col = c("transparent"), point.size = Point_Size,
        x.label = "", y.label = ""
      )
    )),
    msg = sprintf("Rendering %s points\u2026", .point_count)

  )


  }

  message("Define alternating chromosome colours")

  unique_chr <- sort(unique(a$data[[Chromosome_Column]]))

  chr_colour_map <- setNames(
    rep(Chromosome_Colours, length.out = length(unique_chr)),
    unique_chr
  )

  message("Verifying colour scale")

  # Use detected chromosome column in the ggmanh data

  chr_field <- Chromosome_Column

  # Assign alternating colours:

  a$data$colour_group <- chr_colour_map[as.character(a$data[[chr_field]])]

  message("Swapping lead SNPs for diamond index point shape")

  suppressMessages(suppressWarnings({

    # Override for top==TRUE on selected chromosomes

    a$data$colour_group[a$data$top == TRUE & a$data[[chr_field]] %in% Chromosome_Diamond] <- "transparent"

  }))

  if (!is.null(Top_Expand)) {

  # hide the forced top anchor point by its ID

  id_field <- SNP_ID_Column

  if (!is.null(id_field) && id_field %in% names(a$data)) {

    a$data$colour_group[as.character(a$data[[id_field]]) == "FAKETOP70"] <- "transparent"

   }

  }

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
                              chr.order = tight_order,
                              chr.col = Chromosome_Colours, chrlabs = c(1:22, "X", "Y", "M"), rescale = Condense_Scale,
                              signif = Break_Point, rescale.ratio.threshold = 0.0, signif.rel.pos = 0.8,
                              signif.col = c("transparent"),  point.size = Point_Size, x.label = "", y.label = "")

      }))

    }

  }

 if(is.null(Top_Expand))

 {

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


 }else{

   b <- a

   if(Diamond_Index == TRUE)

   {

   b2 <- a2

   }

 }

  message("Formatting Plot")

  # tick length: 0 when called by Regional_Plot(), else 0.8 cm

  tick_len_x <- if (is_called_by_regional_plot) grid::unit(0, "pt") else grid::unit(0.1, "cm")

  p <- b + ggplot2::theme(panel.border =  ggplot2::element_blank(),
                 panel.grid.major =  ggplot2::element_blank(),
                 panel.grid.minor =  ggplot2::element_blank(),
                 panel.background =  ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(hjust = 0.5, size = Title_Size),
                 axis.text.x = ggplot2::element_text(size = Chromosome_Label_Size, vjust = X_Axis_Text_Vjust, colour = "black"),
                 axis.line =  ggplot2::element_line(),
                 axis.title.x =  ggplot2::element_text(size = X_Axis_Title_Size, vjust = X_Axis_Title_Vjust ),
                 axis.title.y =  ggplot2::element_text(size = Y_Axis_Title_Size, vjust = Y_Axis_Title_Vjust),
                 axis.text.y = ggtext::element_markdown(size = Y_Axis_Text_Size, vjust = Y_Axis_Text_Vjust, colour = "black"),

                 axis.ticks.length.x  =  tick_len_x,
                 legend.position = "none",
                 plot.margin =  ggplot2::margin(t = 20,  # Top margin
                                      r = 40,  # Right margin
                                      b = 60,  # Bottom margin
                                      l = 60)) # Left margin


  if(Diamond_Index == TRUE)

  {

  p2 <- b2 + ggplot2::theme(panel.border =  ggplot2::element_blank(),
                          panel.grid.major =  ggplot2::element_blank(),
                          panel.grid.minor =  ggplot2::element_blank(),
                          panel.background =  ggplot2::element_blank(),
                          plot.title = ggplot2::element_text(hjust = 0.5, size = Title_Size),
                          axis.text.x = ggplot2::element_text(size = Chromosome_Label_Size, vjust = X_Axis_Text_Vjust, colour = "black"),
                          axis.line =  ggplot2::element_line(),
                          axis.title.x =  ggplot2::element_text(size = X_Axis_Title_Size, vjust = X_Axis_Title_Vjust ),
                          axis.title.y =  ggplot2::element_text(size = Y_Axis_Title_Size, vjust = Y_Axis_Title_Vjust),
                          axis.text.y = ggtext::element_markdown(size = Y_Axis_Text_Size, vjust = Y_Axis_Text_Vjust, colour = "black"),

                          axis.ticks.length.x  =  ggplot2::unit(0.8,"cm"),
                          legend.position = "none",
                          plot.margin =  ggplot2::margin(t = 20,  # Top margin
                                                         r = 40,  # Right margin
                                                         b = 60,  # Bottom margin
                                                         l = 60)) # Left margin


  }

  message("Formatting axes and labels")

  p <- p +  ggplot2::xlab(X_Axis_Title) + ggplot2::theme(axis.title.x =  ggplot2::element_text(vjust= X_Axis_Title_Vjust, size = X_Axis_Title_Size, hjust = 0.5))
  p <- p +  ggplot2::ylab(Y_Axis_Title) + ggplot2::theme(axis.title.y =  ggplot2::element_text(vjust = Y_Axis_Title_Vjust, size = Y_Axis_Title_Size))


  if(Diamond_Index == TRUE)

  {

  p2 <- p2 +  ggplot2::xlab(X_Axis_Title) + ggplot2::theme(axis.title.x =  ggplot2::element_text(vjust= X_Axis_Title_Vjust, size = X_Axis_Title_Size, hjust = 0.5))
  p2 <- p2 +  ggplot2::ylab(Y_Axis_Title) + ggplot2::theme(axis.title.y =  ggplot2::element_text(vjust= Y_Axis_Title_Vjust, size = Y_Axis_Title_Size))


  }

  if(Interactive == TRUE)

  {

  message("Adjusting axes and labels for interactive mode")

  dot_padding <- "\n \n "  # add more lines if needed

  Y_Axis_Title <- paste0(Y_Axis_Title, dot_padding)
  X_Axis_Title <- paste0(dot_padding, X_Axis_Title)

  p <- p +  ggplot2::labs(y = Y_Axis_Title, x = X_Axis_Title, title = Title)

  if(Diamond_Index == TRUE)

    {

    p2 <- p2 +  ggplot2::labs(y = Y_Axis_Title, x = X_Axis_Title)

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

  message("Anchoring")

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

  if(Colour_Index == TRUE & (Diamond_Index == FALSE | !is.null(Diamond_Index_SNP) ))

  {

    message("Colouring Index circle")

    df <- as.data.frame(p$data)

    special_points <- df %>% dplyr::filter(COLOUR == 4)

    # Also need if only one CHROM passed to Manhattan

    if (is_called_by_regional_plot == TRUE | length(unique(Data[[Chromosome_Column]])) == 1)

    {

      p <- p + ggplot2::geom_point(
        data = special_points,
        ggplot2::aes(x = .data[[Position_Column]], y = log10pval),
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

      chr_field <- Chromosome_Column

      chroms_all_fake <- df %>%
        dplyr::group_by(.data[[chr_field]]) %>%
        dplyr::filter(all(P == 1 & FAKE == 1)) %>%
        dplyr::pull(.data[[chr_field]]) %>%
        unique()

     # Filter out such chromosomes from special_points

      special_points <- special_points[!(special_points[[Chromosome_Column]] %in% chroms_all_fake), ]

 }

    message("Filtering Diamond Index SNPs")

    special_points <- special_points[special_points[[Chromosome_Column]] %in% Chromosome_Diamond, ]

    if(is_called_by_regional_plot == TRUE)

    {

      suppressMessages({

        p <- p + ggplot2::geom_point(
          data = special_points,
          ggplot2::aes(x = .data[[Position_Column]], y = log10pval),
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

  suppressMessages(suppressWarnings({

  p <- p +  ggplot2::coord_cartesian(clip = "off")

  }))

  message("Scaling Axes")

  if (length(unique(Data[[Chromosome_Column]])) > 1) {

  message("Using full scaling")

  df <- as.data.frame(p$data)

  middle_new_pos_values <- numeric()

  for (chrom in c(1:22, "X", "Y", "M")) {

    # Filter for the current chromosome (use detected column)

    df_filtered <- df %>%
      dplyr::filter(.data[[Chromosome_Column]] == chrom)

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

  chroms_in_data <- unique(as.character(Data[[Chromosome_Column]]))

  # Define maximum full label set

  all_chromosomes <- as.character(c(1:22, "X", "Y", "M"))

  # Filter out dropped and missing chromosomes

  chromosome_labels <- setdiff(all_chromosomes, as.character(Chromosome_Label_Drops))
  chromosome_labels <- intersect(chromosome_labels, chroms_in_data)

  # Ensure label order matches middle_new_pos_values

  final_labels <- ifelse(all_chromosomes %in% chromosome_labels, all_chromosomes, "")

  # defaults

  right_pad <- 0.005
  left_pad  <- 0.005

  # override if provided/not NULL

  if (!is.null(Right_Mult) && is.finite(Right_Mult)) {
    right_pad <- Right_Mult
  }

  if (!is.null(Left_Mult) && is.finite(Left_Mult)) {
    left_pad <- Left_Mult
  }


  if(!is.null(Top_Expand))  {

    suppressMessages(suppressWarnings({

    d <- p +  ggplot2::scale_x_continuous(breaks = middle_new_pos_values,
                                          labels = final_labels,
                                          position = "bottom",
                                          expand =  ggplot2::expansion(mult = c(left_pad, right_pad)) )

    }))

    # match axis linewidth from your theme (axis.line)

    axis_lwd <- 0.5

    tmp <- try(ggplot2::calc_element("axis.line", d$theme), silent = TRUE)

    if (!inherits(tmp, "try-error")) {

      if (!is.null(tmp$linewidth)) axis_lwd <- tmp$linewidth

      else if (!is.null(tmp$size)) axis_lwd <- tmp$size

    }

    # remove default x-axis line

    d <- d + ggplot2::theme(axis.line.x = ggplot2::element_blank())

    # compute axis line end INCLUDING the right-pad expansion multiplier

    df_axis <- df

    if ("FAKE" %in% names(df_axis)) {

      df_axis <- df_axis[df_axis$FAKE != 1, , drop = FALSE]

    }

    x_min <- min(df_axis$new_pos, na.rm = TRUE)
    x_max <- max(df_axis$new_pos, na.rm = TRUE)

    # expansion(mult=...) adds (range * mult) padding on each side

    x_span <- x_max - x_min

    if (!is.finite(x_span) || x_span == 0) x_span <- 1

    x_axis_end <- x_max + (x_span * right_pad)

    # draw x-axis line up to the expanded axis end (respects right_pad)

    d <- d + ggplot2::annotate(
      "segment",
      x = -Inf, xend = x_axis_end,
      y = -Inf, yend = -Inf,
      colour = "black",
      linewidth = axis_lwd,
      lineend = "butt"
    )

    d <- d + ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 60, l = 60))

    sig_y <- -log10(Sig_Threshold)

    if (!is.null(Top_Expand)) {

      # Use the computed details

      x_start <- x_min - (x_span * left_pad)
      x_end   <- x_max + (x_span * right_pad)

      d <- d + ggplot2::annotate(
        "segment",
        x = x_start, xend = x_end,
        y = sig_y,   yend = sig_y,
        linetype = Sig_Line_Type,
        colour   = Sig_Line_Colour,
        linewidth = Sig_Line_Width
      )

    }

  }else{

    suppressMessages(suppressWarnings({

    d <- p +  ggplot2::scale_x_continuous(breaks = middle_new_pos_values,
                                          labels = final_labels,
                                          position = "bottom",
                                          expand =  ggplot2::expansion(mult = c(left_pad, right_pad)) )

    }))

  }

    if(!is.null(Top_Expand))

     {

    # Extract CURRENT y-axis breaks from ggplot object `d`

    # get the computed y labels from ggplot object `d` (as shown: "0 | 2 | ...")

    y_idx <- which(vapply(d$scales$scales, function(s) "y" %in% s$aesthetics, logical(1)))[1]

    sc_y  <- d$scales$scales[[y_idx]]

    bld <- ggplot2::ggplot_build(d)
    pp  <- bld$layout$panel_params[[1]]

    y_breaks <- if (!is.null(pp$y$breaks)) pp$y$breaks else if (!is.null(pp$y.major)) pp$y.major else sc_y$get_breaks()

    y_labels <- sc_y$get_labels(y_breaks)

    cat(paste(y_labels, collapse = " | "), "\n")

    expanded_labels <- y_labels

    # attach the vector to the ggplot object (as an attribute)

    attr(d, "expanded_labels") <- expanded_labels

    y_idx <- which(vapply(a_Orig$scales$scales, function(s) "y" %in% s$aesthetics, logical(1)))[1]

    sc_y  <- a_Orig$scales$scales[[y_idx]]

    bld <- ggplot2::ggplot_build(a_Orig)
    pp  <- bld$layout$panel_params[[1]]

    y_breaks <- if (!is.null(pp$y$breaks)) pp$y$breaks else if (!is.null(pp$y.major)) pp$y.major else sc_y$get_breaks()

    y_labels <- sc_y$get_labels(y_breaks)

    cat(paste(y_labels, collapse = " | "), "\n")

    orig_labels <- y_labels

    # numeric max of orig_labels (robust if labels are characters/markdown)

    orig_labels_num <- suppressWarnings(as.numeric(gsub("[^0-9\\.\\-]+", "", as.character(orig_labels))))
    orig_labels_num <- orig_labels_num[is.finite(orig_labels_num)]

    max_orig_label <- if (length(orig_labels_num)) max(orig_labels_num, na.rm = TRUE) else NA_real_

     # remove default y-axis line (it runs the full height)

     d <- d + ggplot2::theme(axis.line.y = ggplot2::element_blank())

     # draw y-axis line ONLY up to y = 50

     y_cut <- max_orig_label + 1

     # start at baseline (0) or data min if it dips below 0

     y0 <- 0

     if ("log10pval" %in% names(df)) {

       y0 <- min(0, min(df$log10pval, na.rm = TRUE))

     }

     # match axis linewidth from the theme (axis.line)

     axis_lwd <- 0.5

     tmp <- try(ggplot2::calc_element("axis.line", d$theme), silent = TRUE)

     if (!inherits(tmp, "try-error")) {

       if (!is.null(tmp$linewidth)) axis_lwd <- tmp$linewidth

       else if (!is.null(tmp$size)) axis_lwd <- tmp$size

     }

     d <- d + ggplot2::annotate(
       "segment",
       x = -Inf, xend = -Inf,
       y = -Inf, yend = y_cut,   # start at -Inf so it meets the x-axis line
       colour = "black",
       linewidth = axis_lwd,
       lineend = "butt"
     )

     # shrink y tick labels above 50, without changing the scale

     y_threshold <- max_orig_label + 1

     small_y <- max(0, round(Y_Axis_Text_Size * 0.0))  # very small above

     y_idx <- which(vapply(d$scales$scales, function(s) "y" %in% s$aesthetics, logical(1)))[1]

     if (!is.na(y_idx) && length(y_idx)) {

       d$scales$scales[[y_idx]]$labels <- function(x) {

         vapply(x, function(v) {

           if (is.na(v)) return(NA_character_)

           txt <- format(v, trim = TRUE, scientific = FALSE)

           if (v > y_threshold) {

             sprintf("<span style='font-size:%dpt'>%s</span>", small_y, txt)

           } else {

             txt

           }

         }, character(1))

       }

     }

     # Force to max original

     y_threshold <- max_orig_label + 1

     # wanted_breaks = the literal numeric values present in orig_labels (no seq/pretty)

     wanted_breaks <- suppressWarnings(as.numeric(gsub("[^0-9\\.\\-]+", "", as.character(orig_labels))))
     wanted_breaks <- wanted_breaks[is.finite(wanted_breaks)]

     # (optional) keep only within 0..y_threshold, and ensure clean ordering/uniqueness

     wanted_breaks <- sort(unique(wanted_breaks[wanted_breaks >= 0 & wanted_breaks <= y_threshold]))

     # Find the y scale object that already exists on the plot

     y_idx <- which(vapply(d$scales$scales, function(s) "y" %in% s$aesthetics, logical(1)))[1]

     if (!is.na(y_idx) && length(y_idx)) {

       d$scales$scales[[y_idx]]$breaks <- wanted_breaks

       d$scales$scales[[y_idx]]$labels <- function(x) format(x, trim = TRUE, scientific = FALSE)

     }

  }

  }else{

    message("Determining X axis optimal breaks and labels")

    pos_range <- range(Data[[Position_Column]], na.rm = TRUE)
    breaks <- pretty(pos_range, n = 8)

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

    # defaults

    right_pad <- 0.01
    left_pad  <- 0.01

    # override if provided

    if (!is.null(Right_Mult) && is.finite(Right_Mult)) {

      right_pad <- Right_Mult

    }

    if (!is.null(Left_Mult) && is.finite(Left_Mult)) {

      left_pad <- Left_Mult

    }

    suppressMessages(suppressWarnings({

    d <- p + ggplot2::scale_x_continuous(
      breaks = breaks,
      labels = labels,
      position = "bottom" ,
      expand =  ggplot2::expansion(mult = c(left_pad, right_pad))
    )

    }))

  }

  # Due to Verbose = FALSE layering issue.

  Plot_Outcome <- d

  if (is_called_by_regional_plot) { # need separate defaults here due to scaling

    y_mult <- c(0.03, 0.05)  # c(lower, upper)

    # override lower if provided

    if (!is.null(Lower_Mult) && is.finite(Lower_Mult)) {

      y_mult[1] <- Lower_Mult

    }

    # override upper if provided

    if (!is.null(Upper_Mult) && is.finite(Upper_Mult)) {

      y_mult[2] <- Upper_Mult

    }

    idx <- which(sapply(Plot_Outcome$scales$scales, function(s) "y" %in% s$aesthetics))

    Plot_Outcome$scales$scales[[idx]]$expand <- ggplot2::expansion(mult = y_mult)

    x_mult <- c(0.015, 0.015)  # c(left, right) - default regional

    if (!is.null(Left_Mult) && is.finite(Left_Mult)) {

      x_mult[1] <- Left_Mult

    }

    if (!is.null(Right_Mult) && is.finite(Right_Mult)) {

      x_mult[2] <- Right_Mult

    }

    # x stays unchanged (default)

    idx <- which(sapply(Plot_Outcome$scales$scales, function(s) "x" %in% s$aesthetics))

    Plot_Outcome$scales$scales[[idx]]$expand <- ggplot2::expansion(mult = x_mult)

  }

   if (!is_called_by_regional_plot) {

    y_mult <- c(0.01, 0.01)  # c(bottom, top)

    # override if provided

    if (!is.null(Lower_Mult) && is.finite(Lower_Mult)) {

      y_mult[1] <- Lower_Mult

    }

    if (!is.null(Upper_Mult) && is.finite(Upper_Mult)) {

      y_mult[2] <- Upper_Mult

    }

  idx <- which(vapply(d$scales$scales, function(s) "y" %in% s$aesthetics, logical(1)))[1]

  if (!is.na(idx)) {

    d$scales$scales[[idx]]$expand <- ggplot2::expansion(mult = y_mult)

   }

  }

  if(Title_On == FALSE) {

    Plot_Outcome <- Plot_Outcome + ggplot2::labs(title = NULL)

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

    # resolve character Data as either file paths OR object names

    .resolve_char_data <- function(x, env) {
      x <- unname(x)

      # If they are all existing files

      if (length(x) >= 1 && all(file.exists(x))) {

        # If just ONE file: return a data.frame (not a list)

        if (length(x) == 1) {

          suppressMessages(

            suppressWarnings(

          df <- vroom::vroom(x[[1]], show_col_types = FALSE, progress = F)

            ))

          return(df)

        }

        # Multiple files: return named list

        suppressMessages(

          suppressWarnings(

        data_list <- lapply(x, function(f) vroom::vroom(f, show_col_types = FALSE, progress = F))

          ))

        names(data_list) <- basename(tools::file_path_sans_ext(x))

        return(data_list)

      }

      # Otherwise treat as object names

      missing <- x[!vapply(x, exists, logical(1), envir = env, inherits = TRUE)]

      if (length(missing)) {

        stop(

          "Data contained character values that are neither existing files nor objects in the calling environment: ",
          paste(missing, collapse = ", "),

          call. = FALSE

        )

      }

      objs <- lapply(x, function(nm) get(nm, envir = env, inherits = TRUE))

      if (!all(vapply(objs, is.data.frame, logical(1)))) {

        bad <- x[!vapply(objs, is.data.frame, logical(1))]

        stop(

          "These Data object names did not resolve to data frames: ",

          paste(bad, collapse = ", "),

          call. = FALSE

        )

      }

      # If just ONE object: return the data.frame directly (not a list)

      if (length(objs) == 1) {

        return(objs[[1]])

      }

      # Multiple objects: return named list

      names(objs) <- x

      objs

    }
    raw_data_expr <- call_expr[["Data"]]

    # special handling

    if (!is.null(raw_data_expr) && is.call(raw_data_expr) && identical(raw_data_expr[[1]], as.name("c"))) {

      arg_exprs <- as.list(raw_data_expr[-1])

      arg_vals  <- lapply(arg_exprs, function(e) eval(e, parent.frame()))

      if (all(vapply(arg_vals, is.data.frame, logical(1)))) {

        names(arg_vals) <- vapply(arg_exprs, function(e) deparse(e), character(1))
        clean_args$Data <- arg_vals

      } else if (all(vapply(arg_vals, is.character, logical(1)))) {

        # if not files, resolve as object names

        clean_args$Data <- .resolve_char_data(unlist(arg_vals), parent.frame())

      } else {

        stop("Mixed or invalid types in c(...) for Data.", call. = FALSE)

      }

    }

    Data <- clean_args$Data

    is_df <- function(x) is.data.frame(x)

    # handle Data = "Object" OR Data = c("Obj1","Obj2") even when NOT written as c(...) in the call

    if (is.character(Data) && length(Data) >= 1) {

      # Existing behavior for multiple file paths (still works), but now also supports object names

      clean_args$Data <- .resolve_char_data(Data, parent.frame())

      Data <- clean_args$Data

    }

    # existing: list to multi-plot

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

    # single df path already handled above; single df object below:

    if (is.data.frame(Data)) {

      expr <- raw_data_expr
      name <- "Data_1"

      if (!is.null(expr)) {

        if (is.character(expr)) {

          # Data = "Intelligence_Sum_Stats" - use the string without adding extra quotes to final plot title

          name <- expr

        } else {

          # Data = Intelligence_Sum_Stats - keep symbol name

          name <- deparse(expr)

        }

      }

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
