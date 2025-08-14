
#' Annotated, Summary Statistics with RSIDs
#'
#' @param Data These are the summary statistics (either file path(s),environmental object(s) or a variable(s) containing such contents) to be plotted; defaults to NULL
#' @param Chromosome_Column Manually specify chromosome column; defaults to NULL, leading to auto-detection
#' @param Position_Column Manually specify chromosome genomic position column; defaults to NULL, leading to auto-detection
#' @param SNP_ID_Column Manually specify SNP ID column; defaults to NULL, leading to auto-detection
#' @param PValue_Column Manually specify P Value column; defaults to NULL, leading to auto-detection
#' @param Reference_Allele_Column Manually specify reference allele column; defaults to NULL, leading to auto-detection
#' @param Effect_Allele_Column Manually specify effect allele column; defaults to NULL, leading to auto-detection
#' @param Verbose Prevent display of progress bar as function is running and instead show key milestone outputs/messages (mainly for debugging purposes); defaults to FALSE.
#' @param Genome_Build Reference genome to base required RSID annotations around; defaults to grch38
#' @param ... Shadow argument to detect inherited arguments passed and modified from the Single_Plot() function; defaults to NULL
#'
#' @return Annotated dataframe with Lab column containing RSIDs is allocated to specified object and the resulting summary statistics can then be passed to other MiamiR functions
#'
#' @export
#'
#' @examples  Annotated_Data <- Annotate_Data(Data = Intelligence_Sum_Stats)


Annotate_Data <- function(Data = NULL,
                          Chromosome_Column = NULL,
                          Position_Column = NULL,
                          SNP_ID_Column = NULL,
                          PValue_Column = NULL,
                          Reference_Allele_Column = NULL,
                          Effect_Allele_Column = NULL,
                          Genome_Build = "grch38",
                          ...,
                          Verbose = FALSE)

{

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

  message("Allocating Index SNPs in specified regions")

  Data <- Data %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(min_P_GENPOS = GENPOS[which.min(P)]) %>%
    dplyr::ungroup()

  message("Clearing current Lab information")

  Data$Lab <- NULL

  message("Creating space for new Lab column")

  Data <- Data %>%
    dplyr::mutate(Lab = ifelse(GENPOS == min_P_GENPOS, ID, ""))

  message("Just taking index SNPs for annotation")

  SNPs <- Data[Data$Lab != "",]

  message("Preparing query dataset")

  Data$Lab <- NULL

  #Get biomaRt search format
  SNPs <- SNPs %>% dplyr::select(CHROM, GENPOS, GENPOS, ALLELE0, ALLELE1)

  SNPs_Back <- SNPs %>%
    dplyr::mutate(
      ref_alt = paste(ALLELE0, ALLELE1, sep = "/"),
      alt_ref = paste(ALLELE1, ALLELE0, sep = "/")
    )

  SNPs$CHR <- SNPs$CHROM
  SNPs$chr_start <- SNPs$GENPOS
  SNPs$chr_end <- SNPs$GENPOS

  #Don't need these in search file
  SNPs$CHROM <- NULL
  SNPs$GENPOS <- NULL
  SNPs$ALLELE0 <- NULL
  SNPs$ALLELE1 <- NULL

  message("Creating SNP mart object depending on genome build")

  if(Genome_Build == "grch37")
  {

  snp_mart <- biomaRt::useEnsembl(biomart="ENSEMBL_MART_SNP",
                                               host="https://grch37.ensembl.org",
                                              dataset="hsapiens_snp")
  }

  if(Genome_Build == "grch38") #stable in 19
  {

    snp_mart <- biomaRt::useEnsembl(biomart="ENSEMBL_MART_SNP",
                                    host="https://grch37.ensembl.org",
                                    dataset="hsapiens_snp")

  }

   message("Formatting search terms")

   position <- apply(SNPs, 1, paste, collapse = ":")

   coords <- position

   coords <- gsub("^23:", "X:", coords)

   message("Creating storage results")

   c <- data.frame()

   # Loop through each row of the coords object
   for (i in 1:length(coords)) {

     # Extract the current row (chromosomal region) from coords
     current_coord <- coords[i]

     message(paste0("Processing: ", current_coord))

     # Query BioMart with the current coordinate
     a <- biomaRt::getBM(
       attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'allele'),
       filters = c('chromosomal_region'),
       values = current_coord,
       mart = snp_mart
     )

     # Convert the result to a dataframe
     b2 <- as.data.frame(a)

     message("Match(s): ")
     message(b2$refsnp_id)

     # Combine the current result with the previously accumulated results
     c <- rbind(c, b2)

     progress <- message(paste0("Obtained Index SNP RS code for ", i, " out of ", length(coords), " Index SNPs"))


   }

   message("Formatting receiving data frame and results")

   Data$CHROM <- as.character(Data$CHROM)
   SNPs_Back$CHROM <- as.character(SNPs_Back$CHROM)

   c <- c %>%
     dplyr::mutate(chr_name = as.character(chr_name))

   message("Joining and ensuring where multiple correct RSID is assigned")

   # Join and apply flexible allele filtering
   joined <- SNPs_Back %>%
     dplyr::inner_join(
       c,
       by = c("CHROM" = "chr_name", "GENPOS" = "chrom_start")
     ) %>%
     dplyr::filter(GENPOS == chrom_end) %>%
     dplyr::rowwise() %>%
     dplyr::filter({
       # Split alleles into parts
       allele_parts <- base::unlist(base::strsplit(base::toupper(allele), "/"))
       ref_parts <- base::unlist(base::strsplit(base::toupper(ref_alt), "/"))
       alt_parts <- base::unlist(base::strsplit(base::toupper(alt_ref), "/"))

       # Check if either direction matches all alleles
       base::all(ref_parts %in% allele_parts) || base::all(alt_parts %in% allele_parts)
     }) %>%
     dplyr::ungroup()

   # Join back the matching refsnp_id values
   SNPs_Back <- SNPs_Back %>%
     dplyr::left_join(
       joined %>% dplyr::select(CHROM, GENPOS, refsnp_id),
       by = c("CHROM", "GENPOS")
     )

    #Join SNPs to summary stats
    Data <- dplyr::left_join(
     Data,
     SNPs_Back %>%
       dplyr::select(GENPOS, CHROM, Lab = refsnp_id),
     by = c("GENPOS", "CHROM")
    )

    message("The following were annotated:")

    message(Data[!is.na(Data$Lab), c("ID", "Lab")])

    message("Returning annotated summary stats")

    #Return the df itself
    return(Data)

}

.Annotate_Data_original <- Annotate_Data

Annotate_Data <- function(..., session = NULL) {

  # Try to capture the name
  mc <- match.call(expand.dots = FALSE)
  data_arg_name <- NULL

  if (!is.null(mc$...)) {

    dots_expr <- as.list(mc$...)
    if ("Data" %in% names(dots_expr)) {

      data_arg_name <- deparse(dots_expr$Data, backtick = TRUE)

    }

  }

  args <- list(...)
  args$session <- session
  args$.dots <- args

  # Fallbacks: if Data is a file path (character), use that
  if (is.null(data_arg_name)) {

    if ("Data" %in% names(args) && is.character(args$Data) && length(args$Data) == 1) {

      data_arg_name <- args$Data

    } else {

      data_arg_name <- "Data"

    }
  }

  message(sprintf("Processing dataset: %s", data_arg_name))

  verbose_mode <- if ("Verbose" %in% names(args)) isTRUE(args$Verbose) else FALSE

  if (verbose_mode) {

    # Direct call with messages
    return(do.call(.Annotate_Data_original, args))

  } else {

    # Quiet run, but the one line above still shows

    return(
      suppressMessages(
        suppressWarnings(
          run_with_counter(
            func    = .Annotate_Data_original,
            args    = args,
            session = session
          )
        )
      )
    )
  }
}
