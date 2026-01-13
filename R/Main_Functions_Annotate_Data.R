
#' Annotated, Summary Statistics with chromosomal index RsIDs
#'
#' @param Data These are the summary statistics (either file path(s),environmental object(s) or a variable(s) containing such contents) to be plotted; defaults to NULL (can be a list for multiple runs)
#' @param Chromosome_Column Manually specify chromosome column; defaults to NULL, leading to auto-detection
#' @param Position_Column Manually specify chromosome genomic position column; defaults to NULL, leading to auto-detection
#' @param SNP_ID_Column Manually specify SNP ID column; defaults to NULL, leading to auto-detection
#' @param PValue_Column Manually specify P Value column; defaults to NULL, leading to auto-detection
#' @param Reference_Allele_Column Manually specify reference allele column; defaults to NULL, leading to auto-detection
#' @param Effect_Allele_Column Manually specify effect allele column; defaults to NULL, leading to auto-detection
#' @param Verbose Prevent display of progress bar as function is running and instead show key milestone outputs/messages (mainly for debugging purposes); defaults to FALSE.
#' @param Genome_Build Reference genome to base required RSID annotations around; defaults to grch38 (grch37 also allowed)
#' @param ... Shadow argument to detect inherited arguments passed and modified from the Single_Plot() function; defaults to NULL
#'
#' @return Annotated data frame with Lab column containing RsIDs for chromosome index SNPs is allocated to specified object and the resulting summary statistics can then be passed to other MiamiR functions
#'
#' @export
#'
#' @examples  Annotated_Data <- Annotate_Data(Data = Intelligence_Sum_Stats,  Genome_Build = "grch37")


Annotate_Data <- function(Data = NULL,
                          Chromosome_Column = NULL,
                          Position_Column = NULL,
                          SNP_ID_Column = NULL,
                          PValue_Column = NULL,
                          Reference_Allele_Column = NULL,
                          Effect_Allele_Column = NULL,
                          Genome_Build = "grch38",
                          Verbose = FALSE,
                          ...)

 {


  # Need P and ID for Indexing

  if (is.character(Data) && length(Data) == 1)

  {

    file_path <- Data

    message("Dataset absent from environment")

    if (file.exists(file_path))

    {

      message("Reading data from file: ", file_path)

      message("Loading data using vroom...")

      Data <- vroom::vroom(file_path, show_col_types = FALSE, progress = F)

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

  # don't need as processing RsID only and may be fed later to plots which do this

  # if (!is.null(PValue_Column) && grepl("log", PValue_Column, ignore.case = TRUE)) {
  #
  #   message("Converting log value to plain P")
  #
  #   x <- Data[[PValue_Column]]
  #
  #   x <- if (is.numeric(x)) x else as.numeric(x)
  #
  #   Data[[PValue_Column]] <- exp(-x * 2.302585092994046)  # 2.30258... = log(10)
  #
  # }

  Data[[Ref_Allele_Column]] <- toupper(Data[[Ref_Allele_Column]])
  Data[[Alt_Allele_Column]] <- toupper(Data[[Alt_Allele_Column]])

  message("Allocating Index SNPs in specified regions")

  Data <- Data %>%
    dplyr::group_by(.data[[Chromosome_Column]]) %>%
    dplyr::mutate(min_P_GENPOS = .data[[Position_Column]][which.min( .data[[PValue_Column]] )]) %>%
    dplyr::ungroup()

  message("Clearing current Lab information")

  Data$Lab <- NULL

  message("Creating space for new Lab column")

  Data <- Data %>%
    dplyr::mutate(Lab = ifelse(.data[[Position_Column]] == min_P_GENPOS, "TOP", ""))  # Can't use ID as sometimes NA ID

  message("Just taking index SNPs for annotation")

  SNPs <- Data[Data$Lab == "TOP",]

  message("Preparing query dataset")

  Data$Lab <- NULL

  # Get biomaRt search format

  SNPs <- SNPs %>% dplyr::select(.data[[Chromosome_Column]], .data[[Position_Column]], .data[[Position_Column]], .data[[Ref_Allele_Column]], .data[[Alt_Allele_Column]])

  # Indel normalization (simple anchor trimming)- rare

  #   Convert anchored GWAS-style alleles like:
  #     insertion:  A -> AC     into   - -> C
  #     insertion:  A -> AAC    into   - -> AC
  #     deletion:   AT -> A     into   T -> -
  #     deletion:   GTC -> G    into   TC -> -


  message("Filtered to biomaRt data required")

  SNPs <- SNPs %>%

    dplyr::mutate(

      !!Ref_Allele_Column := toupper(as.character(.data[[Ref_Allele_Column]])),
      !!Alt_Allele_Column := toupper(as.character(.data[[Alt_Allele_Column]]))

    )%>%

    dplyr::rowwise() %>%
    dplyr::mutate(

      # Keep originals for logic inside this row

      .ref0 = as.character(.data[[Ref_Allele_Column]]),
      .alt0 = as.character(.data[[Alt_Allele_Column]]),

      # Compute normalized REF allele (ALLELE0)

      !!Ref_Allele_Column := {

        ref <- .ref0
        alt <- .alt0

        # Only normalize if it's an indel (different lengths) and both exist

        if (nchar(ref) != nchar(alt) && nchar(ref) > 0 && nchar(alt) > 0) {

          # Trim common leading bases (shared anchor)

          # Example insertion: ref="A", alt="AC"  => trim "A" => ref="", alt="C"

          # Example deletion:  ref="AT", alt="A"  => trim "A" => ref="T", alt=""

          while (nchar(ref) > 0 &&
                 nchar(alt) > 0 &&
                 substr(ref, 1, 1) == substr(alt, 1, 1)) {

            ref <- substr(ref, 2, nchar(ref))
            alt <- substr(alt, 2, nchar(alt))

          }

          # If REF becomes empty, that means nothing is present on REF side

          # Represent empty as "-" to match Ensembl style

          if (nchar(ref) == 0) "-" else ref

        } else {

          # Not an indel (or missing allele) - leave as-is

          ref

        }

      },

      # Compute normalized ALT allele (ALLELE1)

      !!Alt_Allele_Column  := {

        ref <- .ref0
        alt <- .alt0

        if (nchar(ref) != nchar(alt) && nchar(ref) > 0 && nchar(alt) > 0) {

          # Do the exact same trimming

          while (nchar(ref) > 0 &&
                 nchar(alt) > 0 &&
                 substr(ref, 1, 1) == substr(alt, 1, 1)) {

            ref <- substr(ref, 2, nchar(ref))
            alt <- substr(alt, 2, nchar(alt))

          }

          # If ALT becomes empty, that means nothing is present on ALT side

          # Represent empty as "-" (deletion).

          if (nchar(alt) == 0) "-" else alt

        } else {

          alt

        }

      }

    ) %>%

    dplyr::ungroup() %>%

    # Drop temporary columns used only for row wise computation

    dplyr::select(-.ref0, -.alt0)

   SNPs_Back <- SNPs %>%
    dplyr::mutate(
      ref_alt = paste(.data[[Ref_Allele_Column]] , .data[[Alt_Allele_Column]], sep = "/"),
      alt_ref = paste(.data[[Alt_Allele_Column]], .data[[Ref_Allele_Column]], sep = "/")
    )

  message("Creating SNP mart object depending on genome build:")

  if(Genome_Build == "grch37")

  {

  message("grch37")

  snp_mart <- biomaRt::useEnsembl(biomart="ENSEMBL_MART_SNP",
                                               host="https://grch37.ensembl.org",
                                               mirror  = "www",
                                               dataset="hsapiens_snp")
  }

  if(Genome_Build == "grch38")  #

  {

    message("grch38")

    snp_mart <- biomaRt::useEnsembl(biomart="ENSEMBL_MART_SNP",
                             #       host="https://grch38.ensembl.org", - default host is hg38
                                    mirror  = "www",
                                    dataset="hsapiens_snp")

  }

   message("Formatting search terms")

  eff_nchar <- function(x) {

    x <- as.character(x)
    x <- ifelse(is.na(x) | x == "-", "", x)
    nchar(x)

  }

  # Per-row buffer column (W) as indels start and end away from anchor at times.

  # Treat "-" as length 0 so -/C counts as 1 bp insertion, -/AC counts as 2, etc.

  len0 <- eff_nchar(SNPs[[Ref_Allele_Column]])
  len1 <- eff_nchar(SNPs[[Alt_Allele_Column]])

  diff_len <- abs(len1 - len0) # so no neg at one calc

  # Indel if effective lengths differ (includes "-" cases automatically)

  is_indel <- diff_len > 0

  # Create a buffer column on SNPs (W):

  #  - SNPs: 0

  SNPs <- SNPs %>%
    dplyr::mutate(

      W = ifelse(is_indel, pmax(0L, diff_len), 0L),
      chr_start = pmax(1L, as.integer(.data[[Position_Column]]) - as.integer(W)),
      chr_end   = as.integer(.data[[Position_Column]]) + as.integer(W),
      coords    = paste0(.data[[Chromosome_Column]], ":", chr_start, ":", chr_end)

    )

  # Your existing X fix

  SNPs$coords <- gsub("^23:", "X:", SNPs$coords)

  # coords vector for the loop (same as before)

  coords <- SNPs$coords

    message("Creating storage results frame")

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

     # Convert the result to a data frame

     b2 <- as.data.frame(a)

     # Combine the current result with the previously accumulated results

     c <- rbind(c, b2)

     progress <- message(paste0("Searching Index SNP RsID for ", i, " out of ", length(coords), " Index SNPs"))

   }

   message("Formatting receiving data frame and results")

   Data[[Chromosome_Column]] <- as.character(Data[[Chromosome_Column]])
   SNPs_Back[[Chromosome_Column]] <- as.character(SNPs_Back[[Chromosome_Column]])

   c <- c %>%
     dplyr::mutate(chr_name = as.character(chr_name))

   message("Joining and ensuring where multiple correct RSID is assigned")

   # Join and apply flexible allele filtering

   joined <- SNPs_Back %>%
     dplyr::inner_join(c, by = stats::setNames("chr_name", Chromosome_Column)) %>%
     dplyr::filter(.data[[Position_Column]] == chrom_start |
                   .data[[Position_Column]] == chrom_end)  %>%

     # then make sure either start or end - fixed query of one base for single snps, allows flex inv

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
       joined %>% dplyr::select(.data[[Chromosome_Column]], .data[[Position_Column]], refsnp_id),
       by = c(Chromosome_Column, Position_Column)
     )

   missing_lab_snps <- SNPs_Back %>%
     dplyr::filter(is.na(refsnp_id) | refsnp_id == "")

   message(sprintf("SNPs with absent RsID (likely structural variation in complex regions): %d", nrow(missing_lab_snps)))

    # wrangled above not joining back anyway

    # Join SNPs to summary stats

    Data <- dplyr::left_join(

     Data,
     SNPs_Back %>%
       dplyr::select(.data[[Position_Column]], .data[[Chromosome_Column]], Lab = refsnp_id),
     by = c(Position_Column, Chromosome_Column)

    )

    message("Final Clean")

    Data$min_P_GENPOS <- NULL # final clean step

    message("Returning annotated summary stats")

    # Return the df itself

    return(Data)

}

if (!exists("use_wrapper")) use_wrapper <- TRUE # for shiny!

if (use_wrapper == TRUE) {

  .Annotate_Data_original <- Annotate_Data

  Annotate_Data <- function(..., session = NULL) {

    args <- list(...)
    args$session <- session
    args$.dots <- args

    # Where to look up object names supplied as strings

    # (search parent frames first, then global env)

    resolve_data_name <- function(name) {

      if (!is.character(name) || length(name) != 1) return(name)

      # Try parent frames

      pf <- parent.frame()

      if (exists(name, envir = pf, inherits = TRUE)) {

        return(get(name, envir = pf, inherits = TRUE))

      }

      # Then global env

      if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {

        return(get(name, envir = .GlobalEnv, inherits = FALSE))

      }

      # Otherwise keep as string (could be a filepath)

      name

    }

    is_multi_data <- function(x) {

      if (is.null(x)) return(FALSE)

      # data.frame/tibble is ALWAYS a single dataset

      if (is.data.frame(x)) return(FALSE)

      # list of datasets (but not a data.frame)

      if (is.list(x)) return(length(x) > 1)

      # atomic vector of datasets (multiple object names or file paths)

      if (is.atomic(x) && length(x) > 1) return(TRUE)

      FALSE

    }

    as_data_list <- function(x) {

      if (is.data.frame(x)) return(list(x))
      if (is.list(x)) return(x)
      if (is.atomic(x) && length(x) > 1) return(as.list(x))
      list(x)

    }

    # per-run argument picker:

    pick_i <- function(val, i, n) {

      if (is.null(val)) return(NULL)

      if (is.list(val) && !is.data.frame(val)) {

        if (length(val) == 1) return(val[[1]])

        if (length(val) == n) return(val[[i]])

        return(val[[1]])

      }

      if (is.atomic(val) && length(val) > 1) {

        if (length(val) == n) return(val[i])

        return(val[1])

      }

      val

    }

    infer_names <- function(dlist) {

      n <- length(dlist)
      nm <- names(dlist)

      if (!is.null(nm) && any(nzchar(nm))) {

        nm[nm == ""] <- paste0("run_", which(nm == ""))

        return(nm)

      }

      # if list elements are character (names/paths), use basename for pretty labels

      flat <- unlist(dlist, recursive = FALSE, use.names = FALSE)

      if (length(flat) == n && all(vapply(flat, is.character, logical(1)))) {

        out <- basename(flat)
        out[out == "" | is.na(out)] <- paste0("run_", seq_len(n))
        return(out)

      }

      paste0("run_", seq_len(n))

    }

    allowed_patterns <- c(
      "^Processing:\\s",
      "^Searching Index SNP RsID\\s+for\\s+"
    )

    allow_msg <- function(msg) {

      any(vapply(allowed_patterns, grepl, logical(1), x = msg))

    }

    run_one <- function(run_args) {

      # Resolve Data if it's a character object name (single string)

      if ("Data" %in% names(run_args)) {

        run_args$Data <- resolve_data_name(run_args$Data)

      }

      verbose_mode <- if ("Verbose" %in% names(run_args)) isTRUE(run_args$Verbose) else FALSE

      if (verbose_mode) {

        # Full messages, no counter wrapper

        return(do.call(.Annotate_Data_original, run_args))

      }

      # Verbose = FALSE:

      withCallingHandlers(

        suppressWarnings(

          run_with_counter(

            func    = .Annotate_Data_original,
            args    = run_args,
            session = session

          )

        ),

        message = function(m) {

          msg <- conditionMessage(m)

          if (allow_msg(msg)) {

            message(msg)  # allow only these messages through

          }

          invokeRestart("muffleMessage")  # obfuscate everything else (and avoid duplicates)

        }

      )

    }

    # MULTI-RUN

    if ("Data" %in% names(args) && is_multi_data(args$Data)) {

      Data_list <- as_data_list(args$Data)
      n <- length(Data_list)

      # Names based on the original user-supplied entries (strings become names)

      out_names <- infer_names(Data_list)

      results <- vector("list", n)
      names(results) <- out_names

      base_args <- args
      base_args$Data <- NULL

      for (i in seq_len(n)) {

        run_args <- base_args
        run_args$Data <- Data_list[[i]]

        # Per-run vector/list support for ALL args (including ...)

        for (nm in names(run_args)) {

          if (nm %in% c("session", ".dots")) next

          run_args[[nm]] <- pick_i(run_args[[nm]], i, n)

        }

        message(sprintf("Processing dataset %d/%d: %s", i, n, out_names[i]))

        results[[i]] <- run_one(run_args)
      }

      return(results)

    }

    # Resolve Data if it's a character object name (single string)

    if ("Data" %in% names(args)) {

      args$Data <- resolve_data_name(args$Data)

    }

    return(run_one(args))

  }

}
