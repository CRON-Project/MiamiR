
#' Munged summary statistics in preparation for METASOFT-like meta-analysis
#'
#' @param Match_Allele_Study Data set used as reference point for effect allele direction; defaults to NULL - automatically assigned as first set
#' @param Output A selection of columns to isolate, in order, alternating by dataset; defaults to NULL, but automatically assigned c("BETA", "SE") upon function call - standard METASOFT ordered format; GENPOS, P, CHROM, ALLELE0, ALLELE1 (all options)
#' @param Data These are the summary statistics (either file path(s), environmental object(s) or a variable(s) containing such contents) to be plotted; defaults to NULL (must be one list only)
#' @param BETA_Column Manually specify BETA Test Statistic column; defaults to NULL, leading to auto-detection
#' @param Position_Column Manually specify chromosome genomic position column; defaults to NULL, leading to auto-detection
#' @param SNP_ID_Column Manually specify SNP ID column; defaults to NULL, leading to auto-detection
#' @param PValue_Column Manually specify P Value column; defaults to NULL, leading to auto-detection
#' @param Reference_Allele_Column Manually specify reference allele column; defaults to NULL, leading to auto-detection
#' @param Effect_Allele_Column Manually specify effect allele column; defaults to NULL, leading to auto-detection
#' @param SE_Column Manually specify SE of Test Statistic column; defaults to NULL, leading to auto-detection
#' @param Chromosome_Column Manually specify chromosome column; defaults to NULL, leading to auto-detection
#' @param Output_Reference_Identifier Include the supplied internal identifier of the Match_Allele_Study dataset in final output; defaults to TRUE
#' @param Verbose Prevent display of progress bar as function is running and instead show key milestone outputs/messages (mainly for debugging purposes); defaults to FALSE
#'
#' @return Munged dataframe amenable to METASOFT meta-analysis - only including SNPs present in backbone/reference file as other SNPs effects cannot be consistently aligned
#' @export
#'
#' @examples METASOFT <- METASOFT_File_Gen(Data = c("Intelligence_Sum_Stats", "Household_Income_Sum_Stats"))
#'

METASOFT_File_Gen <- function(Data = NULL,
                              BETA_Column = NULL,
                              Chromosome_Column = NULL,
                              Position_Column = NULL,
                              SNP_ID_Column = NULL,
                              PValue_Column = NULL,
                              Reference_Allele_Column = NULL,
                              Effect_Allele_Column = NULL,
                              SE_Column = NULL,
                              Match_Allele_Study = NULL,
                              Output = NULL,
                              Output_Reference_Identifier = TRUE,
                              Verbose = FALSE

)

{


  # Progress-safe messaging (show NADA unless Verbose=TRUE OR inside .progress_break)

  if (!exists(".progress", inherits = TRUE))        .progress        <- function(...) invisible(NULL)

  if (!exists(".progress_done", inherits = TRUE))   .progress_done   <- function(...) invisible(NULL)

  if (!exists(".progress_pause", inherits = TRUE))  .progress_pause  <- function(...) invisible(NULL)

  if (!exists(".progress_resume", inherits = TRUE)) .progress_resume <- function(...) invisible(NULL)

  .allow_msgs <- isTRUE(Verbose)

  # local message() for Verbose

  message <- function(..., domain = NULL, appendLF = TRUE) {

    if (isTRUE(get0(".allow_msgs", ifnotfound = FALSE, inherits = TRUE))) {

      base::message(..., domain = domain, appendLF = appendLF)

    }

    invisible(NULL)

  }

  # hider helper messaging

  .quiet_call <- function(expr) {

    if (isTRUE(.allow_msgs)) return(force(expr))

    withCallingHandlers(

      force(expr),
      message = function(m) invokeRestart("muffleMessage")

    )

  }

  # pause progress bar, temporarily allow messages, run expr, resume bar

  .progress_break <- function(expr, msg = NULL, pad = 0L) {

    .progress_pause()
    cat("\r\n", sep = "")  # end progress line cleanly

    old <- .allow_msgs
    .allow_msgs <<- TRUE

    on.exit({

      .allow_msgs <<- old
      .progress_resume()
    }, add = TRUE)

    if (!is.null(msg)) {

      cat(strrep("\n", pad), msg, "\n", sep = "")
      utils::flush.console()

    }

    force(expr)

  }

     if(is.null(Output))

     {

     Output <- c("BETA", "SE")

     }

      Match_Allele_Direction = TRUE # Always!

      Data_Sets <- Data

      if (!is.null((Data_Sets))) {

        message("Using Data Set Names as study names")

        Names <- (Data_Sets)  # Use names of the Data_Sets vector

      }

      if (is.null(Match_Allele_Study)) {

        Match_Allele_Study <- Names[1]

        message(paste0("Matching allele directions automatically to ", Match_Allele_Study) )

      }

      # force reference study first so we can filter others immediately

      Names <- c(Match_Allele_Study, setdiff(Names, Match_Allele_Study))

      # have to reorder this too

      Data_Sets <- Names

      # Initialise once before loop

      Chromosome_Columns   <- list()
      PValue_Columns       <- list()
      Position_Columns     <- list()
      SNP_ID_Columns       <- list()
      Ref_Allele_Columns   <- list()
      Alt_Allele_Columns   <- list()
      BETA_Columns         <- list()
      SE_Columns           <- list()

      Loaded_Data <- list()

      ref_keys <- NULL  # will hold GENPOS keys from reference study

      for (Set in Names) {

        Spec_Name <- Set

        message(Spec_Name)



        if ((is.character(Spec_Name) && grepl("/", Spec_Name)) ||

            (is.character(Spec_Name) && file.exists(Spec_Name))) {

          # Treat as file path

          message("Loading ", Spec_Name, " from file path")

          Data <- vroom::vroom(Spec_Name, show_col_types = FALSE, progress = FALSE) # don't want interference

        } else if (exists(Spec_Name, envir = .GlobalEnv)) {

          # Treat as object in the environment

          Data <- get(Spec_Name, envir = .GlobalEnv)

        } else {

          stop(sprintf("Spec_Name '%s' is not an object or a valid file path", Spec_Name))

        }

        Chromosome_Columns[[Spec_Name]] <- .quiet_call(detect_chromosome_column(Data, Chromosome_Column))

        message("Addressing variable chromosome X nomenclature")

        col <- Chromosome_Columns[[Spec_Name]]
        v <- Data[[col]]

        # detect whether this dataset needs a 23 -> X

        need_fix <- FALSE
        detected_as <- NULL

        need_fix <- any(v == 23L, na.rm = TRUE)
        detected_as <- if (need_fix) "numeric value 23" else "no value 23"

        # v is always numeric here

        has_23 <- any(v == 23L, na.rm = TRUE)

        if (has_23) {

          message("[CHROM] ", Spec_Name, ": detected numeric value 23 → converting 23 -> X")

         vc <- as.character(v)        # only pay computational cost if needed
         vc[v == 23L] <- "X"
         Data[[col]] <- vc

        } else {

          # numeric columns cannot contain "X"

          message("[CHROM] ", Spec_Name, ": can't find 23 or X → probably no X data")

        }

        message("Assigning other key data columns")

        PValue_Columns[[Spec_Name]]     <- .quiet_call(detect_pvalue_column(Data, PValue_Column))
        Position_Columns[[Spec_Name]]   <- .quiet_call(detect_position_column(Data, Position_Column))

        col_pos <- Position_Columns[[Spec_Name]]

        pos <- Data[[col_pos]]
        pos_i <- suppressWarnings(as.integer(pos))   # keep if needed

        if (Spec_Name == Match_Allele_Study) {

          ref_pos_set <- unique(pos_i)

          message("Reference unique GENPOS: ", length(ref_pos_set))

        } else {

          before <- nrow(Data)

          keep <- !is.na(pos_i) & (pos_i %in% ref_pos_set)
          Data <- Data[keep, , drop = FALSE]

          after <- nrow(Data)

          message("Kept ", after, " / ", before, " rows (dropped ", before - after, ") for lazy filtering")

        }

        SNP_ID_Columns[[Spec_Name]]     <- .quiet_call(detect_snp_column(Data, SNP_ID_Column))
        Ref_Allele_Columns[[Spec_Name]] <- .quiet_call(detect_reference_allele_column(Data, Reference_Allele_Column))
        Alt_Allele_Columns[[Spec_Name]] <- .quiet_call(detect_effect_allele_column(Data, Effect_Allele_Column))
        BETA_Columns[[Spec_Name]]       <- .quiet_call(detect_beta_column(Data, BETA_Column))
        SE_Columns[[Spec_Name]]         <- .quiet_call(detect_se_column(Data, SE_Column))

        Loaded_Data[[Spec_Name]] <- Data

      }

      Combined_Processed_Data <- data.frame()

      for (i in seq_along(Data_Sets)) {

      dataset_name <- Data_Sets[i]  # Get the dataset name

      Data <- Loaded_Data[[dataset_name]]

      corresponding_name <- Names[i]  # Get the corresponding name

      message(paste0("Processing: ", corresponding_name))

      Chromosome_Column <- Chromosome_Columns[[i]]
      Position_Column <- Position_Columns[[i]]
      SNP_ID_Column <- SNP_ID_Columns[[i]]
      PValue_Column <- PValue_Columns[[i]]
      Standard_Error_Column <- SE_Columns[[i]]
      Beta_Column <- BETA_Columns[[i]]
      Reference_Allele <- Ref_Allele_Columns[[i]]
      Effect_Allele <- Alt_Allele_Columns[[i]]

      # Print dataset name and corresponding name

      message(paste("Corresponding name:", corresponding_name))

      # Manually assign columns for ease of use - would normally get directly but easier here as cross-comparing different dfs.

      Data$CHROM <- Data[[Chromosome_Column]]
      Data$GENPOS <- Data[[Position_Column]]
      Data$ID <- Data[[SNP_ID_Column]]
      Data$P <- Data[[PValue_Column]]
      Data$SE <- Data[[Standard_Error_Column]]
      Data$BETA <- Data[[Beta_Column]]
      Data$STUDY <- corresponding_name

      if (!is.null(PValue_Column) && grepl("log", PValue_Column, ignore.case = TRUE)) {

        message("Converting log value to plain P")

        x <- Data[[PValue_Column]]

        x <- if (is.numeric(x)) x else as.numeric(x)

        Data[[PValue_Column]] <- exp(-x * 2.302585092994046)  # 2.30258... = log(10)

      }

      message("Ensure uppercase alleles")

      cols <- c(Reference_Allele, Effect_Allele)

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

     Data$ALLELE0 <- (Data[[Reference_Allele]])
     Data$ALLELE1 <- (Data[[Effect_Allele]])

     message("Subsetting key data")

     Data <- Data %>% dplyr::select(ID, ALLELE0, ALLELE1, CHROM, GENPOS, BETA, SE, P, STUDY)

     message("Setting a reference allele column")

     # Data$COORD_Norm <- stringi::stri_c("chr", Data$CHROM, ":", Data$GENPOS, ":", Data$ALLELE0, ":", Data$ALLELE1)

    n <- nrow(Data)

    ch  <- Data[["CHROM"]]
    pos <- Data[["GENPOS"]]
    a0  <- Data[["ALLELE0"]]
    a1  <- Data[["ALLELE1"]]

    out_norm <- character(n)

    bs <- 2e6
    n_chunks <- ceiling(n / bs)

    .progress_break({
      t0 <- Sys.time()

      for (k in seq_len(n_chunks)) {

        i <- (k - 1L) * bs + 1L
        j <- min(k * bs, n)

        out_norm[i:j] <- stringi::stri_paste(
          "chr", ch[i:j], ":", pos[i:j], ":", a0[i:j], ":", a1[i:j],
          sep = ""
        )

        # your chunk status lines (only visible while bar is paused)

        base::message(sprintf(

          "Creating COORD_Norm for %s | chunk %d/%d (rows %d–%d) | %.1fs",
          corresponding_name,
          k, n_chunks, i, j,
          as.numeric(difftime(Sys.time(), t0, units = "secs"))

        ))

      }

    }, msg = paste0("Building COORD_Norm (standard, reference-direction coordinate) in ", n_chunks, " chunks:"))

    Data$COORD_Norm <- out_norm

    message("Deriving COORD_Alt by allele swap")

    # fast swap of last two fields: :A:B -> :B:A

    Data$COORD_Alt <- sub(

      ":(.):(.)$",
      ":\\2:\\1",
      Data$COORD_Norm,
      perl = TRUE

    )

    message("Finished COORD build in ",
            round(difftime(Sys.time(), t0, units = "secs"), 1), "s")


    message("Setting an effect allele based alternate column for wide join")

    # Data$COORD_Alt <- stringi::stri_c("chr", Data$CHROM, ":", Data$GENPOS, ":", Data$ALLELE1, ":", Data$ALLELE0)

    # Combined_Processed_Data <- rbind(Combined_Processed_Data, Data)

    Combined_Processed_Data <- data.table::rbindlist(list(Combined_Processed_Data, Data), use.names = TRUE, fill = TRUE)

    message(paste0("Processed: ", corresponding_name))

      }


   if(Match_Allele_Direction == T)

   {

    message(paste0("Processing of matching allele directions to: ", Match_Allele_Study))

    Match_Allele_Study_Clean <- Match_Allele_Study

      # Ensure STUDY column has no leading/trailing spaces

      Combined_Processed_Data$STUDY <- trimws(Combined_Processed_Data$STUDY)

      message("Generating reference alleles from aforementioned data...")

        # Extract reference alleles, matching STUDY even if the number prefix differs

        reference <- Combined_Processed_Data %>%

          dplyr::mutate(STUDY_Clean = STUDY) %>%  # Remove number prefix
          dplyr::filter(STUDY_Clean == Match_Allele_Study_Clean) %>%  # Compare cleaned names
          dplyr::select(ID, ALLELE0, ALLELE1, COORD_Norm, COORD_Alt) %>%
          dplyr::rename(Ref_ALLELE0 = ALLELE0, Ref_ALLELE1 = ALLELE1)

        reference <- reference %>%
          dplyr::rename_with(~ paste0("REF_", .), everything())

        message("Aligning reference alleles with data...")

        #  Combined_Processed_Data2 <- Combined_Processed_Data %>%
        #  dplyr::left_join(reference, by = c("COORD_Norm" = "REF_COORD_Norm"))

        idx <- match(Combined_Processed_Data$COORD_Norm, reference$REF_COORD_Norm)

        ref_cols_to_add <- setdiff(names(reference), "REF_COORD_Norm")

        Combined_Processed_Data2 <- Combined_Processed_Data

        for (nm in ref_cols_to_add) {

          Combined_Processed_Data2[[nm]] <- reference[[nm]][idx]

        }

        message("Finding alleles that failed to match at first attempt...")

        # Step 2: Identify rows with no match in COORD_Norm

        unmatched <- Combined_Processed_Data2 %>%

          dplyr::filter(is.na(REF_ID)) %>%

          dplyr::select(names(Combined_Processed_Data))  # Retain only original columns for the fallback join not joined ones too.

          message("Trying to align reference to failures...")

          # Step 3: Perform the second join using COORD_Alt

          # fallback <- unmatched %>%

          # dplyr::left_join(reference, by = c("COORD_Norm" = "REF_COORD_Alt"))

        # match unmatched$COORD_Norm to reference$REF_COORD_Alt

        idx_alt <- match(unmatched$COORD_Norm, reference$REF_COORD_Alt)

        # add all reference columns except the join key (REF_COORD_Alt)

        ref_cols_to_add_alt <- setdiff(names(reference), "REF_COORD_Alt")

        fallback <- unmatched

        for (nm in ref_cols_to_add_alt) {

          fallback[[nm]] <- reference[[nm]][idx_alt]

        }

        message("Removing failures from initial success...and binding second attempt of failures...")

          # Step 4: Combine matched rows from the first and second join

          Combined_Processed_Data2 <- Combined_Processed_Data2 %>%

            dplyr::filter(!is.na(REF_ID)) %>%  # Keep rows matched on COORD_Norm
            dplyr::bind_rows(fallback)  # Add rows matched on COORD_Alt

          Combined_Processed_Data <- Combined_Processed_Data2

          message("Reversing BETA values of alleles where REF and ALT match reference study REF & ALT but are flipped...")

          Combined_Processed_Data <- Combined_Processed_Data %>%

            dplyr::mutate(

              BETA_Flipped = !is.na(REF_Ref_ALLELE0) & !is.na(REF_Ref_ALLELE1) &
                ALLELE0 == REF_Ref_ALLELE1 & ALLELE1 == REF_Ref_ALLELE0,
              BETA = dplyr::if_else(BETA_Flipped, BETA * -1, BETA),
              No_Match = (is.na(REF_COORD_Norm) & is.na(REF_COORD_Alt))

            )

          message("Unified ID dervies from the reference one as alleles matched to this direction")

          # View only the rows where BETA was flipped

          flipped_rows <- Combined_Processed_Data %>%
            dplyr::filter(BETA_Flipped)

          message(paste0("Number of alleles flipped across any study: ", nrow(flipped_rows)))

          # Filter only the rows where BETA was flipped

          unmatched_rows_2 <- Combined_Processed_Data %>%
            dplyr::filter(No_Match)

          message(paste0("Number of alleles which were not present either aligned or flipped from any study in the reference study: ", nrow(unmatched_rows_2)))

          # Make a unified ID showing direction of effect for all based on reference picked for allele flips

          message("Creating unified IDs for Flips, based on reference study direction of effect")

          # Combined_Processed_Data$Dir <- stringi::stri_c("chr", Combined_Processed_Data$CHROM, ":", Combined_Processed_Data$GENPOS, ":", Combined_Processed_Data$REF_Ref_ALLELE0, ":", Combined_Processed_Data$REF_Ref_ALLELE1)

          n <- nrow(Combined_Processed_Data)

          ch  <- Combined_Processed_Data[["CHROM"]]
          pos <- Combined_Processed_Data[["GENPOS"]]
          r0  <- Combined_Processed_Data[["REF_Ref_ALLELE0"]]
          r1  <- Combined_Processed_Data[["REF_Ref_ALLELE1"]]

          out <- character(n)

          bs <- 2e6
          n_chunks <- ceiling(n / bs)

          .progress_break({

            t0 <- Sys.time()

            for (k in seq_len(n_chunks)) {

              i <- (k - 1L) * bs + 1L
              j <- min(k * bs, n)

              out[i:j] <- stringi::stri_paste(

                "chr", ch[i:j], ":", pos[i:j], ":", r0[i:j], ":", r1[i:j],
                sep = ""

              )

              base::message(sprintf(

                "Creating Dir for %s | chunk %d/%d (rows %d–%d) | %.1fs",
               corresponding_name , k, n_chunks, i, j,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))

              ))

            }

          }, msg = paste0("Building Dir (final direction coordinate) in ", n_chunks, " chunks:"))

          Combined_Processed_Data$Dir <- out

          # Safe and makes sure all exact from ref - it is also matched to same allele direction within this already

          Combined_Processed_Data$COORD_Uni <- Combined_Processed_Data$REF_ID # every df has

  }

  # SOME NAs obviously if no match ie. ref like UKB not in other one

  # Remove rows where matched for allele flip, but actually strand change as ID is different

  # Issue if there were RS possibly but okay for now as all HG38 COORDS - RS still fine though for creating custom dir ID

    if(Match_Allele_Study == TRUE)

    {

      Combined_Processed_Data <- Combined_Processed_Data[!is.na(Combined_Processed_Data$ID) & !is.na(Combined_Processed_Data$REF_ID) & Combined_Processed_Data$ID == Combined_Processed_Data$REF_ID, ]

      message("Initial filtering completed")

    # Also if ref ID didn't match/is NA means these IDs are not in the ref study so drop

    Combined_Processed_Data <- Combined_Processed_Data %>%

    dplyr::filter(!is.na(REF_ID))

    }

    message("Splitting to separate dfs")

       # Split dataframe by unique STUDY values

        study_dfs <- lapply(split(Combined_Processed_Data, Combined_Processed_Data$STUDY), function(sub_df) {

        study_name <- unique(sub_df$STUDY)

        cols_to_rename <- !(colnames(sub_df) %in% c("COORD_Uni", "Dir"))

        # Rename only selected columns

        colnames(sub_df)[cols_to_rename] <- paste0(study_name, "_", colnames(sub_df)[cols_to_rename])

        # keep dir too

        return(sub_df)

      })

      message(paste0("Joining all to ", Match_Allele_Study, " ", "by unified ID" ))

      # Extract the first study's dataframe as the reference

      reference_study_name <- Match_Allele_Study  # Take the first study name
      reference_df <- study_dfs[[reference_study_name]]  # Extract first dataframe

      # Loop over remaining studies and left join them to the reference dataframe using COORD_Uni

      for (study_name in setdiff(names(study_dfs), reference_study_name)) {

        message(paste0("Joining ", study_name, " to ", reference_study_name))

        coord_col <- grep("COORD_Uni", colnames(study_dfs[[study_name]]), value = TRUE)

        reference_df <- dplyr::left_join(reference_df, study_dfs[[study_name]], by = coord_col)

      }

     # Store the final merged dataframe

     Combined_Processed_Data_Joined <- reference_df

     genpos_cols <- grep("GENPOS$", colnames(Combined_Processed_Data_Joined), value = TRUE)

     # Filter genpos_cols to only keep those containing Match_Allele_Study

     genpos_cols <- genpos_cols[grepl(Match_Allele_Study, genpos_cols)]

     # Could allow to keep all, remove any blank in future like with the below code originally....

     # Need to adjust in future for non-regenie X, is it 23, X, chr x etc.

     # Could add combo cols below to generate files with all combos of studies to meta.

     # Drop rows where all GENPOS columns are NA

      if (length(genpos_cols) > 0) {

        Combined_Processed_Data_Joined <- Combined_Processed_Data_Joined %>%
          dplyr::filter(!dplyr::if_any(dplyr::all_of(genpos_cols), is.na))

      }

      # Identify columns that end exactly with "BETA" or "SE" plus "COORD_Uni"

      keep_cols <- colnames(Combined_Processed_Data_Joined)[

        grepl(paste0(Output, "$", collapse = "|"), colnames(Combined_Processed_Data_Joined)) |

        colnames(Combined_Processed_Data_Joined) == "COORD_Uni" | colnames(Combined_Processed_Data_Joined) == "Dir.x"

      ]

      Ref_Study_Identifier <- paste0(Match_Allele_Study, "_ID")

      if(Output_Reference_Identifier == TRUE)

      {

      keep_cols <- c(keep_cols, Ref_Study_Identifier) # have to add _ID to pull out

      }

      # ID is always renamed reference file one for original identifier

      message("These column types will be kept from each data frame and married to unified ID in order:")

      # Preserve the original order found in the dataframe

      ordered_cols <- colnames(Combined_Processed_Data_Joined)[colnames(Combined_Processed_Data_Joined) %in% keep_cols]

      # Ensure COORD_Uni is always first while maintaining the original order for the rest

      ordered_cols <- c("COORD_Uni", setdiff(ordered_cols, "COORD_Uni"))

      message("Filtering columns")

      # Select columns in the correct order

      Combined_Processed_Data_Joined <- Combined_Processed_Data_Joined %>%
        dplyr::select(all_of(ordered_cols))

      Combined_Processed_Data_Joined$ID_DIR <- Combined_Processed_Data_Joined$Dir.x
      Combined_Processed_Data_Joined$Dir.x <- NULL
      Combined_Processed_Data_Joined$COORD_Uni <- NULL

      message("Final Output")

      Combined_Processed_Data_Joined <- Combined_Processed_Data_Joined[,

        c("ID_DIR", setdiff(names(Combined_Processed_Data_Joined), "ID_DIR"))

      ]

      return(Combined_Processed_Data_Joined)

}

  if (!exists("use_wrapper")) use_wrapper <- TRUE

  if (isTRUE(use_wrapper)) {

    .METASOFT_File_Gen_original <- METASOFT_File_Gen

    METASOFT_File_Gen <- function(..., session = NULL, Debug = FALSE) {

      args <- list(...)

      if (!is.null(args$Data)) {

        labs <- character(0)

        if (is.character(args$Data)) labs <- unname(args$Data)

        if (length(labs)) message("Processing datasets: ", paste(labs, collapse = ", "))

      }

      verbose_mode <- isTRUE(args$Verbose)

      # Verbose: run the real function directly (no progress wrapper)

      if (verbose_mode) {

        if (Debug) message("[Debug] METASOFT_File_Gen: Verbose=TRUE (no progress wrapper)")

        return(do.call(.METASOFT_File_Gen_original, args))

      }

      # Non-verbose: run with progress counter

      # Pass the 'real' function + args (no runner indirection)

      if (Debug) message("[Debug] METASOFT_File_Gen: Verbose=FALSE (run_with_counter)")

      out <- suppressWarnings(

        run_with_counter(

          func    = .METASOFT_File_Gen_original,
          args    = args,
          session = session

        )

      )

      if (Debug) message("[Debug] METASOFT_File_Gen: finished")

      out

    }

  }


