
#' Munged summary statistics in preparation for METASOFT-like meta-analysis
#'
#' @param Match_Allele_Direction Match the allele directions to avoid effect direction flips; defaults to TRUE
#' @param Match_Allele_Study Data set used as reference point for effect allele direction; defaults to NULL - automatically assigned as first set
#' @param Output A selection of columns to isolate, in order, alternating by dataset; defaults to c("BETA", "SE") - standard METASOFT ordered format
#' @param Data These are the summary statistics (either file path(s), environmental object(s) or a variable(s) containing such contents) to be plotted; defaults to NULL
#' @param BETA_Column Manually specify BETA Test Statistic column; defaults to NULL, leading to auto-detection
#' @param Position_Column Manually specify chromosome genomic position column; defaults to NULL, leading to auto-detection
#' @param SNP_ID_Column Manually specify SNP ID column; defaults to NULL, leading to auto-detection
#' @param PValue_Column Manually specify P Value column; defaults to NULL, leading to auto-detection
#' @param Reference_Allele_Column Manually specify reference allele column; defaults to NULL, leading to auto-detection
#' @param Effect_Allele_Column Manually specify effect allele column; defaults to NULL, leading to auto-detection
#' @param SE_Column Manually specify SE of Test Statistic column; defaults to NULL, leading to auto-detection
#'
#' @return Munged dataframe amenable to METASOFT meta-analysis - only including SNPs present in backboned file
#' @export
#'
#' @examples METASOFT <- METASOFT_File_Gen(Data = c("Intelligence_Sum_Stats", "Household_Income_Sum_Stats"))
#'


METASOFT_File_Gen <- function(Data = NULL,
                              BETA_Column = NULL,
                              Position_Column = NULL,
                              Chromosome_Column = NULL,
                              SNP_ID_Column = NULL,
                              PValue_Column = NULL,
                              Reference_Allele_Column = NULL,
                              Effect_Allele_Column = NULL,
                              SE_Column = NULL,
                              Match_Allele_Direction = TRUE,
                              Match_Allele_Study = NULL,
                              Output = c("BETA", "SE"),
                              Verbose = FALSE

)


{

      Data_Sets <- Data


      if (!is.null((Data_Sets))) {

        message("Using Data Set Names")
        Names <- (Data_Sets)  # Use names of the Data_Sets vector

      }

      message(Names)

      if (is.null(Match_Allele_Study)) {

        Match_Allele_Study <- Names[1]

        message(paste0("Matching allele directions automatically to ", Match_Allele_Study) )

      }

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

      for (Set in Names) {
        message(Set)

        Spec_Name <- Set

        message(Spec_Name)

        if ((is.character(Spec_Name) && grepl("/", Spec_Name)) ||
            (is.character(Spec_Name) && file.exists(Spec_Name))) {

          # Treat as file path
          Data <- vroom::vroom(Spec_Name, show_col_types = FALSE)

        } else if (exists(Spec_Name, envir = .GlobalEnv)) {

          # Treat as object in the environment
          Data <- get(Spec_Name, envir = .GlobalEnv)

        } else {

          stop(sprintf("Spec_Name '%s' is not an object or a valid file path", Spec_Name))
        }

        Loaded_Data[[Spec_Name]] <- Data

        Chromosome_Columns[[Spec_Name]] <- detect_chromosome_column(Data, Chromosome_Column)

        if (!is.null(Chromosome_Columns[[Spec_Name]])) {
          Data <- Data %>%
            dplyr::mutate(!!Chromosome_Columns[[Spec_Name]] :=
                            ifelse(.data[[Chromosome_Columns[[Spec_Name]]]] == "23", "X",
                                   .data[[Chromosome_Columns[[Spec_Name]]]]))
        }

        PValue_Columns[[Spec_Name]]     <- detect_pvalue_column(Data, PValue_Column)
        Position_Columns[[Spec_Name]]   <- detect_position_column(Data, Position_Column)
        SNP_ID_Columns[[Spec_Name]]     <- detect_snp_column(Data, SNP_ID_Column)
        Ref_Allele_Columns[[Spec_Name]] <- detect_reference_allele_column(Data, Reference_Allele_Column)
        Alt_Allele_Columns[[Spec_Name]] <- detect_effect_allele_column(Data, Effect_Allele_Column)
        BETA_Columns[[Spec_Name]]       <- detect_beta_column(Data, BETA_Column)
        SE_Columns[[Spec_Name]]         <- detect_se_column(Data, SE_Column)
      }

    Combined_Processed_Data <- data.frame()

    for (i in seq_along(Data_Sets)) {

    dataset_name <- Data_Sets[i]  # Get the dataset name

    message(dataset_name)

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
    message(paste("Processing dataset:", dataset_name))
    message(paste("Corresponding name:", corresponding_name))

    #Manually assign columns for ease of use
    Data$CHROM <- Data[[Chromosome_Column]]
    Data$GENPOS <- Data[[Position_Column]]
    Data$ID <- Data[[SNP_ID_Column]]
    Data$P <- Data[[PValue_Column]]
    Data$SE <- Data[[Standard_Error_Column]]
    Data$BETA <- Data[[Beta_Column]]
    Data$STUDY <- corresponding_name
    Data$ALLELE0 <- toupper(Data[[Reference_Allele]])
    Data$ALLELE1 <- toupper(Data[[Effect_Allele]])

    Data <- Data %>% dplyr::select(ID, ALLELE0, ALLELE1, CHROM, GENPOS, BETA, SE, P, STUDY)

    Data$COORD_Norm <- stringi::stri_c("chr", Data$CHROM, ":", Data$GENPOS, ":", Data$ALLELE0, ":", Data$ALLELE1)
    Data$COORD_Alt <- stringi::stri_c("chr", Data$CHROM, ":", Data$GENPOS, ":", Data$ALLELE1, ":", Data$ALLELE0)

    Combined_Processed_Data <- rbind(Combined_Processed_Data, Data)

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


        Combined_Processed_Data2 <- Combined_Processed_Data %>%
          dplyr::left_join(reference, by = c("COORD_Norm" = "REF_COORD_Norm"))


        message("Finding alleles that failed to match at first attempt...")

        # Step 2: Identify rows with no match in COORD_Norm
        unmatched <- Combined_Processed_Data2 %>%
          dplyr::filter(is.na(REF_ID)) %>%
          dplyr::select(names(Combined_Processed_Data))  # Retain only original columns for the fallback join not joined ones too.

          message("Trying to align reference to failures...")

          # Step 3: Perform the second join using COORD_Alt
          fallback <- unmatched %>%
            dplyr::left_join(reference, by = c("COORD_Norm" = "REF_COORD_Alt"))

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

          # View only the rows where BETA was flipped
          unmatched_rows_2 <- Combined_Processed_Data %>%
            dplyr::filter(No_Match)

          message(paste0("Number of alleles which were not present either aligned or flipped from any study in the reference study: ", nrow(unmatched_rows_2)))

          #Make a unified ID showing direction of effect for all based on reference picked for allels flips

          message("Creating unified IDs for Flips, based on reference study direction of effect")

          Combined_Processed_Data$Dir <- stringi::stri_c("chr", Combined_Processed_Data$CHROM, ":", Combined_Processed_Data$GENPOS, ":", Combined_Processed_Data$REF_Ref_ALLELE0, ":", Combined_Processed_Data$REF_Ref_ALLELE1)

          #Safe and makes sure all exact from ref - it is also matched to same allele direction within this already.
          Combined_Processed_Data$COORD_Uni <- Combined_Processed_Data$REF_ID

          }

  #SOME NAs obvisouly if no match ie. ref like UKB not in other one
  #Remove rows where matched for allele flip, but actually strand change as ID is different.
  #Issue if there were RS possibly but okay for now as all HG38 COORDS

  #Must be the same
  Combined_Processed_Data <- Combined_Processed_Data[!is.na(Combined_Processed_Data$ID) & !is.na(Combined_Processed_Data$REF_ID) & Combined_Processed_Data$ID == Combined_Processed_Data$REF_ID, ]

  message("Initial filtering completed")

  #Also if ref ID didnt match/is NA means these IDs are not in the ref study so drop

  Combined_Processed_Data <- Combined_Processed_Data %>%
    dplyr::filter(!is.na(REF_ID))



  message("Splitting to separate dfs")

       # Split dataframe by unique STUDY values
        study_dfs <- lapply(split(Combined_Processed_Data, Combined_Processed_Data$STUDY), function(sub_df) {
        study_name <- unique(sub_df$STUDY)

        cols_to_rename <- !(colnames(sub_df) %in% c("COORD_Uni", "Dir"))

        # Rename only selected columns
        colnames(sub_df)[cols_to_rename] <- paste0(study_name, "_", colnames(sub_df)[cols_to_rename])
        #keep dir too

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

     message("Only excluding if missing value for this/backbone:")

     message("Removing SNPs with missing data from backbone/allele flip study...")

     #could allow to keep all, remove any blank in futue like with the below code originally....
     #Need to adjust in future for non-regenie X, is it 23, X, chr x etc.
     #Could add combo cols below to generate files with all combos of studies to meta.

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

      message("These column types will be kept from each data frame and married to unified ID in order:")

      # Preserve the original order found in the dataframe
      ordered_cols <- colnames(Combined_Processed_Data_Joined)[colnames(Combined_Processed_Data_Joined) %in% keep_cols]

      # Ensure COORD_Uni is always first while maintaining the original order for the rest
      ordered_cols <- c("COORD_Uni", setdiff(ordered_cols, "COORD_Uni"))

      message("Therefore these columns will be kept in saved output:")

      message(ordered_cols)

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

if(use_wrapper == TRUE)
{
.METASOFT_File_Gen_original <- METASOFT_File_Gen

METASOFT_File_Gen <- function(..., session = NULL, Debug = FALSE) {
  args <- list(...)

  show_dataset_message <- function(x) {

    labs <- character(0)
    if (!is.null(x)) {

      if (is.character(x)) {

        labs <- unname(x)
      } else if (!is.null(names(x)) && any(nzchar(names(x)))) {

        labs <- unname(names(x))

      }

    }

    if (length(labs)) {

      message("Processing datasets: ", paste(labs, collapse = ", "))

    }
  }
  show_dataset_message(args$Data)

  verbose_mode <- isTRUE(args$Verbose)

  if (verbose_mode) {

    # Direct call, show everything
    return(do.call(.METASOFT_File_Gen_original, args))

  } else {

    # Silent/progress mode
    runner <- function() {

      if (Debug) message("[Debug] METASOFT_File_Gen: start")

      res <- do.call(.METASOFT_File_Gen_original, args)

      if (Debug) message("[Debug] METASOFT_File_Gen: end")

      res

    }

    on.exit({

      if (Debug) message("[Debug] METASOFT_File_Gen: wrapper on.exit() reached")

    }, add = TRUE)

    return(

      suppressMessages(
        suppressWarnings(
          run_with_counter(
            func    = runner,
            args    = list(),
            session = session

          )
        )
      )
    )
  }
}
}

