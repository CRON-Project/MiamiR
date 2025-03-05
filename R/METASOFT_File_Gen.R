
#Combined <- METASOFT_File_Gen(Locations)
#Locations <- list.files("C:/Users/callumon/Downloads", full.names = T)
#Locations <- Locations[grepl("Savage|Hill", Locations)]


#' Title Create datasets of values of interest across multiple GWAS datasets for meta-analysis
#'
#' @param Data_Sets These are a list of GWAS summary statistics to be evaluated; defaults to empty list.
#' @param Chromosome_Columns These are a list of manual chromosome column names for the data sets used in the order specified; defaults to c()
#' @param Reference_Alleles These are a list of manual reference allele column names for the data sets used in the order specified; defaults to c()
#' @param Effect_Alleles These are a list of manual effect allele column names for the data sets used in the order specified; defaults to c()
#' @param Position_Columns These are a list of manual position column names for the data sets used in the order specified; defaults to c()
#' @param SNP_ID_Columns These are a list of manual SNP ID column names for the data sets used in the order specified; defaults to c()
#' @param Beta_Columns These are a list of manual BETA column names for the data sets used in the order specified; defaults to  c()
#' @param Standard_Error_Columns These are a list of manual SE column names for the data sets used in the order specified; defaults to c()
#' @param PValue_Columns These are a list of manual P column names for the data sets used in the order specified; defaults to c()
#' @param Spec_CHROM Particular Chromosome to select for a subsetted analysis
#' @param Match_Allele_Direction Do you want to match the allele directions to avoid flips (T/F); defaults to T
#' @param Match_Allele_Study If the above is T, which data sets should be the reference point for effect allele direction; defaults to "LbDementia_Sum_Stats"; THIS SHOULD BE THE NAME YOU ASSIGN
#' @param Output A selection of columns to isolate, in order, alternating by dataset; defaults to c(BETA, SE)
#' @param File_NameFile name to save data as; defaults to "METASOFT_Format"
#' @param File_Type File type of saved data; defaults to "txt"
#'
#' @return Data frame of information combined across the datasets
#' @export
#'
#' @examples Combined <- METASOFT_File_Gen(Data_Sets = c("Intelligence_Peaks"))
#'
#'
#'
#'


METASOFT_File_Gen <- function(Data_Sets = c(),
                        Chromosome_Columns = c(),
                        Reference_Alleles = c(),
                        Effect_Alleles = c(),
                        Position_Columns = c() , SNP_ID_Columns = c(),
                        Beta_Columns = c(), Standard_Error_Columns = c(),
                        PValue_Columns = c(),
                        Spec_CHROM = NULL,
                        Match_Allele_Direction = TRUE,
                        Match_Allele_Study = NULL,
                        Output = c("BETA", "SE"),
                        File_Name = "METASOFT_Format",
                        File_Type = "txt"
)




{


  # Check if Data_Sets contains file paths
  if (all(file.exists(Data_Sets))) {
    message("Loading datasets from file paths...")

    dataset_names <- c()  # Store dataset names

    for (path in Data_Sets) {
      # Extract filename without extension
      dataset_name <- tools::file_path_sans_ext(basename(path))


      message("Processing file: ", path)
      message("Extracted dataset name: ", dataset_name)

      # Detect if the file is compressed
      is_compressed <- grepl("\\.gz$", path, ignore.case = TRUE)

      # Read only the first few rows to detect column names
      header_df <- tryCatch({
        if (is_compressed) {

          #data.table::fread(cmd = paste0("zcat ", path, " | head -n 5"))
          read.csv(path, sep = "", nrows = 5)

        } else {
          data.table::fread(path, nrows = 5)
        }
      }, error = function(e) {
        message(paste("Skipping file due to error:", path, "\n", e))
        next  # Skip to next file instead of stopping
      })

      # **Ensure column names are standardized**
      colnames(header_df) <- trimws(colnames(header_df))  # Remove spaces

      # Detect chromosome column dynamically
      allowed_names_chromosomes <- c("Chromosome", "chromosome", "chrom", "chr", "CHROM", "CHR", "Chr", "Chrom")

      chrom_col <- NULL  # Initialize as NULL

      for (allowed_name in allowed_names_chromosomes) {
        if (allowed_name %in% colnames(header_df)) {
          chrom_col <- allowed_name  # Ensure exact match
          message(paste("✅ Detected chromosome column:", chrom_col))
          break
        }
      }

      # If no Chromosome column is found, continue without filtering
      if (is.null(chrom_col)) {
        message(paste("⚠️ No chromosome column found in", dataset_name, "- loading full dataset."))
      }

      # **Optimized Direct Filtering While Loading**
      df <- tryCatch({
        if (is_compressed) {
          message(paste0("📂 Detected compressed file, using read.csv() to load and filtering CHROM == ", Spec_CHROM, "."))
          #  temp_df <- data.table::fread(cmd = paste0("zcat ", path, " | grep '^1\\s'"), sep = "\t")
          temp_df <- read.csv(path, sep = "")
      #    print("CHROM Filt...")
          temp_df <- temp_df[temp_df[[chrom_col]] == Spec_CHROM, ]
        } else {
          if (!is.null(Spec_CHROM)) {
            message(paste0("📂 Detected average file, using vroom() to load and filtering CHROM == ", Spec_CHROM, "."))
          } else {
            message("📂 Detected average file, using vroom() to load. All CHROMs maintained.")
          }
          temp_df <- vroom::vroom(path)  # Read everything first (unavoidable)

         # print(temp_df)
        #  print(chrom_col)

          # **Apply filtering right after reading**
          if (!is.null(chrom_col)) {
            # temp_df <- temp_df[temp_df$Chr == 1,]
            if(is.null(Spec_CHROM))
            {
              print("Whole file desired")
            }else{
          #    print("CHROM Filt...")
              temp_df <- temp_df[temp_df[[chrom_col]] == Spec_CHROM, ]
            }
          }
        }

        if (nrow(temp_df) == 0) {
          message(paste("⚠️ No rows found with", chrom_col, "== 1 in", dataset_name, "- skipping dataset."))
        }

        temp_df  # Return filtered dataset

      }, error = function(e) {
        message(paste("Skipping file due to error:", path, "\n", e))
        next  # Skip to next file instead of stopping
      })

      # Store dataset in environment
      assign(dataset_name, df, envir = .GlobalEnv)
      message(paste("✅ Dataset", dataset_name, "loaded into environment."))

      # Print the first few rows of the dataset
   #   message(paste("🔍 Preview of", dataset_name, ":"))
      #print(head(df))  # Show first few rows of the dataset

      # Append dataset name to vector
      dataset_names <- c(dataset_names, dataset_name)
    }

    # Update Data_Sets with dataset names only
    Data_Sets <- dataset_names

  } else {
    message("Using datasets from the R environment...")
    Data_Sets <- lapply(Data_Sets, function(name) {
      if (!is.character(name)) stop("Dataset name must be a character string.")

      if (exists(name, envir = .GlobalEnv)) {
        return(name)  # Keep dataset name in Data_Sets
      } else {
        stop(paste("Dataset", name, "not found in environment."))
      }
    })

    # Convert list to character vector
    Data_Sets <- unlist(Data_Sets)
  }

  message("📌 Final Data_Sets list: ", paste(Data_Sets, collapse = ", "))

  # Return Data_Sets with correct names
  #  return(Data_Sets)


  #z


  if (!is.null((Data_Sets))) {
    print("Using Data Set Names")
    Names <- (Data_Sets)  # Use names of the Data_Sets vector
  } else {
    print("Auto-naming")
    Names <- paste("Dataset", seq_along(Data_Sets))  # Fallback to generic names
  }

  #print(Names)



  #Goes below above as need to match based on STUDY assigned from Name
  if (is.null(Match_Allele_Study)) {

    Match_Allele_Study <- Names[1]
    print(paste0("Matching allele directions automatically to ", Match_Allele_Study) )

   # print(Match_Allele_Study)

    #    if (num_pattern_count > 2) {
    #
    #     Match_Allele_Study <- sub("^[0-9]+_", "", Match_Allele_Study)
    #    }

    # print(Match_Allele_Study)

  }


  Combined_Processed_Data <- data.frame()



  # print(Data_Sets)

  for (i in seq_along(Data_Sets)) {

    #   print(Data_Sets)



    # dataset_name <- Data_Sets[i]
    # Data <- dataset_name

    dataset_name <- Data_Sets[i]  # Get the dataset name

    #  print(dataset_name)

    Data <- get(dataset_name)  # Load the dataset using get()

    #    print(Data)

    corresponding_name <- Names[i]  # Get the corresponding name


    print(paste0("Processing: ", corresponding_name))


    Chromosome_Column <- Chromosome_Columns[i]
    Position_Column <- Position_Columns[i]
    SNP_ID_Column <- SNP_ID_Columns[i]

    PValue_Column <- PValue_Columns[i]
    Standard_Error_Column <- Standard_Error_Columns[i]


    Beta_Column <- Beta_Columns[i]


    Reference_Allele <- Reference_Alleles[i]
    Effect_Allele <- Effect_Alleles[i]



    # Print dataset name and corresponding name
    print(paste("Processing dataset:", dataset_name))
    print(paste("Corresponding name:", corresponding_name))


    if(is.null(Reference_Allele))
    {

      allowed_names_reference_allele <- c("A0", "A2", "REF", "other_allele", "a0", "a2", "ref", "Non_effect_Allele", "Reference", "reference", "allele0", "allele2", "ALLELE0", "ALLELE2")

      for (allowed_name in allowed_names_reference_allele) {
        if (allowed_name %in% colnames(Data)) {
          usable_ref_allele <- colnames(Data)[which(colnames(Data) == allowed_name)]
          break
        }
      }



      print(paste0("Using", " ", allowed_name, " ", "as reference allele (automatic)"))
      Reference_Allele <- allowed_name

    }else{
      print(paste0("Using", " ", Reference_Allele, " ", "as reference allele (manual)"))
    }


    if(is.null(Effect_Allele))
    {

      allowed_names_alternate_allele <- c("A1", "ALT", "a1", "alt", "Alternate", "Effect_Allele", "effect_allele", "alternate", "allele1", "ALLELE1")

      for (allowed_name in allowed_names_alternate_allele) {
        if (allowed_name %in% colnames(Data)) {
          usable_effect_allele <- colnames(Data)[which(colnames(Data) == allowed_name)]
          break
        }
      }



      print(paste0("Using", " ", allowed_name, " ", "as effect allele (automatic)"))
      Effect_Allele <- allowed_name


    }else{
      print(paste0("Using", " ", Effect_Allele, " ", "as effect allele (manual)"))
    }




    if(is.null(Chromosome_Column))
    {

      allowed_names_chromosomes <- c("chromosome", "chrom", "chr", "CHROM", "Chromosome", "CHR", "Chr", "Chrom")

      for (allowed_name in allowed_names_chromosomes) {
        if (allowed_name %in% colnames(Data)) {
          usable_chrom_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
          break
        }
      }


      print(paste0("Using", " ", allowed_name, " ", "as chromosome column (automatic)"))
      Chromosome_Column <- allowed_name


    }else{
      print(paste0("Using", " ", Chromosome_Column, " ", "as chromosome column (manual)"))
    }






    if(is.null(Beta_Column))
    {


      allowed_names_betas <- c("BETA", "Beta", "beta", "B", "stdBeta")



      for (allowed_name in allowed_names_betas) {
        if (allowed_name %in% colnames(Data)) {
          usable_chrom_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
          break
        }
      }


      print(paste0("Using", " ", allowed_name, " ", "as beta column (automatic)"))
      Beta_Column <- allowed_name



    }else{
      print(paste0("Using", " ", Beta_Column, " ", "as beta column (manual)"))
    }




    if(is.null(Standard_Error_Column))
    {


      allowed_names_ses<- c("SE", "se", "Se", "Standard_Error", "standard_error", "Standard_Error_of_Beta")

      for (allowed_name in allowed_names_ses) {
        if (allowed_name %in% colnames(Data)) {
          usable_chrom_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
          break
        }
      }


      print(paste0("Using", " ", allowed_name, " ", "as standard error column (automatic)"))
      Standard_Error_Column <- allowed_name



    }else{
      print(paste0("Using", " ", Standard_Error_Column, " ", "as standard error column (manual)"))
    }



    if(is.null(PValue_Column))
    {
      allowed_names_pvals <- c("P", "p", "Pvalue", "pvalue", "P-Value", "p-value", "p-Value",
                               "P-VALUE", "logp","LogP", "LOGP", "Logp", "log10p","Log10P",
                               "LOG10P", "Log10p",
                               "log10p", "LOG10P", "-LOG10P", "", "p_value")

      for (allowed_name in allowed_names_pvals) {
        if (allowed_name %in% colnames(Data)) {
          usable_p_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
          break
        }
      }


      print(paste0("Using", " ", allowed_name, " ", "as P-Val column (automatic)"))
      PValue_Column <- allowed_name

    }else{
      print(paste0("Using", " ", PValue_Column, " ", "as P-Val column (manual)"))
    }



    if (!("P" %in% colnames(Data)) & any(colnames(Data) %in% c("logp", "LogP", "LOGP", "Logp",
                                                               "log10p", "Log10P", "LOG10P",
                                                               "Log10p", "-LOG10P"))) {

      print("No P Value in first dataset, Calculating from LOG10P column detected")
      Data$P <- 10^-(as.numeric(Data[[PValue_Column]]))
      PValue_Column <- "P"

    }



    if(is.null(Position_Column))
    {

      allowed_names_pos <- c("POS", "pos", "Pos", "BP", "BPos", "bpos", "BPOS", "bPos",
                             "Position", "position", "POSITION", "genpos", "GENPOS",
                             "Genpos", "base_pair_location")


      for (allowed_name in allowed_names_pos) {
        if (allowed_name %in% colnames(Data)) {
          usable_pos_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
          break
        }
      }


      print(paste0("Using", " ", allowed_name, " ", "as position column (automatic)"))
      Position_Column <- allowed_name

    }else{
      print(paste0("Using", " ", Position_Column, " ", "as position column (manual)"))
    }


    if(is.null(SNP_ID_Column))
    {


      allowed_names_SNP <- c("ID", "Id", "ID", "RsID", "RsId","RSID", "snp", "SNP", "Snp",
                             "snv" ,"SNV" , "Snv", "RS", "rs", "variant_id" )

      for (allowed_name in allowed_names_SNP) {
        if (allowed_name %in% colnames(Data)) {
          usable_snp_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
          break
        }
      }


      print(paste0("Using", " ", allowed_name, " ", "as SNP column (automatic)"))
      SNP_ID_Column <- allowed_name

    }else{
      print(paste0("Using", " ", SNP_ID_Column, " ", "as SNP column (manual)"))
    }



    #Manually assign columns for ease of use
    Data$CHROM <- Data[[Chromosome_Column]]
    Data$GENPOS <- Data[[Position_Column]]
    Data$ID <- Data[[SNP_ID_Column]]





    Data$P <- Data[[PValue_Column]]
    Data$SE <- Data[[Standard_Error_Column]]


    Data$BETA <- Data[[Beta_Column]]

    Data$STUDY <- corresponding_name


    #If from non-REGENIE won't match
    # Assign and convert to uppercase
    Data$ALLELE0 <- toupper(Data[[Reference_Allele]])
    Data$ALLELE1 <- toupper(Data[[Effect_Allele]])




    #sill fine if for OR as renaming BETA

    # print(Data)

    Data <- Data %>% dplyr::select(ID, ALLELE0, ALLELE1, CHROM, GENPOS, BETA, SE, P, STUDY)


 #   print("Naming")

   # colnames(Data) <- paste0(corresponding_name, colnames(Data))

 #   print("Named")

    Data$COORD_Norm <- stringi::stri_c("chr", Data$CHROM, ":", Data$GENPOS, ":", Data$ALLELE0, ":", Data$ALLELE1)
    Data$COORD_Alt <- stringi::stri_c("chr", Data$CHROM, ":", Data$GENPOS, ":", Data$ALLELE1, ":", Data$ALLELE0)

    #think if strands in future and how to keep? Also, must have same ID though.


    Combined_Processed_Data <- rbind(Combined_Processed_Data, Data)


    print(paste0("Processed: ", corresponding_name))

 #   print("Combined")

    #No longer need once joined

    print("Removing Dataframes Loaded into the Environment to save RAM")
#
#      # if (exists(corresponding_name)) {
#      #   rm(list = corresponding_name, envir = .GlobalEnv)
#      #   cat("Dataset", corresponding_name, "removed from the environment.\n")
#      # } else {
#      #   cat("Dataset", corresponding_name, "does not exist in the environment.\n")
#      # }
#
#     #

  }

 # print(Combined_Processed_Data)







  if(Match_Allele_Direction == T)

  {


  #    print(Combined_Processed_Data)

   #   print(Combined_Processed_Data$STUDY)

    print(paste0("Processing of matching allele directions to: ", Match_Allele_Study))





        Match_Allele_Study_Clean <- Match_Allele_Study


      # Ensure STUDY column has no leading/trailing spaces
      Combined_Processed_Data$STUDY <- trimws(Combined_Processed_Data$STUDY)



      print("Generating reference alleles from aforementioned data...")

        # Extract reference alleles, matching STUDY even if the number prefix differs
        reference <- Combined_Processed_Data %>%
          dplyr::mutate(STUDY_Clean = STUDY) %>%  # Remove number prefix
          dplyr::filter(STUDY_Clean == Match_Allele_Study_Clean) %>%  # Compare cleaned names
          dplyr::select(ID, ALLELE0, ALLELE1, COORD_Norm, COORD_Alt) %>%
          dplyr::rename(Ref_ALLELE0 = ALLELE0, Ref_ALLELE1 = ALLELE1)




        reference <- reference %>%
          dplyr::rename_with(~ paste0("REF_", .), everything())

#print(reference)


#print(Combined_Processed_Data)


        print("Aligning reference alleles with data...")


        print(Combined_Processed_Data)
        print(reference)



Combined_Processed_Data2 <- Combined_Processed_Data %>%
  dplyr::left_join(reference, by = c("COORD_Norm" = "REF_COORD_Norm"))


#print(sum(is.na(Combined_Processed_Data2$REF_Ref_ALLELE1)))
#print(sum(is.na(Combined_Processed_Data2$REF_Ref_ALLELE0)))


#return(Combined_Processed_Data2)

print("Finding alleles that failed to match at first attempt...")

# Step 2: Identify rows with no match in COORD_Norm
unmatched <- Combined_Processed_Data2 %>%
  dplyr::filter(is.na(REF_ID)) %>%
  dplyr::select(names(Combined_Processed_Data))  # Retain only original columns for the fallback join not joined ones too.





print("Trying to align reference to failures...")

#print(unmatched)


# Step 3: Perform the second join using COORD_Alt
fallback <- unmatched %>%
  dplyr::left_join(reference, by = c("COORD_Norm" = "REF_COORD_Alt"))

#print(fallback)


print("Removing failures from initial success...and binding second attempt of failures...")

# Step 4: Combine matched rows from the first and second join
Combined_Processed_Data2 <- Combined_Processed_Data2 %>%
  dplyr::filter(!is.na(REF_ID)) %>%  # Keep rows matched on COORD_Norm
  dplyr::bind_rows(fallback)  # Add rows matched on COORD_Alt


Combined_Processed_Data <- Combined_Processed_Data2


# Combined_Processed_Data <- Combined_Processed_Data  %>%
#       dplyr::mutate(
#           # Flip BETA only if alleles are swapped
#           BETA = dplyr::if_else(!is.na(Ref_ALLELE0) & !is.na(Ref_ALLELE1) &
#                                   ALLELE0 == Ref_ALLELE1 & ALLELE1 == Ref_ALLELE0,
#                                 BETA * -1, BETA)
#         )


print("Reversing BETA values of alleles where REF and ALT match reference study REF & ALT but are flipped...")

Combined_Processed_Data <- Combined_Processed_Data %>%
  dplyr::mutate(
    BETA_Flipped = !is.na(REF_Ref_ALLELE0) & !is.na(REF_Ref_ALLELE1) &
      ALLELE0 == REF_Ref_ALLELE1 & ALLELE1 == REF_Ref_ALLELE0,
    BETA = dplyr::if_else(BETA_Flipped, BETA * -1, BETA),
    No_Match = (is.na(REF_COORD_Norm) & is.na(REF_COORD_Alt))
  )


print("Unified ID dervies from the reference one as alleles matched to this direction")


# View only the rows where BETA was flipped
flipped_rows <- Combined_Processed_Data %>%
  dplyr::filter(BETA_Flipped)

print(paste0("Number of alleles flipped across any study: ", nrow(flipped_rows)))


#print(table(is.na(Combined_Processed_Data$REF_COORD_Norm) & is.na(Combined_Processed_Data$REF_COORD_Alt)))

#unmatched_rows <- Combined_Processed_Data %>%
 # dplyr::filter(is.na(Combined_Processed_Data$REF_COORD_Norm) & is.na(Combined_Processed_Data$REF_COORD_Alt))


#print(unmatched_rows)



# View only the rows where BETA was flipped
unmatched_rows_2 <- Combined_Processed_Data %>%
  dplyr::filter(No_Match)

#print(unmatched_rows_2)

print(paste0("Number of alleles which were not present either aligned or flipped from any study in the reference study: ", nrow(unmatched_rows_2)))





#Make a unified ID showing direction of effect for all based on reference picked for allels flips

print("Creating unified IDs for Flips, based on reference study direction of effect")

Combined_Processed_Data$COORD_Uni <- stringi::stri_c("chr", Combined_Processed_Data$CHROM, ":", Combined_Processed_Data$GENPOS, ":", Combined_Processed_Data$REF_Ref_ALLELE0, ":", Combined_Processed_Data$REF_Ref_ALLELE1)

}

    #  print(Combined_Processed_Data)

  return(Combined_Processed_Data)

#

      # Print final dataset
   #   print(Combined_Processed_Data)
  #    print(Match_Allele_Study)

  #    print(sum(is.na(Combined_Processed_Data$REF_Ref_ALLELE1)))
  #    print(sum(is.na(Combined_Processed_Data$REF_Ref_ALLELE0)))


      #Sometimes bp may match but alleles are different some blank means they weren't in the reference gwas e.g. blank
      #savage extra cols means that SNP wasn't in Hill used as ref so couldn't bind.
      #CHROM may also be dif

      #ALSO if some studies missing for SNP what can the flip be matched to?

  #    Combined_Processed_Datamiss <- Combined_Processed_Data[is.na(Combined_Processed_Data$REF_Ref_ALLELE1),]


 #     print(sum(is.na(Combined_Processed_Data$COORD_Uni)))


  print("Splitting to separate dfs")

      # Split dataframe by unique STUDY values
      study_dfs <- lapply(split(Combined_Processed_Data, Combined_Processed_Data$STUDY), function(sub_df) {
        study_name <- unique(sub_df$STUDY)

        # Rename columns by prefixing the study name, except for COORD_Uni
        colnames(sub_df) <- c(paste0(study_name, "_", colnames(sub_df)[-ncol(sub_df)]), "COORD_Uni")

        return(sub_df)
      })



   #   print(study_dfs)



    #  print(study_dfs)


      # # Rejoin data by COORD_Uni
      # Combined_Processed_Data_Rejoined <- Reduce(function(x, y) dplyr::full_join(x, y, by = "COORD_Uni"), study_dfs)
      #
      # # Display result
      # print(Combined_Processed_Data_Rejoined)

      print(paste0("Joining all to ", Match_Allele_Study, " ", "by unified ID" ))



      # Extract the first study's dataframe as the reference
      reference_study_name <- Match_Allele_Study  # Take the first study name
      reference_df <- study_dfs[[reference_study_name]]  # Extract first dataframe

    #  return(study_dfs)

      # Loop over remaining studies and left join them to the reference dataframe using COORD_Uni
      for (study_name in setdiff(names(study_dfs), reference_study_name)) {

        print(paste0("Joining ", study_name, " to ", reference_study_name))
        coord_col <- grep("COORD_Uni", colnames(study_dfs[[study_name]]), value = TRUE)

        reference_df <- dplyr::left_join(reference_df, study_dfs[[study_name]], by = coord_col)
      }

      # Store the final merged dataframe
      Combined_Processed_Data_Joined <- reference_df


    #  return(Combined_Processed_Data_Joined)


  #    print(nrow(Combined_Processed_Data_Joined))

  #    print(colnames(Combined_Processed_Data_Joined))
      genpos_cols <- grep("GENPOS$", colnames(Combined_Processed_Data_Joined), value = TRUE)


     print("GENPOS cols")
     print(genpos_cols)


     # Filter genpos_cols to only keep those containing Match_Allele_Study
     genpos_cols <- genpos_cols[grepl(Match_Allele_Study, genpos_cols)]


     print("Only excluding if missing value for this/backbone:")
     # Print the result
     print(genpos_cols)


     #only if not in backbone is excluded as allele flip cant be referenced to something and out method is that now

   #   return(Combined_Processed_Data_Joined)

     print("Removing SNPs with missing data from backbone/allele flip study...")

     #could allow to keep all, remove any blank in futue like with the below code originally....
     #Need to adjust in future for non-regenie X, is it 23, X, chr x etc.
     #Could add combo cols below to generate files with all combos of studies to meta.

      # Drop rows where all GENPOS columns are NA
      if (length(genpos_cols) > 0) {
        Combined_Processed_Data_Joined <- Combined_Processed_Data_Joined %>%
          dplyr::filter(!dplyr::if_any(dplyr::all_of(genpos_cols), is.na))
      }



   #  print("HO")

 #    return(Combined_Processed_Data_Joined)



      #    print(nrow(Combined_Processed_Data_Joined))

      #remove those which aren't present in all studies

      # Display the result
   #   print(Combined_Processed_Data_Joined)

      #ADD REMOVE NO BETA FLIP EVAL COLS AND ADJUST FOR WHICH STUDY TO JOIN TO AND REMPOVE!
      #check if genpos removes one and the other
      #BETA, SEs only and COORD_Uni ID of ref name only.




      # Identify columns that end exactly with "BETA" or "SE" plus "COORD_Uni"
      keep_cols <- colnames(Combined_Processed_Data_Joined)[
        grepl(paste0(Output, "$", collapse = "|"), colnames(Combined_Processed_Data_Joined)) |
          colnames(Combined_Processed_Data_Joined) == "COORD_Uni"
      ]


      print("These column types will be kept from each data frame and married to unified ID in order:")
      print(Output)



      # Preserve the original order found in the dataframe
      ordered_cols <- colnames(Combined_Processed_Data_Joined)[colnames(Combined_Processed_Data_Joined) %in% keep_cols]

      # Ensure COORD_Uni is always first while maintaining the original order for the rest
      ordered_cols <- c("COORD_Uni", setdiff(ordered_cols, "COORD_Uni"))

      print("Therefore these columns will be kept in saved output:")
      print(ordered_cols)





      print("Filtering columns")


      # Select columns in the correct order
      Combined_Processed_Data_Joined <- Combined_Processed_Data_Joined %>%
        dplyr::select(all_of(ordered_cols))



      print("Saving Output")

      filename <- paste0(File_Name, ".", File_Type)

      print(filename)  # Optional: to see the generated filename


    #  write.table(Combined_Processed_Data_Joined, filename, sep = "\t", row.names = FALSE, quote = FALSE)



      return(Combined_Processed_Data_Joined)




}

