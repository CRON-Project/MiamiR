
#INCOME AND INTELLIGENCE BOTH HG37 NEED SAME BUILD FOR SAME SNP UNLESS RSID
#' Title
#'
#' @param Data_Sets These are a list of GWAS summary statistics to be evaluated; defaults to c("LbDementia_Sum_Stats", "Intelligence_Sum_Stats")
#' @param Names These are a list of names for the data sets to be labelled with in the output file - this also specifies the order down the page; defaults to c("Dementia", "Intelligence")
#' @param Chromosome_Columns These are a list of manual chromosome column names for the data sets used in the order specified; defaults to c()
#' @param Reference_Alleles These are a list of manual reference allele column names for the data sets used in the order specified; defaults to c()
#' @param Effect_Alleles These are a list of manual effect allele column names for the data sets used in the order specified; defaults to c()
#' @param Position_Columns These are a list of manual position column names for the data sets used in the order specified; defaults to c()
#' @param SNP_ID_Columns These are a list of manual SNP ID column names for the data sets used in the order specified; defaults to c()
#' @param Beta_Columns These are a list of manual BETA column names for the data sets used in the order specified; defaults to  c()
#' @param Standard_Error_Columns  These are a list of manual SE column names for the data sets used in the order specified; defaults to c()
#' @param PValue_Columns These are a list of manual P column names for the data sets used in the order specified; defaults to c()
#' @param Peak_Separation Distnance in BP necessary from a lead SNP for another region of association to begin; defaults to 100K
#' @param File_Name File name to save data as; defaults to "Peaks"
#' @param File_Type File type of saved data; defaults to "txt"
#'
#' @return
#' @export
#'
#' @examples
#'
#'

#'Peak_Finder_Outcome <- Peak_Finder(Data_Sets = Locations,
#'                                 Peak_Separation = 2000000,
#'                                 File_Name = "Peaks",
#'                                 Spec_CHROM = 1,
#'                                 File_Type = "txt")
#'
#'
#'
#'

# setwd("C:/Users/callumon/Miami_Package_R/MiamiR/inst/extdata")

#USE THIS TO LOAD DATA STORED ELSEWHERE as binary
# load("HillWD_31844048_household_Income.txt.rda")
# load("SavageJansen_2018_intelligence_metaanalysis.rda")
# HillWD_31844048_household_Income.txt <- HillWD_31844048_household_Income.txt[HillWD_31844048_household_Income.txt$Chr == 1,]
# SavageJansen_2018_intelligence_metaanalysis <- SavageJansen_2018_intelligence_metaanalysis[SavageJansen_2018_intelligence_metaanalysis$CHR == 1,]


Peak_Finder <- function(Data_Sets = c("HillWD_31844048_household_Income.txt", "SavageJansen_2018_intelligence_metaanalysis"),
                        Names = NULL,
                        Chromosome_Columns = c(),
                        Reference_Alleles = c(),
                        Effect_Alleles = c(),
                        Position_Columns = c() , SNP_ID_Columns = c(),
                        Beta_Columns = c(), Standard_Error_Columns = c(),
                        PValue_Columns = c(),
                        Spec_CHROM = NULL,
                        Peak_Separation = 1000000,
                        File_Name = "Peaks",
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
          data.table::fread(cmd = paste0("zcat ", path, " | head -n 5"))
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
          message("📂 Detected compressed file, using `fread(cmd = 'grep')` to filter CHROM == 1.")
          temp_df <- data.table::fread(cmd = paste0("zcat ", path, " | grep '^1\\s'"), sep = "\t")
        } else {
          message("🚀 Using `fread()` to load only CHROM == 1.")
          temp_df <- vroom::vroom(path)  # Read everything first (unavoidable)

          print(temp_df)
          print(chrom_col)

          # **Apply filtering right after reading**
          if (!is.null(chrom_col)) {
           # temp_df <- temp_df[temp_df$Chr == 1,]
            if(is.null(Spec_CHROM))
            {
              print("Whole file desired")
            }else{
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
      message(paste("🔍 Preview of", dataset_name, ":"))
      print(head(df))  # Show first few rows of the dataset

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


    Data$COORD_Norm <- stringi::stri_c("chr", Data$CHROM, ":", Data$GENPOS, ":", Data$ALLELE0, ":", Data$ALLELE1)
    Data$COORD_Alt <- stringi::stri_c("chr", Data$CHROM, ":", Data$GENPOS, ":", Data$ALLELE1, ":", Data$ALLELE0)


  Combined_Processed_Data <- rbind(Combined_Processed_Data, Data)


  print("Combined")
#  print(Combined_Processed_Data)

}


#  print(sum(is.na(Combined_Processed_Data$GENPOS)))
#  print(table(Combined_Processed_Data$CHROM))
  #out of loop now

  find_min_p_with_distance <- function(data) {
    # Ensure GENPOS is numeric
    data <- data %>%
      dplyr::mutate(GENPOS = as.numeric(GENPOS))

    # Sort data by P-value (smallest first)
    data <- data[order(data$P), ]

    # Initialize selected rows
    selected_rows <- data.frame()

    # While there are rows left in the data
    while (nrow(data) > 0) {
      # Select the row with the smallest P-value
      best_row <- data[1, ]
      selected_rows <- rbind(selected_rows, best_row)

      # Filter out rows within Peak_Separation distance
      data <- data %>%
        dplyr::filter(abs(GENPOS - best_row$GENPOS) >= Peak_Separation)
    }

    return(selected_rows)
  }





  ####



  result <- Combined_Processed_Data %>%
    dplyr::filter(P < 5e-8) %>%
    dplyr::group_by(STUDY, CHROM) %>%
    dplyr::group_modify(~ find_min_p_with_distance(.x)) %>%
    dplyr::ungroup()



  print("Peak Algorithm Search")
#    return(result)



  print("Study specific peaks per CHROM")
#  print(result)



  #get matching rows from original combined df for comparison

  # Extract unique COORD_Norm and COORD_Alt from the result
  unique_coords <- unique(c(result$COORD_Norm, result$COORD_Alt))

#  print("Unique coords:")
#  print(unique_coords)



  # Filter rows from Combined_Processed_Data where COORD_Norm or COORD_Alt matches any unique_coords
  matching_rows <- Combined_Processed_Data %>%
    dplyr::filter(COORD_Norm %in% unique_coords | COORD_Alt %in% unique_coords)

  # Print the matching rows
#  print(matching_rows)



  matching_rows <- matching_rows %>%
    dplyr::group_by(CHROM, GENPOS) %>%
    dplyr::mutate(
      COORD_Uni = ifelse(
        ALLELE0 == unlist(strsplit(COORD_Norm, ":"))[3] &
          ALLELE1 == unlist(strsplit(COORD_Norm, ":"))[4],
        COORD_Norm, # Use COORD_Norm if order matches
        NA          # Set NA for others initially
      )
    ) %>%
    dplyr::ungroup()

  # Fill in COORD_Uni for duplicates based on the row that matched the order
  matching_rows <- matching_rows %>%
    dplyr::group_by(CHROM, GENPOS) %>%
    dplyr::mutate(
      COORD_Uni = ifelse(
        is.na(COORD_Uni),
        COORD_Uni[!is.na(COORD_Uni)][1], # Assign the matching row's COORD_Uni
        COORD_Uni
      )
    ) %>%
    dplyr::ungroup()

 # print(matching_rows)

#   return(matching_rows)
# zzz
  #Matching rows gets the

 # return(matching_rows)

  # Exclude rows from the same STUDY as in `result`
  #matching_other_studies <- matching_rows %>%
   # dplyr::filter(!(STUDY %in% result$STUDY))

  # Print matching rows for other studies
#  print(matching_rows)


  # print(nrow(result))
  # print(nrow(matching_rows))
  # combined_result <- rbind(result, matching_rows)
  # print(nrow(combined_result))
  # # Print the combined result
  # print(combined_result)
  #
  #
  # print("removing duplicate rows, just want ones from other studies")
  #
  #
  # print(nrow(combined_result))
  #
  # # Remove duplicate rows
  # unique_combined_result <- combined_result %>%
  #   dplyr::distinct()
  #
  # print(nrow(unique_combined_result))
  #

 # filename <- paste0(File_Name, ".", File_Type)
  #print(filename)

  #write.table(matching_rows, filename)



  # Split `matching_rows` into a list of dataframes by `STUDY`
  split_dfs <- split(matching_rows, matching_rows$STUDY)

  # Iterate through each study-specific dataframe and save them
  for (study_name in names(split_dfs)) {


  #  print(study_name)

    # Get the dataframe for the current study
    study_df <- split_dfs[[study_name]]


  #  print(study_df)

    # Drop columns `COORD_Alt`, `ID`, and `STUDY`
    study_df <- study_df %>%
      dplyr::select(-c(COORD_Alt, COORD_Norm, ID, STUDY))

    # Rename `COORD_Norm` to `ID`
    colnames(study_df)[colnames(study_df) == "COORD_Uni"] <- "ID"

    # Generate a filename with the study name
    filename <- paste0(Spec_CHROM, "_", study_name, "_", File_Name, ".", File_Type)
    print("SAVING...")
    print(filename)  # Optional: to see the generated filename

 #   setwd("C:/Users/callumon/Miami_Package_R/MiamiR/data")
    # Save the dataframe to the file
    write.table(study_df, filename, sep = "\t", row.names = FALSE, quote = FALSE)
  }






  }




