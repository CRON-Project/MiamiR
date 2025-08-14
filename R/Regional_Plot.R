

#' Title
#'
#' @inheritParams Single_Plot
#'
#' @param Data
#' @param Chromosome_Column
#' @param Position_Column
#' @param SNP_ID_Column
#' @param PValue_Column
#' @param Chromosome
#' @param MIN_GENE_GAP
#' @param Separation_Distance
#' @param Region_Window
#' @param P_threshold
#' @param Genome_Build
#' @param Population
#' @param Auto_LD_Region_Size
#' @param Auto_LD
#' @param Gene_Tracks
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
#'

Regional_Plot <- function(Data = NULL,
                          Chromosome_Column = NULL,
                          Position_Column = NULL, SNP_ID_Column = NULL,
                          PValue_Column = NULL,
                          Chromosomes =  NULL, # c(1:22, "X", "Y"),
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
                          Gene_Biotype = NULL, #c("coding", "microrna", "lncRNA", "ppseudo"),
                          Gene_Biotype_Colour = NULL, # c("red", "yellow", "brown", "green" ),
                          Gene_Name = NULL, #c("CRIP1P1", "TPRXL"),
                          Auto_LD_Region_Size = 1e6,
                          Auto_LD = FALSE,
                          Reference_Allele_Column = NULL,  Effect_Allele_Column = NULL,
                          Verbose = FALSE,
                          Interactive = FALSE, # full T/F to TRUE/FALSE critical
                          Gene_Tracks = TRUE,  .dots = list()) {





    #Need genetic data
  default_gap_fraction <- 0.3

  MIN_GENE_GAP <- if (!is.null(MIN_GENE_GAP)) {
    MIN_GENE_GAP
  } else {
    default_gap_fraction * Region_Window
  }





 #a <-  table(gene_data$gene_biotype)
 # return(gene_data)

  # Detect chromosome column
  Chromosome_Column <- detect_chromosome_column(Data)

  # Normalize "23" to "X"



  Data[[Chromosome_Column]] <- as.character(Data[[Chromosome_Column]])
#  Data[[Chromosome_Column]][Data[[Chromosome_Column]] == "23"] <- "X"

  #Don't need that here


  PValue_Column     <- detect_pvalue_column(Data, PValue_Column)
  Position_Column   <- detect_position_column(Data, Position_Column)
  SNP_ID_Column     <- detect_snp_column(Data, SNP_ID_Column)
  Ref_Allele_Column     <- detect_reference_allele_column(Data, Reference_Allele_Column)
  Alt_Allele_Column     <- detect_effect_allele_column(Data, Effect_Allele_Column)


  #Manually assign columns for ease of use
  Data$CHROM <- Data[[Chromosome_Column]]
  Data$GENPOS <- Data[[Position_Column]]
  Data$ID <- Data[[SNP_ID_Column]]
  Data$P <- Data[[PValue_Column]]

  Data$ALLELE0 <- Data[[Ref_Allele_Column]]
  Data$ALLELE1 <- Data[[Alt_Allele_Column]]

  Data$ALLELE0 <- toupper(Data$ALLELE0 )
  Data$ALLELE1 <- toupper(Data$ALLELE1 )

  # Filter to selected chromosome

  #make sure here regarless

  filtered_data <- Data  # Default is whole dataset

  # if(!is.null(Chromosome)) # may not be specified ; default is NULL
  #
  # {
  #
  # filtered_data <- Data[Data[[Chromosome_Column]] == Chromosome, , drop = FALSE]
  #
  # }


  #importnat for shiny





  filtered_data$CHROM <- as.numeric(filtered_data$CHROM)






  # if (!is.null(Chromosome) && Chromosome != "" && Chromosome != "NULL") {
  #   Chromosome <- as.character(Chromosome)
  #   Data[[Chromosome_Column]] <- as.character(Data[[Chromosome_Column]])
  #   filtered_data <- Data[Data[[Chromosome_Column]] == Chromosome, , drop = FALSE]
  # } else {
  #   filtered_data <- Data
  # }


#  chrom_list <- c(1,2,3)
#


  if(!is.null(Chromosomes))

  {

  filtered_data <- filtered_data[filtered_data$CHROM %in% Chromosomes,]


  }else{

  #  Chromosomes <- 1
   #filtered_data <- filtered_data[filtered_data$CHROM %in% c(Chromosomes) ,]

   }





  find_peaks <- function(data, Separation_Distance = 1e6, P_threshold = 5e-8) {
    # Group by chromosome (or single group if only one CHROM)
    chrom_list <- if (length(unique(data$CHROM)) == 1) list(data) else base::split(data, data$CHROM)

    results <- base::lapply(chrom_list, function(chrom_data) {
      # Filter to only SNPs below genome-wide threshold
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
          # Only keep if this SNP is at least Separation_Distance away from all previously selected
          too_close <- base::any(base::abs(snp$GENPOS - selected$GENPOS) < Separation_Distance)
          if (!too_close) {
            selected <- dplyr::bind_rows(selected, snp)  # dplyr::bind_rows used explicitly
          }
        }
      }

      selected
    })

    dplyr::bind_rows(results)
  }



  peaks <- find_peaks(filtered_data, Separation_Distance = 1e6, P_threshold = 5e-8)



  # Step 1: Create LD_Query column (chr:pos format)
  peaks <- peaks %>%
   # mutate(
    #  LD_Query = paste0("chr", CHR, ":", POS)
  #  )
    dplyr::mutate(
      LD_Query = paste0("chr", CHROM, ":", GENPOS)
    )






  # Step 2: Loop through each row and query LDproxy()
  ld_results_list <- list()






  if(Auto_LD == TRUE)
  {



  for (i in seq_len(nrow(peaks))) {
    message("Calculating LD for lead SNP from peak ", i, " ...")

    query <- peaks$LD_Query[i]


    token <- "a3a5b2b4d4c5"





    result <- LDlinkR::LDproxy(
      snp = query,      # Query SNP
      pop = Population,           # Population: EUR, AFR, EAS, AMR, SAS
      r2d = "r2",            # LD metric: "r2" or "d"
      token = token,  # Or use get_token() if you've set it
      win_size = Auto_LD_Region_Size,
      genome_build = Genome_Build
    )

    #https://www.rdocumentation.org/packages/LDlinkR/versions/1.4.0/topics/LDproxy
 #   Choose between one of the three options...`grch37` for genome build GRCh37 (hg19), `grch38` for GRCh38 (hg38), or `grch38_high_coverage` for GRCh38 High Coverage (hg38) 1000 Genome Project data sets. Default is GRCh37 (hg19).


      # Store result in the list with identifier
      ld_results_list[[i]] <- result %>%
        dplyr::mutate(Peak_Index = i, Lead_SNP = peaks$SNP[i])


  }


  }





  get_region_around_all_peaks <- function(filtered_data, peaks, Region_Window) {
    if (nrow(peaks) == 0) return(filtered_data[0, , drop = FALSE])  # return empty df if no peaks

    # Store list of regional data frames
    region_list <- lapply(seq_len(nrow(peaks)), function(i) {
      peak <- peaks[i, ]
      chrom <- peak$CHROM
      pos <- peak$GENPOS

      # Subset to same chromosome
      chr_data <- filtered_data[filtered_data$CHROM == chrom, , drop = FALSE]

      # Subset to window around peak
      region_data <- chr_data[chr_data$GENPOS >= (pos - Region_Window) & chr_data$GENPOS <= (pos + Region_Window), ]

      # Add a Peak_ID column
      region_data$Peak_ID <- paste0("peak", i)

      return(region_data)
    })

    # Combine all regions into one data frame
    combined_regions <- do.call(rbind, region_list)

    return(combined_regions)
  }






  region_df <- get_region_around_all_peaks(filtered_data, peaks, Region_Window = Region_Window)



  # Add pretty peak and final plot title
  region_df$pretty_peak_id_singleton <- gsub("peak(\\d+)", "Peak \\1", region_df$Peak_ID)
  region_df$Singleton_Name <- paste0("Chromosome ", region_df$CHROM, " ", region_df$pretty_peak_id_singleton, " Regional Plot")

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





  # Capture user arguments (e.g. Condense_Scale, Title, etc.)
  # Capture user arguments
  user_args <- .dots




  # Pull defaults from Single_Plot
  default_args <- formals(.Single_Plot_original) #changed to incorporate loading bar



  # Evaluate any symbols (like language objects)
  # evaluated_defaults <- lapply(default_args, function(x) {
  #   if (is.symbol(x) || is.language(x)) eval(x) else x
  # })
  #
  evaluated_defaults <- lapply(default_args, function(x) {
    if (missing(x)) return(NULL)  # just in case, though rarely hit
    if (is.symbol(x) || is.language(x)) {
      tryCatch(eval(x), error = function(e) NULL)
    } else {
      x
    }
  })


#


  # Merge user args over defaults (user overrides default)
  args <- modifyList(evaluated_defaults, user_args)



#print(args)
Colour_Of_Diamond <- args$Colour_Of_Diamond



  # Provide defaults only if not supplied
  if (is.null(args$Title)) {
    args$Title <- paste0("Chromosome", Chromosomes, " Regional Plot")
  }





  # Inject required args
  args$Data <- filtered_data
 # args$Point_Size <- 25







  if(Interactive == TRUE)
  {
    args$Interactive <- TRUE
  }else{
    args$Interactive <- FALSE
  }



  # Add Mb-based x-axis
  pos_range <- range(filtered_data$POS, na.rm = TRUE)
  breaks <- pretty(pos_range, n = 6)
  labels <- paste0(round(breaks / 1e6, 1), " Mb")

  # plot <- plot + ggplot2::scale_x_continuous(
  #   breaks = breaks,
  #   labels = labels,
  #   position = "bottom",
  #   expand = ggplot2::expansion(mult = c(0.01, 0.01))
  # )








  plot_all_peak_regions <- function(region_df, user_args, Chromosome = NULL) {
    if (!"Peak_ID" %in% names(region_df)) {
      stop("The dataframe must include a 'Peak_ID' column.")
    }




    region_df$pretty_peak_id <- gsub("peak(\\d+)", "Peak \\1",  region_df$Peak_ID)

    region_df$Title_Name <- paste0("Chromosome ", region_df$CHROM, " ", region_df$pretty_peak_id, " " ,"Regional Plot") # real title with spaces

    #save name with no title spacing
    region_df$Save_Name <- paste0("Chromosome ", region_df$CHROM, " ", region_df$pretty_peak_id, " " ,"Regional Plot") # real title with spaces

    Title_Names <- unique(region_df$Save_Name)





    # Get all unique peak IDs
    peak_ids <- unique(region_df$Peak_ID)

    CHROM_Lab <- unique(region_df$CHROM)




    # Pull default args from Single_Plot
    default_args <- formals(Single_Plot)


 #   evaluated_defaults <- lapply(default_args, function(x) {
#      if (is.symbol(x) || is.language(x)) eval(x) else x
#    })


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



#    user_args <- list(...)
    args_base <- modifyList(evaluated_defaults, user_args)






    # Create list of plots
    plots <- lapply(seq_along(peak_ids), function(i) {


      peak_id <- peak_ids[i]
      filtered_data <- region_df[region_df$Peak_ID == peak_id, , drop = FALSE]

      filtered_data <- filtered_data %>%
        dplyr::mutate(Coord = paste0("chr", CHR, ":", as.integer(POS))) %>%
       dplyr::mutate(
          LD_Query = paste0("chr", CHROM, ":", GENPOS)
        )




      if(Auto_LD == TRUE)

      {

      # Get corresponding LD result
      ld_df <- ld_results_list[[i]]


      }






      if(Auto_LD == TRUE)

      {
        if (!is.null(ld_df)) {



          ld_df_clean <- ld_df
          #  select(Coord, RS_Number, R2) %>%
          ld_df_clean$Proxy_SNP <- ld_df_clean$RS_Number
          ld_df_clean$R2_LD <- ld_df_clean$R2
         #   rename(Proxy_SNP = RS_Number, R2_LD = R2)




          if (!"Coord" %in% colnames(ld_df_clean)) {
            ld_df_clean$Coord <- NA_character_
            ld_df_clean$R2_LD <- NA_character_ #need na values for later to grey
        }

          #sometimes lead snp missing from ld database and coord fails return in that case add NA coord


          # Join with filtered_data
          filtered_data <- filtered_data %>%
            dplyr::left_join(ld_df_clean, by = "Coord")

        } else {
          filtered_data$R2_LD <- NA
          filtered_data$Proxy_SNP <- NA
        }


      }


      if (nrow(filtered_data) == 0) return(NULL)  # skip empty



      #save real args for later injecitons
      args_real <- args


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




      lead_snp <- filtered_data %>%
        dplyr::filter(P == min(P, na.rm = TRUE)) %>%
        dplyr::pull(ID)  # or whatever your SNP ID column is called



      if(Auto_LD == TRUE)

      {

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







      if (Interactive == TRUE) {
        # filtered_data$Hover_Info <- paste0("SNP: ", filtered_data$ID, "\nCHR: ", filtered_data$CHROM, "\nPOS: ", filtered_data$GENPOS,
        #                                    "\nP: ", signif(filtered_data$P, 4),
        #                                    "\nREF: ", filtered_data$ALLELE0, " ALT: ", filtered_data$ALLELE1,
        #                                    "\nLD: ", signif(filtered_data$R2_LD, 2))
        #

        if(Auto_LD == TRUE)

        {


        filtered_data$Hover_Info <- paste0(
          "SNP: ", filtered_data$ID, "\n",
          "CHR: ", filtered_data$CHROM, "\n",
          "POS: ", filtered_data$GENPOS, "\n",
          "P: ", signif(filtered_data$P, 4), "\n",
          "REF: ", filtered_data$ALLELE0, " ALT: ", filtered_data$ALLELE1, "\n",
          "LD: ",  signif(filtered_data$R2_LD,2)
        )


        }else{


          filtered_data$Hover_Info <- paste0(
            "SNP: ", filtered_data$ID, "\n",
            "CHR: ", filtered_data$CHROM, "\n",
            "POS: ", filtered_data$GENPOS, "\n",
            "P: ", signif(filtered_data$P, 4), "\n",
            "REF: ", filtered_data$ALLELE0, " ALT: ", filtered_data$ALLELE1)


        }

        #match single plot format

        filtered_data$Hover_Info <- NA_character_
      }


      #need diamond for each as not always lead chrom peak
      filtered_data <- filtered_data %>%
        dplyr::mutate(top = (P == min(P, na.rm = TRUE)))


    #  return(filtered_data)

      args$Data <- filtered_data



#      print(peak_id)

 #     print(CHROM_Lab)

      pretty_peak_id <- gsub("peak(\\d+)", "Peak \\1", peak_id)



  #    print(pretty_peak_id)

    #  Title_Name <- paste0("Chromosome ", CHROM_Lab, " ",pretty_peak_id, " " ,"Regional Plot", "\n\n")


    #  Title_Name <- "bah"
      Title_Name <- unique(filtered_data$Title_Name)

  #    print(Title_Name)

      if (is.null(args$Title)) {
      #  args$Title <- paste0("Chromosome ", CHROM_Lab, " ",pretty_peak_id, " " ,"Regional Plot\n\n\n") #this is the actual one
      args$Title <- paste0(Title_Name, "\n\n\n")
        }

   #   args$X_Axis_Title <- "Genomic Position (Mb)"


      # args <- modifyList(args_base, list(
      #   Data = filtered_data,
      #   Title = Title_Name,
      #   X_Axis_Title = "Genomic Position (Mb)"
      # ))

      args <- modifyList(args, list(
        Data = filtered_data
     #   Title = Title_Name,
   #     X_Axis_Title = "Genomic Position (Mb)"
      ))



    #  return(filtered_data)
      #args[names(user_args)] <- user_args



      #
      # if(Interactive == TRUE) {
      #   Data$Hover_Info <- paste0(
      #     "CHR: ", Data$CHROM, "\n",
      #     "POS: ", Data$GENPOS, "\n",
      #     "P: ", signif(Data$P, 4)
      #   )
      # }



      # Generate plot safely
   #   plot <- tryCatch({


      #will fail if draft plot included as TRUE!

        #p <- do.call(Single_Plot, args)



      args$session <- NULL

      #X2 defa

      #manually set
  #    args$Point_Size <- 10
       #Don't want alternating
       args$Chromosome_Colours <- Point_Colour

      #line below fail


      if(is.null(Recombination_Axis_Title_Size))
      {
        Recombination_Axis_Title_Size <- args$Y_Axis_Title_Size
      }
       if(is.null(Recombination_Axis_Text_Size))
       {
         Recombination_Axis_Text_Size <- args$Y_Axis_Text_Size
       }


    #   print(args)

      p <- do.call(.Single_Plot_original, args)





    #  print( args$Title )
    #  return(p)




      retrieve <- peak_id # unique(region_df$Save_Name)  # e.g., "Chromosome 22 Peak 1 Regional Plot"


      chrom <- unique(p$data$CHROM)
  #    print(retrieve)
  #   chrom <- unique(p$data$CHROM)
  #    print(Title_Names)


    #  return(retrieve)



#return(p)


#
#         if(Interactive == TRUE)
#
#         {
#           p <- p + ggplot2::aes(text = Hover_Info)
#         }

  #    print("real break determinant")

        # Add Mb-based x-axis
        pos_range <- range(filtered_data$POS, na.rm = TRUE)
        breaks <- pretty(pos_range, n = 5)
        labels <- paste0(round(breaks / 1e6, 2), " Mb")




        #auto scale mb

        # Generate initial breaks/labels
        pos_range <- range(filtered_data$POS, na.rm = TRUE)
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








    #    print(chrom)


        if(Genome_Build == "grch38")
        {
          setwd("C:/Users/callumon/Miami_Package_R/MiamiR/Processed_Gene_Data")
          read_name <- paste0("Chromosome_", chrom, "_", "Gene_Data_HG38_Processed.json")
          gene_data <- jsonlite::read_json(read_name, simplifyVector = TRUE)
          setwd("C:/Users/callumon/Miami_Package_R/MiamiR")
          recomb_data <- read.csv("Processed_AVG_Recomb_HG38.csv")

        }

        if(Genome_Build == "grch37")
        {
          setwd("C:/Users/callumon/Miami_Package_R/MiamiR/Processed_Gene_Data")
          read_name <- paste0("Chromosome_", chrom, "_", "Gene_Data_HG19_Processed.json")
          gene_data <- jsonlite::read_json(read_name, simplifyVector = TRUE)
          setwd("C:/Users/callumon/Miami_Package_R/MiamiR")
          recomb_data <- read.csv("Processed_AVG_Recomb_HG19.csv")
        }






        if(Recombination_Line == TRUE)

        {

          # 1️⃣ Subset recombination data for correct chromosome and region
          chrom <- unique(p$data$CHROM)
          genpos <- p$data$GENPOS
          xrange <- range(genpos, na.rm = TRUE)

          recomb_df <- recomb_data %>%
            filter(CHROM == paste0("chr", chrom),
                   end >= xrange[1],
                   start <= xrange[2]) %>%
            mutate(midpoint = (start + end) / 2)

          #recomb_df$score <- recomb_df$score  / 10

          # 2️⃣ Compute scale factor for nice overlay
          #    max_p <- max(p$data$log10pval, na.rm = TRUE)
          #      max_rr <- max(recomb_df$score, na.rm = TRUE)
          #      scale_factor <- max_p / max_rr


          #print(max(recomb_df$score))




          # No library() calls needed if you fully namespace all functions

          # 1️⃣ Extract chromosome and plot x-range
          chrom <- unique(p$data$CHROM)
          genpos <- p$data$GENPOS
          xrange <- range(genpos, na.rm = TRUE)

          # 2️⃣ Filter recombination data for current chromosome and plot region
          recomb_df <- dplyr::filter(
            recomb_data,
            CHROM == paste0("chr", chrom),
            end >= xrange[1],
            start <= xrange[2]
          ) %>%
            dplyr::mutate(
              midpoint = (start + end) / 2
            )

          # 3️⃣ Compute scale factor to align recombination rate with -log10(p) scale
          max_p <- max(p$data$log10pval, na.rm = TRUE)
          max_rr <- max(recomb_df$score, na.rm = TRUE)
          scale_factorz <- max_p / max_rr

          #print(paste("Scale factor:", scale_factorz))

          # 4️⃣ Add recombination rate track as overlay with secondary y-axis
          p <- p +
            ggplot2::geom_line(
              data = recomb_df,
              mapping = ggplot2::aes(x = midpoint, y = score * scale_factorz),
              color = Recombination_Line_Colour,
              size = Recombination_Line_Thickness,
              linetype = Recombination_Line_Type
            ) +
            ggplot2::scale_y_continuous(
            #  name = expression(-log[10](p)),
              name = args$Y_Axis_Title, #does inherit fine but needs to be explicit here to map if changed from default
              sec.axis = ggplot2::sec_axis(~ . / scale_factorz, name = Recombination_Axis_Title)
            ) +
            ggplot2::theme(
              axis.title.y.right = ggplot2::element_text(color = "black", vjust = 7, size = Recombination_Axis_Title_Size ),
              axis.text.y.right =  ggplot2::element_text(size = Recombination_Axis_Text_Size, colour = "black"),
            )

          # 5️⃣ Print the updated plot
          #return(p2)

          # 3️⃣ Overlay as line with secondary axis

          #
          #         p  <- p +
          #           geom_line(
          #             data = recomb_df,
          #             aes(x = midpoint, y = score * scale_factor),
          #             color = Recombination_Line_Colour,
          #             linetype = Recombination_Line_Type,   # ← line type control
          #             size = Recombination_Line_Thickness             # ← line thickness (default ≈ 1, smaller = thinner)
          #           ) +
          #           scale_y_continuous(
          #             name = "-log10(P)",
          #             sec.axis = sec_axis(~ . / scale_factor, name = "Recombination rate (cM/Mb)")
          #           )


        }






        if (!is.null(Gene_Biotype)) {

          # Alias table (as you wrote it)
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

          #  print(valid_biotypes)

          gene_data <- gene_data %>%
            dplyr::filter(gene_biotype %in% valid_biotypes)

          # Handle colouring logic
          if (!is.null(Gene_Biotype_Colour)) {

            if (length(Gene_Biotype_Colour) == 1) {
              # Single color → apply to all
              gene_data$Gene_Biotype_Colour <- Gene_Biotype_Colour

            } else if (length(Gene_Biotype_Colour) >= length(valid_biotypes)) {
              # Multiple colors → assign based on order of valid_biotypes
              color_map <- setNames(Gene_Biotype_Colour[seq_along(valid_biotypes)], valid_biotypes)
              gene_data$Gene_Biotype_Colour <- color_map[gene_data$gene_biotype]
            } else {
              # Not enough colors provided → fallback (e.g., default single color)
              warning("Not enough colors for Gene_Biotype; using first color for all.")
              gene_data$Gene_Biotype_Colour <- Gene_Biotype_Colour[1]
            }

          } else {
            # If Gene_Biotype_Colour not provided → optional: assign default color
            gene_data$Gene_Biotype_Colour <- "black"
          }

        }


        if(!is.null(Gene_Name))

        {

          gene_data <- gene_data %>%
            dplyr::filter(gene_name %in% Gene_Name)


        }




  #      return(filtered_data)
        # pos_range <- range(filtered_data$POS, na.rm = TRUE)
        # span_bp <- diff(pos_range)
        # span_mb <- span_bp / 1e6
        #
        # decimal_places <- dplyr::case_when(
        #   span_mb < 0.001 ~ 5,
        #   span_mb < 0.01 ~ 4,
        #   span_mb < 0.1 ~ 3,
        #   span_mb < 3 ~ 5,
        #   span_mb < 10 ~ 1,
        #   TRUE ~ 0
        # )
        #
        # # Define desired number of breaks (adjust based on span if needed)
        # n_breaks <- if (span_mb < 0.05) 2 else if (span_mb < 0.1) 3 else 5
        #
        # # Force evenly spaced breaks including min/max
        # break_step <- span_bp / (n_breaks - 1)
        # breaks <- seq(from = pos_range[1], to = pos_range[2], by = break_step)
        #
        # # Ensure last value exactly = pos_range[2]
        # if (tail(breaks, 1) < pos_range[2]) {
        #   breaks <- c(breaks, pos_range[2])
        # } else {
        #   breaks[length(breaks)] <- pos_range[2]
        # }
        #
        # labels <- paste0(formatC(breaks / 1e6, format = "f", digits = decimal_places), " Mb")
        #
        # print(pos_range)
       # return(filtered_data)


        #need to make sure breaks arent outside the pos range due to round down at low end or round up at  high end



      #   if (Gene_Tracks) {
      # #    p <- p + ggplot2::scale_x_continuous(
      # #      breaks = breaks,
      # #      labels = labels,
      # #      position = "bottom",
      # #      expand = ggplot2::expansion(mult = c(0.01, 0.01))
      # #    ) +
      #  p <- p   +   ggplot2::theme(axis.title.x = element_text(size = 0),
      #                      axis.text.x = element_blank(),
      #                      axis.ticks.x = element_blank(),
      #                      plot.margin = margin(t = 0, b = 0, l = 50))
      #   }

        # Filter gene_data to only genes within the region of interest

#        print(gene_data)


        # Assuming pos_range is already defined
        top <- max(pos_range, na.rm = TRUE)
        bottom <- min(pos_range, na.rm = TRUE)

 #       print(top)
  #      print(bottom)

        # Step-by-step filtering:
        gene_data_top <- gene_data[gene_data$start <= top, , drop = FALSE]
        gene_data_filtered <- gene_data_top[gene_data_top$end >= bottom, , drop = FALSE]





        gene_data <- gene_data_filtered



#        gene_data <- gene_data[gene_data$start <= max(pos_range) & gene_data$end >= min(pos_range), drop = FALSE]


   #     return(gene_data)

#
#         label_nudge <- 0.015 * diff(pos_range)
#
#
#         gene_data <- gene_data %>%
#           dplyr::mutate(
#             visible_start = pmax(start, pos_range[1]),
#             visible_end   = pmin(end,   pos_range[2]),
#             visible_span  = visible_end - visible_start,
#             view_width    = diff(pos_range),
#             frac_visible          = visible_span / (end - start),
#             frac_visible_of_view  = visible_span / view_width,
#             mid_visible = (visible_start + visible_end) / 2,
#
#             label_x = dplyr::case_when(
#               frac_visible_of_view > 0.5 ~ mid_visible,                                # Centered
#               frac_visible < 0.5 & start < pos_range[1] ~ visible_end + label_nudge,  # Near left
#               frac_visible < 0.5 & end   > pos_range[2] ~ visible_start - label_nudge,# Near right
#               TRUE ~ mid_visible
#             ),
#
#             label_hjust = dplyr::case_when(
#               frac_visible_of_view > 0.5 ~ 0.5,
#               frac_visible < 0.5 & start < pos_range[1] ~ 1,  # Right-align left-clipped
#               frac_visible < 0.5 & end   > pos_range[2] ~ 0,  # Left-align right-clipped
#               TRUE ~ 0.5
#             )
#           )






        edge_buffer <- 0.2 * diff(pos_range)

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
              # 🛑 If the gene takes up enough space, just center it
              frac_visible_of_view > 0.05 ~ mid_visible,

              # 🔸 Near left edge — regardless of visibility
              start < pos_range[1] + edge_buffer ~ visible_end + label_nudge,

              # 🔸 Near right edge — regardless of visibility
              end > pos_range[2] - edge_buffer ~ visible_start - label_nudge,

              # ✅ Mostly visible and not near edge
              frac_visible_of_view > 0.5 ~ mid_visible,

              # 🔻 Clipped left
      #        frac_visible < 0.5 & start < pos_range[1] ~ visible_end + label_nudge,

              # 🔻 Clipped right
       #       frac_visible < 0.5 & end > pos_range[2] ~ visible_start - label_nudge,

              # 🔘 Fallback
              TRUE ~ mid_visible
            ),

            label_hjust = dplyr::case_when(
              start < pos_range[1] + edge_buffer ~ 1,  # Near left edge → right align
              end > pos_range[2] - edge_buffer ~ 0,    # Near right edge → left align
              frac_visible_of_view > 0.5 ~ 0.5,
              frac_visible < 0.5 & start < pos_range[1] ~ 1,
              frac_visible < 0.5 & end > pos_range[2] ~ 0,
              TRUE ~ 0.5
            )
          )



        # Minimum gap (in base pairs) to allow genes on the same line
        min_gene_gap <- MIN_GENE_GAP  # 100 kb

        # Sort by genomic start
        gene_data <- gene_data[order(gene_data$start), ]

       # print(gene_data)

     #   gene_data$label <- dplyr::if_else(
    #      gene_data$strand == "+",
    #      paste0(gene_data$gene_name, "→"),
    #      paste0(gene_data$gene_name, "←")
    #    )

    #    gene_data$arrow_label <- ifelse(gene_data$strand == "+", "→", "←")



       # return(gene_data)



#        return(gene_data)

        if(  nrow(gene_data) != 0)

        {


          #Don't need for now
          gene_data$label <- gene_data$label |>
            gsub("^\\s*<-\\s*", "", x = _) |>              # Remove `<-` at start with any spacing
            gsub("\\s*(->|<-)\\s*$", "", x = _)            # Remove `->` or `<-` at end




        # Initialize y-tier
        gene_data$y <- NA
        tiers <- list()



        for (i in seq_len(nrow(gene_data))) {
          gene_start <- gene_data$start[i]
          gene_end <- gene_data$end[i]

          placed <- FALSE



          for (tier_idx in seq_along(tiers)) {
            last_end <- tiers[[tier_idx]]

            if (gene_start > (last_end + min_gene_gap)) {
              gene_data$y[i] <- tier_idx * .75  # 👈 increase spacing here
              tiers[[tier_idx]] <- gene_end
              placed <- TRUE
              break

            }
          }

          if (!placed) {
            tier_number <- length(tiers) + 1
            gene_data$y[i] <- tier_number * .75  # 👈 same here
            tiers[[tier_number]] <- gene_end
          }
        }

        # Assign each gene a y-axis position
    #    gene_data$y <- factor(gene_data$gene_name, levels = rev(unique(gene_data$gene_name)))



        # Sort genes by position
        gene_data <- gene_data[order((gene_data$start + gene_data$end) / 2), ]


        }else{

          Gene_Tracks <- FALSE
}


        # Assign tiered y-positions to spread overlapping genes apart
    #    gene_data$y <- rep(seq(1, length.out = nrow(gene_data)), length.out = nrow(gene_data))


  #    }, error = function(e) {
  #      warning(paste("Skipping", peak_id, "due to error:", e$message))
   #     NULL
    #  })







      #edit at end works


      x_range <- diff(pos_range)

      x_expand <- x_range * 0.05  # 5% extra on both sides




        p <- p + ggplot2::labs(x = NULL)   # <-- Add this
      #    aes(
      #      colour = LD_Bin
      #    ) +
          # Plot all points (with LD colours)


        filtered_data_No_Lead <- filtered_data[filtered_data$top == FALSE,]


       # return(filtered_data_No_Lead)

        suppressMessages(suppressWarnings({

        if(Auto_LD == TRUE & Interactive == TRUE)
        {




        p <- p +  ggplot2::geom_point(
            data = filtered_data,
            ggplot2::aes(x = GENPOS, y = -log10(P), colour = LD_Bin,text = Hover_Info),
            size = args$Point_Size # gives mapping space under cusros rsame as point
          )




        }


#^ and below only for mapping

#
#           if(Auto_LD == TRUE & Interactive == FALSE)
#           {
#
#
#
#             p <- p +  ggplot2::geom_point(
#               data = filtered_data,
#               ggplot2::aes(x = GENPOS, y = -log10(P), colour = LD_Bin),
#               size = 0
#             )
#
#
# return(p)
#
#
#
#           }





}))





       #   ggplot2::scale_x_continuous(
#
 #           expand = ggplot2::expansion(add = c(x_expand, x_expand))
  #        ) +

     #   return(p)


     #   if(1==2){

        suppressMessages(suppressWarnings({

        if(nrow(args$Data) > 1){

        if(nrow(args$Data) <= 1){

          p <- p +   ggplot2::scale_x_continuous(
            breaks = NULL,
            labels = NULL,
            position = "bottom")
            #   expand = ggplot2::expansion(mult = c(0.01, 0.01))
       #     expand = ggplot2::expansion(add = c(x_expand, x_expand))) #+

          #visual bug with x_expand when 0 as only one point.

        }else{

   p <- p +   ggplot2::scale_x_continuous(
            breaks = NULL,
            labels = NULL,
            position = "bottom" )#,
            #   expand = ggplot2::expansion(mult = c(0.01, 0.01))

            #   expand = ggplot2::expansion(add = c(x_expand, x_expand))) #+


        }

        }


        }))


          # Plot lead SNP on top (big purple diamond)
     #     ggplot2::geom_point(
    #        data = filtered_data %>% dplyr::filter(SNP == lead_snp),
    #        ggplot2::aes(x = GENPOS, y = -log10(P)),
    #        shape = 18,
    #        size = 4,
    #        colour = "purple4"
    #      )



       #

#        return(filtered_data)



        suppressMessages(suppressWarnings({

        if(Auto_LD == TRUE)
        {







#
#           # Explicit mapping:
#           p <- p + ggplot2::aes(colour = LD_Bin)
#
#           # Consistent LD colour scale:
#           p <- p +  ggplot2::scale_colour_manual(
#             #     name = "R^2~LD",
#             values = c(
#               "NA" = "darkgrey",
#               "<0.2" = "darkblue",
#               "0.2-0.4" = "lightblue",
#               "0.4-0.6" = "green",
#               "0.6-0.8" = "yellow",
#               "0.8-<1" = "red",
#               "Lead" = "transparent"
#             ),
#             na.value = "darkgrey"
#           )+
#             ggplot2::labs(colour = "R² LD\n")+
#             ggplot2::guides(
#               colour = ggplot2::guide_legend(override.aes = list(shape = 16, size = 5)),
#               shape = "none",
#               size = "none"
#             )


          # Ensure LD_Bin is a factor with correct level ordering:
          ld_levels <- c("NA", "<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-<1", "Lead")
          p$data$LD_Bin <- factor(p$data$LD_Bin, levels = ld_levels)

          # Plot color scale where Lead = "transparent" in the plot:
          plot_colors <- c(
            "NA" = "darkgrey",
            "<0.2" = "darkblue",
            "0.2-0.4" = "lightblue",
            "0.4-0.6" = "green",
            "0.6-0.8" = "yellow",
            "0.8-<1" = "red",
            "Lead" = "transparent"  # Plot Lead as invisible
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
     #     ggplot2::scale_x_continuous(
    #        breaks = breaks,
    #        labels = labels,
    #        position = "bottom",
    #        expand = ggplot2::expansion(mult = c(0.01, 0.01))
     #     ) +
       #   ggplot2::geom_point(
      #      data = filtered_data,
      #      aes(x = GENPOS, y = -log10(P), colour = R2_LD),
      #      size = 2
      #    ) +
        #  ggplot2::scale_colour_gradientn(
        #    colours = c("grey", "blue", "green", "orange", "red"),
        #    limits = c(0, 1),
        #    na.value = "grey70",
        #    name = expression(LD~(r^2))
        #  )+

      #   p <- p +
      #     ggplot2::labs(title = Title_Name)+
      #     ggplot2::theme(
      #       legend.key.height = ggplot2::unit(1, "cm"),
      #     #  plot.title = ggplot2::element_text(
      #     #    hjust = 0.5,         # Center horizontally
      #     #    vjust = -8        # Optional: make it look nice
      #     #  ),
      #
      #      # axis.title.x =    ggplot2::element_blank(),
      #     #  axis.text.x =    ggplot2::element_blank(),
      #       legend.position = "right",
      #       legend.box.just = "left",
      #       legend.spacing.y = ggplot2::unit(40, "cm"),
      #            legend.title = ggplot2::element_text(
      #             margin = ggplot2::margin(b = 30),
      #   size = LD_Legend_Title_Size,
      #   face = "bold",
      #   family = "sans" #FiraCode
      # ),
      # legend.text = ggplot2::element_text(
      #   size = LD_Legend_Text_Size,
      #   family = "sans",
      #   margin = ggplot2::margin(l = 8, b = 2)# match font for entries too# match font for entries too
      # ),
      # legend.key = ggplot2::element_blank() # ensures square color blocks]]

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
          legend_pos <- "right"  # restores ggplot2 default behavior outside plot
          legend_just <- NULL    # no justification needed for "right"

        }


        #need separate logic or early trigger
        if (LD_Legend_On == FALSE) {
          legend_pos <- "none"  # restores ggplot2 default behavior outside plot
          legend_just <- NULL    # no justification needed for "right"

        }


        if (Auto_LD == TRUE) {

        p <- p +
  #        ggplot2::labs(title = Title_Name) +
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

            #diabled originally
       #     axis.ticks.length.x  = ggplot2::unit(0.3, "cm"),
           # axis.ticks.x = element_blank(),
  #    axis.line.x = element_blank(),
    #   plot.title = ggplot2::element_blank(),
         #   plot.margin =    ggplot2::margin(t = 0, b = 0, r = 10, l = 50)
       #     panel.background = element_rect(fill = "white", colour = NA),  # keep panel white
          #  plot.background = element_rect(fill = "grey")  # shade around plot (including below axis)
      #     axis.line.x = ggplot2::element_line(color = "black")

      #    )

      # return(p)





        suppressMessages(suppressWarnings({


        if(Gene_Tracks == TRUE & Interactive == FALSE)
        {
          p <- p + ggplot2::theme(axis.line.x = ggplot2::element_blank(), plot.margin =  ggplot2::margin(t = 0, b = 0, r = 10, l = 50),
                         axis.ticks.x = ggplot2::element_blank(),
                         axis.title.x =    ggplot2::element_blank(),
                         axis.text.x =    ggplot2::element_blank())
        }







        if(Gene_Tracks == TRUE & Interactive == TRUE)
        {
          p <- p + ggplot2::theme(plot.margin =  ggplot2::margin(t = 0, b = 0, r = 10, l = 50),
                                  axis.ticks.x = ggplot2::element_blank(),
                                  axis.title.x =    ggplot2::element_blank(),
                                  axis.text.x =    ggplot2::element_blank())
        }



}))



        suppressMessages(suppressWarnings({



        if(Gene_Tracks == FALSE & (nrow(args$Data) > 1))
        {


      #    return(p)

      #    if(nrow(args$Data) <= 1){


   #         p <- p + ggplot2::scale_x_continuous(
    #          breaks = breaks,
     #         labels = labels,
      #        position = "bottom",
         #     #   expand = ggplot2::expansion(mult = c(0.01, 0.01))
       #       expand = ggplot2::expansion(add = c(x_expand, x_expand))
        #    )

        #    return(p)


       #   }else{

            p <- p + ggplot2::scale_x_continuous(
              breaks = breaks,
              labels = labels,
              position = "bottom" ,
              #   expand = ggplot2::expansion(mult = c(0.01, 0.01))
            #  expand = ggplot2::expansion(add = c(x_expand, x_expand))
            expand = c(0,0)
            )



        }





        if(Gene_Tracks == FALSE )
        {

        p <- p +
          ggplot2::labs(x = args_real$X_Axis_Title)+ ggplot2::theme(
          plot.margin =    ggplot2::margin(t = 0, b = 50, r = 10, l = 50),
          axis.title.x = ggplot2::element_text(size = args_real$X_Axis_Title_Size, vjust = -7, colour = "black"),
          axis.text.x  = ggplot2::element_text(size = args_real$Chromosome_Label_Size, vjust = -2, colour = "black" )
          )


        }




        # if (nrow(Data) == 1) {
        #   x_val <- Data[[Position_Column]][1]
        #
        #   p <- p +
        #     scale_x_continuous(
        #       limits = c(x_val - 1, x_val + 1),    # Give axis a range
        #       expand = expansion(add = c(0, 0))
        #     )
        # }
        #
        # p <- p +
        #   labs(x = "Your X-Axis Title") +
        #   theme(
        #   axis.title.x = element_text(),
        #     # Avoid pushing too low
        #     plot.margin = margin(t = 10, r = 10, b = 60, l = 10)  # Ensure enough space for it to show
        #   )
        #
        # p <- p + theme(
        #   plot.background = element_rect(color = "red", fill = NA, linewidth = 1)
        # )
        #



     #   return(p)






      if(Auto_LD == FALSE)
      {
        p <- p + ggplot2::theme(legend.position = "none")
      }



        }))





        p_main_plot <- p






      # p <- p +
      #   annotation_custom(
      #     grob = grid::rectGrob(
      #       x = grid::unit(0.5, "npc"),
      #       y = grid::unit(-0.2, "npc"),  # below the plot area
      #       width = grid::unit(1, "npc"),
      #       height = grid::unit(0.12, "npc"),
      #       gp = grid::gpar(fill = "grey90", col = NA)
      #     )
      #   ) +
      #   ggplot2::coord_cartesian(clip = "off")  # allow drawing outside plot area


        if( nrow(gene_data) != 0){

      exons_df2 <-  suppressMessages( gene_data %>%
        dplyr::select(transcript_id, label, y, exons) %>%
        tidyr::unnest(exons)
      )

      # Unnest introns (with name repair to fix duplicate column name issue)
      introns_df2 <-  suppressMessages( gene_data %>%
        dplyr::select(transcript_id, label, y, introns) %>%
        tidyr::unnest(introns, names_repair = "unique"))  #%>%
#        rename(start = intron_start, end = intron_end)


        introns_df2$start <-  introns_df2$intron_start
        introns_df2$end <-  introns_df2$intron_end

      gene_data$gene_biotype <- factor(gene_data$gene_biotype)

    #  return(gene_data)

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


        #FINE HERE


     #   return(gene_data)


        suppressMessages(suppressWarnings({



        if( nrow(gene_data) != 0)
        {

      p_genes <- ggplot2::ggplot(gene_data) +
        # ggplot2::geom_segment(
        #   ggplot2::aes(x = start, xend = end, y = y, yend = y),
        #   size = 4,
        #   color = "black",
        #   alpha = 0.4  # adjust alpha for transparency
        # )+

        ggplot2::geom_segment(data = exons_df2,
                     ggplot2::aes(x = start, xend = end, y = y, yend = y),
                     size = Gene_Structure_Size_Exons,
                     color = Gene_Structure_Colour_Exons) # +


      n_gene_rows <- max(gene_data$y)




      # Dynamic scaling of total height
      total_height <- (5 + (n_gene_rows * 2)) * (Gene_Structure_Size_Exons/7) * (Gene_Label_Size/5)  # Base height + per-row bonus


      #FINE

  #    return(p)

      #changing to z works conflicting with other scale factro

      #this scale factor below breaks it
      scale_factor <- 0.6 * (Gene_Structure_Size_Exons / 7) + 0.4 * (Gene_Label_Size / 5)

   #   return(p)

      total_height <- (8 + (n_gene_rows * 2)) * scale_factor


  #    return(p)

      gene_data <- gene_data %>%
        dplyr::mutate(
          arrow_char = ifelse(strand == "+", "→", "←"),
          arrow_x = ifelse(strand == "+", end, start),
          arrow_y = (y + 0.013) * (scale_factor/1.8285) # 👈 Small upward nudge to center better
        )




      # p_genes <- p_genes + ggplot2::geom_text(
      #   data = gene_data,
      #   ggplot2::aes(x = arrow_x, y = y, label = arrow_char),
      #   inherit.aes = FALSE,
      #   size = (Gene_Structure_Size_Introns * 5) ,
      #   family = "sans",
      #   vjust = 0.5,
      #   hjust = ifelse(gene_data$strand == "+", 0, 1),  # right-align for →, left-align for ←
      #   color = "black"
      # )







  # print(gene_data)


      if(nrow(introns_df2) > 0)
      {




p_genes <-  p_genes +  ggplot2::geom_segment(data = introns_df2,
                      ggplot2::aes(x = start, xend = end, y = y, yend = y),
                      color = Gene_Structure_Colour_Introns,
                      alpha = 0.1,
                      size = Gene_Structure_Size_Introns)



if(Sense_Arrow == TRUE)
{


# Compute view width for arrow segment and gap
view_width <- diff(pos_range)
arrow_body_length <- Sense_Arrow_Body_Length * view_width     # Total arrow line length
arrow_gap         <- Sense_Arrow_Gene_Gap * view_width    # Gap from gene edge

# Create arrow segment with small starting gap
strand_highlight_df <- gene_data %>%
  dplyr::mutate(
    xstart = dplyr::case_when(
      strand == "+" ~ end + arrow_gap,
      strand == "-" ~ start - arrow_gap
    ),
    xend = dplyr::case_when(
      strand == "+" ~ end + arrow_gap + arrow_body_length,
      strand == "-" ~ start - arrow_gap - arrow_body_length
    )
  )



threshold <- 10  # small buffer near plot edges

plot_xmin <- pos_range[1]
plot_xmax <- pos_range[2]
threshold <- 0.01 * view_width  # e.g., 1% of view width as a threshold near the edges

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
  Sense_Arrow_Body_Thickness <- Gene_Structure_Size_Introns * 0.6
}





# Draw arrowed segments with visible shaft and head
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


#return(p_genes)

      }











      suppressMessages(suppressWarnings({


        if(!is.null(Gene_Biotype_Colour))
        {


          # Ensure factor levels exactly match valid_biotypes in order
     #     gene_data$gene_biotype <- factor(
    #        gene_data$gene_biotype,
     #       levels = valid_biotypes
      #    )


      #    print(unique(valid_biotypes))
       #   print(gene_data$gene_biotype)



          valid_biotypes_pretty <- gsub("_", " ", valid_biotypes)

          gene_data$gene_biotype <- factor(
            gene_data$gene_biotype,
            levels = valid_biotypes_pretty
          )

          full_color_map <- setNames(Gene_Biotype_Colour, valid_biotypes_pretty)

          present_biotypes <- levels(droplevels(gene_data$gene_biotype))

          color_map_present <- full_color_map[present_biotypes]


     #     print(biotype_color_map)
      #    print(present_biotypes)

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


      # p_genes <-  p_genes +  ggplot2::geom_text(
      #     ggplot2::aes(
      #  #     x = (start + end) / 2,
      #       #x= label_x,
      #
      #       x = label_x, hjust = label_hjust,
      #       y = y + 0.375,
      #       label = label #,  # ⬅️ color by biotype
      #
      #     ),
      #     color = Gene_Biotype_Colour,
      #     show.legend = FALSE,
      #     angle = 0,
      #     family = "sans",
      #     fontface = "plain",
      #     vjust = 0.5,
      #     hjust = 0.5,
      #     size = Gene_Label_Size
      #   ) +
      #   ggplot2::coord_cartesian(clip = "off")+
      #   ggplot2::guides(
      #     color = ggplot2::guide_legend(
      #       override.aes = list(
      #         shape = 15,   # square
      #         size = 6,     # good size
      #         stroke = 0,   # clean edge
      #         vjust = 0.5   # force vertical align
      #       )
      #     )
      #   ) +
      #   ggplot2::geom_point(
      #     data = unique(gene_data[c("gene_biotype", "start", "y")]),
      #     ggplot2::aes(x = start, y = y, color = Gene_Biotype_Colour),
      #     shape = 15, size = 0, alpha = 0, show.legend = TRUE, inherit.aes = FALSE)
        # )+
        # ggplot2::guides(
        #   color = ggplot2::guide_legend(
        #     override.aes = list(
        #       shape = 15,
        #       size = 6,
        #       stroke = 0,
        #       vjust = 0.5,
        #       alpha = 1     # 🔥 Critical: ensure legend symbol is visible even if data points alpha=0
        #     )
        #   )
        # )
     #   geom_point(
    #      data = unique(gene_data[c("gene_biotype", "start", "y")]),
    #      aes(x = start, y = y, color = gene_biotype),
    #      shape = 15, size = 0, show.legend = TRUE, inherit.aes = FALSE
    #    )+
        # guides(
        #   color = guide_legend(
        #     override.aes = list(
        #       shape = 15,         # Square blocks
        #       size = 6
        #     )
        #   )
        # )+


      #below breaks





        }else{

          p_genes <-  p_genes +  ggplot2::geom_text(
            ggplot2::aes(
              #     x = (start + end) / 2,
              #x= label_x,

              x = label_x, hjust = label_hjust,
              y = y + 0.375,
              label = label ,  # ⬅️ color by biotype
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
            palette = "Paired",
            name = Gene_Legend_Title,
            labels = levels(gene_data$gene_biotype)
            ) +
            ggplot2::coord_cartesian(clip = "off")+
            ggplot2::guides(
              color = ggplot2::guide_legend(
                override.aes = list(
                  shape = 15,   # square
                  size = 6,     # good size
                  stroke = 0,   # clean edge
                  vjust = 0.5   # force vertical align
                )
              )
            ) +
            ggplot2::geom_point(
              data = unique(gene_data[c("gene_biotype", "start", "y")]),
              ggplot2::aes(x = start, y = y, color = gene_biotype),
              shape = 15, size = 0, alpha = 0, show.legend = TRUE, inherit.aes = FALSE
            )+
            ggplot2::guides(
              color = ggplot2::guide_legend(
                override.aes = list(
                  shape = 15,
                  size = 6,
                  stroke = 0,
                  vjust = 0.5,
                  alpha = 1     # 🔥 Critical: ensure legend symbol is visible even if data points alpha=0
                )
              )
            )





          #   geom_point(
          #      data = unique(gene_data[c("gene_biotype", "start", "y")]),
          #      aes(x = start, y = y, color = gene_biotype),
          #      shape = 15, size = 0, show.legend = TRUE, inherit.aes = FALSE
          #    )+
          # guides(
          #   color = guide_legend(
          #     override.aes = list(
          #       shape = 15,         # Square blocks
          #       size = 6
          #     )
          #   )
          # )+


          #below breaks

        }

}))






   #   return(p_genes)


      # p_genes <- p_genes + ggplot2::geom_text(
      #   data = gene_data,
      #   ggplot2::aes(
      #     x = label_x,
      #     y = y + 0.275,  # slightly below the gene name
      #     label = arrow_label
      #   ),
      #   size = Gene_Label_Size * 0.8,
      #   family = "sans",
      #   vjust = 1,       # top-aligned
      #   hjust = 0.5,     # centered horizontally
      #   color = "black",
      #   show.legend = FALSE
      # )






      p_genes <- p_genes +

        ggplot2::scale_x_continuous(
      #    limits = pos_range,
          breaks = breaks,
          labels = labels,
    #      expand = ggplot2::expansion(add = c(x_expand, x_expand)),
          sec.axis = ggplot2::dup_axis(
            labels = labels,
            name = args_real$X_Axis_Title  # Or NULL if you don’t want a title at top
          )
        )+
        ggplot2::coord_cartesian(xlim = pos_range)

      #exact clip rather than expand


      #fine here



      #above ^ sec
     # return(p_genes)


      p_genes <- p_genes +

    #    ggplot2::scale_y_continuous(limits = c(NULL, max(gene_data$y, na.rm = TRUE) + 1),
     #                               expand = c(0, 0)  )+
      #  ggplot2::theme_minimal(base_size = 10) +
        #   labs(x = NULL, y = NULL) +
        ggplot2::labs(x = args_real$X_Axis_Title, y = NULL) +
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = args_real$X_Axis_Title_Size, vjust = -7, colour = "black"),
          axis.text.x  = ggplot2::element_text(size = args_real$Chromosome_Label_Size, vjust = -2, colour = "black" ),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.title.x.top = ggplot2::element_blank(),   # ⬅️ no title at top
          axis.text.x.top = ggplot2::element_blank(),   # ⬅️ no title at top
      #    axis.ticks.x.top =  ggplot2::element_blank(),
      axis.ticks.length.x  = ggplot2::unit(0.3, "cm"),
      axis.ticks.length.x.top = grid::unit(-0.3, "cm"),
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(),
          plot.title = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(t = 0, r = 10, b = 60, l = 50),  # adds space below
          plot.title.position = "plot",
      axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.5),
      legend.justification = "left",                          # move left
      legend.title = ggplot2::element_text(
        margin = ggplot2::margin(b = 30), #t to control space between legends.
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

  #    p_genes <- p_genes + ggplot2::coord_cartesian(clip = "off")


   #   return(p_genes)


      if(Interactive == FALSE)

      {

      p_genes <- p_genes + ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(add = c(0.5, 0.5))  # bottom, top 0.5, 0.5
        )

      }else{


        y_min <- 1.25 #smaller for bigger gap
        y_max <- max(gene_data$y, na.rm = TRUE)
        pad <- 1

        y_lims <- c(y_min - pad, y_max + pad)


        p_genes <- p_genes + ggplot2::scale_y_continuous(
          limits = y_lims,
          expand = c(0, 0)  # turn off auto-padding
        )


        p_genes <- p_genes + ggplot2::theme(legend.position = "none")

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
        legend_pos_gene <- "right"  # restores ggplot2 default behavior outside plot
        legend_just_gene <- NULL    # no justification needed for "right"

      }


      #need separate logic or early trigger
      if (Gene_Legend_On == FALSE) {
        legend_pos_gene <- "none"  # restores ggplot2 default behavior outside plot
        legend_just_gene <- NULL    # no justification needed for "right"

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
            linewidth = 0.5,    # Thickness of border   # Legend background fill (optional)
          )
        )






        }




        }))



        if (Auto_LD == TRUE) {


        if(LD_Legend_Location == "Top Right")

        {
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

    #    return(p)
#ok
#
#       if(Auto_LD == FALSE)
#       {
#         p_genes <- p_genes + theme(legend.title = ggplot2::element_text(
#           margin = ggplot2::margin(b = 30, t = 350), #t to control space between legends.
#           size = 30,
#           face = "bold",
#           family = "FiraCode"
#         ))
#       }





  #  plot <- p / p_genes   + patchwork::plot_layout(heights = c(3, 2.5))


      if(Auto_LD == TRUE & Gene_Tracks == TRUE & nrow(gene_data) != 0)
      {
      #make sure legends align
    plot <- p / p_genes + patchwork::plot_layout(heights = c(3, 3), guides = "collect")



    plot <- cowplot::plot_grid(
      p, p_genes,
      ncol = 1,
      rel_heights = c(3, 3),
      align = "v",          # vertically align plots
      axis = "lr"           # align left and right axes
    )



      }






    if(Auto_LD == FALSE & Gene_Tracks == TRUE)
    {
   p_genes <-  p_genes + ggplot2::theme(legend.title = ggplot2::element_text(
        margin = ggplot2::margin(b = 10, t = 0), #t to control space between legends.
        size = 30,
        face = "bold",
        family = "sans"
      ))
      plot <- p / p_genes   + patchwork::plot_layout(heights = c(3, 3))

      plot <- cowplot::plot_grid(
        p, p_genes,
        ncol = 1,
        rel_heights = c(4.5, 3),
        align = "v",          # vertically align plots
        axis = "lr"           # align left and right axes
      )




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
    }



      if(Auto_LD == TRUE & Gene_Tracks == FALSE)
      {

        plot <- p  # / p_genes   + patchwork::plot_layout(heights = c(3, 3))
      }


if(Auto_LD == FALSE & Gene_Tracks == FALSE)
{

  plot <- p  # / p_genes   + patchwork::plot_layout(heights = c(3, 3))
}



       # plot <- patchwork::as.ggplot(plot)


    #    return(plot)

###

      if(Interactive == TRUE)

      {



       tick_y_base <- max(gene_data$y, na.rm = TRUE)
       tick_df <- data.frame(
         x = breaks,
         y = tick_y_base + 1, #.5 + current 5 expand
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



   #   p_genes <- p_genes + ggplot2::coord_cartesian(clip = "off")



      }




    if(Gene_Tracks == TRUE)
    {

    n_gene_rows <- max(gene_data$y)


    # Dynamic scaling of total height
    total_height <- (5 + (n_gene_rows * 2)) * (Gene_Structure_Size_Exons/7) * (Gene_Label_Size/5)  # Base height + per-row bonus


    scale_factor <- 0.6 * (Gene_Structure_Size_Exons / 7) + 0.4 * (Gene_Label_Size / 5)
    total_height <- (6 + (n_gene_rows * 2)) * scale_factor

#5


    attr(plot, "dynamic_height") <- total_height

    }









      if (Interactive == TRUE) {

        # Assume you already made p and p_genes separately (static ggplots)
#
#         fig1 <- ggplotly(p, tooltip = "text") %>%
#           layout(margin = list(t = 50, b = 50, l = 50, r = 50)) %>%
#           config(displayModeBar = TRUE, responsive = TRUE)
#
#         fig2 <- ggplotly(p_genes, tooltip = "text") %>%
#           layout(margin = list(t = 100, b = 50, l = 50, r = 50)) %>%
#           config(displayModeBar = TRUE, responsive = TRUE)
#
#
#
#         interactive_full <- subplot(
#           fig1, fig2,
#           nrows = 2,
#           shareY = TRUE,
#           shareX = FALSE,
#           titleY = TRUE,
#           margin = 0,                # <<<<<< KEY CHANGE
#           heights = c(0.5, 0.5)       # <<<<<< shrink gene panel
#         ) %>%
#           layout(
#             hovermode = "closest",
#             margin = list(t = 50, b = 50, l = 50, r = 50)
#           ) %>%
#           config(displayModeBar = TRUE, responsive = TRUE)
#
#         attr(interactive_full, "dynamic_height") <- total_height
#         attr(interactive_full, "interactive_panel") <- interactive_full
#
#         plot <- interactive_full


        # Always do this:
        attr(p, "main_ggplot") <- p
        attr(p, "gene_track_plot") <- p_genes
        attr(p, "dynamic_height") <- total_height
        plot <- p

        # Do NOT create any plotly or subplot logic here
     #   return(plot)


      }



        # if (Interactive == TRUE) {
        #   attr(p_genes, "plot_width") <- 60  # or dynamically calculate based on region width or number of genes
        #   attr(p, "dynamic_height") <- total_height
        #   attr(p, "interactive_panel") <- p
        #   attr(p, "gene_track_panel") <- p_genes
        #
        #   return(p)
        # }
    #  plot <- p





        return(invisible(plot))
    })

    names(plots) <- Title_Names

    for (n in names(plots)) {
      if (!is.null(singleton_positions[[n]]) && !is.na(singleton_positions[[n]])) {
        xpos <- singleton_positions[[n]]
        xpos_mb <- xpos / 1e6
        plots[[n]] <- plots[[n]] +
          ggplot2::scale_x_continuous(
            limits = c(xpos, xpos),
            breaks = xpos,
            labels = paste0(formatC(xpos_mb, format = "f", digits = 2), " Mb")
          )
      }
    }




    return(invisible(plots))
  }





  region_plots <- plot_all_peak_regions(region_df, user_args = user_args)






  return(invisible(region_plots))


  #return(invisible(region_plots))
  #invis to get rid of tables

}
# Save original

.Regional_Plot_original <- Regional_Plot  # Save original

Regional_Plot <- function(..., session = NULL) {
  dots <- list(...)

  # Get formals from original Regional_Plot
  fn_formals <- formals(.Regional_Plot_original)
  valid_args <- names(fn_formals)

  # Extract args that match Regional_Plot formal parameters
  clean_args <- dots[names(dots) %in% valid_args]

  # Fill in any missing defaults
  for (param in setdiff(valid_args, names(clean_args))) {
    default_val <- fn_formals[[param]]
    if (!identical(default_val, quote(expr = ))) {
      clean_args[[param]] <- eval(default_val)
    }
  }

  # Detect regional-only args (i.e., not inherited from Single_Plot)
  single_formals <- names(formals(.Single_Plot_original))
  regional_only <- setdiff(valid_args, c(single_formals, "...", "session"))

  # Evaluate regional-only args from the calling environment
  caller_env <- parent.frame()
  regional_args <- setNames(
    lapply(regional_only, function(arg) {
      if (exists(arg, envir = caller_env, inherits = TRUE)) {
        get(arg, envir = caller_env, inherits = TRUE)
      } else {
        NULL  # not supplied, let defaults handle it later
      }
    }),
    regional_only
  )


  # Inject regional-only args into .dots (user-specified args override)
  clean_args$.dots <- modifyList(dots, regional_args)

  # ✅ Improved missing required args detection:
  dots_names <- names(clean_args$.dots)
  if (is.null(dots_names)) dots_names <- character(0)

  all_args_names <- union(names(clean_args), dots_names)

  required_args <- names(fn_formals)[sapply(fn_formals, function(x) identical(x, quote(expr = )))]
  missing_required <- setdiff(required_args, all_args_names)

  if (length(missing_required) > 0) {
    stop("❌ Missing required argument(s): ", paste(missing_required, collapse = ", "))
  }



  # ✅ Execute
  if (identical(clean_args$Verbose, FALSE)) {
    run_with_counter(.Regional_Plot_original, args = clean_args, session = session)
  } else {
    do.call(.Regional_Plot_original, clean_args)
  }
}


# Save the original
# Save the original function
# Save original function


#CHROM <- Regional_Plot(Data = Intelligence_Sum_Stats, Chromosome = 3)
#CHROM <- Single_Plot(Data = Intelligence_Sum_Stats)


#Regional_Plot(Data = Intelligence_Peaks, Condense_Scale = F)
