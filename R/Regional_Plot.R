

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
                          Chromosome = NULL,
                          MIN_GENE_GAP = 150000,
                          Separation_Distance = 1e6,
                          Region_Window = 1e4,
                          P_threshold = 5e-8,
                          Genome_Build = "grch37",
                          Population = "EUR",
                          Auto_LD_Region_Size = 100,
                          Auto_LD = TRUE,
                          Reference_Allele_Column = NULL,  Effect_Allele_Column = NULL,
                          Interactive = F,
                          Gene_Tracks = TRUE, ...) {




  #Need genetic data


  setwd("C:/Users/callumon/Miami_Package_R/MiamiR")

  gene_data <- jsonlite::read_json("Gene_Data_HG38_Processed.json", simplifyVector = TRUE)




  # Detect chromosome column
  Chromosome_Column <- detect_chromosome_column(Data)

  # Normalize "23" to "X"


  print(Chromosome_Column)
  print("Hi")

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

  print(filtered_data$CHROM)

  print(sum(is.na(filtered_data$CHROM)))




  filtered_data$CHROM <- as.numeric(filtered_data$CHROM)


  message(print(sum(is.na(filtered_data$CHROM))))

  # Print rows where CHROM is NA
 # na_rows <- filtered_data[is.na(filtered_data$CHROM), ]
#  print(na_rows)


  print(filtered_data$CHROM)

  # if (!is.null(Chromosome) && Chromosome != "" && Chromosome != "NULL") {
  #   Chromosome <- as.character(Chromosome)
  #   Data[[Chromosome_Column]] <- as.character(Data[[Chromosome_Column]])
  #   filtered_data <- Data[Data[[Chromosome_Column]] == Chromosome, , drop = FALSE]
  # } else {
  #   filtered_data <- Data
  # }

  filtered_data <- filtered_data[filtered_data$CHROM == 3,]

  print(filtered_data$CHROM)





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



  print(peaks)

  print(filtered_data)


  print(peaks$CHROM)
  print(peaks$GENPOS)

  # Step 1: Create LD_Query column (chr:pos format)
  peaks <- peaks %>%
   # mutate(
    #  LD_Query = paste0("chr", CHR, ":", POS)
  #  )
    dplyr::mutate(
      LD_Query = paste0("chr", CHROM, ":", GENPOS)
    )



  print("a")


  #query <- peaks$LD_Query[1]

#print(query)



  # Step 2: Loop through each row and query LDproxy()
  ld_results_list <- list()


  Auto_LD <- TRUE
  print(Auto_LD)

  if(Auto_LD == TRUE)
  {

    print(nrow(peaks))

  for (i in seq_len(nrow(peaks))) {
    message("Calculating LD for lead SNP from peak ", i, " ...")

    query <- peaks$LD_Query[i]


    token <- "a3a5b2b4d4c5"


    print(query)
    print(Population)
    print(token)
    print(Auto_LD_Region_Size)
    print(Genome_Build)



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
        mutate(Peak_Index = i, Lead_SNP = peaks$SNP[i])


  }


  }


  print("Hi Again")

 # return(result)

  get_region_around_all_peaks <- function(filtered_data, peaks, Region_Window = 1e6) {
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



  region_df <- get_region_around_all_peaks(filtered_data, peaks, Region_Window = 1e6)





  print(region_df$Peak_ID)


#  region_df <- region_df[region_df$Peak_ID == "peak1",]


  filtered_data <- region_df


  print(filtered_data)


  print("here2")

  # Capture user arguments (e.g. Condense_Scale, Title, etc.)
  # Capture user arguments
  user_args <- list(...)  # all user-provided inputs

  print("here3")

  # Remove the '...' key if it's there
  user_args$... <- NULL

  print("here")

  # Pull defaults from Single_Plot
  default_args <- formals(Single_Plot)
  print("hereaoh")

  # Evaluate any symbols (like language objects)
  evaluated_defaults <- lapply(default_args, function(x) {
    if (is.symbol(x) || is.language(x)) eval(x) else x
  })


  print("hereaoh1.2")

  # Merge user args over defaults (user overrides default)
  args <- modifyList(evaluated_defaults, user_args)

  # Provide defaults only if not supplied
  if (is.null(args$Title)) {
    args$Title <- paste0("Chromosome", Chromosome, " Regional Plot")
  }

    args$X_Axis_Title <- "Genomic Position (Mb)"

    print("hereaoh2")

  # Inject required args
  args$Data <- filtered_data
  args$Point_Size <- 5

  print("hereaoh3")

  Interactive <- TRUE

  print(Interactive)



  if(Interactive == TRUE)
  {
    args$Interactive <- TRUE
  }else{
    args$Interactive <- FALSE
  }

  print("Hi2")

  #return(filtered_data)

  # # Generate the plot
  # plot <- do.call(Single_Plot, args)
  #

  print("Hi")

  # Add Mb-based x-axis
  pos_range <- range(filtered_data$POS, na.rm = TRUE)
  breaks <- pretty(pos_range, n = 5)
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

    region_df$Title_Name <- paste0("Chromosome ", region_df$CHROM, " ", region_df$pretty_peak_id, " " ,"Regional Plot")

    Title_Names <- unique(region_df$Title_Name)


   # print(Title_Names)



    # Get all unique peak IDs
    peak_ids <- unique(region_df$Peak_ID)

    CHROM_Lab <- unique(region_df$CHROM)

    # Pull default args from Single_Plot
    default_args <- formals(Single_Plot)
    evaluated_defaults <- lapply(default_args, function(x) {
      if (is.symbol(x) || is.language(x)) eval(x) else x
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


      print(ld_df)


      if(Auto_LD == TRUE)

      {
        if (!is.null(ld_df)) {

          print(head(ld_df, 10))  # 👈 or use View(ld_df) if in RStudio

          ld_df_clean <- ld_df
          #  select(Coord, RS_Number, R2) %>%
          ld_df_clean$Proxy_SNP <- ld_df_clean$RS_Number
          ld_df_clean$R2_LD <- ld_df_clean$R2
         #   rename(Proxy_SNP = RS_Number, R2_LD = R2)

          print(filtered_data)
          print(ld_df_clean)

          # Join with filtered_data
          filtered_data <- filtered_data %>%
            dplyr::left_join(ld_df_clean, by = "Coord")

        } else {
          filtered_data$R2_LD <- NA
          filtered_data$Proxy_SNP <- NA
        }


      }

      print("hi4")

      if (nrow(filtered_data) == 0) return(NULL)  # skip empty


      print("hi5")


      # Merge args and inject filtered data
      args <- args_base

      print("hi6")


      if(Auto_LD == TRUE)

      {

      filtered_data$R2_LD <- as.numeric(filtered_data$R2_LD)

      }

      lead_snp <- filtered_data %>%
        filter(P == min(P, na.rm = TRUE)) %>%
        pull(ID)  # or whatever your SNP ID column is called



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


      print("hi4")


      if (Interactive == TRUE) {
        # filtered_data$Hover_Info <- paste0("SNP: ", filtered_data$ID, "\nCHR: ", filtered_data$CHROM, "\nPOS: ", filtered_data$GENPOS,
        #                                    "\nP: ", signif(filtered_data$P, 4),
        #                                    "\nREF: ", filtered_data$ALLELE0, " ALT: ", filtered_data$ALLELE1,
        #                                    "\nLD: ", signif(filtered_data$R2_LD, 2))
        #

        filtered_data$Hover_Info <- paste0(
          "SNP: ", filtered_data$ID, "\n",
          "CHR: ", filtered_data$CHROM, "\n",
          "POS: ", filtered_data$GENPOS, "\n",
          "P: ", signif(filtered_data$P, 4), "\n",
          "REF: ", filtered_data$ALLELE0, " ALT: ", filtered_data$ALLELE1, "\n",
          "LD: ",  signif(filtered_data$R2_LD,2)
        )

        #match single plot format

      } else if (!Interactive) {
        filtered_data$Hover_Info <- NA_character_
      }

      args$Data <- filtered_data

      pretty_peak_id <- gsub("peak(\\d+)", "Peak \\1", peak_id)

      Title_Name <- paste0("Chromosome ", CHROM_Lab, " ",pretty_peak_id, " " ,"Regional Plot", "\n\n")

      if (is.null(args$Title)) {
        args$Title <- paste0("Chromosome ", CHROM_Lab, " ",pretty_peak_id, " " ,"Regional Plot")
      }

      args$X_Axis_Title <- "Genomic Position (Mb)"



      args <- modifyList(args_base, list(
        Data = filtered_data,
        Title = Title_Name,
        X_Axis_Title = "Genomic Position (Mb)"
      ))

      #args[names(user_args)] <- user_args



      #
      # if(Interactive == TRUE) {
      #   Data$Hover_Info <- paste0(
      #     "CHR: ", Data$CHROM, "\n",
      #     "POS: ", Data$GENPOS, "\n",
      #     "P: ", signif(Data$P, 4)
      #   )
      # }


      print("hi7")

      # Generate plot safely
   #   plot <- tryCatch({


      #will fail if draft plot included as TRUE!

        p <- do.call(Single_Plot, args)

        print("hi8")
#
#         if(Interactive == TRUE)
#
#         {
#           p <- p + ggplot2::aes(text = Hover_Info)
#         }


        # Add Mb-based x-axis
        pos_range <- range(filtered_data$POS, na.rm = TRUE)
        breaks <- pretty(pos_range, n = 5)
        labels <- paste0(round(breaks / 1e6, 1), " Mb")



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
        gene_data <- gene_data[gene_data$start <= max(pos_range) & gene_data$end >= min(pos_range), , drop = FALSE]


        print("Hihalf")

        # Minimum gap (in base pairs) to allow genes on the same line
        min_gene_gap <- MIN_GENE_GAP  # 100 kb

        # Sort by genomic start
        gene_data <- gene_data[order(gene_data$start), ]

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

        print("Hiz")

        # Sort genes by position
        gene_data <- gene_data[order((gene_data$start + gene_data$end) / 2), ]



print("HiA")

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

        if(Auto_LD == TRUE)
        {
        p <- p +  ggplot2::geom_point(
            data = filtered_data,
            ggplot2::aes(x = GENPOS, y = -log10(P), colour = LD_Bin,text = Hover_Info),
            size = 5
          )


        }
       #   ggplot2::scale_x_continuous(
#
 #           expand = ggplot2::expansion(add = c(x_expand, x_expand))
  #        ) +

   p <- p +   ggplot2::scale_x_continuous(
            breaks = NULL,
            labels = NULL,
            position = "bottom",
            #   expand = ggplot2::expansion(mult = c(0.01, 0.01))
            expand = ggplot2::expansion(add = c(x_expand, x_expand)))+

          # Plot lead SNP on top (big purple diamond)
          ggplot2::geom_point(
            data = filtered_data %>% filter(SNP == lead_snp),
            ggplot2::aes(x = GENPOS, y = -log10(P)),
            shape = 18,
            size = 4,
            colour = "purple4"
          )

        if(Auto_LD == TRUE)
        {

        p <- p +  ggplot2::scale_colour_manual(
       #     name = "R^2~LD",
            values = c(
              "NA" = "darkgrey",
              "<0.2" = "darkblue",
              "0.2-0.4" = "lightblue",
              "0.4-0.6" = "green",
              "0.6-0.8" = "yellow",
              "0.8-<1" = "red",
              "1" = "red"
            ),
            na.value = "darkgrey"
          )+
          ggplot2::labs(colour = "R² LD")+
          ggplot2::guides(
            colour = ggplot2::guide_legend(override.aes = list(shape = 16, size = 5)),
            shape = "none",
            size = "none"
          )
        }
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

        p <- p +
          ggplot2::labs(title = Title_Name)+
          ggplot2::theme(
            legend.key.height = ggplot2::unit(1, "cm"),
          #  plot.title = ggplot2::element_text(
          #    hjust = 0.5,         # Center horizontally
          #    vjust = -8        # Optional: make it look nice
          #  ),

           # axis.title.x =    ggplot2::element_blank(),
          #  axis.text.x =    ggplot2::element_blank(),
            legend.position = "right",
            legend.box.just = "left",
            legend.spacing.y = ggplot2::unit(40, "cm"),
                 legend.title = ggplot2::element_text(
                  margin = ggplot2::margin(b = 30),
        size = 30,
        face = "bold",
        family = "FiraCode"
      ),
      legend.text = ggplot2::element_text(
        size = 20,
        family = "FiraCode",
        margin = ggplot2::margin(l = 8, b = 2)# match font for entries too# match font for entries too
      ),
      legend.key = ggplot2::element_blank() # ensures square color blocks


            #diabled originally
       #     axis.ticks.length.x  = ggplot2::unit(0.3, "cm"),
           # axis.ticks.x = element_blank(),
  #    axis.line.x = element_blank(),
    #   plot.title = ggplot2::element_blank(),
         #   plot.margin =    ggplot2::margin(t = 0, b = 0, r = 10, l = 50)
       #     panel.background = element_rect(fill = "white", colour = NA),  # keep panel white
          #  plot.background = element_rect(fill = "grey")  # shade around plot (including below axis)
      #     axis.line.x = ggplot2::element_line(color = "black")

          )

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

        if(Gene_Tracks == FALSE)
        {
        p <- p + ggplot2::scale_x_continuous(
          breaks = breaks,
          labels = labels,
          position = "bottom",
       #   expand = ggplot2::expansion(mult = c(0.01, 0.01))
       expand = ggplot2::expansion(add = c(x_expand, x_expand))
        ) +
          ggplot2::labs(x = args$X_Axis_Title)+ ggplot2::theme(
          plot.margin =    ggplot2::margin(t = 0, b = 50, r = 10, l = 50),
          axis.title.x = ggplot2::element_text(size = args$X_Axis_Title_Size, vjust = -7, colour = "black"),
          axis.text.x  = ggplot2::element_text(size = args$Chromosome_Label_Size, vjust = -2, colour = "black" )
          )
      }



      if(Auto_LD == FALSE)
      {
        p <- p + ggplot2::theme(legend.position = "none")
      }


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


      exons_df2 <- gene_data %>%
        select(transcript_id, label, y, exons) %>%
        tidyr::unnest(exons)

      # Unnest introns (with name repair to fix duplicate column name issue)
      introns_df2 <- gene_data %>%
        select(transcript_id, label, y, introns) %>%
        tidyr::unnest(introns, names_repair = "unique")  #%>%
#        rename(start = intron_start, end = intron_end)


        introns_df2$start <-  introns_df2$intron_start
        introns_df2$end <-  introns_df2$intron_end

      gene_data$gene_biotype <- factor(gene_data$gene_biotype)
      levels(gene_data$gene_biotype) <- gsub("_", " ", levels(gene_data$gene_biotype))



      p_genes <- ggplot2::ggplot(gene_data) +
        # ggplot2::geom_segment(
        #   ggplot2::aes(x = start, xend = end, y = y, yend = y),
        #   size = 4,
        #   color = "black",
        #   alpha = 0.4  # adjust alpha for transparency
        # )+
        ggplot2::geom_segment(data = introns_df2,
                     ggplot2::aes(x = start, xend = end, y = y, yend = y),
                     color = "grey", size = 1) +
        ggplot2::geom_segment(data = exons_df2,
                     ggplot2::aes(x = start, xend = end, y = y, yend = y),
                     size = 6,
                     color = "black") +
        ggplot2::geom_text(
          ggplot2::aes(
            x = (start + end) / 2,
            y = y + 0.375,
            label = label,
            color = gene_biotype #,  # ⬅️ color by biotype

          ),
          show.legend = FALSE,
          angle = 0,
          family = "FiraCode",
          fontface = "plain",
          vjust = 0.5,
          hjust = 0.5,
          size = 6
        ) + ggplot2::scale_color_brewer(
          palette = "Dark2",
          name = "Gene Biotype",
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
          shape = 15, size = 0, alpha = 0, show.legend = FALSE, inherit.aes = FALSE
        )+
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
        ggplot2::scale_x_continuous(
          limits = pos_range,
          breaks = breaks,
          labels = labels,
          expand = ggplot2::expansion(add = c(x_expand, x_expand)),
          sec.axis = ggplot2::dup_axis(
            labels = labels,
            name = args$X_Axis_Title  # Or NULL if you don’t want a title at top
          )
        ) +

    #    ggplot2::scale_y_continuous(limits = c(NULL, max(gene_data$y, na.rm = TRUE) + 1),
     #                               expand = c(0, 0)  )+
      #  ggplot2::theme_minimal(base_size = 10) +
        #   labs(x = NULL, y = NULL) +
        ggplot2::labs(x = args$X_Axis_Title, y = NULL) +
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = args$X_Axis_Title_Size, vjust = -7, colour = "black"),
          axis.text.x  = ggplot2::element_text(size = args$Chromosome_Label_Size, vjust = -2, colour = "black" ),
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
        margin = ggplot2::margin(b = 30, t = 100), #t to control space between legends.
        size = 30,
        face = "bold",
        family = "FiraCode"
      ),
      legend.text = ggplot2::element_text(
        size = 20,
        vjust = 0,
        family = "FiraCode",
        margin = ggplot2::margin(l = 8, b = 8)# match font for entries too
      ),
      legend.text.align = 0,
      legend.key = ggplot2::element_rect(fill = NA)  # ensures square color blocks
        )

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




#print(p_genes)
  #  plot <- p / p_genes   + patchwork::plot_layout(heights = c(3, 2.5))


      if(Auto_LD == TRUE & Gene_Tracks == TRUE)
      {
      #make sure legends align
    plot <- p / p_genes + patchwork::plot_layout(heights = c(3, 3), guides = "collect")
}

    if(Auto_LD == FALSE & Gene_Tracks == TRUE)
    {
   p_genes <-  p_genes + ggplot2::theme(legend.title = ggplot2::element_text(
        margin = ggplot2::margin(b = 10, t = 0), #t to control space between legends.
        size = 30,
        face = "bold",
        family = "FiraCode"
      ))
      plot <- p / p_genes   + patchwork::plot_layout(heights = c(3, 3))
    }

      if(Auto_LD == TRUE & Gene_Tracks == FALSE)
      {
        p_genes <-  p_genes + ggplot2::theme(legend.title = ggplot2::element_text(
          margin = ggplot2::margin(b = 10, t = 0), #t to control space between legends.
          size = 30,
          face = "bold",
          family = "FiraCode"
        ))
        plot <- p / p_genes   + patchwork::plot_layout(heights = c(3, 3))
      }



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

    if(Gene_Tracks == FALSE)
    {
      plot <- p
    }


    if(Gene_Tracks == TRUE)
    {

    n_gene_rows <- max(gene_data$y)

    print(n_gene_rows)

    # Dynamic scaling of total height
    total_height <- 5 + n_gene_rows * 2  # Base height + per-row bonus


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



      return(plot)
    })

    names(plots) <- Title_Names

    return(plots)
  }





  region_plots <- plot_all_peak_regions(region_df, user_args = user_args)


 # return(p_genes)

  return(region_plots)
}




#CHROM <- Regional_Plot(Data = Intelligence_Sum_Stats, Chromosome = 3)
#CHROM <- Single_Plot(Data = Intelligence_Sum_Stats)


#Regional_Plot(Data = Intelligence_Peaks, Condense_Scale = F)
