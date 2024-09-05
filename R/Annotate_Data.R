#' Title Use biomaRt to obtain RSID from COORDS if required before plotting
#'
#' @param Data This is the name of the data frame to annotate; defaults to Intelligence_Sum_Stats
#' @param Chromosome_Column Manually specify chromosome column; defaults to "CHROM"
#' @param Position_Column Manually specify chromosome position column; defaults to "GENPOS"
#' @param SNP_ID_Column  Manually specify SNP ID column; defaults to "ID"
#' @param PValue_Column Manually specify P Value column; defaults to "P"
#' @param Build Genome reference to query biomaRt for ("HG38" or "HG19"); defaults to "HG38"
#' @param Reference_Allele The reference allele in the GWAS summary statistics; defaults to "ALLELE0"
#' @param Effect_Allele The alternate/test/effect allele in the GWAS summary statistics; defaults to "ALLELE1"
#'
#' @return Outputs the input dataframe with a column called Lab with
#'         the corresponding RSIDs of the index SNP per chromosome
#'
#' @export
#'
#' @examples  Labelled_Data <- Annotate_Data(Data = Intelligence_Sum_Stats,
#'            Chromosome_Column = "CHROM",
#'            Position_Column = "GENPOS", SNP_ID_Column = "ID",
#'            PValue_Column = "P", Build = "HG38")
#'
#'

Annotate_Data <- function(Data = Intelligence_Sum_Stats,
                        Chromosome_Column = "CHROM",
                        Position_Column = "GENPOS", SNP_ID_Column = "ID",
                        PValue_Column = "P", Build = "HG38",
                        Reference_Allele = "ALLELE0", Effect_Allele = "ALLELE1")


{


  allowed_names_chromosomes <- c("A0", "A2", "REF", "a0", "a2", "ref", "Reference", "reference", "allele0", "allele2", "ALLELE0", "ALLELE2")

  for (allowed_name in allowed_names_chromosomes) {
    if (allowed_name %in% colnames(Data)) {
      usable_ref_allele <- colnames(Data)[which(colnames(Data) == allowed_name)]
      break
    }
  }



  print(paste0("Using", " ", allowed_name, " ", "as reference allele"))
  Reference_Allele <- allowed_name



  allowed_names_chromosomes <- c("A1", "REF", "a1", "ref", "Reference", "reference", "allele1", "ALLELE1")

  for (allowed_name in allowed_names_chromosomes) {
    if (allowed_name %in% colnames(Data)) {
      usable_effect_allele <- colnames(Data)[which(colnames(Data) == allowed_name)]
      break
    }
  }



  print(paste0("Using", " ", allowed_name, " ", "as effect allele"))
  Effect_Allele <- allowed_name



  allowed_names_chromosomes <- c("chromosome", "chrom", "chr", "CHROM", "Chromosome", "CHR", "Chr", "Chrom")

  for (allowed_name in allowed_names_chromosomes) {
    if (allowed_name %in% colnames(Data)) {
      usable_chrom_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
      break
    }
  }



  print(paste0("Using", " ", allowed_name, " ", "as chromosome column"))
  Chromosome_Column <- allowed_name



  allowed_names_pos <- c("POS", "pos", "Pos", "Position", "position", "POSITION", "genpos", "GENPOS",
                         "Genpos")

  for (allowed_name in allowed_names_pos) {
    if (allowed_name %in% colnames(Data)) {
      usable_pos_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as position column"))
  Position_Column <- allowed_name



  allowed_names_SNP <- c("ID", "Id", "ID", "RsID", "RsId","RSID", "snp", "SNP", "Snp",
                         "snv" ,"SNV" , "Snv", "RS", "rs")

  for (allowed_name in allowed_names_SNP) {
    if (allowed_name %in% colnames(Data)) {
      usable_snp_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as SNP column"))
  SNP_ID_Column <- allowed_name







  allowed_names_pvals <- c("P", "p", "Pvalue", "pvalue", "P-Value", "p-value", "p-Value",
                           "P-VALUE", "logp","LogP", "LOGP", "Logp", "log10p","Log10P",
                           "LOG10P", "Log10p",
                           "log10p", "LOG10P", "-LOG10P", "")

  for (allowed_name in allowed_names_pvals) {
    if (allowed_name %in% colnames(Data)) {
      usable_p_top <- colnames(Data)[which(colnames(Data) == allowed_name)]
      break
    }
  }


  print(paste0("Using", " ", allowed_name, " ", "as P-Val column"))
  PValue_Column <- allowed_name





  if (!("P" %in% colnames(Data)) & any(colnames(Data) %in% c("logp", "LogP", "LOGP", "Logp",
                                                             "log10p", "Log10P", "LOG10P",
                                                             "Log10p", "-LOG10P"))) {

    print("No P Value in first dataset, Calculating from LOG10P column detected")
    Data$P <- 10^-(as.numeric(Data[[PValue_Column]]))
    PValue_Column <- "P"

  }





  #Manually assign columns for ease of use
  Data$CHROM <- Data[[Chromosome_Column]]
  Data$GENPOS <- Data[[Position_Column]]
  Data$ID <- Data[[SNP_ID_Column]]
  Data$P <- Data[[PValue_Column]]
  Data$REF <- Data[[Reference_Allele]]
  Data$ALT <- Data[[Effect_Allele]]
  Data$REF <- toupper(Data[[Reference_Allele]]) #make capital always to match biomart search
  Data$ALT <- toupper(Data[[Effect_Allele]])
  Data$REF_ALT <- paste0(Data$REF, "/", Data$ALT)
  Data$ALT_REF <- paste0(Data$ALT, "/", Data$REF)



  #Get min P row per chromosome - options in plotting to include/exclude
  Data <- Data %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(min_P_GENPOS = GENPOS[which.min(P)]) %>% #P is formed in function
    dplyr::ungroup()




  # Label the minimum P value for each CHROM
  Data <- Data %>%
    dplyr::mutate(Lab = ifelse(GENPOS == min_P_GENPOS, ID, ""))




  #Look for ones with index chromosome IDs
  SNPs <- Data[Data$Lab != "",]


  #Remove this as need this column later
  Data$Lab <- NULL



  #Get biomaRt search format
  SNPs <- SNPs %>% dplyr::select(CHROM, GENPOS, GENPOS)


  SNPs$CHR <- SNPs$CHROM
  SNPs$chr_start <- SNPs$GENPOS
  SNPs$chr_end <- SNPs$GENPOS

  #Don't need these in search file
  SNPs$CHROM <- NULL
  SNPs$GENPOS <- NULL



  #Create SNP mart object depending on genome build

  if(Build == "HG19")
  {

  snp_mart <- biomaRt::useEnsembl(biomart="ENSEMBL_MART_SNP",
                                               host="https://grch37.ensembl.org",
                                              dataset="hsapiens_snp")

  }

  if(Build == "HG38")
  {


   snp_mart <- biomaRt::useEnsembl(biomart="ENSEMBL_MART_SNP", #38 override - make command to spec
                          dataset="hsapiens_snp")



  }



   position <- apply(SNPs, 1, paste, collapse = ":")


   coords <- position

   coords <- gsub("^23:", "X:", coords)




   # Initialize an empty dataframe to store the results
   c <- data.frame()

   # Loop through each row of the coords object
   for (i in 1:length(coords)) {
     # Extract the current row (chromosomal region) from coords
     current_coord <- coords[i]

     print(current_coord)
     # Query BioMart with the current coordinate
     a <- biomaRt::getBM(
       attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'allele'),
       filters = c('chromosomal_region'),
       values = current_coord,
       mart = snp_mart
     )

     # Convert the result to a dataframe
     b <- as.data.frame(a)





     b2 <- dplyr::left_join(b, Data, by = c("chrom_start" = "GENPOS")) #link df alleles

     #only need to do this if multiple matches, as SNPs may be unified under one rs
     if(nrow(b2) > 1)
     {

     print("Multiple matches, filering for SNP by position")


      b2 <- b2[b2$min_P_GENPOS == b2$chrom_start & b2$min_P_GENPOS == b2$chrom_end,  ]
      b2 <- b2[!is.na(b2$refsnp_id),] #NA issue for cols not labelled

      if(nrow(b2) > 1)
      {
        print("Multiple matches still, filering for SNPs by alleles")
      b2 <- b2[b2$ALT_REF == b2$allele | b2$REF_ALT == b2$allele,] #keep the RS matching the alleles if multiple matches still after pos's checked
      }

     }


     print("Match:")
     print(b2$refsnp_id)

     # Combine the current result with the previously accumulated results
     c <- rbind(c, b2)

     print("")
     progress <- paste0("Obtained Index SNP RS code for ", i, " out of ", length(coords), " Index SNPs")
     print(progress)
   }

   #Get the RSID and the coords which are one base to rejoin
   c <- c %>% dplyr::select(refsnp_id, chrom_start) %>%
     #Format to Lab as this is needed in plot functions later
     dplyr::rename(Lab = refsnp_id)



#Join SNPs to summary stats
Data <- dplyr::left_join(Data, c, by = c("GENPOS" = "chrom_start"))

print("The Following RSIDs were assigned to index SNPs")
print("")
print(table(Data$Lab))

print("Returning annotated summary stats")

#remove redundant cols
Data$ID <- NULL
Data$GENPOS <- NULL
Data$CHROM <- NULL
Data$min_P_GENPOS <- NULL
Data$REF <- NULL
Data$ALT <- NULL
Data$REF_ALT <- NULL
Data$ALT_REF <- NULL

#Return the df itself
return(Data)

}


