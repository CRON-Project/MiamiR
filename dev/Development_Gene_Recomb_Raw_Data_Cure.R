

# Load GTF file (offline)
#https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

gtf_path <- "C:/Users/callumon/Downloads/Homo_sapiens.GRCh38.109.gtf.gz"

gtf <- import(gtf_path)

# Convert to data.frame
gtf_df <- as.data.frame(gtf)

colnames(gtf_df)

table(gtf_df$seqnames) #shows chroms/reigons available

#not keeping mt or any scaffolds but could always add in future

CHROMS <- c(seq(1:22), "X", "Y") #Massive, need to loop

for(CHROM in CHROMS)
{

print("Starting...")
print(CHROM)

# Filter to all Chromosome and valid gene names
gene_data <- gtf_df[ gtf_df$seqnames == CHROM &  !is.na(gtf_df$gene_name), ]

#Add strand-aware gene label
gene_data$label <- ifelse(
  gene_data$strand == "+",
  paste0(gene_data$gene_name, " ->"),
  paste0("<- ", gene_data$gene_name)
)

#Filter to exons
exons_all <- gene_data %>%
  filter(type == "exon") %>%
  arrange(gene_id, transcript_id, start)

#Identify longest transcript per gene
longest_transcripts <- exons_all %>%
  group_by(gene_id, transcript_id) %>%
  summarise(tx_start = min(start), tx_end = max(end), .groups = "drop") %>%
  mutate(tx_length = tx_end - tx_start) %>%
  group_by(gene_id) %>%
  slice_max(tx_length, n = 1, with_ties = FALSE) %>%
  ungroup()

#Get exon rows from the longest transcripts only
exons <- exons_all %>%
  semi_join(longest_transcripts, by = c("gene_id", "transcript_id"))

#Calculate introns from exon gaps
introns <- exons %>%
  arrange(transcript_id, start) %>%
  group_by(transcript_id) %>%
  mutate(
    intron_start = lag(end) + 1,
    intron_end = start - 1
  ) %>%
  filter(!is.na(intron_start), intron_end > intron_start) %>%
  ungroup() %>%
  select(transcript_id, intron_start, intron_end)

#Build final summary table
gene_summary <- longest_transcripts %>%
  left_join(
    gene_data %>%
      select(gene_id, gene_name, strand, label, gene_biotype, seqnames) %>%
      distinct(),
    by = "gene_id"
  ) %>%
  mutate(
    exons = map(transcript_id, ~ exons %>% filter(transcript_id == .x) %>% select(start, end)),
    introns = map(transcript_id, ~ introns %>% filter(transcript_id == .x))
  )


gene_summary$start <- gene_summary$tx_start
gene_summary$end <- gene_summary$tx_end

gene_data <- gene_summary

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/inst/extdata/Processed_Gene_Data")

save_name <- paste0("Chromosome_", CHROM, "_", "Gene_Data_HG38_Processed.json")

print("Saving...")
print(CHROM)

write_json(gene_data, save_name, pretty = TRUE)

}


# Load GTF file (offline)
#Newest hg19/37
#https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz


gtf_path <- "C:/Users/callumon/Downloads/Homo_sapiens.GRCh37.75.gtf.gz"
gtf <- import(gtf_path)

# Convert to data.frame
gtf_df <- as.data.frame(gtf)

colnames(gtf_df)

table(gtf_df$seqnames)

#there are other bits not chroms here

CHROMS <- c(seq(1:22), "X", "Y") #Massive, need to loop.

for(CHROM in CHROMS)
{

  print("Starting...")
  print(CHROM)

#Filter to all Chromosome and valid gene names
gene_data <- gtf_df[gtf_df$seqnames == CHROM & !is.na(gtf_df$gene_name), ]

#Add strand-aware gene label
gene_data$label <- ifelse(
  gene_data$strand == "+",
  paste0(gene_data$gene_name, " ->"),
  paste0("<- ", gene_data$gene_name)
)

#Filter to exons
exons_all <- gene_data %>%
  filter(type == "exon") %>%
  arrange(gene_id, transcript_id, start)

#Identify longest transcript per gene
longest_transcripts <- exons_all %>%
  group_by(gene_id, transcript_id) %>%
  summarise(tx_start = min(start), tx_end = max(end), .groups = "drop") %>%
  mutate(tx_length = tx_end - tx_start) %>%
  group_by(gene_id) %>%
  slice_max(tx_length, n = 1, with_ties = FALSE) %>%
  ungroup()

#Get exon rows from the longest transcripts only
exons <- exons_all %>%
  semi_join(longest_transcripts, by = c("gene_id", "transcript_id"))

#Calculate introns from exon gaps
introns <- exons %>%
  arrange(transcript_id, start) %>%
  group_by(transcript_id) %>%
  mutate(
    intron_start = lag(end) + 1,
    intron_end = start - 1
  ) %>%
  filter(!is.na(intron_start), intron_end > intron_start) %>%
  ungroup() %>%
  select(transcript_id, intron_start, intron_end)

#Build final summary table
gene_summary <- longest_transcripts %>%
  left_join(
    gene_data %>%
      select(gene_id, gene_name, strand, label, gene_biotype, seqnames) %>%
      distinct(),
    by = "gene_id"
  ) %>%
  mutate(
    exons = map(transcript_id, ~ exons %>% filter(transcript_id == .x) %>% select(start, end)),
    introns = map(transcript_id, ~ introns %>% filter(transcript_id == .x))
  )

gene_summary$start <- gene_summary$tx_start
gene_summary$end <- gene_summary$tx_end

gene_data <- gene_summary

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/inst/extdata/Processed_Gene_Data")

save_name <- paste0("Chromosome_", CHROM, "_", "Gene_Data_HG19_Processed.json")

print("Saving...")
print(CHROM)

write_json(gene_data, save_name, pretty = TRUE)

}


#LD query API - may max out with users - adapt later.
token <- "a3a5b2b4d4c5"

#https://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recombAvg.bw


setwd("C:/Users/callumon/Downloads")

# Import entire file: backup option beliw
#bw_data <- import("recomb1000GAvg (2).bw")
bw_data <- import("recombAvg.bw")

# Convert GRanges to a tibble
bw_df <- as_tibble(bw_data) %>%
  transmute(
    CHROM = as.character(seqnames),
    start = start,
    end = end,
    score = score
  )

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/inst/extdata/Processed_Recombination_Data")

write.csv(bw_df, file = "Processed_AVG_Recomb_HG38.csv", row.names = FALSE)


#http://hgdownload.soe.ucsc.edu/gbdb/hg19/decode/SexAveraged.bw

setwd("C:/Users/callumon/Downloads")

# Import entire file:
bw_data <- import("SexAveraged.bw")

# Convert GRanges to a tibble
bw_df <- as_tibble(bw_data) %>%
  transmute(
    CHROM = as.character(seqnames),
    start = start,
    end = end,
    score = score
  )

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/inst/extdata/Processed_Recombination_Data")

write.csv(bw_df, file = "Processed_AVG_Recomb_HG19.csv", row.names = FALSE)

#Create backups of RDA dfs in raw format for later demonstration.

# Write Intelligence_Sum_Stats to .txt
write.table(Intelligence_Sum_Stats,
            file = "C:/Users/callumon/Miami_Package_R/MiamiR/inst/extdata/Example_Data_Raw/Intelligence_Sum_Stats.txt",
            sep = "\t",                # Tb-separated
            row.names = FALSE,         # No row names
            quote = FALSE              # No quoting of strings
)

# Write Household_Income_Sum_Stats to .txt
write.table(Household_Income_Sum_Stats,
            file = "C:/Users/callumon/Miami_Package_R/MiamiR/inst/extdata/Example_Data_Raw/Household_Income_Sum_Stats.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE
)

