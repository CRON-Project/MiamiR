
#' METASOFT Wrapper
#'
#' @param Data These are the harmonised summary statistics (either file path(s), environmental object(s) or a variable(s) containing such contents) to be run through the METASFOT wrapper. These must match the METASFOT input file format and can be generated with METASOFT_File_Gen(); defaults to NULL
#' @param options Character vector of extra METASOFT options (e.g. c("-verbose", "mvalue"); defaults to character()/NULL
#' @param REGENIE_Cols Include REGENIE format-style columns; defaults to TRUE
#'
#' @return METASOFT output file is allocated to specified object
#'
#' @export
#'
#' @examples METASOFT_Output <- METASOFT_Run(Data = METASOFT_Example_Input)

METASOFT_Run <- function(Data, options = character(), REGENIE_Cols = TRUE ) {


  #Auto-detect Java path assuming wsl available in terminal

  java_path <- tryCatch(system("wsl which java", intern = TRUE), error = function(e) "")

  if (length(java_path) == 0 || java_path == "") {

    stop("Java not found. Please install Java (>= 8) and ensure it is on PATH.")

  }

  #Locate files in MiamiR package (inst/java)

  jar_path <- system.file("java", "METASOFT.jar", package = "MiamiR")
  p_table  <- system.file("java", "HanEskinPvalueTable.txt", package = "MiamiR")


  #Preview input data (first 10 rows)

  message("Data preview (first 10 rows):")
  print(utils::head(as.data.frame(Data), 10))

  #Create temp input/output

  input_file  <- tempfile(fileext = ".meta")
  output_file <- tempfile(fileext = ".out")

  #Remove header from input DF for METASOFT

  utils::write.table(Data, input_file, quote = FALSE,
                     row.names = FALSE, col.names = FALSE)

  #Convert Windows to WSL paths and quote them

  to_wsl_quoted <- function(path) {
    shQuote(sub("^([A-Za-z]):", "/mnt/\\L\\1", gsub("\\\\", "/", path), perl = TRUE))
  }

  jar_path_wsl     <- to_wsl_quoted(jar_path)
  p_table_wsl      <- to_wsl_quoted(p_table)
  input_file_wsl   <- to_wsl_quoted(input_file)
  output_file_wsl  <- to_wsl_quoted(output_file)

  #Build arguments for WSL call of METASOFT


  args <- c("java", "-jar", jar_path_wsl,
            "-input", input_file_wsl,
            "-output", output_file_wsl,
            "-pvalue_table", p_table_wsl,
            options)

  #Debug

  message("jar_path: ", jar_path)
  message("java_path: ", paste(java_path, collapse = " "))
  message("input_file: ", input_file)
  message("output_file: ", output_file)
  message("args (WSL):"); print(args)

  #Run METASOFT in WSL

  res <- system2("wsl", args = args, stdout = TRUE, stderr = TRUE)

  #Show METASOFT logs first

  message("METASOFT log")

  if (length(res)) cat(paste(res, collapse = "\n"), "\n")

  #Load the METASOFT .out file

  read_metasoft_out <- function(path) {

    #Delimiter

    sep <- "\t"

    #Fallback: pad header to max fields, then read skipping header

    header_line <- tryCatch(readLines(path, n = 1L, warn = FALSE), error = function(e) "")

    header <- if (sep == "\t") {

      strsplit(header_line, "\t", fixed = TRUE)[[1]]

    } else {

      scan(text = header_line, what = character(), quiet = TRUE)

    }

    #Count max fields across lines

    nf <- tryCatch(
      utils::count.fields(path, sep = if (sep == "\t") "\t" else " ",
                          quote = "", blank.lines.skip = TRUE),
      error = function(e) integer()
    )

    if (!length(nf)) return(NULL)

    maxf <- max(nf, na.rm = TRUE)
    need <- maxf - length(header)

    if (need > 0) header <- c(header, paste0("V", seq_len(need)))

    #Read body and set names

    df2 <- tryCatch(
      utils::read.table(path, header = FALSE, sep = sep, quote = "",
                        comment.char = "", stringsAsFactors = FALSE,
                        fill = TRUE, skip = 1),
      error = function(e) NULL
    )

    if (is.null(df2)) return(NULL)

    #Pad columns if still short

    if (ncol(df2) < length(header)) {

      df2[(length(header))] <- NA

    }

    names(df2) <- header[seq_len(ncol(df2))]
    df2

  }

  out_df <- read_metasoft_out(output_file)

  out_df$V1 <- NULL
  out_df$V2 <- NULL
  out_df$V3 <- NULL

  if (isTRUE(REGENIE_Cols)) {


    rs <- as.character(out_df$RSID)

    m  <- regexec("^(?:chr)?([^:]+):([^:]+):([^:]+):([^:]+)$", rs, perl = TRUE)
    mm <- regmatches(rs, m)

    mat <- do.call(rbind, lapply(mm, function(x)

      if (length(x) == 5) x[-1] else rep(NA_character_, 4)

    ))

    colnames(mat) <- c("CHROM","GENPOS","ALLELE0","ALLELE1")

    out_df$CHROM   <- mat[, "CHROM"]
    out_df$GENPOS  <- suppressWarnings(as.numeric(mat[, "GENPOS"]))  # numeric position
    out_df$ALLELE0 <- mat[, "ALLELE0"]
    out_df$ALLELE1 <- mat[, "ALLELE1"]

    out_df$P <- out_df$PVALUE_FE
    out_df$BETA <- out_df$BETA_FE
    out_df$SE <- out_df$STD_FE

  }


  if (!is.null(out_df)) {

    message("Output preview (first 10 rows)")

    print(utils::head(out_df, 10))

  } else {

    message("Output file exists but could not be parsed cleanly. You can inspect it at: ", output_file)

  }

  #Return results

  return(out_df)

}
