
#' METASOFT Wrapper
#'
#' @param Data These are the harmonised summary statistics (either file path(s), environmental object(s) or a variable(s) containing such contents) to be run through the METASFOT wrapper. These must match the METASFOT input file format and can be generated with METASOFT_File_Gen(); defaults to NULL
#' @param options Character vector of extra METASOFT options (e.g. c("-verbose", "mvalue"); defaults to character()/NULL
#' @param REGENIE_Cols Include REGENIE format-style columns in output; defaults to TRUE
#' @param Verbose Prevent display of progress bar as function is running and instead show key milestone outputs/messages (mainly for debugging purposes); defaults to FALSE.
#'
#' @return METASOFT output file is allocated to specified object
#'
#' @export
#'
#' @examples METASOFT_Output <- METASOFT_Run(Data = METASOFT_Example_Input)

  METASOFT_Run <- function(Data = NULL, options = character(), REGENIE_Cols = TRUE, Verbose = FALSE ) {

  # Allow Data to be a file path

     if (is.character(Data) && length(Data) == 1L && grepl("/", Data, fixed = TRUE)) {

      file_path <- Data

      message("Dataset absent from environment")

      if (file.exists(file_path)) {

        message("Reading data from file: ", file_path)

        message("Loading data using vroom...")

        suppressWarnings(

        suppressMessages(

        Data <- vroom::vroom(file_path, show_col_types = FALSE, progress = FALSE)

        )

        )

        message("Finished reading")

      } else {

        stop("The provided character string does not point to an existing file: ", file_path, call. = FALSE)

      }

    }

  # Auto-detect Java path assuming wsl available in terminal

  message("Finding java path")

  if(tolower(Sys.info()[["sysname"]]) == "linux")

  {

    java_path <- tryCatch(system("which java", intern = TRUE), error = function(e) "")

  }else{

    java_path <- tryCatch(system("wsl which java", intern = TRUE), error = function(e) "")

  }

  if (length(java_path) == 0 || java_path == "") {

    stop("Java not found. Please install Java (>= 8) and ensure it is on PATH.")

  }

  message("Located functional path, locating packaged software")

  # Locate files in MiamiR package (inst/java)

  jar_path <- system.file("java", "METASOFT.jar", package = "MiamiR")
  p_table  <- system.file("java", "HanEskinPvalueTable.txt", package = "MiamiR")

  # Preview input data (first 10 rows)

  message("Gathering arguments")

  # Create temp input/output

  input_file  <- tempfile(fileext = ".meta")
  output_file <- tempfile(fileext = ".out")

  # Remove header from input DF for METASOFT

  utils::write.table(Data, input_file, quote = FALSE,
                     row.names = FALSE, col.names = FALSE)

  # Convert Windows to WSL paths and quote them

  # Convert Windows paths to quoted WSL paths (for use with `wsl`)

  to_wsl_quoted <- function(path) {

    path_win <- tryCatch(

      normalizePath(path, winslash = "\\", mustWork = FALSE),
      error = function(e) path

    )

    if (!nzchar(path_win)) path_win <- path

    shQuote(sub("^([A-Za-z]):", "/mnt/\\L\\1", gsub("\\\\", "/", path_win), perl = TRUE))

  }

  if (tolower(Sys.info()[["sysname"]]) == "linux") {

    # Native Linux: just call java directly with normal paths

    args <- c(

      "-jar", jar_path,
      "-input", input_file,
      "-output", output_file,
      "-pvalue_table", p_table,
      options

    )


    message("jar_path: ", jar_path)
    message("java_path: ", paste(java_path, collapse = " "))
    message("input_file: ", input_file)
    message("output_file: ", output_file)

    message("Running METASOFT meta-analysis...")

    res <- system2("java", args = args, stdout = TRUE, stderr = TRUE)


  } else {

    # Windows via WSL: convert all paths to /mnt/c/... and quote them

    jar_path_wsl    <- to_wsl_quoted(jar_path)
    p_table_wsl     <- to_wsl_quoted(p_table)
    input_file_wsl  <- to_wsl_quoted(input_file)
    output_file_wsl <- to_wsl_quoted(output_file)

    args <- c(

      "java", "-jar", jar_path_wsl,
      "-input", input_file_wsl,
      "-output", output_file_wsl,
      "-pvalue_table", p_table_wsl,
      options

    )


    message("jar_path: ", jar_path)
    message("java_path: ", paste(java_path, collapse = " "))
    message("input_file: ", input_file)
    message("output_file: ", output_file)

    message("Running METASOFT meta-analysis...")

    res <- system2("wsl", args = args, stdout = TRUE, stderr = TRUE)

  }

  # Debug

  # Fail early if METASOFT didn't produce an output file

  if (!file.exists(output_file)) {

    stop("METASOFT failed: output file was not created. See log above.")

  }

  # Show METASOFT logs first

  message("METASOFT log of call:")

  if(Verbose == TRUE)

  {

  if (length(res)) cat(paste(res, collapse = "\n"), "\n")

  }

  out_df <- suppressWarnings(

    suppressMessages(

      vroom::vroom(output_file, show_col_types = FALSE, progress = FALSE)

    )

  )

  # Locate MVALUES column

  mcol <- base::grep("^MVALUES_OF_STUDIES", base::names(out_df), value = TRUE)

  if (length(mcol) != 1L) {

    base::stop("Couldn't uniquely find the MVALUES column. Found: ", base::paste(mcol, collapse = ", "))

  }

  # Clean and count tokens

  out_df[[mcol]] <- stringr::str_squish(base::as.character(out_df[[mcol]]))

  max_tokens <- base::max(stringr::str_count(out_df[[mcol]], "\\S+") + 0L, na.rm = TRUE)

  # split into separate columns

  out_df <- tidyr::separate_wider_delim(

    out_df,
    cols   = dplyr::all_of(mcol),
    delim  = " ",
    names  = base::paste0("MVALUE_", base::seq_len(max_tokens)),
    too_few = "align_start"

  )

  if (isTRUE(REGENIE_Cols)) {

    message("Adding REGENIE Style information columns")

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

    message("Outputting file...")

  } else {

    message("Output file exists but could not be parsed cleanly. You can inspect it at: ", output_file)

  }

  # Return results

  return(out_df)

  }

  if (!exists("use_wrapper")) use_wrapper <- TRUE

  if (isTRUE(use_wrapper)) {

    .METASOFT_Run_original <- METASOFT_Run

    METASOFT_Run <- function(..., session = NULL) {

      # Capture the *exact* expressions supplied in ...

      dots_expr <- base::as.list(base::substitute(base::list(...)))[-1]

      args <- base::list(...)
      args$session <- session
      args$.dots   <- args

      # helpers

      resolve_data_name <- function(name) {

        if (!base::is.character(name) || base::length(name) != 1L) return(name)

        pf <- base::parent.frame()

        if (base::exists(name, envir = pf, inherits = TRUE)) {

          return(base::get(name, envir = pf, inherits = TRUE))

        }

        if (base::exists(name, envir = .GlobalEnv, inherits = FALSE)) {

          return(base::get(name, envir = .GlobalEnv, inherits = FALSE))

        }

        name
      }

      is_multi_data <- function(x) {

        if (base::is.null(x)) return(FALSE)

        if (base::is.data.frame(x)) return(FALSE)

        if (base::is.list(x)) return(base::length(x) > 1L)

        if (base::is.atomic(x) && base::length(x) > 1L) return(TRUE)

        FALSE

      }

      as_data_list <- function(x) {

        if (base::is.data.frame(x)) return(base::list(x))

        if (base::is.list(x))       return(x)

        if (base::is.atomic(x) && base::length(x) > 1L) return(base::as.list(x))

        base::list(x)
      }

      pick_i <- function(val, i, n) {

        if (base::is.null(val)) return(NULL)

        if (base::is.list(val) && !base::is.data.frame(val)) {

          if (base::length(val) == 1L) return(val[[1L]])

          if (base::length(val) == n)  return(val[[i]])

          return(val[[1L]])

        }

        if (base::is.atomic(val) && base::length(val) > 1L) {

          if (base::length(val) == n) return(val[i])

          return(val[1L])

        }

        val
      }

      infer_names <- function(dlist) {

        n  <- base::length(dlist)
        nm <- base::names(dlist)

        if (!base::is.null(nm) && base::any(nzchar(nm))) {

          nm[nm == ""] <- base::paste0("run_", base::which(nm == ""))

          return(nm)

        }

        flat <- base::unlist(dlist, recursive = FALSE, use.names = FALSE)

        if (base::length(flat) == n && base::all(base::vapply(flat, base::is.character, logical(1)))) {

          out <- base::basename(flat)
          out[out == "" | base::is.na(out)] <- base::paste0("run_", base::seq_len(n))

          return(out)

        }

        base::paste0("run_", base::seq_len(n))

      }

      label_from_expr <- function(expr) {

        if (base::is.null(expr)) return("run_1")

        lab <- base::paste(base::deparse(expr, width.cutoff = 500L), collapse = "")

        # if user literally typed a quoted string path, show basename

        if (base::grepl('^".*"$', lab)) {

          lab <- base::gsub('^"|"$', "", lab)

          if (nzchar(lab)) lab <- base::basename(lab)

        }

        if (!nzchar(lab)) "run_1" else lab
      }

      send_pause <- function() {

        if (!base::is.null(session)) {

          session$sendCustomMessage(

            "plot_progress",
            base::list(pct = NA, msg = "Running METASOFT meta-analysis")

          )

        }

      }

      # run one file/df

      run_one <- function(run_args) {

        if ("Data" %in% base::names(run_args)) {

          run_args$Data <- resolve_data_name(run_args$Data)

        }

        verbose_mode <- if ("Verbose" %in% base::names(run_args)) isTRUE(run_args$Verbose) else FALSE

        run_args$session <- NULL
        run_args$.dots   <- NULL

        # silent progress bar

        .bar_draw <- function(pct) {

          pct_i <- base::as.integer(base::floor(base::max(0, base::min(100, pct))))
          bar_len <- 30L
          filled  <- base::as.integer(base::round(bar_len * pct_i / 100))

          base::cat(

            "\r[", base::strrep("=", filled), base::strrep(" ", bar_len - filled), "] ",
            base::sprintf("%3d%%", pct_i),
            sep = ""

          )

          utils::flush.console()

        }

        .bar_done <- function() {

          base::cat("\n")
          utils::flush.console()

        }

        if (verbose_mode) {

          return(base::do.call(.METASOFT_Run_original, run_args))

        }

        # Verbose = FALSE: bar + METASOFT stdout/stderr only

        .bar_draw(0)

        out <- withCallingHandlers(

          base::do.call(.METASOFT_Run_original, run_args),

          message = function(m) {

            msg <- base::conditionMessage(m)

            # show + trigger pause only for this milestone

            if (base::grepl("^Running METASOFT meta-analysis", msg)) {

              send_pause()
              .bar_draw(50)

              return()  # DO NOT muffle this need for running pause

            }

            # muffle everything else in Verbose=FALSE

            base::invokeRestart("muffleMessage")

          }


        )

        .bar_draw(100)
        .bar_done()

        out

      }

      # MULTI-RUN

      if ("Data" %in% base::names(args) && is_multi_data(args$Data)) {

        Data_list <- as_data_list(args$Data)
        n <- base::length(Data_list)

        out_names <- infer_names(Data_list)

        results <- base::vector("list", n)
        base::names(results) <- out_names

        base_args <- args
        base_args$Data <- NULL

        for (i in base::seq_len(n)) {

          run_args <- base_args
          run_args$Data <- Data_list[[i]]

          for (nm in base::names(run_args)) {

            if (nm %in% c("session", ".dots")) next

            run_args[[nm]] <- pick_i(run_args[[nm]], i, n)

          }

          base::message(base::sprintf("Processing dataset %d/%d: %s", i, n, out_names[i]))

          results[[i]] <- run_one(run_args)
        }


        return(results)

      }

      # SINGLE

      # Use the *expression user typed* for Data

      supplied_lab <- "run_1"

      if ("Data" %in% base::names(dots_expr)) {

        supplied_lab <- label_from_expr(dots_expr$Data)

      } else if (base::length(dots_expr) >= 1L) {

        # handle unnamed first arg: METASOFT_Run(METASOFT_Example_Input)

        supplied_lab <- label_from_expr(dots_expr[[1L]])

      }

      base::message(base::sprintf("Processing dataset 1/1: %s", supplied_lab))

      if ("Data" %in% base::names(args)) {

        args$Data <- resolve_data_name(args$Data)

      }

      run_one(args)

    }

  }
