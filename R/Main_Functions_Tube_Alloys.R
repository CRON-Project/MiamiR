
#' Visualise GWAS Data
#'
#' @param Data These are the summary statistics (either file path(s), environmental object(s) or a variable(s) containing such contents) to be plotted; defaults to NULL
#' @param Chromosomes Manually specify which chromosomes from the summary statistics to run the chromosome sensitive functions on; defaults to NULL, leading to selective default filtering
#' @param Phenos Specify GWAS phenotype file, when Phenos_Covars is FALSE - used for GWAS model overview Forest_Plot(); defaults to NULL
#' @param Covars Specify GWAS covariate file, when Phenos_Covars is FALSE - used for GWAS model overview Forest_Plot(); defaults to NULL
#' @param Phenos_Covars Specify combined phenotype & covariate file, when Phenos & Covars are both FALSE - used for GWAS model overview Forest_Plot(); defaults to NULL
#' @param Pheno_Name When Phenos_Covars is TRUE, specify the name of the column containing the GWAS Phenotype - used for GWAS model overview Forest_Plot(); defaults to NULL
#'
#' @return Images of Manhattan Plot(s), Regional Plots(s) and Forest Plot(s) are allocated to specified object and the resulting ggplot and grob object(s) can then be saved to an image
#' @export
#'
#' @examples Plots <- Tube_Alloys(Data = Intelligence_Sum_Stats, Phenos = Fake_PHENOS_Binary, Covars = Fake_COVARS, Chromosomes = 1 )
#'

Tube_Alloys <- function(Data = NULL, Chromosomes = NULL, Phenos = NULL, Covars = NULL,
                        Phenos_Covars = NULL, Pheno_Name = NULL)

{

   data_sym  <- substitute(Data)
   data_name <- deparse(data_sym)

   #Works out if multiple data sets are passed
   a <- grepl(".+,.*", data_name)

   parsed <- parse(text = data_name)[[1]]

   # Extract the symbols/strings inside c()
   dataset_names <- vapply(as.list(parsed)[-1L], deparse1, character(1))

   #Regardless gets clean names of dfs/paths
   x_clean <- sub("^c\\((.*)\\)$", "\\1", data_name)

   message("Data object passed: ", x_clean)


   message("Running Single_Plot()")
   SinglePlot <- eval(bquote(Single_Plot(Data = .(data_sym))))

   message("Running Regional_Plot()")
   RegionalPlot <- eval(bquote(Regional_Plot(Data = .(data_sym), Chromosomes = Chromosomes)))

   message("Running Forest_Plot()")

   #a signifies multiple files
   if(a == TRUE)

   {

   ForestPlot <- Forest_Plot(Data = c(dataset_names), Chromosomes = Chromosomes, Model_Reference = F)

   }else{

   ForestPlot <- Forest_Plot(Data = x_clean, Chromosomes = Chromosomes, Model_Reference = F)

   }

   #define load function

   load_input <- function(x) {

     if (is.null(x)) return(NULL)

     if (exists(x, envir = .GlobalEnv)) {

       get(x, envir = .GlobalEnv)

     } else if (file.exists(x)) {

       message("Reading from file path: ", x)
       vroom::vroom(x, show_col_types = FALSE)

     } else {

       stop("The name or file path '", x, "' does not exist.")

     }
   }

   #Define binary-like phenotype detection function

   is_binary <- function(x) {
     u <- unique(x[!is.na(x)])
     length(u) == 2
   }

   message("Checking for accessory GWAS file data")

   Phenos <- if (is.character(Phenos)) load_input(Phenos) else Phenos
   Covars <- if (is.character(Covars)) load_input(Covars) else Covars
   Phenos_Covars <- if (is.character(Phenos_Covars)) load_input(Phenos_Covars) else Phenos_Covars

   if( (!is.null(Phenos) & !is.null(Covars) &  !is.null(Phenos_Covars) ) | (!is.null(Phenos) & !is.null(Covars) ) | (!is.null(Phenos_Covars) ))
   {

   if (is.null(Phenos_Covars) && !is.null(Phenos) && !is.null(Covars)) {

     message("Binding covariate and phenotype data together by identified participant ID columns")

     get_match <- function(df, wanted) {

       cn <- names(df); up <- toupper(cn)
       sapply(wanted, function(w) {

         i <- match(toupper(w), up)
         if (is.na(i)) NA_character_ else cn[i]

       }, USE.NAMES = TRUE)

     }

     #Phenotype column from Phenos (first non-ID col)
     id_cols <- c("FID","IID","ID","SampleID","eid")
     pheno_col <- setdiff(names(Phenos), id_cols)[1]

     #Initially try 2-key join if both sides have FID & IID
     pair <- c("FID","IID")
     x_pair <- get_match(Phenos, pair)
     y_pair <- get_match(Covars, pair)

     if (all(!is.na(x_pair)) && all(!is.na(y_pair))) {

       cov_tmp <- Covars
       names(cov_tmp)[match(y_pair, names(cov_tmp))] <- x_pair
       by_keys <- unname(x_pair)
       Phenos_Covars <- merge(Phenos, cov_tmp, by = by_keys, all = TRUE)
       message("Joined Phenos(", paste(by_keys, collapse = "+"),
               ") <-> Covars(", paste(unname(y_pair), collapse = "+"), ") [2 keys]")
     } else {

       #Else use any single ID on each side (names can differ)
       singles <- id_cols
       x_s <- get_match(Phenos, singles)
       y_s <- get_match(Covars, singles)
       x_key <- unname(x_s[!is.na(x_s)])[1]
       y_key <- unname(y_s[!is.na(y_s)])[1]

       if (length(x_key) && length(y_key)) {

         cov_tmp <- Covars
         names(cov_tmp)[match(y_key, names(cov_tmp))] <- x_key
         Phenos_Covars <- merge(Phenos, cov_tmp, by = x_key, all = TRUE)
         message("Joined Phenos(", x_key, ") <-> Covars(", y_key, ") [1 key]")

       } else {

         message("Could not find ID cols to join (looked for: ", paste(id_cols, collapse=", "), ").")

           }
     }


     if (!is.null(Phenos_Covars) && !is.na(pheno_col)) {

       covar_cols <- setdiff(names(Phenos_Covars), c(id_cols, pheno_col))

       form <- as.formula(paste(pheno_col, "~", paste(covar_cols, collapse = " + ")))

       if (is_binary(Phenos_Covars[[pheno_col]])) {

         message("Running glm binomial model")
         message(form)

         fit <- glm(form, data = Phenos_Covars, family = binomial())

       } else {

         message("Running lm model")
         message(form)

         fit <- lm(form, data = Phenos_Covars)

       }

       message("Model fit for phenotype: ", pheno_col)

     }
   }

   id_cols <- c("FID","IID","ID","SampleID","eid")

   if (!is.null(Phenos_Covars) && is.null(Phenos) && is.null(Covars)) {

     message("Provided GWAS accessory files already have joined covariate and phenotype data")
     message(paste0("Using:"), " ", Phenos, " ", "as phenotype")

     # Resolve phenotype name
     if (!is.null(Pheno_Name)) {

       # case-insensitive match
       nm <- names(Phenos_Covars)
       i <- match(toupper(Pheno_Name), toupper(nm))

       if (is.na(i)) stop("Pheno_Name '", Pheno_Name, "' not found in Phenos_Covars.")

       ph_col <- nm[i]

     } else {

       ph_col <- setdiff(names(Phenos_Covars), id_cols)[1]
       if (is.na(ph_col)) stop("Could not infer phenotype column; please provide Pheno_Name.")

     }

     covar_cols <- setdiff(names(Phenos_Covars), c(id_cols, ph_col))

     if (!length(covar_cols)) stop("No covariate columns found.")

     form <- as.formula(paste(ph_col, "~", paste(covar_cols, collapse = " + ")))
     fit <- if (is_binary(Phenos_Covars[[ph_col]])) {

       message("Running glm binomial model")
       message(form)

       glm(form, data = Phenos_Covars, family = binomial())

     } else {

       message("Running glm binomial model")
       message(form)

       lm(form, data = Phenos_Covars)

     }

   }

   tmp_model_name <- paste0("Covariate_Phenotype_Data")
   assign(tmp_model_name, fit, envir = .GlobalEnv)

   message("Running Covariate File through Forest_Plot()")
   GWASModelPlot <- Forest_Plot(Data = tmp_model_name, Model_Reference = T, Verbose = T)


   }

   message("Collecting valid available outputs")

   out <- mget(
     c("SinglePlot", "RegionalPlot", "ForestPlot", "GWASModelPlot"),
     ifnotfound = list(NULL, NULL, NULL, NULL)
   )

   # Remove NULLs and empty lists
   out <- out[vapply(out, function(x) !is.null(x) && !(is.list(x) && length(x) == 0), logical(1))]

   return(out)

}
