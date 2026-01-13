
  # Minimal GWAS sim: covars + binary / quantitative phenotypes

  simulate_gwas_covariates <- function(

      n = 2000,
      locations = c("UK","IE","PL"),
      loc_probs = c(0.6, 0.25, 0.15),
      pc_loc_means = list(UK=c(0.20,0.10), IE=c(-0.10,-0.05), PL=c(-0.25,0.15)),
      pc_sd_within = 0.08,
      age_by_loc = list(UK=c(mean=52,sd=10), IE=c(mean=49,sd=11), PL=c(mean=47,sd=9)),
      sex_prob_male = 0.48,
      bmi_base = 26, bmi_age_slope = 0.05, bmi_sd = 3.5,
      seed = 123

  ) {

    stopifnot(abs(sum(loc_probs) - 1) < 1e-8)

    set.seed(seed)

    IID <- sprintf("ID%06d", 1:n); FID <- IID
    loc <- sample(locations, n, TRUE, loc_probs)
    mu <- t(vapply(loc, function(L) pc_loc_means[[L]], numeric(2L)))
    PC1 <- rnorm(n, mu[,1], pc_sd_within); PC2 <- rnorm(n, mu[,2], pc_sd_within)
    ap <- t(vapply(loc, function(L) age_by_loc[[L]], numeric(2L)))
    Age <- round(pmax(18, rnorm(n, ap[,1], ap[,2])))
    Sex <- ifelse(runif(n) < sex_prob_male, 1L, 2L)  # 1=male, 2=female
    BMI <- rnorm(n, mean = bmi_base + bmi_age_slope * (Age - mean(Age)), sd = bmi_sd)

    data.frame(

      FID = FID, IID = IID,
      Age = Age, Sex = Sex,
      Location = factor(loc, levels = locations),
      BMI = round(BMI, 2),
      PC1 = PC1, PC2 = PC2,
      check.names = FALSE

    )

  }

  .build_X <- function(covars) model.matrix(~ Age + Sex + BMI + PC1 + PC2 + Location, data = covars)

  simulate_quant_pheno <- function(

      covars,
      betas = c(Age=0.015, Sex=0.20, BMI=0.10, PC1=-0.30, PC2=0.10, LocationIE=0.20, LocationPL=-0.10),
      target_R2 = 0.30,
      seed = 1

  ) {

    stopifnot(target_R2 > 0, target_R2 < 1); set.seed(seed)
    X <- .build_X(covars); cols <- setdiff(colnames(X), "(Intercept)")
    b <- setNames(rep(0, length(cols)), cols); b[names(betas)] <- betas[names(betas)]

    lin <- drop(X[, cols, drop=FALSE] %*% b); vlin <- var(lin)

    sigma <- sqrt(vlin * (1 - target_R2) / target_R2)
    y <- lin + rnorm(nrow(X), 0, sigma)
    data.frame(

      FID = covars$FID, IID = covars$IID,
               PHENO_QUANT = as.numeric(scale(y, center=TRUE, scale=FALSE))

      )

  }

  simulate_binary_pheno <- function(

      covars,
      betas = c(Age=0.015, Sex=-0.15, BMI=0.08, PC1=0.20,
                PC2=-0.05, LocationIE=0.25, LocationPL=-0.10),
      prevalence = 0.25,
      seed = 2

  ) {

    stopifnot(prevalence > 0, prevalence < 1); set.seed(seed)
    X <- .build_X(covars); cols <- setdiff(colnames(X), "(Intercept)")
    b <- setNames(rep(0, length(cols)), cols); b[names(betas)] <- betas[names(betas)]

    lin <- drop(X[, cols, drop=FALSE] %*% b)

    a0 <- uniroot(function(a) mean(plogis(a + lin)) - prevalence, c(-15, 15))$root

    p <- plogis(a0 + lin); y01 <- rbinom(nrow(X), 1, p)
    data.frame(FID = covars$FID, IID = covars$IID, PHENO_BIN = y01)

  }

  # Run once

  covars <- simulate_gwas_covariates(n = 5000, seed = 42)
  ph_q  <- simulate_quant_pheno(covars, target_R2 = 0.30, seed = 100)
  ph_b  <- simulate_binary_pheno(covars, prevalence = 0.25, seed = 200)

  setwd("C:/Users/callumon/Miami_Package_R/MiamiR/inst/extdata/Fake_GWAS_Files_Raw")

  # Separate phenotype files

  write.table(ph_b, "Fake_PHENOS_Binary.tsv", sep="\t", quote=FALSE, row.names=FALSE)
  write.table(ph_q, "Fake_PHENOS_Quant.tsv", sep="\t", quote=FALSE, row.names=FALSE)

  # Covariates only

  write.table(covars, "Fake_COVARS.tsv", sep="\t", quote=FALSE, row.names=FALSE)

  # Covariates with attached phenotypes

  covars_with_bin   <- merge(covars, ph_b, by=c("FID","IID"))
  covars_with_quant <- merge(covars, ph_q, by=c("FID","IID"))

  setwd("C:/Users/callumon/Miami_Package_R/MiamiR/inst/extdata/Fake_GWAS_Files_Raw")

  write.table(covars_with_bin,   "Fake_COVARS_and_Binary.tsv",   sep="\t", quote=FALSE, row.names=FALSE)
  write.table(covars_with_quant, "Fake_COVARS_and_Quant.tsv", sep="\t", quote=FALSE, row.names=FALSE)

  # Simulate a 100K “fake HLA” GWAS-like dataset for Segregate_HLA_Plot()

  set.seed(1)

  simulate_fake_hla_100k <- function(

      n = 100000,
      chr = 6,
      props = c(HLA = 0.25, AA = 0.25, SNPS = 0.25, RSID = 0.25),
      include_other = FALSE,
      other_prop = 0.0

  ) {

    props <- props / sum(props)

    if (isTRUE(include_other) && other_prop > 0) {

      props <- c(props, OTHER = other_prop)
      props <- props / sum(props)

    }

    classes <- sample(names(props), size = n, replace = TRUE, prob = props)

    # ID gen

    hla_loci   <- c("A","B","C","DRB1","DQA1","DQB1","DPB1","DPA1")

    hla_allele <- sprintf("%02d", 1:40)  # "01".."40"

    make_HLA <- function(k) {

      paste0("HLA_", sample(hla_loci, k, TRUE), "*", sample(hla_allele, k, TRUE))

    }

    aa_cat <- c("A","B","C","DRB1","DQA1","DQB1","DPB1","DPA1")
    aa_aa  <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

    make_AA <- function(k) {

      indel <- sample(c(FALSE, TRUE), k, TRUE, prob = c(0.8, 0.2))
      rest  <- paste0("POS", sample(1:350, k, TRUE), "_", sample(aa_aa, k, TRUE))
      base  <- paste0(sample(aa_cat, k, TRUE), "_", rest)
      ifelse(indel, paste0("INDEL_AA_", base), paste0("AA_", base))

    }

    snps_cat <- c("HLA","MHC","ClassI","ClassII","MHCIII")

    make_SNPS <- function(k) {

      indel <- sample(c(FALSE, TRUE), k, TRUE, prob = c(0.8, 0.2))
      rest  <- paste0("rs", sample(1e6:9e6, k, TRUE), "_",
                      sample(c("A","C","G","T"), k, TRUE), "_",
                      sample(c("A","C","G","T"), k, TRUE))
      base  <- paste0(sample(snps_cat, k, TRUE), "_", rest)

      ifelse(indel, paste0("INDEL_SNPS_", base), paste0("SNPS_", base))

    }

    make_RSID  <- function(k) paste0("rs", sample(1e5:9e8, k, TRUE))
    make_OTHER <- function(k) paste0("XYZ_", sample(LETTERS, k, TRUE), sample(1000:9999, k, TRUE))

    SNP <- character(n)

    k <- sum(classes == "HLA");  if (k) SNP[classes == "HLA"]  <- make_HLA(k)

    k <- sum(classes == "AA");   if (k) SNP[classes == "AA"]   <- make_AA(k)

    k <- sum(classes == "SNPS"); if (k) SNP[classes == "SNPS"] <- make_SNPS(k)

    k <- sum(classes == "RSID"); if (k) SNP[classes == "RSID"] <- make_RSID(k)

    if ("OTHER" %in% classes) {

      k <- sum(classes == "OTHER"); SNP[classes == "OTHER"] <- make_OTHER(k)

    }

    # P Vals

    # 90% uniform null, 10% enriched small P

    is_signal <- runif(n) < 0.10
    P <- numeric(n)

    P[!is_signal] <- runif(sum(!is_signal), min = 0, max = 1)

    # enriched tail: Beta(a=0.25,b=1) gives lots of tiny p-values

    P[is_signal]  <- stats::rbeta(sum(is_signal), shape1 = 0.25, shape2 = 1)

    # Avoid zeros

    P[P < 1e-300] <- 1e-300

    POS <- sample(25e6:34e6, n, TRUE)
    A1  <- sample(c("A","C","G","T"), n, TRUE)
    A2  <- sample(c("A","C","G","T"), n, TRUE)

    # ensure A1 != A2

    same <- A1 == A2

    while (any(same)) {

      A2[same] <- sample(c("A","C","G","T"), sum(same), TRUE)
      same <- A1 == A2
    }

    out <- dplyr::tibble(
      SNP = SNP,
      CHR = rep(chr, n),
      POS = POS,
      A1  = A1,
      A2  = A2,
      P   = P
    )

    out
  }

  Fake_HLA_100K <- simulate_fake_hla_100k(n = 10000, chr = 6)

  # Optional sanity check: does it fall into normal categories

  table(

    dplyr::case_when(

      stringr::str_starts(Fake_HLA_100K$SNP, "HLA_") ~ "HLA",
      stringr::str_starts(Fake_HLA_100K$SNP, "AA_") |
        stringr::str_starts(Fake_HLA_100K$SNP, "INDEL_AA_") ~ "AA",
      stringr::str_starts(Fake_HLA_100K$SNP, "SNPS_") |
        stringr::str_starts(Fake_HLA_100K$SNP, "INDEL_SNPS_") ~ "SNPS",
      stringr::str_detect(Fake_HLA_100K$SNP, "^rs\\d+") ~ "RSID",
      TRUE ~ "OTHER"

    )

  )

  use_data(Fake_HLA_100K)

