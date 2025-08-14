# --- Minimal GWAS sim: covars + binary / quantitative phenotypes ---

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
  data.frame(FID = covars$FID, IID = covars$IID, PHENO_QUANT = as.numeric(scale(y, center=TRUE, scale=FALSE)))
}

simulate_binary_pheno <- function(
    covars,
    betas = c(Age=0.015, Sex=-0.15, BMI=0.08, PC1=0.20, PC2=-0.05, LocationIE=0.25, LocationPL=-0.10),
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

#Run once; change n/seed as needed
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
