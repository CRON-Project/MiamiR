
#Need non-CRAN ggmanh first

# .onLoad <- function(libname, pkgname) {
#   if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
#   }
#   if (!requireNamespace("ggmanh", quietly = TRUE)) {
#     message("ggmanh package not found. Installing it from Bioconductor...")
#     BiocManager::install("ggmanh")
#   }
# }
#
# .onLoad <- function(libname, pkgname) {
#   if (!requireNamespace("ggmanh", quietly = TRUE)) {
#     message("ggmanh package not found. Installing it from Bioconductor...")
#     if (!requireNamespace("BiocManager", quietly = TRUE)) {
#       install.packages("BiocManager")
#     }
#     BiocManager::install("ggmanh")
#   }
# }
