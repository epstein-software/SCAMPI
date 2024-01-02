# This gds data set and workflow are borrowed from SISG, 2021, University of
# Washington, Seattle.
# https://uw-gac.github.io/SISG_2021/gds-format.html

options(shiny.maxRequestSize = 10000 * 1024^2)

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!require("SeqArray", quietly = TRUE)) {
  BiocManager::install("SeqArray")
}

library(SeqArray)

if (!require("SeqVarTools", quietly = TRUE)) {
  BiocManager::install("SeqVarTools")
}
library(SeqVarTools)

library(curl)

source_data_repo_path <- "https://github.com/UW-GAC/SISG_2021/raw/master/data/"
local_data_repo_path <- "data-raw/.internal_data/"

if (!dir.exists(local_data_repo_path)) {
  dir.create(local_data_repo_path)
}

source_data_file <- paste(source_data_repo_path, "1KG_phase3_subset_chr1.vcf.gz", sep = "")
local_data_file <- paste(local_data_repo_path, "1KG_phase3_subset_chr1.vcf.gz", sep = "")

if (!file.exists(local_data_file)) {
  curl::curl_download(source_data_file, local_data_file)
}

# convert the VCF to GDS
gdsfile <- paste(local_data_repo_path, "1KG_phase3_subset_chr1.gds", sep = "")
# seqVCF2GDS(local_data_file,
#            gdsfile,
#            fmt.import="GT",
#            storage.option="LZMA_RA")


# open a connection to the GDS file
gds <- seqOpen(gdsfile)
variantInfo(gds)
usethis::use_data(gds, overwrite = TRUE)
seqClose(gds)

# file.copy(gdsfile, "data/")
