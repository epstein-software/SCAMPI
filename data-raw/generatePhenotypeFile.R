# This phenotype data set and workflow are borrowed from SISG, 2021, University of
# Washington, Seattle.
# https://uw-gac.github.io/SISG_2021/gds-format.html

library(tidyverse)
options(shiny.maxRequestSize = 10000 * 1024^2)

repo_path <- "https://github.com/UW-GAC/SISG_2021/raw/master"
local_data_repo_path <- "data-raw/.internal_data/"

if (!dir.exists(local_data_repo_path)) {
  dir.create(local_data_repo_path)
}

source_data_repo_path <- "https://github.com/UW-GAC/SISG_2021/raw/master/data/"
local_data_repo_path <- "data-raw/.internal_data/"

pheno_files <- c("pheno_data_study_1.txt", "pheno_data_study_2.txt", "pheno_data_study_3.txt")
source_data_file <- paste(source_data_repo_path, pheno_files, sep = "")
local_data_file <- paste(local_data_repo_path, pheno_files, sep = "")

for (i in 1:3) {
  if (!file.exists(local_data_file[i])) {
    curl::curl_download(source_data_file[i], local_data_file[i])
  }
}

### Inspect individual study data in R

study_1 <- read.table(local_data_file[1], header = TRUE, sep = "\t", as.is = TRUE)
head(study_1)
dim(study_1)

study_2 <- read.table(local_data_file[2], header = TRUE, sep = "\t", as.is = TRUE)
head(study_2)
dim(study_2)

study_3 <- read.table(local_data_file[3], header = TRUE, sep = "\t", as.is = TRUE)
head(study_3)
dim(study_3)

study_2 <- study_2 %>%
  rename(sex = Sex, age = Age, height = Height)

study_3 <- study_3 %>%
  mutate(height = height * 2.54)

###  Compare study values

study_1$study <- "study_1"
study_2$study <- "study_2"
study_3$study <- "study_3"

all.equal(names(study_1), names(study_2))
all.equal(names(study_1), names(study_3))

phen <- dplyr::bind_rows(study_1, study_2, study_3)

# saveRDS(phen, paste(local_data_repo_path, "phen.rda", sep = ""))
usethis::use_data(phen)
