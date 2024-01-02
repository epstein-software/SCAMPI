# This R is used to generate 6 traits for the 1126 samples randomly. The data
# will only be used for demonstration purposes.

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!require("tidyverse")) {
  install.packages("tidyverse")
}
if (!require("MASS")) {
  install.packages("MASS")
}
if (!require("gdata")) {
  install.packages("gdata")
}
if (!require("CompQuadForm")) {
  install.packages("CompQuadForm")
}
if (!require("dglm")) {
  install.packages("dglm")
}
if (!require("lubridate")) {
  install.packages("lubridate")
}
if (!require("GENESIS")) {
  BiocManager::install("GENESIS")
}
if (!require("SeqArray", quietly = TRUE)) {
  BiocManager::install("SeqArray")
}
if (!require("SNPRelate", quietly = TRUE)) {
  BiocManager::install("SNPRelate")
}
if (!require("mvtnorm", quietly = TRUE)) {
  BiocManager::install("mvtnorm")
}

library(tidyverse)
library(MASS)
library(gdata)
library(CompQuadForm)
library(dglm)
library(lubridate)
library(GENESIS)
library(SeqArray)
library(SNPRelate)
library(mvtnorm)

local_data_repo_path <- "data-raw/.internal_data/"
trait <- readRDS(paste(local_data_repo_path, "phen.rda", sep = ""))

local_data_repo_path <- "data-raw/.internal_data/"
gdsfile <- paste(local_data_repo_path, "1KG_phase3_subset_chr1.gds")
gds <- seqOpen(gdsfile)

# Get the sample.id that will be included
SNV.id <- seqGetData(gds, "variant.id")

### I only want to use 3 SNPs out of the 1120 SNPs to generate the traits
set.seed(123)
target <- sample(1:1120, 3, replace = F) # 415 463 179
Sigma_Collect <-
  list(
    Sigma1 = matrix(c(1, 0.3216, 0.3216, 1), ncol = 2),
    Sigma2 = matrix(c(1, 0.021, 0.021, 1), ncol = 2),
    Sigma3 = matrix(c(1, -0.4856, -0.4856, 1), ncol = 2)
  )

dos <- seqGetData(gds, "$dosage")

set.seed(456)
for (genom in target) {
  Genom <- dos[, genom]

  beta0 <- runif(2, 0, 5)
  beta1 <- rnorm(2, 0, 1)
  beta2 <- c(runif(1, 0, 1), runif(1, -1, 0))
  beta3 <- c(runif(1, 0, 1), runif(1, -1, 0))

  Sigma <- Sigma_Collect[[which(target == genom)]]

  data_with <- data.frame(
    G = double(),
    val1 = double(),
    val2 = double(),
    stringsAsFactors = FALSE
  )

  for (iter in 1:nrow(dos)) {
    G <- Genom[iter]
    W <- rnorm(1, 1, 1)

    ### --- Simulate the data with interaction
    mu_with <- c(
      beta0[1] + beta1[1] * G + beta2[1] * W + beta3[1] * G * W,
      beta0[2] + beta1[2] * G + beta2[2] * W + beta3[2] * G * W
    )

    with_inters <- mvtnorm::rmvnorm(1, mean = mu_with, sigma = Sigma)

    temp_with <- data.frame(
      G = G,
      val1 = with_inters[1],
      val2 = with_inters[2],
      stringsAsFactors = FALSE
    )
    data_with <- rbind(data_with, temp_with)
  }

  names(data_with)[2:3] <- paste(c("trait", "trait"),
    c(2 * which(target == genom) - 1, 2 * which(target == genom)),
    sep = "_"
  )

  trait <- cbind(trait, data_with[, 2:3])
}

seqClose(gds)

# saveRDS(trait, paste(local_data_repo_path, "trait.rda", sep = ""))
# usethis::use_data(trait)
