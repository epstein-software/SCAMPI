################################## DGLM ######################################
#' Double Generalized Linear Models (DGLM) for Multiple Outcomes
#' @description
#' Apply Double Generalized Linear Models (DGLM) to standardize each
#' column in matrix y. The mean sub-model will be standardized by one genotype variable g_var and multiple confounders z_var;
#' The variance sub-model will be standardized by confounders z_var. The confounders z_var is optional.
#'
#' @param y outcome, a matrix with multiple columns, representing multiple traits
#' @param g_var genotype, column matrix G with one column, representing one genotype variable
#' @param z_var optional confounder, matrix Z with multiple columns, representing multiple confounders
#' @returns standardized traits of the original trait matrix y, a matrix with multiple columns
#' @examples
#' DgLm(
#'   y = matrix(rnorm(500), ncol = 5),
#'   z_var = matrix(rnorm(400), ncol = 4),
#'   g_var = rbinom(100, 2, 0.25)
#' )
#' @importFrom dglm dglm
#' @importFrom magrittr %>%
#' @export
DgLm <- function(y, g_var, z_var = NULL) {
  Y_correct <- matrix(, nrow = nrow(y))

  for (col in 1:ncol(y)) {
    # If the confounder exists
    if (!is.null(z_var)) {
      out <- tryCatch(
        dglm(y[, col] ~ factor(g_var) + z_var,
          ~z_var,
          family = stats::gaussian,
          na.action = na.exclude
        ),
        error = function(e) {
          NA
        },
        warning = function(e) {
          NA
        }
      )

      # If the confounder does not exist
    } else {
      out <- tryCatch(
        dglm(y[, col] ~ factor(g_var),
          ~1,
          family = stats::gaussian,
          na.action = na.exclude
        ),
        error = function(e) {
          NA
        },
        warning = function(e) {
          NA
        }
      )
    }

    if (length(out) == 1) {
      if (is.na(out)) {
        scaled <- rep(NA, dim(y)[1])
      }
    } else if (all(!is.na(summary(out)$coefficients))) {
      scaled <- residuals(out) / sqrt(predict.glm(out$dispersion.fit, type = "response"))
    } else {
      scaled <- rep(NA, dim(y)[1])
    }
    Y_correct <- cbind(Y_correct, scaled)
  }

  Y_correct_final <- as.matrix(Y_correct[, -1])

  return(Y_correct_final)
}

########################### SCAMPI ###################################
#' Scalable Cauchy Aggregate test using Multiple Phenotypes to test Interactions (SCAMPI)
#'
#' @description
#' SCAMPI method is applied for detecting the G x E or G x G interaction effects by
#' utilizing both variance and covariance structure of multiple traits.
#' @param y outcome, a matrix with multiple columns, representing multiple traits
#' @param g_var genotype, column matrix G with one column, representing one genotype variable
#' @param z_var optional confounder, matrix Z with multiple columns, representing multiple confounders
#' @returns SCAMPI p-value
#' @examples
#' SCAMPI(
#'   y = matrix(rnorm(500), ncol = 5),
#'   z_var = matrix(rnorm(400), ncol = 4),
#'   g_var = rbinom(100, 2, 0.25)
#' )
#' @export
SCAMPI <- function(y,
                   g_var,
                   z_var = NULL) {
  # 1. Used DGLM to regress out G and Z from the mean effect, and G from the variance effect
  scaled_sim_y <- DgLm(y = y, g_var = g_var, z_var = z_var)

  # 2. Create all pairwise combinations of products
  if (ncol(scaled_sim_y) >= 2) {
    scaled_sim_y_cxp <- apply(
      scaled_sim_y,
      1,
      function(z) {
        c(combn(z, 2, function(z) {
          z[1] * z[2]
        }), z^2)
      }
    )
  } else {
    scaled_sim_y_cxp <- matrix(apply(
      scaled_sim_y,
      1,
      function(z) {
        c(z^2)
      }
    ), nrow = 1)
  }
  scaled_sim_y_cxp_tran <- t(scaled_sim_y_cxp)

  # 3. Regress cross product on genotype

  crox_prod_pvalue <- apply(
    scaled_sim_y_cxp_tran,
    2,
    function(x) {
      lm_alt <- tryCatch(lm(x ~ factor(g_var)),
        error = function(e) {
          NA
        },
        warning = function(e) {
          NA
        }
      )
      lm_nul <- tryCatch(lm(x ~ 1),
        error = function(e) {
          NA
        },
        warning = function(e) {
          NA
        }
      )
      anova_res <- tryCatch(anova(lm_nul, lm_alt),
        error = function(e) {
          NA
        },
        warning = function(e) {
          NA
        }
      )
      tryCatch(anova_res$`Pr(>F)`[2],
        error = function(e) {
          NA
        },
        warning = function(e) {
          NA
        }
      )
    }
  )


  # 4. Use CCT to aggregate the column level p-values and get one final SCAMPI pvalue
  CCT_pvalue <- tryCatch(CCT(crox_prod_pvalue),
    error = function(e) {
      NA
    },
    warning = function(e) {
      NA
    }
  )
  return(list(aggregate_pvalue = CCT_pvalue, single_pvalue = crox_prod_pvalue))
}

################################## CCT #####################################
#' Cauchy Combination Test (CCT)
#' @description
#' CCT is used to aggregate multiple p-value developed by Liu and Xie, 202, JASA, 115:529, 393-402, DOI: 10.1080/01621459.2018.1554485
#' The code used for implementing CCT is adopted from STAAR_v0.9.7 contributed by Xihao Li and Zilin Li at
#' https://github.com/xihaoli/STAAR/blob/master/R/CCT.R
#' @param pvals p-values, an array of p-values.
#' @param weights weights associated with each p-value, SCAMPI used the default weight, 1/len(pvals).
#' @returns a single p-value aggregated from pvals.
#' @examples
#' set.seed(123)
#' pval <- runif(20)
#' CCT(pval)
#' @export
CCT <- function(pvals, weights = NULL) {
  #### check if there is NA
  if (sum(is.na(pvals)) > 0) {
    stop("Cannot have NAs in the p-values!")
  }

  #### check if all p-values are between 0 and 1
  if ((sum(pvals < 0) + sum(pvals > 1)) > 0) {
    stop("All p-values must be between 0 and 1!")
  }

  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals == 0) >= 1)
  is.one <- (sum(pvals == 1) >= 1)
  if (is.zero && is.one) {
    stop("Cannot have both 0 and 1 p-values!")
  }
  if (is.zero) {
    return(0)
  }
  if (is.one) {
    warning("There are p-values that are exactly 1!")
    return(1)
  }

  #### check the validity of weights (default: equal weights) and standardize them.
  if (is.null(weights)) {
    weights <- rep(1 / length(pvals), length(pvals))
  } else if (length(weights) != length(pvals)) {
    stop("The length of weights should be the same as that of the p-values!")
  } else if (sum(weights < 0) > 0) {
    stop("All the weights must be positive!")
  } else {
    weights <- weights / sum(weights)
  }

  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0) {
    cct.stat <- sum(weights * tan((0.5 - pvals) * pi))
  } else {
    cct.stat <- sum((weights[is.small] / pvals[is.small]) / pi)
    cct.stat <- cct.stat + sum(weights[!is.small] * tan((0.5 - pvals[!is.small]) * pi))
  }

  #### check if the test statistic is very large.
  if (cct.stat > 1e+15) {
    pval <- (1 / cct.stat) / pi
  } else {
    pval <- 1 - pcauchy(cct.stat)
  }
  return(pval)
}

########################### Levene test #################################
#' Univariate Levene's Test
#'
#' @description
#' The univariate Levene's test can be applied to a single trait. It is implemented
#' based on the formula in Pare and Cook, 2010, PLoS Genet 6(6): e1000981. doi:10.1371/journal.pgen.1000981
#'
#' @param y outcome, column matrix, representing one single trait
#' @param g_var genotype, column matrix G with one column, representing one genotype variable
#' @returns one single p-value
#' @import dplyr
#' @examples
#' set.seed(123)
#' N <- 1000
#' G <- matrix(rbinom(N, 2, 0.25), ncol = 1)
#' W <- matrix(rnorm(N), ncol = 1)
#' Y <- 0.2 + 0.1 * G + 0.3 * W + 2 * G * W + rnorm(N)
#' levene_test(y = Y, g_var = G)
#' @export
levene_test <- function(y, g_var) {
  if (nrow(y) != nrow(g_var)) {
    stop("y has wrong data slice ")
  }

  n_phe <- dim(y)[2]
  N <- length(g_var)

  N0 <- length(which(g_var == 0))
  N1 <- length(which(g_var == 1))
  N2 <- length(which(g_var == 2))
  N_i <- matrix(c(N0, N1, N2), ncol = 3)

  if (length(which(N_i == 0)) >= 1) {
    zero_mark <- which(N_i == 0)
    N_i <- N_i[, -zero_mark]
  }

  K <- length(unique(g_var))
  W_P <- array(NA, dim = c(n_phe))

  for (w in 1:n_phe) {
    if (dim(y)[1] != length(g_var)) {
      stop("Dimension of the phenotype is different from the length of snp")
    }

    # Column Y and SNP used for the analysis
    temp_df <- data.frame(cbind(y[, w], g_var))
    names(temp_df) <- c("pheno", "g")

    # Average Y per SNP
    Y_bar_g <- temp_df %>%
      group_by(g) %>%
      summarise(Y_bar = mean(pheno)) %>%
      as.data.frame()

    Z <- suppressMessages(temp_df %>%
      left_join(., Y_bar_g) %>%
      mutate(Z = abs(pheno - Y_bar)))

    Z_bar_g <- Z %>%
      group_by(g) %>%
      summarise(Z_bar = mean(Z)) %>%
      arrange(g) %>%
      as.data.frame()

    Z <- suppressMessages(Z %>% left_join(., Z_bar_g))
    Z_bar_overall <- mean(Z$Z)

    # Compute the numerator
    Z_bar_i <- matrix((Z_bar_g$Z_bar - Z_bar_overall)^2, ncol = 1)
    numerator <- (N - K) * (N_i %*% Z_bar_i)

    # Compute the denominator
    denominator <- (K - 1) * sum((Z$Z - Z$Z_bar)^2)

    # Compute W statistics
    W_phe <- numerator / denominator
    W_P_temp <- pf(W_phe, df1 = K - 1, df2 = N - K, lower.tail = FALSE)
    # Compute the p-value for the F statistics
    W_P[w] <- W_P_temp

    rm(
      temp_df, Y_bar_g, Z, Z_bar_g, Z_bar_overall, Z_bar_i, numerator, denominator,
      W_phe, W_P_temp
    )
  }

  return(W_P)
}

########################### Multivariate Levene's Test #################################
#' Multivariate Levene's Test
#'
#' @description
#' Multivariate Levene's Test is extended in our paper from the univariate Levene's Test to
#' accomodate multiple traits.
#' @param y outcome, a matrix with multiple columns, representing multiple traits
#' @param g_var genotype, column matrix G with one column, representing one genotype variable
#' @param z_var optional confounder, matrix Z with multiple columns, representing multiple confounders
#' @returns one aggregated p-value for Multivariate Levene's Test
#' @examples
#' set.seed(123)
#' multivariate_levene_test(
#'   y = matrix(rnorm(500), ncol = 5),
#'   z_var = matrix(rnorm(400), ncol = 4),
#'   g_var = matrix(rbinom(100, 2, 0.25), ncol = 1)
#' )
#' @export
multivariate_levene_test <- function(y, g_var, z_var = NULL) {
  # 1. Used DGLM to regress out G and Z from the mean effect, and G from the variance effect
  scaled_sim_y <- DgLm(y = y, g_var = g_var, z_var = z_var)

  # 2. Compute the Levene's p-value for each column
  single_pvalue <- apply(
    scaled_sim_y,
    2,
    function(x) {
      tryCatch(
        levene_test(
          y = matrix(x, ncol = 1),
          g_var = g_var
        ),
        error = function(e) {
          NA
        },
        warning = function(e) {
          NA
        }
      )
    }
  )

  # 3. Use CCT to aggregate the column level p-values and get one final revised Multivariate Levene's pvalue
  revised_levene_pvalue <- tryCatch(CCT(single_pvalue),
    error = function(e) {
      NA
    },
    warning = function(e) {
      NA
    }
  )
  return(list(aggregate_pvalue = revised_levene_pvalue, single_pvalue = single_pvalue))
}
########################### PC SCAMPI #################################
#' PC SCAMPI
#'
#' @description
#' This is a variation of the SCAMPI method. PC SCAMPI treats the PC of traits as
#' the outcomes.
#' @param y outcome, a matrix with multiple columns, representing multiple traits
#' @param g_var genotype, column matrix G with one column, representing one genotype variable
#' @param z_var optional confounder, matrix Z with multiple columns, representing multiple confounders
#' @returns PC SCAMPI p-value
#' @examples
#' set.seed(123)
#' PC_scampi_test(
#'   y = matrix(rnorm(500), ncol = 5),
#'   z_var = matrix(rnorm(400), ncol = 4),
#'   g_var = matrix(rbinom(100, 2, 0.25), ncol = 1)
#' )
#' @export
PC_scampi_test <- function(y, g_var, z_var = NULL) {
  # 1. Used DGLM to regress out G and Z from the mean effect, and G from the variance effect
  scaled_sim_y <- DgLm(y = y, g_var = g_var, z_var = z_var)

  # 2. Construct the PCs using the data scaled by DGLM
  scaled_sim_y_pca <- tryCatch(prcomp(scaled_sim_y, center = TRUE, scale. = TRUE),
    error = function(e) {
      NA
    },
    warning = function(e) {
      NA
    }
  )
  scaled_sim_y_pca_pred <- tryCatch(predict(scaled_sim_y_pca),
    error = function(e) {
      NA
    },
    warning = function(e) {
      NA
    }
  )

  # 3. Compute the Levene's p-value for each PC
  if (is.null(dim(scaled_sim_y_pca_pred))) {
    single_pca_pvalue <- NA
  } else {
    single_pca_pvalue <- apply(
      scaled_sim_y_pca_pred,
      2,
      function(x) {
        tryCatch(
          levene_test(
            y = matrix(x, ncol = 1),
            g_var = g_var
          ),
          error = function(e) {
            NA
          },
          warning = function(e) {
            NA
          }
        )
      }
    )
  }

  # 4. Use CCT to aggregate the column level p-values and get one final revised Levene's pvalue
  revised_levene_pca_pvalue <- tryCatch(CCT(single_pca_pvalue),
    error = function(e) {
      NA
    },
    warning = function(e) {
      NA
    }
  )
  return(list(aggregate_pvalue = revised_levene_pca_pvalue, single_pvalue = single_pca_pvalue))
}
