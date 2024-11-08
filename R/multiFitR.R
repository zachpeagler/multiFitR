#' Continuous distributions
#'
#' This function returns a list of continuous distributions.
#'
#' @export
cont_distributions <- function() {
  return(c("normal", "lognormal", "gamma", "exponential"))
}

#' Multiple Proportional Density Functions for Continuous Variables
#'
#' This function gets the proportional density functions for selected distributions
#' against continuous, non-negative numbers. Possible distributions include "normal",
#' "lognormal", "gamma", "exponential", "cauchy", "t", "weibull", "logistic",
#' and "all".
#'
#' @param x The variable of which to get the PDF
#' @param seq_length The length of sequence to fit the distribution to
#' @param distributions The distributions to fit x against
#' @returns A dataframe with x, the real density, and the pdf of the desired
#'  distributions with length (nrows) equal to seq_length +1.
#' @export
multiPDF_cont <- function(x, seq_length, distributions){
  # get a sequence from the minimum to maximum of x with length
  #equal to seq_length + 1
  x_seq <- seq(min(x), max(x), length.out = seq_length+1)
  # create real density for x
  x_pdf <- density(x, n=seq_length+1)
  # initialize df of x and the real density
  pdf_df <- as.data.frame(x_seq)
  pdf_df$dens = x_pdf$y
  ## see if "all" is in distributions
  if ("all" %in% distributions) {
    distributions <- c("normal", "lognormal", "gamma", "exponential")
  }

  if ("normal" %in% distributions) {
    x_n <- MASS::fitdistr(x, "normal")
    x_pdf_n <- dnorm(x_seq, mean=x_n$estimate[1],
                     sd = x_n$estimate[2])
    pdf_df$pdf_normal = x_pdf_n
  }
  if ("lognormal" %in% distributions) {
    x_ln <- MASS::fitdistr(x, "lognormal")
    x_pdf_ln <- dlnorm(x_seq, meanlog=x_ln$estimate[1],
                       sdlog = x_ln$estimate[2])
    pdf_df$pdf_lognormal = x_pdf_ln
  }
  if ("gamma" %in% distributions) {
    x_g <- MASS::fitdistr(x, "gamma")
    x_pdf_g <- dgamma(x_seq, shape=x_g$estimate[1],
                      rate=x_g$estimate[2])
    pdf_df$pdf_gamma = x_pdf_g
  }
  if ("exponential" %in% distributions) {
    x_exp <- MASS::fitdistr(x, "exponential")
    x_pdf_exp <- dexp(x_seq, rate = x_exp$estimate)
    pdf_df$pdf_exponential = x_pdf_exp
  }
  ## return dataframe with pdfs
  return(pdf_df)
}

#' Multiple Cumulative Distribution Functions for Continuous Variables
#'
#' This function gets the cumulative distribution function for selected distributions
#' against a continuous, non-negative input variable. Possible distributions include "normal",
#' "lognormal", "gamma", "exponential", "cauchy", "t", "weibull", "logistic",
#' and "all".
#'
#' @param x The variable of which to get the CDF
#' @param seq_length The length of sequence to fit the distribution to
#' @param distributions The distributions to fit x against
#' @returns A dataframe with x, the real density, and the pdf of the desired
#'  distributions with length (nrows) equal to seq_length +1.
#' @export
multiCDF_cont <- function(x, seq_length, distributions){
# get a sequence from the minimum to maximum of x with length
#equal to seq_length + 1
x_seq <- seq(min(x), max(x), length.out = seq_length+1)
# create real cumulative density for x
x_cdf <- ecdf(x)(x_seq)
# initialize df of x and the cumulative density
cdf_df <- as.data.frame(x_seq)
cdf_df$dens = x_cdf
if ("all" %in% distributions) {
  distributions <- c("normal",
                     "lognormal",
                     "gamma",
                     "exponential")
}

if ("normal" %in% distributions) {
  x_n <- MASS::fitdistr(x, "normal")
  x_cdf_n <- pnorm(x_seq, mean=x_n$estimate[1],
                   sd = x_n$estimate[2])
  cdf_df$cdf_normal = x_cdf_n
}
if ("lognormal" %in% distributions) {
  x_ln <- MASS::fitdistr(x, "lognormal")
  x_cdf_ln <- plnorm(x_seq, meanlog=x_ln$estimate[1],
                     sdlog = x_ln$estimate[2])
  cdf_df$cdf_lognormal = x_cdf_ln
}
if ("gamma" %in% distributions) {
  x_g <- MASS::fitdistr(x, "gamma")
  x_cdf_g <- pgamma(x_seq, shape=x_g$estimate[1],
                    rate=x_g$estimate[2])
  cdf_df$cdf_gamma = x_cdf_g
}
if ("exponential" %in% distributions) {
  x_exp <- MASS::fitdistr(x, "exponential")
  x_cdf_exp <- pexp(x_seq, rate = x_exp$estimate)
  cdf_df$cdf_exponential = x_cdf_exp
}

return(cdf_df)
}

#' Multiple Kolmogorov-Smirnov Tests for Continuous Variables
#'
#' This function gets the distance and p-value from a Kolmogorov-smirnov test for selected distributions
#' against a continuous, non-negative input variable. Possible distributions include "normal",
#' "lognormal", "gamma", "exponential", "cauchy", "t", "weibull", "logistic",
#' and "all".
#'
#' @param x The variable to perform KS tests against
#' @param distributions The distributions to test x against
#' @returns A dataframe with the distance and p value for each performed
#' KS test
#' @export
multiKS_cont <- function(x, distributions) {
  # check if "all" was passed to distributions
  if ("all" %in% distributions) {
    distributions <- c("normal",
                       "lognormal",
                       "gamma",
                       "exponential")
  }
  KS_df <- data.frame(matrix(ncol=3, nrow=0))
  colnames(KS_df) <- c("Distribution", "Distance", "P-Value")
  # check normal
  if ("normal" %in% distributions) {
    x_n <- MASS::fitdistr(x, "normal")
    x_KS_n <- ks.test(x, "pnorm", mean=x_n$estimate[1],
                      sd = x_n$estimate[2])
    KS_n <- data.frame(matrix(ncol=0, nrow=1))
    KS_n$Distribution <- "Normal"
    KS_n$Distance <- if (is.null(x_KS_n$statistic)
                         == FALSE) {x_KS_n$statistic}
    else {"NA"}
    KS_n$PValue <- if (is.null(x_KS_n$p.value)
                       == FALSE) {x_KS_n$p.value}
    else {"NA"}
    KS_df <- rbind(KS_df, KS_n)
  }
  if ("lognormal" %in% distributions) {
    x_ln <- MASS::fitdistr(x, "lognormal")
    x_KS_ln <- ks.test(x, "plnorm",
                       meanlog=x_ln$estimate[1],
                       sdlog = x_ln$estimate[2])[c(1, 2)]
    KS_ln <- data.frame(matrix(ncol=0, nrow=1))
    KS_ln$Distribution <- "Lognormal"
    KS_ln$Distance <- if (is.null(x_KS_ln$statistic)
                          == FALSE) {x_KS_ln$statistic}
    else {"NA"}
    KS_ln$PValue <- if (is.null(x_KS_ln$p.value)
                        == FALSE) {x_KS_ln$p.value}
    else {"NA"}
    KS_df <- rbind(KS_df, KS_ln)
  }
  if ("gamma" %in% distributions) {
    x_g <- MASS::fitdistr(x, "gamma")
    x_KS_g <- ks.test(x, "pgamma", shape=x_g$estimate[1],
                      rate=x_g$estimate[2])
    KS_g <- data.frame(matrix(ncol=0, nrow=1))
    KS_g$Distribution <- "Gamma"
    KS_g$Distance <- if (is.null(x_KS_g$statistic)
                         == FALSE) {x_KS_g$statistic}
    else {"NA"}
    KS_g$PValue <- if (is.null(x_KS_g$p.value)
                       == FALSE) {x_KS_g$p.value}
    else {"NA"}
    KS_df <- rbind(KS_df, KS_g)
  }
  if ("exponential" %in% distributions) {
    x_exp <- MASS::fitdistr(x, "exponential")
    x_KS_exp <- ks.test(x, "pexp", rate = x_exp$estimate)
    KS_exp <- data.frame(matrix(ncol=0, nrow=1))
    KS_exp$Distribution <- "Exponential"
    KS_exp$Distance <- if (is.null(x_KS_exp$statistic)
                           == FALSE) {x_KS_exp$statistic}
    else {"NA"}
    KS_exp$PValue <- if (is.null(x_KS_exp$p.value)
                         == FALSE) {x_KS_exp$p.value}
    else {"NA"}
    KS_df <- rbind(KS_df, KS_exp)
  }

  KS_df$Distribution = as.factor(KS_df$Distribution)
  KS_df$Distance = as.numeric(KS_df$Distance)
  KS_df$PValue = as.numeric(format(as.numeric(KS_df$PValue),
                                 scientific = FALSE))
  KS_df$Distance <- round(KS_df$Distance, 3)
  KS_df$PValue <- round(KS_df$PValue, 3)

  return(KS_df)
}