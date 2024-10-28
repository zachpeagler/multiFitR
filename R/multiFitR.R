#' Continuous distributions
#'
#' This function returns a list of continuous distributions.
#'
#' @export
cont_distributions <- function() {
  return(c("normal", "lognormal", "gamma", "exponential",
           "cauchy", "t", "weibull", "logistic"))
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
    distributions <- c("normal", "lognormal", "gamma", "exponential",
                       "cauchy", "t", "weibull", "logistic")
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
  if ("cauchy" %in% distributions) {
    x_cau <- MASS::fitdistr(x, "cauchy")
    x_pdf_cau <- dcauchy(x_seq, location=x_cau$estimate[1],
                         scale = x_cau$estimate[2])
    pdf_df$pdf_cauchy = x_pdf_cau
  }
  if ("t" %in% distributions) {
    x_t <- MASS::fitdistr(x, "t")
    x_pdf_t <- dt(x_seq, df = x_t$estimate[3])
    pdf_df$pdf_t = x_pdf_t
  }
  if ("weibull" %in% distributions) {
    x_wei <- MASS::fitdistr(x, "weibull")
    x_pdf_wei <- dweibull(x_seq, shape = x_wei$estimate[1],
                          scale = x_wei$estimate[2])
    pdf_df$pdf_weibull = x_pdf_wei
  }
  if ("logistic" %in% distributions) {
    x_logis <- MASS::fitdistr(x, "logistic")
    x_pdf_logis <- dlogis(x_seq, x_logis$estimate[1],
                          x_logis$estimate[2])
    pdf_df$pdf_logistic = x_pdf_logis
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
                     "exponential",
                     "cauchy",
                     "t",
                     "weibull",
                     "logistic")
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
if ("cauchy" %in% distributions) {
  x_cau <- MASS::fitdistr(x, "cauchy")
  x_cdf_cau <- pcauchy(x_seq,
                       location = x_cau$estimate[1],
                       scale = x_cau$estimate[2])
  cdf_df$cdf_cauchy = x_cdf_cau
}
if ("t" %in% distributions) {
  x_t <- MASS::fitdistr(x, "t")
  x_cdf_t <- pt(x_seq, df = x_t$estimate[3])
  cdf_df$cdf_t = x_cdf_t
}
if ("weibull" %in% distributions) {
  x_wei <- MASS::fitdistr(x, "weibull")
  x_cdf_wei <- pweibull(x_seq, shape = x_wei$estimate[1],
                        scale = x_wei$estimate[2])
  cdf_df$cdf_weibull = x_cdf_wei
}
if ("logistic" %in% distributions) {
  x_logis <- MASS::fitdistr(x, "logistic")
  x_cdf_logis <- plogis(x_seq, x_logis$estimate[1],
                        x_logis$estimate[2])
  cdf_df$cdf_logistic = x_cdf_logis
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
                       "exponential",
                       "cauchy",
                       "t",
                       "weibull",
                       "logistic")
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
  if ("cauchy" %in% distributions) {
    x_cau <- MASS::fitdistr(x, "cauchy")
    x_KS_cau <- ks.test(x, "pcauchy",
                        location = x_cau$estimate[1],
                        scale = x_cau$estimate[2])
    KS_cau <- data.frame(matrix(ncol=0, nrow=1))
    KS_cau$Distribution <- "Cauchy"
    KS_cau$Distance <- if (is.null(x_KS_cau$statistic)
                           == FALSE) {x_KS_cau$statistic}
    else {"NA"}
    KS_cau$PValue <- if (is.null(x_KS_cau$p.value)
                         == FALSE) {x_KS_cau$p.value}
    else {"NA"}
    KS_df <- rbind(KS_df, KS_cau)
  }
  if ("t" %in% distributions) {
    x_t <- MASS::fitdistr(x, "t")
    x_KS_t <- ks.test(x, "pt", df = x_t$estimate[3])
    KS_t <- data.frame(matrix(ncol=0, nrow=1))
    KS_t$Distribution <- "t"
    KS_t$Distance <- if (is.null(x_KS_t$statistic)
                         == FALSE) {x_KS_t$statistic}
    else {"NA"}
    KS_t$PValue <- if (is.null(x_KS_t$p.value)
                       == FALSE) {x_KS_t$p.value}
    else {"NA"}
    KS_df <- rbind(KS_df, KS_t)
  }
  if ("weibull" %in% distributions) {
    x_wei <- MASS::fitdistr(x, "weibull")
    x_KS_wei <- ks.test(x, "pweibull",
                        shape = x_wei$estimate[1],
                        scale = x_wei$estimate[2])
    KS_wei <- data.frame(matrix(ncol=0, nrow=1))
    KS_wei$Distribution <- "Weibull"
    KS_wei$Distance <- if (is.null(x_KS_wei$statistic)
                           == FALSE) {x_KS_wei$statistic}
    else {"NA"}
    KS_wei$PValue <- if (is.null(x_KS_wei$p.value)
                         == FALSE) {x_KS_wei$p.value}
    else {"NA"}
    KS_df <- rbind(KS_df, KS_wei)
  }
  if ("logistic" %in% distributions) {
    x_logis <- MASS::fitdistr(x, "logistic")
    x_KS_logis <- ks.test(x, "plogis", x_logis$estimate[1],
                          x_logis$estimate[2])
    KS_logis <- data.frame(matrix(ncol=0, nrow=1))
    KS_logis$Distribution <- "Loglogistic"
    KS_logis$Distance <- if (is.null(x_KS_logis$statistic)
                             == FALSE) {x_KS_logis$statistic}
    else {"NA"}
    KS_logis$PValue <- if (is.null(x_KS_logis$p.value)
                           == FALSE) {x_KS_logis$p.value}
    else {"NA"}
    KS_df <- rbind(KS_df, KS_logis)
  }
  KS_df$Distribution = as.factor(KS_df$Distribution)
  KS_df$Distance = as.numeric(KS_df$Distance)
  KS_df$PValue = as.numeric(format(as.numeric(KS_df$PValue),
                                 scientific = FALSE))
  KS_df$Distance <- round(KS_df$Distance, 3)
  KS_df$PValue <- round(KS_df$PValue, 3)

  return(KS_df)
}

#' Multiple PDF Plot For Continuous Variables
#' Using ggplot2
#'
#' This function returns a plotly output showing the PDFs for selected distributions
#' against a continuous, non-negative input variable. Possible distributions include "normal",
#' "lognormal", "gamma", "exponential", "cauchy", "t", "weibull", "logistic",
#' and "all".
#'
#' @param x The variable to for which to plot PDFs
#' @param seq_length The number of points over which to fit x
#' @param distributions The distributions to fit x against
#' @returns A dataframe with x, the real density, and the pdf of the desired
#'  distributions with length(nrow) equal to seq_length +1.
#' @export
multiPDF_plot <- function (x, seq_length, distributions) {
  # check if "all" was passed to distributions
  if ("all" %in% distributions) {
    distributions <- c("normal",
                       "lognormal",
                       "gamma",
                       "exponential",
                       "cauchy",
                       "t",
                       "weibull",
                       "logistic")
  }
  # calculate PDFs
  data <- multiPDF_cont(x, seq_length, distributions)
  # create plot with real density
  p <- ggplot2::ggplot(data) +
    ggplot2::geom_line(aes(x=x_seq, y=dens, color="Real Density"))+
    ggplot2::xlab("x")+
    ggplot2::ylab("PDF")+
    ggplot2::labs(title=paste("PDF plot for x over selected distributions"))+
    ggplot2::guides(color=guide_legend(title="Distribution"))+
    ggplot2::theme_bw()
  # check for each type of distribution in the distributions, and add it if present
  if ("normal" %in% distributions == TRUE) {
    p <- p + ggplot2::geom_line(aes(x=x_seq, y=pdf_normal, color='Normal'))
  }
  if ("lognormal" %in% distributions == TRUE) {
    p <- p + ggplot2::geom_line(aes(x=x_seq, y=pdf_lognormal, color='Lognormal'))
  }
  if ("gamma" %in% distributions == TRUE) {
    p <- p + ggplot2::geom_line(aes(x=x_seq, y=pdf_gamma, color='Gamma'))
  }
  if ("exponential" %in% distributions == TRUE) {
    p <- p + ggplot2::geom_line(aes(x=x_seq, y=pdf_exponential, color='Exponential'))
  }
  if ("cauchy" %in% distributions == TRUE) {
    p <- p + ggplot2::geom_line(aes(x=x_seq, y=pdf_cauchy, color='Cauchy'))
  }
  if ("t" %in% distributions == TRUE) {
    p <- p + ggplot2::geom_line(aes(x=x_seq, y=pdf_t, color='t'))
  }
  if ("weibull" %in% distributions == TRUE) {
    p <- p + ggplot2::geom_line(aes(x=x_seq, y=pdf_weibull, color='Weibull'))
  }
  if ("logistic" %in% distributions == TRUE) {
    p <- p + ggplot2::geom_line(aes(x=x_seq, y=pdf_logistic, color='Loglogistic'))
  }
  return(p)
}

#' Multiple CDF Plot For Continuous Variables
#' Using ggplot2
#'
#' This function returns a plotly output showing the CDFs for selected distributions
#' against a continuous, non-negative input variable. Possible distributions include "normal",
#' "lognormal", "gamma", "exponential", "cauchy", "t", "weibull", "logistic",
#' and "all".
#'
#' @param x The variable to for which to plot PDFs
#' @param seq_length The number of points over which to fit x
#' @param distributions The distributions to fit x against
#' @returns A dataframe with x, the real density, and the pdf of the desired
#'  distributions with length(nrow) equal to seq_length +1.
#' @export
multiCDF_plot <- function (x, seq_length, distributions) {
  # check if "all" was passed to distributions
  if ("all" %in% distributions) {
    distributions <- c("normal",
                       "lognormal",
                       "gamma",
                       "exponential",
                       "cauchy",
                       "t",
                       "weibull",
                       "logistic")
  }
  # calculate PDFs
  data <- multiCDF_cont(x, seq_length, distributions)
  # create plot with real density
  p <- ggplot2::ggplot(data) +
    ggplot2::geom_line(aes(x=x_seq, y=dens, color="Real Distribution"))+
    ggplot2::xlab("x")+
    ggplot2::ylab("CDF")+
    ggplot2::labs(title=paste("CDF plot for x over selected distributions"))+
    ggplot2::guides(color=guide_legend(title="Distribution"))+
    ggplot2::theme_bw()
  # check for each type of distribution in the distributions, and add it if present
  if ("normal" %in% distributions == TRUE) {
    p <- p + ggplot2::geom_line(aes(x=x_seq, y=cdf_normal, color='Normal'))
  }
  if ("lognormal" %in% distributions == TRUE) {
    p <- p + ggplot2::geom_line(aes(x=x_seq, y=cdf_lognormal, color='Lognormal'))
  }
  if ("gamma" %in% distributions == TRUE) {
    p <- p + ggplot2::geom_line(aes(x=x_seq, y=cdf_gamma, color='Gamma'))
  }
  if ("exponential" %in% distributions == TRUE) {
    p <- p + ggplot2::geom_line(aes(x=x_seq, y=cdf_exponential, color='Exponential'))
  }
  if ("cauchy" %in% distributions == TRUE) {
    p <- p + ggplot2::geom_line(aes(x=x_seq, y=cdf_cauchy, color='Cauchy'))
  }
  if ("t" %in% distributions == TRUE) {
    p <- p + ggplot2::geom_line(aes(x=x_seq, y=cdf_t, color='t'))
  }
  if ("weibull" %in% distributions == TRUE) {
    p <- p + ggplot2::geom_line(aes(x=x_seq, y=cdf_weibull, color='Weibull'))
  }
  if ("logistic" %in% distributions == TRUE) {
    p <- p + ggplot2::geom_line(aes(x=x_seq, y=cdf_logistic, color='Loglogistic'))
  }
  return(p)
}

