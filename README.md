# multiFitR
[![License: MIT](https://img.shields.io/badge/License-MIT-lightgrey.svg)](https://opensource.org/license/mit)
![experimental](https://img.shields.io/badge/lifecycle-experimental-orange)
![year](https://img.shields.io/badge/year-2024-blue)

## Description

An R package that helps with fitting multiple proportional density functions (PDFs), cumulative distribution functions (CDFs), and Kolmogorov-Smirnov tests at a time.

## Functions

### cont_distributions()

#### Description
This function returns a list of continuous distributions: "normal",
"lognormal", "gamma", "exponential", "cauchy", "t", "weibull", and "logistic"
#### Usage
> cont_distributions()

### eligible_distributions(x)
---
#### Description
A function that returns the eligible distributions for a given input variable.

#### Usage
> eligible_distributions(x)

#### Arguments
- **x** - The variable for which to find distribution eligibility.

### multiPDF_cont()
---
#### Description
This function gets the proportional density functions for selected distributions
against continuous, non-negative numbers. Possible distributions include "normal",
"lognormal", "gamma", "exponential", "cauchy", "t", "weibull", "logistic",
and "all".
#### Usage
> multiPDF_cont(x, seq_length, distributions)
#### Arguments
- **x** - The variable of which to fit the PDF
- **seq_length** - The length + 1 to which to set x and fit the distribution
- **distributions** - The distributions to fit "x" against. See cont_distributions for options.

### multiCDF_cont()
---
#### Description
This function gets the cumulative distribution function for selected distributions
against a continuous, non-negative input variable. Possible distributions include "normal",
"lognormal", "gamma", "exponential", "cauchy", "t", "weibull", "logistic",
and "all".
#### Usage
> multiCDF_cont(x, seq_length, distributions)
#### Arguments

- **x** - The variable of which to fit the CDF
- **seq_length** - The length + 1 to which to set x and fit the distribution
- **distributions** - The distributions to fit "x" against. See cont_distributions for options.


### multiKS_cont
---
#### Description
This function gets the distance and p-value from a Kolmogorov-smirnov test for selected distributions
against a continuous, non-negative input variable. Possible distributions include "normal",
"lognormal", "gamma", "exponential", "cauchy", "t", "weibull", "logistic",
and "all".
#### Usage
> multiKS_cont(x, distributions)
#### Arguments
- **x** - The variable to perform KS tests on.
- **distributions** - The distributions to fit "x" against. See cont_distributions for options.