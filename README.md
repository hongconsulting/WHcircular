
<!-- README.md is generated from README.Rmd. Please edit that file -->

# WHcircular

<!-- badges: start -->
<!-- badges: end -->

This package re-implements circular-linear regression and other selected
functionality from the `circular` package for improved performance,
particularly by avoiding construction of $O(n^2)$ diagonal matrices.
Additional functions are provided for time-of-day variables in “HH:MM”
format. ‘WH’ code is designed to maximize speed, for use in simulation
studies. As such, input validation is minimal and the user is
responsible for ensuring valid use.

## Installation

You can install the development version of WHcircular from
[GitHub](https://github.com/hongconsulting/WHcircular) with:

``` r
remotes::install_github("hongconsulting/WHcircular")
```

## Example 1

``` r
library(WHcircular)
y <- c("01:00","01:15","01:30","01:45","02:00",
       "23:00","23:15","23:30","23:45","00:00")
X <- c(rep(0, 5), rep(1, 5))
print(HHMM.mean(y))
#> [1] "00:30"
print(HHMM.lm(y, ~X))
#>      [,1]          [,2]                        [,3]     
#> [1,] ""            "Time (95%CI)"              "p"      
#> [2,] "(Intercept)" "01:30 (01:16 to 01:44)"    ""       
#> [3,] "β1"          "−02:00 (−02:18 to −01:41)" "< 0.001"
```

## Example 2

``` r
library(WHcircular)
library(nycflights13)
data <- nycflights13::flights
y <- sprintf("%04d", data$sched_dep_time)
y <- paste0(substr(y, 1, 2), ":", substr(y, 3, 4))
print(HHMM.lm(y = y, X = ~as.factor(data$origin),
              name = c("EWR", "JFK", "LGA")))
#>      [,1]  [,2]                        [,3]     
#> [1,] ""    "Time (95%CI)"              "p"      
#> [2,] "EWR" "13:25 (13:24 to 13:27)"    ""       
#> [3,] "JFK" "01:07 (01:05 to 01:09)"    "< 0.001"
#> [4,] "LGA" "−00:15 (−00:17 to −00:13)" "< 0.001"
```

## Example 3

``` r
library(WHcircular)
y <- c("01:00","01:15","01:30","01:45","02:00",
       "23:00","23:15","23:30","23:45","00:00")
X <- c(rep(0, 5), rep(1, 5))
theta <- WH_HHMM_to_rad(y)
print(WH_rad_to_HHMM(WH_rad_mean(theta)))
#> [1] "00:30"
print(WH_rad_to_HHMM(WH_rad_sd(theta)))
#> [1] "01:04"
fit <- WH_reg_circlinear(theta, as.matrix(X))
print(fit$coef)
#> [1] 117.1153750   0.3926991  -0.2679492
print(fit$SE)
#> [1] 52.26269262  0.03086760  0.02219328
print(WH_rad_to_HHMM((fit$coef[2])))
#> [1] "01:30"
print(WH_rad_to_dHHMM(2 * atan(fit$coef[3])))
#> [1] "−02:00"
print(WH_CordeiroPaulaBotter(theta, X))
#> [1] 0.002683716
```

## Further Reading

1.  Cordeiro, G.M., Paula, G.A. and Botter, D.A., 1994. Improved
    likelihood ratio tests for dispersion models. *International
    Statistical Review*, pp. 257–274.
2.  Fisher, N.I. and Lee, A.J., 1992. Regression models for an angular
    response. *Biometrics*, pp. 665–677.
