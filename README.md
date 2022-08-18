# rate

Calculate the responder average treatment effect for binary, continuous and time-to-event outcomes subject to right censoring.

## Installation

A recent version of lava will be required:

``` r
# install.packages("devtools")
devtools::install_github("kkholst/lava")
```

You can install the development version of rate from [GitHub](https://github.com/) with:

``` r
devtools::install_github("Andreas-Nordland/rate")
```

## Example: Binary Outcome

``` r
library(rate)
library(lava)
```

Simulation with binary outcome:

``` r
# structural model, see ?lava::lvm
m <- lvm() |>
  addvar("a") |>
  distribution(~a, ones.lvm(p = 0.5)) |>
  addvar(var = paste("w", 1:10, sep = "")) |>
  regression("d1", value = function(w1, w2, w3) (w1>0) * 2*sin(2*w2) + exp(w3)) |>
  distribution(~ d1, binomial.lvm()) |>
  intervention("d",  value=function(a, d1) a * d1) |>
  regression("y1", value = function(d1, w1,w4,w5) d1 * (2 * cos(2*w4) - 1) + w1*w5 + log(abs(w5*w4))) |>
  distribution(~ y1, binomial.lvm()) |>
  regression("y0", value = function(w1,w4,w5) w1*w5 + log(abs(w5*w4))) |>
  distribution(~ y0, binomial.lvm()) |>
  intervention("y",  value=function(y1, y0, a)  a * y1 + (1-a)*y0)
```

``` r
n <- 1e3
dat <- sim(m, n)

rate.est <- RATE(
    response=y ~ a + w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9  + w10,
    post.treatment = d ~ a + w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9  + w10,
    treatment=a~1,
    data=dat,
    M=5, # number of folds for cross-fitting
    SL.args.response = list(family = binomial(),
                            SL.library = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.gam")),
    SL.args.post.treatment = list(family = binomial(),
                                  SL.library = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.gam"))
  )
#> Loading required namespace: SuperLearner
#> Loading required namespace: ranger

rate.est
#>      Estimate Std.Err    2.5%    97.5%   P-value
#> a1     0.2254 0.01795  0.1902  0.26060 3.503e-36
#> a0     0.3396 0.02001  0.3003  0.37877 1.367e-64
#> d      0.7319 0.01924  0.6942  0.76963 0.000e+00
#> rate  -0.1559 0.03455 -0.2237 -0.08821 6.400e-06
```

## Simulations

Simulations can be replicated using the following scripts in inst/sim/

``` r
# binary outcome
source(system.file("sim", "rate.R", package="rate"))
# survival outcome
source(system.file("sim", "rate_surv.R", package="rate"))
```
