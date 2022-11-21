library(lava)
library(future.apply)
library(progressr)
library(rate)

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


# True value --------------------------------------------------------------

# tmp <- sim(m, n = 1e7)
# Psi0 <- mean((tmp$y1 - tmp$y0)[tmp$d1 == 1])
# rm(tmp)

# Simulation --------------------------------------------------------------

M <- 5
onerun <- function(n){
  dat <- sim(m, n)

  r <- RATE(response=y~1, post.treatment=d~1, treatment=a~1, data=dat, efficient = FALSE)

  r_eff <- RATE(
    response=y ~ a + w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9  + w10,
    post.treatment = d ~ a + w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9  + w10,
    treatment=a~1,
    data=dat,
    M=M,
    SL.args.response = list(family = binomial(),
                            SL.library = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.gam")),
    SL.args.post.treatment = list(family = binomial(),
                                  SL.library = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.gam"))
  )

  out <- c(
    empir=coef(r)["rate"], empir.se = vcov(r)["rate","rate"]^.5,
    eff=coef(r_eff)["rate"], eff.se=vcov(r_eff)["rate","rate"]^.5
  )

  return(out)
}
# test
# onerun(1e3)

# Simulation
plan("multicore")
handlers(global = TRUE)
R <- 1e3
n <- 1e3
res <- sim(onerun, R = R, args = list(n = n), seed = TRUE)
plan("sequential")
summary(res,estimate=c(1,3),se=c(2,4), true = c(Psi0, Psi0))

save.image(file = "simulation_super.RData")
