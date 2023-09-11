library("future.apply")
library("lava")
library("survival")
library("rate")
library("mets")
library("data.table")

progressr::handlers(global = TRUE)
R0 <- 1000
n.grp0 <- 500

sim_surv <- function(n.grp, beta, zeta, kappa){
  n <- 2 * n.grp

  # id
  id <- 1:n

  # covariate W
  W <- runif(n, min = 1, max = 3)

  # treatment
  A <- c(rep(0, n.grp), rep(1, n.grp))

  # observed exposure
  pd <- lava::expit(kappa[1] + kappa[2] * W)
  D1 <- rbinom(n,1,pd)
  D0 <- rep(0, n)
  D <- D1 * A + D0 * (1-A)

  # simulate T
  TT1 <- c(unlist(rexp(n, 1) / exp(matrix(c(rep(1,n),W, D1, W * D1), ncol = 4) %*% beta)))
  TT0 <- c(unlist(rexp(n, 1) / exp(matrix(c(rep(1,n),W, D0, W * D0), ncol = 4) %*% beta)))
  TT <- TT1 * A + TT0 * (1-A)

  # simulate C
  C <- c(unlist(rexp(n, 1) / exp(matrix(c(rep(1,n), A, D), ncol = 3) %*% zeta)))

  time <- apply(cbind(TT, C), 1, min)
  event <- TT < C

  d <- data.frame(
    id = id,
    W = W,
    A = A,
    D1 = D1,
    D0 = D0,
    D = D,
    TT1 = TT1,
    TT0 = TT0,
    TT = TT,
    C = C,
    time = time,
    event = event
  )
  d <- d[order(time), ]

  return(d)
}

par0 <- list(
  beta = c(-2, 2, -0.2, -0.4),
  zeta = c(1, 1, -1),
  kappa = c(2, -0.5),
  tau = 0.5
)

# true causal target parameter ----------------------------------------------------------

set.seed(1)
d0 <- sim_surv(n.grp = 5e6, beta = par0$beta, zeta = par0$zeta, kappa = par0$kappa)
Psi0_A1 <- mean(((d0$TT1 <= par0$tau)))
Psi0_A0 <- mean((d0$TT0 <= par0$tau))
Psi0_D1 <- mean(d0$D[d0$A == 1])
Psi0 <- mean(((d0$TT1 <= par0$tau) - (d0$TT0 <= par0$tau))[d0$D1 == 1])
rm(d0)


# onerun ------------------------------------------------------------------

onerun <- function(n.grp,
                   treatment = A ~ 1,
                   post.treatment = D ~ A * W,
                   call.response = "phreg",
                   response = Surv(time, event) ~ W * D,
                   args.response = list(),
                   call.censoring = "phreg",
                   censoring = Surv(time, event == 0) ~ A + D,
                   args.censoring = list(),
                   M = 1,
                   return.data = FALSE
){
  dt <- sim_surv(
    n = n.grp,
    beta = par0$beta,
    zeta = par0$zeta,
    kappa = par0$kappa
  )
  if (return.data)
    return(dt
           )

  require("rate")
  require("SuperLearner")
  require("mets")
  require("ranger")
  est <- RATE.surv(
    treatment = treatment,
    post.treatment = post.treatment,
    SL.args.post.treatment = list(family = binomial(),
                                  SL.library = c("SL.glm"),
                                  cvControl = SuperLearner.CV.control(V = 2L)),
    response = response,
    call.response = call.response,
    args.response = args.response,
    censoring = censoring,
    call.censoring = call.censoring,
    args.censoring = args.censoring,
    tau = par0$tau,
    data = dt,
    M = M
  )

  out <- c(
    est$coef,
    setNames(diag(est$vcov)^(0.5), paste(names(est$coef), ".se", sep = ""))
  )
  return(out)
}
# set.seed(1)
# onerun(1e3)

# diagnostics --------------------------------------------------------------------

# tmp <- onerun(n.grp = n.grp0, return.data = TRUE)
# sf <- survfit(Surv(time, event) ~ A, data = tmp)
# plot(sf, mark.time = TRUE)
# table(tmp$event)
#
# plot(x = x = tmp$TT, y = tmp$C)
# abline(a = 0, b = 1)

# simulation --------------------------------------------------------------

start <- Sys.time()

# true cox models:
future::plan(list(tweak("multisession", workers = 6)))
sim.res.cox <- sim(onerun, R = R0, args = list(n.grp = n.grp0), seed = 1)
future::plan("sequential")
summary(sim.res.cox, estimate = 1:4, se = 5:8, true = c(Psi0_A1, Psi0_A0, Psi0_D1, Psi0))

# simple cox models:
future::plan(list(tweak("multisession", workers = 6)))
sim.res.cox.simple <- sim(onerun, R = R0,
                           args = list(n.grp = n.grp0,
                                       post.treatment = D ~ A,
                                       response = Surv(time, event) ~ D,
                                       censoring = Surv(time, event == 0) ~ A+D),
                           seed = 1)
future::plan("sequential")
summary(sim.res.cox.simple, estimate = 1:4, se = 5:8, true = c(Psi0_A1, Psi0_A0, Psi0_D1, Psi0))

# biased cox models:
future::plan(list(tweak("multisession", workers = 6)))
sim.res.cox.bias <- sim(onerun, R = R0,
                        args = list(n.grp = n.grp0,
                                    response = Surv(time, event) ~ W * A,
                                    censoring = Surv(time, event == 0) ~ W * A),
                        seed = 1)
future::plan("sequential")
summary(sim.res.cox.bias, estimate = 1:4, se = 5:8, true = c(Psi0_A1, Psi0_A0, Psi0_D1, Psi0))

end <- Sys.time()


# save --------------------------------------------------------------------

end <- Sys.time()

save.image(file = "rate_surv_2.RData")

# latex output ------------------------------------------------------------

idx <- c(1,2,10, 11, 12, 14)
res.tab <- rbind(
  summary(sim.res.cox, estimate = 4, se = 8, true = Psi0)[idx, ],
  summary(sim.res.cox.simple, estimate = 4, se = 8, true = Psi0)[idx, ],
  summary(sim.res.cox.bias, estimate = 4, se = 8, true = Psi0)[idx, ]
)

dimnames(res.tab)[[1]] <- c("$X = (W, A, D)$", "$X = (A, D)$", "$X = (W,A)$")

library("xtable")
xtab <- xtable(res.tab, method = "compact", label = "tab:surv", digits = 4)
print(xtab, sanitize.text.function = function(x){x})


# kaplan-meier ------------------------------------------------------------

onerun_km <- function(n.grp,
                           treatment = A ~ 1,
                           post.treatment = D ~ A,
                           call.response = "phreg",
                           response = Surv(time, event) ~ A,
                           args.response = list(),
                           call.censoring = "phreg",
                           censoring = Surv(time, event == 0) ~ 1,
                           args.censoring = list(),
                           M = 1,
                           return.data = FALSE
){
  dt <- sim_surv(
    n = n.grp,
    beta = par0$beta,
    zeta = par0$zeta,
    kappa = par0$kappa
  )
  if (return.data)
    return(dt
    )

  require("rate")
  require("SuperLearner")
  require("mets")
  require("ranger")
  est <- RATE.surv(
    treatment = treatment,
    post.treatment = post.treatment,
    SL.args.post.treatment = list(family = binomial(),
                                  SL.library = c("SL.glm"),
                                  cvControl = SuperLearner.CV.control(V = 2L)),
    response = response,
    call.response = call.response,
    args.response = args.response,
    censoring = censoring,
    call.censoring = call.censoring,
    args.censoring = args.censoring,
    tau = par0$tau,
    data = dt,
    M = M
  )

  km <- survfit(response, data = dt)
  skm <- summary(km)
  pd <- data.table(surv = skm$surv, time = skm$time, strata = skm$strata)
  ss <- 1 - pd[time < par0$tau][,.SD[.N], strata]$surv
  te <- ss[2] - ss[1]

  coef <- c(ss[2], ss[1], est$coef[3], te/est$coef[3])
  names(coef) <- names(est$coef)
  out <- c(
    coef,
    setNames(diag(est$vcov)^(0.5), paste(names(est$coef), ".se", sep = ""))
  )
  return(out)
}
# onerun_km(n.grp = 500)

future::plan(list(tweak("multisession", workers = 6)))
sim.res.km <- sim(onerun_km, R = R0, args = list(n.grp = n.grp0), seed = 1)
future::plan("sequential")
summary(sim.res.km, estimate = 1:4, se = 5:8, true = c(Psi0_A1, Psi0_A0, Psi0_D1, Psi0))

idx <- c(1,2,10, 11, 12, 14)
res.tab <- rbind(
  summary(sim.res.km, estimate = 4, se = 8, true = Psi0)[idx, ]
)

library("xtable")
xtab <- xtable(res.tab, method = "compact", label = "tab:surv", digits = 4)
print(xtab, sanitize.text.function = function(x){x})
