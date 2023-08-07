library("future.apply")
library("lava")
library("survival")
library("rate")
library("mets")
library("rate")

progressr::handlers(global = TRUE)
R0 <- 100 #1000

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
  TT1 <- c(unlist(rexp(n, 1) / exp(matrix(c(W, D1, W * D1), ncol = 3) %*% beta)))
  TT0 <- c(unlist(rexp(n, 1) / exp(matrix(c(W, D0, W * D0), ncol = 3) %*% beta)))
  TT <- TT1 * A + TT0 * (1-A)

  # simulate C
  C <- c(unlist(rexp(n, 1) / exp(matrix(c(A, D), ncol = 2) %*% zeta)))

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
  beta = c(0.6, -0.5, -0.3),
  zeta = c(0.2, -0.8),
  kappa = c(2, -0.5),
  tau = 0.5
)

# true values ----------------------------------------------------------

set.seed(1)
d0 <- sim_surv(n.grp = 5e6, beta = par0$beta, zeta = par0$zeta, kappa = par0$kappa)
Psi0_A1 <- mean(((d0$TT1 <= par0$tau)))
Psi0_A0 <- mean((d0$TT0 <= par0$tau))
Psi0_D1 <- mean(d0$D[d0$A == 1])
Psi0 <- mean(((d0$TT1 <= par0$tau) - (d0$TT0 <= par0$tau))[d0$D1 == 1])
rm(d0)

onerun <- function(n.grp,
                   treatment = A ~ 1,
                   post.treatment = D ~ A * W,
                   call.response = "phreg",
                   response = Surv(time, event) ~ W * D,
                   args.response = list(),
                   call.censoring = "phreg",
                   censoring = Surv(time, event == 0) ~ A + D,
                   args.censoring = list(),
                   M = 1
){
  dt <- sim_surv(
    n = n.grp,
    beta = par0$beta,
    zeta = par0$zeta,
    kappa = par0$kappa
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

# true models
future::plan(list(tweak("multisession", workers = 5)))
sim.res.cox <- sim(onerun, R = R0, args = list(n.grp = 500), seed = 1)
future::plan("sequential")
summary(sim.res.cox, estimate = 1:4, se = 5:8, true = c(Psi0_A1, Psi0_A0, Psi0_D1, Psi0))

# random forest models: X = (W,A,D)
future::plan(list(tweak("multisession", workers = 5)))
sim.res.ranger <- sim(onerun, R = R0,
                      args = list(n.grp = 500,
                                  response = Surv(time, event) ~ W + A + D,
                                  call.response = "ranger",
                                  args.response = list(num.threads = 1,
                                                       save.memory = TRUE),
                                  censoring = Surv(time, event == 0) ~ W + A + D,
                                  call.censoring = "ranger",
                                  args.censoring = list(num.threads = 1,
                                                        save.memory = TRUE)),
                      seed = 2)
future::plan("sequential")
summary(sim.res.ranger, estimate = 1:4, se = 5:8, true = c(Psi0_A1, Psi0_A0, Psi0_D1, Psi0))
# 100 * 146 / 60 /  60 ~ 4 hours for 500 replications
saveRDS(sim.res.ranger, file = "sim.res.ranger.rds")

# random forest models: X = (A,D)
future::plan(list(tweak("multisession", workers = 5)))
sim.res.ranger.2 <- sim(onerun, R = R0,
                        args = list(n.grp = 500,
                                    response = Surv(time, event) ~ A + D,
                                    call.response = "ranger",
                                    args.response = list(num.threads = 1,
                                                         save.memory = TRUE),
                                    censoring = Surv(time, event == 0) ~ A + D,
                                    call.censoring = "ranger",
                                    args.censoring = list(num.threads = 1,
                                                          save.memory = TRUE)),
                        seed = 3)
future::plan("sequential")
saveRDS(sim.res.ranger.2, file = "sim.res.ranger.2.rds")

# random forest models: X = (W,A)
future::plan(list(tweak("multisession", workers = 5)))
sim.res.ranger.3 <- sim(onerun, R = R0,
                        args = list(n.grp = 500,
                                    response = Surv(time, event) ~ W + A,
                                    call.response = "ranger",
                                    args.response = list(num.threads = 1,
                                                         save.memory = TRUE),
                                    censoring = Surv(time, event == 0) ~ W + A,
                                    call.censoring = "ranger",
                                    args.censoring = list(num.threads = 1,
                                                          save.memory = TRUE)),
                        seed = 4)
future::plan("sequential")
saveRDS(sim.res.ranger.3, file = "sim.res.ranger.3.rds")

