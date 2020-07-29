suppressPackageStartupMessages({
  require(data.table)
})

.args <- c("threshold_grid.rds")
.args <- commandArgs(trailingOnly = TRUE)

ref  <- data.table(expand.grid(
  R0 = 2.5, #seq(1.5, 3, by=0.5),
  k = c(0.16, 0.54, 2),
  eta = c(0, 1-.46, 1-.35), # actual values from input simulation
  theta = c(0, seq(0.3, 0.7, by=0.2))
))[!(eta == 0 & theta == 0)]

outbreak_constraint <- function(p, R, k) { (1+(R/k)*p)^(-k)-1+p }

outbreak_p <- Vectorize(function(R, k) { 
  uniroot(outbreak_constraint, c(1e-6, 1-1e-6), R=R, k=k)$root
})

Tfun <- function(R0, k, c = 0.5){
  # 1-c the probability of outbreak
  # c the probability of control
  -log(c)/log(R0)*(0.334 + 0.689/k + 0.408/R0 - 
                     0.507/(k*R0) - 0.356/(R0^2) + 0.467/(k*R0^2))
}

ref[, qb := outbreak_p(R0, k) ]
ref[R0*(1-theta) > 1, qprime := outbreak_p(R0*(1-theta),k) ]
ref[R0*(1-theta) <= 1, qprime := 0 ]
ref[, T0b := log(2)/log((1-qb)^-1) ]
ref[, Tf := Tfun(R0, k) ]
ref[, T0prime := log(2)/log((1-qprime)^-1) ]
ref[, Tfp := Tfun(R0*(1-theta), k) ]

ref[is.finite(T0prime)][which.max(T0prime)]

ggplot(ref) + 
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_point(aes(x=T0prime, y=Tfp, color="intervention"), data=function(dt) dt[!is.infinite(T0prime)]) +
  geom_point(aes(x=T0b, y=Tf, color="baseline")) +
  scale_x_continuous("T0 via root finding") + scale_y_continuous("T0 via approximation") +
  scale_color_discrete("R0 Calculation") +
  theme_minimal() + theme(
    legend.position = c(0.2, 0.8)
  )

lvls <- c(lo.lo=0.025, lo=0.25, md=0.5, hi=0.75, hi.hi=0.975)

resolution <- 1e5
q <- runif(resolution)

ref[, names(lvls) := Inf ]

ref[!is.infinite(T0prime), names(lvls) := {
  qs <- quantile(qgamma(q, Tfp, rate = 1-eta) - qgamma(q, Tf, rate=1), probs = lvls)
  names(qs) <- names(lvls)
  as.list(qs)
}, by = .(R0, k, eta, theta)]

saveRDS(ref, tail(.args, 1))
