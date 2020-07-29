suppressPackageStartupMessages({
  require(data.table)
})

.args <- c("sampled_threshold_grid.rds")
.args <- commandArgs(trailingOnly = TRUE)

rg <- .95
ps <- c(0, 1) + (1-rg)/2*c(1,-1)
qs <- c(1.4, 3.9)
qmatcher <- function(ps, qs, qdistro, init) optim(init, function(xpars) {
  return(sum((do.call(qdistro, args = c(list(p=ps), as.list(xpars)))-qs)^2))
})
 
distropars <- qmatcher(ps, qs, qgamma, c(shape=1, scale=1))$par
# consistency check
# do.call(qgamma, args=c(list(p=c(.25, .75)), distropars1))

ref  <- data.table(expand.grid(
  k = c(0.16, 0.54, 2),
  eta = 1-c(0, .35, .46),#c(0, .42, .46),
  theta = c(0, seq(0.3, 0.7, by=0.2))
))[!(eta == 0 & theta == 0)]

samples <- 1e5

outbreak_constraint <- function(p, R, k) { (1+(R/k)*p)^(-k)-1+p }

outbreak_p <- Vectorize(function(R, k) { 
  uniroot(outbreak_constraint, c(1e-6, 1-1e-6), R=R, k=k)$root
})

lvls <- c(lo.lo=0.025, lo=0.25, md=0.5, hi=0.75, hi.hi=0.975)

ret <- ref[,{
  set.seed(42)
  Rs <- sort(do.call(rgamma, c(list(n=samples), distropars)))
  redRs <- Rs*(1-theta)
  ind <- which.min(redRs <= 1)
  qprime <- ifelse(ind == 1, outbreak_p(redRs, k), outbreak_p(redRs[-c(1:(ind-1))], k))
  qb <- ifelse(ind == 1, outbreak_p(Rs, k), outbreak_p(Rs[-c(1:(ind-1))], k))
  T0b <- log(2)/log((1-qb)^-1)
  T0prime <- log(2)/log((1-qprime)^-1)
  qs <- runif(samples-ind+1)
  res <- quantile(
    c(qgamma(qs, T0prime, rate = 1-eta) - qgamma(qs, T0b, rate=1), rep(Inf,ind-1)),
    probs = lvls
  )
  names(res) <- names(lvls)
  as.list(res)
}, by=1:nrow(ref)]

saveRDS(ret, tail(.args, 1))
