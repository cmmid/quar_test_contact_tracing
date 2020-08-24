covid_pal <- c("#e66101", "#5e3c99", "#0571b0")

#extrafont::loadfonts()
pdf.options(useDingbats=FALSE)

pre_board_labels <- c("NA" = "None",
                      "1"  = "-1 day",
                      "4"  = "-4 days",
                      "7"  = "-7 days")

released_labels <- c("Released after mandatory isolation" =
                       "None",
                     "Released after first test" =
                       "One",
                     "Released after second test" =
                       "Two")

test_labeller <- function(x){
  mutate(x,
         tests = case_when(!is.na(second_test_delay) ~ "Two",
                           screening                 ~ "One",
                           TRUE                      ~ "Zero"),
         tests = factor(tests, levels = c("Zero", "One", "Two"),
                        ordered = T),
         stringency = factor(stringency,
                             levels = c("low",
                                        "moderate",
                                        "high",
                                        "maximum"),
                             labels = c("Low",
                                        "Mod.",
                                        "High",
                                        "Max."),
                             ordered = T))
}

type_labels <- c("asymptomatic" =
                   "Asymptomatic",
                 "symptomatic" =
                   "Pre-symptomatic")

# is this not common to many scripts?
# move to utils.R
main_scenarios <-
  list(`low` = 
         crossing(released_test = c("Released after first test",
                                    "Released after mandatory isolation")),
       `moderate` = 
         crossing(released_test = c("Released after first test",
                                    "Released after mandatory isolation")),
       `high` = 
         crossing(released_test = "Released after second test"),
       `maximum` = 
         crossing(released_test = c("Released after first test",
                                    "Released after mandatory isolation"))
  ) %>%
  bind_rows(.id = "stringency") %>%
  mutate(stage_released = "Infectious",
         stringency = fct_inorder(stringency)) 



add_pre_board_labels <- function(x){
  #inner_join(input) %>% # might we do this here?
  mutate(x,
         pre_board_screening_label =
           as.factor(pre_board_screening),
         pre_board_screening_label = 
           fct_explicit_na(pre_board_screening_label, "NA"),
         pre_board_screening_label = factor(pre_board_screening_label,
                                            levels = names(pre_board_labels),
                                            labels = pre_board_labels, ordered = T))}

syndromic_sensitivity <- 0.7


probs        <- c(0.025,0.25,0.5,0.75,0.975)
lshtm_greens <- rev(c("#00BF6F","#0d5257"))

mv2gamma <- function(mean, var){
  list(shape = mean^2/var,
       rate  = mean/var,
       scale = var/mean) 
}

gamma2mv <- function(shape, rate=NULL, scale=NULL){
  if (is.null(rate)){
    rate <- 1/scale
  }
  
  list(mean = shape/rate,
       var  = shape/rate^2)
}


time_to_event <- function(n, mean, var){
  if (var > 0){
    parms <- mv2gamma(mean, var)
    return(rgamma(n, shape = parms$shape, rate = parms$rate))
  } else{
    return(rep(mean, n))
  }
}

time_to_event_lnorm <- function(n, meanlog, sdlog){
  return(rlnorm(n, meanlog = meanlog, sdlog = sdlog))
}

gamma.parms.from.quantiles <- function(q, p=c(0.025,0.975),
                                       precision=0.001, derivative.epsilon=1e-3, start.with.normal.approx=T, start=c(1.1, 1.1), plot=F, plot.xlim=numeric(0))
{
  # Version 1.0.2 (December 2012)
  #
  # Function developed by 
  # Lawrence Joseph and Patrick Bélisle
  # Division of Clinical Epidemiology
  # Montreal General Hospital
  # Montreal, Qc, Can
  #
  # patrick.belisle@rimuhc.ca
  # http://www.medicine.mcgill.ca/epidemiology/Joseph/PBelisle/GammaParmsFromQuantiles.html
  #
  # Please refer to our webpage for details on each argument.
  
  f <- function(x, theta){dgamma(x, shape=theta[1], rate=theta[2])}
  F.inv <- function(x, theta){qgamma(x, shape=theta[1], rate=theta[2])}
  f.cum <- function(x, theta){pgamma(x, shape=theta[1], rate=theta[2])}
  f.mode <- function(theta){shape <- theta[1]; rate <- theta[2]; mode <- ifelse(shape>1,(shape-1)/rate,NA); mode}
  theta.from.moments <- function(m, v){shape <- m*m/v; rate <- m/v; c(shape, rate)}
  dens.label <- 'dgamma'
  parms.names <- c('shape', 'rate')
  
  
  if (length(p) != 2) stop("Vector of probabilities p must be of length 2.")
  if (length(q) != 2) stop("Vector of quantiles q must be of length 2.")
  p <- sort(p); q <- sort(q)
  
  #_____________________________________________________________________________________________________
  
  print.area.text <- function(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim, M=30, M0=50)
  {
    par.usr <- par('usr')
    par.din <- par('din')
    
    p.string <- as.character(round(c(0,1) + c(1,-1)*p.check, digits=4))
    str.width <- strwidth(p.string, cex=cex)
    str.height <- strheight("0", cex=cex)
    
    J <- matrix(1, nrow=M0, ncol=1)
    
    x.units.1in <- diff(par.usr[c(1,2)])/par.din[1]
    y.units.1in <- diff(par.usr[c(3,4)])/par.din[2]
    aspect.ratio <- y.units.1in/x.units.1in
    
    # --- left area  -----------------------------------------------------------
    
    scatter.xlim <- c(max(plot.xlim[1], par.usr[1]), q[1])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from=0, to=p[1], length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-1]*diff(tmp.x))/diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points to the right of the mode, if any
    w <- which(x>mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[1]+str.width[1]) <= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y, fromLast=T))
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[1], adj=c(1,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[1], mean(par.usr[c(3,4)]), labels=p.string[1], col='gray', cex=cex, srt=90, adj=c(1,0))
    }
    
    # --- right area  ----------------------------------------------------------
    
    scatter.xlim <- c(q[2], plot.xlim[2])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from=p[2], to=f.cum(plot.xlim[2], theta), length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-length(h)]*diff(tmp.x))/diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points to the left of the mode, if any
    w <- which(x<mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[2]-str.width[2]) >= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y))
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[2], adj=c(0,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[2], mean(par.usr[c(3,4)]), labels=p.string[2], col='gray', cex=cex, srt=-90, adj=c(1,0))
    }
  }
  
  # ......................................................................................................................................
  
  Newton.Raphson <- function(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start)
  {
    Hessian <- matrix(NA, 2, 2)
    
    if (start.with.normal.approx)
    {
      # Probably not a very good universal choice, but proved good in most cases in practice
      m <-  diff(q)/diff(p)*(0.5-p[1]) + q[1]
      v <- (diff(q)/diff(qnorm(p)))^2
      theta <- theta.from.moments(m, v)
    }
    else theta <- start
    
    
    change <- precision + 1
    niter <- 0
    # Newton-Raphson multivariate algorithm
    while (max(abs(change)) > precision)
    {
      Hessian[,1] <- (f.cum(q, theta) - f.cum(q, theta - c(derivative.epsilon, 0))) / derivative.epsilon
      Hessian[,2] <- (f.cum(q, theta) - f.cum(q, theta - c(0, derivative.epsilon))) / derivative.epsilon
      
      f <- f.cum(q, theta) - p
      change <- solve(Hessian) %*% f
      last.theta <- theta
      theta <- last.theta - change
      
      # If we step out of limits, reduce change
      
      if (any(theta<0))
      {
        k <- min(last.theta/change)
        theta <- last.theta - k/2*change
      }
      
      niter <- niter + 1
    }
    
    list(theta=as.vector(theta), niter=niter, last.change=as.vector(change))
  }
  
  # ...............................................................................................................
  
  plot.density <- function(p, q, f, f.cum, F.inv, mode, theta, plot.xlim, dens.label, parms.names, cex)
  {
    if (length(plot.xlim) == 0)
    {
      plot.xlim <- F.inv(c(0, 1), theta)
      
      if (is.infinite(plot.xlim[1]))
      {
        tmp <- min(c(0.001, p[1]/10))
        plot.xlim[1] <- F.inv(tmp, theta)
      }  
      
      if (is.infinite(plot.xlim[2]))
      {
        tmp <- max(c(0.999, 1 - (1-p[2])/10))
        plot.xlim[2] <- F.inv(tmp, theta)
      }
    }
    plot.xlim <- sort(plot.xlim)
    
    
    x <- seq(from=min(plot.xlim), to=max(plot.xlim), length=1000)
    h <- f(x, theta)
    x0 <- x; f0 <- h
    ylab <- paste(c(dens.label, '(x, ', parms.names[1], ' = ', round(theta[1], digits=5), ', ', parms.names[2], ' = ', round(theta[2], digits=5), ')'), collapse='')
    plot(x, h, type='l', ylab=ylab)
    
    # fill in area on the left side of the distribution
    x <- seq(from=plot.xlim[1], to=q[1], length=1000)
    y <- f(x, theta)
    x <- c(x, q[1], plot.xlim[1]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # fill in area on the right side of the distribution
    x <- seq(from=max(plot.xlim), to=q[2], length=1000)
    y <- f(x, theta)
    x <- c(x, q[2], plot.xlim[2]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # draw distrn again
    points(x0, f0, type='l')
    h <- f(q, theta)
    points(rep(q[1], 2), c(0, h[1]), type='l', col='orange')
    points(rep(q[2], 2), c(0, h[2]), type='l', col='orange')
    # place text on both ends areas
    print.area.text(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim)  
    
    xaxp <- par("xaxp")
    x.ticks <- seq(from=xaxp[1], to=xaxp[2], length=xaxp[3]+1)
    q2print <- as.double(setdiff(as.character(q), as.character(x.ticks)))
    
    mtext(q2print, side=1, col='orange', at=q2print, cex=0.6, line=2.1)
    points(q, rep(par('usr')[3]+0.15*par('cxy')[2], 2), pch=17, col='orange')
  }
  
  #________________________________________________________________________________________________________________
  
  
  parms <- Newton.Raphson(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start=start)
  p.check <- f.cum(q, parms$theta)
  
  if (plot) plot.density(p, q, f, f.cum, F.inv, f.mode(parms$theta), parms$theta, plot.xlim, dens.label, parms.names, 0.8)
  
  list(shape=parms$theta[1], rate=parms$theta[2], scale=1/parms$theta[2], last.change=parms$last.change, niter=parms$niter, q=q, p=p, p.check=p.check)
}

beta.parms.from.quantiles <- function(q, p=c(0.025,0.975),
                                      precision=0.001, derivative.epsilon=1e-3, start.with.normal.approx=T, start=c(1, 1), plot=F)
{
  # Version 1.3 (February 2017)
  #
  # Function developed by 
  # Lawrence Joseph and pbelisle
  # Division of Clinical Epidemiology
  # Montreal General Hospital
  # Montreal, Qc, Can
  #
  # patrick.belisle@rimuhc.ca
  # http://www.medicine.mcgill.ca/epidemiology/Joseph/PBelisle/BetaParmsFromQuantiles.html
  #
  # Please refer to our webpage for details on each argument.
  
  f <- function(x, theta){dbeta(x, shape1=theta[1], shape2=theta[2])}
  F.inv <- function(x, theta){qbeta(x, shape1=theta[1], shape2=theta[2])}
  f.cum <- function(x, theta){pbeta(x, shape1=theta[1], shape2=theta[2])}
  f.mode <- function(theta){a <- theta[1]; b <- theta[2]; mode <- ifelse(a>1, (a-1)/(a+b-2), NA); mode}
  theta.from.moments <- function(m, v){a <- m*m*(1-m)/v-m; b <- a*(1/m-1); c(a, b)}
  plot.xlim <- c(0, 1)
  
  dens.label <- 'dbeta'
  parms.names <- c('a', 'b')
  
  if (length(p) != 2) stop("Vector of probabilities p must be of length 2.")
  if (length(q) != 2) stop("Vector of quantiles q must be of length 2.")
  p <- sort(p); q <- sort(q)
  
  #_____________________________________________________________________________________________________
  
  print.area.text <- function(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim, M=30, M0=50)
  {
    par.usr <- par('usr')
    par.din <- par('din')
    
    p.string <- as.character(round(c(0,1) + c(1,-1)*p.check, digits=4))
    str.width <- strwidth(p.string, cex=cex)
    str.height <- strheight("0", cex=cex)
    
    J <- matrix(1, nrow=M0, ncol=1)
    
    x.units.1in <- diff(par.usr[c(1,2)])/par.din[1]
    y.units.1in <- diff(par.usr[c(3,4)])/par.din[2]
    aspect.ratio <- y.units.1in/x.units.1in
    
    # --- left area  -----------------------------------------------------------
    
    scatter.xlim <- c(max(plot.xlim[1], par.usr[1]), q[1])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from=0, to=p[1], length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-1]*diff(tmp.x))/diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points to the right of the mode, if any
    w <- which(x>mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[1]+str.width[1]) <= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y, fromLast=T))
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[1], adj=c(1,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[1], mean(par.usr[c(3,4)]), labels=p.string[1], col='gray', cex=cex, srt=90, adj=c(1,0))
    }
    
    # --- right area  ----------------------------------------------------------
    
    scatter.xlim <- c(q[2], plot.xlim[2])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from=p[2], to=f.cum(plot.xlim[2], theta), length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-length(h)]*diff(tmp.x))/diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points to the left of the mode, if any
    w <- which(x<mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[2]-str.width[2]) >= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y))
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[2], adj=c(0,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[2], mean(par.usr[c(3,4)]), labels=p.string[2], col='gray', cex=cex, srt=-90, adj=c(1,0))
    }
  }
  
  # ......................................................................................................................................
  
  Newton.Raphson <- function(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start)
  {
    Hessian <- matrix(NA, 2, 2)
    
    if (start.with.normal.approx)
    {
      # Probably not a very good universal choice, but proved good in most cases in practice
      m <-  diff(q)/diff(p)*(0.5-p[1]) + q[1]
      v <- (diff(q)/diff(qnorm(p)))^2
      theta <- theta.from.moments(m, v)
    }
    else theta <- start
    
    
    change <- precision + 1
    niter <- 0
    # Newton-Raphson multivariate algorithm
    while (max(abs(change)) > precision)
    {
      Hessian[,1] <- (f.cum(q, theta) - f.cum(q, theta - c(derivative.epsilon, 0))) / derivative.epsilon
      Hessian[,2] <- (f.cum(q, theta) - f.cum(q, theta - c(0, derivative.epsilon))) / derivative.epsilon
      
      f <- f.cum(q, theta) - p
      change <- solve(Hessian) %*% f
      last.theta <- theta
      theta <- last.theta - change
      
      # If we step out of limits, reduce change
      
      if (any(theta<0))
      {
        w <- which(theta<0)
        k <- min(last.theta[w]/change[w])
        theta <- last.theta - k/2*change
      }
      
      niter <- niter + 1
    }
    
    list(theta=as.vector(theta), niter=niter, last.change=as.vector(change))
  }
  
  # ...............................................................................................................
  
  plot.density <- function(p, q, f, f.cum, F.inv, mode, theta, plot.xlim, dens.label, parms.names, cex)
  {
    if (length(plot.xlim) == 0)
    {
      plot.xlim <- F.inv(c(0, 1), theta)
      
      if (is.infinite(plot.xlim[1]))
      {
        tmp <- min(c(0.001, p[1]/10))
        plot.xlim[1] <- F.inv(tmp, theta)
      }  
      
      if (is.infinite(plot.xlim[2]))
      {
        tmp <- max(c(0.999, 1 - (1-p[2])/10))
        plot.xlim[2] <- F.inv(tmp, theta)
      }
    }
    plot.xlim <- sort(plot.xlim)
    
    
    x <- seq(from=min(plot.xlim), to=max(plot.xlim), length=1000)
    h <- f(x, theta)
    x0 <- x; f0 <- h
    ylab <- paste(c(dens.label, '(x, ', parms.names[1], ' = ', round(theta[1], digits=5), ', ', parms.names[2], ' = ', round(theta[2], digits=5), ')'), collapse='')
    plot(x, h, type='l', ylab=ylab)
    
    # fill in area on the left side of the distribution
    x <- seq(from=plot.xlim[1], to=q[1], length=1000)
    y <- f(x, theta)
    x <- c(x, q[1], plot.xlim[1]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # fill in area on the right side of the distribution
    x <- seq(from=max(plot.xlim), to=q[2], length=1000)
    y <- f(x, theta)
    x <- c(x, q[2], plot.xlim[2]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # draw distrn again
    points(x0, f0, type='l')
    h <- f(q, theta)
    points(rep(q[1], 2), c(0, h[1]), type='l', col='orange')
    points(rep(q[2], 2), c(0, h[2]), type='l', col='orange')
    # place text on both ends areas
    print.area.text(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim)  
    
    xaxp <- par("xaxp")
    x.ticks <- seq(from=xaxp[1], to=xaxp[2], length=xaxp[3]+1)
    q2print <- as.double(setdiff(as.character(q), as.character(x.ticks)))
    
    mtext(q2print, side=1, col='orange', at=q2print, cex=0.6, line=2.1)
    points(q, rep(par('usr')[3]+0.15*par('cxy')[2], 2), pch=17, col='orange')
  }
  
  #________________________________________________________________________________________________________________
  
  
  parms <- Newton.Raphson(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start=start)
  p.check <- f.cum(q, parms$theta)
  
  if (plot) plot.density(p, q, f, f.cum, F.inv, f.mode(parms$theta), parms$theta, plot.xlim, dens.label, parms.names, 0.8)
  
  list(a=parms$theta[1], b=parms$theta[2], last.change=parms$last.change, niter=parms$niter, q=q, p=p, p.check=p.check)
}


#Nishiura et al. 2020 serial interval
mean_si=4.7
sd_si=2.9

si <- distcrete::distcrete("lnorm",
                           meanlog=log(mean_si),
                           sdlog=log(sd_si),
                           interval = 1,
                           w = 0)

#He et al. Infectivity profile
infect_shape=97.18750 
infect_rate=3.71875
infect_shift=25.62500

# list of pathogens that may be worth considering as sensitivity

gamma.parms.from.quantiles(q = c(5.1, 11.5),
                           p = c(0.5, 0.975)) %>%
  {list(shape = .$shape, scale = .$scale)} -> inc_parms

pathogen <- list(
  symptomatic = 
    # review paper Byrne et al. (2020) https://doi.org/10.1101/2020.04.25.20079889
    # define T1 as infection to beginning of presymptomatic infectious period
    
    
    append(
      # https://www.acpjournals.org/doi/10.7326/M20-0504
      rriskDistributions::get.lnorm.par(q = c(5.1, 11.5),
                                        p = c(0.5, 0.975),
                                        plot = F,
                                        show.output = F) %>%
        as.list %>%
        set_names(., c("mu_inc", "sigma_inc")),
      
      # Li et al https://www.nejm.org/doi/full/10.1056/nejmoa2001316
      {c(9.1, 14.7)} %>% 
        set_names(., c("mu_inf", "sigma_inf"))),
  
  
  
  asymptomatic = 
    append(
      # https://www.acpjournals.org/doi/10.7326/M20-0504
      rriskDistributions::get.lnorm.par(q = c(5.1, 11.5),
                                        p = c(0.5, 0.975),
                                        plot = F,
                                        show.output = F) %>%
        
        set_names(., c("mu_inc", "sigma_inc")),
      # https://doi.org/10.1101/2020.04.25.20079889
      list(
        mu_inf    =  6,
        sigma_inf = 12))) %>%
  map(~data.frame(.x), .id = "type")

# https://www.medrxiv.org/content/10.1101/2020.04.25.20079103v3
asymp_fraction <- beta.parms.from.quantiles(q = c(0.24,  0.38),
                                            p = c(0.025, 0.975)) %>%
  {list(shape1 = .$a,
        shape2 = .$b)}


gen_screening_draws <- function(x){
  #browser()
  
  n <- nrow(x)
  
  # generate screening random draws for comparison
  x <- mutate(x, 
              screen_1 = runif(n, 0, 1),  # on arrival
              screen_2 = runif(n, 0, 1))  # follow-up
}

# given infection histories above, what proportion of travellers end up being 
# caught at each step in the screening process?

calc_outcomes <- function(x, dat_gam){
  
  # generate required times for screening events
  x <- mutate(x,
              first_test_t        = traced_t + first_test_delay,
              second_test_delay   = second_test_delay,
              second_test_t       = first_test_t + second_test_delay)
  
  # what's the probability of PCR detection at each test time?
  x <- mutate(x,first_test_p = 
                c(predict(object  = dat_gam,
                          type    = "response",
                          newdata = data.frame(day = first_test_t))),
              second_test_p = 
                c(predict(object  = dat_gam,
                          type    = "response",
                          newdata = data.frame(day = second_test_t)))
  ) 
  
  # asymptotic infections have a lower detectability
  x <- mutate_at(x, 
                 .vars = vars(ends_with("test_p")),
                 .funs = ~ifelse(type == "asymptomatic",
                                 0.62 * .,
                                 .))
  # if pre-flight PCR is before infection, it's impossible to be detectable
  # this assumes 100% specificity
  # if we look at non-infected travellers in future this will change
  x <- mutate(x,
              first_test_p = ifelse(first_test_t < 0,
                                    0,
                                    first_test_p))
  
  # make comparisons of random draws to screening sensitivity
  # conduct symptomatic screening as well
  x <-
    mutate(x,
           first_test_label       = detector(pcr = first_test_p,  u = screen_1),
           second_test_label      = detector(pcr = second_test_p, u = screen_2))
}

when_released <- function(x){
  #browser()
  mutate(x, released_test = case_when(
    
    !screening    ~ 
      "Released after mandatory isolation",
    
    screening & !first_test_label & is.na(second_test_label) ~
      "Released after first test",
    
    screening & !first_test_label & !second_test_label     ~
      "Released after second test",
    
    first_test_label                            ~
      "Released after first test + mandatory quarantine",
    
    !first_test_label  & second_test_label      ~
      "Released after second test + mandatory quarantine",
    
    TRUE                                        ~ 
      "Mandatory quarantine"
  ),
  released_t = case_when(
    
    released_test == "Released after mandatory isolation"                   ~
      traced_t + first_test_delay, # BILLY TO CHECK
    
    released_test == "Released after first test"                           ~ 
      first_test_t + results_delay,
    
    released_test == "Released after second test"                           ~ 
      second_test_t + results_delay,
    
    released_test == "Released after first test + mandatory quarantine"     ~ 
      first_test_t  + max_mip,
    
    released_test == "Released after second test + mandatory quarantine"    ~
      second_test_t + max_mip,
    
    released_test == "Mandatory quarantine"                                 ~
      traced_t + max_mip)) %>% 
    mutate(released_test =  ifelse(type == "symptomatic" & 
                                     onset > traced_t &
                                     onset < released_t,
                                   "Symptomatic during quarantine",
                                   released_test),
           released_t = ifelse(released_test == "Symptomatic during quarantine",
                               traced_t + pmax(onset + post_symptom_window,
                                               symp_end, max_mip),
                               released_t))
}

stage_when_released <- function(x){
  mutate(x, stage_released = as.factor(case_when(
    is.na(released_t)         ~ "Prevented from boarding",
    released_t < inf_end      ~ "Infectious",
    released_t >= inf_end     ~ "Post-infectious"
  ))) %>% 
    mutate(
      # how many days infectious after tracing
      days_released_inf = 
        case_when(
          stage_released == "Post-infectious" ~ 0,
          type           == "asymptomatic"    ~ inf_end - pmax(released_t, inf_start),
          type           == "symptomatic"     ~ onset   - pmax(released_t, inf_start)),
      # how many days infectious before tracing
      days_prior_inf =
        case_when(
          traced_t > inf_end   ~ inf_dur,
          traced_t < inf_start ~ 0,
          TRUE                 ~ traced_t - inf_start))
  
}


detector <- function(pcr, u = NULL, spec = 1){
  
  if (is.null(u)){
    u <- runif(n = length(pcr))
  }
  
  # true positive if the PCR exceeds a random uniform
  # when uninfected, PCR will be 0
  TP <- pcr > u
  
  # false positive if in the top (1-spec) proportion of random draws
  FP <- (pcr == 0)*(runif(n = length(pcr)) > spec)
  
  return(TP | FP)
}


make_plot_labels <- function(x){
  mutate(x, country = paste("Origin:",
                            ifelse(country == "United States of America",
                                   "USA", country))) %>%
    mutate_at(.vars = vars(contains("test_delay"),
                           pre_board_screening,
                           max_quarantine),
              .funs = function(x){ 
                sub(pattern = "1 days",
                    x = paste(x, "days"),
                    replacement = "1 day")
              })  %>%
    mutate(
      pre_board_screening = ifelse(pre_board_screening == "NA days",
                                   "No syndromic screening",
                                   paste("Syndromic screening", 
                                         pre_board_screening,
                                         "prior to departure"))) %>%
    mutate(max_quarantine = ifelse(max_quarantine == "0 days",
                                   "No isolation",
                                   paste("Isolate for",
                                         max_quarantine))) %>%
    mutate_at(.vars = vars(pre_board_screening,
                           matches("(first|second)_test_delay")),
              .funs = function(x){ifelse(grepl(pattern = "NA days",
                                               x = x),
                                         NA_character_,
                                         x)}) %>%
    mutate(stringency = stringr::str_to_title(stringency),
           stringency = paste(stringency, "stringency") ) %>%
    unite(col = "label",
          country,
          stringency,
          pre_board_screening,
          max_quarantine,
          sep = "\n") %>%
    dplyr::select(label) %>%
    {bind_cols(x, .)}
}


make_delay_label <- function(x,s){
  paste(na.omit(x), s)
}


capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

make_incubation_times <- function(n_travellers,
                                  pathogen,
                                  asymp_parms){
  #browser()
  incubation_times <- crossing(i  = 1:n_travellers,
                               type = c("symptomatic",
                                        "asymptomatic") %>%
                                 factor(x = .,
                                        levels = .,
                                        ordered = T)) %>%
    mutate(idx=row_number()) %>% 
    select(-i) %>% 
    split(.$type) %>%
    map2_df(.x = .,
            .y = pathogen,
            ~mutate(.x,
                    exp_to_onset   = time_to_event_lnorm(n = n(),
                                                         meanlog = .y$mu_inc, 
                                                         sdlog   = .y$sigma_inc),
                    onset_to_recov = time_to_event(n = n(),
                                                   mean = .y$mu_inf, 
                                                   var  = .y$sigma_inf))) 
  
  source("wolfel.R")
  source("he.R")
  # infectious period from time of onset to no longer infectious
  incubation_times %<>% 
    mutate(u = runif(n = nrow(.), 0.01, 0.99)) %>%
    mutate(#inf_from_onset = 
      #  approx(x    = wolfel_pred$y, 
      #         y    = wolfel_pred$day, 
      #        xout = u)$y,
      pre_symp_lead  = 
        approx(x    = HE$p,
               y    = HE$delay,
               xout = pmin(1 - 1e-5,
                           pmax(1e-5,
                                pgamma(q = exp_to_onset,
                                       shape = inc_parms$shape,
                                       scale = inc_parms$scale))))$y
    )
  
  incubation_times %<>% 
    mutate(onset     = exp_to_onset,
           #inf_start = onset - pre_symp_lead,
           #inf_end   = ifelse(type == "asymptomatic",
           #exp_to_onset + onset_to_recov,
           #exp_to_onset + inf_from_onset),
           symp_end  = ifelse(type == "asymptomatic",
                              onset, # but really never matters because asymptomatics are never symptomatic!
                              exp_to_onset + onset_to_recov),
           #inf_dur   = inf_end - inf_start,
           symp_dur  = symp_end - onset)
  
  incubation_times %<>% gen_screening_draws
  
  return(incubation_times)
  
}


## just making sure the proportion of cases are secondary or not
make_sec_cases <- function(prop_asy, incubation_times){
  #browser()
  
  #browser()
  props <- c("asymptomatic"=prop_asy,
             "symptomatic"=(1-prop_asy))
  
  split_inc <- split(incubation_times,incubation_times$type)
  
  res <- lapply(seq_along(props), function(x) sample_frac(split_inc[[x]],props[[x]]))
  
  res <- do.call("rbind",res)
  res
}

make_arrival_scenarios <- function(input, 
                                   inf_arrivals, 
                                   incubation_times){
  #source('kucirka_fitting.R', local=T)
  
  arrival_scenarios <- crossing(input, inf_arrivals)
  
  # calculate outcomes of screening
  arrival_scenarios %<>% calc_outcomes(., dat_gam)
  
  return(arrival_scenarios)
  
}

index_test_labeller <- function(x, newline = FALSE){
  paste0("Delay between index case's\nonset and having a test:",
         ifelse(newline, "\n", " "),
         x,
         ifelse(x == 1, " day", " days"))
}

delay_scaling_labeller <- function(x, newline = FALSE){
  paste0("Contact tracing delay",
         ifelse(newline, "\n", " "),
         "scaling factor: ",
         x)
}


waning_labeller <- function(x){
  paste("Waning:",
        dplyr::case_when(x == "waning_canada_total" ~ "Exponential decay",
                         x == "waning_constant"     ~ "Constant",
                         TRUE ~ "Unknown"))
}

percentage <- function(x, ...){
  
  ans <- scales::percent(x, ...)
  ans[x > 1] <- ""
  
  ans
  
}

pretty_percentage <- function(x){
  ans <- pretty(x)
  ans[ans <= 1]
}


make_release_figure <- function(x_summaries,
                                input,
                                xlab = "Days in quarantine",
                                ylab = "",
                                text_size = 2.5,
                                text_angle = 45,
                                h_just = 0,
                                log_scale = FALSE,
                                hline = 0,
                                faceting = NULL,
                                percent = FALSE){
  #browser()
  x_summaries %<>% test_labeller
  
  # how to do presymptomatic
  
  facet_vars <- all.vars(faceting)
  
  if ("type" %in% facet_vars){
    x_summaries %<>% mutate(type = factor(type,
                                          levels = c("asymptomatic",
                                                     "symptomatic"),
                                          labels = c("Asymptomatic",
                                                     "Presymptomatic")))
  }
  
  require(formula.tools)
  
  
  dy <- select(x_summaries, !!!facet_vars,
               `97.5%`, `75%`) %>%
    gather(key, value, -c(!!!facet_vars)) %>%
    group_by_at(.vars = vars(-value)) %>%
    filter(value == max(value)) %>%
    distinct %>%
    spread(key, value) %>%
    mutate(stringency = fct_collapse(stringency, 
                                     "High" = "High",
                                     other_level = "Other")) %>%
    group_by_at(.vars = vars(-`75%`, -`97.5%`)) %>%
    summarise_all(.funs = ~max(.)) %>%
    tidyr::pivot_wider(names_from = "stringency", 
                       names_glue = "{stringency}_{.value}",
                       values_from = c(`75%`, `97.5%`))
  
  
  
  needs_expanding <- dy %>% 
    filter(`High_97.5%` > 0.75*`Other_97.5%`) %>%
    {nrow(.) > 0L}
  
  dy %<>% mutate(ypos = 0.05 * pmax(`High_97.5%`, `Other_97.5%`)) %>%
    select(one_of(all.vars(faceting)), ypos)
  
  
  
  x_summaries %<>% left_join(dy)
  
  figure <-  
    ggplot(data=x_summaries, aes(x = time_in_iso, 
                                 y = `50%`, 
                                 color = tests)) +
    geom_hline(aes(yintercept=1),linetype=hline)+
    geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`,
                       group = delays),
                   position = position_dodge2(width = 0.75),
                   alpha = 0.3,
                   size = 3) +
    geom_linerange(aes(ymin = `25%`, ymax = `75%`,
                       group = delays),
                   position = position_dodge2(width = 0.75),
                   alpha = 0.5,
                   size = 3) +
    geom_point(pch = "-", size = 12,
               position = position_dodge2(width = 0.75),
               aes(y = `50%`,
                   group = delays)
    ) +
    geom_text(data=filter(x_summaries, stringency=="High"),
              aes(x     = time_in_iso,
                  y     = ifelse(percent, 1.05, `97.5%` + ypos),
                  label = delays),
              angle     = text_angle,
              hjust     = h_just,
              vjust     = 0,
              size      = text_size,
              position  = position_dodge2(width = 0.75),
              check_overlap = TRUE,
              show.legend = F)+
    scale_color_manual(name = "Number of negative tests required for release",
                       values = covid_pal)+
    theme_minimal()+
    theme(axis.ticks = element_line(),
          panel.grid.major.x = element_blank(),
          panel.border = element_rect(fill=NA),
          legend.position = "bottom",
          strip.placement = "outside") +
    ggplot2::labs(x = xlab,
                  y = ylab) +
    xlab("Days in quarantine\n(including 1 day delay on testing results)")
  
  figure <- figure + facet_nested(
    nest_line = T,
    facets = faceting,
    labeller = labeller(index_test_delay = index_test_labeller,
                        delay_scaling    = delay_scaling_labeller,
                        waning           = waning_labeller),
    scales = "free_x", space = "free")
  # }
  
  # check if the top of the y axis needs adjustment
  
  
  
  # end check top y
  
  if (log_scale) {
    mult <- c(0.1, ifelse(needs_expanding, 0.5, 0.1))
    figure <- figure +
      scale_y_log10(labels=format_format(scientific=FALSE),
                    expand = expansion(mult = mult)) +
      theme(panel.grid.minor.y = element_blank()) +
      annotation_logticks(sides = "l")
    
  } else {
    mult <- c(0.1, ifelse(needs_expanding, 0.25, 0.1))
    figure <- figure + 
      scale_y_continuous(#limits=c(-dy/5,NA),
        breaks = pretty_percentage,
        expand = expansion(mult = mult),
        #limits = ifelse(percent, c(0,1), NULL),
        labels = ifelse(percent, percentage, scales::number))
  }
  
  
  return(figure)
  
}

plot_data <- function(input, 
                      x_summaries,
                      main_scenarios = NULL){
  
  dat <- x_summaries  %>%
    inner_join(input) %>% # should really carry these through when summarising
    mutate(second_test_delay_ = 
             ifelse(is.na(second_test_delay),
                    0,
                    second_test_delay),
           time_in_iso = 
             first_test_delay + 
             second_test_delay_+
             screening)
  
  
  if (!is.null(main_scenarios)){
    main_scenarios %<>% select(-one_of("released_test")) %>% distinct
    dat <- left_join(dat, main_scenarios)
  }
  
  dat %>%  
    #tidyr::nest(data = -c(first_test_delay, second_test_delay)) %>%
    dplyr::mutate(delays = paste(first_test_delay, "&",
                                 first_test_delay + second_test_delay)) %>%
    #tidyr::unnest(data) %>%
    dplyr::mutate(time_in_iso = factor(time_in_iso, 
                                       levels = sort(unique(.$time_in_iso)),
                                       ordered = T)) %>%
    dplyr::filter(M!=0) %>%  # if the mean is zero, this group is empty
    return
}

make_arrivals_table <- function(x, table_vars = c("country")){
  x %>%
    mutate(sim = factor(sim, levels = 1:n_arrival_sims)) %>%
    group_by_at(.vars = vars(sim, one_of(table_vars))) %>%
    count(.drop = F) %>%
    group_by_at(.vars = vars(-n, -sim)) %>%
    nest() %>%
    mutate(Q = map(data, ~quantile(.x$n, probs = probs))) %>%
    unnest_wider(Q) %>%
    select(-data) %>%
    group_by_at(.vars = vars(one_of(table_vars))) %>%
    transmute(value = sprintf("%0.0f (95%%: %0.0f, %0.0f)", `50%`, `2.5%`, `97.5%`)) %>%
    mutate_at(.vars = vars(-c(country, value)), 
              .funs = str_to_title) %>%
    spread(country, value)
}


make_plots <- function(
  x, 
  input,
  main_scenarios = NULL,
  log_scale = FALSE,
  #fixed = TRUE,
  text_size = 2.5,
  #trav_vol_manual = NULL,
  xlab = "Days in quarantine\n(including 1 day delay on testing results)",
  sum = FALSE,
  y_var = "days_released_inf",
  faceting = NULL){
  
  #browser()
  
  ylabA = "Number of infectious persons\nreleased per index case"
  
  if (sum){
    ylabB = 
      "Total number of person-days infectiousness\nremaining for released secondary case"
  } else {
    ylabB = 
      "Number of days infectiousness\nremaining per released secondary case"
  }
  
  
  all_grouping_vars <- all.vars(faceting)
  
  # ... should be country, type, pre_board_screening, etc.
  x_summaries <- x %>%
    filter(stage_released == "Infectious") %>%
    make_released_quantiles(., all_grouping_vars) %>%
    inner_join(input)
  
  
  # why are the plo medians coming out backwars in Low?
  figA_data <- plot_data(input, 
                         x_summaries,
                         main_scenarios)
  
  
  # deal with the join and stage_released # need to have scenario guide the labelling
  
  # browser()
  figA <- figA_data %>% 
    #filter(pre_board_screening == "None") %>% 
    # # this filter should be done outside the function
    make_release_figure(
      x_summaries = .,
      input     = input,
      text_size = text_size,
      xlab      = xlab,
      ylab      = ylabA, 
      faceting  = faceting) 
  
  ## person-days
  
  x_days_summaries <- 
    make_released_time_quantiles(x, 
                                 y_var = y_var,
                                 all_grouping_vars,
                                 sum = sum)
  
  
  
  figB_data <- plot_data(input = input, 
                         x_summaries = 
                           x_days_summaries,
                         main_scenarios)
  
  figB <- figB_data %>% 
    #filter(pre_board_screening == "None") %>%  
    make_release_figure(
      x = .,
      input=input,
      xlab = "Days in quarantine\n(including 1 day delay on testing results)",
      text_size = text_size,
      ylab = ylabB,
      faceting = faceting) 
  
  
  fig <- figA + figB + plot_layout(ncol = 1, guide = "collect") +
    plot_annotation(tag_levels = "A",
                    theme = theme(legend.position = "bottom"))
  
  return(fig)
  
}



make_released_quantiles <- function(x, vars){
  
  dots1 <- rlang::exprs(sim, scenario)
  dots2 <- lapply(vars, as.name)
  
  dots <- append(dots1, dots2)
  
  x_count <- x %>%
    dplyr::ungroup(.) %>%
    dplyr::select(!!! dots) %>%
    dplyr::group_by_all(.) %>%
    dplyr::count(.) 
  
  x_count %>%
    dplyr::ungroup(.) %>%
    select(-n) %>%
    ungroup %>%
    as.list %>%
    map(unique) %>%
    expand.grid %>%
    left_join(x_count) %>%
    mutate(n = ifelse(is.na(n), 0, n)) %>%
    tidyr::nest(data = c(sim, n)) %>%
    dplyr::mutate(
      Q = purrr::map(
        .x = data,
        .f = ~quantile(.x$n, probs = probs)),
      M = purrr::map_dbl(
        .x = data, 
        .f = ~mean(.x$n))) %>%
    tidyr::unnest_wider(Q) %>%
    dplyr::select(-data) %>%
    dplyr::ungroup(.)
}

make_released_time_quantiles <- function(x, y_var, vars, sum = FALSE){
  #browser()
  
  dots1 <- rlang::exprs(sim, scenario)
  dots2 <- lapply(vars, as.name)
  y_var <- as.name(y_var)
  dots  <- append(dots1, dots2)
  
  if (sum){
    x <- x %>%
      dplyr::select(!!! dots, y_var) %>%
      group_by_at(.vars = vars(-y_var)) %>%
      summarise(y_var = sum(y_var, na.rm=T))
  }
  
  x_days <- x %>%
    dplyr::select(!!! dots, !! y_var) %>%
    dplyr::filter( !!y_var > 0)
  
  x_days %>%
    nest(data = c(!!y_var, sim)) %>%
    mutate(Q = purrr::map(.x = data, ~quantile( .x[[y_var]],
                                                probs = probs)),
           M = map_dbl(.x = data, ~mean(.x[[y_var]]))) %>%
    unnest_wider(Q) %>%
    select(-data)
  
}



save_plot <- function(plot   = ggplot2::last_plot(),
                      prefix = stringi::stri_rand_strings(1, length = 8),
                      base   = NULL, # to help identify which figure in the paper
                      device = NULL, # additional devices, e.g. "png" 
                      width  = 210, 
                      height = 210,
                      dpi    = 300,
                      units  = "mm"){
  
  file <- paste0("results/",
                 paste(prefix, base, sep = "_"),
                 ".pdf")
  
  pdf(file = file,
      width  = measurements::conv_unit(x = width,  units, "inch"),
      height = measurements::conv_unit(x = height, units, "inch"),
      useDingbats = FALSE)
  
  print(plot)
  
  dev.off()
  
  if (length(device) > 0L){
    #img <- pdftools::pdf_render_page(file, dpi = dpi)
    
    img <- magick::image_read_pdf(file, density = dpi)
    
    purrr::map(.x = device,
               .f = 
                 ~magick::image_convert(image = img, format = .x) %>%
                 magick::image_write(., path = sub(pattern = "pdf",
                                                   replacement = .x,
                                                   x = file)))
  }
}


delay_to_gamma <- function(x){
  ans <- dplyr::transmute(x, left = t - 0.5, right = t + 0.5) %>%
    dplyr::mutate(right = ifelse(right == max(right), Inf, right)) %>%
    {fitdistcens(censdata = data.frame(.),
                 distr = "gamma", 
                 start = list(shape = 1, rate = 1))} 
  
  return(gamma2mv(ans$estimate[["shape"]],
                  ans$estimate[["rate"]]))
}

run_analysis <- 
  function(n_sims          = 1000,
           n_sec_cases     = 1000, # this shouldn't matter. just needs to be Big Enough
           n_ind_cases     = 10000,
           input,
           seed            = 145,
           P_c, P_r, P_t,
           dat_gam,
           asymp_parms){       # a list with shape parameters for a Beta
    
    #browser()
    
    message(sprintf("\n%s == SCENARIO %d ======", Sys.time(), input$scenario))
    
    #browser()
    set.seed(seed)
    
    message("Generating incubation times")
    
    #Generate incubation periods to sample
    incubation_times <- make_incubation_times(
      n_travellers = n_ind_cases,
      pathogen     = pathogen,
      asymp_parms  = asymp_parms)
    
    message("Generating asymptomatic fractions")
    inf <- data.frame(prop_asy = rbeta(n = n_sims,
                                       shape1 = asymp_parms$shape1,
                                       shape2 = asymp_parms$shape2)) 
    
    
    message("Generating index cases' transmissions")
    # Generate index cases' inc times
    ind_inc <- incubation_times %>% 
      filter(type=="symptomatic") %>% 
      sample_n(n_sims) %>% 
      mutate(sim = seq(1L, n_sims, by = 1L)) %>% 
      bind_cols(inf) %>% 
      #sample test result delay
      ## sample uniformly between 0 and 1 when 0.5...
      mutate(index_result_delay = time_to_event(n = n(),
                                                mean = P_r[["mean"]],
                                                var  = P_r[["var"]])) %>% 
      #sample contact info delay
      mutate(contact_info_delay = time_to_event(n = n(),
                                                mean = P_c[["mean"]],
                                                var  = P_c[["var"]])) %>% 
      #sample tracing delay
      mutate(tracing_delay      = time_to_event(n = n(),
                                                mean = P_t[["mean"]],
                                                var  = P_t[["var"]])) %>% 
      #add index test delay (assume 2 days post onset)
      #crossing(distinct(input, index_test_delay, delay_scaling, waning)) %>%     
      crossing(distinct(input, index_test_delay, delay_scaling, waning)) %>%     
      mutate_at(.vars = vars(tracing_delay, contact_info_delay, index_result_delay),
                .funs = ~(. * delay_scaling)) %>%
      rename("index_onset" = onset) %>% 
      mutate(index_testing_t    = index_onset + index_test_delay,
             index_result_t     = index_onset + index_test_delay + index_result_delay,
             traced_t           = index_onset + index_test_delay + index_result_delay +
               contact_info_delay + tracing_delay)
    
    rm(list = c("P_t", "P_r", "P_c", "inf"))
    
    message("Generating secondary cases' incubation times")
    
    #Generate secondary cases
    sec_cases <- make_incubation_times(
      n_travellers = n_sec_cases,
      pathogen     = pathogen,
      asymp_parms  = asymp_parms)
    
    #browser()
    message("Generating secondary cases' exposure times")
    ind_inc %<>% 
      nest(data = -c(sim,
                     prop_asy,
                     index_onset,
                     index_test_delay,
                     index_result_delay,
                     contact_info_delay,
                     tracing_delay,
                     index_testing_t,
                     traced_t,
                     delay_scaling))
    
    ind_inc %<>% 
      mutate(prop_asy    = as.list(prop_asy)) %>%
      mutate(sec_cases   = map(.x = prop_asy, 
                               .f  = ~make_sec_cases(as.numeric(.x),
                                                     sec_cases)
      ))
    
    rm(sec_cases)
    
    ind_inc %<>%
      unnest(prop_asy) %>%
      unnest(sec_cases) %>% 
      ungroup() 
    
    ind_inc %<>%
      select(-data)
    
    ind_inc %<>% 
      #rowwise %>%
      mutate(exposed_t = index_onset + (
        rtgamma(n     = n(),
                b     = infect_shift + index_testing_t,
                shape = infect_shape,
                rate  = infect_rate) - infect_shift)
      ) #%>% ungroup
    
    
    message("Shifting secondary cases' times relative to index cases' times")
    #exposure date relative to index cases exposure
    incubation_times_out <- ind_inc %>% 
      mutate(onset     = onset     + exposed_t,
             symp_end  = symp_end  + exposed_t) 
    
    ## need to ditch dead columns
    rm(ind_inc)
    
    incubation_times_out <- left_join(input, incubation_times_out,
                                      by = c("index_test_delay", "delay_scaling"))
    
    
    #source('kucirka_fitting.R',local=T)  
    
    #calc outcomes 
    message("Calculating outcomes for each traveller")
    incubation_times_out %<>% calc_outcomes(., dat_gam)
    
    #when released
    message("Calculating when travellers released")
    incubation_times_out %<>% when_released()
    
    #stage of infection when released
    #message("Calculating infection status on release")
    # incubation_times_out %<>% stage_when_released()
    
    message("Transmission potential of released travellers")
    incubation_times_out %<>% transmission_potential
    
    return(incubation_times_out)
    
  }

run_rr_analysis <- function(
  released_times,
  main_scenarios,
  baseline_scenario,
  text_size = 2.5,
  log_scale=TRUE,
  y_var = infectivity_post,
  faceting = ~ stringency){
  set.seed(145)
  
  #browser()
  #Parameters
  
  y_var <- enquo(y_var)
  
  baseline <- inner_join(baseline_scenario, input )
  
  stringencies <- distinct(released_times, stringency, scenario)
  
  
  released_times_summaries <- 
    mutate(released_times, 
           time_in_iso = released_t - traced_t) %>% 
    mutate(infectivity=!!y_var) %>% 
    select(scenario,stringency,sim,idx,released_test,infectivity, 
           index_test_delay, -contains("delay"),screening)
  
  
  baseline_summaries <- 
    inner_join(released_times_summaries,
               baseline) %>% 
    rename("baseline_infectivity"   = infectivity,
           #"baseline_released_test" = released_test,
           "baseline_scenario"      = scenario,
           "baseline_stringency"    = stringency) %>%
    select(sim,idx, contains("baseline"), -screening,-pathogen,-contains("delay") ,index_test_delay)
  
  
  n_risk_ratios <- released_times_summaries %>% 
    filter(grepl(x = released_test,
                 pattern ="Released after"),
           !grepl(x = released_test,
                  pattern = "\\+")) %>%
    inner_join(baseline_summaries) %>%   
    mutate(ratio=(infectivity)/(baseline_infectivity)) %>% 
    replace_na(list(ratio=1)) %>% 
    nest(data = -c(scenario,baseline_scenario,index_test_delay)) %>%
    mutate(Q = map(.x = data, ~quantile(.x$ratio, probs = probs)),
           M = map_dbl(.x = data, ~mean(.x$ratio))) %>%
    unnest_wider(Q) %>%
    select(-data) %>%
    inner_join(stringencies) %>%
    inner_join(input)
  
  ylabA <- sprintf("Rate ratio of infectivity in comparison to\n%s stringency, %i day quarantine, %s scenario",
                   baseline_scenario$stringency,
                   with(baseline_scenario,
                        first_test_delay + screening + 
                          ifelse(is.na(second_test_delay), 0, second_test_delay)),"no testing")
  
  rr_fig_data <-
    plot_data(input, n_risk_ratios, main_scenarios = NULL) 
  
  rr_fig <- rr_fig_data %>%
    make_release_figure(
      x         = .,
      input     = input,
      text_size = text_size,
      xlab      = xlab,
      ylab      = ylabA, 
      log_scale = log_scale,
      hline     = "dashed",
      faceting  = faceting)  
  
  
  file <- paste(names(baseline_scenario),
                baseline_scenario, sep = "-", 
                collapse = " ")
  
  list("png", "pdf") %>%
    map(~ggsave(filename = paste0("results/rr_figs_baseline_",
                                  file,".",.x),
                plot=rr_fig,
                width = 260, 
                height = 80*nrow(distinct(ungroup(n_risk_ratios),
                                          !!lhs(faceting))), units="mm",
                dpi = 320,
                device = ifelse(.x=="pdf",cairo_pdf,
                                "png")))
  
  
  return(rr_fig_data)
}

make_days_plots <- 
  function(x, 
           #input,
           main_scenarios = NULL,
           log_scale = FALSE,
           #fixed = TRUE,
           text_size = 2.5,
           #trav_vol_manual = NULL,
           xlab = "Days in quarantine\n(including 1 day delay on testing results)",
           sum = FALSE,
           y_labels = NULL, # MUST BE PASSED IN!!!
           # pass in y_vars as a named list
           faceting = NULL,
           base = stringi::stri_rand_strings(1, 8)){
    
    #browser()
    all_grouping_vars <- all.vars(faceting)
    
    if (sum){
      y_labels <- sub("^Average", "Total", y_labels)
    } 
    
    #browser()
    
    ## need to do summaries over y_labels and x
    ## need a list length(y_labels) long, so perhaps an inner function mapping over?
    
    
    
    x_days_summaries <-
      as.list(names(y_labels)) %>%
      lapply(X = ., FUN = function(y){
        map_df(.x = x,
               ~make_released_time_quantiles(.x,
                                             y_var = y, sum = sum,
                                             vars = all_grouping_vars))})
    ## end summaries
    
    fig_data <- x_days_summaries %>% 
      map(~plot_data(input = input, # should we still pass it in? recover?
                     x_summaries = 
                       .x,
                     main_scenarios))
    
    figs <- map2(
      .x = fig_data,
      .y = y_labels,
      .f = ~make_release_figure(
        x         = .x,
        #input     = input,
        xlab      = "Days in quarantine\n(including 1 day delay on testing results)",
        text_size = text_size,
        ylab      = .y,
        faceting  = faceting,
        percent   = TRUE) )
    
    
    fig <-  wrap_plots(figs, nrow=1,
                       guides = "collect")
    
    if (length(y_labels) > 1L) {
      fig <- fig + 
        plot_annotation(tag_levels = "A")
    }
    
    fig <- fig & theme(legend.position = "bottom")
    
    
    list("png", "pdf") %>%
      map(~ggsave(filename = paste0("results/days_plots_",base,".",.x),
                  plot=fig,
                  width  = 60*nrow(distinct(fig_data[[1]][,get.vars(rhs(faceting))]))*
                    length(fig_data), 
                  height = 80*nrow(distinct(fig_data[[1]][,get.vars(lhs(faceting))])), 
                  dpi = 300,
                  units = "mm",
                  device = ifelse(.x == "pdf",
                                  cairo_pdf,
                                  "png")))
    
    return(
      set_names(fig_data, sub(pattern = "infectivity_", 
                              replacement = "", x = names(y_labels)))
      
    )
    
  }


summarise_results <- function(x, reduction = TRUE){
  if (!is.logical(reduction)){
    stop("reduction must be logical")
  }
  x <- mutate_at(x, 
                 .vars = vars(contains("%")),
                 .funs = function(x){
                   reduction*(1 - x) +
                     (1 - reduction)*x})
  if (reduction){
    percentages <- grep(x = names(x), pattern = "%")
    x[,rev(percentages)] <- x[,percentages]  
  }
  
  x
  
}

show_results <- function(x, reduction = TRUE){
  select(x, delays,
         one_of(all.vars(faceting)),
         screening, time_in_iso,
         contains("%")) %>%
    group_by(stringency, index_test_delay) %>%
    group_split %>%
    map(summarise_results, reduction = reduction) 
}

transmission_potential <- function(x){
  #browser()
  
  x %<>% 
    mutate(
      onset_sec   = onset      - exposed_t,
      release_sec = released_t - exposed_t, 
      b = (onset + 10) - onset + infect_shift,
      q_release = released_t - onset + infect_shift,
      q_traced  = traced_t   - onset + infect_shift,
      a = 0) %>%
    mutate(
      post_untruncated        = pgamma(q     = q_release, 
                                       shape = infect_shape,
                                       rate  = infect_rate, 
                                       lower.tail = F),
      infectivity_denominator = pgamma(q     = b,
                                       shape = infect_shape,
                                       rate  = infect_rate),
      time_from_b_to_Inf      = 1 - infectivity_denominator,
      infectivity_post        = ifelse(b < q_release,
                                       0, 
                                       (post_untruncated - time_from_b_to_Inf)/
                                         infectivity_denominator)
    ) %>%
    mutate(
      pre_untruncated         = pgamma(q     = q_traced,
                                       shape = infect_shape,
                                       rate  = infect_rate),
      infectivity_pre         = pmin(1, 
                                     pre_untruncated/infectivity_denominator)
    ) 
  
  #browser()
  
  x %<>%
    
    mutate(
      quar_untruncated    =
        pmap_dbl(.l = list(q_traced,
                           q_release, 
                           waning),
                 .f = ~integrate(
                   f = function(x){
                     dgamma(x, shape = infect_shape, rate  = infect_rate) * 
                       (1 - get(..3)(x - ..1))
                   },
                   lower = ..1,
                   upper = ..2)$value),
      infectivity_quar    =
        pmax(0,pmin(1,quar_untruncated/infectivity_denominator)),
      infectivity_total   =
        (infectivity_quar + infectivity_post + infectivity_pre),
      infectivity_averted = 1 - infectivity_total)
  
  # , 
  #               infectivity_pre = ptrunc(spec="gamma",
  #                                        a= traced_t -onset +infect_shift - 10,
  #                                        b= traced_t -onset +infect_shift + 10,
  #                                        q = traced_t - onset + infect_shift,
  #                                         shape = infect_shape,
  #                                         rate  = infect_rate,
  #                                         lower.tail = TRUE))
  
  return(x)
}

#function with if statements removed
ptrunc_v <- function (q, spec, a = -Inf, b = Inf, ...) 
{
  #browser()
  tt <- q
  aa <- rep(a, length(q))
  bb <- rep(b, length(q))
  G <- get(paste("p", spec, sep = ""), mode = "function")
  tt <- G(apply(cbind(apply(cbind(q, bb), 1, min), aa), 1, 
                max), ...)
  tt <- tt - G(aa, ...)
  G.a <- G(aa, ...)
  G.b <- G(bb, ...)
  result <- tt/(G(bb, ...) - G(aa, ...))
  return(result)
}


rtgamma <- function(n = 1, a = 0, b = Inf, shape, rate = 1, scale = 1/rate){
  
  p_b <- pgamma(q = b, shape = shape, rate = rate)
  p_a <- pgamma(q = a, shape = shape, rate = rate)
  
  u   <- runif(n = n, min = p_a, max = p_b)
  q   <- qgamma(p = u, shape = shape, rate = rate)
  
  return(q)
}


check_unique_values <- function(df, vars){
  # given a data frame and a vector of variables to be used to facet or group, 
  # which ones have length < 1?
  
  l <- lapply(X = vars, 
              FUN =function(x){
                length(unique(df[, x]))
              })
  
  vars[l > 1]
  
}


waning_piecewise_linear <- function(x, ymax, ymin, k, xmax){
  
  if (ymin == ymax){
    Beta = c(0, ymin)
  } else {
    
    Beta <- solve(a = matrix(data = c(xmax, 1,
                                      k,    1),    ncol = 2, byrow = T),
                  b = matrix(data = c(ymin, ymax), ncol = 1))
  }
  
  (x >= 0)*pmin(ymax, pmax(0, Beta[2] + Beta[1]*x))
  
}

waning_points <- function(x, X, Y, log = FALSE){
  
  if (length(X) != length(Y)){
    stop("X and Y must be same length")
  }
  
  if (length(Y) == 1){
    return(rep(Y, length(x)))
  }
  
  if (log){
    Y <- log(Y)
  }
  
  Beta <- solve(a = cbind(X, 1), b = matrix(Y,ncol=1))
  
  Mu <- Beta[2] + Beta[1]*x
  if (log){
    Mu <- exp(Mu)
  }
  (x >= 0)*pmax(0, Mu)
  
}

