covid_pal <- c("#e66101", "#5e3c99", "#0571b0")

extrafont::loadfonts()
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
                           post_flight_screening     ~ "One",
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
                                    "Released after mandatory isolation"),
                  pre_board_screening = c(NA,1,4,7)),
       `moderate` = 
         crossing(released_test = c("Released after first test",
                                    "Released after mandatory isolation"),
                  pre_board_screening = c(NA,1,4,7)),
       `high` = 
         crossing(released_test = "Released after second test",
                  pre_board_screening = c(NA,1,4,7)),
       `maximum` = 
         crossing(released_test = c("Released after first test",
                                    "Released after mandatory isolation"),
                  pre_board_screening = c(NA,1,4,7))
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

input <- 
  tibble(pathogen = "SARS-CoV-2") %>%
  mutate(syndromic_sensitivity = syndromic_sensitivity)  %>%
  bind_cols(., list(
    `low` = 
      crossing(pre_board_screening = c(NA,1,4,7),
               post_flight_screening = c(TRUE,FALSE),
               first_test_delay = 0,
               second_test_delay = NA), 
    `moderate` = 
      crossing(pre_board_screening = c(NA,1,4,7),
               post_flight_screening = c(TRUE,FALSE),
               first_test_delay = c(3,5,7),
               second_test_delay = NA),
    `high` = 
      crossing(pre_board_screening = c(NA,1,4,7),
               post_flight_screening = TRUE,
               first_test_delay =  c(0:3),
               second_test_delay = c(2,4,6)),
    `maximum` = 
      crossing(pre_board_screening = c(NA,1,4,7),
               post_flight_screening = c(TRUE,FALSE),
               first_test_delay = 14,
               second_test_delay = NA)) %>%
      bind_rows(.id = "stringency")) %>% 
  mutate(scenario=row_number()) 

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
      gamma.parms.from.quantiles(q = c(5.1, 11.5),
                                 p = c(0.5, 0.975)) %>%
        {list(shape = .$shape, scale = .$scale)}  %>% 
        {gamma2mv(.$shape, scale = .$scale)} %>% 
        set_names(., c("mu_inc", "sigma_inc")),
      
      # Li et al https://www.nejm.org/doi/full/10.1056/nejmoa2001316
      {c(9.1, 14.7)} %>% 
        set_names(., c("mu_inf", "sigma_inf"))),
  
  
  
  asymptomatic = 
    append(
      # https://www.acpjournals.org/doi/10.7326/M20-0504
      gamma.parms.from.quantiles(q = c(5.1, 11.5),
                                 p = c(0.5,0.975)) %>%
        {list(shape = .$shape, scale = .$scale)}  %>% 
        {gamma2mv(.$shape, scale = .$scale)} %>% 
        set_names(., c("mu_inc", "sigma_inc")),
      # https://doi.org/10.1101/2020.04.25.20079889
      list(
        mu_inf    =  6,
        sigma_inf = 12))) %>%
  map(~data.frame(.x), .id = "type")

asymp_fraction <- beta.parms.from.quantiles(q = c(0.03,  0.55),
                                            p = c(0.025, 0.975)) %>%
  {list(shape1 = .$a,
        shape2 = .$b)}





arrivals_fun <- function(x, trav_vol, sims=1){
  
  weekly_arrivals <- rbinom(n = sims, size = trav_vol, prob = 7/30)  
  
  infected_arrivals <- rmultinom(n = sims, 
                                 size = weekly_arrivals,
                                 prob = x$prev) %>%
    unlist %>%
    {data.frame(arrivals = ., type = x$type)} %>%
    filter(type != "uninfected")
  
  #rbinom(n = sims, prob=prev, size = weekly_arrivals)
  
  infected_arrivals
}

# red zones are countries with a higher prevalence of active cases

prev_est <- read_csv("data/currentPrevalenceEstimates_20_07_2020.csv")

zones <- 
  dplyr::select(prev_est, country, propCurrentlyInfMid) %>%
  {set_names(pull(., propCurrentlyInfMid), pull(., country))} %>%
  outer(X = ., Y = ., FUN = "-") %>%
  as.data.frame %>%
  rownames_to_column(var = "destination") %>%
  gather(origin, value, -destination) %>%
  mutate(zone = ifelse(value < 0,
                       "red",
                       "green")) %>%
  filter(origin != destination) %>%
  left_join(dplyr::select(prev_est, 
                          origin = country, 
                          origin_prev = propCurrentlyInfMid),
            by = "origin")

# this might be better done with a test of equal proportions 
# requires reverse-engineering the standard errors
# countries with prevalences not significantly different would be in each others' green zones

## now we need to get the prevalence in the EU

prev_est_region <-
  prev_est %>%
  mutate(region = countrycode::countrycode(country,
                                           "country.name.en",
                                           "eu28", warn = FALSE),
         region = ifelse(country == "United States of America",
                         "USA",
                         region)) %>%
  filter(!is.na(region),
         country != "United Kingdom") %>%
  rename(country.name.en = country,
         country = region) %>%
  nest(q = c(propCurrentlyInfMid, propCurrentlyInfLow)) %>% 
  mutate(gamma_parms = map(.x = q,
                           .f = ~gamma.parms.from.quantiles(q = unlist(.x), 
                                                            #start.with.normal.approx = T,
                                                            p = c(0.5, 0.025))
  ),
  gamma_parms_ = map(gamma_parms, ~list(shape = .x$shape,
                                        rate  = .x$rate,
                                        scale = .x$scale))) %>%
  unnest_wider(gamma_parms_)


prev_est_eu <- 
  prev_est %>%
  inner_join(select(eurostat::eu_countries, country = name),
             by = "country") %>%
  filter(country != "United Kingdom") %>%
  #select(country, propCurrentlyInfMid, population) %>%
  ungroup %>%
  {bind_cols(
    summarise(., propCurrentlyInfMid = weighted.mean(x = propCurrentlyInfMid,
                                                     w = population)),
    summarise_at(., .vars = vars(population, totalCases, totalNewCases),
                 .funs = sum))
  } %>%
  mutate(country = "EU")

prev_est %<>% bind_rows(prev_est_eu)

## flight volumes
## source: CAA tables 10.1, 12.1 for July 2019, with May year on year change
US_flight_vol <- read.csv("data/Table_12_1_Intl_Air_Pax_Traffic_Route_Analysis.csv") %>% 
  filter(foreign_country=="USA") %>% 
  summarise(total=sum(total_pax_this_period))

flight_vols <-
  bind_rows(
    data.frame(year = 2019,
               EU  = 18186680,
               USA =  US_flight_vol$total),
    data.frame(year = 2020,
               EU  = 18186680*0.01,
               USA = US_flight_vol$total*0.01)) 


flight_times <- data.frame(country    = c("EU", "USA"),
                           dur_flight = c(2/24, 8/24))

make_proportions <- function(prev_est_region,
                             origin_country = "United Kingdom",
                             asymp_parms,
                             n = 1){
  
  prev <- prev_est_region %>% 
    filter(country == origin_country) %>%
    mutate(prev = rgamma(n= n(),
                         shape = shape,
                         rate = rate)) %>%
    ungroup %>%
    
    summarise(., prev = weighted.mean(x = prev,
                                      w = population)) %>% unlist
  
  
  # sample prevalences from the country
  
  prop.asy <- rbeta(n = n, 
                    shape1 = asymp_parms$shape1,
                    shape2 = asymp_parms$shape2)
  
  #prev <- unlist(filter(prev_est, country == origin_country)$propCurrentlyInfMid)
  
  data.frame(symptomatic = (1-prop.asy)*prev,  # symptomatic
             asymptomatic = prop.asy*prev) %>%
    mutate(sim = 1:n())# will never have symptoms
  
}

make_prevalence <- function(prev_est_region,
                            origin_country = "United Kingdom",
                            n = 1){
  
  prev_est_region %>% 
    filter(country == origin_country) %>%
    split(.$country.name.en) %>%
    map_df(~data.frame(prev = rgamma(n = n, shape = .x$shape, rate = .x$rate),
                       id = 1:n),
           .id = "country.name.en") %>%
    inner_join(., select(prev_est_region, country.name.en, population),
               by = "country.name.en") %>%
    group_by(id) %>%
    dplyr::summarise(prev = weighted.mean(x = prev,
                                          w = population),
                     .groups = "drop") %$% prev
  
  
}


gen_screening_draws <- function(x){
  #browser()
  
  n <- nrow(x)
  
  # generate screening random draws for comparison
  x <- mutate(x, 
              screen_0 = runif(n, 0, 1),  # pre-departure
              screen_1 = runif(n, 0, 1),  # on arrival
              screen_2 = runif(n, 0, 1),  # follow-up
              screen_s = runif(n, 0, 1))  # syndromic screening
}

# given infection histories above, what proportion of travellers end up being 
# caught at each step in the screening process?

calc_outcomes <- function(x, dat_gam){
  
  # generate required times for screening events
  x <- mutate(x,
              symp_screen_t       = flight_departure,
              flight_arrival      = flight_departure + dur_flight,
              first_test_t        = flight_arrival + first_test_delay,
              second_test_delay   = second_test_delay,
              second_test_t       = first_test_t + second_test_delay,
              zeroth_test_t       = flight_departure - pre_board_screening)
  
  # what's the probability of PCR detection at each test time?
  x <- mutate(x,
              zeroth_test_p =
                c(predict(object  = dat_gam,
                          type    = "response",
                          newdata = data.frame(day = zeroth_test_t))),
              first_test_p = 
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
              zeroth_test_p = ifelse(zeroth_test_t < 0,
                                     0,
                                     zeroth_test_p))
  
  # make comparisons of random draws to screening sensitivity
  # conduct symptomatic screening as well
  x <-
    mutate(x,
           zeroth_test_label      = detector(pcr = zeroth_test_p, u = screen_0),
           first_test_label       = detector(pcr = first_test_p,  u = screen_1),
           second_test_label      = detector(pcr = second_test_p, u = screen_2),
           symp_screen_label      = 
             symp_screen_t > onset &                     # after onset of symptoms
             symp_screen_t < symp_end &                  # still symptomatic
             runif(n = n(), 0, 1) < syndromic_sensitivity) # detected
  
}

when_released <- function(x){
  
  mutate(x, released_test = case_when(
    
    zeroth_test_label                           ~ 
      "Positive at zeroth test, prevented from boarding",
    
    symp_screen_label  & type != "asymptomatic" ~ 
      "Symptomatic at departure, prevented from boarding",
    
    !post_flight_screening    ~ 
      "Released after mandatory isolation",
    
    
    post_flight_screening & !first_test_label & is.na(second_test_label) ~
      "Released after first test",
    
    post_flight_screening & !first_test_label & !second_test_label     ~
      "Released after second test",
    
    first_test_label                            ~
      "Released after first test + mandatory quarantine",
    
    !first_test_label  & second_test_label      ~
      "Released after second test + mandatory quarantine",
    
    TRUE                                        ~ 
      "Mandatory quarantine"
  ),
  released_t = case_when(
    
    released_test == "Symptomatic at departure, prevented from boarding"    ~
      NA_real_, #
    
    released_test == "Positive at zeroth test, prevented from boarding"     ~
      NA_real_, # if they never board, they don't contribute to quarantine
    
    released_test == "Released after mandatory isolation"                   ~
      flight_arrival + first_test_delay, # BILLY TO CHECK
    
    released_test == "Released after first test"                           ~ 
      first_test_t + 1,
    
    released_test == "Released after second test"                           ~ 
      second_test_t + 1,
    
    released_test == "Released after first test + mandatory quarantine"     ~ 
      first_test_t  + 14,
    
    released_test == "Released after second test + mandatory quarantine"    ~
      second_test_t + 14,
    
    released_test == "Mandatory quarantine"                                 ~
      flight_arrival + 14)) %>% 
    mutate(released_test =  ifelse(type == "symptomatic" & 
                                     onset > flight_arrival &
                                     onset < released_t,
                                   "Symptomatic during quarantine",
                                   released_test),
           released_t = ifelse(released_test == "Symptomatic during quarantine",
                               flight_arrival + pmax(onset + 7, symp_end, 14), released_t))
}

stage_when_released <- function(x){
  
  mutate(x, stage_released = as.factor(case_when(
    is.na(released_t)         ~ "Prevented from boarding",
    released_t < inf_end      ~ "Infectious",
    released_t >= inf_end     ~ "Post-infectious"
  ))) %>% mutate(
    days_released_inf = 
      case_when(stage_released == "Post-infectious" ~ 
                  0,
                stage_released == "Prevented from boarding" ~ 
                  NA_real_, 
                type           == "asymptomatic" ~ 
                  inf_end - pmax(released_t, inf_start),
                type           == "symptomatic" ~ 
                  onset   - pmax(released_t, inf_start)))
  
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
                           contains("test_delay")),
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

bivariate_color_scale <- tibble(
  "3 - 6" = "#574249", 
  "2 - 6" = "#5B6570",
  "1 - 6" = "#608997", 
  "0 - 6" = "#64ACBE", 
  "3 - 4" = "#985356",
  "2 - 4" = "#A07E84", 
  "1 - 4" = "#A8AAB1",
  "0 - 4" = "#B0D5DF",
  "3 - 2" = "#C85A5A", 
  "2 - 2" = "#D38989",
  "1 - 2" = "#DDB9B9",
  "0 - 2" = "#E8E8E8"
) %>%
  gather("delays_minus", "colour") %>% 
  mutate(delays=str_replace_all(delays_minus," - ", " + "))

bivariate_color_scale_leg <- bivariate_color_scale %>%
  separate(delays_minus, into = c("first_test_delay", "second_test_delay"), sep = " - ")

legend <- ggplot() +
  geom_tile(
    data = bivariate_color_scale_leg,
    mapping = aes(
      x = first_test_delay,
      y = second_test_delay,
      fill = colour)
  ) +
  scale_fill_identity() +
  labs(x = "Days until first test",
       y = "Days after first test\n until second test") +
  theme_minimal()+
  labs(subtitle = "Two test regimen")+
  theme(plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 12)
  ) +
  # quadratic tiles
  coord_fixed()

make_incubation_times <- function(
  n_travellers,
  pathogen,
  syndromic_sensitivity = 0.7){
  
  incubation_times <- crossing(idx  = 1:n_travellers,
                               type = c("symptomatic",
                                        "asymptomatic") %>%
                                 factor(x = .,
                                        levels = .,
                                        ordered = T)) %>%
    split(.$type) %>%
    map2_df(.x = .,
            .y = pathogen,
            ~mutate(.x,
                    exp_to_onset   = time_to_event(n = n(),
                                                   mean = .y$mu_inc, 
                                                   var  = .y$sigma_inc),
                    onset_to_recov = time_to_event(n = n(),
                                                   mean = .y$mu_inf, 
                                                   var  = .y$sigma_inf))) 
  
  source("wolfel.R")
  source("he.R")
  # infectious period from time of onset to no longer infectious
  incubation_times %<>% 
    mutate(u = runif(n = nrow(.), 0.01, 0.99)) %>%
    mutate(inf_from_onset = 
             approx(x    = wolfel_pred$y, 
                    y    = wolfel_pred$day, 
                    xout = u)$y,
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
           inf_start = onset - pre_symp_lead,
           inf_end   = ifelse(type == "asymptomatic",
                              exp_to_onset + onset_to_recov,
                              exp_to_onset + inf_from_onset),
           symp_end  = ifelse(type == "asymptomatic",
                              onset, # but really never matters because asymptomatics are never symptomatic!
                              exp_to_onset + onset_to_recov),
           inf_dur   = inf_end - inf_start,
           symp_dur  = symp_end - onset)
  
  # add flight
  incubation_times %<>% 
    mutate(flight_departure = runif(n = nrow(.),
                                    min = 0,
                                    max = onset + onset_to_recov)) %>%
    mutate(symp_screen_label      = 
             flight_departure > onset &
             flight_departure < symp_end &
             runif(n = nrow(.)) < syndromic_sensitivity)
  
  incubation_times %<>% gen_screening_draws
  
  return(incubation_times)
  
}

make_inf_arrivals <- function(countries,
                              prev_est_region,
                              n_arrival_sims,
                              asymp_fraction,
                              flight_vols = NULL,
                              trav_vol_manual = NULL,
                              trav_vol_p = 1,
                              flight_times,
                              incubation_times,
                              fixed = FALSE){
  
  inf_arrivals <- as.list(countries) %>%
    purrr::set_names(., .) %>%
    purrr::map(~make_prevalence(prev_est_region = prev_est_region,
                                origin_country = .x, 
                                n = n_arrival_sims)) %>%
    purrr::map_dfr(.id = "country", ~data.frame(pi = .x)) %>%
    dplyr::mutate(alpha = rbeta(n = nrow(.),
                                shape1 = asymp_fraction$shape1,
                                shape2 = asymp_fraction$shape2)) %>%
    tidyr::nest(data = -c(country)) 
  
  if (!fixed){
    inf_arrivals <- dplyr::inner_join(inf_arrivals,
                                      dplyr::filter(flight_vols, year == 2020) %>%
                                        tidyr::gather(country, trav_vol),
                                      by = "country") %>%
      dplyr::mutate(trav_vol = trav_vol/2)
  } else {
    inf_arrivals <- dplyr::mutate(inf_arrivals, trav_vol = trav_vol_manual)
  }
  
  inf_arrivals <-  
    mutate(inf_arrivals,
           travellers = purrr::map2(.x = data,
                                    .y = trav_vol,
                                    .f = 
                                      ~make_travellers(
                                        x = .x,
                                        incubation_times = incubation_times,
                                        trav_vol = .y,
                                        trav_vol_p = trav_vol_p,
                                        fixed = fixed))) 
  
  inf_arrivals  <- 
    mutate(inf_arrivals,
           individuals = furrr::future_map(.x = travellers,
                                           .f = travellers_to_individuals,
                                           incubation_times = incubation_times))
  
  inf_arrivals <- tidyr::unnest(inf_arrivals, individuals)
  
  # add flight duration
  inf_arrivals <- dplyr::inner_join(inf_arrivals, flight_times)
  
  return(inf_arrivals)
  
}


make_travellers <- function(x, # contains relevant parameters
                            incubation_times,
                            trav_vol,
                            trav_vol_p = 7/30,
                            xi = NULL,
                            fixed = TRUE){
  
  # incubation_times should be big
  # sims: number of simulations being performed. should be 1.
  # trav_vol: should already be scaled
  
  sims <- nrow(x)
  
  weekly_arrivals <- rbinom(n = sims, size = trav_vol, prob = trav_vol_p)  
  
  if (is.null(xi)){
    xi <- incubation_times %>%
      filter(type == "symptomatic") %>% # all ever-symptomatics
      summarise(xi = mean(symp_screen_label)) %$% xi # what proportion are screened
  }
  
  if (fixed){
    weekly_intended <- weekly_arrivals
  } else {
    weekly_intended <- weekly_arrivals + 
      rnbinom(n    = sims,
              size = weekly_arrivals,
              prob = 1 - x$pi*(1-x$alpha)*syndromic_sensitivity*xi)
  }
  
  uninfected     <- rbinom(sims, weekly_intended, 1 - x$pi)
  infected       <- weekly_intended - uninfected
  asymptomatic   <- rbinom(sims, infected, x$alpha)
  ever_symp      <- infected - asymptomatic
  currently_symp <- rbinom(sims, ever_symp, xi)
  #symp_removed   <- rbinom(sims, currently_symp, syndromic_sensitivity)
  not_symp       <- ever_symp - currently_symp
  
  #total_travellers <- uninfected + asymptomatic + symp_permitted
  
  data.frame(sim = 1:sims,
             uninfected,
             asymptomatic,
             currently_symp,
             #total_travellers,
             not_symp)
}


travellers_to_individuals <- function(x, incubation_times){
  
  keys <- list(asymptomatic   =
                 data.frame(type = 'asymptomatic',
                            symp_screen_label = FALSE),
               not_symp = 
                 data.frame(type = 'symptomatic',
                            symp_screen_label = FALSE),
               currently_symp   = 
                 data.frame(type = 'symptomatic',
                            symp_screen_label = TRUE))
  
  y <- gather(select(x,-uninfected),#,-total_travellers),
              key, value, -sim)
  
  y <- split(y, y$sim)
  
  l <- lapply(y, FUN = function(z){
    
    map2_df(.x = as.list(set_names(z$value, z$key)),
            .y = keys, 
            .f = ~sample_n(inner_join(incubation_times,
                                      .y,
                                      by = c("type",
                                             "symp_screen_label")), 
                           size = .x, replace = F))
  }) 
  
  r <- bind_rows(l, .id = "sim") 
  
  r <- mutate(r, sim = parse_number(sim))
  
  return(r)
  
}


make_arrival_scenarios <- function(input, inf_arrivals, incubation_times){
  source('kucirka_fitting.R', local=T)
  
  arrival_scenarios <- crossing(input, inf_arrivals)
  
  # calculate outcomes of screening
  arrival_scenarios %<>% calc_outcomes(., dat_gam)
  
  return(arrival_scenarios)
  
}

make_release_figure <- function(x,
                                input,
                                xlab = "Days in quarantine",
                                ylab = "",
                                text_size = 2.5,
                                text_angle = 45,
                                h_just = 0,
                                log_scale = FALSE,
                                hline = 0,
                                faceting = NULL){
  
  x %<>% test_labeller
  
  # how to do presymptomatic
  
  facet_vars <- all.vars(faceting)
  
  if ("type" %in% facet_vars){
    x %<>% mutate(type = factor(type,
                                levels = c("asymptomatic",
                                           "symptomatic"),
                                labels = c("Asymptomatic",
                                           "Presymptomatic")))
  }
  
  require(formula.tools)
  
  dy <- select(x, !!!facet_vars,
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
    summarise_all(.funs = max) %>%
    tidyr::pivot_wider(names_from = "stringency", 
                       names_glue = "{stringency}_{.value}",
                       values_from = c(`75%`, `97.5%`))
  
  
  
  needs_expanding <- dy %>% 
    filter(`High_97.5%` > 0.75*`Other_97.5%`) %>%
    {nrow(.) > 0L}
  
  dy %<>% mutate(ypos = 0.05 * pmax(`High_97.5%`, `Other_97.5%`)) %>%
    select(one_of(all.vars(faceting)), ypos)
  
  
  
  x %<>% inner_join(dy)
  
  figure <-  
    ggplot(data=x, aes(x = time_in_iso, 
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
    geom_point(pch = "-", size = 6, 
               position = position_dodge2(width = 0.75),
               aes(y = `50%`,
                   group = delays)
    ) +
    geom_text(data=filter(x, stringency=="High"),
              aes(x     = time_in_iso,
                  y     = `97.5%` + ypos,
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
    facets = faceting,
    nest_line = T,
    scales = "free", space = "free_x")
  # }
  
  # check if the top of the y axis needs adjustment
  
  
  
  # end check top y
  
  mult <- c(0.1, ifelse(needs_expanding, 0.5, 0.1))
  
  if (log_scale) {
    figure <- figure +
      scale_y_log10(labels=format_format(scientific=FALSE),
                    expand = expansion(mult = mult)) +
      theme(panel.grid.minor.y = element_blank()) +
      annotation_logticks(sides = "l")
    
  } else {
    figure <- figure + 
      scale_y_continuous(#limits=c(-dy/5,NA),
        breaks = pretty_breaks(), 
        expand = expansion(mult = mult))
  }
  
  
  return(figure)
  
}




plot_data <- function(input, arrival_released_times_summaries,
                      main_scenarios = NULL){
  
  dat <- input %>%
    mutate(second_test_delay_ = 
             ifelse(is.na(second_test_delay),
                    0,
                    second_test_delay),
           time_in_iso = 
             first_test_delay + 
             second_test_delay_+
             post_flight_screening) %>% 
    # mutate(pre_board_screening = as.factor(pre_board_screening)) %>% 
    # mutate(pre_board_screening = fct_explicit_na(pre_board_screening, "NA")) %>% 
    # mutate(pre_board_screening = factor(pre_board_screening,
    #                                     levels = names(pre_board_labels),
    #                                     labels = pre_board_labels, ordered = T)) %>% 
    
    inner_join(arrival_released_times_summaries) 
  
  if (!is.null(main_scenarios)){
    main_scenarios %<>% select(-released_test) %>% distinct
    dat <- left_join(dat, main_scenarios)
  }
  
  dat %>%  
    tidyr::nest(data = -c(first_test_delay, second_test_delay)) %>%
    tidyr::unite(col = "delays",
                 first_test_delay, second_test_delay,
                 sep = " + ", remove = FALSE) %>%
    tidyr::unnest(data) %>%
    dplyr::mutate(time_in_iso = factor(time_in_iso, 
                                       levels = unique(.$time_in_iso),
                                       ordered = T)) %>%
    # dplyr::mutate(stringency = factor(stringency,
    #                                   levels = c("low",
    #                                              "moderate",
    #                                              "high",
    #                                              "maximum"),
    #                                   labels = c("Low",
    #                                              "Mod.",
    #                                              "High",
    #                                              "Max."),
    #                                   ordered = T)) %>%
    dplyr::filter(M!=0) %>%  # if the mean is zero, this group is empty
    #dplyr::mutate(released_test = factor(released_test,
    #                                     levels = names(released_labels),
    #                                     labels = released_labels, ordered = T)) 
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
  arrival_released_times, 
  input,
  main_scenarios = NULL,
  log_scale = FALSE,
  fixed = FALSE,
  text_size = 2.5,
  trav_vol_manual = NULL,
  xlab = "Days in quarantine\n(including 1 day delay on testing results)",
  sum = FALSE,
  pre_board_screening = FALSE,
  faceting = NULL){
  
  ylabA = "Number of infectious persons\nreleased per"
  if (fixed){
    ylabA <- paste(ylabA, sprintf("%i travellers", trav_vol_manual))
  } else {
    ylabA <- paste(ylabA, "week")
  }
  
  if (sum){
    ylabB = 
      "Total number of person-days infectiousness\nremaining for released travellers"
  } else {
    ylabB = 
      "Number of days infectiousness\nremaining per released traveller"
  }
  
  
  if (!pre_board_screening){
    arrival_released_times <- dplyr::filter(arrival_released_times,
                                            is.na(pre_board_screening))
  }
  
  all_grouping_vars <- all.vars(faceting)
  
  # ... should be country, type, pre_board_screening, etc.
  arrival_released_times_summaries <- arrival_released_times %>%
    filter(stage_released == "Infectious") %>%
    make_arrival_released_quantiles(., 
                                    all_grouping_vars) %>%
    inner_join(input)
  
  # why are the plo medians coming out backwars in Low?
  figS3A_data <- plot_data(input, 
                           arrival_released_times_summaries,
                           main_scenarios)
  
  
  # deal with the join and stage_released # need to have scenario guide the labelling
  
  
  figS3A <- figS3A_data %>% 
    #filter(pre_board_screening == "None") %>% 
    # # this filter should be done outside the function
    make_release_figure(
      x         = .,
      input     = input,
      text_size = text_size,
      xlab      = xlab,
      ylab      = ylabA, 
      faceting  = faceting) 
  
  ## person-days
  
  arrival_released_days_summaries <- 
    make_arrival_released_time_quantiles(arrival_released_times, 
                                         all_grouping_vars,
                                         sum = sum)
  
  
  
  figS3B_data <- plot_data(input = input, 
                           arrival_released_times_summaries = 
                             arrival_released_days_summaries,
                           main_scenarios)
  
  figS3B <- figS3B_data %>% 
    #filter(pre_board_screening == "None") %>%  
    make_release_figure(
      x = .,
      xlab = "Days in quarantine\n(including 1 day delay on testing results)",
      text_size = text_size,
      ylab = ylabB,
      faceting = faceting) 
  
  
  figS3 <- figS3A + figS3B + plot_layout(ncol = 1, guide = "collect") +
    plot_annotation(tag_levels = "A",
                    theme = theme(legend.position = "bottom"))
  
  return(figS3)
  
}



make_arrival_released_quantiles <- function(x, vars){
  
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

make_arrival_released_time_quantiles <- function(x, vars, sum = FALSE){
  
  
  dots1 <- rlang::exprs(sim, scenario)
  dots2 <- lapply(vars, as.name)
  
  dots  <- append(dots1, dots2)
  
  if (sum){
    x <- x %>%
      dplyr::select(!!! dots, days_released_inf) %>%
      group_by_at(.vars = vars(-days_released_inf)) %>%
      summarise(days_released_inf = sum(days_released_inf, na.rm=T))
  }
  
  x_days <- x %>%
    dplyr::select(!!! dots, days_released_inf) %>%
    dplyr::filter(!is.na(days_released_inf), days_released_inf > 0)
  
  x_days %>%
    nest(data = c(days_released_inf, sim)) %>%
    mutate(Q = purrr::map(.x = data, ~quantile( .x$days_released_inf,
                                                probs = probs)),
           M = map_dbl(.x = data, ~mean(.x$days_released_inf))) %>%
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
    
    img <- magick::image_read(file, density = dpi)
    
    purrr::map(.x = device,
               .f = 
                 ~magick::image_convert(image = img, format = .x) %>%
                 magick::image_write(., path = sub(pattern = "pdf",
                                                   replacement = .x,
                                                   x = file)))
  }
}
