# parameters for simulation

# Ashcroft et al. infectivity profile
# https://smw.ch/article/doi/smw.2020.20336
infect_shape = 97.18750 
infect_rate  =  3.71875
infect_shift = 25.62500

# # McAloon et al. incubation period meta-analysis
#https://bmjopen.bmj.com/content/10/8/e039652
inc_parms <- list(mu_inc = 1.63,
                  sigma_inc = 0.5)

gen_shape <- 9.89
gen_rate <- 2.06

pathogen <- list(
  symptomatic = 
    # review paper Byrne et al. (2020) https://doi.org/10.1101/2020.04.25.20079889
    # define T1 as infection to beginning of presymptomatic infectious period
    append(
      # https://www.acpjournals.org/doi/10.7326/M20-0504
      inc_parms,
      
      # Li et al https://www.nejm.org/doi/full/10.1056/nejmoa2001316
      # variance calculated by inverting confidence interval
      list(mu_inf    = 9.1,
           sigma_inf = 14.7)),
  
  asymptomatic = 
    append(
      # https://www.acpjournals.org/doi/10.7326/M20-0504
      inc_parms,
      # https://doi.org/10.1101/2020.04.25.20079889
      list(
        mu_inf    =  6,
        sigma_inf = 12))) %>%
  map(~data.frame(.x), .id = "type")

# https://www.medrxiv.org/content/10.1101/2020.04.25.20079103v3
asymp_fraction <- rriskDistributions::get.beta.par(
  q = c(0.24,  0.38),
  p = c(0.025, 0.975), 
  show.output = F, plot = F) %>%
  as.list

waning_none <- function(x){
  waning_points(x, X = 0, Y = 1)
}

adhere_10 <- function(x){
  waning_points(x, X = 0, Y = 0.1)
}
adhere_20 <- function(x){
  waning_points(x, X = 0, Y = 0.2)
}
adhere_30 <- function(x){
  waning_points(x, X = 0, Y = 0.3)
}
adhere_40 <- function(x){
  waning_points(x, X = 0, Y = 0.4)
}
adhere_50 <- function(x){
  waning_points(x, X = 0, Y = 0.5)
}
adhere_60 <- function(x){
  waning_points(x, X = 0, Y = 0.6)
}
adhere_70 <- function(x){
  waning_points(x, X = 0, Y = 0.7)
}
adhere_80 <- function(x){
  waning_points(x, X = 0, Y = 0.8)
}
adhere_90 <- function(x){
  waning_points(x, X = 0, Y = 0.9)
}
adhere_100 <- function(x){
  waning_points(x, X = 0, Y = 1)
}

waning_constant <- function(x){
  waning_points(x, X = 0, Y = 0.75)
}

# waning_drop <- function(x){
#   waning_piecewise_linear(x, 0.75, 0.25, 7, 14)
# }

# waning_linear <- function(x){
#   waning_piecewise_linear(x, ymax = 0.75, .16, 0, 8.3)
# }

# waning_canada_community <- function(x){
#   waning_points(x, X = c(0, 30), Y = c(1, 0.541), log = T)
# }

waning_canada_total <- function(x){
  waning_points(x, X = c(0, 30), Y = c(1, 0.158), log = T)
}

smith_uk <- function(x){
  waning_points(x, X = 0, Y = 0.109)
}

default_testing <- c(1, 3, 5, 7, 10, 14)