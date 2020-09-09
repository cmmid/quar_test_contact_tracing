if (!require(pacman)){
  install.packages("pacman")
}

pacman::p_load(char = c("tidyverse",
                        "data.table",
                        "grid",
                        "ggh4x",
                        "stringi",
                        "scales",
                        "mgcv",
                        "furrr",
                        "tictoc",
                        "countrycode",
                        "eurostat",
                        "conflicted",
                        "magrittr",
                        "epitools",
                        "rlang",
                        "cowplot",
                        "ggnewscale",
                        "formula.tools",
                        "extrafont",
                        "measurements",
                        "distcrete",
                        "fitdistrplus",
                        "magick",
                        "pdftools",
                        "truncdist",
                        "rriskDistributions",
                        "ggridges",
                        "readODS"))

conflicted::conflict_prefer("set_names", "purrr")
conflicted::conflict_prefer("melt", "reshape2")
map(.x = c("mutate", "select", "filter"), 
    .f = function(x){conflicted::conflict_prefer(name = x, "dplyr")})

