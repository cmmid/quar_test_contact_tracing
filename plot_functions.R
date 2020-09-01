## functions used for ploting


covid_pal <- c("#e66101", "#5e3c99", "#0571b0")

#extrafont::loadfonts()
pdf.options(useDingbats=FALSE)


test_labeller <- function(x){
  mutate(x,
         stringency = factor(stringency,
                             levels = c("none",
                                        "one",
                                        "two"),
                             labels = c("None",
                                        "One",
                                        "Two"),
                             ordered = T))
}

type_labels <- c("asymptomatic" =
                   "Asymptomatic",
                 "symptomatic" =
                   "Pre-symptomatic")