###### use this file to prepend the appropriate yaml and knit a set of slides
###### sam wrote this

library(magrittr)
library(purrr)
library(rmarkdown)
library(here)

fileIn <- "slides.Rmd"

commonheader <- "Infrastructure/commonlectureheader.yml"

process_Rmd <- function(fileIn,
                        commonheader, 
                        output_format = "all",
                        delete_Rmd = FALSE){
    
    commonheader <- readLines(con = commonheader)
    
    textRmd <- readLines(con = fileIn) 
    
    fileOut <- sub(pattern = ".Rmd",
                   replacement = "_processed.Rmd",
                   x = fileIn,
                   fixed = T)
    
    textHead_limits <- grep(pattern = "---",
                            x = textRmd)
    
    textHead_indices <- seq(textHead_limits[1], 
                            textHead_limits[2])
    
    textHead <- textRmd[textHead_indices]
    textBody <- textRmd[-textHead_indices]
    
    textHeadToKeep <-
        grep(value = T,
             pattern = "(^((T|t)itle)|(A|a)uthor|(D|d)ate):", 
             x = textHead)
    
    writeLines(text = c("---",
                        textHeadToKeep,
                        commonheader,
                        "---",
                        textBody),
               con  = fileOut)
    
    rmarkdown::render(input = fileOut, 
                      output_format = output_format)
    
    if (delete_Rmd){
        file.remove(fileOut)
    }
    
    
}

fileIn %>%
    as.list %>%
    map(~process_Rmd(.x,
                     commonheader,
                     output_format = "powerpoint_presentation",
                     delete_Rmd = TRUE)
    )

