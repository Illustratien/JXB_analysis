rm(list=ls())
library(dplyr)
conflicted::conflict_prefer("filter", "dplyr") %>% suppressMessages()
conflicted::conflict_prefer("select", "dplyr") %>% suppressMessages()

dir.create("result/", showWarnings = FALSE)
folder.vec <- c("plot","table")
for (i in folder.vec) {
  if (!dir.exists(paste0("result/", i))) {
    dir.create(paste0("result/", i), showWarnings = FALSE)
  }
}

source("src/Fig1.R")
source("src/Fig2.R")
source("src/Fig3_prep.R")
source("src/Fig3.R")
source("src/Fig4_5.R")
source("src/Fig6_7.R")
source("src/FigS3.R")
source("src/FigS4.R")
