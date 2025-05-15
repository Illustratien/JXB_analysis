rm(list=ls())
packages_to_check <- c('dplyr','ggplot2','magick','ggcorrplot','toolPhD','tibble',
                       'scales','ggh4x','rlang','gridExtra','conflicted','agricolae',
                       'tidyr','ggpmisc','parallel','doParallel',
                       'pacman','plyr','purrr','knitr','toolStability','cowplot','ggrepel','corrplot')
check_and_install <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    install.packages(package_name, repos = "https://cran.r-project.org")
    #  requireNamespace(package_name, quietly = TRUE) # No need to load here, the user will load.
  } else {
    message(paste(package_name, "is already installed.")) # Inform the user that the package was already installed
  }
}

# Use lapply to apply the function to each package in the vector
lapply(packages_to_check, check_and_install)


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
source("src/FigS5_6.R")
