install.packages()
# dplyr,ggplot2,toolPhD
list.files("src",pattern = ".R$",full.names = T)

dir.create("./result", showWarnings = FALSE)

source("src/Fig1.R")
source("src/Fig2.R")
source("src/Fig3_prep.R")
source("src/Fig3.R")
source("src/Fig4_5.R")
source("src/Fig6_7.R")