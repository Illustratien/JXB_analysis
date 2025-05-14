rm(list=ls())
pacman::p_load(dplyr,ggplot2,ggforce,ggbiplot,purrr,gtable,grid,
               toolPhD,factoextra,toolStability)
options(dplyr.summarise.inform = FALSE)

source("src/fun/pca_fun.R")

lookup <- data.frame(
  ori =c("plot_grain_number","plot_tkw","S61_87","S77_87",'S71_77',"gSPAD_61","kernel_Nitrogen87",
         "GCD61" , "straw_CHO87",
         'tiller_61','tiller_31',"FE","mobile61_DM","gSPAD_duration","mobile61_N","S23_31","S31_41"),
  new=c('GN','TGW','TT[61-87]',
        'TT[77-87]','TT[71-77]','SPAD[61]','"["*N*"]"["grain, 87"]','GCD','"["*WSC*"]"["87,straw"]',
        'tiller[61]','tiller[31]',"FE","DM[mobile]","GFD","N[mobile]","tiller","SE"
  ))

BASE <- c("plot_grain_number","plot_tkw")


tiff(filename='result/plot/Fig6.tiff',
     units="cm",
     width=21,#21.6
     height=15,#16.1
     compression = "lzw",
     pointsize=3,
     res=500,# dpi)
     family="Arial")
source("src/fun/pca_fun.R")
pca.check(c(BASE,
            'tiller_61',
            'S71_77',
            'gSPAD_61',"GCD61","straw_CHO87")) %>% print()

dev.off()

tiff(filename='result/plot/Fig7.tiff',
     units="cm",
     width=21,#21.6
     height=15,#16.1
     compression = "lzw",
     pointsize=3,
     res=500,# dpi)
     family="Arial")
source("src/fun/pca_fun.R")
pca.check(c(BASE,
            'tiller_61',
            'S71_77',
            'gSPAD_61',"GCD61","straw_CHO87"),typ = "mana") %>% print()

dev.off()
