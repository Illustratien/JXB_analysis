rm(list=ls())
pacman::p_load(dplyr,ggplot2,ggforce,ggbiplot,purrr,gtable,grid,
               toolPhD,factoextra,toolStability)
options(dplyr.summarise.inform = FALSE)

lookup <- data.frame(
  ori =c("plot_grain_number","plot_tkw","S61_87","S77_87","gSPAD_61","kernel_Nitrogen87",
         "GCD61" , "straw_CHO87",
         'tiller_61','tiller_31',"FE","mobile61_DM","gSPAD_duration","mobile61_N","S23_31","S31_41"),
  new=c('GN','TGW','TT[61-87]',
        'TT[77-87]','SPAD[61]','"["*N*"]"["grain, 87"]','GCD','"["*WSC*"]"["87,straw"]',
        'tiller[61]','tiller[31]',"FE","DM[mobile]","GFD","N[mobile]","tiller","SE"
  ))
para_vec <- c("plot_grain_number",
          "plot_tkw", 'tiller_61',
          'S71_77',
          'gSPAD_61',"GCD61","straw_CHO87"
 )

tiff(filename='result/Fig.6.tiff',
     units="cm",
     width=21,#21.6
     height=15,#16.1
     compression = "lzw",
     pointsize=3,
     res=500,# dpi)
     family="Arial")
source("src/fun/pca_fun.R")
pca.check(para_vec)
dev.off()

tiff(filename='result/Fig.7.tiff',
     units="cm",
     width=21,#21.6
     height=15,#16.1
     compression = "lzw",
     pointsize=3,
     res=500,# dpi)
     family="Arial")
source("src/fun/pca_fun.R")
pca.check(para_vec,typ = "mana")
dev.off()
