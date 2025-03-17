# barplot -----------------------------------------------------------------
rm(list=ls())
pacman::p_load(dplyr,ggplot2,toolPhD)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
centroid_dot_pattern <- function(params, boundary_df, aspect_ratio, legend) {
  # https://coolbutuseless.github.io/package/ggpattern/articles/developing-patterns-2.html
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Convert the simple `boundary_df` polygon information into a 
  # simple features polygon object i.e. {sf}
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  boundary_sf <- ggpattern::convert_polygon_df_to_polygon_sf(boundary_df)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Now that we have the boundary as an {sf} object, we can use a simple
  # features' function to find the centroid
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  centroid    <- sf::st_centroid(boundary_sf)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Create a single character at the cenroid
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  grid::pointsGrob(
    x    = centroid[1],
    y    = centroid[2],
    pch  = params$pattern_shape,
    size = unit(params$pattern_size, 'char'),
    gp   = grid::gpar(
      col = ggplot2::alpha(params$pattern_fill, params$pattern_alpha)
    )
  )
}
smtble<- function(df,gf){
  list(ecovalence(df,"Trait",gf,c(setdiff(c('g','var'),gf),"DFG_year"),unit.correct = T),
       genotypic_superiority_measure(df,"Trait",gf,c(setdiff(c('g','var'),gf),"DFG_year"),
                                     unit.correct = T)
  ) %>% Reduce('merge',.)
}
options(ggpattern_geometry_funcs = list(centroid = centroid_dot_pattern))




plot_yield <- read.csv('data/merge_profile.csv') %>% 
  filter(namCombine%in%c("plot_yield"),!DFG_year=="DFG2019")

stable_merg <- bind_rows(
  # management 
  smtble(plot_yield,'g') %>% 
    tidyr::pivot_longer(ecovalence:genotypic.superiority.measure,
                        values_to = 'SI',names_to = 'SI_name') %>% 
    mutate(trait_name="plot_yield",
           Typ="management"),
  # genotype 
  smtble(plot_yield,'var') %>% 
    tidyr::pivot_longer(ecovalence:genotypic.superiority.measure,
                        values_to = 'SI',names_to = 'SI_name') %>% 
    mutate(
      trait_name="plot_yield",
      Typ="genotype")
) %>% rename(mean=Mean.Trait)

Geno <- c("Capone","Pionier","Patras","Apertus","Torrild","Alves","Potenzial","Esket")
col.clr <- c("#003f5c","#58508d","#8a508f","#bc5090", "#de5a79", "#ff6361","#ff8531", "#ffa600")

names(col.clr) <- Geno

col.nitro <- c("#B3AD00","#1F9600")
names(col.nitro) <- c("176","220")
shape_geno <- c(0,1,2,5,6,8,12,13)
names(shape_geno) <- Geno

unit.pi <-quote('t/ha')

sivec <- stable_merg$SI_name %>% unique()

# -------------------------------------------------------------------------
for (si in sivec){
  stable_tabley <- dplyr::filter(stable_merg,trait_name=="plot_yield",SI_name==si)
  
  si.txt <- ifelse(si=="ecovalence","W","P")
  pg.txt <-bquote(bolditalic(.(si.txt))[bold('g,GY')]~bold('('*.(unit.pi)*')'))
  pm.txt <-bquote(bolditalic(.(si.txt))[bold('m,GY')]~bold('('*.(unit.pi)*')'))
  
  manage_st <- stable_tabley%>%
    dplyr::filter(Typ=="management")%>%
    dplyr::mutate(g=Genotype,
                  Genotype=gsub("[_]",'\n',Genotype))%>%
    tidyr::separate(col = g,
                    into=c('application','nitrogen',
                           'sowing_date'),
                    sep = '[_]')%>%
    dplyr::mutate(treatment=interaction(application,sowing_date,sep=" ")%>%
                    factor(.,levels=c("Split Early","Split Late",
                                      "Combined Early","Combined Late")),
                  SI_l=round(SI,2),
                  SI_y=case_when(SI<.1~SI*1.5,
                                 SI<.5~SI*1.15,
                                 SI>.5~SI*1.08,
                                 SI>1~SI*1.01),
                  Shp_y=case_when(SI_y<.1~SI_y*1.5,
                                  SI_y<.5~SI_y*1.15,
                                  SI_y>.5~SI_y*1.08,
                                  SI_y>1~SI_y*1.01))
  
  geno_st <- stable_tabley%>%dplyr::filter(Typ=="genotype")%>%
    mutate(SI_l=round(SI,2),
           SI_y=case_when(SI<.1~SI*1.3,SI<.5~SI*1.15,SI>.5~SI*1.08,SI>1~SI*1.01))
  
  # combined plot of geno and management -------------------------------------------------------------------------
  p_geno2 <- ggplot(geno_st,
                    aes(x=reorder(Genotype,SI),y=SI,label=SI_l)) +
    geom_bar(stat="identity",alpha=.5)+
    ggplot2::scale_y_continuous(
      limits = c(-.18,max(stable_tabley$SI)*1.2))+
    geom_point(aes(y=-.12,shape=Genotype,color=Genotype),
               size=3,stroke=1.5,show.legend = F)+
    scale_shape_manual(values = shape_geno)+
    scale_color_manual(values = col.clr)+
    xlab("Genotype")+
    ylab(pg.txt)+
    theme_phd_facet(ax.txt.siz = 6,ax.tit.siz = 8,t=20)+
    geom_text(aes(y=SI_y),size=3)+
    theme(axis.title.x = element_text(vjust=-12),
          axis.title.x.bottom = element_blank())
  
  p_manage2 <- ggplot(manage_st,
                      aes(x=reorder(Genotype,SI),y=SI,label=SI_l,
                          color=nitrogen)) +
    geom_bar(stat="identity",aes(fill=nitrogen),alpha=.5)+
    theme_phd_facet(ax.txt.siz = 6,ax.tit.siz = 8,t=20)+
    geom_point(aes(y=-.12,shape=treatment),size=3,stroke=1.5)+
    ylab(pm.txt)+
    scale_shape_manual(values = 
                         c("Split Early" = 1, "Combined Early" = 16,
                           "Split Late"=0,"Combined Late"=15))+
    geom_text(aes(y=SI_y),show.legend = F,color="black",size=3,)+
    ggplot2::scale_y_continuous(
      limits = c(-.18,max(stable_tabley$SI)*1.2))+
    ggplot2::theme(
      axis.title.x=element_blank(),
      legend.position = "none")
  
  legend_b <- cowplot::get_plot_component(
    ggplot(manage_st,
           aes(x=Genotype,y=SI,
               color=nitrogen,
               shape=treatment)) +
      geom_point(size=3,alpha=.8,stroke = 1)+
      scale_shape_manual(name = "Application x Sowing Day",
                         values=
                           c("Split Early" = 1, "Combined Early" = 16,
                             "Split Late"=0,"Combined Late"=15))+
      theme_phd_facet(lgd.tit.siz = 6,lgd.txt.siz = 6)+
      theme(plot.margin =margin(1,1,1,5,"cm"),
            legend.spacing.x = unit(.5, 'cm'),
            legend.position = "bottom",
            legend.title =element_text(hjust=-30),
            legend.key.size = unit(.8, 'cm'))+
      guides(color=guide_legend(nrow=1, byrow=TRUE,
                                title.position="top", 
                                title.hjust = 0.5),
             shape=guide_legend(nrow=1, byrow=TRUE,
                                title.position="top", 
                                title.hjust = 0.5)),
    'guide-box-top', return_all = TRUE)
  
  prow2 <-    cowplot::plot_grid(
    cowplot::plot_grid(
      p_geno2,p_manage2,align = 'hv',
      labels =c("A Genotype SI","B Management SI"),
      label_size=10,rel_widths = c(.9,1.1),
      hjust = -1,nrow = 1),
    legend_b,rel_heights = c(.9,.1),ncol=1)
  
}

# separated bar plot -------------------------------------------------------------------------
legend_a <- cowplot::get_plot_component(
  dplyr::filter(stable_merg,trait_name=="plot_yield",Typ=="genotype") %>% 
    select(-mean) %>% 
    tidyr::pivot_wider(values_from = SI,names_from = SI_name)%>%
    ggplot(.,
           aes(x=ecovalence,
               y=genotypic.superiority.measure,
               color=Genotype,
               shape=Genotype)) +
    geom_point(size=3,stroke=1.5)+
    scale_shape_manual(values = shape_geno)+
    guides(color=guide_legend(nrow=1, byrow=TRUE,
                              title.position="top", 
                              title.hjust = 0.5),
           shape=guide_legend(nrow=1, byrow=TRUE,
                              title.position="top", 
                              title.hjust = 0.5))+
    theme_phd_facet(lgd.tit.siz = 6,lgd.txt.siz = 6)+
    theme(plot.margin =margin(1,1,1,5,"cm"),
          legend.spacing.x = unit(.5, 'cm'),
          legend.position = "bottom",
          legend.title =element_text(hjust=-30),
          legend.key.size = unit(.8, 'cm'))+
    scale_color_manual(values = col.clr) ,'guide-box-bottom',return_all=T)

legend_b <- cowplot::get_plot_component(
  ggplot(manage_st,
         aes(x=Genotype,y=SI,
             color=nitrogen,
             shape=treatment)) +
    geom_point(size=3,alpha=.8,stroke = 1)+
    scale_shape_manual(name = "Application x Sowing Day",
                       values=
                         c("Split Early" = 1, "Combined Early" = 16,
                           "Split Late"=0,"Combined Late"=15))+
    theme_phd_facet(lgd.tit.siz = 6,lgd.txt.siz = 6)+
    theme(plot.margin =margin(1,1,1,5,"cm"),
          legend.spacing.x = unit(.5, 'cm'),
          legend.position = "bottom",
          legend.title =element_text(hjust=-30),
          legend.key.size = unit(.8, 'cm'))+
    scale_color_manual(values=col.nitro)+
    guides(color=guide_legend(nrow=1, byrow=TRUE,
                              title.position="top", 
                              title.hjust = 0.5),
           shape=guide_legend(nrow=1, byrow=TRUE,
                              title.position="top", 
                              title.hjust = 0.5)),
  'guide-box-bottom',return_all=T
)

plist <- purrr::map(c("genotype","management"),
                    function(typ){
                      subdata <-   dplyr::filter(stable_merg,trait_name=="plot_yield",Typ==typ) %>% 
                        select(-mean) %>% 
                        tidyr::pivot_wider(values_from = SI,names_from = SI_name)
                      
                      si.txt <- ifelse(typ=="management","m","g")
                      x.txt <-bquote(bolditalic("W")[bold(.(si.txt)*',GY')]~bold('('*.(unit.pi)*')'))
                      y.txt <-bquote(bolditalic("P")[bold(.(si.txt)*',GY')]~bold('('*.(unit.pi)*')'))
                      if(typ=="management"){
                        
                        p0 <- subdata %>%
                          dplyr::mutate(g=Genotype,
                                        Genotype=gsub("_",'\n',Genotype)) %>% 
                          tidyr::separate(col = g,
                                          into=c('application','nitrogen',
                                                 'sowing_date'),
                                          sep = '[_]') %>% 
                          mutate(treatment=interaction(application,sowing_date,sep=" ")%>%
                                   factor(.,levels=c("Split Early","Split Late",
                                                     "Combined Early","Combined Late"))) %>% 
                          ggplot(.,
                                 aes(x=ecovalence,
                                     y=genotypic.superiority.measure,
                                     shape=treatment,
                                     color=nitrogen))+
                          scale_color_manual(values=col.nitro)+
                          geom_hline(yintercept = subdata$genotypic.superiority.measure %>% mean(),lty=2)+
                          geom_vline(xintercept = subdata$ecovalence %>% mean(),lty=2)+
                          ggrepel::geom_text_repel(aes(label=Genotype),size=2,show.legend = F)
                        
                        shp <- c("Split Early" = 1, "Combined Early" = 16,
                                 "Split Late"=0,"Combined Late"=15)
                      }else{
                        p0 <-  subdata %>%ggplot(.,
                                                 aes(x=ecovalence,
                                                     y=genotypic.superiority.measure,
                                                     color=Genotype,
                                                     shape=Genotype)) +
                          scale_color_manual(values = col.clr)+
                          geom_hline(yintercept = subdata$genotypic.superiority.measure %>% mean(),lty=2)+
                          geom_vline(xintercept = subdata$ecovalence %>% mean(),lty=2)+
                          ggrepel::geom_text_repel(aes(label=Genotype),size=2.5,show.legend = F)
                        shp <- shape_geno
                      }
                      
                      pl <- p0 +
                        geom_point(size=3,stroke=1.5,show.legend = F)+
                        scale_shape_manual(values = shp)+
                        xlab(x.txt)+
                        ylab(y.txt)+
                        theme_phd_facet(ax.txt.siz = 6,ax.tit.siz = 8,b=10,t=25)
                      
                      print(pl)
                      return(pl)
                    })


prow2 <-    cowplot::plot_grid(
  cowplot::plot_grid(
    plotlist = plist,align = 'hv',
    labels =c("A Genotype SI","B Management SI"),
    label_size=10,rel_widths = c(1,1),
    hjust = -1,nrow = 1),
  legend_a,legend_b,
  rel_heights = c(.9,.1,.1),ncol=1)

tiff(filename=paste0("result/Fig.1.tiff"),
     units="cm",
     width=17.4,
     height=12.5,
     compression = "lzw",
     pointsize=12,
     res=400)# dpi
print(prow2)
dev.off()

