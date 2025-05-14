rm(list=ls())
pacman::p_load(toolStability,ggplot2,dplyr)
# color setting -------------------------------------------------------------------------
Geno  <- c("Capone","Pionier","Patras","Apertus","Torrild","Alves","Potenzial","Esket")
col.clr <- c("#003f5c","#58508d","#8a508f","#bc5090", "#de5a79", "#ff6361","#ff8531", "#ffa600")
names(col.clr) <- Geno
# shapes for genotype
shapes <- c(0,1,2,5,6,8,12,13)
names(shapes) <- Geno

# for single figure
th_lgd <- theme(legend.title = element_text(size = 6,face='bold'),# legend size
                legend.key.height = unit(.6,"cm"),
                legend.key.width = unit(.6,"cm"),
                legend.text = element_text(size =6,face='bold',
                                           margin = margin(t = 6,b=6)))
th_rm <- theme(legend.position = "none")
pi.vec <- c("pi","ecovalence")
# add straw -------------------------------------------------------------------------

datm <- read.csv('data/merge_profile.csv') %>%
  filter(!DFG_year=='DFG2019')
g_df <- datm %>% 
  dplyr::filter(trait=='yield',part=='plot') %>% 
  genotypic_superiority_measure(.,'Trait',genotype = 'var',
                                environment = c('appl','DFG_year',
                                                'timeid','nitrogen'),
                                unit.correct = T) %>% 
  dplyr::select(-Mean.Trait)  %>%  
  arrange(genotypic.superiority.measure) %>% 
  rename('pg_rank'='genotypic.superiority.measure','var'='Genotype')
trait.vec<- c("plot_yield","plot_tkw",
              "plot_grain_number","HI",
              "spike_number_87","plot_grain_per_spike",
              "grain_length","grain_width",
              "straw_DM87",
              "post61_DM","mobile61_DM")

stable_table <- datm %>% 
  dplyr::filter(namCombine%in%trait.vec) %>%
  group_by(trait) %>% group_split() %>% 
  purrr::map_dfr(.,~{
    genotypic_superiority_measure(.x,'Trait',genotype = 'var',
                                  environment = c('appl','DFG_year',
                                                  'timeid','nitrogen'),
                                  unit.correct = T) %>% 
      dplyr::select(-Mean.Trait) %>% 
      mutate(trait=.x$namCombine[1])
  }) %>%
  rename('pg_rank'='genotypic.superiority.measure','var'='Genotype')

pgrank <- stable_table %>% 
  tidyr::pivot_wider(names_from = trait,values_from = pg_rank) %>% 
  dplyr::select(-c(var)) %>% 
  relocate("plot_yield") 

# matrix of the p-value of the correlation
p.mat <- corrplot::cor.mtest(pgrank,
                             method = "spearman")$p


M <- pgrank%>% 
  cor(.,method = "spearman")

lk <- data.frame(ori=colnames(M),
                 new=c("GY",'DM["87,straw"]',"HI","GL","GN","GpS","GW",'DM["61-87,straw"]',
                       'DM["87-61,all"]',"SN","TGW"
                 )) %>% 
  mutate(bn=paste0('bold(italic(P)[g*","~',new,"])"))

colnames(M) <-
  with(lk,bn[match(colnames(M),ori)])
rownames(M) <- with(lk,bn[match(rownames(M),ori)])

colnames(p.mat) <-
  with(lk,bn[match(colnames(p.mat),ori)])
rownames(p.mat) <- with(lk,bn[match(rownames(p.mat),ori)])
# single plot -------------------------------------------------------------
si.vec <- "genotypic.superiority.measure"
udf <- read.csv('data/unit.csv') %>% 
  mutate(
    unit=gsub('deg',"°",unit) %>%
      gsub("(\\{|\\})","",.) %>%
      gsub("\\*","x",.))

unit.df<- purrr::map_dfr(trait.vec,~{
  id <- .x %>% gsub("(plot_|_87)?", "", .)%>% 
    grepl(.,udf$trait,perl = T)%>% which()
  udf[id,] %>% mutate(namCombine=.x)
})

# get corresponding unit based on trait
unit0 <- unit.df$unit %>% 
  gsub('(\\/ha|\\/ ha)',"~ha^-1",.) %>% 
  gsub('(tha|t ha)','t~ha',.) %>% 
  gsub('\\/ 1000 grain',"1000~grain~^-1",.) %>% 
  gsub(' ','~',.)
# unit0
lsp.vec <- rep(5,length(unit0))

# -------------------------------------------------------------------------
stable_table<- datm %>% 
  dplyr::filter(namCombine%in%trait.vec) %>%
  group_by(trait) %>% group_split() %>% 
  purrr::map_dfr(.,~{
    genotypic_superiority_measure(.x,'Trait',genotype = 'var',
                                  environment = c('appl','DFG_year',
                                                  'timeid','nitrogen'),
                                  unit.correct = T) %>% 
      mutate(trait=.x$namCombine[1]) %>% 
      rename(value=Mean.Trait)
  }) %>%
  rename('SI_value'='genotypic.superiority.measure')


pll2 <- purrr::map(si.vec,function(si.nam){
  # and for each trait
  pl <- purrr::map(1:3,function(i){
    # length(trait.vec)
    trait.name <- trait.vec[i]
    lsp <- lsp.vec[i]
    trnam <-   with(lk,new[match(trait.vec[i],ori)])
    xlab1<- bquote(bold(.(trnam)~.(unit0[[i]])))
    
    Si.txt <- ifelse(si.nam=="Ecovalence","W","P")
    ylab1 <- bquote(bolditalic(.(Si.txt))[bold("g, "*.(trnam))]~bold(.(unit0[[i]])))
    if(grepl("\\^",unit0[[i]])){
      ylab1 <- bquote(bolditalic(.(Si.txt))[bold("g, "*.(trnam))]~bold("x"*"("*"10"^"6"~"ha"^-1*")"))
      xlab1<- bquote(bold(.(trnam)~bold("x"*"("*"10"^"6"~"ha"^-1*")")))
    }
    
    sub <- dplyr::filter(stable_table,
                         trait==trait.name)%>% 
      mutate(Genotype=factor(Genotype,levels=g_df$var))
    p <- ggplot(sub,
                aes(x     =value,y=SI_value,
                    color =Genotype,
                    shape =Genotype,label=Genotype))+
      toolPhD::theme_phd_facet(ax.tit.siz = 8,ax.txt.siz = 8,l = 5,t=10,b=3)+
      geom_point(size=1,alpha=.8,stroke = 1.2)+
      scale_color_manual(name='Genotype',values=col.clr)+
      scale_shape_manual(name = "Genotype",values=shapes)+
      scale_x_continuous(limits=c(range(sub$value)[1]*.98,
                                  range(sub$value)[2]*1.02))+
      scale_y_continuous(limits=c(range(sub$SI_value)[1]*.96,
                                  range(sub$SI_value)[2]*1.1))+
      ggrepel::geom_text_repel(seed = 195,show.legend = F,size=3,nudge_x = .1,box.padding =.3 )+
      xlab(xlab1)+
      ylab(ylab1)
    return(p)
    
  })
  names(pl) <- trait.vec[1:3]
  return(pl)
})


p <- ggcorrplot::ggcorrplot(M, 
                            type="lower",lab_col = "white",lab_size=1,
                            p.mat = p.mat, sig.level = 0.05,
                            insig = "pch",tl.cex =5,legend.title = parse(text="Spearman~italic(r)"),
                            color=c( "#CB7000","white", "#007878"),
                            lab = TRUE)+
  theme(legend.title     = element_text(size=7),
        legend.text      = element_text(size=7),
        legend.key.width = unit(dev.size()[1] / 100, "inches"))+
  scale_x_discrete(labels = scales::parse_format())+
  scale_y_discrete(labels = scales::parse_format())

# -------------------------------------------------------------------------
stable_table2 <- datm %>% 
  dplyr::filter(namCombine%in%trait.vec) %>%
  mutate(M=interaction(appl,timeid,nitrogen)) %>% 
  group_by(trait) %>% group_split() %>% 
  purrr::map_dfr(.,~{
    genotypic_superiority_measure(as.data.frame(.x),'Trait',genotype ="M" ,
                                  environment = c('var','DFG_year'),
                                  unit.correct = T) %>% 
      dplyr::select(-Mean.Trait) %>% 
      mutate(trait=.x$namCombine[1])
  }) %>%
  rename('pg_rank'='genotypic.superiority.measure','var'='Genotype')

pgrank2 <- stable_table2 %>% 
  tidyr::pivot_wider(names_from = trait,values_from = pg_rank) %>% 
  dplyr::select(-c(var)) %>% 
  relocate("plot_yield") 

# matrix of the p-value of the correlation
p.mat2 <- corrplot::cor.mtest(pgrank2,
                              method = "spearman")$p

M2 <- pgrank2%>% 
  cor(.,method = "spearman")

lk2 <- data.frame(ori=colnames(M2),
                  new=c("GY",'DM["87,straw"]',"HI","GL","GN","GpS","GW",'DM["61-87,straw"]',
                        'DM["87-61,all"]',"SN","TGW"
                  )) %>% 
  mutate(bn=paste0('bold(italic(P)[m*","~',new,"])"))

colnames(M2) <-
  with(lk2,bn[match(colnames(M2),ori)])
rownames(M2) <- with(lk2,bn[match(rownames(M2),ori)])
colnames(p.mat2) <-
  with(lk2,bn[match(colnames(p.mat2),ori)])
rownames(p.mat2) <- with(lk2,bn[match(rownames(p.mat2),ori)])
p2 <- ggcorrplot::ggcorrplot(M2, 
                             type="lower",lab_col = "white",lab_size=1,
                             p.mat = p.mat2, sig.level = 0.05,
                             insig = "pch",tl.cex =5,
                             legend.title = parse(text="Spearman~italic(r)"),
                             color=c( "#CB7000","white", "#007878"),
                             lab = TRUE)+
  theme(legend.title     = element_text(size=7),
        legend.text      = element_text(size=7),
        legend.key.width = unit(dev.size()[1] / 100, "inches"))+
  scale_x_discrete(labels = scales::parse_format())+
  scale_y_discrete(labels = scales::parse_format())

lgd2 <- cowplot::get_legend(p2)

lgd <- cowplot::get_legend(pll2[[1]][["plot_yield"]]+
                             th_lgd+
                             guides(color=guide_legend(ncol=2,
                                                       title.position = "top",
                                                       override.aes = list(size = 2)),
                                    shape=guide_legend(ncol=2,
                                                       title.position = "top",
                                                       override.aes = list(size = 2))))

# -------------------------------------------------------------------------
p3 <- cowplot::plot_grid(pll2[[1]][["plot_yield"]]+th_rm,
                         pll2[[1]][["plot_tkw"]]+th_rm,
                         pll2[[1]][["plot_grain_number"]]+th_rm,
                         p+th_rm,p2+th_rm, 
                         cowplot::plot_grid(lgd2,lgd,ncol=2,rel_widths = c(.3,1),
                                            align="v"
                         ),
                         ncol = 3,labels = c("A","B","C","D","E",""),
                         rel_widths = c(1,1,1,
                                        1,1,.5),
                         align = "hv", axis="lrtb",
                         label_x = .05,label_y = .99,
                         label_size = 10)

tiff(filename='result/plot/Fig2.tiff',
     units="cm",
     width=20,#21.6
     height=14,#16.1
     compression = "lzw",
     pointsize=12,
     res=500,# dpi)
     family="Arial")

p3 %>% print()
dev.off()
