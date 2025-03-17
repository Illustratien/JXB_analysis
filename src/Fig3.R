rm(list=ls())
options(dplyr.summarise.inform = FALSE)
pacman::p_load(toolStability,dplyr,toolPhD,foreach,ggplot2,cowplot,purrr)

hex <- c("#f6ab49","#97383c","#4b87a1")

ldf <- data.frame(
  namCombine=c("ear_DM61",'mobile61_DM','post61_DM'),
  fn=c('DM["61,spike"]',
       'Delta*DM["61-87,straw"]',
       'Delta*DM["87-61,all"]'))
names(hex) <- ldf$fn

# -------------------------------------------------------------------------
datm <- read.csv('data/merge_profile.csv')%>%
  filter(!DFG_year=='DFG2019')

g_rank <- datm %>% 
  dplyr::filter(trait=='yield',part=='plot') %>% 
  genotypic_superiority_measure(.,'Trait',genotype = 'var',
                                environment = c('appl','DFG_year',
                                                'timeid','nitrogen'),
                                unit.correct = T) %>% 
  dplyr::select(-Mean.Trait)  %>%  
  arrange(genotypic.superiority.measure) %>% 
  rename('pg_rank'='genotypic.superiority.measure','var'='Genotype')

hsd_tbl <- read.csv("result/HSDtable.csv") %>% 
  dplyr::select(groups,trait,treatment,Year,combi) %>% 
  rename(namCombine=trait,DFG_year=Year,var=treatment) %>% 
  filter(namCombine%in%c('post61_DM','mobile61_DM',"ear_DM61"),
         combi=="var") %>%
  dplyr::select(-combi) 

# grain source filling -------------------------------------------------------------------------

filling_strategy <- read.csv('data/merge_profile.csv')%>% 
  filter(namCombine%in%c('post61_DM','mobile61_DM',"ear_DM61")) %>% 
  group_by(var,DFG_year,namCombine) %>% 
  summarise(Tra=mean(Trait,na.rm=T),Trasd=sd(Trait,na.rm=T)) %>%
  left_join(.,hsd_tbl) %>% 
  rename(Genotype=var) %>% 
  merge(.,ldf) %>% 
  mutate(namCombine=factor(fn,levels=ldf$fn %>% rev()),
         Genotype=factor(Genotype,
                         levels=g_rank$var%>% 
                           rev() ),
         l=case_when(is.na(groups)~paste0('bold(',round(Tra,1),')'),
                     T~ paste0('bold(',round(Tra,1),'^',groups,')')),
         Y=case_when(DFG_year=="DFG2019"~"(A) 2019",
                     DFG_year=="DFG2020"~"(B) 2020",
                     T~"(C) 2021"))

p2 <- ggplot(filling_strategy,
             aes(y=Tra,fill=namCombine, x=Genotype,label=l)) + 
  geom_bar(
    stat="identity"
  )+
  theme_phd_facet(b=3,t=3,
                  lgd.tit.siz = 10,lgd.txt.siz = 10)+ # spacing of x and axis)+
  guides(fill=guide_legend(title="grain filling source"))+
  facet_grid(~Y)+
  scale_fill_manual(values=hex,
 
                    labels = scales::parse_format()
  )+
  ylab(parse(text='bold(Contribution~to~DM["87,spike"]~(t/ha))'))+

  geom_text(
    size = 3.6, parse=T,position = position_stack(vjust = 0.5),
    color="white",fontface="bold")+
  theme(
    legend.justification='left',
    legend.key.height = unit(1, "cm"),
    axis.title.x = element_text(margin = margin(b=8),vjust=-5))+
  coord_flip()+
  xlab("cultivar")

# -------------------------------------------------------------------------

tr <- c('mobile_effciency_WSC','post61_DM',"plot_yield","GCD61")

look.vec <- c('"[WSC]"["straw"]~frac("["*61-87*"]","["*61*"]")',
              'Delta*DM["87-61,all"]~"("*t/ha*")"','GY~"("*t/ha*")"',
              'GCD61~"("*d*degree*C*")"'
) %>% map_chr(.,~{paste0('bold(',.x,')')})

v <- read.csv("data/merge_profile.csv") %>% 
  filter(namCombine%in%tr,
         ) %>% 
  group_by(DFG_year,var,namCombine) %>% 
  summarise(Trait=mean(Trait)) %>% 
  tidyr::pivot_wider(names_from = 'namCombine',values_from = 'Trait') %>% 
  mutate(post61DM_per_mobile_WSC=post61_DM/mobile_effciency_WSC)

sv <- v %>% 
  group_by(var) %>% 
  dplyr::summarise(across(where(is.numeric),~mean(.x,na.rm = T)))

xy_ls<- list(c(2,3))

dr <- bind_rows(v,sv) %>% 
  ungroup() %>% 
  select(-DFG_year) %>%
  tidyr::pivot_longer(-var,names_to = "trait",values_to = "Trait") %>% 
  group_by(trait) %>% 
  summarise(min=min(Trait)*.98,max=max(Trait)*1.05)

rangedf <- data.frame(t(dr[,-1])) %>% 
  setNames(., dr[,1] %>% unlist())

# -------------------------------------------------------------------------
hex <- scales::hue_pal()(3) 

names(hex) <-c("2019","2020","2021")
ltype <- c("dotted","solid","dashed")
names(ltype) <- c("2019","2020","2021")


qp4<- function(datta,xy_vec){
  tr_vec <- tr[xy_vec]
  xvar <- tr_vec[1] %>% sym()
  yvar <- tr_vec[2]%>% sym()
  
  datta %>%
    mutate(
      Year=gsub("DFG","",DFG_year),
      Year=factor(Year,levels=c("2019","2020","2021"))) %>%
    ggplot(aes(x=!!xvar,y=!!yvar))+
    theme_phd_facet(b=3,t=3,l=20,
                    lgd.tit.siz = 10,lgd.txt.siz = 10)+#ax.txt.siz =6,ax.tit.siz = 6,l=0,b=2,t=1,r=1
    scale_y_continuous(limits = rangedf[[tr_vec[2]]],labels = scales::label_number())+
    scale_x_continuous(limits = rangedf[[tr_vec[1]]],labels = scales::label_number())+
    ggpmisc::stat_poly_line(formula = y ~ x,
                            aes(group = Year,linetype=Year,color=Year),
                            # color="darkgray",
                            se=F,linewidth=.5) +
    ggpmisc::stat_poly_eq(formula =  y ~x,
                          aes(
                            group = Year,color=Year,
                            label = paste(
                              after_stat(eq.label),
                              after_stat(rr.label),
                              sep = "*\", \"*")),
                          label.x = c("left","right","right"),
                          label.y = c("top","bottom","bottom"),
                          size = 2.8)+
    scale_color_manual(values=hex)+
    scale_linetype_manual(values=ltype)+
    labs(x=parse(text=look.vec[xy_vec[1]]),
         y=parse(text=look.vec[xy_vec[2]]))+
    geom_point(size=4,shape=1,aes(color=Year))+
    facet_grid(~source)+
    ggrepel::geom_text_repel(aes(label=var),size=2.7,show.legend = F)+
    theme(axis.title.x = element_text(margin = margin(b=8),vjust=-2))
}


tiff(filename='result/Fig3_post.tiff',
     units="cm",
     width=21,#21.6
     height=18,#16.1
     compression = "lzw",
     pointsize=3,
     res=500,# dpi)
     family="Arial")
cowplot::plot_grid(
  p2+theme(legend.position = "none"),
  cowplot::plot_grid(qp4(
    v %>% mutate(source="(D) post-anthesis assimilates"),
    xy_ls[[1]]),
    get_plot_component(p2,pattern = 'guide-box-right',return_all = T) %>% 
      ggdraw(),nrow=1,align="h",rel_widths = c(1,.65)),
  align="h",ncol=1)
dev.off()
