rm(list=ls())
pacman::p_load(dplyr,purrr,toolPhD,ggplot2,toolStability,ggnewscale)

b0 <- read.csv("data/merge_profile.csv") %>% 
  dplyr::filter(namCombine%in%c("spike_number_87","tiller_31")) %>%
  dplyr::select(-c(Acc_Radiation,group,TT,Das,BBCH,part,trait)) %>% 
  tidyr::pivot_wider(values_from = "Trait",names_from = "namCombine") %>% 
  mutate(tiller_fertility=spike_number_87/tiller_31) %>% 
  select(timeid:appl,tiller_fertility) %>% 
  rename(Trait=tiller_fertility) %>% 
  mutate(namCombine="tiller_fertility")

b <- read.csv("data/merge_profile.csv") %>% 
  select(timeid:appl,Trait,namCombine) %>% 
  bind_rows(.,b0)

col.nitro <- c("#B3AD00","#1F9600")
names(col.nitro) <- c("176","220")
hex <- scales::hue_pal()(3) 
names(hex) <-c("2019","2020","2021")
Geno <- c("Capone","Pionier","Patras","Apertus","Torrild","Alves","Potenzial","Esket")
col.clr <- c("#003f5c","#58508d","#8a508f","#bc5090", "#de5a79", "#ff6361","#ff8531", "#ffa600")
names(col.clr) <- Geno
shape_geno <- c(0,1,2,5,6,8,12,13)
names(shape_geno) <- Geno

g_rank <- b %>% 
  dplyr::filter(namCombine=="plot_yield") %>% 
  genotypic_superiority_measure(.,'Trait',genotype = 'var',
                                environment = c('appl','DFG_year',
                                                'timeid','nitrogen'),
                                unit.correct = T) %>% 
  dplyr::select(-Mean.Trait)  %>%  
  arrange(genotypic.superiority.measure) %>% 
  rename('pg_rank'='genotypic.superiority.measure','var'='Genotype')
# -------------------------------------------------------------------------
glance.lm <- function(x, ...) {
  #https://www.tidymodels.org/learn/develop/broom/
  sx <- summary(x)
  result <- sx$coefficients %>%
    tibble::as_tibble(rownames = "term") %>%
    dplyr::rename(slope = Estimate,
                  p.value = `Pr(>|t|)`) %>% 
    cbind(.,
          with(sx,
               tibble::tibble(
                 r.squared = r.squared,
                 adj.r.squared = adj.r.squared))) %>% 
    .[,!grepl(" ",names(.))]
  return(result)  
}
sc_ppp <- function(xv,yv,xlabt=NULL,lgd=NULL,...){
  k <- b %>% 
    filter(namCombine%in%c(xv,yv),
           !DFG_year=="DFG2019"
    ) %>% 
    group_by(var,namCombine) %>% 
    summarise(Trait=mean(Trait,na.rm=T),.groups = "drop") %>% 
    # change the nitrogen from t/ha kg/ha
    mutate(Trait=case_when(grepl("_(N|N87)\\b",namCombine)~Trait*1000,
                           T~Trait)) %>% 
    tidyr::pivot_wider(names_from="namCombine",values_from="Trait") %>% 
    tidyr::pivot_longer(-c(xv,"var"),values_to = "Trait",names_to = "trait")
  
  lmdf<-
    k %>%
    group_by(trait) %>% 
    tidyr::nest() %>% 
    mutate(model = map(data, ~sprintf("Trait~%s",xv) %>% 
                         as.formula() %>% lm(.,data=.x) %>% 
                         glance.lm(.))) %>% 
    tidyr::unnest(model) %>% 
    dplyr::select(-data) %>% 
    filter(term==xv) %>% 
    mutate(sign=ifelse(slope>0,1,-1)) %>% 
    arrange(desc(sign),desc(r.squared)) %>% 
    left_join(look.up) %>% ungroup() %>% 
    mutate(n=1:n(),
           fn=paste(LETTERS[n],"~",fn) %>% 
             map_chr(.,~paste0('bold(',.,')')))
  
  p1 <-
    k %>% 
    merge(.,lmdf %>% select(trait,sign,fn,p.value)) %>% 
    # merge(.,look.up) %>% 
    mutate(fn=factor(fn,levels=lmdf$fn),
           var=factor(var,levels=Geno)) %>% 
    ggplot(aes(.data[[xv]],Trait))+
    toolPhD::theme_phd_facet(b=10,
                             strp.txt.siz = 9)+
    
    theme(axis.title.y=element_blank(),
          legend.key.height = unit(1, "cm"),
          ...)+
    geom_point(aes(shape=var,color=var),
               stroke =2,alpha=.5,size=2.5)+
    scale_color_manual(values=col.clr)+
    scale_shape_manual(values=shape_geno)+
    ggpmisc::stat_poly_line(
      data = . %>% group_by(fn) %>% filter(any(p.value<0.05)),
      formula = y ~ x,
      aes(group = fn),
      se=F,linewidth=.5) +
    ggpmisc::stat_poly_eq(
      data = . %>% group_by(fn) %>% filter(any(p.value<0.05)),
      formula =  y ~x,
      aes(
        group=fn,
        label = paste(
          after_stat(rr.label),
          sep = "*\", \"*")),
      label.x = .9,
      label.y=.1,
      size = 2.8)+
    facet_wrap(~fn,ncol = 3,labeller = label_parsed,
               scales = "free_y")+
    guides(
      shape=guide_legend(title="genotype",
      ),
      color=guide_legend(title="genotype",
      ))
  if(is.null(lgd)){
    # put it on the lower right corner
    p1 <- p1+theme(
      legend.direction ="horizontal",
      
      legend.position=   c(.85,.1)
    )+
      guides(
        shape=guide_legend(ncol=3,
                           title="genotype",
                           direction = "horizontal",
                           title.position = "top"),
        color=guide_legend(ncol=3,
                           title="genotype",
                           direction = "horizontal",
                           title.position = "top"))
    
  }
  if(!is.null(xlabt)){
    p1 <- p1+
      xlab( parse(text=xlabt))
  }
  print(p1)
}

yv<- c("mobile61_WSC",  # negative #21
       "mobile61_N", # negative, w21
       "mobile61_DM",# negative
       "kernel_total_N87", # positive , w21
       "straw_DM87", # positive
       "straw_CHO87", # positive    
       "post61_DM", # positive
       "plot_yield"
)

look.up <- data.frame(
  trait=c("mobile61_WSC",  # negative #21
          "mobile61_N", # negative, w21
          "mobile61_DM",# negative
          "kernel_total_N87", # positive , w21
          "straw_DM87", # positive
          "straw_CHO87", # positive    
          "post61_DM", # positive
          "plot_yield"
  ),
  fn= c(
    'Delta*WSC["61-87, straw"]~"("*t/ha*")"',
    'Delta*N["61-87, straw"]~"("*kg/ha*")"',
    'Delta*DM["61-87, straw"]~"("*t/ha*")"',
    'N["87, grain"]~"("*kg/ha*")"',
    'DM["87, straw"]~"("*t/ha*")"',
    'WSC["87, straw"]~"("*mg/g*")"',
    'Delta*DM["87-61, all"]~"("*t/ha*")"',
    'GY~"("*t/ha*")"'
  ) %>% map_chr(.,~{paste0('bold(',.x,')')})
  
)

# -------------------------------------------------------------------------
look.up <- data.frame(
  trait=c("mobile61_WSC",  # negative #21
          "mobile61_N",
          "kernel_Nitrogen87",
          "straw_CHO87", # negative, w21
          "GCD61",
          "S61_87"),
  
  fn=c('Delta*WSC["61-87,straw"]~"("*t~ha^-1*")"',
       'Delta*N["61-87,straw"]~"("*kg~ha^-1*")"',
       '"["*N*"]"["87,grain"]~"("*mg~g^-1*")"',
       '"["*WSC*"]"["87,straw"]~"("*mg~g^-1*")"',
       'GCD~"("*degree*C*d*")"',
       'TT["61-87"]~"("*degree*C*d*")"')
  
) 

tiff(filename='result/plot/Fig4.tiff',
     units="cm",
     width=21,#21.6
     height=15,#16.1
     compression = "lzw",
     pointsize=3,
     res=500,# dpi)
     family="Arial")
sc_ppp("post61_DM",look.up$trait,
       'Delta*DM["87-61,all"]~"("*t~ha^-1*")"',lgd="right") %>% print()
dev.off()

# -------------------------------------------------------------------------
hex <- scales::hue_pal()(3) 
names(hex) <-c("2019","2020","2021")
ltype <- c("solid","solid","dashed")
names(ltype) <- c("2019","2020","2021")
glk <- data.frame(var=g_rank$var) %>% 
  mutate(n=1:n(),
         genotype=paste0(LETTERS[n]," ",var)
  )

pp<- function(yv,xv){
  # individual year
  k <-  b %>% filter(namCombine%in%c(xv,yv,"plot_yield")) %>% 
    tidyr::pivot_wider(names_from="namCombine",values_from="Trait") %>% 
    tidyr::pivot_longer(-c(xv,"var","DFG_year","timeid","appl","nitrogen"),
                        values_to = "Trait",names_to = "trait")%>% 
    filter(trait==yv) 
  
  lmdf<-
    k %>%
    group_by(DFG_year,var) %>% 
    tidyr::nest() %>% 
    mutate(model = map(data, ~sprintf("Trait~%s",xv) %>% 
                         as.formula() %>% lm(.,data=.x) %>% 
                         glance.lm(.))) %>% 
    tidyr::unnest(model) %>% 
    dplyr::select(-data) %>% 
    filter(term==xv) %>% 
    mutate(sign=ifelse(slope>0,1,-1)) %>% 
    arrange(desc(sign),desc(r.squared))
  
  
  st <- k%>% 
    dplyr::mutate(treatment=interaction(appl,timeid,sep=" ")%>%
                    factor(.,levels=c("Split Early","Split Late",
                                      "Combined Early","Combined Late")),
                  nitrogen=as.factor(nitrogen),
                  year=gsub("DFG","",DFG_year)
    ) %>% 
    merge(.,lmdf) %>% 
    merge(.,glk)
  
  p <-  st %>% 
    ggplot(aes(.data[[xv]],Trait))+
    toolPhD::theme_phd_facet(legend.position="bottom")+
    geom_point(aes(
      color=year,
    ),shape=1,
    size=2)+
    
    ggpmisc::stat_poly_line(
      data=st,
      formula = y ~ x,
      mapping=aes(group = year,linetype=year,color=year),
      se=F,linewidth=.5) +
    scale_color_manual(values=hex)+
    scale_linetype_manual(values=ltype)+
    ggpmisc::stat_poly_eq(
      data=. %>% group_by(var,DFG_year) %>% filter(any(p.value<0.05)),
      formula =  y ~x,
      aes(
        color=year,
        group = year,
        label = paste(
          after_stat(rr.label),
          sep = "*\", \"*")),
      size = 2.8)+
    ggh4x::facet_nested_wrap(
      
      genotype~.,ncol=4) +
    xlab(parse(text='bold(GCD~"("*degree*C*d*")")'))+
    ylab(parse(text='bold("["*WSC*"]"["87, straw"]~"("*mg~g^-1*")")'))+
    theme(strip.placement = "right")
  print(p)
}

# -------------------------------------------------------------------------
tiff(filename='result/plot/Fig5.tiff',
     units="cm",
     width=18,#21.6
     height=12,#16.1
     compression = "lzw",
     pointsize=3,
     res=500,# dpi)
     family="Arial")
pp("straw_CHO87","GCD61") %>% print()
dev.off()

