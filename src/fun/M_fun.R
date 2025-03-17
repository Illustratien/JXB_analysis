# fun ---------------------------------------------------------------------
# round2<- function(vec){
#   map_dbl(vec,function(x){
#     if(x<=.1){
#       round(x,3)
#     }else if(.1<x&x<1){
#       round(x,2)
#     }else{
#       round(x,1)
#     }
#   }) %>% as.character()
# }
table_fun1 <- function(df,x){
  # input:
  # df: data frame
  # x: character, specifiying the group for summary
  #    could be either "var" or "management"
  # 
  # output:
  # formated summarised data frame 
  if(x=='g'){
    subhsd <- hsd_tbl %>% 
      mutate(treatment=gsub("_","\n",treatment))
  }else{
    subhsd <- hsd_tbl
  }
  subhsd <-subhsd %>% rename(!!quo_name(x):="treatment")
  
  
  res <- df %>% 
    group_by_at(c(x,"trait_name","DFG_year")) %>% 
    summarise(
      # Trait=mean(Trait,na.rm=T) %>% round2(.)
      Trait=mean(Trait,na.rm=T) %>% toolPhD::round_scale(.),
      # trait_sd=sd(Trait,na.rm=T) %>% toolPhD::round_scale(.)
    ) %>% 
    left_join(.,subhsd, c(x, "trait_name", "DFG_year")) %>% 
    # tidyr::unite("Trait",c(trait_m,trait_sd),sep='±') %>%
    mutate(Trait=case_when(!is.na(groups)~paste(Trait,groups,sep="^"),
                           T~Trait)) %>% 
    tidyr::pivot_wider(.,names_from = trait_name,values_from = "Trait") 
  return(res)
}
table_fun2 <- function(df,x){
  # input:
  # df: data frame
  # x: character, specifiying the group for summary
  #    could be either "var" or "management"
  # 
  # output:
  # formated summarised data frame 
  if(x=='g'){
    subhsd <- hsd_tbl %>% 
      mutate(treatment=gsub("_","\n",treatment))
  }else{
    subhsd <- hsd_tbl
  }
  subhsd <-subhsd %>% rename(!!quo_name(x):="treatment")
  
  
  res <- df %>% 
    group_by_at(c(x,"trait_name","DFG_year")) %>% 
    summarise(
      # Trait=mean(Trait,na.rm=T) %>% round2(.)
      trait_m=mean(Trait,na.rm=T) %>% toolPhD::round_scale(.),
      trait_sd=sd(Trait,na.rm=T) %>% toolPhD::round_scale(.)
    ) %>% 
    left_join(.,subhsd, c(x, "trait_name", "DFG_year")) %>% 
    tidyr::unite("Trait",c(trait_m,trait_sd),sep='±') %>%
    mutate(Trait=case_when(!is.na(groups)~paste(Trait,groups,sep=" "),
                           T~Trait)) %>% 
    tidyr::pivot_wider(.,names_from = trait_name,values_from = "Trait") 
  return(res)
}

get_mean <- function(vec){
  # transform column of output formated text 
  # into plot ready numeric
  gsub("±.*","",vec) %>% as.numeric()
}
get_sd <- function(vec){
  # transform column of output formated text 
  # into plot ready numeric
  gsub(".+?±","",vec) %>%
    gsub("[a-z]+?","",.) %>%  as.numeric()
}

layout_para<- function(nplot,  N.per.page=4){
  
  # for the display purpose,
  # calculate the position where legend should be shown
  # nplot: number of plots in total,
  # N.per.page: how many plot to dsiplay per page.
  
  if(nplot>N.per.page){
    n_iter <- nplot%/%N.per.page
    
    pind <- map(1:n_iter,~{seq(N.per.page*(.x-1)+1,.x*N.per.page)})
    
    if(nplot%%N.per.page>0){
      pind[[n_iter+1]] <- seq(N.per.page*n_iter+1,
                              N.per.page*n_iter+nplot%%N.per.page)
      n_iter <- n_iter+1
    }
    lgd_pos_ind <- pind %>% map_dbl(.,~{.x[length(.x)]})
  }else{
    lgd_pos_ind <-nplot
  }
  return(list(lgd_pos_ind,n_iter,pind))
}

condi1 <- rlang::quo(trait%in%c("CHO","total_C","Nitrogen","total_N")&BBCH%in%c(61,87))

condi2 <- rlang::quo(trait%in%(trait[grep('mob|pre|post|TE$',trait)])|
                       (trait_name%in%trait_name[grep("ear_biomass61",trait_name)]))

condi3 <- rlang::quo((trait%in%(trait[grep('tkw|grain_number|yield|grain_per_spike|grain_(width|length)',trait)])&part=="plot")|
                       (trait_name%in%trait_name[grep("spike_number|biomass_ear|FE",trait_name)]))
condi4 <- rlang::quo((trait%in%(trait[grep('GCD61|S$|RUE|gL(A)?I_int$|gSPAD_int|(^g.*61$)',trait)])|(trait%in%trait[grep("Delta",trait)]&BBCH=="61_87")))

condi5 <- rlang::quo((trait_name%in%( unique(hsd_tbl$trait_name))))

condi6 <- rlang::quo(trait_name%in%c("post61_DM","pre61_DM","mobile61_DM",
              "ear61_DM_rel" ,"ear_DM61","post61_DM_rel","pre61_DM_rel","mobile61_DM_rel",
              "plot_yield"))

condilist <- list(condi1,condi2,condi3,condi4,condi6,condi5)
condiName<- c("nutri","translocation","ycomponent","postanthesis","grainfill",'sig_trait')

# eval <- function(con){
#   sub_df %>% dplyr::filter(!!con) %>%
#     arrange(BBCH,part,trait) %>% 
#     .$trait_name %>% unique()
# }
# merge_dat$trait_name %>% unique() %>% .[grep("mobile",.)]
# -------------------------------------------------------------------------
con1 <- rlang::quo(grepl("DM",trait)&BBCH%in%c(51,61,87)|
                     grepl("total_C|WSC|total_N|Nitrogen|CHO",trait)&BBCH%in%c(51,61,87)|
                     grepl("mobile61|(pre|post)",trait_name))

con2 <- rlang::quo(grepl("FE|spike|grain|tkw|plot_yield",trait_name)| grepl("tiller",trait)&BBCH%in%c(51,61,87))

con3<- rlang::quo((trait%in%(trait[grep('GCD61|S$|RUE|gL(A)?I_int$|gSPAD_int|(^g.*61$)',trait)])|
                     (trait%in%trait[grep("Delta",trait)]&BBCH%in%c("61_87","51_61"))))
condils <- list(con1,con2,con3,condi6,condi5)
condiName<- c("DM_trans","ycomponent","postanthesis","grainfill",'sig_trait')

# 
# rlang::quo(trait%in%(trait[grep('mob|post',trait)])|
#              (trait_name%in%trait_name[grep("(ear_biomass61|ear_total_C61|TE)",trait_name)]))
# eval(con2)

range01 <- function(x){
  funs <- c(min,max)
  res<- lapply(funs, function(f) f(x, na.rm = TRUE))
  (x-res[[1]])/(res[[2]]-res[[1]])
}

gg_color_hue <- function(n) {
  # https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

table_fun3 <- function(df,x){
  # input:
  # df: data frame
  # x: character, specifiying the group for summary
  #    could be either "var" or "management"
  # 
  # output:
  # formated summarised data frame 
  # subhsd <-subhsd %>% rename(!!quo_name(x):="treatment")
  
  
  res <- df %>% 
    group_by_at(c(x,"cl","trait_name","DFG_year")) %>% 
    summarise(
      Trait=mean(Trait,na.rm=T) %>% toolPhD::round_scale(.),
      # trait_m=mean(Trait,na.rm=T) %>% round2(.),
      # trait_sd=sd(Trait,na.rm=T) %>% round2(.)
      ) %>% 
    left_join(.,hsd_sub, c(x, "trait_name", "DFG_year")) %>% 
    # tidyr::unite("Trait",c(trait_m,trait_sd),sep='±') %>%
    mutate(Trait=case_when(!is.na(groups)~paste(Trait,groups,sep=" "),
                           T~Trait))
  return(res)
}

have_common_characters <- function(string1, string2) {
  common_chars <- intersect(strsplit(string1, "")[[1]], strsplit(string2, "")[[1]])
  return(length(common_chars) > 0)
}

overlap<- function(x,y){
  map2_lgl(x,y,~{
    if(any(is.na(c(.x,.y)))){
      NA
    }else{ 
      have_common_characters(.x,.y)
    }
    
  })
}
shapev <- c(0,2,1)
names(shapev) <- paste0("DFG",2019:2021)
plot_fun <- function(df){
  p <- df %>% 
    ggplot(aes(trait_name,r,color=var))+
    geom_point(aes(shape=DFG_year))+
    ggrepel::geom_text_repel(data=df %>% filter(d2==T),
                             mapping=aes(label=Trait),show.legend = F,
                             # nudge_y=.04,
                             size=3,
                             box.padding = .05)+
    scale_shape_manual(values=shapev)+
    geom_hline(yintercept = .9,linetype=2,alpha=.3,color='darkgrey')+
    geom_hline(yintercept = 1.1,linetype=2,alpha=.3,color='darkgrey')+
    geom_hline(yintercept = 1,linetype=2,alpha=.5,color='darkgrey')+
    theme_bw()+
    scale_color_manual(values=colors_border)+
    # ggh4x::facet_grid2(DFG_year~.,scales = "free_x",space="free")+
    ylab(paste0(var_tar,collapse="/"))+
    theme(legend.position = c(.85,.08),
          legend.direction = "horizontal",
          strip.background = element_blank(),
          panel.grid = element_blank(),
          legend.background = element_blank())
  return(p)
}
y <- merge_dat %>% dplyr::filter(namCombine=="plot_yield")

g_rank <- genotypic_superiority_measure(y,'Trait',genotype = 'var',
                                        environment = c('appl','DFG_year',
                                                        'timeid','nitrogen'),
                                        unit.correct = T) %>% 
  dplyr::select(-Mean.Trait)  %>%  arrange(genotypic.superiority.measure) %>% 
  rename('pg_rank'='genotypic.superiority.measure','var'='Genotype')

m_rank <-genotypic_superiority_measure(y,'Trait',genotype = 'g',
                                       environment = c('DFG_year',
                                                       'var'),
                                       unit.correct = T) %>% 
  
  dplyr::select(-Mean.Trait)  %>%  
  arrange(genotypic.superiority.measure) %>% 
  rename('pm_rank'='genotypic.superiority.measure','g'='Genotype')
