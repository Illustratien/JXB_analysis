rm(list=ls())
pacman::p_load(dplyr,purrr,ggplot2,stickylabeller,ggh4x,
               toolPhD,toolStability)

options(dplyr.summarise.inform = FALSE)
# -------------------------------------------------------------------------

hsd_tbl <- read.csv("result/Table/sssp/HSDtable.csv") %>% 
  dplyr::select(groups,trait,treatment,Year,combi) %>% 
  rename(trait_name=trait,DFG_year=Year)
resdf <- read.csv("result/Table/sssp/HSDtable_phenology_gcd.csv") %>% 
  dplyr::select(groups,trait,treatment,Year,combi,BBCH) %>% 
  mutate(
    BBCH=gsub("\\\n","-",BBCH),
    trait=
      case_when(grepl("\\-",BBCH)~gsub("Acc_Temperature","S",trait),
                
                T~gsub("Acc_Temperature","TT",trait)),
    
    BBCH=gsub("\\-","_",BBCH),
    BBCH=case_when(grepl("GCD61|^S$","trait")~"",T~BBCH),
    Year=paste0("DFG",Year),
    
  ) %>% 
  rename(DFG_year=Year) %>% 
  
  mutate(
    trait_name=case_when(!trait==BBCH&(!BBCH=="")~paste(trait,BBCH,sep="-"),
                         T~trait),
    trait_name=gsub("TT-","TT",trait_name) %>% 
      gsub("S-","S",.) )%>% 
  select(-trait,-BBCH)
hsd_tbl <- bind_rows(hsd_tbl,resdf) %>% distinct()

# -------------------------------------------------------------------------
periodic <- read.csv(file = "result/Climate/periodic_BBCH_thermaltime.csv") %>% 
  rename("BBCH"="B",trait=Climate)
bbc <- read.csv("result/Climate/cumulative_BBCH_thermaltime.csv")%>% 
  dplyr::select(Year,var,nitrogen,appl,timeid,BBCH,any_of(starts_with("Acc"))) %>% 
  tidyr::pivot_longer(starts_with("Acc"),names_to="trait",values_to ="value") %>% 
  mutate(BBCH=as.character(BBCH))
GCD <-read.csv("data/GCD_merge.csv") %>% 
  dplyr::select(DFG_year:S,any_of(ends_with("61")),GCD87) %>% 
  tidyr::pivot_longer(GLA50:GCD87,names_to = "trait",values_to = "Trait") %>% 
  mutate(DFG_year=gsub("DFG","",DFG_year) %>% as.numeric(),
         BBCH="")

data_list <- bind_rows(periodic,bbc) %>%
  rename(DFG_year=Year,Trait=value) %>%
  bind_rows(.,GCD) %>%
  tidyr::unite("namCombine",c(trait,BBCH),sep="-") %>%
  mutate(DFG_year=paste0("DFG",DFG_year),
         namCombine=
           case_when(
             grepl("\\\n",namCombine)~
               gsub("Acc_Temperature-","S",namCombine) %>% 
               gsub("\\\n","_",.),
             T~gsub("Acc_Temperature-","TT",namCombine)
           ))

merge_dat <- read.csv("data/merge_profile.csv") %>% 
  select(timeid:appl,Trait,namCombine) %>% 
  filter(namCombine%in%c("TT61","TT87","tiller_61","spike_number_87")) %>% 
  tidyr::pivot_wider(names_from=namCombine,values_from=Trait) %>% 
  mutate(T61_T87=TT61/TT87,
         SN_T61=spike_number_87/tiller_61) %>% 
  select(-c(tiller_61:TT87)) %>%
  tidyr::pivot_longer(names_to="namCombine",values_to="Trait",cols=c(T61_T87,SN_T61)) %>% 
  bind_rows(.,read.csv("data/merge_profile.csv") %>% 
              select(timeid:appl,Trait,namCombine) 
  ) %>%
  bind_rows(.,data_list) %>% 
  filter(!DFG_year=="DFG2019") %>% 
  mutate(trait_name=namCombine,
         g=paste(nitrogen,appl,timeid,sep="_"),
         Year=gsub("DFG","",DFG_year) %>% factor()) %>% distinct()
# %>% 
# select(-namCombine)


# -------------------------------------------------------------------------

y <- read.csv("data/merge_profile.csv") %>%
  dplyr::filter(trait=='yield',part=='plot',!DFG_year=="DFG2019")

g_rank <- genotypic_superiority_measure(y,'Trait',genotype = 'var',
                                        environment = c('appl','DFG_year',
                                                        'timeid','nitrogen'),
                                        unit.correct = T) %>% 
  dplyr::select(-Mean.Trait)  %>%  
  arrange(genotypic.superiority.measure) %>% 
  rename('pg_rank'='genotypic.superiority.measure','var'='Genotype')

col_pal <- scales::hue_pal()(3)
names(col_pal) <- as.character(2019:2021)
source("src/fun/M_fun.R")

# uni <- read.csv("data/unit.csv")

# -------------------------------------------------------------------------
plot_trait<- function(tv,ltb){
  typ <- "var"
  rank_level<- g_rank[[typ]]
  
  label_table <-  merge_dat %>% 
    dplyr::filter(namCombine%in%tv) %>% 
    # out put ready format
    table_fun1(.,typ)%>% 
    mutate({{typ}}:=factor(.data[[typ]],rank_level))%>% 
    # remove NA rows
    dplyr::select(-groups,-combi) %>%
    tidyr::pivot_longer(-c(var,DFG_year),
                        names_to = "namCombine",
                        values_to = "groups") %>%
    na.omit() %>%
    left_join(ltb) %>% 
    arrange(letters) %>% 
    mutate( cvnam=paste0('~"("*',letters,'*")"~',nam)%>% factor(),
            unit=case_when(unit==""~unit,
                           T~paste0('(',unit,')')))
  
  p <- merge_dat%>% 
    left_join(label_table,.) %>%
    mutate(
      {{typ}}:=factor(.data[[typ]],rank_level),
      Year=factor(Year,levels=c("2020","2021")))%>% 
    arrange(var,Year) %>% 
    ggplot(.,aes(x=groups,y=Trait,color=Year,groups=Year))+
    geom_boxplot(outlier.shape=NA)+
    scale_x_discrete(labels = function(l) parse(text=l))+
    geom_point(position=position_jitterdodge(.4),
               alpha=.6,size=2,shape=1)+
    theme_classic()+
    ggh4x::facet_nested(var~cvnam+unit,independent = "y",
                        nest_line=T, 
                        scales="free",switch = "y",
                        labeller = label_parsed)+
    scale_color_manual(values = col_pal)+
    theme( 
      axis.text.x = element_text(size=10,angle=90),
      # ggtext element_markdown(), ... failed
      axis.text.y = element_text(size=10,# colour =color.vec
      ),
      axis.title.x=element_blank(),
      # nested_line
      strip.text.y =  element_text(size=8),
      strip.text.x =  element_text(size=8),
      strip.background = element_blank(),
      legend.position = "bottom",
      strip.placement = "outside",
      axis.title.y=element_blank())+
    coord_flip()
}

# individual input -------------------------------------------------------------------------
# tv <- c("plot_yield","plot_grain_number",
#         "plot_tkw","plot_grain_per_spike",
#         "spike_number_87","straw_DM87")
# 
# ltb <- data.frame(namCombine=tv,
#                   nam=c('GY','GN','TGW','GpS','SN','Straw'),
#                   letters =LETTERS[1:6],
#                   unit= c("t~ha^-1",'x10^6~ha^-1',
#                           'g/1000~grain','x10^6~ha^-1~spike^-1',
#                           'x10^6~ha^-1',"t~ha^-1"))
toolPhD::df_ue(merge_dat,"trait_name","g")

# quality -------------------------------------------------------------------------


# toolPhD::df_ue(merge_dat,"namCombine","LAI")
tv <- c( 
  # "straw_CHO61","straw_CHO87","mobile61_WSC",
  # "S61_87",
  # "SN_T61","T61_T87",
  "S", "GCD61",
"S31_41",
 "S61_71","S71_77", "S61_87",
  "TT61"
  # ,"TT87"
)
# tv <- c(
#   # "grain_area",
#   # "grain_bulk_density",
#         # "kernel_Nitrogen87",
#         # "kernel_total_N87",
#   # "straw_CHO87" ,"mobile61_WSC",
#   # "TT61",  "TT87",
#   "Height_61","Height_83","S31_41","S31_61",
#       # "S","GCD61",
#   "height_61"
#   # "straw_CHO61","straw_DM61","ear_DM61"  
#         )

ltb <- data.frame(namCombine=tv,
                  # nam=tv %>% gsub("\\_","~",.),
                  nam=c('GCD.S','GCD',
                        'TT[31-41]',
                        'TT[61-71]',
                         'TT[71-77]','TT[61-87]',
                        'TT61'),
                  unit=c("",'d*degree*C','d*degree*C', 'd*degree*C', 
                                'd*degree*C', 'd*degree*C', 'd*degree*C'#kg/t
                                ),
                  # unit=NA,
                  letters =LETTERS[1:length(tv)])

tiff(filename=paste0("paper_fig/FigS6.tiff"),
     units="cm",
     width=21,
     height=14,
     compression = "lzw",
     pointsize=3,
     res=600)# dpi
plot_trait(tv,ltb) %>% print()
dev.off()



# 61 trait ----------------------------------------------------------------


toolPhD::df_ue(merge_dat,"namCombine","LAI")

tv <- c("GCD61", "LAImax61",
        "straw_CHO61",
        "tiller_61","straw_CHO87" ,"Height_61")

ltb <- data.frame(namCombine=tv,
                  nam=tv,unit=NA,
                  letters =LETTERS[1:length(tv)])
plot_trait(tv,ltb) %>% print()

# -------------------------------------------------------------------------
# tv <- c("mobile61_WSC",
#         "WSCTE","straw_CHO87",
#         "straw_CHO61","ear_CHO61")
# 
# ltb <- data.frame(namCombine=tv,
#                   nam=c('Delta~WSC',
#                         'WSCTE',#mobile/all61
#                         'CHO["Straw,87"]','CHO["Straw,61"]',
#                         'CHO["Ear,61"]'),
#                   letters =LETTERS[1:5],
#                   unit= c("kg/ha",
#                           '10^-3',#kg/t
#                           'mg/g',
#                           'mg/g',"mg/g"))

# -------------------------------------------------------------------------
tv <- c("rachis_DM87","rachis_ratio_DM87",
        "grain_area","grain_bulk_density",
        "grain_width" ,"grain_length")

ltb <- data.frame(namCombine=tv,
                  nam=c('rachis','rachis/spike',
                        'area',#mobile/all61
                        'density','width',
                        'length'),
                  letters =LETTERS[1:6],
                  unit= c('t/ha','"%"',
                          'mm^2',#kg/t
                          'kg/hl',
                          'mm',"mm"))

# -------------------------------------------------------------------------

tv <- c("rachis_DM87","rachis_ratio_DM87",
        "grain_area","grain_bulk_density",
        "grain_width" ,"grain_length")

ltb <- data.frame(namCombine=tv,
                  nam=c('rachis','rachis/spike',
                        'area',#mobile/all61
                        'density','width',
                        'length'),
                  letters =LETTERS[1:6],
                  unit= c('t/ha','"%"',
                          'mm^2',#kg/t
                          'kg/hl',
                          'mm',"mm"))

# -------------------------------------------------------------------------
tv <- c("GCD61","S71_77",
        "S61_71" ,"S61_87",
        "gSPAD_61","gLI_int")

ltb <- data.frame(namCombine=tv,
                  nam=c('GCD','grain~filling',
                        'flowering~time',#mobile/all61
                        'post-anthesis~time','SPAD~61',
                        'LI~integral'),
                  letters =LETTERS[1:6],
                  unit= c('d*degree*C','d*degree*C',
                          'd*degree*C',#kg/t
                          'd*degree*C',
                          'a','a'))



# -------------------------------------------------------------------------

tv <- c("plot_yield","plot_grain_number",
        "plot_tkw","plot_grain_per_spike",
        "spike_number_87","grain_width" ,"grain_length")

ltb <- data.frame(namCombine=tv,
                  nam=c('GY','GN','TGW','GpS','SN','GW',"GL"),
                  letters =LETTERS[1:7],
                  unit= c("t~ha^-1",'x10^6~ha^-1',
                          'g/1000~grain','x10^6~ha^-1~spike^-1',
                          'x10^6~ha^-1','mm',"mm"))
# -------------------------------------------------------------------------
# get means for the non-significant traits 

p <- plot_trait(tv,ltb)
p
tiff(filename=paste0("paper_fig/FigS5.tiff"),
     units="cm",
     width=21,
     height=14,
     compression = "lzw",
     pointsize=3,
     res=600)# dpi
p %>% print()
dev.off()

# try to control the label text 
# color.vec <- merge_dat%>% 
#   dplyr::filter(namCombine%in%ltb$namCombine) %>%
#   mutate({{typ}}:=factor(.data[[typ]],rank_level))%>% 
#   left_join(label_table) %>%
#   group_by(letters,var,Year) %>%
#   summarise(m=mean(Trait)) %>% 
#   ungroup() %>% 
#   arrange(letters,var,desc(m)) %>% 
#   .$Year %>% as.character() %>% 
#   col_pal[.]
# 
# Table2<- merge_dat %>% 
#   # subset by condition 
#   dplyr::filter(!!condi) %>% 
#   # out put ready format
#   table_fun2(.,typ)%>% 
#   mutate({{typ}}:=factor(.data[[typ]],rank_level)) %>% 
#   select(-groups,-combi) %>%
#   tidyr::pivot_longer(-c(var,DFG_year),
#                       names_to = "namCombine",
#                       values_to = "groups") %>%
#   na.omit() %>%
#   tidyr::pivot_wider(names_from = "namCombine",values_from ="groups" )
