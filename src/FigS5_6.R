rm(list=ls())
pacman::p_load(dplyr,purrr,ggplot2,stickylabeller,ggh4x,
               toolPhD,toolStability)

options(dplyr.summarise.inform = FALSE)
# -------------------------------------------------------------------------

hsd_tbl <- read.csv("result/HSDtable.csv") %>% 
  dplyr::select(groups,trait,treatment,Year,combi) %>% 
  rename(trait_name=trait,DFG_year=Year)
resdf <- read.csv("result/HSDtable_phenology_gcd.csv") %>% 
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
periodic <- read.csv(file = "data/periodic_BBCH_thermaltime.csv") %>% 
  rename("BBCH"="B",trait=Climate)
bbc <- read.csv("data/cumulative_BBCH_thermaltime.csv")%>% 
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

# ------------------------------------------------------------------------
tv <- c("plot_yield","plot_grain_number",
        "plot_tkw","plot_grain_per_spike",
        "spike_number_87","grain_width" ,"grain_length")

ltb <- data.frame(namCombine=tv,
                  nam=c('GY','GN','TGW','GpS','SN','GW',"GL"),
                  letters =LETTERS[1:7],
                  unit= c("t~ha^-1",'x10^6~ha^-1',
                          'g/1000~grain','x10^6~ha^-1~spike^-1',
                          'x10^6~ha^-1','mm',"mm"))

p <- plot_trait(tv,ltb)

tiff(filename=paste0("result/Fig.S5.tiff"),
     units="cm",
     width=21,
     height=14,
     compression = "lzw",
     pointsize=3,
     res=600)# dpi
p %>% print()
dev.off()
# quality -------------------------------------------------------------------------

tv <- c( 
  "S", "GCD61",
  "S31_41",
  "S61_71","S71_77", "S61_87",
  "TT61"
)


ltb <- data.frame(namCombine=tv,
                  nam=c('GCD.S','GCD',
                        'TT[31-41]',
                        'TT[61-71]',
                        'TT[71-77]','TT[61-87]',
                        'TT61'),
                  unit=c("",'degree*C*d','degree*C*d', 'degree*C*d', 
                         'degree*C*d', 'degree*C*d', 'degree*C*d'#kg/t
                  ),
  
                  letters =LETTERS[1:length(tv)])

tiff(filename=paste0("result/Fig.S6.tiff"),
     units="cm",
     width=21,
     height=14,
     compression = "lzw",
     pointsize=3,
     res=600)# dpi
plot_trait(tv,ltb) %>% print()
dev.off()
