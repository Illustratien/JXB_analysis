rm(list = ls())
pacman::p_load(dplyr,purrr,ggplot2,magrittr);options(dplyr.summarise.inform = FALSE)

management <-  read.csv("data/management_time.csv") %>% 
  mutate(date=as.Date(date,format = "%d.%m.%Y",tz = "CET"),
         DFG_year=sub('DFG','',DFG_year)) %>% 
  rename('Year'='DFG_year')


climate <- read.csv("data/climate.csv")%>%
  mutate(
    DayTime=as.Date(DayTime,format="%Y-%m-%d",tz = "CET"),
    sowing_date=stringr::str_to_title(sowing_date),
    DFG_year=sub('DFG','',DFG_year))%>%
  rename("timeid"="sowing_date",'Year'='DFG_year') %>% 
  mutate(Acc_Radiation=Acc_Radiation/1000)# W to kw  
# select three important stages for highlighting
# note that this is the harvest biomass cut time not the evaluating time
three_stages<- management %>% 
  dplyr::filter(stage%in%paste0("BBCH",c(51,61,87))) %>%
  dplyr::select(date,Year,stage) %>% arrange(date) %>% 
  left_join(.,climate[,c("DayTime","Year","timeid","Das")] %>%
              rename("date"="DayTime"))

#For displaying purpose, select only one level to show 
sowingExample <- 'Early'

# shape vector for different years
shape.vec <- c(1,2,0)
names(shape.vec) <- c(2019:2021)
# function-------------------------------------------------------------------------

sfun <- function(df){
  # small function for calculation of periodic subtraction 
  # input a dataframe
  # output a dataframe
  # column name
  coln <-paste0("B",df$BBCH[1]) 
  df %>% 
    # change the column name to BBCH stage
    dplyr::rename(!!coln:="value") %>%
    # remove disturbing columns for pivot_longer
    dplyr::select(-contains(c("BBCH","Date")))
}

# # Daily climate plot -----------------
# try to label key stages with either rectangle or points
reproductive_stage_df <- three_stages %>%
  dplyr::filter(timeid=="Early",stage%in%paste0("BBCH",c(51,61,87)))
# label_df <- reproductive_stage_df[rep(1:nrow(reproductive_stage_df),6),] %>% mutate(value=-5)

res.l <- tidyr::gather(climate,"climate","value", DailyMean_Temperature:DailySum_Radiation) %>% 
  mutate(climate=factor(climate,levels = c('DailyMean_Temperature','DailySum_Radiation',
                                           'DailySum_Percipitation','Acc_Temperature',
                                           'Acc_Radiation','Acc_Percipitation')),
         # You can use recode() directly with factors; 
         # it will preserve the existing order of levels while changing the values.
         Climate=recode(climate,
                        'DailyMean_Temperature'= 'Daily~mean~temperature~(degree*C)',
                        "DailySum_Radiation"="Daily~irradiation~(W~m^-2)",
                        "DailySum_Percipitation"="Daily~precipitation~(mm)",
                        "Acc_Percipitation"="Accumulated~precipitation~(mm)",
                        "Acc_Radiation"="Accumulated~Irradiation~(kW~m^{-2})",
                        "Acc_Temperature"="Accumulated~thermal~sum~(degree*C*d)"))
#data
daily_dat <- res.l%>% filter(timeid==sowingExample,!grepl("Acc",climate)) %>% 
  dplyr::select(-c(DayTime:timeid,Climate)) %>% 
  tidyr::pivot_wider(.,names_from = 'climate',values_from = 'value')
#label df
label_day <- reproductive_stage_df%>% 
  mutate(stage=sub("BBCH","",stage)) %>% 
  left_join(.,daily_dat )
#look up table for title
lookup <- res.l %>% dplyr::select(Climate,climate) %>% distinct()
#legend position
lgd.pos <- list('none','none',c(.83,.10))
# prefix for subplot
pre.vec <- c('A.','B.','C.')
# List of plots
daily_plist <- purrr::map(1:3,~{
  y <- daily_dat %>% names() %>% .[length(.)-.x+1]
  yt <- paste0(pre.vec[.x],'~',
               with(lookup,lookup[match(y,climate),'Climate']) %>% as.character())
  # anno <- data.frame(Das=150,y=daily_dat[[y]] %>% max(),label=)
  
  p <- ggplot()+
    geom_line(data=daily_dat,
              alpha=.5,aes_string(x='Das',y=y,color='Year',linetype="Year"),
              size=.6,show.legend = F)+
    # facet_wrap(~Climate,scales="free",nrow=3,labeller = label_parsed)+
    theme_bw()+
    ggtitle(parse(text=yt))+
    scale_linetype_manual(values=c("solid","dotted","dashed"))+
    scale_y_continuous(
      labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
    theme(legend.position = "none",
          legend.direction = 'horizontal',
          plot.margin = margin(0.5,r=0.5,0,l= .5, "cm"),
          axis.title.y=element_blank(),
          legend.text = element_text(size=8),
          legend.title = element_text(size=8),
          strip.background = element_blank(),
          legend.background = element_blank(),
          plot.title=element_text(family='',hjust = .05, face='bold', size=14, margin=margin(l=2,t=5,b=-30)),
          strip.text = element_text(size=8))+
    # ylab( parse(text = yt))+
    geom_point(data=label_day ,
               aes_string('Das',y,color='Year',shape='stage'),size=1.5,alpha=.8)
  # p
  if(.x==3){
    p+xlab('Days after sowing')+
      theme(
        # legend.spacing.y = unit(-0.2, "cm"),
        legend.position = "bottom")+
      guides(color=guide_legend(
        override.aes = list(size=3),
        keywidth=0.2,
        keyheight=0.2,
        default.unit="cm",title = 'Year'),
        shape=guide_legend(
          override.aes = list(size=3),
          keywidth=0.2,
          keyheight=0.2,
          default.unit="cm",title = 'Biomass cut'),nrow=1)
  }else{
    p+theme(axis.title.x = element_blank(),
            axis.text.x= element_blank())
  }
})

tiff(filename='result/plot/FigS3.tiff',
     units="cm",
     width=17.4,
     height=15,
     compression = "lzw",
     pointsize=12,
     res=400)# dpi)
cowplot::plot_grid(plotlist =daily_plist,nrow = 3,
                   rel_heights = c(.7,.7,1) ) %>% print()
dev.off()
