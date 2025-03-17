rm(list=ls())
pacman::p_load(dplyr,magrittr,colorRamps,ggplot2)
interp_df <- function(df){
  # generate equal spaced xy interpolated df
  # using interpolation function from akima
  # https://www.rdocumentation.org/packages/akima/versions/0.6-2.1/topics/interp
  interp.list <- with(df,akima::interp(Date%>%
                                         strftime(.,format = "%j")%>%as.numeric(),
                                       depth, moisture))
  # length of the vector
  vec.l <- interp.list[[1]]%>%length()
  # rotate the resulted matrix
  interp.list[[3]] <-  t(apply(interp.list[[3]], 2, rev))
  
  interp.df <- purrr::map_dfr(1:vec.l,~{
    data.frame(x= as.Date(interp.list[[1]]%>%rev%>%.[.x],
                          origin = paste0(df$DFG_year[1] %>%substr(.,4,8),"-01-01"))-1,
               y=interp.list[[2]],
               z=interp.list[[3]][,.x])
    
  })
  return( interp.df)
}

#For displaying purpose, select only one level to show 
sowingExample <- 'Early'
management <-  read.csv("data/management_time.csv") %>% 
  mutate(date=as.Date(date,format = "%d.%m.%Y",tz = "CET"),
         DFG_year=sub('DFG','',DFG_year)) %>% 
  rename('Year'='DFG_year')
exp  <- map_dfr(2019:2021,~{
  read.csv(paste0("data/Exp_DFG",.x,"_design.csv")) %>% 
    mutate(DFG_year=paste0("DFG",.x))
})
# filter the sowing date dataframe
sow <-  management %>% dplyr::filter(activity=="sowing")%>%
  dplyr::select(date,Year,treatment)
# Diviner -----------------------------------------------------------------
dir_diviner <-'raw_data/Diviner'
fname <- list.files(dir_diviner)
# read diviner
diviner<- purrr::map_dfr(fname,~{read.csv(paste0(dir_diviner,"/",.x))}) %>% 
  group_by(Date)%>%
  mutate(
    Date=as.Date(Date,format="%Y-%m-%d",tz = "CET"),
    # Date=case_when(grepl("\\.",Date)~as.Date(Date,format="%d.%m.%Y"),
    #                T~as.Date(Date,format="%m/%d/%Y")
    depth=sub('X','',depth)%>%
      as.numeric()) %>%
  # merge sowing date to the diviner 
  # use later for filter the smallest range
  # there is 2018 data in diviner 11 in 2019 from Caroline
  merge(.,sow %>% dplyr::filter(treatment=="Early") %>% 
          mutate(DFG_year=paste0("DFG",Year)) %>% 
          dplyr::select(DFG_year,date))

diviner_range <- diviner %>% group_by(DFG_year) %>% 
  summarise(min=min(Date[Date>date]),max=max(Date))

climate <- read.csv(paste0("data/climate.csv"))%>%
  mutate(
    DayTime=as.Date(DayTime,format="%Y-%m-%d",tz = "CET"),
    sowing_date=stringr::str_to_title(sowing_date),
    DFG_year=sub('DFG','',DFG_year))%>%
  rename("timeid"="sowing_date",'Year'='DFG_year') 

# select early as example 
climate_sub<- climate %>% 
  dplyr::filter(timeid=="Early") %>% 
  dplyr::select(grep("(Percipitation|DayTime|Das)",
                     climate %>% names)) %>%
  rename("Date"="DayTime") 

#merge diviner with climate 
diviner_climate <- purrr::map_dfr(1:nrow(diviner_range),~{
  climate_sub %>% 
    # adjust the lower limit of precipitation to ten days before diviner
    dplyr::filter(between(Date,diviner_range[.x,]$min-10,
                          #doesn't change the maximum
                          diviner_range[.x,]$max)) %>% 
    mutate(Year=diviner_range[.x,"DFG_year"] %>% sub("DFG","",.))
})


diviner_sub <- diviner %>%
  group_by(DFG_year,Date) %>% 
  # summarise for three layers
  summarise(m30=mean(value[depth<=30]),
            m60=mean(value[depth<=60&depth>30]),
            m100=mean(value[depth>60])) %>%
  tidyr::pivot_longer(starts_with("m"),
                      values_to = "humidity",names_to = "depth") %>% 
  mutate(depth=factor(depth,levels=c("m30","m60","m100")),
         DFG_year=sub('DFG','',DFG_year)) %>% 
  left_join(.,diviner_climate[,c("Date","Das")]) %>% 
  rename('Year'='DFG_year') %>% 
  na.omit()#remove the na from 2018 diviner
df_lst<- diviner %>% 
  group_by(DFG_year)%>% 
  filter(Date>date) %>% # remove the 2018 winter diviner
  group_by(Date,depth,DFG_year)%>%
  summarise(moisture=mean(value))%>%
  mutate(Date=Date,
         depth=sub('X','',depth)%>%
           as.numeric()) %>% 
  ungroup() %>% 
  group_by(DFG_year) %>% group_split()
df_interpolate <- map(df_lst,~{
  interp_df(.x) %>% 
    mutate(Year=.x$DFG_year %>% sub("DFG","",.) %>% .[1]) %>% 
    left_join(.,sow %>% dplyr::filter(treatment==sowingExample)) %>% 
    mutate(Das=x-date)
})

# interpolate df,not necessarily better??
dfI<- df_interpolate %>% map_dfr(.,~{.x}) %>% 
  rename("Date"="x","depth"="y","value"="z") %>% 
  arrange(Date) %>% 
  dplyr::select(-date) %>% 
  distinct() %>% 
  group_by(Year,Date) %>% 
  # summarise for three layers
  summarise(m30=mean(value[depth<=30]),
            m60=mean(value[depth<=60&depth>30]),
            m100=mean(value[depth>60]),
            Das=Das %>% as.numeric()) %>% 
  tidyr::pivot_longer(starts_with("m"),
                      values_to = "humidity",names_to = "depth") %>% 
  mutate(depth=factor(depth,levels=c("m30","m60","m100")))

# select the BBCH stage to focus on
peri_vec <- c(9,13,23,31,41,51,61,71,77,87) %>% stringr::str_pad(string = .,width = 2,pad = "0")#c(23,41,51,61,71,77,87)


diviner_bbc <- map_dfr(2019:2021,~{
  read.csv(paste0("raw_data/Field_Measurement/BBCH_",.x,".csv"))
}) %>%
  tidyr::pivot_longer(starts_with("BBCH"),names_to = "BBCH",values_to="Date") %>% 
  mutate(Date=case_when(grepl("\\.",Date)~as.Date(Date,format="%d.%m.%Y"),
                        T~as.Date(Date,format="%m/%d/%Y")),
         BBCH=sub("BBCH","",BBCH)) %>% 
  merge(.,exp ,c("plot_id","DFG_year")) %>% 
  dplyr::filter(BBCH%in%peri_vec,timeid==sowingExample) %>% 
  dplyr::select(DFG_year:Date) %>% distinct() %>% 
  mutate(DFG_year=sub('DFG','',DFG_year)) %>% 
  rename("Year"="DFG_year") %>% 
  left_join(.,diviner_sub) %>% na.omit()
# Diviner and percipitation-------------------------------------------------------------------------
diviner_names <- c( 'm30'="(A) 0-30 cm",
                    'm60'="(B) 30-60 cm",
                    'm100'="(C) 60-100 cm")

p <- ggplot() +
  geom_line( diviner_sub,
             mapping=aes(x=Das,y=humidity,
                         color=Year),
             linetype='twodash',
             size=1.2) +
  geom_line( diviner_climate, 
             mapping=aes(x=Das,color=Year,
                         y=DailySum_Percipitation),
             size=.8,
             alpha=.3,show.legend = F) + 
  
  facet_grid(~depth,labeller=as_labeller(diviner_names))+
  # Divide by 10 to get the same range as the temperature
  scale_y_continuous(
    # Features of the first axis (left)
    name = bquote(atop("Volumetric soil water content ("*theta[v]*") %",
                       "\ndashed line")),
    # Add a second axis and specify its features (right)
    sec.axis = sec_axis(trans=~.*4/10,
                        name="Precipitation (mm), solid line"))+
  theme_bw()+
  theme(legend.position = c(.27,.86),
        legend.background = element_blank(),
        strip.background = element_blank(),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6),
        axis.title = element_text(size=9),
        plot.title = element_text(size=16),
        strip.text = element_text(size=10))+
  scale_x_continuous(labels = seq(150,270,30),breaks = seq(150,270,30))+
  xlab("Days after sowing")+
  # ggtitle('Precipitation v.s. soil moisture')+
  guides(color=guide_legend(
    override.aes = list(size=1),
    keywidth=0.2,
    keyheight=0.3,
    default.unit="cm",title = 'Year'))
# p
tiff(filename='paper_fig/FigS4.tiff',
     units="cm",
     width=17.4,
     height=8,
     compression = "lzw",
     pointsize=12,
     res=400)# dpi)
p
dev.off()


