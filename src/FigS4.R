rm(list=ls())
pacman::p_load(dplyr,magrittr,colorRamps,ggplot2)
df <- readRDS("data/diviner.rds")
# Diviner and percipitation-------------------------------------------------------------------------
diviner_names <- c( 'm30'="A 0-30 cm",
                    'm60'="B 30-60 cm",
                    'm100'="C 60-100 cm")

p <- ggplot() +
  geom_line( df$diviner_sub,
             mapping=aes(x=Das,y=humidity,
                         color=Year),
             linetype='twodash',
             size=1.2) +
  geom_line( df$diviner_climate, 
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
  guides(color=guide_legend(
    override.aes = list(size=1),
    keywidth=0.2,
    keyheight=0.3,
    default.unit="cm",title = 'Year'))
# p
tiff(filename='result/plot/FigS4.tiff',
     units="cm",
     width=17.4,
     height=8,
     compression = "lzw",
     pointsize=12,
     res=400)# dpi)
p %>% print()
dev.off()


