
#***************** Up in the Air: The Effects of Air Pollution on Conflict in West Africa**********************

#        MASTER THESIS at the Georg-August University of Goettingen 

#                       Ankhbayar Delgerchuluun

#                               RSCRIPT
#                          ------------------
#                          Data Visualization

####____________________________README_____________________________________________####

#This R script carries out all the figures that included in the thesis.



####_________Required packages_____________________________________________________#####
#packages                    
library(dplyr, warn.conflicts = FALSE)
library("raster")
library(lubridate)
library(ggplot2)
library(hrbrthemes)
library("raster")
library("terra")
library("ggthemes")
library(viridis)
library("patchwork")
library("rgdal")
library("ggpubr")
library("maps")
library("mapproj")

####_________Colors_____________________________________________________#####
###colors
color<-"#e7e2dc"
color1<-"#a6611a"
color2<-"#dfc27d"
color3<-"#80cdc1"
color4<-"#018571"
color5<-"#555577"
color6<-"#4a9143"
color7<-"#020E24"
color8<-"#234178"
color9<-"#CFDDF5"
color10<-"#779147"
color11<-"#ea573d"
color12<-"#2a251f"
color13<-"#e6c69d"
color14<-"#e59d37"
color15<-"#72a278"

West<-vect('West.shp')
wa<-readOGR("West.shp")
season<-c("Harmattan", "Non-Harmattan")
dates_month <- seq(as.Date("2003-01-01"), by = "month", length.out = 76)
dates_month<-as.Date(dates_month)
month_single<-format(as.Date(dates_month), "%m")
month_single<-as.numeric(month_single)
harmattan<-ifelse(month_single<4 | month_single>10, 1,0)

data<-read.csv("datasets/data_01.csv")
summary(data)
#### FIGURE 1: Air pollution over time  ####
figure1<-data %>%
  group_by(date)%>% 
  summarise(pm25=mean(pm25_mg,na.rm=TRUE), 
            midas=mean(midas,na.rm=TRUE))
summary(figure1)

coef<-0.1
figure1<-mutate(figure1, date=as.Date(date))
ggplot(figure1, aes(x=date))+
  geom_line(aes(y=pm25, colour="PM2.5 mg/m3 (LHS)"),  size=1.25)+
  geom_line(aes(y=midas*coef,  colour="DoD Midas  at 550nm (RHS)"), size=1.25)+
  scale_colour_manual(name="", values=c("PM2.5 mg/m3 (LHS)"=color11,"DoD Midas  at 550nm (RHS)"=color12))+
  scale_y_continuous(
    name="Air pollution PM2.5",
    sec.axis=sec_axis(~./coef, name="Dust measurment"))+
  theme_ipsum()+
  xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.grid.minor.y = element_line(), axis.line = element_line(colour = "black"))+
  theme(legend.position="bottom",legend.text=element_text(size=14))


#### FIGURE 2: Spatial distribution of Air pollution PM2.5 ####

pm25 <- rast('datasets/final_temp/pm25_rast_2.tif')
pm25_g<-tapp(pm25, harmattan, fun=mean, rm.na=TRUE)

names(pm25_g)<-season
pm25_df<-as.data.frame(pm25_g, xy=TRUE)
names(pm25_df)[4]<-"Nonharmattan"
summary(pm25_df)
pm<-c(0, 100, 200, 300)
harm_pm25<-ggplot()+
  geom_tile(pm25_df, mapping= aes(x=x, y=y, fill=Harmattan), alpha=0.8)+
  geom_polygon(data = wa, aes(x = long, y = lat, group = group), 
               colour = "black", fill = NA)+
  ggtitle("Harmattan period")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank(),legend.position="bottom", legend.box = "horizontal", plot.title = element_text(hjust = 0.5), legend.title = element_blank())+
  scale_fill_gradient(low=color6, high=color11, limits=c(0, 300))+
  #scale_fill_viridis(scale_color_gradientn, alpha=0.9)+
  coord_equal()

nonh_pm25<-ggplot()+
  geom_tile(pm25_df, mapping= aes(x=x, y=y, fill=Nonharmattan), alpha=0.8)+
  geom_polygon(data = wa, aes(x = long, y = lat, group = group), 
               colour = "black", fill = NA)+
  ggtitle("Non-Harmattan period")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank(),legend.position="bottom", legend.box = "horizontal", plot.title = element_text(hjust = 0.5), legend.title = element_blank())+
  scale_fill_gradient(low=color6, high=color11, limits=c(0, 300))+
  #scale_fill_viridis(scale_color_gradientn, alpha=0.9)+
  coord_equal()

ggarrange(harm_pm25, nonh_pm25,legend=c("bottom"),common.legend = TRUE)




#### FIGURE 3: Conflict over time ####
conflict_g<-read.csv("datasets/conflict/ACLED/conflict_g.csv")
summary(conflict_g)
barc1<-ggplot(conflict_g, aes(x=year,  y=death,  fill=type)) + 
  geom_bar(position="stack", stat="identity", width=0.8)+
  ggtitle("Number of fatalities")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_blank(),legend.position="bottom", legend.box = "horizontal", 
        plot.title = element_text(hjust = 0.5), panel.grid.minor.y = element_line(), legend.title = element_blank())+
  scale_fill_manual(values=c(color10, color11, color12, color13, color14, color15))+
  xlab("")

barc2<-ggplot(conflict_g, aes(x=year, y=conflict, fill=type)) + 
  geom_bar(position="stack", stat="identity", width=0.8)+
  ggtitle("Number of conflict")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_blank(),legend.position="bottom", legend.box = "horizontal", 
        plot.title = element_text(hjust = 0.5), panel.grid.minor.y = element_line(), legend.title = element_blank())+
  scale_fill_manual(values=c(color10, color11, color12, color13, color14, color15))+
  xlab("")

ggarrange( barc2,barc1, legend=c("bottom"),common.legend = TRUE)
#scale_fill_viridis(scale_color_gradientn, alpha=0.9)+


#### FIGURE 4: Spatial distribution of conflict incidences ####

figure3<-data %>%
  group_by(x, y)%>% 
  summarise(cfatal=sum(cdeath,na.rm=TRUE), conflict=mean(conflict,na.rm=TRUE), 
            ccount=sum(ccount,na.rm=TRUE))       
figure3<-fortify(figure3)


co1<-ggplot() +
  geom_tile(data = figure3, aes(fill = ccount, x = x, y = y)) +
  binned_scale(aesthetics = "fill",
               scale_name = "stepsn", 
               palette = function(x) c("white",color4, 
                          color6, color11),
               breaks = c(0.5, 1, 5, 25),
               #limits = c(0, 100),
               #guide = guide_coloursteps(even.steps = TRUE,
                #                       show.limits = TRUE),
               position=c("bottom"),
               name=NULL,
               guide="legend",
               labels=c("","1-5","5-25","25+")
  )+
  ggtitle("Number of conflict")+
  geom_polygon(data = wa, aes(x = long, y = lat, group = group), 
               colour = "black", fill = NA, size=0.8)+
  theme_void() +
  theme(plot.title = element_text(hjust=0.5))+
  coord_equal()

co2<-ggplot() +
  geom_tile(data = figure3, aes(fill = cfatal, x = x, y = y)) +
  binned_scale(aesthetics = "fill",
               scale_name = "stepsn", 
               palette = function(x) c("white",color4, 
                                       color6, color11),
               breaks = c(0.5,10, 99, 101),
               limits = c(0, 1000),
               #guide = guide_coloursteps(even.steps = TRUE,
               #                       show.limits = TRUE),
               position=c("bottom"),
               name=NULL,
               guide="legend",
               labels=c("","1-10","10-99","100+")
  )+
  ggtitle("Number of fatalities")+
  geom_polygon(data = wa, aes(x = long, y = lat, group = group), 
               colour = "black", fill = NA, size=0.8)+
  theme_void() +
  theme(plot.title = element_text(hjust=0.5))+
  coord_equal()

ggarrange( co1,co2, legend=c("bottom"))


#### FIGURE 5: Spatial distribution of DoD Midas #####

midas <- rast('datasets/final_temp/midas.tif')
midas<-crop(midas, West, snap="in", mask=T)
midas_g<-tapp(midas, harmattan, fun=mean, na.rm=TRUE)
names(midas_g)<-season
midas_df<-as.data.frame(midas_g, xy=TRUE)
names(midas_df)[4]<-"nonh"
summary(midas_df)
summary(midas_df)

harm_midas<-ggplot()+
  geom_tile(midas_df, mapping= aes(x=x, y=y, fill=Harmattan), alpha=0.8)+
  geom_polygon(data = wa, aes(x = long, y = lat, group = group), 
               colour = "black", fill = NA)+
  ggtitle("Harmattan period")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank(),legend.position="bottom", legend.box = "horizontal", plot.title = element_text(hjust = 0.5), legend.title = element_blank())+
  scale_fill_gradient(low=color6, high=color11, limits=c(0, 1))+
  #scale_fill_viridis(scale_color_gradientn, alpha=0.9)+
  coord_equal()

nonh_midas<-ggplot()+
  geom_tile(midas_df, mapping= aes(x=x, y=y, fill=nonh), alpha=0.8)+
  geom_polygon(data = wa, aes(x = long, y = lat, group = group), 
               colour = "black", fill = NA)+
  ggtitle("Non-Harmattan period")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank(),legend.position="bottom", legend.box = "horizontal", plot.title = element_text(hjust = 0.5), legend.title = element_blank())+
  scale_fill_gradient(low=color6, high=color11, limits=c(0, 1))+
  #scale_fill_viridis(scale_color_gradientn, alpha=0.9)+
  coord_equal()

ggarrange(harm_midas, nonh_midas)
