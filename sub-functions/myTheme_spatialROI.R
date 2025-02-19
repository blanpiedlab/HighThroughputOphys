

scalefactor=0.75


my.theme_spatialROI<- theme(plot.title=element_blank(),#element_text(size=36, family="sans"), # 
                 plot.subtitle=element_blank(),#element_text(size=36, family="sans"),
                 axis.line = element_blank(),#element_line(colour = "black",linewidth=3),
                 axis.text.x=element_blank(),#element_text(colour="black", size=28, family="sans"),
                 axis.text.y=element_blank(),#element_text(colour="black", size=28, family="sans"),
                 axis.title=element_blank(),#element_text(colour="black", size=36, family="sans"), 
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.position="right",
                 legend.title = element_text(colour="black", size=scalefactor*30, family="sans"), 
                 legend.text=element_text(colour="black", size=scalefactor*28, family="sans"),
                 legend.key.size = unit(scalefactor*1.5, "cm"),
                 legend.spacing.y = unit(scalefactor*3, 'cm'),
                 #legend.margin = 
                 legend.box.margin=margin(r=scalefactor*1,l=scalefactor*0.5,t=scalefactor*0.5,b=scalefactor*0.5, unit='cm'),
                 strip.text= element_blank(),#element_text(colour="black", size=26, family="sans"),
                 axis.ticks=element_blank(),#element_line(colour="black",linewidth=2),
                 #axis.ticks.length=unit(0.5, "cm"),
                 panel.spacing = unit(0.3, "cm"),
                 #plot.tag = element_blank() #element_text(colour="black", size=48, family="sans",face="bold")
                 plot.tag.position = c(0.015,0.915),
                 plot.tag = element_text(colour="white", family="sans", face="bold", size=scalefactor*48,vjust = 0, hjust = 0)                                 

)



