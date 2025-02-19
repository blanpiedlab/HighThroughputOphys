
scalefactor=0.75


my.theme<- theme(plot.title=element_text(size=scalefactor*38, family="sans",hjust=0.5), 
                 plot.subtitle=element_text(size=scalefactor*38, family="sans",hjust=0.5),
                 axis.line = element_line(colour = "black",linewidth=scalefactor*1.5),
                 axis.text.x=element_text(colour="black", size=scalefactor*34, family="sans"), #angle=45, hjust=1),
                 axis.text.y=element_text(colour="black", size=scalefactor*34, family="sans"),
                 axis.title=element_text(colour="black", size=scalefactor*38, family="sans"), 
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.position="right",
                 legend.title = element_blank(),#element_text(colour="black", size=22, family="sans"), 
                 legend.text=element_text(colour="black", size=scalefactor*34, family="sans"),
                 strip.text= element_text(colour="black", size=scalefactor*36, family="sans"),
                 axis.ticks=element_line(colour="black",linewidth=scalefactor*1.5),
                 axis.ticks.length=unit(0.25, "cm"),
                 panel.spacing = unit(0.5, "cm"),
                 plot.tag = element_text(colour="black", size=scalefactor*48, family="sans",face="bold")#,#vjust = 0.5, hjust = -0.5)                                       

)



