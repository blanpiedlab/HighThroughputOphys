my.theme<- theme(plot.title=element_text(size=40, family="sans"), 
                 plot.subtitle=element_text(size=30, family="sans"),
                 axis.line = element_line(colour = "black",size=3),
                 axis.text.x=element_text(colour="black", size=40, family="sans"),
                 axis.text.y=element_text(colour="black", size=40, family="sans"),
                 axis.title=element_text(colour="black", size=42, family="sans"), 
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.position="none",
                 legend.title = element_blank(),#element_text(colour="black", size=22, family="sans"), 
                 legend.text=element_text(colour="black", size=36, family="sans"),
                 strip.text= element_text(colour="black", size=20, family="sans"),
                 axis.ticks=element_line(colour="black",size=2),
                 axis.ticks.length=unit(0.75, "cm"),
                 panel.spacing = unit(4, "lines")                                      

)



