# schmittTrig_v4.R
# Schmitt Trigger Function call with plotting to iterate over a dataframe by group.
# Identifies every peak in the dataset on an ROI by ROI basis (assuming called within grouped dplyr pipe) and scores as a 0 or 1. 
# Plots what the Schmitt trigger is doing during the call. 

# Variable calls as of 4.24.22
# x=.$absoluteTime, y=.$dFF,
# a=.$normIntensity, b=.$normBaseline,
# c=.$baseline_idx, 
# d=.$normStDev, e=.$knot_color,
# f=.$normSplineBaseline,
# threshold=threshold, std=.$dFF_stdev,thresh_stdev=thresh_stdv

library(gsignal)


schmittTrig_v4<- function(dataframe, time, dFF, normIntensity, normBaseline, outlier_indices, threshold, std, ymax, ymin, thresh_stdev, smoothed_dFF) {
    lvl = mean(std, na.rm = TRUE)*threshold
    trig = schtrig(dFF, lvl) # function from gsignal which defines a schmitt and applies it to column y (dFF)
    signal = trig$v
    ymax = mean(ymax, na.rm=TRUE)
    ymin = mean(ymin, na.rm=TRUE)
    #colors = c("#DBDBDB","#1E90FF") # potentially useful for defining colors of signals in plot
    if (showGraphs == TRUE) {
    par(mfrow = c(5,1),oma = c(2,2,0,0) + 0.1,mar = c(3.7,3.7,1,1) + 0.2) #tiling 3 by 1 plots
    
    #normIntensity vs time with normBaseline overlay
    plot(time,normIntensity, main = 'Normalized intensity signal with outlier indices', xlab="",col = "white", type = "l")    
    lines(time,normIntensity, col = "black", type="l",lwd=1)
    points(time,normIntensity, col = outlier_indices, pch = 21, cex=1, lwd=1)


    #normIntensity vs time with normBaseline overlay
    plot(time,normIntensity, main = 'Normalized signal with rolling median baseline (window = 100)', xlab="",col = "black", type = "l")    
    lines(time,normBaseline, col = "red", lwd=2)
    #points(time,normBaseline, col = outlier_indices, pch = 21, cex=3, lwd=3)

    #lines(time,(normBaseline+(normStDev*thresh_stdev)), type='l',col = 'green', lwd=2)
    #lines(time,(normBaseline-(normStDev*thresh_stdev)), type='l',col = 'green', lwd=2)
    
    
    #Show dFF vs time with Schmitt Trigger levels
    plot(time,dFF, ylim = c(ymin,ymax), type="l",main = "dFF with peaks detected", xlab="",col="black",lwd=1)
    abline(h = lvl, col = "blue")

    #Show smoothed dFF vs time with Schmitt Trigger levels
    plot(time,smoothed_dFF, ylim = c(ymin,ymax), type="l",main = "dFF with peaks detected", xlab="",col="black",lwd=1)
    abline(h = lvl, col = "blue")

    
    #binary state plot to show schmittTrig behavior
    plot(time,signal,type="S",main="schmittTrig behavior",col="red",ylab="",xlab="Time (s)",ylim=c(-1.5,1.5),lwd=2)
    
    cbind(dataframe,signal)
    } else {

    cbind(dataframe,signal)
    

    }
}




    #legacy
    #compare old baseline, knots, and splineBaseline
    #plot(time,b, col = "black", type="l")
    #curve(time,f, col = "red", type = "l")
    #points(time,b, col = e, pch = 21, cex=3, lwd=3)
    #lines(time,b, col = "red")    
    #knots go here too