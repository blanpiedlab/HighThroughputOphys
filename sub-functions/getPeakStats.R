#getPeakStats.R


tic()
source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3 - PTPsigma/peakStats.R") 


peakStats<-peakStats(df = fitsCleaned,dfHalfWidths = halfWidths, groupers = groupers,plotBy = plotBy,varList=varList, levels = levels, secondAxis=secondAxis, color_override=color_override)


print("Finished getting peak stats.")
toc()


