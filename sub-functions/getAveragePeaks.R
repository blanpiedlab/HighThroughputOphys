#getAveragePeaks.R

tic()
source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3 - PTPsigma/averagePeaks.R") 

averagePeaks(df = fitsCleaned, groupers = groupers_plus, plotBy = plotBy, levels = levels, secondAxis = secondAxis,color_override = color_override)


print("Finished plotting averaged peaks.")
toc()

