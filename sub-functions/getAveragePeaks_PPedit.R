#getAveragePeaks_PPedit.R

tic()
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v4 - markerSegmentation/averagePeaks_PPedit.R") 

averagePeaks(df = fitsCleaned, groupers = groupers_plus, plotBy = plotBy, levels = levels, secondAxis=secondAxis, color_override = color_override,  thirdAxis = thirdAxis,thirdAxis_levels = thirdAxis_levels)
				


print("Finished plotting averaged peaks.")
toc()

