#get_peak_halfwidths.R

tic()
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/def_peak_halfwidths_rise_decaytimes_v1.R")

halfWidths = getHalfWidths(df = fitsCleaned, groupers = groupers_plus)

print("Finished getting full-width half-maximums.")
toc()
