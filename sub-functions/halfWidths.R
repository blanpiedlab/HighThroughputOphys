#HalfWidths.R

tic()
source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/getHalfWidths.R") 

halfWidths = getHalfWidths(df = fitsCleaned, groupers = groupers_plus)

print("Finished getting full-width half-maximums.")
toc()
