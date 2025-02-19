#getTable.R
tic()
source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3 - PTPsigma/table.R") 



table = getTable( df = peakStats, grouper1 = grouper1, grouper2 = grouper2)

print("Finished printing a data table.")


toc()