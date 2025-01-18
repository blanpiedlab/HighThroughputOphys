
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/lighten_up_df.R") 

#break_me""
			
tic()

light_df<-lighten_up(df = df,remove_cols=drops)

print("finished lightening things up.")

toc()

