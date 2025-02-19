
# get_mutate_condition.R

source(paste0(path,"peakAnalysis/mutate_condition_v2.R"))


df_mutated<- mutate_condition(df=df,conditions=some_conditions,conditions2=other_conditions,varname1=.varname1,varname2=NULL)
