#Try to lighten up!

#lighten_up_df.R


lighten_up<-function(df,remove_cols=drops){



			lighter_df <- df[,!(names(df) %in% remove_cols)]

			lighter_df

}
                        