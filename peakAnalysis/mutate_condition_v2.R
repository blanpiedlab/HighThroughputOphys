
# mutate_condition_v2.R


mutate_condition<- function(df=df, conditions = some_conditions,conditions2=other_conditions,varname1 = .varname1, varname2 = .varname2) {

				.condition_vars =  rlang::syms(conditions)

				if (!is.null(conditions2)) {
					switch_it = TRUE
					.more_condition_vars =  rlang::syms(conditions2)					

				} else {
					switch_it = FALSE
				}

					

				df_mutated <- df %>% dplyr::mutate(!!.varname1 := paste(!!!.condition_vars, sep="-"))
				
				if (switch_it) {
					df_mutated <- df_mutated %>% dplyr::mutate(!!.varname2 := paste(!!!.more_condition_vars, sep="-"))
				}

				df_mutated
}
