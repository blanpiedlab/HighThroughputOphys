
run_anova <- function(var1, var2, df) {
  fo <- reformulate(var2, var1)
  test<-do.call("aov", list(fo, substitute(df)))
  print(summary(test))
}

run_stats <- function(parent_A,comparison_A, parent_B,comparison_B){#,list_store) {
    print(paste0("wilcox test for : ", var_list[i]))
    test<-wilcox.test(comparison_A,comparison_B,paired=FALSE)
    print(test)
    comparison = paste(unique(parent_A$Ca),unique(parent_B$Ca),sep='-')
    data_wilcox<- data.frame(comparison = comparison,test = 'wilcoxon',p_value = test$p.value)
    print(paste0("ks.test for : ", var_list[i]))
    test<-ks.test(comparison_A,comparison_B)
    print(test)
    data_ks<- data.frame(comparison = comparison,test = 'ks_test',p_value = test$p.value)
    data_bind <- bind_rows(data_wilcox,data_ks)

    data_bind
    


}


library(ggpubr)
A <- peak_stats %>% dplyr::filter(Ca == "0pt5Ca",segmentation == "marker")
B <- peak_stats %>% dplyr::filter(Ca == "1Ca",segmentation == "marker")
C <- peak_stats %>% dplyr::filter(Ca == "2Ca",segmentation == "marker")
D <- peak_stats %>% dplyr::filter(Ca == "4Ca",segmentation == "marker")



var_list<- c('amplitude','tau_decay_ms','t_rise','t_decay','t_half','interSpike_ms')

p_list<- list()



for(i in 1:length(var_list)) {


    #for k in comparisons
    #select(get(var_list[i]))


    tmp_A <- A %>% ungroup() %>% select_(var_list[i])
    tmp_A <- as.matrix(tmp_A)

    tmp_B <- B %>% ungroup() %>% select_(var_list[i])
    tmp_B <- as.matrix(tmp_B)

    tmp_C <- C %>% ungroup() %>% select_(var_list[i])
    tmp_C <- as.matrix(tmp_C)

    tmp_D <- D %>% ungroup() %>% select_(var_list[i])
    tmp_D <- as.matrix(tmp_D)


    #A-B
    #A-C
    #A-D
    #B-C
    #B-D
    #C-D

    print(paste0("Generating a statistical comparison for the following variable: ", var_list[i] ) )
    print("Comparing 0pt5Ca, 1Ca")
    output<-run_stats(A,tmp_A,B,tmp_B)#,p_list)
    #assign(paste0("stats_",var_list[i],"_05-1Ca"),output)

    print("Comparing 0pt5Ca, 2Ca")
    output<-run_stats(A,tmp_A,C,tmp_C)#,p_list)
    #assign(paste0("stats_",var_list[i],"_05-2Ca"),output)

    print("Comparing 0pt5Ca, 4Ca")
    output<-run_stats(A,tmp_A,D,tmp_D)#,p_list)
    #assign(paste0("stats_",var_list[i],"_05-4Ca"),output)
    
    print("Comparing 1Ca, 2Ca")
    output<-run_stats(B,tmp_B,C,tmp_C)#,p_list)
    #assign(paste0("stats_",var_list[i],"_1-2Ca"),output)
    print("Comparing 1Ca, 4Ca")
    
    output<-run_stats(B,tmp_B,D,tmp_D)#,p_list)
    assign(paste0("stats_",var_list[i],"_1-4Ca"),output)

    print("Comparing 2Ca, 4Ca")
  
    output<-run_stats(C,tmp_C,D,tmp_D)#,p_list)
    #assign(paste0("stats_",var_list[i],"_2-4Ca"),output)
    
    check<- peak_stats %>% select_("Ca",var_list[i]) %>% ungroup() 
   
    check$Ca <- factor(check$Ca, levels = c('0pt5Ca','1Ca','2Ca','4Ca'))
    print("TESTING ANOVA")
    print(paste0("current variable = ", var_list[i]))
    run_anova(var_list[i], "Ca", check)

}
