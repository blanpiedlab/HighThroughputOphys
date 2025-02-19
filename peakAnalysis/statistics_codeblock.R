
run_anova <- function(var1, var2, df) {
  fo <- reformulate(var2, var1)
  test<-do.call("aov", list(fo, substitute(df)))
  print(summary(test))
}

run_stats <- function(comparison_A, comparison_B) {
    print(paste0("wilcox test for : ", var_list[i]))
    test<-wilcox.test(comparison_A,comparison_B,paired=FALSE)
    print(test)

    print(paste0("ks.test for : ", var_list[i]))
    test<-ks.test(comparison_A,comparison_B)
    print(test)

    


}


library(ggpubr)
A <- peak_stats %>% dplyr::filter(Ca == "0pt5Ca")
B <- peak_stats %>% dplyr::filter(Ca == "1Ca")
C <- peak_stats %>% dplyr::filter(Ca == "2Ca")
D <- peak_stats %>% dplyr::filter(Ca == "4Ca")



var_list<- c('amplitude','tau_decay_ms','t_rise','t_decay','t_half')



for(i in 1:length(var_list)) {



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
    run_stats(tmp_A,tmp_B)
    

    print("Comparing 0pt5Ca, 2Ca")
    run_stats(tmp_A,tmp_C)
    

    print("Comparing 0pt5Ca, 4Ca")
    run_stats(tmp_A,tmp_D)
    
    print("Comparing 1Ca, 2Ca")
    run_stats(tmp_B,tmp_C)
  
    print("Comparing 1Ca, 4Ca")
    
    run_stats(tmp_B,tmp_D)
    
    print("Comparing 2Ca, 4Ca")
  
    run_stats(tmp_C,tmp_D)
    

    check<- peak_stats %>% select_("Ca",var_list[i]) %>% ungroup()
    print("TESTING ANOVA")
    print(paste0("current variable = ", var_list[i]))
    run_anova(var_list[i], "Ca", check)

}
