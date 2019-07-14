build_data <- function(probA = 0.1, probB = 0.1, size = 1000){
  contingency <- data.frame()
  disease_A <- rbinom(n = 1, size = size, prob = probA)
  no_disease_A <- size - disease_A
  contingency <- rbind(contingency, data.frame("disease" = disease_A, 
                                               "no_disease" = no_disease_A))
  disease_B <- rbinom(n = 1, size = size, prob = probB)
  no_disease_B <- size - disease_B
  contingency <- rbind(contingency, data.frame("disease" = disease_B, 
                                               "no_disease" = no_disease_B))
  rownames(contingency) <- c("A", "B")
  return(contingency)
}

extract_pvalue <- function(contingency){
  return(chisq.test(contingency)$p.value)
}
  
run_many_times <- function(probA = 0.1, probB = 0.2, size = 1000, 
                            num_replicates = 500, alpha = 0.05){
  all_pvalues <- replicate(n = num_replicates,
                           extract_pvalue(build_data(probA, probB, size)))
  return(data.frame(pvalue = all_pvalues, significant = all_pvalues < alpha))
}  
  
  