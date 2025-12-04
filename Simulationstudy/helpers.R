coverage <- function(my, my_fit, proc = 0.95){
  if(length(my)==1){
    out <- ifelse(my > quantile(my_fit, (1 - proc)/2) & 
                    my < quantile(my_fit, 1 - (1 - proc)/2), 1, 0)
  } else {
    out <- mean(ifelse(my > quantile(my_fit, (1 - proc)/2) & 
                         my < quantile(my_fit, 1 - (1 - proc)/2), 1, 0))
  }
  
  return(out)
}  


