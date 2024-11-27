# Optimisation methods results 
## Multiple datastreams on different resplications 


countNegative <- function(x){
  
 if ( sum(x < 0) > 0){
   
   return(F)
   
 } 
  
  return(T)
  
}



countZero <- function(x){
  
 if ( sum(x == 0) > 0){
   
   return(F)
   
 } 
  
  return(T)
  
}


removal_not_converged <- function(data1, data2, data3){
  
  data <- data.frame('1run' = data1$convergence, '5run' = data2$convergence, '10run' = data3$convergence)
  
  values <- apply(data, FUN = function(x) if(sum(x != 0) > 0){return(F)}else{return(T)}, MARGIN = 1)
  index <- which(values == T)
  
  
  return(index)
  
}



removal_hessian <- function(data1, data2, data3){
  
  hessians <- data.frame('1run' = data1$hessian, '5run' = data2$hessian, '10run' = data3$hessian)
  
  values1 <- apply(hessians, FUN = countNegative, MARGIN = 1)
  values2 <- apply(hessians, FUN = countZero, MARGIN = 1)
  
  index <- which(values1 == T & values2 == T) #good for both 
  index_negative <- which(values1 == F) 
  index_zero <- which(values2 == F)
  
  new_hessians <- hessians[index, ] #good for both
  
  return(list('indeces' = index, 'negative_indeces' = index_negative, 'zero_indeces' = index_zero,  'hessian' = new_hessians))
  
}

#assume a standard normal distribution as a limiting distribution
confidence_interval <- function(data, variance, prob = 0.95){
  
  prob_final <- prob + ((1 - prob)/2)
  
  CI_lower <- data - qnorm(prob_final, 0, 1)*sqrt(variance)
  CI_upper <- data + qnorm(prob_final, 0, 1)*sqrt(variance)
  
  return(list('lower' = CI_lower, 'upper' = CI_upper))
  
}

compute_analysis <- function(data1, data2, data3, true_parm = 0.1){
  
  true_parm <- logit(true_parm)
  
  ids <- data1$ID
  
  #sum convergence 
  indeces <- removal_not_converged(data1, data2, data3)
  data1 <- data1[indeces, ]
  data2 <- data2[indeces, ]
  data3 <- data3[indeces, ]
  
  nc <- nrow(data1) - length(indeces)
  
  
  #removal of bad hessians 
  track <- removal_hessian(data1, data2, data3)
  
  negative <- length(track$negative_indeces)
  zero <- length(track$zero_indeces)
  
  parms <- cbind(data1[track$indeces, ]$input.parm, data2[track$indeces, ]$input.parm, data3[track$indeces, ]$input.parm)
  

  #mean and empirical SE
  means <- apply(parms, FUN = mean, MARGIN = 2)
  var <- apply(parms, FUN = var, MARGIN = 2)
  
  #theoretical mean standard error (one parm only, attention for the matrix ). use 1/hessians as each entry is the hessian of the parm and this is the same as inverting 1x1 matrix
  hessians <- track$hessian
  msde <- sqrt(apply(1/(hessians), MARGIN = 2, FUN = mean))
  
  #95% confidence intervals (normal)
  CI_95_hdata1 <- confidence_interval(parms[, 1], 1/hessians[, 1])
  CI_95_hdata2 <- confidence_interval(parms[, 2], 1/hessians[, 2])
  CI_95_hdata3 <- confidence_interval(parms[, 3], 1/hessians[, 3])
  
  CI_95_data1 <- confidence_interval(parms[, 1], var[1])
  CI_95_data2 <- confidence_interval(parms[, 2], var[2])
  CI_95_data3 <- confidence_interval(parms[, 3], var[3])
  
  CI_95h <- data.frame('Lower1' = CI_95_hdata1$lower, 'Upper1' = CI_95_hdata1$upper,'Lower5' = CI_95_hdata2$lower, 'Upper5' = CI_95_hdata2$upper, 'Lower10' = CI_95_hdata3$lower, 'Upper10' = CI_95_hdata3$upper )
  CI_95 <- data.frame('Lower1' = CI_95_data1$lower, 'Upper1' = CI_95_data1$upper,'Lower5' = CI_95_data2$lower, 'Upper5' = CI_95_data2$upper, 'Lower10' = CI_95_data3$lower, 'Upper10' = CI_95_data3$upper )
  
  #90% confidence intervals (normal)
  CI_90_hdata1 <- confidence_interval(parms[, 1], 1/hessians[, 1], prob = 0.9)
  CI_90_hdata2 <- confidence_interval(parms[, 2], 1/hessians[, 2], prob = 0.9)
  CI_90_hdata3 <- confidence_interval(parms[, 3], 1/hessians[, 3], prob = 0.9)
  
  CI_90_data1 <- confidence_interval(parms[, 1], var[1], prob = 0.9)
  CI_90_data2 <- confidence_interval(parms[, 2], var[2], prob = 0.9)
  CI_90_data3 <- confidence_interval(parms[, 3], var[3], prob = 0.9)
  
  CI_90h <- data.frame('Lower1' = CI_90_hdata1$lower, 'Upper1' = CI_90_hdata1$upper,'Lower5' = CI_90_hdata2$lower, 'Upper5' = CI_90_hdata2$upper, 'Lower10' = CI_90_hdata3$lower, 'Upper10' = CI_90_hdata3$upper )
  CI_90 <- data.frame('Lower1' = CI_90_data1$lower, 'Upper1' = CI_90_data1$upper,'Lower5' = CI_90_data2$lower, 'Upper5' = CI_90_data2$upper, 'Lower10' = CI_90_data3$lower, 'Upper10' = CI_90_data3$upper )
  
  #coverage: how many times is the true value within the CI? is the nominal coverage reached? with mean scross the simulations alias the data_streams 
  coverage_1_95 <- mean(true_parm >= CI_95_hdata1$lower & true_parm <= CI_95_hdata1$upper)
  coverage_5_95 <- mean(true_parm >= CI_95_hdata2$lower & true_parm <= CI_95_hdata2$upper)
  coverage_10_95 <- mean(true_parm >= CI_95_hdata3$lower & true_parm <= CI_95_hdata3$upper)
  coverage_95_h <- c(coverage_1_95, coverage_5_95, coverage_10_95)
  
  coverage_1_90 <- mean(true_parm >= CI_90_hdata1$lower & true_parm <= CI_90_hdata1$upper)
  coverage_5_90 <- mean(true_parm >= CI_90_hdata2$lower & true_parm <= CI_90_hdata2$upper)
  coverage_10_90 <- mean(true_parm >= CI_90_hdata3$lower & true_parm <= CI_90_hdata3$upper)
  coverage_90_h <- c(coverage_1_90, coverage_5_90, coverage_10_90)
  
  
  coverage_1_95 <- mean(true_parm >= CI_95_data1$lower & true_parm <= CI_95_data1$upper)
  coverage_5_95 <- mean(true_parm >= CI_95_data2$lower & true_parm <= CI_95_data2$upper)
  coverage_10_95 <- mean(true_parm >= CI_95_data3$lower & true_parm <= CI_95_data3$upper)
  coverage_95 <- c(coverage_1_95, coverage_5_95, coverage_10_95)
  
  coverage_1_90 <- mean(true_parm >= CI_90_data1$lower & true_parm <= CI_90_data1$upper)
  coverage_5_90 <- mean(true_parm >= CI_90_data2$lower & true_parm <= CI_90_data2$upper)
  coverage_10_90 <- mean(true_parm >= CI_90_data3$lower & true_parm <= CI_90_data3$upper)
  coverage_90 <- c(coverage_1_90, coverage_5_90, coverage_10_90)
  
  
  #Mean-squared-error
  mse <- apply(parms, MARGIN = 2, FUN = function(x, true_parm) mean((x - true_parm)^2), true_parm = true_parm)
  
  
  stats <- data.frame('mean' = means, 'ese' = sqrt(var), 'tse' =  msde, 'mse' = mse, 'coverage95h' = coverage_95_h ,'coverage90h' = coverage_90_h, 'coverage95' = coverage_95, 'coverage90' = coverage_90 )
  
  output <- list('stats' = stats, 'not_converged' = nc, 'negative_hessian' = ids[track$negative_indeces], 'zero_hessian' =  ids[track$zero_indeces], 'CI_95h' = CI_95h, 'CI_95' = CI_95, 'CI_90h' = CI_90h, 'CI_90' = CI_90)
  
  return(output)
  
}

