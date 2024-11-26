#Load libraries 

library(rje)
library(parallel)
library(doParallel)
library(foreach)

#Import functions 
source('ibm_functions_modified.R')

# General functions

## creta_incidence_matrx: function to create the incidence matrix from data 
create_incidence_matrix <- function(data){
  
  tot_pop <- data[1, "S"] + data[1, "I"] + data[1, "R"]
  s <- c(data["S"])
  i <- c(data[1, "I"], rep(0, (nrow(data) - 1)))
  r <- c(data[1, "R"], rep(0, (nrow(data) - 1)))
  kids <- c(data["Kids"])
  adults <- c(data["Adults"])
  incidence_matrix <- data.frame("S" = s, "I" = i, "R" = r, "Kids" = kids, "Adults" = adults)
  
  for (i in 2: nrow(data)){
    
    incidence_matrix[i, "I"] <- -1 *(data[i, "S"] - data[(i - 1), "S"])
    incidence_matrix[i, "Kids"] <- -1 *(data[i, "Kids"] - data[(i - 1), "Kids"])
    incidence_matrix[i, "Adults"] <- -1 *(data[i, "Adults"] - data[(i - 1), "Adults"])
    incidence_matrix[i, "R"] <- (data[i, "R"] - data[(i - 1), "R"])
    
  }
  
  return(incidence_matrix)
  
}


## run_IBM_seed: function to run_ibm_location() at different seeds

run_IBM_seed <- function(seed, parm){
  
  data <- run_ibm_location(rng_seed = seed, transmission_prob = parm)

## compute_likelihood: function to compute the likelihood from data. For now just Poisson and Negative Binomial

compute_likelihood <- function (nuissance_parms = NULL, data, type = 'pois', parms) {
  
  if (type == 'pois'){
    
    data_lik <- dpois(data, lambda = parms)
    
  } else {
    
    data_lik <- dnbinom(data, mu = parms, size = nuissance_parms)
    
  }
  
  ll <- sum(log(data_lik + 0.0001))
  #print(ll)
  ll <- exp(ll)
    
  return(ll)
  
}
  return(data)
}

## overdispersion_nb: function to compute the nuissance overdispersion parameter for the NB 

overdispersion_nb <- function(size = 1, mu, data, upper_nb = 20){
  
  #call the compute_likelihood while rewritten to have the additional parms as first argument so you can use that direclty in the optimisation 
  
  #optim with the data as the original data and the mu as the simulation ones 
  values <- optim(par = size, compute_likelihood, parms = mu, data = data, method = 'Brent', lower = 0, upper = upper_nb, hessian = T)
  
  #return the solution and the overdispersion parameter
  return(list('value' = values$par, 'convergence' = values$convergence, 'hessian' = values$hessian))
  
}


## likleihood_call: compute the IBM and call the likleihood function with the appropriate function parameters. Can return the negative loglikelihood or the likelihood. Function runs in parallel

likelihood_call <- function(parms, data, average_likelihood = T, negativell = T, transformed = T, type = 'pois', additional_parms = NULL, seeds, cl_num = 5) {

  print(parms)
  
  if (transformed){
    
    parms <- 1/(1 + exp(-parms))
    
  }

  #(windows based)
  #create the cluster and run the simulations 
  require(parallel)
  cl <- parallel::makeCluster(cl_num)
  
  #share the environment to the cores (Necessary for Windows based)
  clusterExport(cl, c("create_population_matrix", "run_ibm_location"))
  
  simulated_data <- parLapply(cl, seeds, run_IBM_seed, parm = parms)
  
  stopCluster(cl)

  #(linux based)

  # insert
  
  #compute the likelihoods and the functions of them 
  
 if (average_likelihood){

    #store lieklihood values
    likelihoods <- c()
    
    for (i in seeds){
  
      simulated_incidence <- create_incidence_matrix(simulated_data[[i]])

      #check the type of the likelihood 
      
      if(type != 'poisson' & !is.null(additional_parms)){
        
        overdispersion <- overdispersion_nb(additional_parms, mu = simulated_incidence$I, data = data$I, upper_nb = additional_parms)
        
        likelihoods <- c(likelihoods, compute_likelihood( data = data$I, parms = simulated_incidence$I, type = type, nuissance_parms = overdispersion$value))
        
      } else {
      
        #compute the likelihood for this run 
      likelihoods <- c(likelihoods, compute_likelihood( data = data$I, parms = simulated_incidence$I, type = type))
        
      }
      
    }
    
    output <- mean(likelihoods)
    
  } else {
    
    #a row a time point, a column an Incidence record for a run of the IBM 
    data_records <- as.data.frame(matrix(nrow = nrow(data), ncol = length(seeds)))
    
    for (i in seeds){
  
    simulated_incidence <- create_incidence_matrix(simulated_data[[i]])
    
    data_records[, i] <- simulated_incidence$I
  
    }
    
    average_incidence <- colMeans(t(data_records))
    
    if(type != 'poisson' & !is.NULL(additional_parms)){
        
        overdispersion <- overdispersion_nb(additional_parms, mu = average_incidence, data = data$I, upper_nb = additional_parms)
        
        likelihoods <- c(likelihoods, compute_likelihood( data = data$I, parms = average_incidence, type = type, nuissance_parms = overdispersion$value))
        
      } else {
        
        #compute the likelihood with the mean data 
        output <- compute_likelihood(data = data$I, parms = average_incidence, type = type)
        
      }
    
  }
  
  if (negativell){

     output <- -1*log(output)
  }
  
  return(output)
  
}

# Optimisation 

## design_optimiser: call the optimisation procedure over the likelihood_call

design_optimiser <- function(value, data, seeds, average_likelihood  = T, negativell = T, type = 'pois', additional_parms = NULL, cl_num = 5,  method = "Nelder-Mead", lower = -Inf, upper = Inf){
  
  require(rje)
  
  init <- Sys.time()
  results <- optim(c(logit(value)), likelihood_call, data = data, seeds = seeds, average_likelihood = average_likelihood, negativell = negativell, transformed = T, type = type, additional_parms = additional_parms, cl_num = cl_num,  method = c("Nelder-Mead"), hessian = T, lower = lower, upper = upper)
  end <- Sys.time()
  
  data_collect <- c(length(seeds), results$convergence, results$par, results$hessian, difftime(end, init, unit = "mins"))
  
  return(data_collect)
  
}

## call_optimiser: invokes the design_optimser in a parallelised way or not 
call_optimiser <- function(data, ids, values, method = 'Nelder-Mead', seeds, parallel = T, main_cluster = 2, average_likelihood = T, negativell = T, type = 'pois',  additional_parms = NULL,  cl_num = 3, lower = -Inf, upper = Inf){
  
  if(parallel){
    
    #export the environment 
    export = c('optim', 'likelihood_call', 'compute_likelihood', 'run_IBM_seed', 'run_ibm_location', 'create_population_matrix', 'design_optimiser', 'create_incidence_matrix', 'overdispersion_nb')
    
    my.cluster <- makeCluster(2)
    registerDoParallel(cl = my.cluster)
    clusterExport(my.cluster, export)
    
    values_returned <- foreach(i = 1:length(ids), .combine = "rbind") %dopar% {c(ids[i], design_optimiser(values[i], data = data[data$ID == ids[i], ], seeds = seeds, average_likelihood  = average_likelihood, negativell = negativell, type = type, additional_parms = additional_parms, cl_num = cl_num,  method = method, lower = lower, upper = upper)) }
  
    stopCluster(my.cluster)
    
  } else {
    
    values_returned <- c(ids[1], design_optimiser(value = values[1], data = data[data$ID == ids[1], ], seeds = seeds, average_likelihood  = average_likelihood, negativell = negativell, type = type, additional_parms = additional_parms,  cl_num = cl_num,  method = method, lower = lower, upper = upper))
    
    for (i in 2:length(ids)){
      
     values_returned  <- rbind(values_returned, c(ids[i], design_optimiser(value = values[i], data = data[data$ID == ids[i], ], seeds = seeds, average_likelihood  = average_likelihood, negativell = negativell, type = type, additional_parms = additional_parms, cl_num = cl_num,  method = method, lower = lower, upper = upper)))
         
    }
    
  }
  
  colnames(values_returned) <- c('ID', 'nruns', 'convergence', 'input parm', 'hessian', 'mins')
  
  return(values_returned)
  
}

 

# Other function callers

## subset_likleihood_call: directly calls likelihood_call in a parallelised way 

subset_likleihood_call <- function(parms, data, seeds, cl_num, main_cluster = 2, average_likelihood = T, negativell = T, transformed = T, type = 'pois', additional_parms = NULL){
  
  export = c('likelihood_call', 'compute_likelihood', 'run_IBM_seed', 'run_ibm_location', 'create_population_matrix', 'create_incidence_matrix', 'overdispersion_nb')
    
    my.cluster <- makeCluster(main_cluster)
    registerDoParallel(cl = my.cluster)
    clusterExport(my.cluster, export)
    
    values <- foreach(i = 1:length(parms), .combine = "rbind") %dopar% {c(parms[i], likelihood_call(parms[i], data, average_likelihood = average_likelihood, negativell = negativell, transformed = transformed, type = type, additional_parms = additional_parms,  seeds = seeds, cl_num = cl_num)) }
  
    stopCluster(my.cluster)
    
    return(values)
  
}

## compute_hessian: invoke the numDeriv::hessian to approxiamte the hessian at the value provided by the optim function (or whatever value provided through the parms parmater). Parallelised. 

compute_hessian <- function(inputdata, parms, seeds, cl_num, main_cluster, average_likelihood = T, negativell = T, type = 'pois',  additional_parms = NULL){
  
  export = c('hessian', 'likelihood_call', 'compute_likelihood', 'run_IBM_seed', 'run_ibm_location', 'create_population_matrix', 'create_incidence_matrix')
    
    my.cluster <- makeCluster(main_cluster)
    registerDoParallel(cl = my.cluster)
    clusterExport(my.cluster, export)

  ## chaneg tje +99 doen just for experiments 
    values <- foreach(i = 1:length(parms), .combine = "rbind") %dopar% {c(parms[i], unique(inputdata[inputdata$ID == (99+i), ]$ID), hessian(likelihood_call, parms[i], data = inputdata[inputdata$ID == (99+i), ], average_likelihood = T, negativell = T, transformed = T, type = 'pois', additional_parms = NULL, seeds = seeds, cl_num = cl_num)) }
  
    stopCluster(my.cluster)
    
    return(values)
  
}

# Analysis 
