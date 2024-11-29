# IBM

## IBM_function_methods.R

This code is aimed at running a desired optimisation procedure invoking the R optim() function to obtain a MLE for a Poisson or a Negative Binomial random variable (extensions to be done). The IBM is provided throught the source 'ibm_functions_modified.R'. The entire process can be parallelised in multiple points to speed up the computations.

Note that to obtain the data from 'ibm_functions_modified.R' and convert to the incidence matrix by 'create_incidence_matrix' it is possible to run

```R
data <- run_ibm_location(rng_seed = 2, transmission_prob = 0.7) #default transmission_prob = 0.1 
incidence_data <- create_incidence_matrix(data)
```

The output of 'call_optimiser' is a dataframe with the number of runs/times the IBM has been called, the convergence status, the optimal parmater found, the hessian at that parameter and the time needed. 

The output of 'subset_likleihood' is a matrix with the parameter and the objective function value at that parameter.

The output of 'compute_hessian' is a matrix with the parmater, the ID of the data set against which the IBM is run and the objective function value at that parameter. This function is used to contrast the hessian reported by optim() with the hessian reported by numDeriv::hessian()

Note that each time the IBM is called, the seed can be specified. For the dataset available the seeds were 1 (for 1 run), 1 to 5 (for 5 runs) and 1 to 10 (for 10 runs). 

## IBM_reading_results.R

This script contains the function 'compute_analysis' that given 3 datasets in input gives in output the mean value of the solution found, the empirical standard error, the theoretical standard error (the sd of the mean of the information matrix computed as the inverse of the hessian), the mean squared error, the coverage at 90% and 95% using the empirical or the theoretical standard error and the relative confidence intervals for each of the datasets. It also returns the indeces of the negative or zero values hessians and the not converged runs. Note that the statistics are computed after filtering for the appropriate runs (converged and positive hessian) across all three datasets. 

Here an example of the call to the function. 

```R
# Example 
AvLik_NM_1_runs <- read.csv('./Data/AvLik_Pois_NM_rep1.csv')
AvLik_NM_5_runs <- read.csv( './Data/AvLik_Pois_NM_rep5.csv')
AvLik_NM_10_runs <- read.csv( './Data/AvLik_Pois_NM_rep10.csv')
AvLik_NM_Pois <- compute_analysis(AvLik_NM_1_runs, AvLik_NM_5_runs, AvLik_NM_10_runs)
```

## Data 

Brief description of the available datasets: 
- incidence_data_stream_01_100runs100199.csv is the incidence datset for datastreams with seeds 100 to 199 and transmission prpbabiltiy of 0.1
- AvLik_collect_Poisson_Data_stream_ID_144OptLHS10015_1runs.csv is the dataset with the values of the Average Likelihood on 1 run from Latin Hypercube sampling (LHS) datapoints for each dataset (100 to 144). LHS done with minimisation of the S metric. LHS seed 12. Variables: X (call for computation), Parm (parameter), Average_Likelihood, ID (ID dataset)
-   AvLik_collect_Poisson_Data_stream_ID_144OptLHS10015_1runs_optvalues.csv is the datset with the optimal value among the LHS for each dataset. Variables: X (call for comoutation), ID (dataset), optvalue (optimal value)
- AvLik_Pois_method_numrep or LikAv_Pois_method_numrep.csv are the datsets with the optimal values found with the 'method' on an Average Likleihood or Likelihood of the Average run on 'numrep' IBM calls. Variables: ID (dataset), nruns (number of runs the IBM was called on), convergence (convergence status, should be 0), input.parm (optimal parameter), hessian (hessian value as reported by optim), mins (time needed to compute)
- AvLik_Pois_method_hessian_numrep and LikAv_Pois_method_hessian_numrep.csv is the dataset that for each ID and optimal parm (as reported in the AvLik_Pois_method_numrep) reports the hessian computed with numDervi::hessian function. Variables: X (call for computation), Input.parm (optimal parameter), ID (dataset ID), Hessian (hessian value). For now, only available for Nelder-Mead method on 1 and 5 runs (Poisson)
- AvLik_Pois_method_comparativehessian_finergrid_numrep.csv is the dataset with the values of the Average Likelihood on 5 runs around the optimal solution (as found by optim) for those ID datasets that have discordant hessian between optim and numDeriv::hessian. Variables: Call (call for function), input.parm (parameter for the IBM), Avlik (Average Likelihood across 5 runs), Code (F for finer grid points and L for less finer grid points), ID (dataset ID)
- AvLik_Pois_NM_comparativehessian_numrep.csv is a dataset similar to AvLik_Pois_method_comparativehessian_finergrid_numrep.csv but on a broader range but not so fine grid points

## utils 

An R script with some code used to check some outputs. 
