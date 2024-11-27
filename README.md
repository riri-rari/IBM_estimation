# IBM

## IBM_function_methods.R

This code is aimed at running a desired optimisation procedure invoking the R optim() function to obtain a MLE for a Poisson or a Negative Binomial random variable (extensions to be done). The IBM is provided throught the source 'ibm_functions_modified.R'. The entire process can be parallelised in multiple points to speed up the computations.

The output of 'call_optimiser' is a dataframe with the number of runs/times the IBM has been called, the convergence status, the optimal parmater found, the hessian at that parameter and the time needed. 

The output of 'subset_likleihood' is a matrix with the parameter and the objective function value at that parameter.

The output of 'compute_hessian' is a matrix with the parmater, the ID of the data set against which the IBM is run and the objective function value at that parameter. This function is used to contrast the hessian reported by optim() with the hessian reported by numDeriv::hessian()

## IBM_reading_results.R

```R
# Example 
AvLik_NM_1_runs <- read.csv('./Data/AvLik_Pois_NM_rep1.csv')
AvLik_NM_5_runs <- read.csv( './Data/AvLik_Pois_NM_rep5.csv')
AvLik_NM_10_runs <- read.csv( './Data/AvLik_Pois_NM_rep10.csv')
AvLik_NM_Pois <- compute_analysis(AvLik_NM_1_runs, AvLik_NM_5_runs, AvLik_NM_10_runs)
```

## Data 

## Analyse data 
