#' @title Running the BayesFit optimisation
#' @description "runBayesFit" defines the prior function and runs the BayesFit estimation

#' @param opts list with entries as explained below. 
#'    Options set -- defines the problem and sets some
#'    parameters to control the MCMC algorithm.
#'  model: List of model parameters - to estimate. 
#'    The parameter objects must each have a
#'    'value' attribute containing the parameter's numerical value.
#'  estimate_params: list. 
#'    List of parameters to estimate, all of which must also be listed
#'    in 'options$model$parameters'.
#'  initial_values: list of float, optional. 
#'    Starting values for parameters to estimate. If omitted, will use
#'    the nominal values from 'options$model$parameters'
#'  step_fn: callable f(output), optional. 
#'    User callback, called on every MCMC iteration.
#'  likelihood_fn: callable f(output, position). 
#'    User likelihood function.
#'  prior_fn: callable f(output, position), optional. 
#'    User prior function. If omitted, a flat prior will be used.
#'  nsteps: int. 
#'    Number of MCMC iterations to perform.
#'  use_hessian: logical, optional. 
#'    Wheter to use the Hessian to guide the walk. Defaults to FALSE.
#'  rtol: float or list of float, optional. 
#'    Relative tolerance for ode solver.
#'  atol: float or list of float, optional. 
#'    Absolute tolerance for ode solver.    
#'  norm_step_size: float, optional. 
#'    MCMC step size. Defaults to a reasonable value.
#'  hessian_period: int, optional. 
#'    Number of MCMC steps between Hessian recalculations. Defaults
#'    to a reasonable but fairly large value, as Hessian calculation is expensive.
#'  hessian_scale: float, optional. 
#'    Scaling factor used in generating Hessian-guided steps. Defaults to a
#'    reasonable value.
#'  sigma_adj_interval: int, optional. 
#'    How often to adjust 'output$sig_value' while annealing to meet
#'    'accept_rate_target'. Defaults to a reasonable value.
#'  anneal_length: int, optional. 
#'    Length of initial "burn-in" annealing period. Defaults to 10% of
#'    'nsteps', or if 'use_hessian' is TRUE, to 'hessian_period' (i.e. 
#'    anneal until first hessian is calculated)
#'  T_init: float, optional. 
#'    Initial temperature for annealing. Defaults to a resonable value.
#'  accept_rate_target: float, optional. 
#'    Desired acceptance rate during annealing. Defaults to a reasonable value.
#'    See also 'sigma_adj_interval' above.
#'  sigma_max: float, optional. 
#'    Maximum value for 'output$sig_value'. Defaults to a resonable value.
#'  sigma_min: float, optional. 
#'    Minimum value for 'output$sig_value'. Defaults to a resonable value. 
#'  sigma_step: float, optional. 
#'    Increment for 'output$sig_value' adjustments. Defaults to a resonable 
#'    value. To eliminate adaptive step size, set sigma_step to 1.
#'  thermo_temp: float in the range [0,1], optional. 
#'    Temperature for thermodynamic integration support. Used to scale
#'    likelihood when calculating the posterior value. Defaults to 1,
#'    i.e. no effect.
#' @return The output after the optimisation is finished - a list with entries as explained in 'Arguments'.
runBayesFit <- function(opts) {
  
  print_fit = function(position) {
    new_values = 10^position
    cat("param\t", "actual\t", "fitted\t", "% error\t","\n")
    for(a in 1:length(new_values)) {
      error = abs(1-opts$model$parameters$value[a]/new_values[a])*100
      values = c(opts$model$parameters$name[a],"\t",round(opts$model$parameters$value[a],2),"\t",round(new_values[a],2),"\t",round(error,2))
      cat(values,"\n")
    }
  }
  
  step = function(output) {
    if (output$iter%%20==0) {
      cat("iter=", output$iter, " sigma=", output$sig_value, " T=", output$T,
          " acc=", round(output$acceptance/(output$iter+1),3),
          " lkl=", round(output$accept_likelihood,2),
          " prior=", round(output$accept_prior,3),
          " post=", round(output$accept_posterior,2), "\n", sep="")
    }
  }
  
  prior = function(position) {
    s = sum((position-prior_mean)^2/(2*prior_var)) 
    return(s)
  }
  
  cat("Running BayesFit", "\n")
  cat("=========================", "\n")
  
  opts$prior_fn=prior
  opts$step_fn = step
  opts = validate(options=opts)
  out = list()
  out = initialize(output=out, options=opts)
  out = estimate(output=out, options=opts)
  print_fit(out$positions[ opts$nsteps-1,])
  return(out)
}