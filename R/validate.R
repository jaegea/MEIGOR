#' @title Validate options.
#' 
#' @description Validate the given options.
#' 
#' 
#' @param options: list with entries as explained below. 
#'    Options set -- defines the problem and sets some
#'    parameters to control the MCMC algorithm.
#' @param model: List of model parameters - to estimate. 
#'    The parameter objects must each have a
#'    'value' attribute containing the parameter's numerical value.
#' @param estimate_params: list. 
#'    List of parameters to estimate, all of which must also be listed
#'    in 'options$model$parameters'.
#' @param initial_values: list of float, optional. 
#'    Starting values for parameters to estimate. If omitted, will use
#'    the nominal values from 'options$model$parameters'
#' @param step_fn: callable f(output), optional. 
#'    User callback, called on every MCMC iteration.
#' @param likelihood_fn: callable f(output, position). 
#'    User likelihood function.
#' @param prior_fn: callable f(output, position), optional. 
#'    User prior function. If omitted, a flat prior will be used.
#' @param nsteps: int. 
#'    Number of MCMC iterations to perform.
#' @param use_hessian: logical, optional. 
#'    Wheter to use the Hessian to guide the walk. Defaults to FALSE.
#' @param rtol: float or list of float, optional. 
#'    Relative tolerance for ode solver.
#' @param atol: float or list of float, optional. 
#'    Absolute tolerance for ode solver.    
#' @param norm_step_size: float, optional. 
#'    MCMC step size. Defaults to a reasonable value.
#' @param hessian_period: int, optional. 
#'    Number of MCMC steps between Hessian recalculations. Defaults
#'    to a reasonable but fairly large value, as Hessian calculation is expensive.
#' @param hessian_scale: float, optional. 
#'    Scaling factor used in generating Hessian-guided steps. Defaults to a
#'    reasonable value.
#' @param sigma_adj_interval: int, optional. 
#'    How often to adjust 'output$sig_value' while annealing to meet
#'    'accept_rate_target'. Defaults to a reasonable value.
#' @param anneal_length: int, optional. 
#'    Length of initial "burn-in" annealing period. Defaults to 10% of
#'    'nsteps', or if 'use_hessian' is TRUE, to 'hessian_period' (i.e. 
#'    anneal until first hessian is calculated)
#' @param T_init: float, optional. 
#'    Initial temperature for annealing. Defaults to a resonable value.
#' @param accept_rate_target: float, optional. 
#'    Desired acceptance rate during annealing. Defaults to a reasonable value.
#'    See also 'sigma_adj_interval' above.
#' @param sigma_max: float, optional. 
#'    Maximum value for 'output$sig_value'. Defaults to a resonable value.
#' @param sigma_min: float, optional. 
#'    Minimum value for 'output$sig_value'. Defaults to a resonable value. 
#' @param sigma_step: float, optional. 
#'    Increment for 'output$sig_value' adjustments. Defaults to a resonable 
#'    value. To eliminate adaptive step size, set sigma_step to 1.
#' @param thermo_temp: float in the range [0,1], optional. 
#'    Temperature for thermodynamic integration support. Used to scale
#'    likelihood when calculating the posterior value. Defaults to 1,
#'    i.e. no effect.
#'    
#' 
#' 
#' @return options: A validated copy of options with defaults applied -- list with entries as explained in 'Arguments'.
#' 
validate = function(options){
  
  if (is.null(options$model)){
    stop("model not defined")
  }
  if ((is.null(options$estimate_params)) || (length(options$estimate_params)==0)){
    stop("estimate_params must contain a list of parameters")
  }
  
  # clamp hessian_period to actual number of steps
  if (options$use_hessian == TRUE){
    options$hessian_period = pmin(options$hessian_period, options$nsteps)
  } else {
    options$hessian_period = Inf
  }
  
  if (is.null(options$anneal_length)){
    # default for anneal length is unspecified
    if (options$use_hessian == TRUE){
      # if using hessian, anneal until we start using it
      options$anneal_length = options$hessian_period
    } else{
      # otherwise, anneal for 10% of the run
      options$anneal_length = options$nsteps*0.10
    }
  } else {
    # clamp it to actual number of steps
    options$anneal_length = pmin(options$anneal_length, options$nsteps)
  }
  
  # default for sigma_adj_interval if unspecified
  if (is.null(options$sigma_adj_interval)){
    # default to 10 adjustments throughout the annealing phase
    options$sigma_adj_interval = pmax((options$anneal_length%/%10), 1)
  }
  return (options)
}