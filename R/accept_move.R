#' @title Accept the current proposed move
#' @description "accept_move" updates the output list after accepting the current proposed move.
#' 
#' @param output: List of output values with entries as explained below. 
#' @param num_estimate: int. 
#'    Number of parameters to estimate.
#' @param estimate_idx: list of int.  
#'    Indices of parameters to estimate in the model's full parameter list.
#' @param initial_values: list of float.  
#'    Starting values for parameters to estimate, taken from the parameters'
#'    nominal values in the model or explicitly specified in 'options'.
#' @param initial_position: list of float. 
#'    Starting position of the MCMC walk in parameter space (log10 of 'initial_values').
#' @param position: list of float. 
#'    Current position of MCMC walk in parameter space, i.e. the most
#'    recently accepted move.
#' @param test_position: list of float. 
#'    Proposed MCMC mmove.
#' @param acceptance: int. 
#'    Number of accepted moves.
#' @param T: float.  
#'    Current value of the simulated annealing temperature.
#' @param T_decay: float.  
#'    Constant for exponential decay of 'T', automatically calculated such 
#'    that T will decay from 'options$T_init' down to 1 over the first
#'    'options$anneal_length' setps.
#' @param sig_value: float.  
#'    Current value of sigma, the scaling factor for the proposal distribution.
#'    The MCMC algorithm dynamically tunes this to maintain the aaceptance
#'    rate specified in 'options$accept_rate_target'.
#' @param iter: int.  
#'    Current MCMC step number.
#' @param start_iter: int.  
#'    Starting MCMC step number.
#' @param ode_options: list.  
#'    Options for the ODE integrator, currently just 'rtol' for relative
#'    tolerance and 'atol' for absolute tolerance.
#' @param initial_prior: float. 
#'    Starting prior value, i.e. the value at 'initial_position'.
#' @param initial_likelihood: float.  
#'    Starting likelihood value, i.e. the value at 'initial_position'.
#' @param initial_posterior: float.  
#'    Starting posterior value, i.e. the value at 'initial_position'.
#' @param accept_prior: float.  
#'    Current prior value i.e. the value at 'position'.
#' @param accept_likelihood: float.  
#'    Current likelihood value i.e. the value at 'position'.
#' @param accept_posterior: float.  
#'    Current posterior value i.e. the value at 'position'.
#' @param test_prior: float.  
#'    Prior value at 'test_position'.
#' @param test_likelihood: float.  
#'    Likelihood value at 'test_position'.
#' @param test_posterior: float.  
#'    Posterior value at 'test_position'.
#' @param hessian: array of float. 
#'    Current hessian of the posterior landscape. Size is 
#'    'num_estimate' x 'num_estimate'.
#' @param positions: array of float.  
#'    Trace of all proposed moves. Size is 'num_estimate' x 'nsteps'.
#' @param priors: array of float. 
#'    Trace of all priors corresponding to 'positions'. Length is 'nsteps'.
#' @param likelihoods: array of float. 
#'    Trace of all likelihoods corresponding to 'positions'. Length is 'nsteps'.
#' @param posteriors: array of float. 
#'    Trace of all posteriors corresponding to 'positions'. Length is 'nsteps'.
#' @param alphas: array of float. Trace of 'alpha' parameter and calculated values. Length is 'nsteps'.
#' @param sigmas: array of float. Trace of 'sigma' parameter and calculated values. Length is 'nsteps'.
#' @param delta_posteriors: array of float. Trace of 'delta_posterior' parameter and calculated values. Length is 'nsteps'.
#' @param ts: array of float. Trace of 'T' parameter and calculated values. Length is 'nsteps'.
#' @param accepts: logical array. 
#'    Trace of wheter each proposed move was accepted or not.
#'    Length is 'nsteps'.
#' @param rejects: logical array. Trace of wheter each proposed move was rejected or not. Length is 'nsteps'.
#' @param hessians: array of float. 
#'    Trace of all hessians. Size is 'num_estimate' x 'num_estimate' x 'num_hessians'
#'    where 'num_hessians' is the actual number of hessians to be calculated.
#' 
#' 
#' @return The updated output after the move is accepted -- list with entries as explained in 'Arguments'.
#' 
#' 
accept_move = function(output){
  output$accept_prior = output$test_prior
  output$accept_likelihood = output$test_likelihood
  output$accept_posterior = output$test_posterior
  output$position = output$test_position
  output$acceptance = output$acceptance + 1
  output$accepts[output$iter] = TRUE
  return(output)
}