vns_defaults <- function(...){
  #Assings default values for all the options  
  maxeval                   <-       1000;            #Maximum number of function evaluations
  maxtime                   <-       60;              #Maximum CPU time
  maxdist                   <-       0.5;             #Percentage of the problem dimension which will be perturbed in the furthest neighborhood (vary between 0-1)
  use_local                 <-       1;               #Uses local search (1) or not (0)
  #The following options only apply if use_local is active
  aggr                    <-       0;               #Aggressive search. The local search is only applied when the best solution has been improved (1=aggressive search, 0=non-aggressive search)
  local_search_type	      <-       1;               #Applies a first (=1) or a best (=2) improvement scheme for the local search
  decomp                  <-       1;               #Decompose the local search (=1) using only the variables perturbed in the global phase
  return(list(maxeval=maxeval, maxtime=maxtime, maxdist=maxdist,use_local=use_local,aggr=aggr,local_search_type=local_search_type,
              decomp=decomp))
}

