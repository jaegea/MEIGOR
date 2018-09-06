solnp_eq <-
function(x, ...){
	res<-do.call(fobj_global, list(x, ...));
	n_fun_eval <<- n_fun_eval + 1;
	res_eq<-res[[2]][1:neq_global];
	return(res_eq)
}

