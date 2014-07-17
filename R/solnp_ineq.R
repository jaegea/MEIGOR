solnp_ineq <-
function(x, ...){
	res<-do.call(fobj_global, list(x, ...));
	n_fun_eval <<- n_fun_eval + 1;
	res_ineq<-res[[2]][(neq_global+1):nconst_global]
	return(res_ineq)
}

