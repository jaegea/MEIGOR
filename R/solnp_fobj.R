solnp_fobj <-
function(x, ...){
	res<-do.call(fobj_global, list(x, ...));
	n_fun_eval <<- n_fun_eval + 1;
	res_f<-res[[1]];
	return(res_f)
}

