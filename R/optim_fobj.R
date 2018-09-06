optim_fobj <-
function(x){
  #Trick to avoid passing lists to nls
  extra_args[[1]]<-x
  f<-do.call(fobj_global, extra_args);
  n_fun_eval <<- n_fun_eval + 1;  
	return(f[1])
}

