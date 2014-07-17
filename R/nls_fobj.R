nls_fobj <-
function(x){
  
  #Trick to avoid passing lists to nls
  extra_args[[1]]<-x
	f<-do.call(fobj_global, extra_args);
	n_fun_eval <<- n_fun_eval + 1;
	
  #Return the 3rd argument which should be the vector of residuals
  
  residuals<-as.vector(f[[3]])
  #return(residuals)
}

