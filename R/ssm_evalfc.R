ssm_evalfc <-
function(x,x_L,x_U,fobj,nconst,c_L,c_U,tolc,weight,int_var,bin_var,nvar,...){
	includ <- 1
	
	if (int_var | bin_var){
		x <- ssm_round_int(x,int_var+bin_var,x_L,nvar);
	}
	
	
#Adjust to bounds
	lower <- which(x<x_L);
	upper <- which(x>x_U);
	
	x[lower]=x_L[lower];
	x[upper]=x_U[upper];
	
	
	if (nconst){
		output <- do.call(fobj,list(x,...));
		value <- output[[1]];
		nlc <- output[[2]];
		
		pena <- ssm_penalty_function(nlc,c_L,c_U,tolc);
		value_penalty <- value+weight*pena;
		
	} else{
		output<-do.call(fobj,list(x,...));
		value <- output[[1]];
		nlc <- 0;
		pena <- 0;
		value_penalty <- value;
	}
	
	if (is.infinite(value) | is.nan(value)){
		includ <- 0;
	}
	return(list(value,value_penalty,pena,nlc,includ,x));
}

