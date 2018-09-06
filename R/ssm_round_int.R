ssm_round_int <-
function(x,index,x_L,nvar){

	index <- nvar-index+1;
	index <- index:nvar;
	
	temp=dim(x);
	nsol=temp[1];

	if (is.vector(x)){
		x<-matrix(x,1,nvar);
		nsol=1;
	}
	
	for (i in 1:nsol) {
		 x[i,index]=x_L[index]+floor(0.5+x[i,index]-x_L[index]);
	}
	
	xrounded <- x;
	return(xrounded)
}

