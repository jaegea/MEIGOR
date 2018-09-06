rosen10<-function(x){
	f<-0;
	n=10;
	for (i in 1:(n-1)){
		 f <- f + 100*(x[i]^2 - x[i+1])^2 + (x[i]-1)^2;
	}
	return(f)
}
