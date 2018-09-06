ssm_penalty_function <-
function(nlc,c_L,c_U,tolc){
	P <- numeric(0);
	
	if (length(nlc)>0){
		a=which(nlc < c_L)
		b=which(nlc > c_U)
		
		P=rep(0,length(a)+length(b));
		counter <- 1;
		
		if (length(a)>0){
			for (i in 1:length(a)) {
				P[counter] <- (c_L[a[i]]-nlc[a[i]]);
#P=c(P, (c_L[a[i]]-nlc[a[i]]));#/(1+abs(c_L[a[i]])));
				counter=counter+1;
			}
		}
		
		
		if (length(b)>0){
			
			
			for (i in 1:length(b)) {
				P[counter] <- (nlc[b[i]]-c_U[b[i]]);
#P=c(P, (nlc[b[i]]-c_U[b[i]]));#/(1+abs(c_U[b[i]])));
				counter=counter+1;
			}
		}
	}
	
#maxP=max(P);
	
	
	if (length(P)>0 && max(P)>tolc){
#The penalty is the maximum constraint violation
		fx <- max(P);
	} else{
		fx <- 0;		
	}
	
	return(fx)	
}

