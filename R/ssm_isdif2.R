ssm_isdif2 <-
function(x,group,tol,flag){
	f<-0;
	ind<-numeric(0);
	ind2<-numeric(0);
	
	for (i in 1:nrow(group)){
		num<-abs(x-group[i,]);
		denom<-min(c(abs(x),abs(group[i,])));
		denom[abs(denom)<1e-16]<-1;
		diference<-num/denom;
		aaa<-which(diference>=tol);
		
		if (flag==1){
			if (!length(aaa)){
				f<-f+1;
				ind<-c(ind,i);
			} 
		} else if (flag==2){
			if (length(aaa)!=length(x)){
				ind<-c(ind,i);
			} else{
				ind2<-c(ind2,i);
			}
		}
	}
	return(list(f=f,ind=ind,ind2=ind2))
}

