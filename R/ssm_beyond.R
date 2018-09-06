ssm_beyond <-
function(z1,z2,z2_val,fobj,nrand,tolc,weight,x_L,x_U,c_L,c_U,nconst,int_var,bin_var,nfuneval,prob_bound,nvar,...){
	continuar <- 1;
	denom <- 1;
	n_improve <- 1;
	
	
	vector <- numeric(0);
	vector_value <- numeric(0);
	vector_value_penalty <- numeric(0);
	vector_penalty <- numeric(0);
	vector_nlc <- numeric(0);
	
	new_child <- numeric(0);
	new_child_value <- numeric(0);
	new_child_value_penalty <- numeric(0);
	new_child_nlc <- numeric(0);
	new_child_penalty <- numeric(0);
	
	while (continuar){
#Continue
		nvar<-length(x_L);
		zv<-matrix(0,2,nvar);
		d <- (z2-z1)/denom;
		zv[1,] <- z2;
		zv[2,] <- z2+d;
		
		aaa=which(zv[2,]<x_L);
		bbb=which(zv[2,]>x_U);
		
		if (length(aaa)>0){
			if (runif(1)>prob_bound){
				zv[2,aaa]=x_L[aaa];
			}
		}
        
		if (length(bbb)>0){
			if (runif(1)>prob_bound){
				zv[2,bbb]=x_U[bbb];
			}
		}
		
		xnew <- zv[1,]+(zv[2,]-zv[1,])*runif(nrand);   #non convex
		
#Evaluate
		output <- ssm_evalfc(xnew,x_L,x_U,fobj,nconst,c_L,c_U,tolc,weight,int_var,bin_var,nvar,...);
		val <- output[[1]];
		val_penalty <- output[[2]];
		pena <- output[[3]];
		nlc <- output[[4]];
		includ <- output[[5]];
		x <- output[[6]];
		
		nfuneval <- nfuneval+1;
		
		if (includ){
        new_child <- rbind(new_child,x);
        new_child_value <- c(new_child_value,val);
        new_child_value_penalty <- c(new_child_value_penalty,val_penalty);
        new_child_nlc <- rbind(new_child_nlc,nlc);
        new_child_penalty <- c(new_child_penalty,pena);
		
			if (val_penalty<z2_val){
				z1 <- z2;
				z2 <- xnew;
				z2_val <-val_penalty;
				
				vector <- x;
				vector_value <- val;
				vector_value_penalty <- val_penalty;
				vector_penalty <- pena;
				vector_nlc <- nlc;
				n_improve <- n_improve+1;
				
				if (n_improve==2){
					denom <- denom/2;
					n_improve <- 0;
				}
			} else{
#If it does not improve, break the loop
				break
			}
		} else{
			continuar <- 0;
		}
	}
	return(list(vector,vector_value,vector_penalty,vector_value_penalty,vector_nlc,new_child,new_child_value,new_child_penalty,new_child_value_penalty,new_child_nlc,nfuneval));
}

