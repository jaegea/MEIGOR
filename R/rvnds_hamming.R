rvnds_hamming<-function(problem,opts,...){

  if (!is.numeric(opts$maxeval) && !is.numeric(opts$maxtime)){
    cat("WARNING:Either opts.maxeval or opts.maxtime must be defined as a stop criterion \n")
    cat("Define at least one of these two options and rerun \n")
    Results <- numeric(0);
    return(Results)
    stop()
  } else{
    if (!is.numeric(opts$maxeval)){
      maxeval <- 1e12;
    } else{
      maxeval <- opts$maxeval;
    }
    if (!is.numeric(opts$maxtime)){
      maxtime <- 1e12;
    } else{
      maxtime <- opts$maxtime;
    }
    if(!(maxtime>0)){
      Results <- numeric(0);
      return(Results);
    }
    if(!(maxeval>0)){
      Results <- numeric(0);
      return(Results);
    }
  }
  
  
#Initialize time
	cpu_time<-proc.time()[3];
	
  
  #Load default values for the options
  
  default <- vns_defaults();
  
  nargin <- length(as.list(match.call())) -1;
  
  if (nargin<2){opts <- numeric(0)};
  
  #Set all optiions
  opts<-vns_optset(default,opts);
    
  maxdist<-opts$maxdist;
  use_local<-opts$use_local;
  aggr<-opts$aggr;
  local_search_type<-opts$local_search_type;
  decomp<-opts$decomp;
  
  fobj<-problem$f;
  

	x_U<-problem$x_U;
	x_L<-problem$x_L;
  
	
  #Check if the bounds have the same dimension
  if (length(x_U)!= length(x_L)){
    cat("Upper and lower bounds have different dimensions!!! \n")
    cat("EXITING")
    Results<-numeric(0);
    return(Results)
    stop()
  } else{
    #Number of decision variables
    nvar<-length(x_L);
  }
  
  
  prob_names<-names(problem);
  
  temp<-match("x_0",prob_names);
  if (is.na(temp)){x0<-round(runif(nvar)*(x_U-x_L+1)+(x_L-0.5))} else{x0<-problem$x_0}

	xcurrent<-x0;
	xbest<-x0;
	
	fbest<-do.call(fobj,list(x0,...));
	
	nfuneval<-1;
	
	fcurrent<-fbest;
	
	cat("Init.point - Bestf:", fbest, "CPUTime:", proc.time()[3]-cpu_time,"\n")

	
	Results<-list(fbest=fbest, xbest=x0, func=fbest, neval=1, time=proc.time()[3]-cpu_time, numeval=1, x=xbest);
	
#Initial neighborhood is k=1
	kcurrent<-1
	
#Maximum neighborhood
	kmax<-floor(maxdist*nvar);
  
  #Perturb at least one dimension
  if (!kmax){kmax=1}
	
 	improve<-0;
	
	while (1){
		xnew<-xcurrent;
		shaked_vars<-sample.int(nvar,kcurrent);
		    
		for (i in 1:kcurrent){
			continuar<-0;
			while (!continuar){
				xnew[shaked_vars[i]]<-round(runif(1)*(x_U[shaked_vars[i]]-x_L[shaked_vars[i]]+1)+(x_L[shaked_vars[i]]-0.5));
				if (xnew[shaked_vars[i]]!=xcurrent[shaked_vars[i]]){
					continuar<-1;				
				}
			}
		}
		
		fnew<-do.call(fobj,list(xnew,...));
		nfuneval<-nfuneval+1;
		
		
		
		if (fnew<fbest){
			fbest<-fnew;
			xbest<-xnew;
			xcurrent<-xnew;
			fcurrent<-fnew;
			improve<-1;
			Results$func<-c(Results$func, fbest);
			Results$time<-c(Results$time, proc.time()[3]-cpu_time);
			Results$neval<-c(Results$neval, nfuneval);
            Results$x<-rbind(Results$x,xbest);
		}
		
		
		if (nfuneval>=maxeval || (proc.time()[3]-cpu_time)>=maxtime){
			cat("NFunEvals:", nfuneval, "Bestf:", fbest, "CPUTime:", proc.time()[3]-cpu_time,"\n")
			cat("*************************", "\n")
			cat("END OF THE OPTIMIZATION", "\n")
			cat("Best solution value", fbest, "\n");
            cat("Decision vector", xbest, "\n");
            cat("CPU time:",proc.time()[3]-cpu_time,"\n");
			cat("Number of function evaluations", nfuneval,"\n");
					
			Results$xbest<-xbest;
			Results$fbest<-fbest;
			Results$func<-c(Results$func, fbest);
			Results$x<-rbind(Results$x,xbest);
			Results$time<-c(Results$time, proc.time()[3]-cpu_time);
			Results$neval<-c(Results$neval, nfuneval);
			Results$numeval<-nfuneval;
			Results.cpu_time=proc.time()[3]-cpu_time;
#save Results in a file
			save(opts, problem, Results, file = "VNSR_report.RData")
			return(Results)
		}
		
		
#Start the local phase
		if (use_local){
			if (!aggr){
#In the non-aggressive scheme we apply local search over the new point even if it does no outperformed xbest
				f0<-fnew;
				x0<-xnew;
			} else if (aggr && improve){
				f0<-fcurrent;
				x0<-xcurrent;
			}
			
			if (!aggr || improve){			
#In the aggressive scheme, it does not make sense to use the local search if xnew did not improve xbest
				res<-rvnds_local(x0,f0,fobj,x_L,x_U,local_search_type,shaked_vars,decomp,nvar,maxeval,maxtime,cpu_time,fbest,xbest,nfuneval,Results,...);
				xout<-res[[1]];
				fout<-res[[2]];
				improve_local<-res[[3]];
				evals_local<-res[[4]];
                Results<-res[[5]];
				nfuneval<-nfuneval+evals_local;
				
				if (improve_local==2){
#This means that the local search improved the best value
					improve<-1;
					xbest<-xout;
					fbest<-fout;
				}
			}
		}
		if (improve){
			improve<-0;
			xcurrent<-xbest;
			fcurrent<-fbest;
			kcurrent<-1;
			if (!use_local || improve_local<2){
				cat("NFunEvals:", nfuneval, "Bestf:", fbest, "CPUTime:", proc.time()[3]-cpu_time,"\n")
			}
		} else{
			kcurrent<-kcurrent+1;
			if (kcurrent>kmax){
				kcurrent<-1;
			}
		}
	}
}




