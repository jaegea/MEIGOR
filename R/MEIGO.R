MEIGO<-function(problem,opts,algorithm,...){
  algorithm<-toupper(algorithm)

  if (algorithm=="VNS"){
    Results<-rvnds_hamming(problem,opts,...);
  }
  else if (algorithm=="ESS"){
    Results<-essR(problem,opts,...);
  }
  else if (algorithm=="MULTISTART"){
    Results<-essR_multistart(problem,opts,...);
  }
  else if (algorithm=="CESSR" || algorithm=="CEVNSR"){
    
    
    opts_names<-names(opts);
    temp<-match("ce_maxeval",opts_names);
    if (is.na(temp)){max_eval=Inf} else{max_eval=opts$ce_maxeval}
    
    temp<-match("ce_maxtime",opts_names);
    if (is.na(temp)){max_time=Inf} else{max_time=opts$ce_maxtime}
    
    temp<-match("ce_niter",opts_names);
    if (is.na(temp)){n_iter=1} else{n_iter=opts$ce_niter}
  
    temp<-match("ce_isparallel",opts_names);
    if (is.na(temp)){
		is_parallel=TRUE;
	}else{
		if(opts$ce_isparallel){
			is_parallel=TRUE
		}else{
			is_parallel=FALSE
		}
	}
    
    temp<-match("ce_type",opts_names);
    if (is.na(temp)){type="SOCKS"} else{type=opts$ce_type}
    
    temp<-match("global_save_list",opts_names);
    if (is.na(temp)){global_save_list=NULL} else{global_save_list=opts$global_save_list}

    
    if (algorithm=="CESSR"){
     Results<-CeSSR(problem,opts,max_eval,max_time, n_iter, is_parallel, type,global_save_list,...)
    }
    
    else if (algorithm=="CEVNSR"){
      Results<-CeVNSR(problem, opts, max_eval, max_time, n_iter, is_parallel, type,...)
    }
  }
  
  else {
  
    cat("The method defined in 'algorithm' is not valid \n")
    cat("Define a valid method (VNS or eSS) and rerun \n")
	
    Results <- numeric(0);}
}
