ssm_localsolver <-
  function(x0,x_L,x_U,c_L,c_U,neq,int_var,bin_var,fobj,
           local_solver,local_iterprint,local_tol,weight,nconst,tolc,...){  
    n_fun_eval<<-0;
    fobj_global<<-fobj;
    neq_global<<-neq;
    nconst_global<<-nconst;
    #trick to avoid passing lists to nls. Right now, this function does not admit extra parameters (...)
    extra_args<<-list(NULL,...)
    
    #http://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html
    
    
    if (match(local_solver,"NM", nomatch=0) |match(local_solver,"BFGS", nomatch=0) | match(local_solver,"CG", nomatch=0) |
      match(local_solver,"LBFGSB", nomatch=0) | match(local_solver,"SA", nomatch=0)){
      
      
      ndeps_o<-tolc*(x_U-x_L);
      
      meth<-switch(local_solver, "NM" = "Nelder-Mead",
                   "BFGS" = "BFGS",
                   "CG"="CG",
                   "LBFGSB" = "L-BFGS-B",
                   "SA" = "SANN")
      
      #new_list<-list(fobj,list(...));
      results<-optim(x0, optim_fobj, gr = NULL,
                     method = meth,
                     lower = x_L, upper = x_U,
                     control = list(ndeps=ndeps_o), hessian = FALSE)
      
      res<-list(results$par, results$value, n_fun_eval);
      
    } else if(match(local_solver,"SOLNP", nomatch=0)){
      
      if (neq){
        eqfun_solnp=solnp_eq;
        eqB_solnp=rep(0,neq);
        
      } else{
        eqfun_solnp=NULL;
        eqB_solnp=NULL;
        
      }
      
      if ((nconst-neq)>0){
        ineqfun_solnp=solnp_ineq;
        ineqLB_solnp=c_L[(neq+1):nconst];
        ineqUB_solnp=c_U[(neq+1):nconst];
      } else{
        ineqfun_solnp=NULL;
        ineqLB_solnp=NULL;
        ineqUB_solnp=NULL;
      }
      
      results<-solnp(x0, solnp_fobj, eqfun = eqfun_solnp, eqB = eqB_solnp, ineqfun = ineqfun_solnp, ineqLB = ineqLB_solnp, ineqUB = ineqUB_solnp, 
                     LB = x_L, UB = x_U, control = list(), ...)
      n_val=length(results$values);
      res<-list(results$pars, results$values[n_val], n_fun_eval);
      
    } else if(match(local_solver,"DHC", nomatch=0)){
      #dhc
      nvar<-length(x0);
      x0<-(x0-x_L)/(x_U-x_L);
      #initsize=max((x_U-x_L)/2);
      initsize<-0.1;
      
      
      if (local_tol==1){thres=1e-6;}
      else if (local_tol==2){thres=1e-8;}
      else if (local_tol==3){thres=1e-10;}
      
      
      results<-dhc(fobj,x0,initsize,thres,100*nvar,x_L,x_U,weight,c_L,c_U,local_iterprint,tolc,...);
      
      res<-list(results[[2]], results[[1]], results[[3]]);
      
      
      
    } else if(match(local_solver,"NL2SOL", nomatch=0)){
      
      
      
      
      if (local_iterprint){disp_iter=TRUE}else{disp_iter=FALSE}
      
      if (local_tol==1){tol=1e-4;}
      else if (local_tol==2){tol=1e-5;}
      else if (local_tol==3){tol=1e-6;}
      

      results<- nls( ~ nls_fobj(x), start = list(x=x0), trace=disp_iter, algorithm="port",lower=x_L, upper=x_U, control=list(warnOnly = TRUE, maxiter = 20, tol = 1e-05, minFactor = 1/1024), extra_args)
      
      parameters<-as.vector(coef(results));
      
      res<-list(parameters, NULL, n_fun_eval);
    }
    
    return(res)
  }

