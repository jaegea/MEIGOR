essR_multistart <-
  function(problem,opts=list(maxeval=NULL,maxtime=NULL),...){
    # Function   : essR_multistart
    # Created on : 26/07/2012
    # Last Update: 26/07/2012
    # Email      : josea.egea@upct.es
    #
    # Script to do multistart optimization (e.g. Multiple local optimization
    # starting from different initial points uniformly distributed withing the bounds)
    # 
    #   essR_multistart uses the same syntax as essR and solves the same type of problems
    # 
    #   USAGE: Results_multistart <- essR_multistart(problem,opts,p1,p2,....,pn);
    # 
    # INPUT PARAMETERS:
    # ****************
    #   problem - List containing problem settings
    #       problem$f   = Name of the file containing the objective
    #                     function (String)
    #       problem$x_L = Lower bounds of decision variables (vector)
    #       problem$x_U = Upper bounds of decision variables (vector)
    #       problem$x_0 = Initial point(s) (optional; vector or matrix)
    #
    # 
    #       Additionally, fill the following fields if your problem has
    #       non-linear constraints
    #       problem$neq     = Number of equality constraints (Integer; do not define it
    #                         if there are no equality constraints)
    #       problem$c_L     = Lower bounds of nonlinear inequality constraints
    #                        (vector)
    #       problem$c_U     = Upper bounds of nonlinear inequality constraints
    #                         (vector)
    #       problem$int_var = Number of integer variables (Integer)
    #       problem$bin_var = Number of binary variables  (Integer)
    #
    # 
    # NOTE: The order of decision variables is x <- c(cont, int, bin)
    # 
    #   opts - List containing options (if set as opts <- numeric(0) default options
    #          will be loaded). See the script of default options to know their values
    #
    #           opts$ndiverse           =         Number of uniformly distributed points within
    #                                             the bounds to do the multistart optimization (default 100)
    #           opts$local_solver       =         Choose local solver
    #                                             0: Local search deactivated (Default),
    #                                             "NM", "BFGS", "CG", "LBFGSB",
    #                                             "SA","SOLNP", "DHC", "NLS2SOL"
    #           opts$iterprint            =       Print basic information on screen after 
    #                                             every local search (Binary, default=1);
    #           opts$local_tol                  = Level of tolerance in local
    #                                             search
    #           opts$local_iterprint            = Print each iteration of local
    #                                             solver on screen (Binary, default=0);
    # 
    # 
    #   p1,p2... :  optional input parameters to be passed to the objective
    #   function
    # 
    # 
    # OUTPUT PARAMETERS:
    # *****************
    # A file called "eSSR__multistart_report.Rdata" is generated containing.
    # 
    #   Results_multistart - Data frame containing results
    #   problem            - Data frame containing problem settings
    #   opts              - Data frame containing all options
    #
    # 
    # Fields in Results_multistart
    #       Results_multistart$fbest        = Best objective function value found
    #                                         after the multistart optimization
    #       Results_multistart$xbest        = Vector providing the best
    #                                         function value
    #       Results_multistart$x0           = Matrix containing the vectors used for
    #                                         the multistart optimization (in rows)
    #       Results_multistart$f0           = Vector containing the objective function
    #                                         values of the initial solutions used in 
    #                                         the multistart optimization
    #                                         iteration
    #       Results_multistart$xxx          = Matrix containing the solutions provided
    #                                         by the local optimization (in rows)
    #       Results_multistart$func         = Vector containing the objective function
    #                                         values obtained after every local search
    #       Results_multistart$no_conv      = Matrix containing the initial points that
    #                                         did not provide any solution when the local
    #                                         search was applied (in rows)
    #       Results_multistart$nfuneval     = Vector containing the number of function
    #                                         evaluations in every optimisation
    #       Results_multistart$time         = Total CPU time to carry out the multistart optimization
    #
    #
    # NOTE: To plot an histogram of the results: hist(Results_multistart$func)
    
    
    
    
    # Initialize time
    cpu_time=proc.time()[3];
    
    Results_multistart <- numeric(0);
    x_U<-problem$x_U;
    x_L<-problem$x_L;
    
    prob_names<-names(problem);
    
    temp<-match("neq",prob_names);
    if (is.na(temp)){neq<-0} else{neq<-problem$neq}
    temp<-match("c_U",prob_names);
    if (is.na(temp)){c_L<-numeric(0); c_U<-numeric(0)} else{c_U<-problem$c_U; c_L<-problem$c_L}
    temp<-match("int_var",prob_names);
    if (is.na(temp)){int_var<-0} else{int_var<-problem$int_var}
    temp<-match("bin_var",prob_names);
    if (is.na(temp)){bin_var<-0} else{bin_var<-problem$bin_var}
    
    #Load default values for the options
    default <- ssm_defaults();
    nargin <- length(as.list(match.call())) -1;
    if (nargin<2){opts <- numeric(0)};
    
    #Set all optiions
    opts<-ssm_optset(default,opts);
    
    weight=opts$weight;
    tolc=opts$tolc;
    local_solver=opts$local_solver;
    iterprint=opts$iterprint;
    
    
    #Local options
    if(is.character(local_solver)){local_solver<-toupper(local_solver)}
    local_tol=opts$local_tol;
    local_iterprint=opts$local_iterprint;
    
    
    
    #Check if the bounds have the same dimension
    if (length(x_U)!= length(x_L)){
      cat("Upper and lower bounds have different dimensions!!! \n")
      cat("EXITING")
      Results_multistart<-numeric(0)
      return(Results_multistart)
      stop()
    } else{
      #Number of decision variables
      nvar<-length(x_L);
    }
    
    
    fobj<-problem$f;
    

    
    
    if (int_var+bin_var){
      cat("No local solvers available to handle integer/binary variables at the moment \n");
      cat("EXITING \n");
      Results_multistart<-numeric(0);
      return(Results_multistart)
      stop()      
    } 
    
    
    
    #Check if the objective function has 3 output arguments when NL2SOL is invoked
    if(match(local_solver,"NL2SOL", nomatch=0)){
      n_out_f<-do.call(fobj,list(x_L,...));
      if (length(n_out_f)<3){
        cat("NL2SOL requires 3 output arguments: f, g, and the vector of residuals \n");
        cat("EXITING \n");
        Results_multistart<-numeric(0);
        return(Results_multistart)
        stop()
      }
    }
    
    #If there are equality constraints
    if (neq){
      #Set the new bounds for all the constraints
      c_L=c(-tolc*rep(1,neq), c_L);
      c_U=c(tolc*rep(1,neq), c_U);
    }
    
    nconst<-length(c_U);
    
    #Check if the objective function has 2 output arguments in constrained problems
    if (nconst){
      n_out_f<-do.call(fobj,list(x_L,...));
      if (length(n_out_f)<2){
        cat("For constrained problems the objective function must have at least 2 output arguments \n");
        cat("EXITING \n");
        Results_multistart<-numeric(0);
        return(Results_multistart)
        stop()
      }
    }
    
    func<-numeric(0);
    xxx<-numeric(0);
    no_conv<-numeric(0);
    x0<-numeric(0);
    nfuneval<-numeric(0);
    
    
    if (opts$ndiverse){
      multix<-matrix(runif(opts$ndiverse*nvar),opts$ndiverse,nvar);
      #Put the variables inside the bounds
      multix<-multix*kronecker(matrix(1,opts$ndiverse,1),t((x_U-x_L)))+kronecker(matrix(1,opts$ndiverse,1),t(x_L));
    }else{multix<-numeric(0)}
    
    
    multix<-rbind(x0, multix);
    
    if ((int_var) || (bin_var)){multix<-ssm_round_int(multix,int_var+bin_var,x_L,nvar)};
    
    ndiverse<-nrow(multix) 
    
    f0<-numeric(0)
    for (i in 1:ndiverse){
      temp<-do.call(fobj,list(multix[i,],...));
      f0[i]<-temp[[1]];
      x0<-multix[i,];
      
      if (opts$iterprint){
        cat("\n");
        cat("Local search number:", i, "\n");
        cat("Call local solver:", toupper(local_solver),  "\n")
        cat("Initial point function value:",f0[i] ,"\n");
        tic<-proc.time()[3];
      }
      
      res<-ssm_localsolver(x0,x_L,x_U,c_L,c_U,neq,int_var,bin_var,fobj,local_solver,local_iterprint,local_tol,weight,nconst,tolc,...);
      x<-res[[1]];
      value<-res[[2]];
      numeval<-res[[3]];
      
      feasible<-1;
      
      if (nconst){
        output <- do.call(fobj,list(x,...));
        value <- output[[1]];
        nlc <- output[[2]];
        
        if (any(nlc<c_L) || any(nlc>c_U)){
          feasible<-0
        }
      }       
      
      if (opts$iterprint){
        if (feasible){
          cat("Local solution function value:", value, "\n");
        }else{
          cat("UNFEASIBLE Local Solution! \n");  
        }
        
        cat("Number of function evaluations in the local search:", numeval, "\n");
        cat("CPU Time of the local search:",proc.time()[3]-tic,"seconds \n\n");
      }
      
      if (!feasible | is.infinite(value) | is.nan(value)){
        no_conv<-rbind(no_conv,x0);
      }else{
        func<-c(func, value);
        xxx<-rbind(xxx,x);
      }
      
      nfuneval<-c(nfuneval,numeval);
      
   
    }
    
    cpu_time<-proc.time()[3]-cpu_time;
    
    Results_multistart<-list(fbest=numeric(0), xxx=numeric(0), xbest=numeric(0), func=numeric(0));
    
    
    hist(func,xlab="Objective function value");
    iii<-which.min(func);
    Results_multistart$fbest<-func[iii];
    Results_multistart$xbest<-xxx[iii,];
    Results_multistart$x0<-multix;
    Results_multistart$f0<-f0;
    Results_multistart$func<-func;
    Results_multistart$xxx<-xxx;
    Results_multistart$no_conv<-no_conv;
    Results_multistart$nfuneval<-nfuneval;
    Results_multistart$time<-cpu_time;
    save(opts, problem, Results_multistart, file = "eSSR_multistart.RData") 
  }
