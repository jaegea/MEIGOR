rvnds_local<-function(x0,f0,fobj,x_L,x_U,search_type,shaked_vars,decomp,nvar,maxeval,maxtime,cpu_time,fbest,xbest,nfuneval,Results,...){
  evals_local<-0;
  last_var_used<-0;
  
  fin<-0;
  
  xcurrloc<-x0;
  fcurrloc<-f0;
  
  if (decomp){
    u_vars<-shaked_vars;      
    n_u_vars<-length(shaked_vars);
  } else{
    u_vars<-1:nvar;
    n_u_vars<-nvar;
  }
  
  improve_solution<-0;
  
  while (!fin){
    x_improvements<-NULL;
    f_improvements<-NULL;
    vars_improving<-NULL;
    exit_local<-0;
    
    possible_directions<-kronecker(matrix(1,n_u_vars,1),t(c(-1,1)));
    
    #Find variables touching bounds
    aaa<-which(xcurrloc[u_vars]==x_L[u_vars]);
    bbb<-which(xcurrloc[u_vars]==x_U[u_vars]);
    
    
    #Adjust the possible directions. i.e., variables touching bounds can only go in one direction
    possible_directions[aaa,1]<-0;
    possible_directions[bbb,2]<-0;
    
    #Choose randomly the order of the variables to explore
    ccc<-sample.int(n_u_vars);
    possible_directions<-rbind(NULL,possible_directions[ccc,]);
    #Do the trick above to avoid transformation of a matrix into a vector
    
    local_improvement<-0;
    
    for (i in 1:n_u_vars){
      if (u_vars[ccc[i]]!=last_var_used){
        for (j in 1:2){
          #2 possible directions: +1 and -1 in the best case
          exit_local<-0
          if (possible_directions[i,j]){
            #If we can move to that direction
            xtemp<-xcurrloc;
            xtemp[u_vars[ccc[i]]]<-xcurrloc[u_vars[ccc[i]]]+possible_directions[i,j];
            ftemp<-do.call(fobj,list(xtemp,...));
            evals_local<-evals_local+1;
            
            if (ftemp<fcurrloc){
              local_improvement<-1;
              improve_solution<-max(improve_solution,1);
              if (ftemp<fbest){
                fbest<-ftemp;
                xbest<-xtemp;
                improve_solution<-2;
                    
                Results$xbest<-xbest;
                Results$fbest<-fbest;
                Results$func<-c(Results$func, fbest);
                Results$x<-rbind(Results$x,xbest);
                Results$time<-c(Results$time, proc.time()[3]-cpu_time);
                Results$neval<-c(Results$neval, nfuneval+evals_local);
                Results$numeval<-nfuneval+evals_local;
                Results$cpu_time=proc.time()[3]-cpu_time;
                
                cat("NFunEvals:", nfuneval+evals_local, "Bestf:", fbest, "CPUTime:", proc.time()[3]-cpu_time,"\n")
                
              }
              
              
              #go beyond
              parar<-0;
              while (!parar){
                if (xtemp[u_vars[ccc[i]]]==x_L[u_vars[ccc[i]]] || xtemp[u_vars[ccc[i]]]==x_U[u_vars[ccc[i]]]){
                  #touching bounds!!!
                  parar<-1;
                } else{
                  #Continue searching in the same direction
                  xtemp2<-xtemp;
                  xtemp2[u_vars[ccc[i]]]<-xtemp[u_vars[ccc[i]]]+possible_directions[i,j];
                  ftemp2<-do.call(fobj,list(xtemp2,...));
                  evals_local<-evals_local+1;
                  if (ftemp2<fbest){
                    fbest<-ftemp2;
                    xbest<-xtemp2;
                    improve_solution<-2;
                    
                    Results$xbest<-xbest;
                    Results$fbest<-fbest;
                    Results$func<-c(Results$func, fbest);
                    Results$x<-rbind(Results$x,xbest);
                    Results$time<-c(Results$time, proc.time()[3]-cpu_time);
                    Results$neval<-c(Results$neval, nfuneval+evals_local);
                    Results$numeval<-nfuneval+evals_local;
                    Results.cpu_time=proc.time()[3]-cpu_time;
                    
                    cat("NFunEvals:", nfuneval+evals_local, "Bestf:", fbest, "CPUTime:", proc.time()[3]-cpu_time,"\n")
                  }
                  if (ftemp2<ftemp){
                    xtemp<-xtemp2;
                    ftemp<-ftemp2;
                  } else{
                    parar<-1;
                  }
                }
              }
              
              #end of go_beyond
              if (search_type==1){
                #First improvement
                xcurrloc<-xtemp;
                fcurrloc<-ftemp;
                last_var_used=u_vars[ccc[i]];
                exit_local<-1;
                break 
                #VER SI ESTE BREAK SALE DE LOS DOS FOR
              } else{
                #Best improvement
                x_improvements<-rbind(x_improvements,xtemp);
                f_improvements<-c(f_improvements,ftemp);
                vars_improving<-rbind(vars_improving,u_vars[ccc[i]]);
              }
            }
            if ((nfuneval+evals_local)>=maxeval || (proc.time()[3]-cpu_time)>=maxtime){
              improve<-improve_solution;
              if (search_type==2 && local_improvement){
                best_child<-min(f_improvements);
                iii<-which.min(f_improvements);
                xout<-x_improvements[iii,];
                fout<-best_child;
              } else{
                xout<-xcurrloc;
                fout<-fcurrloc;
              }
              return(list(xout,fout,improve,evals_local, Results))
            }
          }
        }
        if (exit_local){break}
      }
      
    }
    #Get out if the local search did not improve the current solution
    if (!local_improvement){fin<-1}       
    
    if (local_improvement && search_type==2){
      #Once we got out from the 2nd for loop, we check if we explore all the possible u_vars
      best_child<-min(f_improvements);
      iii<-which.min(f_improvements);
      xcurrloc<-x_improvements[iii,];
      fcurrloc<-best_child;
      last_var_used<-vars_improving[iii];
    }
  }
  xout<-xcurrloc;
  fout<-fcurrloc;
  
  improve<-improve_solution;
  
  return(list(xout,fout,improve,evals_local, Results))
}