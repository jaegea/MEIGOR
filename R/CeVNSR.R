ith_VNSR<-function(i ,problem, opts, ...){
	
	t0=Sys.time();
	
	res<-rvnds_hamming(problem[[i]],opts[[i]],...);
   
	res$t0=t0;
	res$tf=Sys.time();
	
    return(res);   
}

CeVNSR <-function(
			problem,		opts,					max_eval=Inf,	
			max_time=Inf, 	n_iter=1, 				is_parallel= TRUE, 		
			type="SOCKS",	global_save_list=NULL,	...
){
	
	results=c();
	total_evals=0;
    fbest=Inf;
    xbest=NULL;
	n_nodes=length(opts$hosts);
	
	temp_problem=problem;
	problem=list();
	
	for(i in 1:n_nodes){
		problem[[i]]=temp_problem;
	}

	#initializations for snowfall, for now we can run SOCK only at our cluster
	#socketHosts corresponds to the names of the nodes or network addres.
	if(type=="SOCKS"){
		sfSetMaxCPUs(number=n_nodes);
		sfInit( parallel = is_parallel, socketHosts = opts$hosts);
	}else if(type=="MPI"){
		#sfSetMaxCPUs(number=n_nodes);
		sfInit( parallel = is_parallel,cpus=n_nodes, type="MPI");
	}else{
		sfInit();
	}
	
	#if we succeed to spawn slaves
	if(sfParallel()){
		#If not the program still runs sequentially... I hope :)
		cat( "Running in parallel mode on", sfCpus(), "nodes.\n" )
		
		#Export global vars
		sfExport(list=global_save_list);
		
		libpaths=.libPaths();
	
		res=sfClusterEval(.libPaths(.libPaths()));

		
	}else{
		cat( "Running in sequential mode.\n" );
	}

	
	#Try catch. THis makes debuging hard but prevents socket connection from staying open
    tryCatch({   
	
		libpaths=.libPaths();
		#If type is not SGE is snowfall. R sessions must be opened and the libraries loaded
		#We use our own library path as a library path for these sessions which are in the nodes
		#Install R libraries at /home/user/R...
		
		res=sfClusterEval(.libPaths(.libPaths()));
	
		res=sfClusterEval(getwd());
		
		sfClusterEval(library("MEIGOR"));
		
        cpu_time_0=proc.time()[3];
		
		#Call VNS in a node
		res=sfLapply(1:n_nodes,ith_VNSR,problem,opts,...);

        results$time=c();
        results$fbest=c();
        results$xbest=list();
        results$f_mean=c();
        results$numeval=c();
		results$iteration_res=list();
		results$iteration_res[[1]]=res;
        
        matF_0=c();
        matX_0=c();

        for(j in 1:n_nodes){
		
			#Put all xbest and fbest into a single matrix
            matX_0=rbind(matX_0,res[[j]]$xbest);
            matF_0=c(matF_0,res[[j]]$fbest);
			
            total_evals=total_evals+res[[j]]$numeval;
			#update fbest and xbest
            if(res[[j]]$fbest<=fbest){
                fbest=res[[j]]$fbest;
                xbest=res[[j]]$xbest;
            }
			
        }
		
		#Find duplicated solutions in refset and remove them
		index_duplicated=which(duplicated(matX_0));
		
		if(length(index_duplicated)>0){
		
			matX_0=matX_0[-index_duplicated,];
			matF_0=matF_0[-index_duplicated];
			
		}
		
		min_index=which(matF_0==min(matF_0));
	
		#Store the initial solutions for each problem
		for(j in 1:n_nodes){
		
			n_elements=dim(matX_0);
			if(!is.null(n_elements) && n_elements>1){
				problem[[j]]$x_0=matX_0[min_index,];
				problem[[j]]$f_0=matF_0[min_index];
			}else{
				problem[[j]]$x_0=matX_0;
				problem[[j]]$f_0=matF_0;
			}
			
		}
		
        #update cputime
        cpu_time_1=proc.time()[3];
        results$numeval[1]=total_evals;
        results$time[1]=cpu_time_1-cpu_time_0;
		
		#compute and print iteration stats
        results$fbest[1]=fbest;
        results$xbest[[1]]=xbest;
        results$f_mean[1]=mean(matF_0);
		results$iteration_res[[1]]=res;
		
		#Its should not be mandatory to print stuff but is OK for now
        print(paste("iteration:",0,"cputime:",results$time[1],'num_evals;',results$numeval[1],
			"best",results$fbest[1],'mean',results$f_mean[1],sep="   "));
		
		save(file="Results_iter_0",results);
        
		#If you only nee
		if(n_iter>0){
		
            for(i in 1:n_iter){
			
                f_0=list();
                x_0=list();
                count=0;
				
				#Break if we pass max_eval
				if(total_evals>max_eval)return(results);
		
				#Call VNS in a node
				res=sfLapply(1:n_nodes,ith_VNSR,problem,opts,...);
				
                matF_0=c();
                matX_0=c();

                for(j in 1:n_nodes){
					#concatenate all the refsets
                    matX_0=rbind(matX_0,res[[j]]$xbest);
                    matF_0=c(matF_0,res[[j]]$fbest);
					
					#global xbest and fbest
                    total_evals=total_evals+res[[j]]$numeval;
                    if(res[[j]]$fbest<=fbest){
						fbest=res[[j]]$fbest;
						xbest=res[[j]]$xbest;
                    }
                }
				
				index_duplicated=which(duplicated(matX_0));
				
				if(length(index_duplicated)>0){
					matX_0=matX_0[-index_duplicated,];
					matF_0=matF_0[-index_duplicated];
				}
				
				min_index=which(matF_0==min(matF_0));
	
				#Store the initial solutions for each problem
				for(j in 1:n_nodes){
				
					n_elements=dim(matX_0);
					if(!is.null(n_elements) && n_elements>1){
						problem[[j]]$x_0=matX_0[min_index,];
						problem[[j]]$f_0=matF_0[min_index];
					}else{
						problem[[j]]$x_0=matX_0;
						problem[[j]]$f_0=matF_0;
					}
					
				}
		
                results$numeval[i+1]=total_evals;
                cpu_time_1=proc.time()[3];
                results$time[i+1]=cpu_time_1-cpu_time_0;
                results$fbest[i+1]=fbest;
                results$xbest[[i+1]]=xbest;
                results$f_mean[i+1]=mean(matF_0);
				results$iteration_res[[i+1]]=res;
	
                print(paste("iteration:",i,"cputime:",results$time[i+1],'num_evals;',results$numeval[i+1],"best_f",results$fbest[i+1],'mean',
results$f_mean[i+1],sep="   "))
                
				#In case it craches
				save(file=paste("Results_iter_",i,sep=""),results);
				
            }                     
		}
	}
	#terrible identation ARGHHH....
	,error=function(e){print("something happened")},finally={
																
		sfStop();
	
	})
																 
	try({
		
			sfStop();
	});
	 
	 return(results);
}
