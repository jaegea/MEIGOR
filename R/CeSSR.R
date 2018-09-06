ith_essR<-function(i ,problem,opts_par,...){
	
	t0=Sys.time();
	
	res=essR(problem[[i]],opts_par[[i]]); 
	
	res$t0=t0;
	res$tf=Sys.time();
	
    return(res);   
}

CeSSR <-function(
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
		cat( "Running in parallel mode on", sfCpus(), "cpus.\n" )
		#Export global vars
		sfExport(list=global_save_list);
	
		res=sfClusterEval(.libPaths(.libPaths()));

		res=sfClusterEval(library(MEIGOR));
		
	}else{
		cat( "Running in sequential mode.\n" );
	}
	
	cpu_time_0=proc.time()[3];

	res=sfLapply(1:n_nodes,ith_essR,problem,opts);
	#If it is SGE a new session is opened every time so pass libpaths
	
	results$time=c();
	results$f_sd=c();
	results$fbest=c();
	results$xbest=list();
	results$f_mean=c();
	results$x_sd=list();
	results$numeval=c();
	results$iteration_res=list();
	results$iteration_res[[1]]=res;
	
	matF_0=c();
	matX_0=c();

	for(j in 1:n_nodes){
		#Put all xbest and fbest into a single matrix/vector
		matX_0=rbind(matX_0,res[[j]]$Refset$x);
		matF_0=c(matF_0,res[[j]]$Refset$f);
		
		total_evals=total_evals+res[[j]]$numeval;
		
		#Xbest is not necessarily in the reference set
		matX_0=rbind(matX_0,res[[j]]$xbest);
		matF_0=c(matF_0,res[[j]]$fbest);
		
		#update global fbest and xbest
		if(res[[j]]$fbest<=fbest){
			fbest=res[[j]]$fbest;
			xbest=res[[j]]$xbest;
		}
	}
	
	matX_0=rbind(matX_0,xbest);
	matF_0=c(matF_0,fbest);
	
	#Find duplicated solutions in refset and remove them
	index_duplicated=which(duplicated(matX_0));
	
	if(length(index_duplicated)>0){
		matX_0=matX_0[-index_duplicated,];
		matF_0=matF_0[-index_duplicated];
	}
	
	#Store the initial solutions for each problem
	for(j in 1:n_nodes){
		problem[[j]]$x_0=matX_0;
		problem[[j]]$f_0=matF_0;
	}
  
	#update cputime
	cpu_time_1=proc.time()[3];
	results$numeval[1]=total_evals;
	results$time[1]=cpu_time_1-cpu_time_0;
	
	#compute and print iteration stats
	results$f_sd[1]=sd(matF_0);
	results$fbest[1]=fbest;
	results$xbest[[1]]=xbest;
	results$f_mean[1]=mean(matF_0);
	results$iteration_res[[1]]=res;
	
	#Its should not be mandatory to print stuff but is OK for now
	print(paste("iteration:",0,"cputime:",results$time[1],'num_evals;',results$numeval[1],
		"best",results$fbest[1],'mean',results$f_mean[1],'sd_F',results$f_sd[1],sep="   "));
	results$x_sd[[1]]=apply(matX_0,2,sd);
	
	#neither so save results...We should add some inter_save option
	save(file="Results_iter_0",results);
	
	#If you only nee
	if(n_iter>0){
		for(i in 1:n_iter)
		{
			f_0=list();
			x_0=list();
			count=0;
			
			#Break if we pass max_eval
			if(total_evals>max_eval)return(results);
			
			problem$x_0=matX_0;
			problem$f_0=matF_0;
			
			#Call eSS in a node
			res=sfLapply(1:n_nodes,ith_essR,problem,opts,...);
			
			matF_0=c();
			matX_0=c();

			for(j in 1:n_nodes)
			{
				#concatenate all the refsets
				matX_0=rbind(matX_0,res[[j]]$Refset$x);
				matF_0=c(matF_0,res[[j]]$Refset$f);
				
				total_evals=total_evals+res[[j]]$numeval;
				
				#Xbest is not necessarily in the reference set
				matX_0=rbind(matX_0,res[[j]]$xbest);
				matF_0=c(matF_0,res[[j]]$fbest);
				
				#global xbest and fbest
				if(res[[j]]$fbest<=fbest){
					fbest=res[[j]]$fbest;
					xbest=res[[j]]$xbest;
				}
			}

			matX_0=rbind(matX_0,xbest);
			matF_0=c(matF_0,fbest);
			
			#Find and remove duplicated solutions
			index_duplicated=which(duplicated(matX_0));
			
			if(length(index_duplicated)>0){
				matX_0=matX_0[-index_duplicated,];
				matF_0=matF_0[-index_duplicated];
			}
			
			for(j in 1:n_nodes){
				problem[[j]]$x_0=matX_0;
				problem[[j]]$f_0=matF_0;
			}

			results$numeval[i+1]=total_evals;
			cpu_time_1=proc.time()[3];
			results$time[i+1]=cpu_time_1-cpu_time_0;
			results$f_sd[i+1]=sd(matF_0);
			results$fbest[i+1]=fbest;
			results$xbest[[i+1]]=xbest;
			results$f_mean[i+1]=mean(matF_0);
			results$x_sd[[i+1]]=apply(matX_0,2,sd);
			results$iteration_res[[i+1]]=res;
			
			print(paste("iteration:",i,"cputime:",results$time[i+1],'num_evals;',results$numeval[i+1],"best_f",results$fbest[i+1],'mean',results$f_mean[i+1],'sd_F',results$f_sd[i+1],sep="   "))
			
			#neither so save results...We should add some inter_save option
			save(file=paste("Results_iter_",i,sep=""),results);
			
		}                     
	}
	 
	 return(results);
}
