library(MEIGOR)

source(system.file("/benchmarks/f20.R",package="MEIGOR"));

problem=list();
problem$f<-f20;
problem$x_L = rep(-100,1000);
problem$x_U = rep(100,1000);

#Set 1 nodes and 2 cpu's per node
n_nodes=1;
n_cpus_per_node=2;
max_time_per_iter=10;

#Set different values for dim_refset, bal and n2 for each of the 10 cpu's to be used
dim1 = 23;     bal1 = 0;     n2_1 = 0;
dim2 = 33;     bal2 = 0;     n2_2 = 0;
dim3 = 46;     bal3 = 0;     n2_3 = 2;
dim4 = 56;     bal4 = 0;     n2_4 = 4;
dim5 = 72;     bal5 = 0.25;  n2_5 = 7;
dim6 = 72;     bal6 = 0.25;  n2_6 = 10;
dim7 = 88;     bal7 = 0.25;  n2_7 = 15;
dim8 = 101;    bal8 = 0.5;   n2_8 = 20;
dim9 = 111;    bal9 = 0.25;  n2_9 = 50;
dim10 = 123;   bal10 = 0.25; n2_10 = 100;

opts_dim=c(dim1,dim2,dim3,dim4,dim5,dim6,dim7,dim8,dim9,dim10);
opts_bal=c(bal1,bal2,bal3,bal4,bal5,bal6,bal7,bal8,bal9,bal10);
opts_n2=c(n2_1,n2_2,n2_3,n2_4,n2_5,n2_6,n2_7,n2_8,n2_9,n2_10);
D=1000;

#Initialize counter and options
counter=0;
opts=list();
hosts=c();

for(i in 1:n_nodes){
  for(j in 1:n_cpus_per_node){
    
    counter=counter+1;
    
    #Set the name of every thread
    if(i<10)hosts=c(hosts,paste('node0',i,sep=""));
    if(i>=10 && i<100)hosts=c(hosts,paste('node',i,sep=""));	
    
    opts[[counter]]=list();
    
    #Set specific options for each thread
    opts[[counter]]$local_balance  	=	opts_bal[counter];
    opts[[counter]]$dim_refset     	= 	opts_dim[counter];
    opts[[counter]]$local_n2		=	opts_n2[counter];
    
    #Set common options for each thread
    opts[[counter]]$maxtime  		=	max_time_per_iter;
    opts[[counter]]$maxeval			=	Inf;
    opts[[counter]]$local_solver	=	'DHC';
    
    #Options not set will take default values for every thread
    
  }
}

#Set the address of each machine, defined inside the 'for' loop
opts$hosts=c('localhost','localhost');
#opts$hosts=hosts;

#Do not define the additional options for cooperative methods (e.g., ce_maxtime, ce_isparallel, etc..)
#They will take their default values
opts$ce_niter=4;
opts$ce_type="SOCKS";
opts$ce_isparallel=TRUE;

#Call the solver
Results<-MEIGO(problem,opts,algorithm="CeSSR")
