\name{CeSSR}
\alias{CeSSR}

\title{
	Global optimization algorithm for MINLPs based on Scatter Search using a Cooperative Strategy
}

\description{
	  CeSSR attempts to solve problems of the form:

	  \eqn{min f(x,p_1,p_2,...,p_n)}

subject to:

	 \eqn{c_e=0}\cr
	 \eqn{c_L \le c(x) \le c_U}\cr
	\eqn{x_L \le x \le x_U}\cr


%	\deqn{ min F(x)  subject to:
%	x\
%		ceq(x) = 0 (equality constraints)\cr
%		c_L <= c(x) <= c_U (inequality constraints)\cr
%		x_L <= x <= x_U (bounds on the decision variables)
}

\usage{
CeSSR(problem, opts, max_eval = Inf, max_time = Inf,
	 n_iter, is_parallel = TRUE, type = "SOCKS", global_save_list = NULL, ...)
}

\arguments{
  \item{problem}{
			List containing problem settings.
}
  \item{opts}{
	A list of n_threads lists containing options for each cooperative instance of essR.
}
  \item{max_eval}{
		Maximum bumber of evaluations. Default is Inf.
}
  \item{max_time}{
	Maximum time, default is Inf.
}
  \item{n_iter}{
	Number of cooperative iterations. Default is 0 which is the same as running multiple single thread (as many as  n_cpus)
	optimization runs.
}
  \item{is_parallel}{
		Default is TRUE. Sometimes this it is useful to use as FALSE for debugging.
}
%  \item{n_cpus}{
%	Number of cpus/cooperative threads per iteration.
%}
  \item{type}{
	Choose between "SOCKS" and "MPI". Default is "SOCKS" (socket-connection). If you are using "SOCKS" option and you
	want to run multiple cpus in different machines you must specify the adress of each machine in hosts.
	
	"MPI" mode requires you to have Rmpi installed.
}

  \item{global_save_list}{
	Specify the names of global variables to be exported.
}
  \item{\dots}{
	Additional variables.
}
}

\details{
	Check essR documentation for more information about the input arguments.
}

\value{
 \item{f_mean}{Vector with size of n_iter+1 containing the mean value of the objective function in each iteration.}
 \item{f_sd}{Vector with size of n_iter+1 containing the standard deviation value of the objective function in each iteration.}
 \item{fbest}{Vector with size of n_iter+1 containing the best value of the objective function in each iteration.}
 \item{iteration_res}{A list containing the results from every eSSR instance initialized. It follows the format: results$iteration_res[[iteration+1]][[thread_number]]. See also \code{\link{essR}}}
 \item{numeval}{Vector with size of n_iter+1 containing the number objective function evaluations at the end of each iteration.}
 \item{time}{Vector with size of n_iter+1 containing the time spent at the end of an iteration.}
 \item{x_sd}{A list containing the standard deviation of decision each variable at the end of an iteration. It follows the format: results$iteration_res[[iteration+1]][[thread_number]]}
 \item{xbest}{A list containing the best set of decision variables found and the end of each iteration.}
}

\seealso{
	\code{\link{essR}}
}

\examples{

rosen10<-function(x){
	f<-0;
	n=length(x);
	for (i in 1:(n-1)){
		 f <- f + 100*(x[i]^2 - x[i+1])^2 + (x[i]-1)^2;
	}
	return(f)
}

nvar=20;
problem<-list(f=rosen10, x_L=rep(-1000,nvar), x_U=rep(1000,nvar));

#Set 1 nodes and 2 cpu's per node
n_nodes=1;
n_cpus_per_node=3;

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
D=10;

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
    
    opts[[counter]]$maxeval			=	10000;
    opts[[counter]]$local_solver	=	"dhc";
    
    #Options not set will take default values for every thread
    
  }
}

#Set the address of each machine, defined inside the 'for' loop
opts$hosts=c('localhost','localhost','localhost');

#Do not define the additional options for cooperative methods (e.g., ce_maxtime, ce_isparallel, etc..)
#They will take their default values
opts$ce_niter=2;
opts$ce_type="SOCKS";
opts$ce_isparallel=TRUE;

#Call the solver
Results<-MEIGO(problem,opts,algorithm="CeSSR")

}

\keyword{cooperative}
\keyword{strategies}
\keyword{scatter}
\keyword{search}
\keyword{metaheuristic}
