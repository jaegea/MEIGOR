\name{CeVNSR}
\alias{CeVNSR}

\title{
	Global optimization algorithm for MINLPs based on VNS using a Cooperative Strategy
}
\description{
	Solves optimization problems with intenger variables. Using several cooperative instances of VNS.
}
\usage{
	CeVNSR(	problem, opts, max_eval = Inf, max_time = Inf,
			n_iter = 1, is_parallel = TRUE, type = "SOCKS",
			global_save_list = NULL, ...)
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
	
	problem[[ith_thread]]=VNS_problem;
	opts[[ith_thread]]=VNS_opts;
	
	VNS_problem and VNS_opts correspond to lists as seen in the \link{rvnds_hamming} documentation.  
	
}
\value{
\item{f_mean}{Vector with size of n_iter+1 containing the mean value of the objective function in each iteration.}
 \item{f_sd}{Vector with size of n_iter+1 containing the standard deviation value of the objective function in each iteration.}
 \item{fbest}{Vector with size of n_iter+1 containing the best value of the objective function in each iteration.}
 \item{iteration_res}{A list containing the results from every VNS instance initialized. It follows the format: results$iteration_res[[iteration+1]][[thread_number]].}
 \item{numeval}{Vector with size of n_iter+1 containing the number objective function evaluations at the end of each iteration.}
 \item{time}{Vector with size of n_iter+1 containing the time spent at the end of an iteration.}
 \item{x_sd}{A list containing the standard deviation of decision each variable at the end of an iteration. It follows the format: results$iteration_res[[iteration+1]][[thread_number]]}
 \item{xbest}{A list containing the best set of decision variables found and the end of each iteration.}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\link{rvnds_hamming}
	\link{MEIGO}
}
%\examples{
	%
	%rosen10<-function(x){
	%	f<-0;
	%	n=length(x);
	%	for (i in 1:(n-1)){
	%		f <- f + 100*(x[i]^2 - x[i+1])^2 + (x[i]-1)^2;
	%	}
	%	return(f)
	%}
%
%	nvar=20;
%
%	sfStop()
%	
%	problem<-list(f=rosen10, x_L=rep(-1000,nvar), x_U=rep(1000,nvar))
%
%	opts=list();
%	opts[[1]]=list(use_local=1,aggr=1,local_search=1,decomp=1,maxdist=0.8,maxeval=2000);
%	opts[[2]]=list(use_local=1,aggr=0,local_search=2,decomp=0,maxdist=0.5,maxeval=2000);
%	opts[[3]]=list(use_local=1,aggr=0,local_search=2,decomp=0,maxdist=0.5,maxeval=2000);
%	opts[[4]]=list(use_local=1,aggr=0,local_search=2,decomp=0,maxdist=0.5,maxeval=2000);
%
%	opts$hosts=c('localhost','localhost','localhost','localhost');
%
%	opts$ce_niter=2;
%	opts$ce_type="SOCKS";
%	opts$ce_isparallel= TRUE;
%
%	Results=MEIGO(problem,opts, algorithm="CeVNSR");
%	sfStop()
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cooperative}
\keyword{strategies}
\keyword{variable}
\keyword{Neighbourhood}
\keyword{search}

