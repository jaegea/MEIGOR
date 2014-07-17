\name{rvnds_hamming}
\alias{rvnds_hamming}

\title{
	Main VNS function
}
\description{
	VNS Kernel function
}
\usage{
rvnds_hamming(problem, opts, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{problem}{
	List containing problem settings definition.
}
  \item{opts}{
	List containing options (if set as opts <- numeric(0) default options will be loaded).  
}
  \item{\dots}{
	Additional variables passed to the objective function
}
}
\details{
		problem$f:	Name of the file containing the objective function (String).\cr
		problem$x_L: Lower bounds of decision variables (vector).\cr
		problem$x_U: Upper bounds of decision variables (vector).\cr
		problem$x_0: Initial point(s) (optional; vector or matrix).\cr
		problem$f_0: Function values of initial point(s) (optional). These values MUST correspond to feasible points.\cr
		
		User options:

		opts$maxeval: Maximum number of function evaluations  (Default 1000).\cr
		opts$maxtime: Maximum CPU time in seconds (Default 60).\cr
		opts$maxdist: Percentage of the problem dimension which will be perturbed in the furthest neighborhood (varies between 0 and1, default is 0.5).\cr
		opts$use_local: Uses local search (1) or not (0). The default is 1. \cr
		
		The following options only apply when the local search is activated:
		
		opts$use_aggr: Aggressive search. The local search is only applied when the best solution has been improved (1=aggressive search, 0=non-aggressive search, default:0).\cr
		opts$local search type: Applies a first (=1) or a best (=2) improvement scheme for the local search (Default: 1).\cr
		opts$decomp: Decompose the local search (=1) using only the variables perturbed in the global phase. Default: 1.\cr
}
\value{
	 \item{fbest}{Best objective function value found after the optimization}
	 \item{xbest}{Vector providing the best function value}
     \item{cpu_time}{Time in seconds consumed in the optimization}
     \item{func}{Vector containing the best objective function value after each iteration}
     \item{x}{Matrix containing the best vector after each iteration}
     \item{time}{Vector containing the cpu time consumed after each iteration}
     \item{neval}{Vector containing the number of function evaluations after each iteration}
     \item{numeval}{Number of function evaluations}
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

	nvar<-10;

	problem<-list(f="rosen10", x_L=rep(-5,nvar), x_U=rep(1,nvar))

	opts<-list(maxeval=2000, maxtime=3600*69, use_local=1, aggr=0, local_search_type=1, decomp=1, maxdist=0.5)

	algorithm<-"VNS";

	Results<-MEIGO(problem,opts,algorithm);
}


