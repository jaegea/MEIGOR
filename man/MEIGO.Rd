\name{MEIGO}
\alias{MEIGO}

\title{
	MEIGO main function
}

\description{
	Wrapper around the different optimisation methods
}

\usage{
	MEIGO(problem, opts, algorithm, ...)
}

\arguments{
  \item{problem}{
			List containing problem settings.
}
 \item{opts}{
	A list of n_threads lists containing options for each cooperative instance of essR.
}
 \item{algorithm}{ One of 	VNS, ESS, MULTISTART, CESSR, CEVNSR. 
  Check the documentation of each algorithm for more information. 
}
\item{...}{
	Additional input arguments.
}
}

\seealso{
	\code{\link{essR}}
	\code{\link{rvnds_hamming}}
	\code{\link{CeVNSR}}
	\code{\link{CeSSR}}
}

\examples{
#global optimum

#x*=[0.0898, -0.7127];
# or    
#x*=[-0.0898, 0.7127];
#
#f(x*)= -1.03163;

ex1 <- function(x){
	y<-4*x[1]*x[1]-2.1*x[1]^4+1/3*x[1]^6+x[1]*x[2]-4*x[2]*x[2]+4*x[2]^4;
	return(y)
}

#========================= PROBLEM SPECIFICATIONS ===========================
problem<-list(f=ex1,x_L=rep(-1,2),x_U=rep(1,2))
opts<-list(maxeval=500, ndiverse=40, local_solver='DHC', local_finish='LBFGSB', local_iterprint=1)
#========================= END OF PROBLEM SPECIFICATIONS =====================

Results<-MEIGO(problem,opts,algorithm="ESS");

}
