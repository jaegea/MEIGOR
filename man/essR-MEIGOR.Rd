\name{essR-package}
\alias{essR-package}
\docType{package}

\title{
	Global optimization algorithm for MINLPs based on Scatter Search
}

\description{
	  essR attempts to solve problems of the form:
 
	\verb{
	min F(x)  subject to: 
	x
		ceq(x) = 0 (equality constraints)
		c_L <= c(x) <= c_U (inequality constraints)
		x_L <= x <= x_U (bounds on the decision variables)
	}
          Constraint functions, if applicable, must be declared in the same script as the objective function as a second output argument, e.g.:
	\verb{
 	 	 myfunction <- function(x){
				calculate fx - scalar containing the objective function value
				calculate gx - vector (or empty) containing the constraints values
				return(list(fx,gx))
		}
	}
}

\details{
	\tabular{ll}{
		Package: \tab MEIGOR\cr
		Type: \tab Package\cr
		Version: \tab 0.99.6\cr
		Date: \tab 2012-02-04\cr
		License: \tab GPL-3\cr
		LazyLoad: \tab yes\cr
	}
}

\author{

	Jose Egea
	David Henriques
	Thomas Cokelaer
	Alejandro F. Villaverde
	Julio R. Banga
	Julio Saez-Rodriguez

	Maintainer:Jose Egea <josea.egea@upct.es>

}
\references{
 Egea, J.A., Maria, R., Banga, J.R. (2010) An evolutionary method for
 complex-process optimization. Computers & Operations Research 37(2):
 315-324.
 
 Egea, J.A., Balsa-Canto, E., Garcia, M.S.G., Banga, J.R. (2009) Dynamic
 optimization of nonlinear processes with an enhanced scatter search
 method. Industrial & Engineering Chemistry Research 49(9): 4388-4401. 
}
\keyword{ package }
\keyword{scatter}
\keyword{search}% __ONLY ONE__ keyword per line
\keyword{optimization}

\seealso{
	\code{\link{essR}}
}
