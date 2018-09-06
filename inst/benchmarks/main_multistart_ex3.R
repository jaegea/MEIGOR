#global optimum
#
#x*=[0.77152
# 	0.516994
# 	0.204189
# 	0.388811
# 	3.0355
# 	5.0973];
# 
#f(x*)= -0.388811;

library(MEIGOR)
source(system.file("/benchmarks/ex3.R",package="MEIGOR"));

#========================= PROBLEM SPECIFICATIONS ===========================
problem<-list(f="ex3",x_L=rep(0,6),x_U=c(rep(1,4),16,16), neq=4, c_L=-Inf, c_U=4)
opts<-list(ndiverse=25, local_solver='solnp', local_tol=3)
#========================= END OF PROBLEM SPECIFICATIONS =====================
k1=0.09755988;
k3=0.0391908;
k2=0.99*k1;
k4=0.9*k3;
Results<-MEIGO(problem,opts,algorithm="multistart",k1,k2,k3,k4);