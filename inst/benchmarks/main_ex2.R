# global optimum
#x*=[2.32952, 3.17849];    
#f(x*)=-5.50801

library(MEIGOR);
source(system.file("/benchmarks/ex2.R",package="MEIGOR"));

#========================= PROBLEM SPECIFICATIONS ===========================
problem<-list(f=ex2,x_L=rep(0,2),x_U=c(3,4), c_L=rep(-Inf,2), c_U=c(2,36))
opts<-list(maxeval=750, local_solver="DHC", local_n1=2, local_n2=3)
#========================= END OF PROBLEM SPECIFICATIONS =====================

Results<-MEIGO(problem,opts,algorithm="ESS");