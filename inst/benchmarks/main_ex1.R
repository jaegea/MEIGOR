#global optimum

#x*=[0.0898, -0.7127];
# or    
#x*=[-0.0898, 0.7127];
#
#f(x*)= -1.03163;

library(MEIGOR)
source(system.file("/benchmarks/ex1.R",package="MEIGOR"));

#========================= PROBLEM SPECIFICATIONS ===========================
problem<-list(f=ex1,x_L=rep(-1,2),x_U=rep(1,2))
opts<-list(maxeval=500, ndiverse=40, local_solver='DHC', local_finish='LBFGSB', local_iterprint=1)
#========================= END OF PROBLEM SPECIFICATIONS =====================

Results<-MEIGO(problem,opts,algorithm="ESS");
