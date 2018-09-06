# global optimum
#
#x*=[2.23607, 0, 1, 0];
#f(x*)=-40.9575;

library(MEIGOR);
source(system.file("/benchmarks/ex4.R",package="MEIGOR"));
#========================= PROBLEM SPECIFICATIONS ===========================
problem<-list(f="ex4", x_L=rep(0,4), x_U=rep(10,4), x_0=c(3,4,5,1),int_var=3, c_L=rep(-Inf,3), c_U=c(8,10,5))
opts<-list(maxtime=2)
#========================= END OF PROBLEM SPECIFICATIONS =====================
Results<-MEIGO(problem,opts,algorithm="ESS");