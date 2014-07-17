library(MEIGOR)
library(deSolve);
source(system.file("/benchmarks/ex5.R",package="MEIGOR"));

#========================= PROBLEM SPECIFICATIONS ===========================
problem<-list(f=ex5, x_L=rep(0,5), x_U=rep(1,5), x_0=rep(0.5,5))
opts<-list(maxeval=3000, log_var=1:5, local_solver="bierdi")
#========================= END OF PROBLEM SPECIFICATIONS =====================
#time intervals

texp<-c(0, 1230, 3060, 4920, 7800, 10680, 15030, 22620, 36420)
 
#Distribution of species concentration
#	     y(1)    y(2)    y(3)    y(4)    y(5)
 
 yexp<-rbind(c(100.0, 0.0, 0.0, 0.0, 0.0),
 c(88.35, 7.3, 2.3, 0.4, 1.75),
 c( 76.4  , 15.6,    4.5 ,   0.7 ,   2.8),
 c(65.1 ,   23.1  ,   5.3 ,    1.1 ,    5.8),
 c(50.4  ,  32.9   ,  6.0   ,  1.5   ,  9.3),
 c(37.5  ,  42.7 ,    6.0 ,    1.9  ,  12.0),
 c(25.9  ,  49.1  ,   5.9  ,   2.2 ,   17.0),
 c( 14.0  ,  57.4    , 5.1  ,   2.6   , 21.0),
 c(4.5 ,   63.1  ,   3.8   ,  2.9  ,  25.7));
  

x0<-rep(0.001,5)
res<-solnp(x0, fun=ex5, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL, ineqUB = NULL, LB=rep(0,5), UB=rep(1,5),  control = list(),texp,yexp)