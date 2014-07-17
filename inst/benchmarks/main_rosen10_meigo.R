library(MEIGOR)
source(system.file("/benchmarks/rosen10.R",package="MEIGOR"));

nvar<-10;

problem<-list(f="rosen10", x_L=rep(-5,nvar), x_U=rep(1,nvar))

opts<-list(maxeval=2000, maxtime=3600*69, use_local=1, aggr=0, local_search_type=1, decomp=1, maxdist=0.5)

algorithm<-"VNS";

Results<-MEIGO(problem,opts,algorithm);