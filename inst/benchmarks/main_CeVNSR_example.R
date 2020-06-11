library('MEIGOR');
source(system.file("/benchmarks/rosen10.R",package="MEIGOR"));

nvar=10;
problem<-list(f=rosen10, x_L=rep(-5,nvar), x_U=rep(1,nvar))

opts=list();
opts[[1]]=list();
opts[[2]]=list();

#Set options for each thread
opts[[1]]$use_local=1;
opts[[1]]$aggr=1;
opts[[1]]$local_search=1;
opts[[1]]$decomp=1;
opts[[1]]$maxdist=0.8;
opts[[1]]$maxeval=2000;

opts[[2]]$use_local=1;
opts[[2]]$aggr=0;
opts[[2]]$local_search=2;
opts[[2]]$decomp=0;
opts[[2]]$maxdist=0.5;
opts[[2]]$maxeval=2000;

opts$hosts=c('localhost','localhost');

opts$ce_niter=4;
opts$ce_type="SOCKS";
opts$ce_isparallel= TRUE;

Results=MEIGO(problem,opts, algorithm="CeVNSR")
