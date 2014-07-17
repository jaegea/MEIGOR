ex5_dynamics<-function(t,y,p){
	dy<-rep(0,5)
	
	dy[1]<--(p[1]+p[2])*y[1];
	dy[2]<-p[1]*y[1];
	dy[3]<-p[2]*y[1]-(p[3]+p[4])*y[3]+p[5]*y[5];
	dy[4]<-p[3]*y[3];
	dy[5]<-p[4]*y[3]-p[5]*y[5];
	return(list(dy))

}