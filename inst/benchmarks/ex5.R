ex5<-function(x,texp,yexp){
	yini<-c(100, 0, 0, 0 ,0);
	times<-texp;
	out   <- lsodes(yini, times, ex5_dynamics, parms = x)
	
	tout<-out[,1];
	yout<-out[,-1];
	J <- sum((yout-yexp)^2);
	g<-0;
	residuals<-(yout-yexp);

	return(list(J,g,residuals))
}
#***************************************************
#Function of the dynamic system
ex5_dynamics<-function(t,y,p){
  dy<-rep(0,5)
  
  dy[1]<--(p[1]+p[2])*y[1];
  dy[2]<-p[1]*y[1];
  dy[3]<-p[2]*y[1]-(p[3]+p[4])*y[3]+p[5]*y[5];
  dy[4]<-p[3]*y[3];
  dy[5]<-p[4]*y[3]-p[5]*y[5];
  return(list(dy))
}
#***************************************************
  