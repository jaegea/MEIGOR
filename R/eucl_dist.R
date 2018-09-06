eucl_dist <-
function(m1,m2){
	#euclidean distance between the rows of two different matrizes
	#the matrices must have the same number of columns
	
	n_a<-nrow(m1);
	n_b<-nrow(m2);
	
	ddd<-kronecker(matrix(1,1,n_b),apply(m1^2,1,sum))+kronecker(matrix(1,n_a,1),t(apply(m2^2,1,sum)))
	ddd<-ddd-2*m1%*%t(m2);
	ddd<-sqrt(ddd);
	
	#The result is a matrix whose element i,j is the euclidean distance between row i of matrix a and row j of matrix b
	
	return(ddd)
}

