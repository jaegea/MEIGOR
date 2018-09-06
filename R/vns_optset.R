vns_optset <-
function(default,opts){
	
	if (length(opts)){
		opts_names <- names(opts);
		default_names <- names(default);
		
		low_opts_names = tolower(opts_names);
		low_default_names = tolower(default_names);
			
		for (i in 1:length(opts_names)){
			j <- match(low_opts_names[i],low_default_names)
			if (is.na(j)){
				cat("Option '",low_opts_names[i],"'is not defined in VNS. It will be ignored \n")
			
			} else{
				default[[j]]<-opts[[i]];
			}
		}
	}
	opts=default;
	return(opts)
}

