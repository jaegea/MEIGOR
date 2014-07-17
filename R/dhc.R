dhc<-
  function(fobj,x,initsize,thres,budget,x_L,x_U,weight,c_L,c_U,iterprint,tolc,...){

    #n_out=nargout(fobj);
    
    if (length(c_L) || length(c_U)){n_out<-2} else{n_out<-1};
    
    NDIM<-length(x);
    THRESHOLD<-thres;
    INIT_SIZE<-initsize;
    
    numeval<-0;
    

    #Initialize vi
    v<-rep(0,NDIM);        
    u<-v;
    vi=-1;
    vvec=1;
    vr=-INIT_SIZE;
    xreal<-x*(x_U-x_L)+x_L;
    
    
    output <- do.call(fobj,list(xreal,...));
    fx <- output[[1]];
    if (n_out>1){
      cx <- output[[2]];
      penalty <- ssm_penalty_function(cx,c_L,c_U,tolc);
      fx<-fx+weight*penalty;
      }
    else{
      cx<-0;
    }

    numeval<-numeval+1;
    fxv<-1e30;
        
    #JRB
    nnn<-0;
    while(abs(vr)>=THRESHOLD){
      if (abs(vr)<2*THRESHOLD){
        maxiter<-2*NDIM;}
      else{maxiter=2;}
      iter<-0;

      while (fxv>=fx && iter<maxiter){
        if(iter==0){xv<-x;}
        else{xv[vi+1]<-xv[vi+1]-vr;}
        if (vvec){vvec=0;}
        
        vr<--vr;
        if (vr>0){vi<-(vi+1)%%NDIM;}
        xv[vi+1]<-xv[vi+1]+vr;      
        
        #Put solution inside bounds
        aaa<-which(xv<0);
        bbb<-which(xv>1);
      
        xv[aaa]<-0;
        xv[bbb]<-1;
        
        xvreal<-xv*(x_U-x_L)+x_L;
        
        
        output <- do.call(fobj,list(xvreal,...));
        fxv<-output[[1]];
        if (n_out>1){
          cxv<-output[[2]];
          penaltyv<-ssm_penalty_function(cxv,c_L,c_U,tolc);
          fxv<-fxv+weight*penaltyv;}
        else {cxv=0;}
        
        
        pen2<-0;
        aaa<-which(xvreal<x_L);
        bbb<-which(xvreal>x_U);
        
        if (length(aaa)){pen2<-pen2+sum((x_L[aaa]-xvreal[aaa]));}        
        if (length(bbb)){pen2<-pen2+sum((xvreal[bbb]-x_U[bbb]));}
      
        fxv<-fxv+weight*pen2;
        
        numeval<-numeval+1;
        iter<-iter+1;
        
        if (numeval>=budget){
          vr<-thres/10;
          break
          }
        }
      
      
      if (fxv>=fx | is.nan(fxv)){
        fxv<-1e30;
        vr<-vr/2;}
      else{
        fx<-fxv;
        x<-xv;
        if (iterprint){
          #nnn=nnn+1;
          #plot(nnn,fx,'o');hold on;drawnow;
          cat("NEvals:", numeval, "Bestf:",fx,"\n");
          }
        
        
        if (iter==0){
          if (vvec){
            u<-u+v;
            v<-v*2;
            xv<-xv+v;
            vr<-vr*2;}
          else{
            u[vi+1]<-u[vi+1]+vr;
            vr<-vr*2;
            xv[vi+1]<-xv[vi+1]+vr;
            }
          
          #Put inside bounds
          aaa<-which(xv<0);
          bbb<-which(xv>1);
          xv[aaa]<-0;
          xv[bbb]<-1;
          
          xvreal<-xv*(x_U-x_L)+x_L;
          
          
          output <- do.call(fobj,list(xvreal,...));
          fxv<-output[[1]];
          if (n_out>1){
            cxv<-output[[2]];
            penaltyv<-ssm_penalty_function(cxv,c_L,c_U,tolc);
            fxv<-fxv+weight*penaltyv;}
          else {cxv=0;}
          
          pen2<-0;
          aaa<-which(xvreal<x_L);
          bbb<-which(xvreal>x_U);
          
          if (length(aaa)){pen2<-pen2+sum((x_L[aaa]-xvreal[aaa]));}
          if (length(bbb)){pen2<-pen2+sum((xvreal[bbb]-x_U[bbb]));}
          
          fxv<-fxv+weight*pen2;
          numeval<-numeval+1;
          
          if (numeval>=budget){
            vr=thres/10;
            break
            }
          }
        else{
          xv<-xv+u;
          xv[vi+1]<-xv[vi+1]+vr;
          #Put inside bounds
          aaa<-which(xv<0);
          bbb<-which(xv>1);
          xv[aaa]<-0;
          xv[bbb]<-1;
          xvreal<-xv*(x_U-x_L)+x_L;
           
          #The next 22 lines are repeated several times within the code. It could be optimized
          output <- do.call(fobj,list(xvreal,...));
          fxv<-output[[1]];
          if (n_out>1){
            cxv<-output[[2]];
            penaltyv<-ssm_penalty_function(cxv,c_L,c_U,tolc);
            fxv<-fxv+weight*penaltyv;}
          else {cxv=0;}
          
          pen2<-0;
          aaa<-which(xvreal<x_L);
          bbb<-which(xvreal>x_U);
          
          if (length(aaa)){pen2<-pen2+sum((x_L[aaa]-xvreal[aaa]));}
          if (length(bbb)){pen2<-pen2+sum((xvreal[bbb]-x_U[bbb]));}
          
          fxv<-fxv+weight*pen2;
          numeval<-numeval+1;
          
          if (numeval>=budget){
            vr=thres/10;
            break
          }
            
          if (fxv>=fx | is.nan(fxv)){
            u<-rep(0,NDIM);
            xv<-x;
            u[vi+1]<-vr;
            vr<-vr*2;
            xv[vi+1]<-xv[vi+1]+vr;
            
            
            #Put inside bounds
            aaa<-which(xv<0);
            bbb<-which(xv>1);
            xv[aaa]<-0;
            xv[bbb]<-1;
            
            xvreal<-xv*(x_U-x_L)+x_L;
            
            output <- do.call(fobj,list(xvreal,...));
            fxv<-output[[1]];
            if (n_out>1){
              cxv<-output[[2]];
              penaltyv<-ssm_penalty_function(cxv,c_L,c_U,tolc);
              fxv<-fxv+weight*penaltyv;}
            else {cxv=0;}
            
            pen2<-0;
            aaa<-which(xvreal<x_L);
            bbb<-which(xvreal>x_U);
            
            if (length(aaa)){pen2<-pen2+sum((x_L[aaa]-xvreal[aaa]));}
            if (length(bbb)){pen2<-pen2+sum((xvreal[bbb]-x_U[bbb]));}
            
            fxv<-fxv+weight*pen2;
            numeval<-numeval+1;
            
            if (numeval>=budget){
              vr=thres/10;
              break
            }
          }
          
          else{
            x<-xv;
            fx<-fxv;
            u[vi+1]<-u[vi+1]+vr;
            v<-2*u;
            vvec<-1;
            xv<-xv+v;
            
            #Put inside bounds
            aaa<-which(xv<0);
            bbb<-which(xv>1);
            xv[aaa]<-0;
            xv[bbb]<-1;
            
            xvreal<-xv*(x_U-x_L)+x_L;
            
            output <- do.call(fobj,list(xvreal,...));
            fxv<-output[[1]];
            if (n_out>1){
              cxv<-output[[2]];
              penaltyv<-ssm_penalty_function(cxv,c_L,c_U,tolc);
              fxv<-fxv+weight*penaltyv;}
            else {cxv=0;}
            
            pen2<-0;
            aaa<-which(xvreal<x_L);
            bbb<-which(xvreal>x_U);
            
            if (length(aaa)){pen2<-pen2+sum((x_L[aaa]-xvreal[aaa]));}
            if (length(bbb)){pen2<-pen2+sum((xvreal[bbb]-x_U[bbb]));}
            
            fxv<-fxv+weight*pen2;
            numeval<-numeval+1;
            
            if (numeval>=budget){
              vr=thres/10;
              break
            }
            
            vr<-0;
            vr<-sum(v^2);
            vr<-sqrt(vr);
          }
        }
      }
    }
    
    if (fxv<fx & !is.nan(fxv)){
      fx<-fxv;
      x<-xv;
    }
    x<-x*(x_U-x_L)+x_L;
    
    return(list(fx,x,numeval));
  }












