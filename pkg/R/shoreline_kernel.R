shoreline_kernel<-function(x){
  kernell<-matrix(rep(-1,x*x),nrow=x)
  if ( x %% 2 == 0){
    y<-x/2
    idx<-matrix(c(y,y,y+1,y+1,y,y+1,y,y+1),ncol=2)
    kernell[idx]<-(x*x-4)/4
  }else{
    kernell[ceiling(x/2),ceiling(x/2)]<-x*x-1
  }
  kernell[kernell==0]<-1
  return(kernell)
}
