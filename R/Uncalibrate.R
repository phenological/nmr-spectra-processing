#'  Uncalibrate
#'  
#'  @param da spectra dataElements (use for dire and pgpe mainly)
#'  @return da.Data after SR = 0
#'  @example 
#'  da<-local(get(load("~/Downloads/gemma_C2_SER_GMAr08_zggppezf.ivdr@PROF_PLASMA_PGPE256_3mm.daE")))
#'  spectra_new<-uncalibrate(da)



uncalibrate<-function(da){
  if (class(da)[1] != "dataElement"){
    cat(crayon::red("Input must be dataElement\n"))
    stop()
  }
  if (!grepl("PROF_",da@method)){
    cat(crayon::red("Input must be spectra data daE\n"))
    stop()
  }  
  if (length(which(grepl("procs",names(da@obsDescr))==TRUE))<1){
    cat(crayon::red("procs information is required \n"))
    stop()
  }
  
  sf<-da@obsDescr$procs$SF
  sr<-da@obsDescr$procs$SR
  shift<-round((sr/sf)/(as.numeric(da@varName[2])-as.numeric(da@varName[1])))
  for(i in 1: nrow(da@.Data)){
    s<-shift[i]
    if(s>0){
      idx<-seq(1,length(da@varName),1)-s
      idx<-idx[which(idx>0)]
      da@.Data[i,]<-c(rep(0,s),da@.Data[i,idx])
    }
    if(s<0){
      idx<-seq(-s,length(da@varName),1)
      da@.Data[i,]<-c(da@.Data[i,idx],rep(0,abs(s)-1))
    }
    if(s==0){
      da@.Data[i,]<-da@.Data[i,]
    }
  }
  colnames(da@.Data)<-as.numeric(da@varName)
  return(da@.Data)
}





