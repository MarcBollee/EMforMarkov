BSCIforMarkov<-function(C.all, S.all,#list of matrices for the one and two-step transition counts all of the same square dimension
                         M1, Counts,
                         Forbidden,          #matrix of forbidden transitions
                         BSsamples=5000,     #Number of samples for the CIs
                         alpha=0.05,
                         w,
                         precision=0.001, maxiteration=20, verbose=FALSE, plot=F, TrueM1,
                         returnVar=F, returnLCI=F, returnBS=F){
  #Set Up and Check S.all first
  S<-length(S.all)
  if(is.list(S.all)){unlist(S.all)}
  #Check all S.all are increasing integers >0
  if(S>1){
    if(sum(S.all!=as.integer(S.all))>0){stop("S.all must be a list of positive integers.")}
    for(t in 2:S){if(S.all[t-1]>=S.all[t]){stop("S.all must be strictly increasing.")}}
  }
  
  S.all=as.integer(S.all)
  if(!missing(M1)){n<-nrow(M1)}
  if(missing(M1)){
    if(missing(C.all)){stop("One of M1 or C.all must be provided.")}
    #Turn C.all into array for speed
    if(is.list(C.all)){
      n<-nrow(C.all[[1]])
      C<-array(0,c(length(C.all),n,n))
      for(t in 1:length(C.all)){C[t,,]=C.all[[t]]}
    }else{
      n<-nrow(C.all[1,,])
      C<-C.all
    }
    #make sure there is enough information for each row (at least one observation)
    TotalrowSums<-apply(C,c(1,3),sum)
    for(i in 1:n){
      if(TotalrowSums[i]==0){stop(paste("Not enough information about row",i))}
    }
    
    if(missing(Counts)){Counts<-apply(C,1,rowSums)}
    
    #Check C.all and S.all are of equal length
    if(S!=length(C[,1,1])){stop("C.all and S.all must be of equal length.")}
  }else{#M1 provided
    n<-nrow(M1)
    if(missing(Counts)){Counts<-matrix(100,n,S)}
  }
  
  #Make sure Forbidden is matrix of TRUE/FALSE values and define Forbidden if not provided
  if(missing(Forbidden)){
    Forbidden<-matrix(F,n,n)
  }else{
    if(nrow(Forbidden)!=n || ncol(Forbidden)!=n){stop("Forbidden matrix must be of the same size as the transition count matrices.")}
    for(i in 1:n){
      for(j in 1:n){
        if(Forbidden[i,j]!=0 && Forbidden[i,j]!=1){
          stop(paste("Forbidden must be a matrix of binary TRUE/FALSE indicators. Forbidden[", i, ",", j, "] is not binary."))
        }
      }
    }
  }
  Forbidden=(Forbidden==1)
  if(missing(w)){w<-rep(1/n,n)}
  ################################################################## Helpful functions
  #Define a function that normalises but does not break on a string of 0's
  carefulNormalise<-function(x){ 
    if(sum(x)==0){return(rep(1/length(x),length(x))) #take into account that vector can be 0, in which case the probability is assumed uniform accross states
    }else{return(x/sum(x))}
  }
  #Define a function that row-normalises a matrix using this careful normalise function defined above
  rowNormalise<-function(x){
    return(t(apply(x,1,carefulNormalise)))
  }
  ######################################################## Bootstrapped CIs
  #Get estimated matrix
  if(missing(M1)){
    M1<-GEMforMarkov7(C.all=C.all,S.all=S.all,Forbidden=Forbidden,precision=precision,maxiteration=maxiteration,verbose=F,w=w)
  }
  
  #Generate a bunch of 1&2-step transition count matrices with as many counts as were originally there in the data matrix
  M1.BS<-array(0,dim=c(BSsamples,n,n)) #Stores bootstrapped transition probability matrices obtained from bootstrapped count matrices
  C.BS<-array(0,c(n,n,S))
  M1.naive<-matrix(1/n,n,n)
  if(verbose==TRUE){message("Starting Bootstrap:")}
  for(b in 1:BSsamples){
    if(verbose==TRUE && b%%round(BSsamples/10)==0){message(paste(100*round(b/BSsamples,digits=2),"%"))}
    
    #Generate a dataset from true transition probability matrix M1
    for(t in 1:S){
      C.BS[,,t]=rMarkovCount(M1,S.all[t],Counts[,t])
    }
    
    #Calculate a naive estimate of M1 as a "good" and computationally cheap starting point
    M1.naive=matrix(1/n,n,n)
    
    if(S.all[1]==1){
      M1.naive=rowNormalise(C.BS[,,1])
    }else{
      for(t in 1:S){
        eigen<-eigen(rowNormalise(C.BS[,,t]))
        if((sum(is.complex(eigen$values))==0) && (sum(eigen$values<0)==0 || (S.all[t]%%2==0)) && (det(eigen$vectors)>=1e-10)){
          M1.naive=eigen$vectors%*%diag((eigen$values)^(1/S.all[t]))%*%solve(eigen$vectors) #Nth root estimate
          #Remove problematic values
          M1.naive[is.na(M1.naive)]=0.1/n
          M1.naive[is.nan(M1.naive)]=0.1/n
          M1.naive[which("&"(M1.naive<0,M1.naive>1))]=0.1/n
          #Row normalise before passing on the EM algorithm as the starting point
          M1.naive=rowNormalise(M1.naive)
          #Once a suitable candidate has been found from the smallest-step matrix, end the loop
          break
        }
      }
    }
    
    #Run the EM algorithm on each of the generated two-step transition count matrices
    M1.BS[b,,]=round(CoreEM(t(apply(C.BS,1,c)),Sall=S.all,n=n,w=w,M1naive=M1.naive,Forbidden=Forbidden,precision=precision,maxiteration=maxiteration),-log10(precision))
  }
  
  #Obtain a value for the alpha and 1-alpha quantiles
  M1.lCI<-matrix(0,n,n)
  M1.uCI<-matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      M1.lCI[i,j]=quantile(M1.BS[,i,j],alpha/2,na.rm=T)
      M1.uCI[i,j]=quantile(M1.BS[,i,j],1-(alpha/2),na.rm=T)
    }
  }
  
  #Plot
  if(plot==T){
    for(i in 1:n){
      for(j in 1:n){
        #print(paste(min(M1.BS[,i,j]),max(M1.BS[,i,j])))
        hist(M1.BS[which("&"(M1.BS[,i,j]>=0,M1.BS[,i,j]<=1)),i,j],
             freq=F, breaks=seq(from=-0.01,to=1.01,by=0.01), 
             main=paste("Estimate of M1[",i,",",j,"]   Using", paste(S.all),"step"), 
             xlim=c(0,1), xlab = "Transition probability")
        abline(v=mean(M1.BS[,i,j]),col="blue",lwd=2)
        abline(v=M1.lCI[i,j],col="red")
        abline(v=M1.uCI[i,j],col="red")
        abline(v=M1[i,j],col="black",lwd=2)
        if(!missing(TrueM1)){abline(v=TrueM1[i,j],col="gold",lwd=2)}
      }
    }
  }
  
  #Obtain Sum of Variances of elements
  if(returnVar==T){
    Var<-0
    for(i in 1:n){
      for(j in 1:n){
        Var=Var+var(M1.BS[,i,j])
      }
    }
  }
  if(returnLCI==T){
    LCI<-sum(M1.uCI-M1.lCI)
  }
  
  #Round to the level of precision specified
  M1.lCI=round(M1.lCI,ceiling(-log10(precision)))
  M1.uCI=round(M1.uCI,ceiling(-log10(precision)))
  
  if(returnLCI==T){
    return(list("Lower Bound"=M1.lCI,"Estimate"=signif(M1,digits=3),"Upper Bound"=M1.uCI,"Length of CI"=LCI))
  }else if(returnVar==T){
    return(list("Lower Bound"=M1.lCI,"Estimate"=signif(M1,digits=3),"Upper Bound"=M1.uCI,"Variance Sum"=Var))
  }else if(returnBS==T){
    return(M1.BS)
  }else{
    return(list("Lower Bound"=M1.lCI,"Estimate"=signif(M1,digits=3),"Upper Bound"=M1.uCI))
  }
}
