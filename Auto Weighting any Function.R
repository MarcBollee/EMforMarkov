GEMforMarkov6AutoW2<-function(Onesteps=list(), Twosteps=list(), Threesteps=list(), #list of matrices for the one and two-step transition counts all of the same square dimension
                             w, Counts,
                             Forbidden,          #matrix of forbidden transitions
                             precision=0.01, maxiteration=20, BSsamples=5000,
                             ObjFunc,
                             verbose=FALSE){
  
  ################################################################## Setup and Errors
  #If the matrices are not not provided as a list, make it into a list
  if(!is.list(Onesteps)){Onesteps<-list(Onesteps)}
  if(!is.list(Twosteps)){Twosteps<-list(Twosteps)}
  if(!is.list(Threesteps)){Threesteps<-list(Threesteps)}
  
  #Check that all the matrices are square and of the same dimensions
  numOnesteps<-length(Onesteps)
  numTwosteps<-length(Twosteps)
  numThreesteps<-length(Threesteps)
  
  if(numOnesteps!=0){
    n<-nrow(Onesteps[[1]])
    ntakenfrom<-"One"
  }else if(numTwosteps!=0){
    n<-nrow(Twosteps[[1]])
    ntakenfrom<-"Two"
  }else if(numThreesteps!=0){
    n<-nrow(Threesteps[[1]])
    ntakenfrom<-"Three"
  }else{stop("All lists of matrices are empty.")}
  
  TotalOnesteps<-matrix(0,n,n)
  TotalTwosteps<-matrix(0,n,n)
  TotalThreesteps<-matrix(0,n,n)
  
  TotalOnestepsRowSums<-rep(0,n)
  TotalTwostepsRowSums<-rep(0,n)
  TotalThreestepsRowSums<-rep(0,n)
  
  if(numOnesteps!=0){
    for(t in 1:numOnesteps){
      #Check dimensions
      if(nrow(Onesteps[[t]])!=ncol(Onesteps[[t]])){stop(paste("One-Step matrix ", t, " is not square."))}
      if(nrow(Onesteps[[t]])!=n){stop(paste("One-Step matrix", t, "has different number of rows to", ntakenfrom,"-Step matrix 1."))}
      #Check that there are no negative values in the counts
      if(sum(Onesteps[[t]]<0)!=0){stop(paste("One-Step matrix", t, "contains negative values"))}
      #Take any NA values as 0
      Onesteps[[t]][which(is.na(Onesteps[[t]]))]<-0
      #Check that there is enough information for each row (later)
      TotalOnestepsRowSums=TotalOnestepsRowSums+rowSums(Onesteps[[t]])
      #Sum up the 1-step matrices, which becomes the new count matrix
      TotalOnesteps=TotalOnesteps+Onesteps[[t]]
    }
  }
  if(numTwosteps!=0){
    for(t in 1:numTwosteps){
      #check dimensions
      if(nrow(Twosteps[[t]])!=ncol(Twosteps[[t]])){stop(paste("Two-Step matrix", t, "is not square."))}
      if(nrow(Twosteps[[t]])!=n){stop(paste("Two-Step matrix", t, "has different number of rows to ", ntakenfrom,"-Step matrix 1."))}
      #Check that there are no negative values in the counts
      if(sum(Twosteps[[t]]<0)!=0){stop(paste("Two-Step matrix", t, "contains negative values"))}
      #Take any NA values as 0
      Twosteps[[t]][which(is.na(Twosteps[[t]]))]<-0
      #Check that there is enough information for each row (later)
      TotalTwostepsRowSums=TotalTwostepsRowSums+rowSums(Twosteps[[t]])
      #Sum up the 1-step matrices, which becomes the new count matrix
      TotalTwosteps=TotalTwosteps+Twosteps[[t]]
    }
  }
  if(numThreesteps!=0){
    for(t in 1:numThreesteps){
      #check dimensions
      if(nrow(Threesteps[[t]])!=ncol(Threesteps[[t]])){stop(paste("Two-Step matrix ", t, " is not square."))}
      if(nrow(Threesteps[[t]])!=n){stop(paste("Two-Step matrix ", t, " has different number of rows to ", ntakenfrom,"-Step matrix 1."))}
      #Check that there are no negative values in the counts
      if(sum(Threesteps[[t]]<0)!=0){stop(paste("Two-Step matrix ", t, " contains negative values"))}
      #Take any NA values as 0
      Threesteps[[t]][which(is.na(Threesteps[[t]]))]<-0
      #Check that there is enough information for each row (later)
      TotalThreestepsRowSums=TotalThreestepsRowSums+rowSums(Threesteps[[t]])
      #Sum up the 1-step matrices, which becomes the new count matrix
      TotalThreesteps=TotalThreesteps+Threesteps[[t]]
    }
  }
  
  #Check that there is enough information for each row
  for(i in 1:n){
    if(TotalOnestepsRowSums[i]+TotalTwostepsRowSums[i]+TotalThreestepsRowSums[i]==0){
      stop(paste("Not enough information for row ", i))
    }
  }
  
  #Make sure Forbidden is matrix of TRUE/FALSE values and define Forbidden if not provided
  if(missing(Forbidden)){
    Forbidden<-matrix(0,n,n)
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
  
  if(missing(w)){w<-rep(1/n,n)}
  if(missing(Counts)){Counts<-matrix(c(TotalOnestepsRowSums,TotalTwostepsRowSums,TotalThreestepsRowSums),length(S.all),n)}
  M1<-GEMforMarkov7(list(TotalOnesteps,TotalTwosteps,TotalThreesteps),S.all = c(1,2,3),w=w)
  
  ################################################################## Finding Optimal Weight
  M1.BS<-BSCIforMarkovCpp(M1=M1,S.all=c(1,2,3),returnBS=T,w=w, BSsamples=BSsamples, Counts=Counts)
  Obj<-rep(0,BSsamples)
  Obj=apply(M1.BS, MARGIN = 1, ObjFunc)
  VarObj<-var(Obj)
  seed<-0
  
  #Recursive narrowing starting with a step size of 0.1
  step<-0.1
  while(step>=precision){
    #Recompute M1 with new weight
    M1<-GEMforMarkov7(list(TotalOnesteps,TotalTwosteps,TotalThreesteps),S.all = c(1,2,3),w=w)
    
    #Propose changes to the weight set
    Testw.north<-w+c(step,-step/2,-step/2)
    Testw.south<-w+c(-step,step/2,step/2)
    Testw.east<-w+c(-step/2,step,-step/2)
    Testw.west<-w+c(step/2,-step,step/2)
    
    #Store proposed weights in a list
    Testw<-list(Testw.north, Testw.south, Testw.east, Testw.west)
    
    #Sample a random integer to set as seed s.t. all bootstraps generate the same bootstrap data
    seed=sample(1:1000000,1)
    
    #Compute the bootstrapped var of the Objective function given for each proposed weight
    set.seed(seed)
    M1.BS<-BSCIforMarkovCpp(M1=M1,S.all=c(1,2,3),BSsamples=BSsamples,returnBS=T,w=Testw.north, Counts=Counts)
    Obj=apply(M1.BS, MARGIN = 1, ObjFunc)
    VarObj.north<-var(Obj)
    
    set.seed(seed)
    M1.BS<-BSCIforMarkovCpp(M1=M1,S.all=c(1,2,3),BSsamples=BSsamples,returnBS=T,w=Testw.south, Counts=Counts)
    Obj=apply(M1.BS, MARGIN = 1, ObjFunc)
    VarObj.south<-var(Obj)
    
    set.seed(seed)
    M1.BS<-BSCIforMarkovCpp(M1=M1,S.all=c(1,2,3),BSsamples=BSsamples,returnBS=T,w=Testw.east, Counts=Counts)
    Obj=apply(M1.BS, MARGIN = 1, ObjFunc)
    VarObj.east<-var(Obj)
    
    set.seed(seed)
    M1.BS<-BSCIforMarkovCpp(M1=M1,S.all=c(1,2,3),BSsamples=BSsamples,returnBS=T,w=Testw.west, Counts=Counts)
    Obj=apply(M1.BS, MARGIN = 1, ObjFunc)
    VarObj.west<-var(Obj)
    
    set.seed(seed)
    M1.BS<-BSCIforMarkovCpp(M1=M1,S.all=c(1,2,3),BSsamples=BSsamples,returnBS=T,w=w, Counts=Counts)
    Obj=apply(M1.BS, MARGIN = 1, ObjFunc)
    VarObj<-var(Obj)
    
    #Store the variances in a vector
    TestVarObj<-c(VarObj.north, VarObj.south, VarObj.east, VarObj.west)
    
    if(VarObj.north>VarObj && VarObj.south>VarObj && VarObj.east>VarObj && VarObj.west>VarObj){
      step=step/2 #If no change is accepted, divide the step size by 2
    }else{
      w=Testw[[which.min(TestVarObj)]] #Else replace the current weights with ones that give lower var
      VarObj=min(TestVarObj)
    }
    if(verbose==T){ #message for user
      message(paste("Weight:",round(w,3)))
      message(paste("Objective:",round(VarObj,3)))
      message(paste("Step:", step))
    }
  }
  
  return(list("Estimate"=round(M1,ceiling(-log10(precision))),"Weight"=w))
}
