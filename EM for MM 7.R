GEMforMarkov7<-function(C.all, #list of matrices for the one and two-step transition counts all of the same square dimension
                        S.all, #for each matrix in C.all[s], S.all[s] indicates that it is an s-step matrix
                        w,     #weights
                        Forbidden,          #matrix of forbidden transitions: Forbidden[i,j]==T implies i->j is forbidden
                        precision=0.001, maxiteration=20, 
                        verbose=FALSE){
  
  ################################################################## Setup and Errors
  #Turn C.all into array for speed
  if(is.list(C.all)){
    n<-nrow(C.all[[1]])
    C<-array(0,c(length(C.all),n,n))
    for(t in 1:length(C.all)){C[t,,]=C.all[[t]]}
  }else{
    n<-nrow(C.all[1,,])
    C<-C.all
  }

  TotalrowSums<-apply(C,c(1,3),sum)

  if(is.list(S.all)){unlist(S.all)}
  S<-length(S.all)
  
  #Check all S.all are increasing integers >0
  if(sum(S.all!=as.integer(S.all))>0){stop("S.all must be a list of positive integers.")}
  for(t in 2:S){if(S.all[t-1]>=S.all[t]){stop("S.all must be strictly increasing.")}}
  
  #Check C.all and S.all are of equal length
  if(S!=dim(C)[1]){stop("C.all and S.all must be of equal length.")}
  
  for(t in 1:S){
    #Check dimensions of C.all
    if((nrow(C[t,,])!=n) || (ncol(C[t,,])!=n)){stop("All matrices in C.all must be of the same square dimention.")}
  }
  #Check there is enough information in each row
  for(i in 1:n){
    if(TotalrowSums[i]==0){stop(paste("Not enough information about row",i))}
  }
  
  #Make sure Forbidden is matrix of TRUE/FALSE values and define Forbidden if not provided
  if(missing(Forbidden)){
    Forbidden<-matrix(F,n,n)
  }else{
    if(nrow(Forbidden)!=n || ncol(Forbidden)!=n){stop("Forbidden matrix must be of the same size as the transition count matrices.")}
    for(i in 1:n){
      for(j in 1:n){
        if(Forbidden[i,j]!=0 && Forbidden[i,j]!=1){
          stop(paste("Forbidden must be a matrix of binary TRUE/FALSE indicators. Forbidden[",i,",",j,"] is not binary."))
        }
      }
    }
  }
  Forbidden=(Forbidden==1)#Ensure values are boolean
  
  #Define w as uniform if missing
  if(missing(w)){w<-rep(1/S,S)}
  
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
  #Define a function for quickly exponentiating a matrix
  pow = \(x, n) Reduce(`%*%`, replicate(n, x, simplify = FALSE))
  
  #Define a function that increments a path
  nplusone<-n+1
  incrementPath<-function(path){
    len<-length(path)
    path[2]=path[2]+1
    for(ind in 2:len){
      if(path[ind]==nplusone){
        path[ind]=1
        path[ind+1]=path[ind+1]+1
      }else{
        break
      }
    }
    return(path)
  }

  ###############################################################Starting Naive Estimate
  #Calculate a naive estimate of M1 as a "good" and computationally cheap starting point
  M1.naive<-array(0,c(S,n,n))
  
  if(S.all[1]==1){
    M1.naive[1,,]=rowNormalise(C[1,,])
    if(S==1){
      if(verbose==TRUE){warning("Returned MLE based on the 1-step transition matrices provided.")}
      return(M1.naive[1,,])
    } #If only a 1-step matrix is provided then we are done
  }
  
  for(t in (1+(S.all[1]==1)):S){
    #Get eigen vectors and eigen values
    eigen<-eigen(rowNormalise(C[t,,]))
    #Only attempt to calculate the square root of the matrix if all eigenvalues are non-negative and the determinant is not close to 0
    if((sum(is.complex(eigen$values))==0) && (sum(eigen$values<0)==0 || (S.all[t]%%2==0)) && (det(eigen$vectors)>=1e-10)){
      M1.naive[t,,]=eigen$vectors%*%diag((eigen$values)^(1/S.all[t]))%*%solve(eigen$vectors) #Nth root estimate
    }
  }
  
  ############################################################# EM Algorithm Initialisation
  #Initialise loop variables
  delta<-1                          #absolute change in M1.EM since last iteration
  iteration<-1                      #keeps track of while-loop iteration
  
  if(S.all[1]==1){                  #peel off C[1,,] as it behaves differently to the rest
    N<-array(0,c(S-1,n,n))          #stores all estimates of the 1-step count
    N1<-C[1,,]
    C=C[-1,,]
    S.EM<-S.all[-1]
    S=S-1
  }else{
    N<-array(0,c(S,n,n))
    S.EM<-S.all
  }
  N.total<-matrix(0,n,n)            #stores weighted average of estimated 1-step counts
  
  Increment<-0                      #convenience variable, used for sums
  path<-rep(1,S-1)                  #contains current path as a vector of intermediate states between k and l 
  pathProb<-0                       #convenience variable, used for path probability calculation
  Mpows<-matrix(0,n,n)              #contains M1.EM to the power of the current number of steps
  
  #M1.EM initialised as the mean of the naive transition probability matrices
  M1.EM<-apply(M1.naive,c(2,3),sum)
  M1.EM[which(Forbidden)]=0        #Set forbidden transitions to 0 and renormalise
  M1.EM[which(M1.EM<0)]=0.01
  M1.EM=rowNormalise(M1.EM)
  M1.EMold<-M1.EM                  #used for computing delta
  
  #If any component of M1.EM is initialised as 0, it will mathematically remain so throughout the algorithm,
  #so warn user about this if not specified in "Forbidden"
  for(i in 1:n){
    for(j in 1:n){
      if(M1.EM[i,j]==0 && !Forbidden[i,j]){
        if(verbose==TRUE){warning(paste("Transition (", i, "->", j, ") is forbidden by default eventhough not specified as such"))}
        Forbidden[i,j]=T #set to 1 so that the EM loop automatically skips the [i,j] entry and does not waste time computing all 0's
      }
    }
  }
  
  ############################################################# Main EM Algorithm Loop
  #Print start for user
  if(verbose==TRUE){message("Starting EM Algorithm")}
  
  while(delta>precision && iteration<maxiteration){
    
    #E Step: Estimate and impute the number of one-step transition counts contained inside C
    #Compute N
    for(t in 1:S){
      s=S.EM[t]
      #Compute powers of M1.EM
      Mpows=pow(M1.EM,s)
      
      Numofpaths=n^(s-1) #number of paths which must be explored
      
      for(k in 1:n){
        for(l in 1:n){
          
          if(Mpows[k,l]==0){ #will result in 0/0 otherwise
            next
          }
          
          path=c(k,rep(1,s-1),l) #reset to base path
          for(pathIndex in 1:Numofpaths){
            #Compute the path probability
            pathProb=1
            for(m in 1:s){
              pathProb=pathProb*M1.EM[path[m],path[m+1]]
            }
            pathProb=pathProb/Mpows[k,l]
            
            #Compute increment for each transition count
            Increment=C[t,k,l]*pathProb
            
            #Increment every transition in this path
            for(m in 1:s){
              N[t,path[m],path[m+1]]=N[t,path[m],path[m+1]]+Increment
            }

            path=incrementPath(path) #go to next path
          }
        }
      }
      if(verbose==T){message(paste("Completed step",s))}
    }

    #Take a weighted sum of the estimated 1-step counts from 1,2&3-step matrices
    if(S.all[1]==1){
      N.total=w[1]*N1
      for(t in 1:S){
        N.total=N.total+(w[t+1]*N[t,,])
      }
    }else{
      N.total=0
      for(t in 1:S){
        N.total=N.total+(w[t]*N[t,,])
      }
    }
    
    #End of E Step
    
    #M Step: Find MLE of 1-step transition probability based on imputed transition counts (i.e. row-normalise transition count matrix)
    M1.EMold=M1.EM
    M1.EM=N.total
    M1.EM[which(Forbidden)]=0
    M1.EM=rowNormalise(M1.EM)
    
    #Observe change in M1.EM since last iteration to compute delta
    delta=max(abs(M1.EM-M1.EMold))
    if(verbose==TRUE){
      message(paste("Iteration: ",iteration))
      message(paste("Delta: ",delta))
    }
    iteration=iteration+1
    #End of M Step
  }#End of EM Loop
  
  #Sometimes components of M1.EM can be negative values close to machine epsilon
  M1.EM[which(M1.EM<0)]=0
  
  
  if(verbose==TRUE){message("Done!")}
  return(round(M1.EM,ceiling(-log10(precision))))
}
