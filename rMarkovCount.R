rMarkovCount<-function(M,s,rowTotals){
  n<-nrow(M)
  C<-matrix(0,n,n)
  
  Ms<-Reduce(`%*%`, replicate(s, M, simplify = FALSE))
  #Every row of C is just a multinomial
  for(i in 1:n){C[i,]=rmultinom(1, rowTotals[i],Ms[i,])}
  return(C)
}

#OLD DEPRICATED VERSION FOR TESTING
# rMarkovCount<-function(M,s,rowTotals){
#   n<-nrow(M)
#   C<-matrix(0,n,n)
#   #If s==1 then every row of C is just a multinomial
#   if(s==1){
#     for(i in 1:n){C[i,]=rmultinom(1, rowTotals[i],M[i,])}
#     return(C)
#   }
#   #Int is an n x (s+1) matrix. The first column holds the starting population. 
#   #For each entry i of the first column(i.e. starting state) the population in i advances through time.
#   #Row i of C is the final population vector after s steps.
#   Int<-matrix(0,n,s+1)
#   Int[,1]=rowTotals
#   for(i in 1:n){
#     Int[,-1]=0 #reset matrix right of the first column
#     
#     Int[,2]=rmultinom(1,Int[i,1],M[i,])
#     for(j in 3:(s+1)){
#       for(i2 in 1:n){
#         Int[,j]=Int[,j]+rmultinom(1,Int[i2,j-1],M[i2,])
#       }
#     }
#     C[i,]=Int[,s+1]
#   }
#   return(C)
# }



