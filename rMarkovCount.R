rMarkovCount<-function(M,s,rowTotals){
  n<-nrow(M)
  C<-matrix(0,n,n)
  
  Ms<-Reduce(`%*%`, replicate(s, M, simplify = FALSE))
  #Every row of C is just a multinomial
  for(i in 1:n){C[i,]=rmultinom(1, rowTotals[i],Ms[i,])}
  return(C)
}
