#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <math.h>
using namespace Rcpp;
using std::vector;

// Function to perform matrix multiplication
NumericMatrix mat_mult(NumericMatrix A, NumericMatrix B) {
  int n = A.nrow();
  NumericMatrix C(n, n);
  
  // Matrix multiplication: C = A * B
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      C(i, j) = 0;
      for (int k = 0; k < n; k++) {
        C(i, j) += A(i, k) * B(k, j);
      }
    }
  }
  return C;
}

// Function to calculate matrix to the power of an integer
NumericMatrix matrix_power(NumericMatrix mat, int power) {
  
  // Get the number of rows and columns (assuming a square matrix)
  int n = mat.nrow();
  
  // If power is 0, return the identity matrix
  if (power == 0) {
    NumericMatrix identity(n, n);
    for (int i = 0; i < n; i++) {
      identity(i, i) = 1;
    }
    return identity;
  }
  
  // Make a copy of the matrix to modify
  NumericMatrix result = clone(mat);
  
  // For positive powers, repeatedly multiply the matrix
  if (power > 0) {
    for (int p = 1; p < power; p++) {
      result = mat_mult(result, mat);  // Use the mat_mult function for matrix multiplication
    }
  } 
  // For negative powers, throw an error (you would need to invert the matrix here if needed)
  else if (power < 0) {
    Rcpp::stop("Negative powers are not supported without matrix inversion.");
  }
  
  return result;
}

//Function to increment the path by One
IntegerVector nextpath(IntegerVector path, int n){
  path(1)+=1;
  for(int ind=1;ind<(path.size()-1);ind++){
    if(path(ind)==n){
      path(ind)=0;
      path(ind+1)+=1;
    }else{
      break;
    }
  }
  return path;
}

//Function to row normalise a matrix, which does not break on a row of zeros
NumericMatrix rowNormalize(NumericMatrix M){
  int n = M.nrow();
  double rowsum = 0.0;
  
  
  for(int i=0;i<n;i++){
    rowsum=0.0;
    for(int j=0;j<n;j++){//compute rowsum of row i
      rowsum+=M(i,j);
    }
    if(rowsum==0.0){//if 0 then normalise to uniform vector
      for(int j=0;j<n;j++){
        double deflt = 1.0/n;
        M(i,j)=deflt;
      }//j
    }else{//else divide it by the row sum s.t. the new row sum is 1
      for(int j=0;j<n;j++){
        M(i,j)/=rowsum;
      }//j
    }//else
  }//i
  
  return M;
}

//Central EM algorithm
// [[Rcpp::export]]
NumericMatrix CoreEM(NumericMatrix Cmatrix, IntegerVector Sall, int n, NumericVector w, 
                     NumericMatrix M1naive, IntegerMatrix Forbidden, 
                     double precision, int maxiteration){
  double delta = 1.0;
  int iteration = 0;
  int S = Sall.size();
  int s =0;
  
  //Turn Cmatrix back into array
  vector<vector<vector<int>>> C(S, vector<vector<int>>(n, vector<int>(n)));
  for(int t=0; t<S; t++){
    for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
        C[t][i][j]=Cmatrix(i, ((t*n)+j));
      }
    }
  }
  
  //Define a variable to store the estimated 1-step transition count matrices
  vector<vector<vector<double>>> N(S, vector<vector<double>>(n, vector<double>(n)));
  NumericMatrix Ntotal(n,n);
  
  //Define Helpful variables
  double Increment =0.0;
  IntegerVector path(Sall(S-1)+2);
  double pathProb =0.0;
  NumericMatrix Mpows(n,n);
  int NumofPaths =0;
  
  //Define M1 related variables
  NumericMatrix M1(n,n);
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if(Forbidden(i,j)==1){
        M1(i,j)=0.0;
      }else{
        M1(i,j)=M1naive(i,j);
      }
      if(M1(i,j)==0 && Forbidden(i,j)==0){
        Forbidden(i,j)=1;
      }
    }
  }
  M1=rowNormalize(M1);
  NumericMatrix M1old(n,n);
  M1old=M1;
  NumericMatrix logM1(n,n);
  
  //Main EM Loop
  while(delta>precision && iteration<maxiteration){
    //E Step
    //Reset N and compute log of M1
    for(int k=0; k<n; k++){
      for(int l=0; l<n; l++){
        for(int t=0; t<S; t++){
          N[t][k][l]=0.0;
        }
        //logM1(k,l) = log(M1(k,l));
      }
    }
    
    
    for(int t=0; t<S; t++){//for each s-step matrix
      s=Sall(t);
      Mpows = matrix_power(M1, s);
      
      NumofPaths = pow(n,s-1);
      
      for(int k=0; k<n; k++){// for each element of this matrix
        for(int l=0; l<n; l++){
          if(Mpows(k,l)==0.0){ //will result in 0/0 otherwise
            continue;
          }
          //Reinitialize path to start with k and end with l and have only 1's in the middle
          path(0)=k;
          for(int i=1; i<s;i++){
            path(i)=0;
          }
          path(s)=l;
          
          for(int pathindex=0; pathindex<NumofPaths; pathindex++){// for each possible path
            //Compute the path probability
            pathProb=1.0;
            for(int m=0;m<s;m++){
              pathProb *= M1(path(m),path(m+1));
            }
            pathProb /= Mpows(k,l);
            
            //Compute increment for each transition count
            Increment = pathProb*((double)C[t][k][l]);
            for(int m=0;m<s;m++){
              N[t][path(m)][path(m+1)] += Increment;
            }
            
            //Increment Path by one
            path=nextpath(path,n);
          }
          
        }//l
      }//k
    }//t
    
    //Reset Ntotal
    Ntotal = NumericMatrix(n,n);
    //Take a weighted sum of the estimated 1-step counts from 1,2&3-step matrices
    for(int t=0;t<S;t++){
      for(int k=0;k<n;k++){
        for(int l=0;l<n;l++){
          Ntotal(k,l)+=w[t]*N[t][k][l];
        }
      }
    }
    //End of E step
    
    //M Step
    M1old=M1;
    for(int k=0;k<n;k++){
      for(int l=0;l<n;l++){
        if((Forbidden(k,l)==1)||(Ntotal(k,l)<0.0)){
          Ntotal(k,l)=0.0;
        }
      }
    }
    M1=rowNormalize(Ntotal);
    //End of M Step
    
    //Observe change in M1.EM since last iteration to compute delta
    delta=0.0;
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        if(abs(M1(i,j)-M1old(i,j))>delta){delta=abs(M1(i,j)-M1old(i,j));}
      }
    }
    iteration++;
    
  }//End of EM Algorithm
  
  return M1;
}//End of function