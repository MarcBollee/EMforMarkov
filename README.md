# **Introduction**  
This project comprises code which is used to fit a discrete Markov model with missing data. 
A weighted expectation-maximation (EM) algorithm is used to obtain an estimated transition probability matrix.
To obtain confidence intervals, the parametric bootstrap is used.

This set of function was developed in [RStudio 2024.09.1 Build 394](https://posit.co/download/rstudio-desktop/), using [R 4.4.2](https://cran.r-project.org/bin/windows/base/).  
In addition, the following libraries from the CRAN repository were used:  
[Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html): For integrating C++ code.  
[LaplacesDemon](https://cran.r-project.org/web/packages/LaplacesDemon/index.html): For its rdirichlet function which allows sampling of Dirichlet random variables. This function is not used in the scripts, but can be useful for generating random transition matrices to test the method on.

# **File Descriptions**  
## **rMarkovCount.R**  
This contains a function which generates a random transition count matrix given an initial population vector, a transition probability matrix, and the number of time steps ellapsed.

## **EMforMarkov.R**  
This file contains the function which performs the weigthed EM algorithm in R. The function takes a input a list (or array) of count matrices "C.all", and list of labels "S.all" for how many time steps is associated with each count matrix.  
In addition, a set of weights "w" may be specified, otherwise uniform weights will be used.  
"Forbidden" is matrix of binary values of the same size are the count matrices, where Forbidden[i,j]=TRUE means that the transition i->j is not allowed in 1 time step, and therefore the estimated transition probability will automatically be set to 0.  
"precision" and "maxiteration" define the stoping criteria of the while loop that performs the EM algorithm.  
The function returns the estimated 1-step transition probability matrix. 

## **CoreEM.cpp**  
This file contains the code used for performing the weighted EM algorithm in C++. The logic is exactly the same as in EMforMarkov.R, but omits all the initial checks of the inputs. It is recommended that this function be used only within an environment where the input variables have already gone through the necessary checks. It is much faster than the equivalent R version of the code, and is therefore used for the computationally intensive bootstrap. In order to load this function to the R Studio environment, one must install the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) package and run the following lines of code:  
library(Rcpp)  
sourceCpp("./CoreEM.cpp") or whichever path contains the CoreEM.cpp file

## **BSCIforMarkov.R**  
_(prerequisites: rMarkovCount, EMforMarkov, CoreEM)_  
This function is used to generate bootstrapped confidence intervals for the estimated transition probabiltiy matrix. 
The function admits as input either the same input data as EMforMarkov, or a specific transition probability matrix with associated transition counts. 
The function by default returns a list of matrices with the lower bounds, estimates, and upper bounds respecively. 
returnVar=TRUE => function in addition returns the sum of the variances of all components of the transition probability matrix.  
returnLCI=TRUE => function in addition returns the sum of the length of the (1-alpha)% confidence intervals of all components of the transition probability matrix. 
  
returnBS=TRUE => function only returns all sampled bootstrapped transition probability matrices. 

## **Auto Weighting any Function.R**  
_(prerequisites: rMarkovCount, EMforMarkov, CoreEM, BSCIforMarkov)_  
This function takes as input transition count matrices for 1, 2 and 3-time step count matrices, and an arbitrary function of the transition probability matrix. The function attempts to find the set of weights to minimises the vriance of the output of that function given all bootstrapped estimates. It finds these weights using a besic grid search startegy over the 3 simplex. Due to the probabilitistic nature of the bootstrapped variance estimate, this algorithm is not guarenteed to converge.  
Note the use of specific seeding before each use of "BSCIforMarkov". This is to ensure that all 4 new proposed weights are matched up with the same bootstrapped data and have an even playing field.

## **ICER CIs.R**  
_(prerequisites: rMarkovCount, EMforMarkov, CoreEM, BSCIforMarkov, AutoWforMarkov)_  
This file is a script to reproduce the results given in the paper. A ficticious ICER is computed from a set of count matrices, and the bootstrapped variance of the estimate is obtained. First a set of count matrices is generated from the true transition probability matrix. Then the optimal set of weights is found using AutoWforMarkov. Finally, the variance of the ICER between weighted and non-weighted estimates are compared and plotted.  
There is an extra section at the end of the script which graphs the variance of the ICER over the entire 3-simplex of possible weights. This can be used in addition to observe that the variance is convex over this space, and there is a unique minimum. 
