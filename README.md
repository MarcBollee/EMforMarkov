**Introduction**
This project comprises code which is used to fit a discrete Markov model with missing data. 
A weighted expectation-maximation (EM) algorithm is used to obtain an estimated transition probability matrix.
To obtain confidence intervals, the parametric bootstrap is used.

This set of function was developed in RStudio.

**File Descriptions**
**rMarkovCount.R**
This contains a function which generates a random transition count matrix given an initial population vector, a transition probability matrix, and the number of time steps ellapsed.

**EM for MM 7.R**
This file contains the function which performs the weigthed EM algorithm in R. The function takes a input a list (or array) of count matrices "C.all", and list of labels "S.all" for how many time steps is associated with each count matrix. 
In addition, a set of weights "w" may be specified, otherwise uniform weights will be used.
"Forbidden" is matrix of binary values of the same size are the count matrices, where Forbidden[i,j]=TRUE means that the transition i->j is not allowed in 1 time step, and therefore the estimated transition probability will automatically be set to 0.
"precision" and "maxiteration" define the stoping criteria of the while loop that performs the EM algorithm.
The function returns the estimated 1-step transition probability matrix. 

**Core EM.cpp**
This file contains the code used for performing the weighted EM algorithm in C++. The logic is exactly the same as in EM for MM 7.R, but omits all the initial checks of the inputs. It is recommended that this function be used only within an environment where the input variables have already gone through the necessary checks. It is much faster than the equivalent R version of the code, and is therefore used for the intensive bootstrap function. In order to load this function to the R Studio environment, one must install the "Rcpp" package and run the following lines of code:
library(Rcpp)
sourceCpp("./CoreEM.cpp")

**BSCI for MarkovCpp.R** _(prerequisites: rMarkovCount, GEMforMarkov, CoreEM)_
This function is used to generate bootstrapped confidence intervals for the estimated transition probabiltiy matrix. 
The function admits as input either the same input data as GEMforMarkov7, or a specific transition probability matrix with associated transition counts. 
The function by default returns a list of matrices with the lower bounds, estimates, and upper bounds respecively. 
returnVar=TRUE => function returns the sum of the variances of all components of the transition probability matrix. 
returnLCI=TRUE => function returns the sum of the length of the (1-alpha)% confidence intervals of all components of the transition probability matrix. 
returnBS=TRUE => function returns all sampled bootstrapped transition probability matrices. 

**Auto Weighting any Function.R** _(prerequisites: rMarkovCount, GEMforMarkov, CoreEM, BSCIforMarkov)_
This function takes as input transition count matrices for 1, 2 and 3-time step count matrices, and an arbitrary function of the transition probability matrix. The function attempts to find the set of weights to minimises the vriance of the output of that function given all bootstrapped estimates. It finds these weights using a besic grid search startegy over the 3 simplex. Due to the probabilitistic nature of the bootstrapped variance estimate, this algorithm is not guarenteed to converge. 

**ICER CIs.R** _(prerequisites: rMarkovCount, GEMforMarkov, CoreEM, BSCIforMarkov,GEMforMarkov6AutoW2)_
This file is a script to reproduce the results given in the paper. A ficticious ICER is computed from a set of count matrices, and the bootstrapped variance of the estimate is obtained. First a set of count matrices is generated from the true transition probability matrix. Then the optimal set of weights is found using GEMforMarkov6AutoW2. Finally, the variance of the ICER between weighted and non-weighted estimates are compared and plotted. Finally, there is an extra section at the end which graphs the variance of the ICER over the entire 3-simplex. 
