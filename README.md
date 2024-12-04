This project comprises code which is used to fit a discrete Markov model with missing data. 
The expectation-maximation (EM) is used to obtain an estimated transition probability matrix, as can be seen in "EM for MM 7.R" and CoreEM.cpp".
To obtain confidence intervals, the parametric bootstrap is used, as can be seen in "BSCI for MarkovCpp.R".

The file "rMarkovCount.R" is a function which generates count data given an initial population vector, a transition probability matrix, and the number of time steps ellapsed.

To illustrate the utility of weighting the estimates, an example is provided via a ficticious ICER calculation in "ICER CIs.R". 
The analysis uses "Auto Weighting any Function.R" which automatically finds the optimal weights to minimise the variance of any target function which is solely dependent on the transition probability matrix.
The results of the ICER analysis can be seen on the ICER plane, in a cost-effectiveness acceptibility curve, and finally a look at the ICER variance over the entire simplex of possible weights.
