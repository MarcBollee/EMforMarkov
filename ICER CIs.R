library(Rcpp)
sourceCpp("./CoreEM.cpp")
#The purpose of this script is to test the weighted EM method to obtain a tighter confidence interval around a
#fictional ICER value than the non-weighted EM method (Craig and Sendi)

ICER<-function(M1.treatment){
  if(nrow(M1.treatment)!=3 || nrow(M1.treatment)!=3){
    print(M1.treatment)
    return(0)
  }
  if(sum(is.na(M1.treatment))!=0){
    print(M1.treatment)
    return(0)
  }
  M1.placebo<-matrix(c(0.6,0.3,0.1,
                       0.1,0.7,0.2,
                       0, 0,1),nrow=3,byrow=TRUE)
  
  #Define Costs
  costs.placebo<-c(100,200,0)
  costs.treatment<-c(500,600,0)
  utility<-c(1,0.7,0)
  maxTime<-240
  
  #Outcomes
  cost.placebo<-0
  cost.treatment<-0
  #ly.placebo<-0
  #ly.treatment<-0
  qaly.placebo<-0
  qaly.treatment<-0
  
  #Engine
  pop.placebo<-matrix(0,maxTime,3)
  pop.treatment<-matrix(0,maxTime,3)
  
  pop.placebo[1,1]=1
  pop.treatment[1,1]=1
  
  for(t in 2:maxTime){
    cost.placebo=cost.placebo+(pop.placebo[t-1,]%*%costs.placebo)
    cost.treatment=cost.treatment+(pop.treatment[t-1,]%*%costs.treatment)
    #ly.placebo=ly.placebo+((pop.placebo[t-1]%*%c(1,1,0))/12)
    #ly.treatment=ly.treatment+((pop.treatment[t-1]%*%c(1,1,0))/12)
    qaly.placebo=qaly.placebo+(pop.placebo[t-1,]%*%utility/12)
    qaly.treatment=qaly.treatment+(pop.treatment[t-1,]%*%utility/12)
    
    pop.placebo[t,]=pop.placebo[t-1,]%*%M1.placebo
    pop.treatment[t,]=pop.treatment[t-1,]%*%M1.treatment
  }
  
  ICER<-(cost.treatment-cost.placebo)/(qaly.treatment-qaly.placebo)
  return(ICER)
}

#Same as ICER but return the incremental cost and incremental QALY separately
CostBenefit<-function(M1.treatment){
  if(nrow(M1.treatment)!=3 || nrow(M1.treatment)!=3){
    print(M1.treatment)
    return(0)
  }
  if(sum(is.na(M1.treatment))!=0){
    print(M1.treatment)
    return(0)
  }
  M1.placebo<-matrix(c(0.6,0.3,0.1,
                       0.1,0.7,0.2,
                       0, 0,1),nrow=3,byrow=TRUE)
  
  #Define Costs
  costs.placebo<-c(100,200,0)
  costs.treatment<-c(500,600,0)
  utility<-c(1,0.7,0)
  maxTime<-240
  
  #Outcomes
  cost.placebo<-0
  cost.treatment<-0
  #ly.placebo<-0
  #ly.treatment<-0
  qaly.placebo<-0
  qaly.treatment<-0
  
  #Engine
  pop.placebo<-matrix(0,maxTime,3)
  pop.treatment<-matrix(0,maxTime,3)
  
  pop.placebo[1,1]=1
  pop.treatment[1,1]=1
  
  for(t in 2:maxTime){
    cost.placebo=cost.placebo+(pop.placebo[t-1,]%*%costs.placebo)
    cost.treatment=cost.treatment+(pop.treatment[t-1,]%*%costs.treatment)
    #ly.placebo=ly.placebo+((pop.placebo[t-1]%*%c(1,1,0))/12)
    #ly.treatment=ly.treatment+((pop.treatment[t-1]%*%c(1,1,0))/12)
    qaly.placebo=qaly.placebo+(pop.placebo[t-1,]%*%utility/12)
    qaly.treatment=qaly.treatment+(pop.treatment[t-1,]%*%utility/12)
    
    pop.placebo[t,]=pop.placebo[t-1,]%*%M1.placebo
    pop.treatment[t,]=pop.treatment[t-1,]%*%M1.treatment
  }
  
  return(c(cost.treatment-cost.placebo,qaly.treatment-qaly.placebo))
}
#####################################################################Start of var(ICER) analysis
#Declaring variables
BSsamples<-5000
M1.BS<-array(0,dim=c(BSsamples,3,3)) #bootstrapped transition probability matrices from bootstrapped count matrices without weighting
M1w.BS<-array(0,dim=c(BSsamples,3,3))#with weighting
ICER.dist<-matrix(0,BSsamples,2)          #bootstrapped samples of the ICER as c(cost, QALY) without weighting
ICERw.dist<-matrix(0,BSsamples,2)         #with weighting
WTP<-20000

#Define the true 1-step transition probability matrix
M1<-matrix(c(0.75,0.2,0.05,
             0.05,0.8,0.15,
             0, 0,1),nrow=3,byrow=TRUE)
Counts<-matrix(c(270,107,138,
                 21, 32 ,4,
                 54, 43 ,7),3,3,byrow=T)

set.seed(31415)

#Generate a sample of 1,2,3-step count matrices from the true matrix
Trans1step<-rMarkovCount(M=M1,s=1,rowTotals=Counts[1,])
Trans2step<-rMarkovCount(M=M1,s=2,rowTotals=Counts[2,])
Trans3step<-rMarkovCount(M=M1,s=3,rowTotals=Counts[3,])

#Use bootstrapping to find the weighting that minimises the variance of the ICER
w<-GEMforMarkov6AutoW2(Trans1step,Trans2step,Trans3step,ObjFunc = ICER, verbose=T, Forbidden = matrix(c(F,F,F,F,F,F,T,T,F),3,3,byrow=T), Counts=Counts, BSsamples=BSsamples)[[2]]

#Do a final bootstrap comparison between weighted and non-weighted
set.seed(271828)
M1.BS<-BSCIforMarkovCpp(C.all=list(Trans1step,Trans2step,Trans3step), S.all=c(1,2,3),BSsamples=BSsamples,returnBS=T,w=rep(1/3,3), Counts=Counts)
set.seed(271828)
M1w.BS<-BSCIforMarkovCpp(C.all=list(Trans1step,Trans2step,Trans3step), S.all=c(1,2,3),BSsamples=BSsamples,returnBS=T,w=w, Counts=Counts)

for(b in 1:BSsamples){
  #Run a cost benefit analysis with using the bootstrapped transition probability matrices
  ICER.dist[b,]=CostBenefit(M1.BS[b,,])
  ICERw.dist[b,]=CostBenefit(M1w.BS[b,,])
}

#Graph spread on the ICER plane
plot(x=ICER.dist[,2],y=ICER.dist[,1],
     xlab="Incremental QALY",ylab="Incremental Cost (£)", pch=20, cex=0.5, col="blue", xaxs="i",yaxs="i", las=1, cex.axis=1, cex.lab=1.3,cex.main=1.3,
     xlim=c(0,1),ylim=c(0,9000))
points(x=ICERw.dist[,2],y=ICERw.dist[,1],pch=20, cex=0.5,col="red")
lines(x=c(-0.5,1.5),y=c(WTP*(-0.5),WTP*(1.5)),col='black', lwd=2)
legend(x="topright",legend=c("Unweighted","Weighted","WTP: £20,000/QALY"),col=c("blue","red","black"),pch = c(20,20,NA),lty=c(NA,NA,1), lwd=c(NA,NA,2))

#Variance and CIs
var(ICER.dist[,1]/ICER.dist[,2])
var(ICERw.dist[,1]/ICERw.dist[,2])
quantile(ICER.dist[,1]/ICER.dist[,2],c(0.025,0.975))
quantile(ICERw.dist[,1]/ICERw.dist[,2],c(0.025,0.975))

#Graph Cost-Effectiveness Acceptability Curve
WTP.seq<-seq(0,50000,1000)
CEA<-matrix(0,length(WTP.seq),2)
for(t in 1:length(WTP.seq)){
  CEA[t,1]=sum((ICER.dist[,1]/ICER.dist[,2])<=WTP.seq[t])/BSsamples
  CEA[t,2]=sum((ICERw.dist[,1]/ICERw.dist[,2])<=WTP.seq[t])/BSsamples
}
plot(x=WTP.seq, y=CEA[,1], type='l',col="blue", lwd=2, xaxs="i",yaxs="i", las=1, cex.axis=1.3, cex.lab=1.3,cex.main=1.3,
     main="Cost-Effectiveness Acceptibility Curve", xlab="Willingness to Pay Threshold (£/QALY)", ylab="Probability of being Cost-Effective")
lines(x=WTP.seq, y=CEA[,2], col="red", lwd=2)
abline(v=WTP,lty=2)
legend(x="topleft",legend=c("Unweighted","Weighted"),col=c("blue","red"),lty=c(1,1),lwd=c(2,2))
lines(x=c(0,WTP),y=rep(CEA[which(WTP.seq==WTP),1],2),lty=2,col="blue")
lines(x=c(0,WTP),y=rep(CEA[which(WTP.seq==WTP),2],2),lty=2,col="red")
axis(side=2,at=round(CEA[which(WTP.seq==WTP),1],2),col="blue",col.axis="blue",las=1,cex.axis=1.1)
axis(side=2,at=round(CEA[which(WTP.seq==WTP),2],2),col="red",col.axis="red",las=1,cex.axis=1.1)
################################################### 
#Graph of the variance of the ICER over the 3-simplex of possible weights
set.seed(31415)
M1<-matrix(c(0.75,0.2,0.05,
             0.05,0.8,0.15,
             0, 0,1),nrow=3,byrow=TRUE)
stepsize<-0.05 #choose something that divides evenly into 1
steps<-(1/stepsize)+1
w<-rep(0,3)
varICER<-matrix(0,steps,steps)

for(i in 1:steps){
  w[1]=1-(stepsize*(i-1))
  for(j in 1:i){
    w[2]=1-w[1]-(stepsize*(j-1))
    w[3]=1-w[1]-w[2]
    print(w)
    M1.BS=BSCIforMarkovCpp(M1=M1,S.all=c(1,2,3),w=w,returnBS=T, Counts=Counts)
    ICER.dist=apply(M1.BS,1,ICER)
    varICER[i,j]=var(ICER.dist)
  }
}
#Plot the plot
startplot<-3
plot(x=c(0,stepsize,2*stepsize,3*stepsize,4*stepsize),y=varICER[5,1:5],ylim=c(0,max(varICER[5:steps])),xlim=c(0,1),xlab="Weight: w[2]=1-w[1]-w[3]",ylab="Variance of ICER",col='blue',type="l")
for(i in 5:steps){
  lines(x=seq(0,stepsize*(i-1),stepsize),y=varICER[i,(1:i)],col="blue",type='l')
}
# for(i in 1:steps){
#   lines(x=seq(0,stepsize*(i-1),stepsize),y=varICER[i,(1:i)],col="blue",type='l')
# }











