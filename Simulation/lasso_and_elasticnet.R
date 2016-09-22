
#### inverse function of logit####
inv_logit_func<-function(x){
  return( 1/(1+exp(-x)) )
}

##==================================================
##    Simulation for lasso and elastic-net model
##==================================================

#### with different level of imbalanceͬ ####
gmeansum_b=NULL
gmeanTPR=NULL
gmeanTNR=NULL
gmeanERROR=NULL
gmeanCUT=NULL
gmeansum_b2=NULL
gmeanTPR2=NULL
gmeanTNR2=NULL
gmeanERROR2=NULL
gmeanCUT2=NULL

overallsum_b=NULL
overallTPR=NULL
overallTNR=NULL
overallERROR=NULL
overallCUT=NULL
overallsum_b2=NULL
overallTPR2=NULL
overallTNR2=NULL
overallERROR2=NULL
overallCUT2=NULL

## sample size is n=300, the minority is 20%

## 48 continuous variables into 12 groups
n=300
inner_size=c(3,3,3,3,4,4,4,4,5,5,5,5)
column_size=0
beta_0p=c(0,rep(-100,3),rep(0,3),rep(0,3),rep(0,3),rep(-100,4),rep(0,4),rep(0,4),rep(0,4),rep(-100,5),rep(0,5),rep(0,5),rep(0,5))
m=c(0.6,0,0,0,0.6,0,0,0,0.6,0,0,0)
length.b=length(beta_0p) 

## trial-1 ##
cutoff = seq(0.005, 0.5, 0.005)

for(iteration in 1:200){
  for(i in 1:10000){
    data=suppressWarnings(GenData(n,m,cor=0.6,inner_size,column_size,beta_0p))
    y=data[,dim(data)[2]]
    if((length(y[y==1])==n*0.1)&(length(which(is.na(y)))==0))
      break;
  }
  
  dimdata=dim(data)[2]
  samp_size=dim(data)[1]
  samp1=sample(which(y==1),0.7*length(y[y==1]))
  samp0=sample(which(y==0),0.7*length(y[y==0]))
  samp=sort(c(samp1,samp0))
  ## training set
  train=data[samp,]
  trainx=train[,-dimdata]
  trainy=train[,dimdata]
  ## test set
  test=data[-samp,]
  testx=test[,-dimdata]
  testy=test[,dimdata]
  
  trainlasso1=glmnet(scale(trainx[,-1]),trainy,family="binomial")
  lambda=c(trainlasso1$lambda,0)
  trainlasso=glmnet(scale(trainx[,-1]),trainy,family="binomial",lambda=lambda)
  cvsim=cvlasso(trainx[,-1],trainy,lambda,blocks=5,sig=0)   
  gmean=cvsim$gmeanmse
  overall=cvsim$overallmse
  
  ## maxgmean and maxoverall: the larger the better
  maxgmean=max(cvsim$gmeanmse,na.rm=TRUE)   
  maxoverall=max(cvsim$overallmse,na.rm=TRUE)
  
  ## the lambda and cutoff correspond to the maximum MSE
  for(i in 1:length(lambda)){
    if(max(cvsim$gmeanmse[i,],na.rm=TRUE)==maxgmean){
      gmean.lam=i
      gmean.cut=which.max(cvsim$gmeanmse[i,])
      break
    }
  }
  for(i in 1:length(lambda)){
    if(max(cvsim$overallmse[i,],na.rm=TRUE)==maxoverall){
      overall.lam=i
      overall.cut=which.max(cvsim$overallmse[i,])
      break
    }
  }
  
  ## cutoff=0.5
  gmean.cut2=which(cutoff==0.5)
  gmean.lam2=which.max(cvsim$gmeanmse[,gmean.cut2])
  overall.cut2=which(cutoff==0.5)
  overall.lam2=which.max(cvsim$overallmse[,overall.cut2])
  
  ## coefficients when cutoff varies
  gmeancoef=coef(trainlasso)[,gmean.lam]
  overallcoef=coef(trainlasso)[,overall.lam]
  ## coefficients when cutoff=0.5
  gmeancoef2=coef(trainlasso)[,gmean.lam2]
  overallcoef2=coef(trainlasso)[,overall.lam2]
  
  testx[,-1]=scale(testx[,-1])
  
  ## the prediction by G-Means
  # cutoff varies
  testpredgmean=inv_logit_func(testx%*%gmeancoef)
  gmeanresult=testpredgmean>=cutoff[gmean.cut]
  gmeantable=matrix(0, 2, 2)
  gmeantable[1,1]=sum(gmeanresult&(testy==1))
  gmeantable[1,2]=sum(gmeanresult&(testy==0))
  gmeantable[2,1]=sum(!gmeanresult&(testy==1))
  gmeantable[2,2]=sum(!gmeanresult&(testy==0))
  gmeanCUT=c(gmeanCUT,cutoff[gmean.cut])
  gmeanTPR=c(gmeanTPR,gmeantable[1,1]/(gmeantable[1,1]+gmeantable[2,1]))
  gmeanTNR=c(gmeanTNR,gmeantable[2,2]/(gmeantable[1,2]+gmeantable[2,2]))
  gmeanERROR=c(gmeanERROR,(gmeantable[2,1]+gmeantable[1,2])/sum(gmeantable))
  gmeansum_b=rbind(gmeansum_b,as.numeric(gmeancoef!=0))
  
  # cutoff is fixed
  testpredgmean2=inv_logit_func(testx%*%gmeancoef2)
  gmeanresult2=testpredgmean2>=cutoff[gmean.cut2]
  gmeantable2=matrix(0, 2, 2)
  gmeantable2[1,1]=sum(gmeanresult2&(testy==1))
  gmeantable2[1,2]=sum(gmeanresult2&(testy==0))
  gmeantable2[2,1]=sum(!gmeanresult2&(testy==1))
  gmeantable2[2,2]=sum(!gmeanresult2&(testy==0))
  gmeanCUT2=c(gmeanCUT2,cutoff[gmean.cut2])
  gmeanTPR2=c(gmeanTPR2,gmeantable2[1,1]/(gmeantable2[1,1]+gmeantable2[2,1]))
  gmeanTNR2=c(gmeanTNR2,gmeantable2[2,2]/(gmeantable2[1,2]+gmeantable2[2,2]))
  gmeanERROR2=c(gmeanERROR2,(gmeantable2[2,1]+gmeantable2[1,2])/sum(gmeantable2))
  gmeansum_b2=rbind(gmeansum_b2,as.numeric(gmeancoef2!=0))
  
  ## the prediction by Overall error
  # cutoff varies
  testpredoverall=inv_logit_func(testx%*%overallcoef)
  overallresult=testpredoverall>=cutoff[overall.cut]
  overalltable=matrix(0,2,2)
  overalltable[1,1]=sum(overallresult&(testy==1))
  overalltable[1,2]=sum(overallresult&(testy==0))
  overalltable[2,1]=sum(!overallresult&(testy==1))
  overalltable[2,2]=sum(!overallresult&(testy==0))
  overallCUT=c(overallCUT,cutoff[overall.cut])
  overallTPR=c(overallTPR,overalltable[1,1]/(overalltable[1,1]+overalltable[2,1]))
  overallTNR=c(overallTNR,overalltable[2,2]/(overalltable[1,2]+overalltable[2,2]))
  overallERROR=c(overallERROR,(overalltable[2,1]+overalltable[1,2])/sum(overalltable))
  overallsum_b=rbind(overallsum_b,as.numeric(overallcoef!=0))
  # cutoff is fixed
  testpredoverall2=inv_logit_func(testx%*%overallcoef2)
  overallresult2=testpredoverall2>=cutoff[overall.cut2]
  overalltable2=matrix(0,2,2)
  overalltable2[1,1]=sum(overallresult2&(testy==1))
  overalltable2[1,2]=sum(overallresult2&(testy==0))
  overalltable2[2,1]=sum(!overallresult2&(testy==1))
  overalltable2[2,2]=sum(!overallresult2&(testy==0))
  overallCUT2=c(overallCUT2,cutoff[overall.cut2])
  overallTPR2=c(overallTPR2,overalltable2[1,1]/(overalltable2[1,1]+overalltable2[2,1]))
  overallTNR2=c(overallTNR2,overalltable2[2,2]/(overalltable2[1,2]+overalltable2[2,2]))
  overallERROR2=c(overallERROR2,(overalltable2[2,1]+overalltable2[1,2])/sum(overalltable2))
  overallsum_b2=rbind(overallsum_b2,as.numeric(overallcoef2!=0))
}
gmeanb=apply(gmeansum_b,2,sum)
overallb=apply(overallsum_b,2,sum)
gmeanb2=apply(gmeansum_b2,2,sum)
overallb2=apply(overallsum_b2,2,sum)

pdf("prop(1:1).pdf")
barplot(gmeanb[-1],ylim=c(0,200),
        xlab="variables selected by G-Mean",ylab="Selection Frequency")
barplot(overallb[-1],ylim=c(0,200),
        xlab="variables selected by overall error",ylab="Selection Frequency")

barplot(gmeanb2[-1],ylim=c(0,200),
        xlab="variables selected by G-Mean(cutoff=0.5)",ylab="Selection Frequency")
barplot(overallb2[-1],ylim=c(0,200),
        xlab="variables selected by overall error(cutoff=0.5)",ylab="Selection Frequency")
dev.off()

write.table(gmeanCUT,"gmeanCUT.txt")
write.table(gmeanCUT2,"gmeanCUT2.txt")
write.table(gmeanTPR,"gmeanTPR.txt")
write.table(gmeanTPR2,"gmeanTPR2.txt")
write.table(gmeanTNR,"gmeanTNR.txt")
write.table(gmeanTNR2,"gmeanTNR2.txt")
write.table(gmeanERROR,"gmeanERROR.txt")
write.table(gmeanERROR2,"gmeanERROR2.txt")
write.table(gmeansum_b,"gmeansum_b.txt")
write.table(gmeansum_b2,"gmeansum_b2.txt")

write.table(overallCUT,"overallCUT.txt")
write.table(overallCUT2,"overallCUT2.txt")
write.table(overallTPR,"overallTPR.txt")
write.table(overallTPR2,"overallTPR2.txt")
write.table(overallTNR,"overallTNR.txt")
write.table(overallTNR2,"overallTNR2.txt")
write.table(overallERROR,"overallERROR.txt")
write.table(overallERROR2,"overallERROR2.txt")
write.table(overallsum_b,"overallsum_b.txt")
write.table(overallsum_b2,"overallsum_b2.txt")
