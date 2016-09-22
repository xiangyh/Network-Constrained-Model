
##==================================
## Generate Continuous Variables
##==================================

GenData_Cont= function (n,m,cor,var_size)
{ 
  X=NULL
  # Specify the number of groups
  # var_size is the number of groups in each group
  group_num=length(var_size)  
  if(group_num==0){
    return(X);
  }
  else{
    for(i in 1:group_num){
      mu=rep(m[i],var_size[i])
      cormatrix = matrix(rep(cor,var_size[i]*var_size[i]),var_size[i],var_size[i])  ##Э????????
      diag(cormatrix)=rep(1,var_size[i])
      X = cbind(X, rmvnorm(n,mu,cormatrix) ) ## multinormal density, mean 0
    }
    return(X);
  }
}

##==================================
## Generate Discrete Variables
##==================================

GenData_Discrete<-function(n,cor,column_size)
{
  # column size is the level of each varaible
  # the number of dummy variable is then column_size-1
  column_number=length(column_size) 
  means=rep(0,column_number)
  cormatrix=matrix(rep(cor,column_number*column_number),column_number,column_number)
  diag(cormatrix)=rep(1,column_number)
  X=rmvnorm(n,means,cormatrix) 
  XX=NULL
  for(i in 1:column_number){
    sdx=sqrt(var(X[,i],X[,i]))
    meanx=mean(X[,i])
    for(j in 1:n){
      # 
      X[j,i]=ceiling(pnorm(X[j,i],meanx,sdx)*column_size[i]) 
    }
    this_X=1*(matrix(X[,i],n,(column_size[i]-1))==matrix(2:column_size[i],n
                                                         ,(column_size[i]-1),byrow=TRUE))   
    XX=cbind(XX,this_X)
  }
  XX
}


##==================================
## Generate Data including Y
##==================================

GenData<-function(n,m,cor,inner_size,column_size,beta_0p) ## beta_0p is artificially specified
{
  X1=GenData_Cont(n,m,cor,inner_size) ## continuous variable
  if(column_size==0){
    X=cbind(1,X1)
  }else{
    X2=GenData_Discrete(n,cor,column_size) ## discrete variable
    X=cbind(1,X1,X2) ## the final design matrix
  }
  output=X
  ## randomly generate Y
  Y=rbinom(n,1,exp(X%*%beta_0p)/(1+exp(X%*%beta_0p))) 
  output=cbind(X,Y) ## the final simulated dataset
}


