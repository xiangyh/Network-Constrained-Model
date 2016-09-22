##==================================
##    compute the edge matrix
##==================================

computeEdgeMatrix=function(subdata,sig=0){  
  p =dim(subdata)[2]
  edges = cor(subdata,subdata)
  edges=abs(edges);
  for (i in 1:p){
    for (j in i:p){
      
      if(i!=j & abs(edges[i,j])<sig){
        edges[i,j] = 0
        edges[j,i] = 0 
      }
    }
  }  
  #edges[edges>0.5] = 1 
  #edges[edges<=0.5] = 0
  diag(edges) = 0
  edges
}

##===========================================
##  convert edge matrix to laplacian matrix
##===========================================

Edge2Laplacian<-function(edges){
  np<-dim(edges)
  p<-np[1]
  if(p<=1){
    return;
  }
  Laplacian<-matrix(rep(0,p*p),p,p)
  vertex=apply(edges,1,sum)
  for(i in 1:p){
    for(j in i:p){
      if(i==j){
        Laplacian[i,j]=1
      }
      else{
        if(vertex[i]==0||vertex[j]==0){
          Laplacian[i,j]=0;
          Laplacian[j,i]=0;
        }
        else{
          Laplacian[i,j]= - edges[i,j]/(sqrt(vertex[i])*sqrt(vertex[j]))
          Laplacian[j,i]=Laplacian[i,j] 
        }
      }
    }
  }
  return(Laplacian)
}

##===========================================
##   Extend (X,Y) to (X*, Y*) for lasso
##===========================================

argumentMatrxIndependent = function(pData,lambda2,sig)
{ 
  dimData=dim(pData)[2]
  if( dimData>= 2 ){
    edges = computeEdgeMatrix(pData,sig)
    laplacian = Edge2Laplacian(edges) 
    ev <- eigen(laplacian)    
    s = sqrt(abs(ev$values))
    s = ev$vectors %*% diag(s)
    s = sqrt(lambda2) * s
    augsubmatrix = matrix(0,dimData, dimData)
    augsubmatrix=t(s) 
    #  data = rbind(data,augsubmatrix)
    Extend = augsubmatrix
    # Data = 1/sqrt(1+lambda2) * Data
  } 
  return(Extend)
}
