gaussKern <- function(x,b=10){
  K <- (1/((sqrt(2*pi))))*exp(-0.5 *(x/b)^2)
  return(K)
}

Resid = function(x,Z=c(),D=NULL){
  # Compute the residuals from regressing X onto Z
  if(length(Z)==0){
    return(x)
  }
  else if(ncol(as.matrix(Z))==1){
    resid = NULL
    for(i in 1:length(x)){
      sum = 0
      weight = 0
      for (j in 1:length(x)){
        d = sqrt(sum((Z[i]-Z[j])^2))
        k = gaussKern(d)
        sum = sum + (k*x[j])
        weight = weight + k 
      }
      resid[i] = x[i] - (sum / weight)
    }
    return(resid)
  }
  else{
    resid = NULL
    for(i in 1:length(x)){
      sum = 0
      weight = 0
      for (j in 1:length(x)){
        d = sqrt(sum((Z[i,]-Z[j,])^2))
        k = gaussKern(d)
        sum = sum + (k*x[j])
        weight = weight + k 
      }
      resid[i] = x[i] - (sum / weight)
    }
    return(resid)
  }
}

Independence = function(x,y,funclist,funcSet,alpha,D=NULL){
  pList = c()
  for(i in 1:21){
    f = funclist[[funcSet[,i][1]]]
    g = funclist[[funcSet[,i][2]]]
    r = cor(f(x),g(x))
    z = 0.5 * log((1+r)/(1-r))
    xp = (f(x)-mean(f(x))/sd(f(x)))
    yp = (g(y)-mean(g(y))/sd(g(y)))
    tau = mean((xp^2)*(yp^2))
    p = 2*(1-abs(pnorm(q = sqrt(length(x))*z,mean = 0,sd = tau)))
    pList = append(pList,p)
  }
  if(any(pList<alpha)){
    return(F)
  }else{
    return(T)
  }
}

CCI = function(x,y,Z,funclist,alpha,D=NULL){
  rx = Resid(x,Z)
  ry = Resid(y,Z)
  Independence(rx,ry,funclist,funcSet,alpha,D=NULL)
}

funclist = list()
funclist[[1]] = function(x){return(x)}
funclist[[2]] = function(x){return(x^2)}
funclist[[3]] = function(x){return(x^3)}
funclist[[4]] = function(x){return(x^4)}
funclist[[5]] = function(x){return(x^5)}
funclist[[6]] = function(x){return(x^6)}
funclist[[7]] = function(x){return(x^7)}

funcSet = combn(seq(length(funclist)),2)

##

Data = read.csv('.../Updated_Data.csv',header=F)

randsample = sample(1:nrow(Data),200,replace = T)
data = Data[randsample,]

A = matrix(T,8,8)
diag(A) = F
seq = seq_len(8)
sepset = lapply(seq, function(.) vector("list", 8))
A

edgesRemaining = (8*8)-8
ntests=0

#K=0
for (i in 1:8){
  for (j in 1:8){
    if (i != j){
      ntests = ntests+1
      pval=cor.test(data[,i],data[,j],method = 'kendall')$p.value
      cat(i,'_|_',j,'P-value:',pval,"\n")
      if(pval > 0.05){
        A[i,j] = F 
        edgesRemaining = edgesRemaining-1
        cat("Edges remaining:",edgesRemaining,"\n")
        sepset[[i]][[j]] = 'NULL'
        sepset[[j]][[i]] = 'NULL'
      }
    }
  }
}
A

# K = 1,...,7.
for(K in 1:7){
  k=K
  for (i in 1:8){
    for (j in 1:8){
      if(i!=j & A[i,j]==T){
        nbrsBool = A[i,]
        nbrsBool[j] = F
        nbrs = seq(8)[nbrsBool]
        if (length(nbrs) <= k){
          next
        }
        else{
          sep_set = combn(nbrs,m=k)
          for (l in 1:ncol(sep_set)){
            ntests = ntests+1
            cat("Number of tests:",ntests,"\n")
            if(CCI(data[,i],data[,j],Z = data[,sep_set[,l]],funclist = funclist,funcSet,alpha = 0.05)==T){
              cat("Pruning edge",i,"-",j,'\n')
              A[i,j]=F
              sepset[[i]][[j]] = sep_set[,l]
              sepset[[j]][[i]] = sep_set[,l]
              edgesRemaining=edgesRemaining-1
              cat("Edges remaining:",edgesRemaining,"\n")
              break        
            }
          }
        }
      }
    }
  }
}


