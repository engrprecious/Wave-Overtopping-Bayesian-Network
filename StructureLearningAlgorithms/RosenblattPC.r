library('CDVine')
library(rvinecopulib)
library(VineCopula)

# Rosenblatt based CI test 

rosenblattCI = function(x,y,S,data){
  
  # Test conditional independence: Ui _|_ Uj | Uk
  # As inputs the function takes the uniform margins"""
  data1 = as.matrix(cbind(data[,S],data[,x]))
  data2 = as.matrix(cbind(data[,S],data[,y]))
  
  # Using the function rosenblatt to produce a dataframe of U1,...,Ud from inputs V1,...,Vd,
  # where U1=F(U1|0),..., Ud = F(Ud|U1,...,Ud-1) to produce the conditional distributions
  rosen1 = rosenblatt(data1, vinecop(data1,family_set = c("indep","gaussian","clayton","frank","student","gumbel","bb8")))
  rosen2 = rosenblatt(data2, vinecop(data2,family_set = c("indep","gaussian","clayton","frank","student","gumbel","bb8")))
  u1 = rosen1[,ncol(rosen1)]
  u2 = rosen2[,ncol(rosen2)]
  
  # Return BiCopIndTest with the input values F(Ui|Uk) and F(Uj|Uk) to test if 
  # C(F(Ui|Uk), F(Uj|Uk)) = F(Ui|Uk)*F(Uj|Uk) i.e. the independence copula
  return(BiCopIndTest(u1,u2)$p.value) 
}


# Begin implementation


Data = read.csv(".../Updated_Data.csv",header=F)

# Produce random subsample

randsample = sample(1:nrow(Data),400,replace = T)
data = Data[randsample,]

# Initialise Adjacency Matrix

A = matrix(T,8,8)
diag(A) = F
seq = seq_len(8)
sepset = lapply(seq, function(.) vector("list", 8))
A

# Unconditional Independence tests


edgesRemaining = (8*8)-8
ntests=0
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

# Conditional independence tests for k=1,..,7.

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
            pval = rosenblattCI(i,j,sep_set[,l],data)
            cat(i,'_|_',j,"|",sep_set[,l],'P-value:',pval,"\n")
            if(pval>0.05){
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



