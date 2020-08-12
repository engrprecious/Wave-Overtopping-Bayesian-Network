##
library(CondIndTests)

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
            pval = CondIndTest(data[,i],data[,j],data[,sep_set[,l]])$pvalue
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

# Save data, adjacency matrix and separation set to use for orientation and copula fitting

saveRDS(A,".../Adjacency_Matrx.rds")
saveRDS(sepset,".../Sepset.rds")
saveRDS(data,".../PCBN_Data.rds")



