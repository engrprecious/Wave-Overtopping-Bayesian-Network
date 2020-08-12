# Gaussian BN 

gaussian = read.csv('.../Gaussian_data.csv',header=T)

# colnames(gaussian) = c('x1','x2','x3','x4','x5','x6','x7','x8')

library(bnlearn)
library("Rgraphviz")
library(pcalg)

# Creating the BN using pcalg to produce a score 

p=length(gaussian)
n=nrow(gaussian)
suffStat=list(C=cor(gaussian),n=n)
alpha=0.01
indepTest = gaussCItest
z=colnames(gaussian)

pc.fit=pc(suffStat,indepTest,alpha, labels = z)

if ( require ( Rgraphviz )) {
  plot (pc.fit , main = "",labels =z)
}

# Note: edges in [Tm.deep <-> Tm.1.0.deep] changes to [Tm.deep -> Tm.1.0.deep]
# and edge [Hm0.toe -> Rc] is reoriented to [Rc -> Hm0.toe]

dag1=model2network("[Ac][Rc|Ac][Tp.deep|Ac:Rc][Tm.deep|Tp.deep][Tm.1.0.deep|Tp.deep:Tm.deep][Hm0.deep|Tm.deep:Rc][Hm0.toe|Ac:Rc:Hm0.deep][q|Rc:Hm0.toe:Tm.1.0.deep]")
bnlearn::score(dag1,gaussian,type='loglik-g')
bnlearn::score(dag1,gaussian,type='bic-g')
graphviz.plot(dag1,shape='ellipse')

# Compare with fit using iamb algorithm

pdag2 = iamb(gaussian)
undirected.arcs(pdag2)
dag2 = set.arc(pdag2,"Tm.deep","Tm.1.0.deep")
dag2 = set.arc(dag2,"Tp.deep","Tm.deep")
dag2 = set.arc(dag2,"Hm0.toe","q")
dag2= set.arc(dag2,"Rc","q")
dag2 = set.arc(dag2,"Ac","Rc")
dag2 = set.arc(dag2,"Tm.1.0.deep","q")
bnlearn::score(dag2,gaussian,type='loglik-g')
bnlearn::score(dag2,gaussian,type='bic-g')
graphviz.plot(dag2,shape='ellipse')

# Compare with fit using hc algorithm

pdag3 = hc(gaussian)
fitHC = bn.fit(pdag3,gaussian)
undirected.arcs(pdag3)
bnlearn::score(pdag3,gaussian,type='loglik-g')
bnlearn::score(pdag3,gaussian,type='bic-g')
graphviz.plot(pdag3,shape='ellipse')


## Hybrid Bayesian Network with discharge discretised

clgaussian = read.csv('.../HybridData.csv',header=T)
dag = hc(clgaussian)
fit = bn.fit(dag,clgaussian)
fit 
graphviz.plot(fit,shape='ellipse')
bnlearn::score(dag,clgaussian,'loglik-cg')
bnlearn::score(dag,clgaussian,'bic-cg')




