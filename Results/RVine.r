# R-Vine 
library(VineCopula)

data = read.csv('.../Updated_Data.csv',header=F)

rvine = RVineStructureSelect(data)

rvine$Matrix

