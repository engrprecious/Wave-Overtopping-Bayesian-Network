# R-Vine 
library(VineCopula)

data = read.csv('/Users/jamie/Desktop/Diss Resources/Jamie/Updated_Data.csv',header=F)

rvine = RVineStructureSelect(data)

rvine$Matrix

