PCCBN_Data <- readRDS(".../PCCBN_Data.rds")

library(VineCopula)

u1 = PCCBN_Data[,1]
u2 = PCCBN_Data[,2]
u3 = PCCBN_Data[,3]
u4 = PCCBN_Data[,4]
u5 = PCCBN_Data[,5]
u6 = PCCBN_Data[,6]
u7 = PCCBN_Data[,7]
u8 = PCCBN_Data[,8]

# Construct Copulas 

c15 = BiCopSelect(u1,u5)
c32 = BiCopSelect(u3,u2)
c71 = BiCopSelect(u7,u1)
c42 = BiCopSelect(u4,u2)
c85 = BiCopSelect(u8,u5)

#

h15 = BiCopHfunc(u1,u5,c15)
c65 = BiCopSelect(u6,u5)
h65 = BiCopHfunc(u6,u5,c65)
c16g5 = BiCopSelect(h15$hfunc2,h65$hfunc2)

h32 = BiCopHfunc(u3,u2,c32)
c52 = BiCopSelect(u5,u2)
h52 = BiCopHfunc(u5,u2,c52)
c35g2 = BiCopSelect(h32$hfunc2,h52$hfunc2)


h71 = BiCopHfunc(u7,u1,c71)
c51 = BiCopSelect(u5,u1)
h51 = BiCopHfunc(u5,u1,c51)
c75g1 = BiCopSelect(h71$hfunc2,h51$hfunc2)
  

h42 = BiCopHfunc(u4,u2,c42)
c62 = BiCopSelect(u6,u2)
h62 = BiCopHfunc(u6,u2,c62)
c46g2 = BiCopSelect(h42$hfunc2,h62$hfunc2)

h85 = BiCopHfunc(u8,u5,c85)
c35 = BiCopSelect(u3,u5)
h35 = BiCopHfunc(u3,u5,c35)
c83g5 = BiCopSelect(h85$hfunc2,h35$hfunc2)


#

h7g51 = BiCopHfunc(h71$hfunc2,h51$hfunc2,c75g1)
c61g5 = BiCopSelect(h65$hfunc2,h15$hfunc2)
h6g51 = BiCopHfunc(h65$hfunc2,h15$hfunc2,c61g5)
c76g51 = BiCopSelect(h7g51$hfunc2,h6g51$hfunc2)
  
h4g62 = BiCopHfunc(h42$hfunc2,h62$hfunc2,c46g2)
c12 = BiCopSelect(u1,u2)
h12 = BiCopHfunc(u1,u2,c12)
c16g2 = BiCopSelect(h12$hfunc2,h62$hfunc2)
h1g62 = BiCopHfunc(h12$hfunc2,h62$hfunc2,c16g2)
c41g62 = BiCopSelect(h4g62$hfunc2,h1g62$hfunc2)
  
h8g35 = BiCopHfunc(h85$hfunc2,h35$hfunc2,c83g5)
h25 = BiCopHfunc(u2,u5,c52)
h35 = BiCopHfunc(u3,u5,c35)
c23g5 = BiCopSelect(h25$hfunc2,h35$hfunc2)
h2g35 = BiCopHfunc(h25$hfunc2,h35$hfunc2,c23g5)
c82g35 = BiCopSelect(h8g35$hfunc2,h2g35$hfunc2)

# 


c45 = BiCopSelect(u4,u5)
h45 = BiCopHfunc(u4,u5,c45)
c43g5 = BiCopSelect(h45$hfunc2,h35$hfunc2)
h4g35 = BiCopHfunc(h45$hfunc2,h35$hfunc2,c43g5)
h2g35 = BiCopHfunc(h25$hfunc2,h35$hfunc2,c23g5)
c42g35 = BiCopSelect(h4g35$hfunc2,h2g35$hfunc2)
h4g235 = BiCopHfunc(h4g35$hfunc2,h2g35$hfunc2,c42g35)
c82g35 = BiCopSelect(h8g35$hfunc2,h2g35$hfunc2)
h8g235 = BiCopHfunc(h8g35$hfunc2,h2g35$hfunc2,c82g35)
c84g235 = BiCopSelect(h8g235$hfunc2,h4g235$hfunc2)

#

h8g4235 = BiCopHfunc(h8g235$hfunc2,h4g235$hfunc2,c84g235)
c13g5 = BiCopSelect(h15$hfunc2,h35$hfunc2)
h1g35 = BiCopHfunc(h15$hfunc2,h35$hfunc2,c13g5)
c12g35 = BiCopSelect(h1g35$hfunc2, h2g35$hfunc2)
h1g235 = BiCopHfunc(h1g35$hfunc2,h2g35$hfunc2,c12g35)
c14g235 = BiCopSelect(h1g235$hfunc2,h4g235$hfunc2)
h1g4235 = BiCopHfunc(h1g235$hfunc2,h4g235$hfunc2,c14g235)
c81g4235 = BiCopSelect(h8g4235$hfunc2,h1g4235$hfunc2)

## Observe Params

c15 
c32 
c71
c42 
c85
c16g5
c35g2
c75g1
c46g2
c83g5
c76g51
c41g62
c82g35
c84g235
c81g4235

# LL 

LL15 = sum(log(BiCopPDF(u1,u5,c15)))
LL32 = sum(log(BiCopPDF(u3,u2,c32))) 
LL71 = sum(log(BiCopPDF(u7,u1,c71)))
LL42 = sum(log(BiCopPDF(u4,u2,c42)))
LL85 = sum(log(BiCopPDF(u8,u5,c85)))
LL16g5 = sum(log(BiCopPDF(h15$hfunc2,h65$hfunc2,c16g5)))
LL35g2 = sum(log(BiCopPDF(h32$hfunc2,h52$hfunc2,c35g2)))
LL75g1 = sum(log(BiCopPDF(h71$hfunc2,h51$hfunc2,c75g1)))
LL46g2 = sum(log(BiCopPDF(h42$hfunc2,h62$hfunc2,c46g2)))
LL83g5 = sum(log(BiCopPDF(h85$hfunc2,h35$hfunc2,c83g5)))
LL76g51 = sum(log(BiCopPDF(h7g51$hfunc2,h6g51$hfunc2,c76g51)))
LL41g62 = sum(log(BiCopPDF(h4g62$hfunc2,h1g62$hfunc2,c41g62)))
LL82g35 = sum(log(BiCopPDF(h8g35$hfunc2,h2g35$hfunc2,c82g35)))
LL84g235 = sum(log(BiCopPDF(h8g235$hfunc2,h4g235$hfunc2,c84g235)))
LL81g4235 = sum(log(BiCopPDF(h8g4235$hfunc2,h1g4235$hfunc2,c81g4235)))

LL15
LL32
LL71 
LL42 
LL85
LL16g5 
LL35g2 
LL75g1 
LL46g2
LL83g5 
LL76g51 
LL41g62 
LL82g35 
LL84g235 
LL81g4235 

## Gaussian Copulas

c15g = BiCopSelect(u1,u5,familyset = c(1))
c32g = BiCopSelect(u3,u2,familyset = c(1))
c71g = BiCopSelect(u7,u1,familyset = c(1))
c42g = BiCopSelect(u4,u2,familyset = c(1))
c85g = BiCopSelect(u8,u5,familyset = c(1))
c16g5g = BiCopSelect(h15$hfunc2,h65$hfunc2,familyset = c(1))
c35g2g = BiCopSelect(h32$hfunc2,h52$hfunc2,familyset = c(1))
c75g1g = BiCopSelect(h71$hfunc2,h51$hfunc2,familyset = c(1))
c46g2g = BiCopSelect(h42$hfunc2,h62$hfunc2,familyset = c(1))
c83g5g = BiCopSelect(h85$hfunc2,h35$hfunc2,familyset = c(1))
c76g51g = BiCopSelect(h7g51$hfunc2,h6g51$hfunc2,familyset = c(1))
c41g62g = BiCopSelect(h4g62$hfunc2,h1g62$hfunc2,familyset = c(1))
c82g35g = BiCopSelect(h8g35$hfunc2,h2g35$hfunc2,familyset = c(1))
c84g235g = BiCopSelect(h8g235$hfunc2,h4g235$hfunc2,familyset = c(1))
c81g4235g = BiCopSelect(h8g4235$hfunc2,h1g4235$hfunc2,familyset = c(1))

LL15g = sum(log(BiCopPDF(u1,u5,c15g)))
LL32g = sum(log(BiCopPDF(u3,u2,c32g))) 
LL71g = sum(log(BiCopPDF(u7,u1,c71g)))
LL42g = sum(log(BiCopPDF(u4,u2,c42g)))
LL85g = sum(log(BiCopPDF(u8,u5,c85g)))
LL16g5g = sum(log(BiCopPDF(h15$hfunc2,h65$hfunc2,c16g5g)))
LL35g2g = sum(log(BiCopPDF(h32$hfunc2,h52$hfunc2,c35g2g)))
LL75g1g = sum(log(BiCopPDF(h71$hfunc2,h51$hfunc2,c75g1g)))
LL46g2g = sum(log(BiCopPDF(h42$hfunc2,h62$hfunc2,c46g2g)))
LL83g5g = sum(log(BiCopPDF(h85$hfunc2,h35$hfunc2,c83g5g)))
LL76g51g = sum(log(BiCopPDF(h7g51$hfunc2,h6g51$hfunc2,c76g51g)))
LL41g62g = sum(log(BiCopPDF(h4g62$hfunc2,h1g62$hfunc2,c41g62g)))
LL82g35g = sum(log(BiCopPDF(h8g35$hfunc2,h2g35$hfunc2,c82g35g)))
LL84g235g = sum(log(BiCopPDF(h8g235$hfunc2,h4g235$hfunc2,c84g235g)))
LL81g4235g = sum(log(BiCopPDF(h8g4235$hfunc2,h1g4235$hfunc2,c81g4235g)))

LL15g
LL32g
LL71g
LL42g
LL85g
LL16g5g 
LL35g2g 
LL75g1g 
LL46g2g
LL83g5g 
LL76g51g 
LL41g62g 
LL82g35g 
LL84g235g 
LL81g4235g 


