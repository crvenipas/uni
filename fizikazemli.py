import matplotlib.pyplot as plt
plt.ion()
from scipy.constants import pi
import numpy as np
import pandas as pd 
import numpy as np 

G=6.67408e-11

po0=13000
po1=12700
po2=12140
po3=9900
po4=5550
po5=4377
po6=4071
po7=3768
po8=3553
po9=3305
po10=2850
po11=2850
po12=1500
po13=1500
po14=1030
po15=1030

r0=0
r1=1217000
r2=3485000
r3=5701000
r4=5951000
r5=6360000
r6=6366000
r7=6367000
r8=6371000

k0= (po1-po0)/(r1-r0)
k1= (po3-po2)/(r2-r1)
k2= (po5-po4)/(r3-r2)
k3= (po7-po6)/(r4-r3)
k4= (po9-po8)/(r5-r4)
k5= (po11-po10)/(r6-r5)
k6= (po13-po12)/(r7-r6)
k7= (po15-po14)/(r8-r7)
b0=po1-k0*r1
b1=po3-k1*r2
b2=po5-k2*r3
b3=po7-k3*r4
b4=po9-k4*r5
b5=po11-k5*r6
b6=po13-k6*r7
b7=po15-k7*r8




def f(x,k,b):
    return (4*pi*x*x*(k*x+b)*25000)
    
def I(x,k,b):
    return (8*pi*x*x*x*x*(k*x+b)*25000/(3*5.93e24*6370000*6370000))

M=0
h=0
Mm=np.array([])
gm=np.array([])

H=np.arange(0, 6400000, 25000)
 
for h in H:
 if (h<r1): 
  print(h)
  M=f(h,k0,b0)+M
  g=M*G/(h*h)
  gm = np.append(gm,g)
  Mm = np.append(Mm,M)
 
 if (h>=r1 and h<r2):
  print(h)
  M=f(h,k1,b1)+M
  g=M*G/(h*h)
  gm = np.append(gm,g)
  Mm = np.append(Mm,M)
 
 if (h>=r2 and h<r3):
  print(h)
  M=f(h,k2,b2)+M
  g=M*G/(h*h)
  gm = np.append(gm,g)
  Mm = np.append(Mm,M)
  
 if (h>=r3 and h<r4):
  print(h)
  M=f(h,k3,b3)+M
  g=M*G/(h*h)
  gm = np.append(gm,g)
  Mm = np.append(Mm,M)
 if (h>=r4 and h<r5):
  print(h)
  M=f(h,k4,b4)+M
  g=M*G/(h*h)
  gm = np.append(gm,g)
  Mm = np.append(Mm,M)
  
 if (h>=r5 and h<r6):
  print(h)
  M=f(h,k5,b5)+M
  g=M*G/(h*h)
  gm = np.append(gm,g)
  Mm = np.append(Mm,M)
  
 if (h>=r6 and h<r7):
  print(h)
  M=f(h,k6,b6)+M
  g=M*G/(h*h)
  gm = np.append(gm,g)
  Mm = np.append(Mm,M)
  
 if (h>=r7): 
  print(h)
  M=f(h,k7,b7)+M
  g=M*G/(h*h)
  gm = np.append(gm,g)
  Mm = np.append(Mm,M)

x=25000
H=np.flip(H,0)
gm=np.flip(gm,0)
p=0
pm=np.array([])
i=0
for i in range(0,len(H)):
 if (H[i]<r1): 
  p=(k0*H[i]+b0)*gm[i]*x+p
  pm = np.append(pm,p)
 
 if (H[i]>=r1 and H[i]<r2):
  p=(k1*H[i]+b1)*gm[i]*x+p
  pm = np.append(pm,p)
  
 if (H[i]>=r2 and H[i]<r3):
  p=(k2*H[i]+b2)*gm[i]*x+p
  pm = np.append(pm,p)
  
 if (H[i]>=r3 and H[i]<r4):
  p=(k3*H[i]+b3)*gm[i]*x+p
  pm = np.append(pm,p) 

 if (H[i]>=r4 and H[i]<r5):
  p=(k4*H[i]+b4)*gm[i]*x+p
  pm = np.append(pm,p)
  
 if (H[i]>=r5 and H[i]<r6):
  p=(k5*H[i]+b5)*gm[i]*x+p
  pm = np.append(pm,p)  

 if (H[i]>=r6 and H[i]<r7):
  p=(k6*H[i]+b6)*gm[i]*x+p
  pm = np.append(pm,p)
  
 if (H[i]>=r7): 
  p=(k7*H[i]+b7)*gm[i]*x+p
  pm = np.append(pm,p)
  
Ii=0
h=0
H=np.arange(0, 6400000, 25000)
Im=np.array([])
for h in H:
 if (h<r1): 
  print(h)
  Ii=I(h,k0,b0)+Ii
  Im = np.append(Im,Ii)
 
 if (h>=r1 and h<r2):
  print(h)
  Ii=I(h,k1,b1)+Ii
  Im = np.append(Im,Ii)

 if (h>=r2 and h<r3):
  print(h)
  Ii=I(h,k2,b2)+Ii
  Im = np.append(Im,Ii)

 if (h>=r3 and h<r4):
  print(h)
  Ii=I(h,k3,b3)+Ii
  Im = np.append(Im,Ii)

 if (h>=r4 and h<r5):
  print(h)
  Ii=I(h,k4,b4)+Ii
  Im = np.append(Im,Ii)

 if (h>=r5 and h<r6):
  print(h)
  Ii=I(h,k5,b5)+Ii
  Im = np.append(Im,Ii)

 if (h>=r6 and h<r7):
  print(h)
  Ii=I(h,k6,b6)+Ii
  Im = np.append(Im,Ii)
 
 if (h>=r7): 
  print(h)
  Ii=I(h,k7,b7)+Ii
  Im = np.append(Im,Ii)
  
  

