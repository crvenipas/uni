import math
from matplotlib import pyplot as plt

t=10.

Tnach=350. #темп перед нагревом
T0=500. #х=0
Tn=400. #х=Н

C0=381. #теплоемк
lam=384. #теплопров
ro=8800. #плотность

H=0.1 #метра
n=10 #по коорд
k=100 #по врем

h=H/(n-1) #шаг
tay=t/k

A=lam/(h*h)
B=2*lam/(h*h)+ro*C0/tay
C=lam/(h*h)
D=[0]*(n+1)


ttec=tay

Tf=[]
T=[T0,Tnach,Tnach,Tnach,Tn]
T1=[0,0,0,0,0]

s=0
Hh=H+h
T=T1

while ttec<=t:
    Tf.append(T1)
    for i in reversed(range(1,n)):
        T1[i]=T1[i] + ((T1[i+1]-2*T1[i]+T1[i-1])*tay*lam/(ro*C0*h*h))
    for i in range(1,n): 
      print(T1[i])
    ttec=tay+ttec
x=[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1]

fig, ax = plt.subplots(figsize=(5, 3))
ax.plot(x, T1, linestyle='solid')
plt.show()
