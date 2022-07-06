import math
from numpy import *
import matplotlib.pyplot as plt
import numpy as np
c=2.99792e10
h=6.626e-27
e=4.802e-10
m=9.109e-28
k=1.38e-16
t1 = 6000
t2 = 10000
t3 = 30000
i = 2
ni1 = c / 1400e-7
ni2 = c / 820e-7
ni3 = c / 364.6e-7


ne=1e15
npl=[1,1/1.6,1/2.5] 
hi1 = (2 * math.pi**2 * m * e**4) / (h*h)

def alpha(T, n, ni,npl):
    a = ne*npl*((2**4 * math.pi**2 * e**6 * k * T) / (10*n**3 * c * h * (6 * math.pi * m * k * T)**(3/2))) * (1 + (2 * hi1) / (k * T) * (math.exp(h * ni / (k * T)) / (i*i*i)))*(1-(exp(-(h*n/(k*T)))))
    return (a)



n1 = linspace(ni1, ni2, 100)
a1 = alpha(t1, n1, ni1,npl[0])
a2 = alpha(t2, n1, ni1,npl[0])
a3 = alpha(t3, n1, ni1,npl[0])
plt.plot(n1, a1/10, 'r')
plt.plot(n1, a2*1e2, 'g')
plt.plot(n1, a3*0.7e6, 'b')

n2 = linspace(ni2, ni3, 100)
b1 = alpha(t1, n2, ni2,npl[1])
b2 = alpha(t2, n2, ni2,npl[1])
b3 = alpha(t3, n2, ni2,npl[1])
plt.plot(n2, b1/10, 'r')
plt.plot(n2, b2*1e2, 'g')
plt.plot(n2, b3*1e6, 'b')

n3 = linspace(ni3, c/91e-7, 100)
c1 = alpha(t1, n3, ni3,npl[2])
c2 = alpha(t2, n3, ni3,npl[2])
c3 = alpha(t3, n3, ni3,npl[2])
plt.plot(n3, c1/10, 'r')
plt.plot(n3, c2*1e2, 'g')
plt.plot(n3, c3*1e6, 'b')


n4 = linspace( c/91e-7,  c/70e-7, 100)

d1 = alpha(t1, n4, c/200e-7,npl[2])
plt.plot(n4, d1*1e3, 'r')


plt.axvline(x = ni1, color = 'y')
plt.axvline(x = ni2, color = 'y')
plt.axvline(x = c/91e-7, color = 'y')
plt.axvline(x = ni3, color = 'y')
plt.text(ni1, 0, 'Ha', fontsize = 15)
plt.text(ni2, 0, 'Hb', fontsize = 15)
plt.text(ni3, 0, 'Hg', fontsize = 15)
plt.legend(['6000K','10000K','30000K'])
plt.xlabel('Частота,Гц')
plt.ylabel('Alpha')
plt.semilogy()




plt.show()




def good(nu,alfa):
	Teff=5778
	m=h*nu/(k*Teff)
	alfa_sr=np.mean(alfa)
	n=alfa/alfa_sr
	f = m*2*h*nu*nu*nu*(exp(np.power(2*n,1/4))/(m*exp(np.power(2*n,1/4)) - m) - (exp(np.power(2*n,1/2)/np.power(8,1/4)) / (m*exp(np.power(2*n,1/2)/np.power(8,1/4)+2*m)- m*exp(np.power(2*n,1/2)))))/(c*c)
	
	return f
	
Teff=5778
ni1 = c / 1400e-7
ni2 = c / 820e-7
ni3 = c / 364.6e-7
n1 = linspace(ni1, ni2, 100)
n2 = linspace(ni2, ni3, 100)
n3 = linspace(ni3, c/200e-7, 100)
#n4 = linspace( c/91e-7,  c/70e-7, 100)

a1 = alpha(Teff, n1, ni1,npl[0])
b1 = alpha(Teff, n2, ni2,npl[1])
c1 = alpha(Teff, n3, ni3,npl[2])
#d1 = alpha(t1, n4, c/200e-7,npl[2])

I1=good(n1,a1)
plt.plot(n1, I1, 'r')
I2=good(n2,b1)
plt.plot(n2, I2, 'r')
I3=good(n3,c1)
plt.plot(n3, I3, 'r')
#I4=good(n4,d1)
#plt.plot(n4, I4, 'r')

#plt.semilogy()
plt.axvline(x = ni1, color = 'y')
plt.axvline(x = ni2, color = 'y')
#plt.axvline(x = c/91e-7, color = 'y')
plt.axvline(x = ni3, color = 'y')
plt.xlabel('Частота,Гц')
plt.ylabel('I, эрг см^-2 c^-1')
plt.show()



plt.close()
