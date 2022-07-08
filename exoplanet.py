import numpy as np
from scipy. integrate import odeint
from scipy. integrate import quad
import matplotlib.pyplot as plt
import math
from scipy import integrate
from astropy.io import fits
import pysynphot as S

import rebound
Rj=69911000
R=np.array([0.096,0.103,0.09,0.115])*Rj
ro=5515
a=np.array([0.04,0.078,0.05,0.106])*1.496e8
T=np.array([5.577212,15.189,7.76,24.15])
e=np.array([0.11,0.11,0.11,0.11])
ic=np.array([89.15,89.54,89.26,89.68])
omega=np.array([7.49,16.83,14.28,11.39])

TIMES = np.array([2457010.376,2456989.465,2456984.788,2456987.054  ])
M = np.array([0,0,0])
M[1]= 360 + (TIMES[3]-TIMES[1])*360/T[1] +omega[3]
M[0]= 5*360 + (TIMES[3]-TIMES[0])*360/T[0] +omega[3]
M[2] = 360 - (TIMES[3]-TIMES[2])*360/T[2] -omega[3]

M_dens = 4*math.pi*R*R*R*ro*5.0279e-31/3

#M_per= 4*math.pi*math.pi*a*a*a/
OMEGA= np.array([0,90,180,270])

for g in range(0,1):
  for h in range(0,4):
    for n in range(0,3):
      sim = rebound.Simulation()
      sim.units = ('km','day','Msun')
      sim.add(m=0.21)
      sim.add(m=M_dens[0],a=a[0], e = e[0], inc=90-ic[0],omega=omega[0],Omega=OMEGA[g],M=M[0])
      sim.add(m=M_dens[1],a=a[1], e = e[1], inc=90-ic[1],omega=omega[1],Omega=OMEGA[h],M=M[1])
      sim.add(m=M_dens[2],a=a[2], e = e[2], inc=90-ic[2],omega=omega[2],Omega=OMEGA[n],M=M[2])
      sim.add(m=M_dens[3],a=a[3], e = e[3], inc=90-ic[3],omega=omega[3],Omega=0, M= 360-omega[3])

      sim.status()
#sim.units = ('AU','day','Msun')
      sim.init_megno()

      sim.integrator = "whfast"
      sim.dt = 1
      times = np.linspace(0, 365.25e6, 10000)
      
      e_w = np.zeros((4,10000))
      inc_w =np.zeros((4,10000))
      Omega_w = np.zeros((4,10000))
      omega_w  = np.zeros((4,10000))
      megno_w = np.zeros(10000)
      lyap_w = np.zeros(10000)
      
      
      for i,time in enumerate(times):
        sim.integrate(time, exact_finish_time=0)
        sim.calculate_orbits()
        lyap_w[i] = sim.calculate_lyapunov()
        e_w[0,i] = sim.particles[1].e
        e_w[1,i] = sim.particles[2].e
        e_w[2,i] = sim.particles[3].e
        e_w[3,i] = sim.particles[4].e
        inc_w[0, i] = sim.particles[1].inc
        inc_w[1, i] = sim.particles[2].inc
        inc_w[2, i] = sim.particles[3].inc
        inc_w[3, i] = sim.particles[4].inc
        Omega_w[0,i] = sim.particles[1].Omega
        Omega_w[1,i] = sim.particles[2].Omega
        Omega_w[2,i] = sim.particles[3].Omega
        Omega_w[3,i] = sim.particles[4].Omega
        omega_w[0,i] = sim.particles[1].omega
        omega_w[1,i] = sim.particles[2].omega
        omega_w[2,i] = sim.particles[3].omega
        omega_w[3,i] = sim.particles[4].omega
        megno_w[i] = sim.calculate_megno()
        print(time)
      if n<=3:
        plt.rcParams.update({'font.size': 12})
        fig, ((axs1,axs2),(axs3,axs4)) = plt.subplots(2, 2)
        axs1.set_xlabel("time")
        axs1.set_ylabel("i")
        axs1.plot(times, inc_w[0])
        axs1.plot(times, inc_w[1])
        axs1.plot(times, inc_w[2])
        axs1.plot(times, inc_w[3])
#        axs1.set_title('Юпитер')


        axs2.set_xlabel("time")
        axs2.set_ylabel("eccentricity")
        axs2.plot(times, e_w[0])
        axs2.plot(times, e_w[1])
        axs2.plot(times, e_w[2])
        axs2.plot(times, e_w[3])
#        axs2.set_title('Юпитер')

        axs3.set_xlabel("time")
        axs3.set_ylabel("omega")
        axs3.plot(times, omega_w[0])
        axs3.plot(times, omega_w[1])
        axs3.plot(times, omega_w[2])
        axs3.plot(times, omega_w[3])
#        axs3.set_title('Сатурн')


        axs4.set_xlabel("time")
        axs4.set_ylabel("OMEGA")
        axs4.plot(times,  Omega_w[0])
        axs4.plot(times, Omega_w[1])
        axs4.plot(times, Omega_w[2])
        axs4.plot(times, Omega_w[3])
#        axs4.set_title('Сатурн')
        plt.show()
      
      
