
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from astropy import constants as const



#Якобиан 
def integ_jacobi (k,n,u,v):
	r=np.sqrt(k*k + n*n)
	return ((u*u + v*v)/2 - (G*m_2)/((mu*mu) *r) + (alpha*k*k)/2)



#ODE система
def func (t,y):
	k,n,u,v =y
	r = np.sqrt(k*k + n*n)
	dk_dt = u
	dn_dt = v
	du_dt = 2*w*v - ((G*m_2)/(mu*mu*r*r*r) + alpha) *k
	dv_dt = -2*w*u - (G*m_2)/(mu*mu*r*r*r)*n
	g = [dk_dt, dn_dt, du_dt, dv_dt]
	return (g)


# функция для пересчета начальной скорости
def n0 (C_1):
	return (np.sqrt(2*(eps_1_0 + (G*m_2)/(mu*mu * abs((-ksi_1t + C_1)))) - alpha*(-ksi_1t + C_1)*(-ksi_1t + C_1)))

def eps_control(eps_1):
	return (abs((eps_1-eps_1_0)/eps_1_0))
	


############################Константы!!!!!!!##############################
m_1 =1 # масса звезды
m_2= 499 # масса скопления
mu = 1 + m_1/m_2
R=8200
G= 0.004535 
alpha = -0.001976995 
w= 0.02829174 

ksi_1t = 10.44666776 
eps_1t = -(3*G*m_2)/(mu**2*ksi_1t*2) 


t0=0
eps_1_0 = -0.6*eps_1t
C_1 =-7.49875 #-6.7##-3.8317
step = -0.00001


####################начальные параметры#############################
eps_1 = eps_1_0 # нач энергия
t = t0 # нач время
y0 = [0, 0, 0, 0] # массив с начальными значениями
y0[0] = -ksi_1t + C_1
y0[1] = 0
y0[2] = 0
y0[3] = n0(C_1)
####################################################################
while C_1 >= -7.49876 : # Основной цикл, пробег по С1
	y0[0] = -ksi_1t + C_1
	y0[1] = 0
	y0[2] = 0
	y0[3] = n0(C_1)
	solution = integrate.RK45(func, t0, y0, t_bound=3000, max_step = 0.001)
	t_values = [] # массивы координаты, их будем заполнять
	k_values = []
	n_values = []
	u_values = []
	v_values = []
	eps_mas = [] # энергия для контроля ошибки
	while solution.t<=155: # Откидываем все кроме одного вращения (по картинке и времени)
		solution.step()
		t_values.append(solution.t)
		k_values.append(solution.y[0])
		n_values.append(solution.y[1])
		u_values.append(solution.y[2])
		v_values.append(solution.y[3])
		eps_1 = integ_jacobi(solution.y[0], solution.y[1], solution.y[2], solution.y[3])
		eps_mas.append(eps_1)
	print (eps_control(max(eps_mas)), C_1) #ошибка + С1
	fig, ax = plt.subplots() # строим орбиту
	ax.plot(k_values, n_values, color = 'red', linewidth = '2')
	ax.grid(color = 'grey', linewidth = '1', linestyle = '--')
	ax.set_xlabel ('xi')
	ax.set_ylabel('eta')
	fig.set_figwidth (5)
	plt. show() 
	qqq= (np.array(eps_mas)-eps_1_0)/eps_1_0
	fig, ax = plt.subplots() # строим орбиту
	ax.plot(t_values,np.absolute(qqq), color = 'red', linewidth = '2')
	ax.grid(color = 'grey', linewidth = '1', linestyle = '--')
	plt. show() 
	k = []
	n = []
	time = []
	for i in range (len(n_values)): #ищем период
		if i!=0:
			if n_values[i-1]*n_values[i]>0:
				k.append(k_values[i])
				n.append(n_values[i])
				time.append(t_values[i])
				p=2*i
				I=i
			else:
				half_period = max(time) # Период ееее
				print (half_period)
				print(i)
				break

	print (abs((k_values[p-2]+ksi_1t - C_1)/(-ksi_1t + C_1))) # Ищем круглость орбиты
	C_1 = C_1 + step # Следующее С1

#########################################################################################################


B = G * m_2 / mu**2

def sys (t, y):
	k, n, u, v, dk_dt, dn_dt, du_dt, dv_dt = y
	r = np.sqrt (k**2 + n**2)
	b = (G * m_2) / (mu**2 * r**3)
	dk_dt = y[2]
	dn_dt = y[3]
	du_dt = 2*w*y[3] - (b+alpha)*y[0]
	dv_dt = -2*w*y[2] - b*y[1]
	vk_dt = y[6]
	vn_dt = y[7]
	vu_dt = 2*w*y[7] - (b * (1 - 3*y[0]**2/r**2) + alpha)*y[4] + (3*b*y[0]*y[1]/r**2)*y[5]
	vv_dt = - 2*w*y[6] + (3*b*y[1]/r**2)*y[0]*y[4] - (b * (1 - 3*y[1]**2/r**2))*y[5]
	return (dk_dt, dn_dt, du_dt, dv_dt, vk_dt, vn_dt, vu_dt, vv_dt)
	
def const_count (R, y2, yб, y3, y7, y0, y4, y1, y5):
	return (y2*yб + y3*y7 + (alpha + B/R**3)*y0*y4 + (B/R**3)*y1*y5)

#Выставляем нач условия
C_1 = -7.49875
ksi_0 = -ksi_1t + C_1
in_cond_0 = [ksi_0, 0, 0, n0(C_1), 1, 0, 0, 0]
in_cond_1 = [ksi_0, 0, 0, n0(C_1), 0, 1, 0, 0]
in_cond_2 = [ksi_0, 0, 0, n0(C_1), 0, 0, 1, 0]
in_cond_3 = [ksi_0, 0, 0, n0(C_1), 0, 0, 0, 1]
in_cond = in_cond_0, in_cond_1, in_cond_2, in_cond_3
n = 0 
array = [[0,0,0,0], [0,0,0,0],[0,0,0,0],[0,0,0,0]]
# Находим матрицу
for i in in_cond:
	print (i)
	solution = integrate.RK45(sys, 0, i, t_bound=3000, max_step = 0.001)
	dk_dt, dn_dt, du_dt, dv_dt, vk_dt, vn_dt, vu_dt, vv_dt = [], [], [], [],[],[], [], []
	eps_mas, eps_delta, t, C, C_delta = [], [], [], [], []
	eps_mas.append(eps_1_0)
	eps_delta.append(eps_control(eps_mas[0]))
	r = np.sqrt (solution.y[0]**2 + solution.y[1]**2)
	const_0 = const_count(r, solution.y[2], solution.y[6], solution.y[3], solution.y[7], solution.y[0], solution.y[4],solution.y[1],solution.y[5])
	C_delta.append(0)
	C.append(const_0)
	t.append(0)

	while solution.t <= 155:
		dk_dt.append(solution.y[0])
		dn_dt.append(solution.y[1])
		du_dt.append(solution.y[2])
		dv_dt.append(solution.y[3])
		vk_dt.append(solution.y[4])
		vn_dt.append(solution.y[5])
		vu_dt.append(solution.y[6])
		vv_dt.append(solution.y[7])
		t. append(solution.t)
		r = np.sqrt(solution.y[0]**2 + solution.y[1]**2)
		Eps = integ_jacobi(solution.y[0], solution.y[1], solution.y[2], solution.y[3])
		eps_mas.append(Eps)
		eps_delta.append(eps_control(Eps))
		const = const_count(r, solution.y[2], solution.y[6], solution.y[3], solution.y[7], solution.y[0], solution.y[4],solution.y[1],solution.y[5])
		C.append(const)
		c_delt = abs((const-const_0)/const_0)
		C_delta.append(c_delt)
		solution. step()


	plt.plot(t, C_delta,color = 'black')
	plt.title ('C delta')
	plt.xlabel ('time')
	
	plt. show()
	
	
	plt. plot(t, C,color = 'black')
	plt.title ('C')
	plt.xlabel ('time')
	plt. show()
	
	array[0][n] = vk_dt[146250]
	array[1][n]= vn_dt[146250]
	array[2][n] = vu_dt[146250]
	array[3][n] = vv_dt[146250]
	print(array)
	print(vk_dt[146250], vn_dt[146250], vu_dt[146250], vv_dt[146250])
	n=n+1




#matrix
print ('matrix', array)
#Находим лямбда
e_vals = np.linalg.eig(array)

print ('Eigenvalues: ' , e_vals[0][0], e_vals[0][1], e_vals[0][2], e_vals[0][3])
#Модуль
L = [abs(e_vals[0][0]), abs(e_vals[0][1]), abs(e_vals [0][2]), abs(e_vals [0][3])]
print (L)



