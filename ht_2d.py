from sympy import linear_eq_to_matrix, MatrixSymbol, Eq, solve_linear_system
from sympy import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

############## Variables Fijas ##############

Q = 15000. # Calor generado en cuadrante simulado
delta = 0.05 # espaciamiento entre nodos
k = 5. # Conductividad de la pieza

# Exterior
T_e = 25. # Temperatura exterior (aire)
u_e = 1.55e-5 # Viscocidad cinematica del aire
k_e = 0.024 # Conductividad del aire
r_e = 1.09 # Densidad del aire
Cp_e = 1005. # Capacidad calorífica del aire

# Refrigerante
T_f = 25. # Temperatura interior (refrigerante)
u_f = 0.295e-6 # Viscocidad cinematica del refrigerante
k_f = 0.58 # Conductividad del refrigerante
r_f = 997. # Densidad del refrigerante
Cp_f = 4181. # Capacidad calorífica del refrigerante

def solver(l,L,v_e,v_f,caso): #caso = True 'Sin aletas', False 'Con aletas'

	print('Caudal del refrigerante = {:.3f} m^3/s'.format(v_f*(0.5*L)**2), '\n')

	############## Calculo de h* ##############

	# Exterior
	alpha_e = k_e/(r_e*Cp_e) # Difusividad termica del aire
	Pr_e = u_e/alpha_e # Prandtl para el aire
	Re_e = v_e*l/u_e # Reynolds para el aire
	if Re_e < 10**5:
		Nu_e = 0.664*(Re_e**(1/2))*(Pr_e**(1/3))
	else:
		Nu_e = 0.037*((Re_e**(4/5))-871)*(Pr_e**(1/3))
	h_e = Nu_e*k_e/l 

	# Refrigerante
	D = L/2.
	alpha_f = k_f/(r_f*Cp_f) # Difusividad termica del refrigerante
	Pr_f = u_f/alpha_f # Prandtl para el refrigerante
	Re_f = v_f*D/u_f # Reynolds para el refrigerante
	if Re_f < 10**4:
		Nu_f = 4.36
	else:
		Nu_f = 0.023*(Re_f**(4/5))*(Pr_f**(1/3))
	h_f = Nu_f*k_f/D 

	# Aletas
	N = 6 # numero de aletas por un lado
	L_a = 0.1 # Longitud de las aletas
	t = L/(2*N) # Espesor de las aletas
	A_o = l*t # Area entre aletas
	A_b = l*t # Area de la base de la aleta
	A_f = 2*L_a*l + t*l # Area de la aleta
	P = 2*l+2*t # Perimetro de la aleta
	m = np.sqrt(h_e*P/(k*A_b)) # Parametro de la aleta
	L_c = L_a+(t/2) # Longitud corregida
	eta = np.tanh(m*L_c)/(m*L_c) # Eficiencia
	H_e = h_e*(A_o + eta*A_f)/(A_o + A_b)
	############## Adimensionalización ##############

	h_l = 0. # frontera aislada
	h_e = (h_e*L/k)*delta # frontera externa
	h_f = (h_f*L/k)*delta # frontera interna
	H_e = (H_e*L/k)*delta # aletas

	q = Q/(l*3*(0.25*L)**2) # Calor por unidad de volumen
	q = (q*L**2/(k*T_e))*(delta**2) # generación de calor adimensional

	if caso:
		print('h_f* = {:.3f}, h_e* = {:.3f}'.format(h_f/delta,h_e/delta), '\n')
	else:
		print('h_f* = {:.3f}, H_e* = {:.3f}'.format(h_f/delta,H_e/delta), '\n')
		h_e = H_e

	############## Solución numérica ##############

	n = int(1/delta) + 1 # numero de nodos
	t = symarray('t',(n,n)) # nodos
	equations = [] # array de ecuaciones
	counter = 0 # conteo para verificación

	for i in range(0,n):
		for j in range(0,n):
			# Fronteras
			if j == 0: # frontera vertical izq
				if i == 0: # punto superior
					equations.append(Eq(t[i+1,j] + t[i,j+1] - (t[i,j]*(2+h_e+h_l)),0))
					counter+=1

				elif i > 0 and i*delta < 0.25:
					equations.append(Eq(2*t[i,j+1] + t[i+1,j] + t[i-1,j] - 2*t[i,j]*(2+h_l),0))
					counter+=1

				elif i*delta >= 0.25 and i*delta <= 0.5:
					equations.append(Eq(2*t[i,j+1] + t[i+1,j] + t[i-1,j] - 2*t[i,j]*(2+h_l),-q))
					counter+=1

				elif i*delta > 0.5 and i*delta < 0.75:
					equations.append(Eq(2*t[i,j+1] + t[i+1,j] + t[i-1,j] - 2*t[i,j]*(2+h_l),0))
					counter+=1

				elif i*delta == 0.75: # frontera interna
					equations.append(Eq(t[i-1,j] + t[i,j+1] - (t[i,j]*(2+h_l+h_f)),0))
					counter+=1

				elif i*delta > 0.75: # punto inferior
					equations.append(Eq(t[i,j],0))
					counter+=1

			if j == n-1: # frontera vertical der
				if i == 0: # punto superior
					equations.append(Eq(t[i+1,j] + t[i,j-1] - (2*t[i,j]*(1+h_e)),0))
					counter+=1

				elif i > 0 and i < n-1: # demas
					equations.append(Eq(2*t[i,j-1] + t[i+1,j] + t[i-1,j] - 2*t[i,j]*(2+h_e),0))
					counter+=1

				elif i == n-1: # punto inferior
					equations.append(Eq(t[i-1,j] + t[i,j-1] - (t[i,j]*(2+h_e+h_l)),0))
					counter+=1

			if i == 0: # frontera horizontal sup
				if j > 0 and j < n-1: # demas
					equations.append(Eq(2*t[i+1,j] + t[i,j-1] + t[i,j+1] - 2*t[i,j]*(2+h_e),0))
					counter+=1

			if i == n-1: # frontera horizontal inf
				if j > 0 and j*delta < 0.25: # puntos en el fluido
					equations.append(Eq(t[i,j],0))
					counter+=1

				elif j*delta == 0.25: # frontera interna
					equations.append(Eq(t[i-1,j] + t[i,j+1] - (t[i,j]*(2+h_l+h_f)),0))
					counter+=1

				elif j*delta > 0.25 and j*delta < 0.5: # demas
					equations.append(Eq(2*t[i-1,j] + t[i,j-1] + t[i,j+1] - 2*t[i,j]*(2+h_l),0))
					counter+=1

				elif j*delta >= 0.5 and j*delta <= 0.75:
					equations.append(Eq(2*t[i-1,j] + t[i,j-1] + t[i,j+1] - 2*t[i,j]*(2+h_l),-q))
					counter+=1

				elif j*delta > 0.75 and j < n-1:
					equations.append(Eq(2*t[i-1,j] + t[i,j-1] + t[i,j+1] - 2*t[i,j]*(2+h_l),0))
					counter+=1

			#Frontera inversa
			if i*delta == 0.75 and j > 0 and j*delta < 0.25:
				equations.append(Eq(2*t[i-1,j] + t[i,j-1] + t[i,j+1] - 2*t[i,j]*(2+h_f),0))
				counter+=1

			if j*delta == 0.25 and i*delta > 0.75 and i < n-1:
				equations.append(Eq(2*t[i,j+1] + t[i+1,j] + t[i-1,j] - (2*t[i,j]*(2+h_f)),0))
				counter+=1

			if i*delta == 0.75 and j*delta == 0.25:
				equations.append(Eq(2*t[i-1,j] + 2*t[i,j+1] + t[i,j-1] + t[i+1,j] - (2*t[i,j]*(3+h_f)),0))
				counter+=1


			# Nodos internos
			if i > 0 and i*delta < 0.75 and j > 0 and j*delta <= 0.25: 
				if i*delta < 0.25 or i*delta > 0.5:
					equations.append(Eq(t[i+1,j] + t[i-1,j] + t[i,j+1] + t[i,j-1] - 4*t[i,j],0))
					counter+=1

				else:
					equations.append(Eq(t[i+1,j] + t[i-1,j] + t[i,j+1] + t[i,j-1] - 4*t[i,j],-q))
					counter+=1

			if i < n-1 and i*delta >= 0.75 and j < n-1 and j*delta > 0.25:
				if j*delta < 0.5 or j*delta > 0.75:
					equations.append(Eq(t[i+1,j] + t[i-1,j] + t[i,j+1] + t[i,j-1] - 4*t[i,j],0))
					counter+=1

				else:
					equations.append(Eq(t[i+1,j] + t[i-1,j] + t[i,j+1] + t[i,j-1] - 4*t[i,j],-q))
					counter+=1

			if i > 0 and i*delta < 0.75 and j*delta > 0.25 and j < n-1:
				if i*delta < 0.25 or i*delta > 0.5 or j*delta < 0.5 or j*delta > 0.75:
					equations.append(Eq(t[i+1,j] + t[i-1,j] + t[i,j+1] + t[i,j-1] - 4*t[i,j],0))
					counter+=1
				else:
					equations.append(Eq(t[i+1,j] + t[i-1,j] + t[i,j+1] + t[i,j-1] - 4*t[i,j],-q))
					counter+=1

			# Puntos restantes en el fluido
			if i*delta > 0.75 and i < n-1 and j*delta < 0.25 and j > 0:
				equations.append(Eq(t[i,j],0))
				counter+=1

	variables = []
	for i in range(0,n):
		for j in range(0,n):
			variables.append(t[i,j])

	a, b = linear_eq_to_matrix(equations, variables)

	a = np.array(a).astype(np.float64)
	b = np.array(b).astype(np.float64)

	t = np.around(np.linalg.solve(a,b), decimals = 3)
	t = np.ma.masked_where(t == 0, t)
	x = np.linspace(0,1,n)
	y = np.linspace(1,0,n)
	x,y = np.meshgrid(x,y)

	if caso:
		print('Sin aletas', '\n')
	else:
		print('Con aletas', '\n')

	print('Temperatura máxima en la pieza: {:.2f} ºC'.format((np.max(t)*T_e) + T_e))
	print('Temperatura máxima en la pared del conducto: {:.2f} ºC'.format((np.min(t)*T_e) + T_e + (4*Q/(r_f*((0.5*L)**2)*v_f*Cp_f))), '\n')


	fig = plt.figure(figsize = (10,5))
	ax1 = fig.add_subplot(1,2,1)
	ax1.imshow(t.reshape(n,n), cmap = cm.viridis, extent = [0,1,0,1], aspect = 'auto')
	ax1.set_xticks([0,0.25,0.5,0.75,1])
	ax1.set_yticks([0,0.25,0.5,0.75,1])
	ax1.set_title('Temperatura por nodo')
	ax1.set_xlabel('$x$ (m)')
	ax1.set_ylabel('$y$ (m)')

	ax2 = fig.add_subplot(1,2,2)
	CS = ax2.contour(x, y, t.reshape(n,n), 25)
	ax2.clabel(CS, inline=1, fontsize=6)
	ax2.set_xticks([0,0.25,0.5,0.75,1])
	ax2.set_yticks([0,0.25,0.5,0.75,1])
	ax2.set_title('Líneas de contorno')
	ax2.set_xlabel('$x$ (m)')
	ax2.set_ylabel('$y$ (m)')


	plt.tight_layout()
	plt.show()
	#if caso:
	#	plt.savefig('sin_aletas.png')
	#else:
	#	plt.savefig('con_aletas.png')

############## Variables para cambiar ##############

solver(2.7,0.3,5.,0.2,True)
solver(2.0,0.3,5.,0.2,False)
