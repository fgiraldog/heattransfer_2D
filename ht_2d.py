from sympy import linear_eq_to_matrix, MatrixSymbol, Eq, solve_linear_system
from sympy import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


Q = 15000. # Calor generado en cuadrante simulado
delta = 0.05 # espaciamiento entre nodos
l = 1. # Longitud de la pieza
L = 0.3 # Media longitud del lado del perfil de la pieza
k = 5. # Conductividad de la pieza

T_e = 25. # Temperatura exterior (aire)
v_e = 10. # Velocidad del aire
u_e = 1. # Viscocidad cinematica del aire
k_e = 1. # Conductividad del aire
r_e = 1. # Densidad del aire
Cp_e = 1. # Capacidad calorífica del aire
alpha_e = k_e/(r_e*Cp_e) # Difusividad termica del aire
Pr_e = u_e/alpha_e # Prandtl para el aire
Re_e = v_e*l/u_e # Reynolds para el aire
if Re_e < 10**5:
	Nu_e = 0.664*(Re_e**(1/2))*(Pr_e**(1/3))
else:
	Nu_e = 0.037*((Re_e**(4/5))-871)*(Pr_e**(1/3))
h_e = Nu_e*k_e/l 

D = L/2
T_f = 25. # Temperatura interior (refrigerante)
v_f = 10. # Velocidad del refrigerante
u_f = 1. # Viscocidad cinematica del refrigerante
k_f = 1. # Conductividad del refrigerante
r_f = 1. # Densidad del refrigerante
Cp_f = 1. # Capacidad calorífica del refrigerante
alpha_f = k_f/(r_f*Cp_f) # Difusividad termica del refrigerante
Pr_f = u_f/alpha_f # Prandtl para el refrigerante
Re_f = v_f*D/u_f # Reynolds para el refrigerante
if Re_f < 10**4:
	Nu_f = 4.36
else:
	Nu_f = 0.023*(Re_f**(4/5))*(Pr_f**(1/3))
h_f = Nu_f*k_f/D 

h_l = 0. # frontera aislada
h_e = (h_e*L/k)*delta # frontera externa
h_f = (h_f*L/k)*delta # frontera interna

q = Q/(l*3*(0.25*L)**2) # Calor por unidad de volumen
q = (q*L**2/(k*T_e))*(delta**2) # generación de calor adimensional

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

fig = plt.figure(figsize = (5,10))
ax1 = fig.add_subplot(2,1,1)
ax1.imshow(t.reshape(n,n), cmap = cm.viridis, extent = [0,1,0,1], aspect = 'auto')
ax1.set_xticks([0,0.25,0.5,0.75,1])
ax1.set_yticks([0,0.25,0.5,0.75,1])
ax1.set_title('Temperatura por nodo')
ax1.set_xlabel('$x$ (m)')
ax1.set_ylabel('$y$ (m)')

ax2 = fig.add_subplot(2,1,2)
CS = ax2.contour(x, y, t.reshape(n,n), 25)
ax2.clabel(CS, inline=1, fontsize=6)
ax2.set_xticks([0,0.25,0.5,0.75,1])
ax2.set_yticks([0,0.25,0.5,0.75,1])
ax2.set_title('Líneas de contorno')
ax2.set_xlabel('$x$ (m)')
ax2.set_ylabel('$y$ (m)')


plt.tight_layout()
plt.show()
