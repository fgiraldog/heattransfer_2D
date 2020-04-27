from sympy import linear_eq_to_matrix, MatrixSymbol, Eq, solve_linear_system
from sympy import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

h_l = 0.
delta = 0.05
h = 1.*delta
h_f = 1.*delta
q = 1.*(delta**2)
n = int(1/delta) + 1 # numero de nodos

t = symarray('t',(n,n)) # nodos
equations = []
counter = 0

for i in range(0,n):
	for j in range(0,n):
		# Fronteras
		if j == 0: # frontera vertical izq
			if i == 0: # punto superior
				equations.append(Eq(t[i+1,j] + t[i,j+1] - (t[i,j]*(2+h+h_l)),0))
				counter+=1

			elif i*delta > 0.75: # punto inferior
				equations.append(Eq(t[i,j],0))
				counter+=1

			elif i*delta == 0.75:
				equations.append(Eq(t[i-1,j] + t[i,j+1] - (t[i,j]*(2+h+h_f)),0))
				counter+=1

			else: # demas
				equations.append(Eq(2*t[i,j+1] + t[i+1,j] + t[i-1,j] - 2*t[i,j]*(2+h_l),0)) #REVISAR
				counter+=1

		if j == n-1: # frontera vertical der
			if i == 0: # punto superior
				equations.append(Eq(t[i+1,j] + t[i,j-1] - (2*t[i,j]*(1+h)),0))
				counter+=1

			elif i == n-1: # punto inferior
				equations.append(Eq(t[i-1,j] + t[i,j-1] - (t[i,j]*(2+h+h_l)),0))
				counter+=1

			else: # demas
				equations.append(Eq(2*t[i,j-1] + t[i+1,j] + t[i-1,j] - 2*t[i,j]*(2+h),0))
				counter+=1

		if i == 0: # frontera horizontal sup
			if j > 0 and j < n-1: # demas
				equations.append(Eq(2*t[i+1,j] + t[i,j-1] + t[i,j+1] - 2*t[i,j]*(2+h),0))
				counter+=1

		if i == n-1: # frontera horizontal inf
			if j*delta < 0.25 and j*delta > 0:
				equations.append(Eq(t[i,j],0))
				counter+=1

			elif j*delta > 0.25 and j < n-1: # demas
				equations.append(Eq(2*t[i-1,j] + t[i,j-1] + t[i,j+1] - 2*t[i,j]*(2+h_l),0)) #REVISAR
				counter+=1

			elif j*delta == 0.25:
				equations.append(Eq(t[i-1,j] + t[i,j+1] - (t[i,j]*(2+h_l+h_f)),0))
				counter+=1


		#Frontera rara
		if i*delta == 0.75 and j > 0 and j*delta < 0.25:
			equations.append(Eq(2*t[i-1,j] + t[i,j-1] + t[i,j+1] - 2*t[i,j]*(2+h_f),0))
			counter+=1

		if j*delta == 0.25 and i < n-1 and i*delta > 0.75:
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

fig = plt.figure(figsize = (16,5))
ax1 = fig.add_subplot(1,3,1)
ax1.imshow(t.reshape(n,n), cmap = cm.viridis, extent = [0,1,0,1], aspect = 'auto')
ax1.set_xticks([0,0.25,0.5,0.75,1])
ax1.set_yticks([0,0.25,0.5,0.75,1])
ax1.set_title('Temperatura por nodo')
ax1.set_xlabel('$x$ (m)')
ax1.set_ylabel('$y$ (m)')

ax2 = fig.add_subplot(1,3,2)
CS = ax2.contour(x, y, t.reshape(n,n), 25)
ax2.clabel(CS, inline=1, fontsize=6)
ax2.set_xticks([0,0.25,0.5,0.75,1])
ax2.set_yticks([0,0.25,0.5,0.75,1])
ax2.set_title('Líneas de contorno')
ax2.set_xlabel('$x$ (m)')
ax2.set_ylabel('$y$ (m)')

ax3 = fig.add_subplot(1,3,3, projection='3d', aspect = 'auto')
surf = ax3.contourf(x, y, t.reshape(n,n), levels = 25, cmap=cm.viridis)
ax3.set_xticks([0,0.25,0.5,0.75,1])
ax3.set_yticks([0,0.25,0.5,0.75,1])
ax3.set_title('Proyección en 3D')
ax3.set_xlabel('$x$ (m)')
ax3.set_ylabel('$y$ (m)')

plt.tight_layout()
plt.show()
