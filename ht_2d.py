from sympy import linear_eq_to_matrix, MatrixSymbol, Eq, solve_linear_system
from sympy import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

q = 1.
h = 1.
t_f = 5.
h_f = 1.
h_l = 0.
delta = 0.05
h = 1.*delta
n = int(1/delta) + 1 # numero de nodos

t = symarray('t',(n,n)) # nodos

equations = []

counter = 0
for i in range(0,n):
	for j in range(0,n):

		# Fronteras
		if j == 0: # frontera vertical izq
			if i == 0: # punto superior
				equations.append(Eq(t[i+1,j] + t[i,j+1] - (2*t[i,j]*(1+h_l)),0))
				counter+=1

			elif i*delta > 0.75: # punto inferior
				equations.append(Eq(t[i,j],t_f))
				counter+=1

			elif i*delta == 0.75:
				equations.append(Eq(t[i-1,j] + t[i,j+1] - (t[i,j]*(2+h+h_f)),0))
				counter+=1

			else: # demas
				equations.append(Eq(2*t[i,j+1] + t[i+1,j] + t[i-1,j] - 2*t[i,j]*(2+h_l),0))
				counter+=1

		if j == n-1: # frontera vertical der
			if i == 0: # punto superior
				equations.append(Eq(t[i+1,j] + t[i,j-1] - (2*t[i,j]*(1+h)),0))
				counter+=1

			elif i == n-1: # punto inferior
				equations.append(Eq(t[i-1,j] + t[i,j-1] - (2*t[i,j]*(1+h)),0))
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
				equations.append(Eq(t[i,j],t_f))
				counter+=1

			elif j*delta > 0.25 and j < n-1: # demas
				equations.append(Eq(2*t[i-1,j] + t[i,j-1] + t[i,j+1] - 2*t[i,j]*(2+h_l),0))
				counter+=1

			elif j*delta == 0.25:
				equations.append(Eq(t[i-1,j] + t[i,j+1] - (t[i,j]*(2+h+h_f)),0))
				counter+=1


		#Frontera rara
		if i*delta == 0.75 and j > 0 and j*delta < 0.25:
			equations.append(Eq(2*t[i-1,j] + t[i,j-1] + t[i,j+1] - 2*t[i,j]*(2+h_f),0))
			counter+=1

		if j*delta == 0.25 and i < n-1 and i*delta > 0.75:
			equations.append(Eq(t[i+1,j] + t[i,j+1] - (2*t[i,j]*(1+h_f)),0))
			counter+=1

		if i*delta == 0.75 and j*delta == 0.25:
			equations.append(Eq(t[i-1,j] + t[i,j+1] - (2*t[i,j]*(1+h)),0))
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
			equations.append(Eq(t[i,j],t_f))
			counter+=1


print(counter)
variables = []
for i in range(0,n):
	for j in range(0,n):
		variables.append(t[i,j])

a, b = linear_eq_to_matrix(equations, variables)

print(a.shape)
a = np.array(a).astype(np.float64)
b = np.array(b).astype(np.float64)

t = np.around(np.linalg.solve(a,b), decimals = 3)

plt.figure()
plt.imshow(t.reshape(n,n), cmap = cm.viridis)
clb = plt.colorbar()
clb.ax.set_title('$^{\circ}$ C')
plt.axis('off')
plt.show()

fig = plt.figure()
ax = fig.gca(projection='3d')
x = np.linspace(0,1,n)
y = np.linspace(0,1,n)
x,y = np.meshgrid(x,y)
surf = ax.plot_surface(x, y, t.reshape(n,n), cmap=cm.viridis, linewidth=0, antialiased=False)
clb = plt.colorbar(surf)
clb.ax.set_title('$^{\circ}$ C')
plt.show()