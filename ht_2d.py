from sympy import linear_eq_to_matrix, MatrixSymbol, Eq, solve_linear_system
from sympy import *
import numpy as np

q = 1.
h = 1.
t_f = 30.
delta = 0.05
h = 1.*delta
n = int(1/delta) + 1 # numero de nodos
t = symarray('t',(n,n)) # nodos

equations = []

for i in range(0,n):
	for j in range(0,n):

		# Fronteras
		if j == 0: # frontera vertical izq
			if i == 0: # punto superior
				equations.append(Eq(t[i+1,j] + t[i,j+1] - (2*t[i,j]*(1+h)),0))

			if i*delta > 0.75: # punto inferior
				equations.append(Eq(t[i,j],t_f))

			else: # demas
				equations.append(Eq(2*t[i,j+1] + t[i+1,j] + t[i-1,j] - 2*t[i,j]*(2+h),0))

		if j == n-1: # frontera vertical der
			if i == 0: # punto superior
				equations.append(Eq(t[i+1,j] + t[i,j-1] - (2*t[i,j]*(1+h)),0))

			if i == n-1: # punto inferior
				equations.append(Eq(t[i-1,j] + t[i,j-1] - (2*t[i,j]*(1+h)),0))

			else: # demas
				equations.append(Eq(2*t[i,j-1] + t[i+1,j] + t[i-1,j] - 2*t[i,j]*(2+h),0))

		if i == 0: # frontera horizontal sup
			if j > 0 and j < n-1: # demas
				equations.append(Eq(2*t[i+1,j] + t[i,j-1] + t[i,j+1] - 2*t[i,j]*(2+h),0))

		if i == n-1:
			if j*delta < 0.25:
				equations.append(Eq(t[i,j],t_f))

			if j*delta >= 0.25 and j < n-1: # demas
				equations.append(Eq(2*t[i-1,j] + t[i,j-1] + t[i,j+1] - 2*t[i,j]*(2+h),0))

		# Nodos internos
		if i < 0 and i*delta < 0.75 and j < 0 and j*delta <= 0.25: 
			if i*delta < 0.25 or i*delta > 0.5:
				equations.append(Eq(t[i+1,j] + t[i-1,j] + t[i,j+1] + t[i,j-1] - 4*t[i,j],0))

			else:
				equations.append(Eq(t[i+1,j] + t[i-1,j] + t[i,j+1] + t[i,j-1] - 4*t[i,j],-q))

		if i > n-1 and i*delta >= 0.75 and j > n-1 and j*delta > 0.25:
			if j*delta < 0.5 or j*delta > 0.75:
				equations.append(Eq(t[i+1,j] + t[i-1,j] + t[i,j+1] + t[i,j-1] - 4*t[i,j],0))

			else:
				equations.append(Eq(t[i+1,j] + t[i-1,j] + t[i,j+1] + t[i,j-1] - 4*t[i,j],-q))

		if i < 0 and i*delta < 0.75 and j*delta > 0.25 and j < n-1:
			if i*delta < 0.25 or i*delta > 0.5 or j*delta < 0.5 or j*delta > 0.75:
				equations.append(Eq(t[i+1,j] + t[i-1,j] + t[i,j+1] + t[i,j-1] - 4*t[i,j],0))
			else:
				equations.append(Eq(t[i+1,j] + t[i-1,j] + t[i,j+1] + t[i,j-1] - 4*t[i,j],-q))

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
plt.imshow(t.reshape(n,n))
plt.show()