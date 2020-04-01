from sympy import linear_eq_to_matrix, symbols, Eq, solve_linear_system
from sympy import *
import numpy as np

t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9, t_10, t_11, t_12 = symbols('t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9, t_10, t_11, t_12')
delta = 1/3. 
h = 1.*delta
q = 1.*(delta**2)


eq1 = Eq(t_2 + t_5 - (2*t_1*(1+h)),0)
eq2 = Eq(2*t_6 + t_1 + t_3 - (2*t_2*(2+h)),0)
eq3 = Eq(2*t_7 + t_4 + t_2 - (2*t_3*(2+h)),0)
eq4 = Eq(t_8 + t_3 - (2*t_4*(1+h)),0)
eq5 = Eq(2*t_6 + t_1 + t_9 - (2*t_5*(2+h)),0)
eq6 = Eq(t_2 + t_10 + t_7 + t_5 - 4*t_6 + q,0)
eq7 = Eq(t_3 + t_11 + t_6 + t_8 - 4*t_7 + q,0)
eq8 = Eq(2*t_7 + t_4 + t_12 - 4*t_8,0)
eq9 = Eq(t_5 + t_10  - (t_9*2),0)
eq10 = Eq(2*t_6 + t_9 + t_11 - 4*t_10,0)
eq11 = Eq(2*t_7 + t_10 + t_12 - 4*t_11,0)
eq12 = Eq(t_8 + t_11 - 2*t_12,0)

variables = [t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9, t_10, t_11, t_12]
equations = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12]

a, b = linear_eq_to_matrix(equations, variables)

a = np.array(a).astype(np.float64)
b = np.array(b).astype(np.float64)

t = np.around(np.linalg.solve(a,b), decimals = 3)

print(t[:4].tolist())
print(t[4:8].tolist())
print(t[8:].tolist())