import math
import numpy as np
import matplotlib.pyplot as plt
from proj_1_module import dynamics_solve, hamiltonian_solve

def population_model_exact_sol(b, d):
	return lambda p0, t: p0 * math.exp((b - d) * t)

def population_model_veloc_function(b, d):
	return lambda t, s: (b - d) * s

def exact_solve(f, d = 1, t_0 = 0.0, s_0 = 1, h = 0.1, N = 100):
	T = np.array([t_0 + n * h for n in range(N + 1)])
	
	if d == 1:
		S = np.zeros(N + 1)
	
	if d > 1:
		S = np.zeros((N + 1, d))
		
	S[0] = s_0

	for n in range(N + 1):
		S[n] = f(s_0, T[n])
	
	return T, S

b = 2
d = 1
step_size = 0.1
n = 100

p_exact = population_model_exact_sol(b, d)
t_exact, s_exact = exact_solve(p_exact, h = step_size, N = n)

p_model = population_model_veloc_function(b, d)
t_1, s_1 = dynamics_solve(p_model, h = step_size, N = n, method = "Euler")

t_2, s_2 = dynamics_solve(p_model, h = step_size, N = n, method = "RK2")

t_4, s_4 = dynamics_solve(p_model, h = step_size, N = n, method = "RK4")

plt.plot(t_1, s_1, color = "blue")
plt.plot(t_2, s_2, color = "green")
plt.plot(t_4, s_4, color = "orange")
plt.plot(t_exact, s_exact, color = "red", linestyle = "--")
plt.show()