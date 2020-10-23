import math
import numpy as np
import matplotlib.pyplot as plt
from proj_1_module import dynamics_solve, hamiltonian_solve

def population_model_exact_sol(b, d):
	"""
	Return exact population model solution for the population at a given time.

	Args:
		b (float): Birth Rate
		d (float): Death Rate

	Returns:
		function: function(s_0, t_0, t) that returns the population at time t, given the initial conditions s_0 and t_0.
	"""
	return lambda s_0, t_0, t: s_0 * math.exp((b - d) * (t - t_0))

def population_model_veloc_function(b, d):
	"""
	Population model velocity function. 

	Args:
		b (float): Birth Rate
		d (float): Death Rate

	Returns:
		function: function(t, s) that returns the result of the first order population value given t and s. 
	"""
	return lambda t, s: (b - d) * s

def dynamics_exact_solve(f, d = 1, t_0 = 0.0, s_0 = 1.0, h = 0.1, N = 100):
	"""
	Produce time and value array for a given exact solution function(s_0, t_0, t)

	Args:
		f (function): function(s_0, t_0, t) to produce exact values for

	Kwargs:
		d: Phase space dimension (int) set to 1 as default
		t_0: Initial time (float) set to 0.0 as default
		s_0: Initial state (float for D=1, ndarray for D>1) set to 1.0 as default
		h: Step size (float) set to 0.1 as default
		N: Number of steps (int) set to 100 as default

	Returns:
		np.array, np.array: Returns numpy array of times, and states at the given times
	"""
	T = np.array([t_0 + n * h for n in range(N + 1)])
	
	if d == 1:
		S = np.zeros(N + 1)
	
	if d > 1:
		S = np.zeros((N + 1, d))
		
	S[0] = s_0

	for n in range(N + 1):
		S[n] = f(s_0, t_0, T[n])
	
	return T, S

b = 2
d = 1
step_size = 0.1
n = 100

p_exact = population_model_exact_sol(b, d)
t_exact, s_exact = dynamics_exact_solve(p_exact, h = step_size, N = n)

p_model = population_model_veloc_function(b, d)
t_1, s_1 = dynamics_solve(p_model, h = step_size, N = n, method = "Euler")

t_2, s_2 = dynamics_solve(p_model, h = step_size, N = n, method = "RK2")

t_4, s_4 = dynamics_solve(p_model, h = step_size, N = n, method = "RK4")

plt.plot(t_1, s_1, color = "blue")
plt.plot(t_2, s_2, color = "green")
plt.plot(t_4, s_4, color = "orange")
plt.plot(t_exact, s_exact, color = "red", linestyle = "--")
plt.show()