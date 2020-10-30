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
	Produce time and value array for a given exact solution function(s_0, t_0, t).

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

def dynamics_display_step_size_plot(method):
	# Define Constants
	b = 1.5
	d = 1
	p_0 = 1000
	t_total = 10
	step_size = [1, 0.5, 0.1, 0.05, 0.01]
	colors = ["#0000FF", "#0040C0", "#008080", "#00C040", "#00FF00", "#00FF00"]

	# Create numpy arrays containing the exact solutions
	p_exact = population_model_exact_sol(b, d)
	t_exact, s_exact = dynamics_exact_solve(p_exact, s_0 = p_0)

	# Define the population model velocity function
	p_model = population_model_veloc_function(b, d)

	# Create numpy arrays containing numerical solutions for different step sizes. 
	t_euler = {}
	s_euler = {}
	for h in step_size:
		t, s = dynamics_solve(p_model, s_0 = p_0, h = h, N = int(t_total/h), method = method)
		t_euler[h] = t
		s_euler[h] = s

	# Plot all of the numerical model solutions
	for i, h in enumerate(step_size):
		plt.plot(t_euler[h], s_euler[h], color = colors[i])

	# Plot exact solution
	plt.plot(t_exact, s_exact, color = "red", linestyle = "--")

	# Define Labels
	plt.xlabel('Time')
	plt.ylabel('Population')
	plt.title('Population Model Solutions for Different Step Sizes')

	# Display Plot
	plt.show()

def dynamics_display_plot_error(method):
	# Define Constants
	b = 1.5
	d = 1
	p_0 = 1000
	t_total = 100
	h = 0.01

	# Create numpy arrays containing the exact solutions
	p_exact = population_model_exact_sol(b, d)
	t_exact, s_exact = dynamics_exact_solve(p_exact, s_0 = p_0, h = h, N = int(t_total/h))

	# Define the population model velocity function
	p_model = population_model_veloc_function(b, d)

	# Create numpy arrays containing numerical solutions for different methods. 
	t_euler, s_euler = dynamics_solve(p_model, s_0 = p_0, h = h, N = int(t_total/h), method = method)

	# Plot all of the model solutions
	plt.plot(t_euler, s_euler, color = "blue")
	plt.plot(t_exact, s_exact, color = "red", linestyle = "--")

	# Define Labels
	plt.xlabel('Time')
	plt.ylabel('Population')
	plt.title('Numerical vs. Exact Solutions for Long Time Period')

	# Display Plot
	plt.show()

	# Calculate percent error after t = 100
	print(f"{round(100 * (s_exact[-1] - s_euler[-1]) / s_exact[-1], 2)}% Error after t = 100.")