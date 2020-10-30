import math
import numpy as np
import matplotlib.pyplot as plt
from proj_1_module import dynamics_solve, hamiltonian_solve

def harmonic_oscillator_d_hamiltonian(m, w):
	"""
	Returns functions for the partial derivatives of the harmonic oscillator hamiltonian.

	Args:
		m (float): Mass
		w (float): Angular Frequency

	Returns:
		function, function: Function for dH/dx, and dH/dp
	"""
	return lambda x, p: m * w**2 * x, lambda x, p: p / m

def harmonic_oscillator_exact_sol(m, w):
	"""
	Return exact 1D simple harmonic oscillator functions.

	Args:
		m (float): Mass
		w (float): Angular Frequency

	Returns:
		function, function: Returns exact solution for x and p with arguments (x_0, p_0, t_0, t)
	"""
	x = lambda x_0, p_0, t_0, t: x_0 * math.cos(w * (t - t_0)) + (p_0 / (m * w)) * math.sin(w * (t - t_0))
	p = lambda x_0, p_0, t_0, t: p_0 * math.cos(w * (t - t_0)) - (m * w * x_0) * math.sin(w * (t - t_0))
	return x, p

def hamiltonian_exact_solve(x, p, d = 1, t_0 = 0.0, q_0 = 0.0, p_0 = 1.0, h = 0.1, N = 100):
	"""
	Produce time and value array for a given exact solution function x and p.

	Args:
		x (function): x(q_0, p_0, t_0, t) that returns the exact position at a given time and initial conditions
		p (function): p(q_0, p_0, t_0, t) that returns the exact momentum at a given time and initial conditions

	Kwargs:
		d: Spatial dimension (int) set to 1 as default
		t_0: Initial time (float) set to 0.0 as default
		q_0: Initial position (float for d=1, ndarray for d>1) set to 0.0 as default
		p_0: Initial momentum (float for d=1, ndarray for d>1) set to 1.0 as default
		h: Step size (float) set to 0.1 as default
		N: Number of steps (int) set to 100 as default

	Returns:
		T: Numpy array of times
		Q: Numpy array of positions at the times given in T
		P: Numpy array of momenta at the times given in T
	"""
	T = np.array([t_0 + n * h for n in range(N + 1)])
	
	if d == 1:
		P = np.zeros(N + 1)
		Q = np.zeros(N + 1)
	
	if d > 1:
		P = np.zeros((N + 1, d))
		Q = np.zeros((N + 1, d))
	
	Q[0] = q_0
	P[0] = p_0

	for n in range(N + 1):
		Q[n] = x(q_0, p_0, t_0, T[n])
		P[n] = p(q_0, p_0, t_0, T[n])
	
	return T, Q, P

def hamiltonian_display_step_size_plot(method, step_size = [0.10, 0.08, 0.06, 0.04, 0.02, 0.01]):
	# Define Constants
	m, w = 1, 1
	t_total = 10 * 2 * math.pi
	colors = ["#0000FF", "#0040C0", "#008080", "#00C040", "#00FF00", "#00FF00"]

	# Create numpy arrays containing the exact solutions
	x_t, p_t = harmonic_oscillator_exact_sol(m, w)
	t_exact, x_exact, p_exact = hamiltonian_exact_solve(x_t, p_t, N = int(t_total / 0.1))

	# Define the Hamiltonian partial derivatives model velocity function
	d_qH, d_pH = harmonic_oscillator_d_hamiltonian(m, w)

	# Create numpy arrays containing numerical solutions for different step sizes. 
	t_sol = {}
	q_sol = {}
	p_sol = {}
	for h in step_size:
		t, q, p = hamiltonian_solve(d_qH, d_pH, h = h, N = int(t_total / h), method = method)
		t_sol[h] = t
		q_sol[h] = q
		p_sol[h] = p

	# Plot all of the numerical model solutions
	for i, h in enumerate(step_size):
		plt.plot(t_sol[h], q_sol[h], color = colors[i])

	# Plot exact solution
	plt.plot(t_exact, x_exact, color = "red", linestyle = "--")

	# Define Labels
	plt.xlabel('Time')
	plt.ylabel('Position')
	plt.title('Position vs. Time for Different Step Sizes')

	# Display Plot
	plt.show()

def hamiltonian_display_plot(methods, step_size, t_total = None):
	# Define Constants
	m, w = 1, 1
	if t_total is None:
		t_total = 10 * 2 * math.pi
	if isinstance(step_size, float) or isinstance(step_size, int):
		step_size = [step_size for i in range(len(methods))]

	# Create numpy arrays containing the exact solutions
	x_t, p_t = harmonic_oscillator_exact_sol(m, w)
	t_exact, x_exact, p_exact = hamiltonian_exact_solve(x_t, p_t, N = int(t_total / 0.1))

	# Define the Hamiltonian partial derivatives model velocity function
	d_qH, d_pH = harmonic_oscillator_d_hamiltonian(m, w)
	
	# Create numpy arrays containing numerical solutions for different methods. 
	t_sol = {}
	q_sol = {}
	p_sol = {}
	for i, method in enumerate(methods):
		h = step_size[i]
		t, q, p = hamiltonian_solve(d_qH, d_pH, h = h, N = int(t_total / h), method = method)
		t_sol[method] = t
		q_sol[method] = q
		p_sol[method] = p

	# Plot all of the numerical model solutions
	for method in methods:
		plt.plot(t_sol[method], q_sol[method])

	# Plot exact solution
	plt.plot(t_exact, x_exact, color = "red", linestyle = "--")

	# Define Labels
	plt.xlabel('Time')
	plt.ylabel('Position')
	plt.title('Position vs. Time for Different Methods')

	# Display Plot
	plt.show()
	
