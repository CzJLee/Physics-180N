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