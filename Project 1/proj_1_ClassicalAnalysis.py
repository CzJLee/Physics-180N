import math
import numpy as np
import matplotlib.pyplot as plt
from proj_1_module import dynamics_solve, hamiltonian_solve

def kepler_d_hamiltonian(m_1, m_2):
	"""
	Returns functions for the partial derivatives of the kepler hamiltonian.

	Args:
		m_1 (float): Mass 1
		m_1 (float): Mass 1

	Returns:
		function, function: Function for dH/dx, and dH/dp. The inputs (x, p) should each be a numpy array.
	"""

	# Reduced mass 
	m = m_1 * m_2 / (m_1 + m_2)
	
	# Gravitational Constant
	g = 6.67408e-11
	k = g * m_1 * m_2

	return lambda x, p: np.array([k * x[0] / ((x[0] ** 2 + x[1] ** 2) ** 1.5), k * x[1] / ((x[0] ** 2 + x[1] ** 2) ** 1.5)]), lambda x, p: p / m, m, k

def plot_orbit(m_1, m_2, q_0, p_0, h = 86400, N = 365, method = "SE", display_plot = True, verbose = False):
	"""
	Return coordinate solutions for the given parameters, and plots the solution.

	Args:
		m_1 (float): Mass 1
		m_2 (float): Mass 2
		q_0 (numpy array): Initial Position
		p_0 (numpy array): Initial Momentum
		h (int, optional): Step Size. Defaults to 86400.
		N (int, optional): Number of steps. Defaults to 365.
		method (str, optional): Numerical solution method to use. Defaults to "SE".
		display_plot (bool, optional): Display the plot of the orbit. Defaults to True.

	Returns:
		(numpy array, numpy array): Return the x and y coordinate solutions.
	"""
	d_qH, d_pH, m, k = kepler_d_hamiltonian(m_1, m_2)

	# Solve system
	t, q, p = hamiltonian_solve(d_qH, d_pH, d = 2, t_0 = 0.0, q_0 = q_0, p_0 = p_0, h = h, N = N, method = method, verbose = verbose)

	# Unpack (x, y) coordinates to be plotted using matplotlib
	x, y = q.T
	px, py = p.T

	# Plot
	if display_plot:
		plt.plot(x, y)
		plt.show()

	# Return coordinate points
	return x, y, px, py

def find_major_axis(m_1, m_2, q_0, p_0, h = 86400, N = 365):
	"""
	Determine the period of orbit and semi major axis distance for a set of coordinate solutions.

	Args:
		m_1 (float): Mass 1
		m_2 (float): Mass 2
		q_0 (numpy array): Initial Position, Expect q_0 = (q_0, 0) where q_0 is positive.
		p_0 (numpy array): Initial Momentum, Expect p_0 = (0, p_0), where p_0 is positive.
		h (int, optional): Step Size. Defaults to 86400.
		N (int, optional): Number of steps. Defaults to 365.

	Returns:
		(float, float): Return the period as the number of steps, and the semi major axis distance. 
	"""
	# Expect q_0 = (q_0, 0), p_0 = (0, p_0), where q_0 and p_0 are both positive. If this is not the case, this function may return bogus. 
	x, y, *p = plot_orbit(m_1, m_2, q_0, p_0, h = h, N = N, method = "SE", display_plot = False)

	major_axis = None
	full_orbit = None
	for i in range(1, len(x)):
		if y[i] <= 0:
			# Coordinate has completed first half of orbit. 
			major_axis = i
			break
	for i in range(major_axis + 1, len(x)):
		if y[i] >= 0:
			# Coordinate has completed full orbit. 
			full_orbit = i
			break

	if major_axis is None:
		raise Exception("Did not complete half an orbit")
	elif full_orbit is None:
		raise Exception("Did not complete full orbit")
	
	# Calculate period, in seconds
	period = full_orbit * h

	# Calculate semi major axis, in meters (or what ever unit q_0 is in.)
	semi_major_length = (x[0] - x[major_axis]) / 2

	return period, semi_major_length