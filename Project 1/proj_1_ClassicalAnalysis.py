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
	# g = 0.0001
	k = g * m_1 * m_2

	return lambda x, p: np.array([k * x[0] / ((x[0] ** 2 + x[1] ** 2) ** 1.5), k * x[1] / ((x[0] ** 2 + x[1] ** 2) ** 1.5)]), lambda x, p: p / m, m, k

# Define using earth masses
m_sun = 1.989e30
m_earth = 5.972e24
d_qH, d_pH, m, k = kepler_d_hamiltonian(m_sun, m_earth)

# Set initial conditions
# Consider coordinates at Aphelion
# https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
q_0 = np.array([152.10e9, 0])
p_0 = np.array([0, 29290 * m])

# Verify elliptical motion
energy = (p_0[0] ** 2 + p_0[1] ** 2) / (2 * m) - k / ((q_0[0] ** 2 + q_0[1] ** 2) ** 0.5)
print(f"Energy is : {energy}")
l = 152.10e9 * 29290 * m
eccentricity = (1 + (2 * energy * l ** 2)/(m * k ** 2)) ** 0.5
print(f"Eccentricity is : {eccentricity}")

t, q, p = hamiltonian_solve(d_qH, d_pH, d = 2, t_0 = 0.0, q_0 = q_0, p_0 = p_0, h = 1000, N = 1000000, method = "Euler")

r_earth = (- m_sun / (m_earth + m_sun)) * q

# Unpack (x, y) coordinates to be plotted using matplotlib
x, y = r_earth.T

plt.plot(x, y)
plt.show()