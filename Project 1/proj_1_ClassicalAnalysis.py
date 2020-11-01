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
	# g = 6.67408e-11
	k = 1

	return lambda x, p: np.array([k * x[0] / ((x[0] ** 2 + x[1] ** 2) ** 1.5), k * x[1] / ((x[0] ** 2 + x[1] ** 2) ** 1.5)]), lambda x, p: p / m

# Define using earth masses
m_sun = 2
m_earth = 1
d_qH, d_pH = kepler_d_hamiltonian(m_sun, m_earth)

# https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
# Consider coordinates at Aphelion
q_0 = np.array([1, 0])
p_0 = np.array([0, 1])

t, q, p = hamiltonian_solve(d_qH, d_pH, d = 2, t_0 = 0.0, q_0 = q_0, p_0 = p_0, h = 0.01, N = 10000, method = "Euler")

r_earth = (- m_sun / (m_earth + m_sun)) * q

# Unpack (x, y) coordinates to be plotted using matplotlib
x, y = r_earth.T

plt.plot(x, y)
plt.show()

print(p_0)
print(q_0)

print(d_pH(q_0, p_0))
print(d_qH(q_0, p_0))