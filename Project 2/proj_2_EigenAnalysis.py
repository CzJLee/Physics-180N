from proj_2_EigenTest import gen_rand_herm, npprint, is_hermitian, is_eigenvector, verify_eigenvectors
from proj_2_module import jacobi, hermitian_matrix_split, sort_eigen_pair, hermitian_eigensystem
import numpy as np
from numpy.polynomial.hermite import hermval
from math import sqrt, factorial, pi, exp
import matplotlib.pyplot as plt

def psi(v, x):
	"""
	Calculate the real value for psi at a given x.

	v is an eigenvector with elements [c_0, c_1, ... , c_n-1] and length n
	This function will calculate and return c_0 * phi_0(x) + c_1 * phi_1(x) + ... + c_n-1 * phi_n-1(x)
	Where phi is (2**n * factorial(n) * sqrt(pi))**(-0.5) * exp(-x**2 / 2) * h_n(x)
	Where h_n(x) is the n'th hermitian polynomail evaluated at x. 

	Use numpy.polynomial.hermite.hermval to compute h_n(x) values. 
	https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermval.html

	We must create a list of coefficients for hermval. The coefficients will be 
	c_n * (2**n * factorial(n) * sqrt(pi))**(-0.5) * exp(-x**2 / 2)

	Args:
		v (vector): Eigenvector
		x (float): Position to calculate psi at

	Returns:
		float: Wave function value at position x
	"""
	c = np.zeros(len(v))

	for n in range(len(v)):
		c[n] = np.real(v[n]) * (2**n * factorial(n) * sqrt(pi))**(-0.5) * exp(-x**2 / 2)

	return hermval(x, c)

def plot_eigenfunction(l = 0, energy_level = 0, dim = 30):
	"""
	Plot the wave function for the anharmonic oscillator. 

	Use the eigenvectors to plot the wave function of the anharmonic oscillator with x^4 perturbation coefficient l at a given energy level. Diagonalize the Hamiltonain matrix to obtain the eigenvectors. 

	Args:
		l (float, optional): Lambda in H(lambda) = H + lambda * x^4. Defaults to 0.
		energy_level (int, optional): Energy Level to plot. Ranges from [0, Infinity). Defaults to 0.
		dim (int, optional): Dimension of matrix to use. Must be larger than energy_level. Defaults to 30.
	"""
	### Create H Matrix ###
	# Create a list of diagonals for H
	h_diag = [m + 0.5 for m in range(1, dim + 1)]
	# Create a matrix of zeros, and fill the diagonals with the diagonal elements. 
	h = np.zeros((dim, dim))
	np.fill_diagonal(h, h_diag)

	### Create x^4 Matrix ###
	x = np.zeros((dim, dim))

	# Iterate over every element in the matrix, and check for deltas = 1.
	for n in range(dim):
		for m in range(dim):
			if n == m + 1:
				x[n, m] = sqrt(m + 1)
			elif n == m - 1:
				x[n, m] = sqrt(m)
				
	x_4 = np.linalg.multi_dot([x, x, x, x])

	# Check if l and energy_levels are floats or lists. Use lists for multiple overlaid plots. 
	if not hasattr(l, "__iter__"):
		l = [l]
	if not hasattr(energy_level, "__iter__"):
		energy_level = [energy_level]

	for i in range(len(l)):
		### Create H(lambda) Matrix ###
		# Add to make H(lambda) Matrix
		h_l = h + l[i] * x_4

		### Diagonalize Hamiltonain Matrix ###
		w, v = np.linalg.eigh(h_l)

		### Plot Wavefunction ###
		p = []

		# Plot from -5 to 5. 
		x = np.arange(-5, 5, 0.01)
		for j, _ in enumerate(x):
			# Use psi function to calculate value at a given point x. 
			p.append(psi(v.T[energy_level[i]], x[j]))

		plt.plot(x, p)
	
	plt.show()