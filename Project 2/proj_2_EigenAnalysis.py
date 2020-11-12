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

dim = 20

# Create a list of diagonals for H
h_diag = [m + 0.5 for m in range(1, dim + 1)]
# Create a matrix of zeros, and fill the diagonals with the diagonal elements. 
h = np.zeros((dim, dim))
np.fill_diagonal(h, h_diag)

x = np.zeros((dim, dim))

# Iterate over every element in the matrix, and check for deltas = 1.
for n in range(dim):
	for m in range(dim):
		if n == m + 1:
			x[n, m] = sqrt(m + 1)
		elif n == m - 1:
			x[n, m] = sqrt(m)
			
x_4 = np.linalg.multi_dot([x, x, x, x])

h_l = h + x_4

w, v = np.linalg.eigh(h_l)
p = []

x = [i / 100 for i in range(-500, 501)]
for i, _ in enumerate(x):
	p.append(psi(v.T[1], x[i]))

plt.plot(x, p)
plt.show()