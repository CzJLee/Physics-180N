import numpy as np
from math import sqrt
# https://numpy.org/doc/stable/reference/generated/numpy.linalg.eigh.html
from numpy.linalg import eigh
from scipy.stats import unitary_group

def npprint(a, text = None):
	if text:
		print(text)
	with np.printoptions(precision=4, suppress=True, formatter={'float': '{:0.3f}'.format}, linewidth=100):
		print(a)

def gen_rand_herm(dim = 3, eigenvalues = None):
	if eigenvalues is None:
		# Generate an array of random integers to use as random eigenvalues
		# I will sample integers from the range [-10, 11). 
		eigenvalues = np.random.randint(low = -10, high = 11, size = dim)

	# Now make a diagonal matrix with the eigenvalues as its diagonal elements. 
	a = np.zeros((dim, dim))
	np.fill_diagonal(a, eigenvalues)

	# Apply a similarity transformation by a unitary matrix to get a non diagonal matrix with the same eigenvalues
	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.unitary_group.html
	x = unitary_group.rvs(dim)
	a = np.linalg.multi_dot((x, a, x.conj().T))

	return eigenvalues, a

def is_hermitian(h):
	# Check if h is hermitian.
	return np.allclose(h, h.conj().T)

def is_eigenvector(a, l, u):
	# Return True if l is an eigenvalue and u is an eigenvector within tolerance. 
	# Verify that A * u = l * u
	return np.allclose(np.dot(a, u), l * u, rtol=1e-4, atol=1e-7)

def verify_eigenvectors(a, w, v):
	# Return True if all eigenvectors in v have eigenvalues in w.
	# Return False if a pair does not match.
	for i in range(len(a)):
		if not is_eigenvector(a, w[i], v.T[i]):
			return False
	return True