from proj_2_EigenTest import gen_rand_herm, npprint, is_hermitian, is_eigenvector, verify_eigenvectors
import numpy as np
from math import sqrt
from proj_2_module import jacobi, hermitian_matrix_split, sort_eigen_pair, hermitian_eigensystem
from scipy.stats import ortho_group

def gen_rand_sym(dim = 3, eigenvalues = None):
	if eigenvalues is None:
		# Generate an array of random integers to use as random eigenvalues
		# I will sample integers from the range [-10, 11). 
		eigenvalues = np.random.randint(low = -10, high = 11, size = dim)

	# Now make a diagonal matrix with the eigenvalues as its diagonal elements. 
	a = np.zeros((dim, dim))
	np.fill_diagonal(a, eigenvalues)

	# Apply a similarity transformation by a unitary matrix to get a non diagonal matrix with the same eigenvalues
	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.unitary_group.html
	x = ortho_group.rvs(dim)
	a = np.linalg.multi_dot((x, a, x.conj().T))

	return eigenvalues, a

# for i in range(2, 31):
# 	e, a = gen_rand_sym(i)

# 	# npprint(e)
# 	# npprint(a)

# 	np_w, np_v = np.linalg.eigh(a)
# 	# npprint(np_w)

# 	my_w, my_v = jacobi(a)
# 	my_w, my_v = sort_eigen_pair(my_w, my_v)
# 	# npprint(my_w)

# 	print(f"Eigenvalues: {np.allclose(np_w, my_w)}", end=", ")
# 	print(f"Eigenvectors: {verify_eigenvectors(a, my_w, my_v)}", f" For {i}x{i}")

for n in range(2, 25):
	e, h = gen_rand_herm(n)
	# npprint(np.sort(e))
	split = hermitian_matrix_split(h)

	w, v = jacobi(split)
	savew = w
	firsthalf = w[:n]
	firsthalf = np.sort(firsthalf)
	secondhalf = w[n:2*n]
	secondhalf = np.sort(secondhalf)
	everyother = np.sort(w)[::2]
	jv = v
	vt = v[:n, :n]
	vb = v[n:2*n, n:2*n]

	# Slice w, v
	w = w[:n]
	v_real = v[:n, :n]
	v_imag = v[n:2*n, :n]
	v = v_real + v_imag * 1.0j

	w, v = sort_eigen_pair(w, v)

	# npprint(np.sort(w))
	# npprint(v)

	print(f"Eigenvalues: {np.allclose(np.sort(e), np.sort(w))}", end=", ")
	if not np.allclose(np.sort(e), np.sort(w)):
		print("Source Eigen")
		npprint(np.sort(e))
		print("Solution Eigen")
		npprint(np.sort(w))
		print("np Eigen")
		np_w, np_v = np.linalg.eigh(h)
		npprint(np_w)
		print("Save w")
		npprint(savew)
		npprint(savew[:n])
		print("Are halves equal?")
		print(np.allclose(firsthalf, secondhalf))
		npprint(firsthalf)
		npprint(secondhalf)
		npprint(everyother)
		print("Is every other equal to solution?")
		print(np.allclose(np.sort(e), everyother))
		print("Jacobi V")
		npprint(jv)
		print(is_hermitian(jv))
		print("Are top and bottom V equal?")
		print(np.allclose(vt, vb))
		break
	print(f"Eigenvectors: {verify_eigenvectors(h, w, v)}", f" For {n}x{n}")