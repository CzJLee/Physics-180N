from proj_2_EigenTest import gen_rand_herm, npprint, is_hermitian, is_eigenvector, verify_eigenvectors
import numpy as np
from math import sqrt
from proj_2_module import jacobi, hermitian_matrix_split, sort_eigen_pair, hermitian_eigensystem

import time


for dim in range(20, 22):
	e, a = gen_rand_herm(dim)
	tstart = time.time()
	w, v = hermitian_eigensystem(a)
	tend = time.time()
	eigenvalues_match = np.allclose(w, np.sort(e))
	eigenvectors_match = verify_eigenvectors(a, w, v)
	if eigenvalues_match and eigenvectors_match:
		dtime = tend - tstart
		print(f"Diagonalized {dim}x{dim} Hermitian matrix in {dtime:.2f} seconds.")
	else:
		print(f"Failed to diagonalize {dim}x{dim} Hermitian matrix.")
		print("Do eigenvalues match?")
		print(np.allclose(np.sort(e), w))
		npprint(np.sort(e))
		npprint(w)