import numpy as np
from numpy.linalg import norm
from math import sqrt
from copy import copy

def jacobi(a, tol = 1.0e-9):
	"""
	Jacobi Method to Diagonalize Matrices.

	This method takes in a real symmetric matrix, a, and diagonalizes it. 
	I am implementing the method as described in "Numerical Recipes in C: The Art of Scientific Computing, Second Edition", ISBN-10: 0521431085, Chapter 11.1. 
	This section can be found here:
		http://phys.uri.edu/nigh/NumRec/bookfpdf/f11-1.pdf
	I will also adapt techniques used in this implementation to improve the code for computation. 
		https://github.com/mateuv/MetodosNumericos/blob/master/python/NumericalMethodsInEngineeringWithPython/jacobi.py
		https://www.patnauniversity.ac.in/e-content/science/physics/MScPhy58.pdf
	I will try to use the same variables and equations used in the reference, with some exceptions. 
	Equation 11.1.10 sets t to the invese of the quadratic solution. Somehow this black magic works, but I would prefer to use the standard solution from the quadratic as described in equation 28 here:
		http://fourier.eng.hmc.edu/e176/lectures/ch1/node2.html
	Also, while the "cyclic Jacobi method" may be slightly faster computationally, I will implement the original method which searches for the largest off diagonal element. Doing so will make checking against tol easier, and I find it easier to code, at the expense of computational time. Besides, if we were going for computaitonal efficiency, there are better algorithms. 

	Args:
		a (numpy square array): Square symmetric numpy array of real numbers.
		tol (float, optional): Tolerance. Defaults to 1.0e-9.

	Returns: 
		w (numpy (n,  ) array): Vector containing eigenvalues
		v (numpy (n, n) array): Rotation matrix and matrix of eigenvectors
	"""
	def max_off_diag(a):
		# Find largest off-diag. element a[k,l]
		# Search upper right triangle of matrix a for the largest off diagonal element. 
		# Return its value, and position. 
		n = len(a)
		a_max = 0.0
		for i in range(n-1):
			for j in range(i+1, n):
				if abs(a[i,j]) >= a_max:
					a_max = abs(a[i,j])
					p = i
					q = j
		return a_max, p, q

	def rotate(a, v, p, q):
		# Rotate to make a[p, q] = 0
		# Let n be the dimension of matrix a. 
		n = len(a)

		# (11.1.8)
		theta = (a[q,q] - a[p,p]) / (2 * a[p,q])

		# (11.1.10) Modified
		# See Eq (28) here: http://fourier.eng.hmc.edu/e176/lectures/ch1/node2.html
		if theta > 0:
			t = -theta + sqrt(theta**2 + 1.0)
		else:
			t = -theta - sqrt(theta**2 + 1.0)

		# Assign cosine (11.1.11) and sine (11.1.12) values
		c = 1.0/sqrt(t**2 + 1.0)
		s = t*c

		# (11.1.18)
		tau = s / (1 + c)

		# Assign new values to the matrix a (11.1.2). Equivalent to a' = V^T * a * V.
		# For visual reference, see "rotated_matrix.jpg"

		# (11.1.14) Purple
		a[p,p] = a[p,p] - t * a[p,q]
		# (11.1.15) Pink
		a[q,q] = a[q,q] + t * a[p,q]
		# Apply rotation to make a[p, q] = 0
		# (11.1.13) Yellow
		a[p,q] = 0

		# To reduce calculations, only update elements in upper right triangle of matrix a. 
		# Also note that a[i, p] == a[p, i]. Index adjustment is made to only update upper right triangle. 
		# (11.1.16)
		# (11.1.17)
		for i in range(p):
			# Orange
			a_ip = a[i, p]
			a_iq = a[i, q]
			a[i, p] = a_ip - s * (a_iq + tau * a_ip)
			a[i, q] = a_iq + s * (a_ip - tau * a_iq)
		for i in range(p+1, q):
			# Green
			a_pi = a[p, i]
			a_iq = a[i, q]
			a[p, i] = a_pi - s * (a_iq + tau * a_pi)
			a[i, q] = a_iq + s * (a_pi - tau * a_iq)
		for i in range(q+1, n):
			# Blue
			a_pi = a[p, i]
			a_qi = a[q, i]
			a[p, i] = a_pi - s * (a_qi + tau * a_pi)
			a[q, i] = a_qi + s * (a_pi - tau * a_qi)

		# Update the transformation matrix V
		# The same equations as (11.1.16) and (11.1.17) are used, and only columns p and q are updated.
		# (11.1.24)
		for i in range(n):
			v_ip = v[i, p]
			v_iq = v[i, q]
			v[i,p] = v_ip - s * (v_iq + tau * v_ip)
			v[i,q] = v_iq + s * (v_ip - tau * v_iq)

	# Since the rotate method mutates matrix a, make a copy to work with. 
	a_copy = copy(a)
		
	# Let n be the dimension of matrix a. 
	n = len(a_copy)

	# Define max number of rotations. If computation needs more than this, break. 
	max_num_rotations = 5 * (n**2)
	num_rotations = 0

	# Initialize transformation matrix V (11.1.23)
	v = np.identity(n)

	# Iterate Jacobi rotation until matrix a is diagonalized 
	for i in range(max_num_rotations): # Jacobi rotation loop 
		# Find the largest off diagonal element, and the indices of the element. 
		a_max, p, q = max_off_diag(a_copy)

		# If the largest off diag element found is smaller than the tolerance given, then we can say that we are close enough to diagonal and return our current matrix. 
		# if a_max < tol: 
		if check_tolerance(a_copy, tol):
			# diagonal(a) will convert the diagonals into a vector. These will be the eigenvalues. 
			# V will be the rotation matrix, and its columns will be eigenvectors. 
			# print(f"Completed in {num_rotations} rotations")
			return np.diagonal(a_copy), v
		
		# Rotate the matrix, providing the matrix, p, and the indices of the largest element. 
		rotate(a_copy, v, p, q)
		num_rotations += 1
	
	# If execution flow reaches this point, matrix a did not diagonalize to the specified tolerance within the maximum number of rotations. 
	raise Exception(f"Jacobi method did not converge in {max_num_rotations} rotations.")

def check_tolerance(a, tol = 1.0e-9):
	"""
	We want to specify a more complex method of checking tolerance.

	Instead of checking if the largest off diagonal element is smaller than the tolerance, we compare the off diagonal elements to the norm. 

	Args:
		a (Numpy Matrix): Numpy Matrix
		tol (float): Tolerance

	Returns:
		bool: Return True if the tolerance is met. 
	"""
	# Numpy has a built in method to calculate the Frobenius Norm, however, since we are only updating the upper right triangle of matrix a, it would not return the right result. 

	# Find sum of off diagonals squared. 
	off = 0
	n = len(a)
	for i in range(n-1):
		for j in range(i+1, n):
			# Iterate through all off diagonals in upper right triangle 
			off += a[i, j] ** 2
	off = 2 * off

	# Find sum of squares of diagonal elements
	diag = np.sum(np.diagonal(a) ** 2)

	# fnorm is sum of all elements squared
	fnorm = off + diag

	# Check if sqrt(off) <= tol * sqrt(fnorm)
	return sqrt(off) <= tol * sqrt(fnorm)


def hermitian_matrix_split(h):
	# I don't know if this method has an official name. 
	# Split hermitian matrix h into a real symmetric matrix s. 
	a = np.real(h)
	b = np.imag(h)
	s_top = np.concatenate((a, -b), axis = 1)
	s_bot = np.concatenate((b, a), axis = 1)
	s = np.concatenate((s_top, s_bot), axis = 0)
	return s

def sort_eigen_pair(w, v):
	# Get order to sort the eigenvalues and eigenvectors in
	order = w.argsort()
	# Sort eigenvalues
	w = w[order]
	# Sort eigenvectors by column
	v = v.T[order].T
	return w, v

def hermitian_eigensystem(H, tolerance = 1.0e-9):
	"""
	Solves for the eigenvalues and eigenvectors of a hermitian matrix
	
	Args:
		H: Hermitian matrix for which we want to compute eigenvalues and eigenvectors
		
		tolerance: A number that sets the tolerance for the accuracy of the computation.  This number
		is multiplied by the norm of the matrix H to obtain a number delta.  The algorithm successively
		applies (via similarity transformation) Jacobi rotations to the matrix H until the sum of the
		squares of the off-diagonal elements are less than delta.
	
	Returns:
		d: Numpy array containing eigenvalues of H in non-decreasing order
		
		U: A 2d numpy array whose columns are the eigenvectors corresponding to the computed
		eigenvalues.
		
	Checks you might need to do:
		
		H * U[:,k] = d[k] *　U[:,k]      k=0,1,2,...,(n-1)
		
		d[0] <= d[1] <= ... <= d[n-1]     (where n is the dimension of H)
	   
		np.transpose(U) * U = U * np.transpose(U) = np.eye(n)
		
	"""
	n = len(H)

	# Split the Hermitian matrix into a real symmetric matrix
	a = hermitian_matrix_split(H)

	# Apply Jacobi on the real matrix
	w, v = jacobi(a, tolerance)

	# So it turns out that applying the Jacobi method to diagonalize a matrix will shuffle the order of the eigenvalues. This means that even though matrix a started with separate sections, the result will have a shuffled order. If we consider all of the eigenvalues in w, we see that each eigenvalue has multiplicity of two. So we need to pick one of each. The way that I will do this is to sort them, and then choose every other eigenvalue. 
	# I also need to find the eigenvectors that correspond to the eigenvalues. I can split the rotation matrix into its four quadrants, recombine them into complex matrices, and then pair them with the eigenvalues. After doing so, I both lists as normal, and choose every other column to be the corresponding eigenvector. 

	# Split rotation matrix into separate parts. 
	v_real_left = v[:n, :n]
	v_imag_left = v[n:2*n, :n]
	v_left = v_real_left + v_imag_left * 1.0j

	v_real_right = v[n:2*n, n:2*n]
	v_imag_right = -v[:n, n:2*n]
	v_right = v_real_right + v_imag_right * 1.0j

	# Concatenate the two eigenvector matrices to have as many columns as there are eigenvalues in w. 
	v = np.concatenate((v_left, v_right), axis = 1)
	w, v = sort_eigen_pair(w, v)

	# Choose every other column in the sorted matrices. 
	w = w[::2]
	v = v.T[::2].T

	return w, v