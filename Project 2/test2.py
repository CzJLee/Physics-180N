from proj_2_EigenTest import gen_rand_herm, npprint, is_hermitian
import numpy as np
from math import sqrt
from proj_2_module import jacobi, hermitian_matrix_split, sort_eigen_pair, hermitian_eigensystem



test_array = np.array([[1,sqrt(2),4],[sqrt(2),3,sqrt(2)],[4,sqrt(2),1]])
aj, p = jacobi(test_array)
# npprint(aj)
# npprint(p)

e, a = gen_rand_herm(3)
npprint(e)
npprint(a)

w, v = hermitian_eigensystem(a)
print("My Function w")
npprint(w)
print("My Function v")
npprint(v)
my_v = v

print()
w, v = np.linalg.eigh(a)
print("Numpy w")
npprint(w)
print("Numpy v")
npprint(v)

np_v = v

x = np_v.T[0]
npprint(x)
s = 0
for i in x:
	s += np.abs(i)
print(s)