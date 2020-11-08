from proj_2_EigenTest import gen_rand_herm, npprint, is_hermitian
import numpy as np
from math import sqrt
from proj_2_module import jacobi, hermitian_matrix_split



test_array = np.array([[1,sqrt(2),4],[sqrt(2),3,sqrt(2)],[4,sqrt(2),1]])
aj, p = jacobi(test_array)
# npprint(aj)
# npprint(p)

e, a = gen_rand_herm(3)
npprint(e)
npprint(a)
npprint(a.conj().T)
print(is_hermitian(a))

s = hermitian_matrix_split(a)
npprint(s)
print(is_hermitian(s))