from proj_2_EigenTest import gen_rand_herm, npprint
import numpy as np
from math import sqrt
from test import jacobi

a = np.array([[1,sqrt(2),2],[sqrt(2),3,sqrt(2)],[2,sqrt(2),1]])

e, a = gen_rand_herm(3)

aj, p = jacobi(a)

npprint(aj)
npprint(p)