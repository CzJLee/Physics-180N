import math
import numpy as np
import matplotlib.pyplot as plt
from proj_1_module import dynamics_solve, hamiltonian_solve

def d_hamiltonian(m, w):
	return lambda x, p: m * w**2 * x, lambda x, p: p / m

d_qH, d_pH = d_hamiltonian(1, 1)
h = 0.01
N = int(100/h)

t, q, p = hamiltonian_solve(d_qH, d_pH, h=h, N=N, method="SE")

plt.plot(t, q)
plt.show()