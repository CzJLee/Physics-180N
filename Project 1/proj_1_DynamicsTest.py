import math
import numpy as np
import matplotlib.pyplot as plt
from proj_1_module import dynamics_solve, hamiltonian_solve

def exact_sol(b, d):
	return lambda p0, t: p0 * math.exp((b - d) * t)

def population_model(b, d):
	return lambda t, s: (b - d) * s

def exact_solve(f, d = 1, t_0 = 0.0, s_0 = 1, h = 0.1, N = 100):
	T = np.array([t_0 + n * h for n in range(N + 1)])
	
	if d == 1:
		S = np.zeros(N + 1)
	
	if d > 1:
		S = np.zeros((N + 1, d))
		
	S[0] = s_0

	for n in range(N + 1):
		S[n] = f(s_0, T[n])
	
	return T, S

p = exact_sol(1.1, 1)
t, s = exact_solve(p)

print(np.shape(t))
print(np.shape(s))

plt.plot(t, s)
plt.show()