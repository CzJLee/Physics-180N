from proj_3 import Ising_2d
import numpy as np

# model = Ising_2d(4, 2)
# model.set_rand_state()
# print(model)
# print(model.energy())

l = np.zeros((16, 16))

a = np.random.randint(16, size = 2)
print(a)
print(l[tuple(a)])