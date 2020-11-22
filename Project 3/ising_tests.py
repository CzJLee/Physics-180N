from proj_3 import Ising_2d, weighted_die
import numpy as np
import matplotlib.pyplot as plt

num_steps = 500000
L = 64
time_step_size = num_steps / (L ** 2)

# Set temperature to be 3. This is larger than the critical temperature. 
temp = 3

# Create our Ising model using L = 16, temp = 3.
model = Ising_2d(L, temp)
print(model.energy())
# model.set_rand_state()

# Let's run the MCMC for multiple steps, and calculate values for U and M for each step. 
vals = model.mcmcm(num_steps, calculate_vals = True, verbose = True)

# Now lets create our list of time steps to plot
time_steps = [time_step_size * i for i in range(1, num_steps + 1)]

plt.plot(time_steps, vals["M"])
plt.show()

print()
print(model.energy())