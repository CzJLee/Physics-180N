from proj_3 import Ising_2d, weighted_die
import numpy as np
import matplotlib.pyplot as plt
import math

# temps = np.arange(4, 0, -0.05)

# def onsager(t):
# 	if t >= 2 / (math.log(1 + math.sqrt(2))):
# 		return 0
# 	else:
# 		return ((1 - (math.sinh(2 / t) ** (-4))) ** (1/8))

# # Define parameters
# L = 16
# num_steps = L * 10000
# m_t_8 = []
# init_temp = 4

# # Create our Ising model.
# model = Ising_2d(L, init_temp)
# model.set_rand_state()

# for temp in temps:
# 	print()
# 	print(f"Calculating T = {temp:.2f}")
# 	# Set temp
# 	model.set_temp(temp)

# 	# Run mcmc algorithm until convergence. 
# 	vals = model.mcmcm(num_steps, calculate_vals = True, converge_stop = True, converge_value = "M", converge_threshold = 0.01, verbose = True)

# 	# Save current value for M(temp) as the average of 25% most recent values
# 	segment = vals["M"][int((1-0.25)*len(vals["M"]))::]
# 	avg = sum(segment) / len(segment)
# 	m_t_8.append(avg)

# m_t = [onsager(t) for t in temps]

# plt.xlabel("Temperature")
# plt.ylabel('Magnetization')
# plt.title('Magnetization for different temperatures for L = 8')
# plt.plot(temps, m_t, color = "green", linestyle = "dashed")
# plt.scatter(temps, m_t_8, color = "red", marker = ".")
# plt.show()






# # Define parameters

# L = 8
# num_steps = L * 10000
# time_step_size = num_steps / (L ** 2)

# # Set temperature to be 3. This is larger than the critical temperature. 
# temp = 4

# # Create our Ising model.
# model = Ising_2d(L, temp)
# # model.set_rand_state()

# # Let's run the MCMC for multiple steps, and calculate values for U and M for each step. 
# vals = model.mcmcm(num_steps, calculate_vals = True, verbose = True, converge_stop = True, converge_value= "M", converge_threshold=0.02)

# # Now lets create our list of time steps to plot
# time_steps = [time_step_size * i for i in range(1, len(vals["U"]) + 1)]

# if len(vals['U']) < num_steps:
#     print()
#     print(f"System has converged in {len(vals['U']) * time_step_size:.2e} steps.")

# plt.xlabel("Scaled Time")
# plt.ylabel('Mean Internal Energy')
# plt.title('Mean Internal Energy convergence over time')
# plt.plot(time_steps, vals["U"])
# plt.show()

# a = np.random.randint(0, 2, (256, 256))
# plt.imshow(a, cmap='Greys', interpolation='nearest')
# plt.show()



# Define parameters
num_steps = 500000
L = 256
time_step_size = num_steps / (L ** 2)

# Set temperature to be 3. This is larger than the critical temperature. 
temp = 1.8

# Create our Ising model.
model = Ising_2d(L, temp)
# model.set_rand_state()

# Let's run the MCMC for multiple steps, and calculate values for U and M for each step. 
vals = model.mcmcm(num_steps, calculate_vals = True, verbose = True)

# Now lets create our list of time steps to plot
time_steps = [time_step_size * i for i in range(1, num_steps + 1)]

plt.xlabel("Scaled Time")
plt.ylabel('Magnetization')
plt.title('Magnetization convergence over time')
plt.plot(time_steps, vals["M"])
plt.show()

plt.xlabel("Scaled Time")
plt.ylabel('Mean Internal Energy')
plt.title('Mean Internal Energy convergence over time')
plt.plot(time_steps, vals["U"])
plt.show()

model.show_state()