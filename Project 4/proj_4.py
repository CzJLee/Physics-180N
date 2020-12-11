import numpy as np
import matplotlib.pyplot as plt

class Diffusion_1D():
	def __init__(self, a = (1.10407 ** 0.5), l = 20, b = 1/6, h = 1, t_max = 120):
		# Diffusivity ⍺ (cm * s^-1/2)
		self.a = a 

		# Length L (cm)
		self.l = l

		# Length Interval h (cm)
		self.h = h

		# Beta parameter (unitless)
		# β = 𝚫t ⍺^2 / h^2
		self.b = b

		# Delta time interval (s)
		# Found by using given β, ⍺, and h.
		# 𝚫t = β * h^2 / ⍺^2
		self.dt = self.b * self.h**2 / self.a**2

		# Total time to run simulation for (s)
		self.t_max = t_max

		# Max integer number of length steps (+1 to account for both ends)
		self.m_max = int((self.l // self.h) + 1)

		# Max integer number of time steps
		self.n_max = int(self.t_max // self.dt)

		# Create array of time elements
		self.time = np.arange(0, self.t_max, self.dt)

		# Create array of length elements
		self.length = np.arange(self.m_max) * self.h

	def set_boundary_conditions(self, temp = 100, bc_left = 0, bc_right = 0):
		self.bc_left = bc_left
		self.bc_right = bc_right

		# Create an array with length = (l // h) + 1 
		# (+1 to account for both ends)
		# Element of index m in the array will represent the temperature at position x = m * h
		self.temp = np.ones((self.m_max, 1)) * temp
		self.temp[0, 0] = self.bc_left
		self.temp[-1, 0] = self.bc_right
		

	def simulate_diffusion(self):
		# Array and code is formatted to match notation of equation (5)
		# For each time step, apply equation (5) and store results as a new entry in self.temp
		for n in range(self.n_max):
			# Iterate over all time steps
			# Append an empty column to right of temp for next time step
			self.temp = np.append(self.temp, np.empty((self.m_max, 1)), axis = 1)

			if self.bc_left is not None:
				# If there is a left boundary condition
				# Set that position to be equal to the bc
				self.temp[0, n + 1] = self.bc_left
			else:
				# BC is set to None, no BC, apply modified equation 5
				self.temp[0, n + 1] = self.temp[0, n] + self.b * (self.temp[1, n] - self.temp[0, n])

			if self.bc_right is not None:
				# If there is a right boundary condition
				# Set that position to be equal to the bc
				self.temp[-1, n + 1] = self.bc_right
			else:
				# BC is set to None, no BC, apply modified equation 5
				self.temp[-1, n + 1] = self.temp[-1, n] + self.b * (self.temp[-2, n] - self.temp[-1, n])

			# Apply equation (5) to each non-edge length segment.
			for m in range(1, self.m_max - 1):
				# Iterate over all inner length segments
				self.temp[m, n + 1] = self.temp[m, n] + self.b * (self.temp[m + 1, n] + self.temp[m - 1, n] - 2 * self.temp[m, n])

	def get_temp_at_time(self, n):
		return self.temp[:, n]

	def plot_temp_at_time(self, time):
		# Get the index of the given time
		times = [int(t // self.dt) for t in time]

		# Get the column of temp data at the given time index
		temp_at_time = [self.get_temp_at_time(n) for n in times]

		# Plot
		for temp in temp_at_time:
			plt.plot(model.length, temp)

		plt.xlabel("Length (cm)")
		plt.ylabel('Temperature (°C)')
		plt.title('Temperature distribution at various times')

		# Convert times into list of times
		time_str = [f"{t} s" for t in time]
		plt.legend(time_str, loc = "upper right")
		plt.show()

model = Diffusion_1D()
model.set_boundary_conditions()
model.simulate_diffusion()
model.plot_temp_at_time([0, 10, 30, 60, 75])