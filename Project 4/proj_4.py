import numpy as np
import matplotlib.pyplot as plt

class Diffusion_1D():
	def __init__(self, a = (1.10407 ** 0.5), l = 20, b = 1/6, h = 1):
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

	def set_boundary_conditions(self, temp = 100, bc_left = 0, bc_right = 0):
		# Create an array with length = (l // h) + 1 
		# (+1 to account for both ends)
		# Element of index m in the array will represent the temperature at position x = m * h
		self.temp = np.ones(((self.l // self.h) + 1, 1)) * temp
		self.temp[0, 0] = bc_left
		self.temp[-1, 0] = bc_right
		self.bc_left = bc_left
		self.bc_right = bc_right

	def simulate_diffusion(self, t_max = 120):
		# Create array of time elements
		self.time = np.arange(0, t_max, self.dt)

		# For each time step, apply equation (5) and store results as a new entry in self.temp
		for n in range(int(t_max // self.dt)):
			# Append a column to right for next time step
			self.temp = np.append(self.temp, np.empty(((self.l // self.h) + 1, 1)), axis = 1)


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

			for m in range(1, int(self.l // self.h)):
				self.temp[m, n + 1] = self.temp[m, n] + self.b * (self.temp[m + 1, n] + self.temp[m - 1, n] - 2 * self.temp[m, n])

	def get_temp_at_time(self, n):
		return self.temp[:, n]

model = Diffusion_1D()
model.set_boundary_conditions()
model.simulate_diffusion()
# print(model.dt)
print(model.temp)
# print(model.time)

n = int(75 // model.dt)
blah = model.get_temp_at_time(n)
print(blah)

lengths = np.arange(0, 21, 1)
print(lengths)


plt.plot(lengths, blah)
plt.show()