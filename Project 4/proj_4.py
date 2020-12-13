import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
# plt.style.use('ggplot')

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
		self.temp[0, 0] = self.bc_left if self.bc_left is not None else temp
		self.temp[-1, 0] = self.bc_right if self.bc_right is not None else temp

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
			plt.plot(self.length, temp)

		plt.xlabel("Length (cm)")
		plt.ylabel("Temperature (°C)")
		plt.title("Temperature distribution of 1D bar at various times")

		# Convert times into list of times
		time_str = [f"{t} s" for t in time]
		plt.legend(time_str, loc = "upper right")
		plt.show()

	def plot_animate(self):
		# Draw animated plot over the time period
		fig, ax = plt.subplots()
		# Set axis bounds
		ax.set(xlim=(0, self.l), ylim=(0, 100))

		T = self.temp
		x = self.length

		line = ax.plot(x, T[:, 0], lw=2)[0]
		ax.set_xlabel("Length (cm)")
		ax.set_ylabel('Temperature (°C)')

		def animate(i):
			line.set_ydata(T[:, i])
			ax.set_title(f"Temperature distribution of 1D bar at various times\nTime: {i * self.dt :.2f} s")

		anim = FuncAnimation(
			fig, animate, interval=10, frames=self.n_max)

		plt.draw()
		plt.show()

class Diffusion_3D(Diffusion_1D):
	def __init__(self, a = (1.10407 ** 0.5), l = 25, b = 1/6, h = 1, t_max = 180):
		super().__init__(a = a, l = l, b = b, h = h, t_max = t_max)

	def set_boundary_conditions(self, temp = 100, bc_right = 0):
		self.bc_left = 0
		self.bc_right = bc_right

		# Create an array with length = (l // h) + 1 
		# (+1 to account for both ends)
		# Element of index m in the array will represent the temperature at position x = m * h
		self.temp = np.ones((self.m_max, 1)) * temp
		# Multiple by r to transform T -> V
		for i, m in enumerate(self.temp):
			m[0] = m[0] * self.length[i]
		# Apply BC
		self.temp[-1, 0] = self.bc_right if self.bc_right is not None else temp

	def calc_parabola(self, x1, y1, x2, y2, x3, y3):
		# Now we need to find the temp at the center where r = 0. The paper says that they use a three point interpolation scheme. I am assuming that they fit the three nearest points as a parabola to find the value for r = 0. I will do the same. 
		# First we need to calculate the equation for the parabola, then determine the value from that equation. 
		# From http://chris35wills.github.io/parabola_python/
		'''
		Adapted and modifed to get the unknowns for defining a parabola:
		http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
		'''

		denom = (x1-x2) * (x1-x3) * (x2-x3)
		A     = (x3 * (y2-y1) + x2 * (y1-y3) + x1 * (y3-y2)) / denom
		B     = (x3*x3 * (y1-y2) + x2*x2 * (y3-y1) + x1*x1 * (y2-y3)) / denom
		# We are only really interested in C, since we will be plugging in x = 0 for ax^2 + bx + c
		C     = (x2 * x3 * (x2-x3) * y1+x3 * x1 * (x3-x1) * y2+x1 * x2 * (x1-x2) * y3) / denom

		return A,B,C

	def radial_transform(self, temps):
		# As described below equation (7), to convert from 1D back to 3D sphere
		# All points except r = 0 can be divided by r. 
		# For given array of temps, divide by self.length
		# Use numpy.divide, but slice out first element to avoid dividing by zero. 
		temps_radial = [temps[i] / self.length[i] for i in range(1, self.m_max)]

		# Get values for the coordinate points for the three closest values
		y1, y2, y3 = temps_radial[:3]
		x1, x2, x3 = self.length[1:4]
		a, b, c = self.calc_parabola(x1, y1, x2, y2, x3, y3)

		# The temp we want for r = 0 is the value for C
		# Insert this into the temp_radial parabola
		return np.insert(temps_radial, 0, c)

	def convert_temp_to_radial(self):
		# Apply the transformation as described below equation (7) to convert all temps from 1D to Radial
		conv_temp = []
		# Transpose matrix so that rows are length segments for each time
		for segment in self.temp.T:
			conv_temp.append(self.radial_transform(segment))
		
		# Apply change to temp variable
		self.temp = np.array(conv_temp).T