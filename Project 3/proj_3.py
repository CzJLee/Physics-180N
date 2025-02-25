import numpy as np
import random
import math
import matplotlib.pyplot as plt

def weighted_die(num_steps):
	"""
	Perform a MCMC simulation on a weighted dice, and return the expected earnings for a predetermined weighted dice and payout. 

	This method runs an Markov Chain Monte Carlo (MCMC) simulation to simulate the results of rolling a weighted dice with six sides. The dice is weighted such that sides 1 and 2 are each three times more likely to land face up than sides 3, 4, 5, and 6. The payout is such that every time the die lands on side 1 or side 2, you win one dollar, and every time the die lands on side 3, 4, 5, or 6, you lose one dollar. 
	
	The specific algorithm used is the Metropolis Algorithm, and the argument num_steps is the number of times this algorithm is performed. In other words, num_steps is the number of dice rolls simulated. 

	Args:
		num_steps (int): Length of MCMC sequence. Number of dice rolls. 

	Returns:
		float: This method returns the expected winnings after num_steps, which is an integer value.
	"""
	winnings = 0
	results = []

	# Define simple payout function that returns 1 if the dice roll is 1 or 2. 
	payout = lambda roll: 1 if roll in {1, 2} else -1
	
	# For the purpose of the Metropolis algorthim, we need to know the ratio of each dice roll compared to another. Since the number of states (6) is small, lets store this data in an array. Let p[m][n] be ratio of the probability of m to the probability of n. 
	p = np.array([
		[1, 1, 3, 3, 3, 3],
		[1, 1, 3, 3, 3, 3],
		[1/3, 1/3, 1, 1, 1, 1],
		[1/3, 1/3, 1, 1, 1, 1],
		[1/3, 1/3, 1, 1, 1, 1],
		[1/3, 1/3, 1, 1, 1, 1]
	])

	# Let s be the list of all possible states. 
	s = {1, 2, 3, 4, 5, 6}

	# Expect num_steps >= 1
	if num_steps < 1:
		raise ValueError("num_steps must be >= 1")

	# Set the first state.
	current_state = random.choice(tuple(s))

	# Run Metropolis Algorithm for num_steps
	while num_steps > 0:
		# Propose a new state that is different than current state. 
		# As a side note, I thought for a while how to represent the different dice states. I wanted to use a set for the easy use of symmetric_difference. But doing so requires converting the set to a tuple each time we choose a random element. I decided to go with this original idea, rather than looking for something more efficient, as this makes sense. 
		# Uniformly pick a random state from the set of states, excluding the current state. 
		proposed_state = random.choice(tuple(s.symmetric_difference({current_state})))

		# Now calculate acceptance probability. This is the chance that the proposed_state is accepted. Else, choose the current_state. 
		# Remember to modify index to match dice rolls with array indexing from 0. 
		acceptance = min(1, p[proposed_state - 1][current_state - 1])

		# Random float in range [0, 1). 
		if random.random() < acceptance:
			# We accept the proposed_state. 
			current_state = proposed_state
		else:
			# We reject the current state. Current state stays the same. 
			pass

		# Determine winnings
		winnings += payout(current_state)

		# Add dice roll and payout to results list.
		results.append([current_state, payout(current_state)])

		num_steps -= 1

	return winnings, results

class Ising_2d:
	def __init__(self, L, temp, H = 0):
		# Critical temperature = 2.2692
		self.temp_crit = 2 / (math.log(1 + math.sqrt(2)))

		# Side length of lattice.
		self.L = L
		self.N = L ** 2

		# Temperature
		self.T = temp

		# External Magnetic Field
		self.H = H

		# Init lattice matrix. Lets call this a. Set all spins to be aligned. 
		self.a = np.ones((L, L))

		# Init Energy (E) and Net Magnetization (S)
		self.E = self.energy()
		self.S = self.spin()

		# Init expectaton values 
		self.E_exp = 0
		self.S_exp = 0
		self.E_squared_exp = 0
		self.S_squared_exp = 0
	
	def __str__(self):
		return str(self.a)

	def set_rand_state(self):
		# Fill the lattice matrix with +1 or -1 spin states randomly with a uniform bias. 
		# Fill an LxL array with 0s and 1s. 
		neg_spins = np.random.randint(2, size = (self.L, self.L))
		# Set a new array with either 1 or -1 by subtracting 2 based on neg_spins. 
		self.a = np.ones((self.L, self.L)) - 2 * neg_spins

		# Init Energy (E) and Net Magnetization (S)
		self.E = self.energy()
		self.S = self.spin()

	def set_temp(self, temp):
		self.T = temp

	def mcmcm(self, num_steps, calculate_vals = False, verbose = False, converge_stop = False, converge_value = "U", converge_range = 0.25, converge_threshold = 0.01):
		"""
		Metropolis Algorithm for Ising Model on a square lattice.

		Args:
			num_steps (int): Number of spin configurations. num_steps >= 1.
			calculate_vals (bool): If True, calculate and return a dict of vals with keys "E", "S", "U", "M", "MS", and "C".
			verbose (bool): If True, print the current step progress
			converge_stop (bool): If True, end simulation after converge_value is within converge_threshold. 
			converge_value (string): The value to check convergence on. Options are "E", "S", "U", "M", "MS", and "C".
			converge_range (float): Value between 0 and 1. Percent of most recent calculated values to check convergence against. 
			converge_threshold (float): Min and Max values within converge_range must be within threshold to declare convergence. 

		Returns: 
			dict: Returns a dict with keys as the strings in calculate, and values as a list of values for each step. The list has length num_steps, where the n'th index is the state AFTER the n'th step. 
		"""

		if calculate_vals:
			vals = {}
			# Init vals dict with empty lists.
			vals["E"] = []
			vals["S"] = []
			vals["U"] = []
			vals["M"] = []
			vals["MS"] = []
			vals["C"] = []

		n = 1
		while n <= num_steps:
			if verbose and n%(num_steps/1000) == 0:
				print(f"Calculating step {n}/{num_steps}", end = "\r")

			# Pick a random site i on the 2D lattice and compute the energy change 𝚫E due to the change of sign in s_i
			rand_site = tuple(np.random.randint(self.L, size = 2))
			del_energy = self.H
			for nn_index in self.nn(rand_site):
				del_energy += self.a[nn_index]
			del_energy *= 2 * self.a[rand_site]

			# If 𝚫E <= 0 then accept the move. Else, accept the move with probability A = exp(-𝚫E/T)
			if del_energy <= 0:
				# Spin flip accepted, flip spin of rand_site
				self.a[rand_site] = -1 * self.a[rand_site]

				# Only update E and S if the state changes. 
				if calculate_vals:
					self.update_E(rand_site)
					self.update_S(rand_site)
			else:
				# Accept the move with probability A = exp(-𝚫E/T)
				P_A = math.exp(-del_energy / self.T)
				# Generate a random float [0, 1)
				if random.random() < P_A:
					# Spin flip accepted, flip spin of rand_site
					self.a[rand_site] = -1 * self.a[rand_site]

					# Only update E and S if the state changes. 
					if calculate_vals:
						self.update_E(rand_site)
						self.update_S(rand_site)
				else:
					# Spin flip is rejected, nothing changes. 
					pass

			# Update values. These update regardless of whether spin state changed. 
			if calculate_vals:
					self.update_E_exp(n)
					self.update_S_exp(n)
					vals["E"].append(self.E)
					vals["S"].append(self.S)
					vals["U"].append(self.U())
					vals["M"].append(self.M())
					vals["MS"].append(self.MS())
					vals["C"].append(self.C())

			# Check convergence condition. Only check every 1% of num_steps to speed up simulation.
			if converge_stop and n%(num_steps/100) == 0:
				# Get range of values
				check_vals = vals[converge_value][int( (1 - converge_range) * n )::]
				# Find min and max
				vals_max = max(check_vals)
				vals_min = min(check_vals)
				diff = abs(vals_max - vals_min)
				if diff < converge_threshold:
					break

			n += 1

		if calculate_vals:
			return vals

	# Use class methods to calculate values

	def nn(self, index):
		"""
		Return list of nearest neighbor indices of index

		Args:
			index (list): Expect list of form [a, b] which are indices of an element in self.a

		Returns:
			list: Return the list of four nearest neighbor indices. Top, Right, Bottom, Left.
		"""
		# Unpack indices to make my life easier. 
		a, b = index

		# Since nearest neighbors wrap around to other side of array, use mod to calculate wrapped values
		top = ((a - 1)%self.L, b)
		bottom = ((a + 1)%self.L, b)
		right = (a, (b + 1)%self.L)
		left = (a, (b - 1)%self.L)

		# Return indices in order top, right, bottom, left
		return top, right, bottom, left

	def energy(self):
		# H is the external magnetic field.
		# For this lab, we assume H = 0.

		# The energy is the negative sum of all nearest neighbor products. 
		# i.e. for every element s_i, 
		# sum: s_i * s_i_right + s_i * s_i_left + s_i * s_i_top + s_i * s_i_bottom
		# Where s_i_right is the right neighbor of s_i. 

		total_energy = 0

		# Iterate over every element in the array.
		for i in range(self.L):
			for j in range(self.L):
				for nn_index in self.nn([i, j]):
					# For each nearest neighbor of index [i, j]
					total_energy += self.a[i, j] * self.a[nn_index]

		# Add the external magnetic field influence 
		if self.H:
			total_energy += self.H * np.sum(self.a)
		
		# Return negative sum
		return -total_energy

	def set_E(self):
		self.E = self.energy()
	
	def update_E(self, index):
		"""
		Update self.E by applying 𝚫E

		Returns:
			float: Energy.
		"""
		del_energy = self.H
		for nn_index in self.nn(index):
			del_energy += self.a[nn_index]
		del_energy *= 2 * self.a[index]

		# Update E
		self.E += del_energy

	def update_E_exp(self, n):
		"""
		Update energy expectation value.

		Args:
			n (int): Current step number.

		Returns:
			float, float: Energy expectation value. Energy squared expectaton value
		"""
		# Based on <O>_n+1 equation in lab proj_3_instructions
		self.E_exp += (1 / n) * (self.E - self.E_exp)
		self.E_squared_exp += (1 / n) * ((self.E) ** 2 - self.E_squared_exp)

		return self.E_exp, self.E_squared_exp

	def spin(self):
		"""
		Net Magnetization (S)

		Sum of all spin states in array.
		"""
		return np.sum(self.a)
	
	def set_S(self):
		self.S = self.spin()

	def update_S(self, index):
		"""
		Update self.S by applying 𝚫S

		Returns:
			float: Energy.
		"""
		del_S = 2 * self.a[index]

		# Update S
		self.S += del_S

	def update_S_exp(self, n):
		"""
		Update Magnetization expectation value. 

		Args:
			n (int): Current step number.

		Returns:
			float, float: Magnetization expectation value. Magnetization squared expectaton value
		"""
		# Based on <O>_n+1 equation in lab proj_3_instructions
		self.S_exp += (1 / n) * (self.S - self.S_exp)
		self.S_squared_exp += (1 / n) * ((self.S) ** 2 - self.S_squared_exp)

		return self.S_exp, self.S_squared_exp

	def U(self):
		"""
		Mean Internal Energy
		"""
		return self.E_exp / self.N

	def M(self):
		"""
		Magnetization
		"""
		return abs(self.S_exp) / self.N

	def MS(self):
		"""
		Magnetic Susceptibility
		"""
		return (self.S_squared_exp - (self.S_exp) ** 2) / (self.N * self.T)

	def C(self):
		"""
		Specific Heat
		"""
		return (self.E_squared_exp - (self.E_exp) ** 2) / (self.N * (self.T ** 2))

	def show_state(self):
		"""
		Show heatmap plot of current spin states
		"""
		plt.imshow(self.a, cmap='Greys', interpolation='nearest')
		plt.show()
