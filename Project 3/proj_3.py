import numpy as np
import random

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