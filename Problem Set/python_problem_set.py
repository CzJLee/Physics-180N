# Physics 180N Problem Set Functions

def test(function, args, expected_output):
	"""
	Simple unit testing function for functions.

	Args:
		function (function): The function to be tested
		args (tuple): Arguments passed to the given function
		expected_output: Expected output of function given the input arguments
	"""
	if not isinstance(args, tuple):
		# If args only has one element, it is not a tuple. Make it a tuple.
		args = (args,)
	arg_string = str(args) if len(args) > 1 else f"({args[0]})"
	print(f"Expect {function.__name__}{arg_string} to be {expected_output}")
	print(f">>> {function.__name__}{arg_string}")
	actual_output = function(*args)
	print(actual_output, "\n")
	if actual_output == expected_output:
		print("PASS: Output matches expected output.")
	else:
		print("ERROR: Output does NOT match expected output.")

################################################################################

# Problem Group 1 - Recursion

## Problem 1.1 - Factorial
def recursive_factorial(n):
	"""
	Recursive function to calculate n!

	Args:
		n (int): Positive integer

	Raises:
		ValueError: Raise error if input n is negative. 

	Returns:
		int: Returns n! = n * (n-1) * (n-2) * ...
	"""
	if n < 0:
		raise ValueError("n must not be negative")
	if n == 0:
		return 1
	else:
		return n * recursive_factorial(n-1)

## Problem 1.2 - Memoization
def recursive_factorial_memo(n, memo_dict):
	"""
	Recursive function to calculate n factorial, using a memoization dictionary.

	Args:
		n (int): Positive integer
		memo_dict (dict): Memoization dictionary

	Raises:
		ValueError: Raises error if input n is negative.

	Returns:
		int: Returns n! = n * (n-1) * (n-2) * ...
	"""
	if n < 0:
		raise ValueError("n must not be negative")
	if n == 0:
		return 1
	elif n not in memo_dict :
		memo_dict[n] = n * recursive_factorial_memo(n-1, memo_dict)
	return memo_dict[n]

## Problem 1.3 - Binomial Coefficients
def binomial_factorial(n, k):
	"""
	Calculate the Binomial Coefficient of (n, k) using the factorial definition.

	Args:
		n (int): Binomial coefficient n, where n >= k
		k (int): Binomial coefficient k, where n >= k

	Returns:
		int: Binomial Coefficient of (n, k) using the factorial definition
	"""
	n_fac = recursive_factorial(n)
	k_fac = recursive_factorial(k)
	nk_fac = recursive_factorial(n-k)
	return int(n_fac / (k_fac * nk_fac))

def binomial_recursive(n, k):
	"""
	Calculate the Binomial Coefficient of (n, k) using the recurrence relation.

	Args:
		n (int): Binomial coefficient n, where n >= k
		k (int): Binomial coefficient k, where n >= k

	Returns:
		int: Binomial Coefficient of (n, k) using the recurrence relation.
	"""
	if k == 0:
		return 1
	elif k == n:
		return 1
	else:
		return binomial_recursive(n-1, k-1) + binomial_recursive(n-1, k)


## Problem 1.4 - The Logistic Map
def logistic(n, r, x0):
	"""
	Calculate the n'th value in the logistic map function x_(n+1) = r * x_n * (1 - x_n).

	Args:
		n (int): The n'th value of the logistic map function
		r (float): parameter value r
		x0 (float): initial value

	Returns:
		float: Returns the n'th value of the logistic map function for the given parameter r and initial value x0.
	"""
	if n == 0:
		return x0
	else:
		x_step = logistic(n-1, r, x0)
		return r * x_step * (1 - x_step)

def logistic_gen(n, r, x0):
	"""
	Generator function to calculate the n'th value in the logistic map function.

	Args:
		n (int): The n'th value of the logistic map function
		r (float): parameter value r
		x0 (float): initial value

	Yields:
		float: Yields the n'th value of the logistic map function for the given parameter r and initial value x0.
	"""
	yield x0
	x_n = x0
	for i in range(n-1):
		x_n = r * x_n * (1 - x_n)
		yield x_n


################################################################################

# Problem Group 2 - Searching for stuff

## Problem 2.1 - Linear Search
def linear_search(l, n):
	"""
	Simple lineary search for a list of numbers.

	Args:
		l (list): List of positive integers
		n (int): Value to search for

	Returns:
		int: Returns the index of the first match found, or None if no matches found.
	"""
	for i in range(len(l)):
		if n == l[i]:
			return i
	print("The specific number is not in the list.")
	return None

## Problem 2.2 - Bisection Search
def bisection_search(l, n):
	"""
	Wrapper for binary_search algorithm.

	Using a simple binary_search algorithm I have written before.

	Args:
		l (list): Sorted list in ascending order. Elements must be comparable
		n (float): Value to search for

	Returns:
		int, None: Returns the index of the first match found, or None if no matches found. 
	"""
	index = binary_search(l, n)
	if index is not None:
		return index
	else: 
		print("The specific number is not in the list.")
		return None

def binary_search(sorted_list, target):
	"""
	Simple Binary Search using pointers. Iterative approach

	Time Complexity: O(n). Space Complexity: O(1)

	Args:
		sorted_list (list): Sorted list in ascending order. Elements must be comparable. 
		target (float): Value to search for

	Returns:
		int, None: Returns the index of the first match found, or None if no matches found. 
	"""
	# Set the left pointer to the beginning of the array, and the right pointer to the end of the array.
	left_pointer = 0
	right_pointer = len(sorted_list)
	
	while left_pointer < right_pointer:
		# Calculate the middle index using the two pointers
		middle_index = (right_pointer + left_pointer) // 2
		middle_value = sorted_list[middle_index]
		if middle_value == target:
			return middle_index
		if target < middle_value:
			right_pointer = middle_index
		if target > middle_value:
			left_pointer = middle_index + 1
	
	return None

## Problem 2.3 - Bisection Root Finding
def bisection_root(f, x_left, x_right, epsilon):
	"""
	Bisection method to find a root of a continuous function f.

	The algorithm applies to any continuous function f on an interval x_left to x_right where the value of the function f changes sign from x_left to x_right.

	Args:
		f (function): Continuous real-valued function of a single real variable
		x_left (float): Estimated value to the left of the root
		x_right (float): Estimated value to the right of the root
		epsilon (float): Tolerance error to find the actual root within

	Returns:
		float: Root of the function f between the given intervals within tolerance epsilon
	"""
	# Make sure x_left and x_right are entered in the right order
	if x_left > x_right:
		# Swap them
		x_left, x_right = x_right, x_left

	# Check if either of the intervals are a zero themselves.
	if f(x_left) == 0:
		return x_left
	elif f(x_right) == 0:
		return x_right
	elif f(x_left) * f(x_right) > 0:
		# If this is True, then the value of the function does not change sign over the interval.
		# The root finding algorithm is not guaranteed to work. 
		print("The function does not pass through zero for the given interval.")
		return None

	num_iterations = 100000
	for i in range(num_iterations):
		middle_interval = (x_left + x_right) / 2
		middle_value = f(middle_interval)
		if abs(middle_value) <= epsilon:
			# middle_value is within range of epsilon
			print(f"Found root in {i+1} iterations.")
			return middle_interval
		elif f(x_left) * f(middle_interval) < 0:
			# f(x_left) and f(middle_interval) have opposite signs. 
			# There must be a zero between them
			# Set right value to the middle value
			x_right = middle_interval
		elif f(x_right) * f(middle_interval) < 0:
			# f(x_right) and f(middle_interval) have opposite signs. 
			# There must be a zero between them
			# Set left value to the middle value
			x_left = middle_interval
		else:
			# All signs are equal, something went wrong
			return None

	# Did not find a root within the number of iterations
	print(f"Did not find a root within {num_iterations} iterations")
	return middle_interval

## Problem 2.4 - A Physical Application: Projectile Range Maximization

################################################################################

# Problem Group 3 - Fun with primes

## Problem 3.1 - Sieve of Eratosthenes
## Fun Fact, I wrote a sieve function in Physics 18L. I will use that same implementation here. 
def sieve(limit):
	"""
	Prime number generator using Sieve of Erathosthenes.

	Args:
		limit (int): Generates prime numbers less than the limit

	Yields:
		int: n'th Prime Number, starting from 2
	"""
	#Implement Sieve of Erathosthenes
	#https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
	nums = [True] * (limit + 1) #Create a list of "True" items corresponding to the first "limit" integers.  
	nums[0] = nums[1] = False #0 and 1 are not prime numbers. 

	for (i, is_prime) in enumerate(nums):
		if is_prime:
			yield i #If the i'th number is prime, yield that number to the returned list. 
			for n in range(i * i, limit + 1, i): 
				nums[n] = False 
				#Sift out all multiples of i, as those can never be prime. 
				#All numbers less than i**2 have already been checked. 

def sieve_Eratosthenes(limit):
	"""
	Return a list of all prime numbers less than limit

	Args:
		limit (int): Maximum value in the list

	Returns:
		list: Returns list of all prime numbers less than limit
	"""
	return list(sieve(limit))

## Problem 3.2 - Prime Factorization
from math import ceil, sqrt
def prime_factors(n):
	"""
	Return a list of all prime factors of n.

	Uses sieve generator to generate prime numbers to check divisibility against. 

	Args:
		n (int): Number to find prime factors of

	Returns:
		list: Returns a list of all prime factors of n. If n has no prime factors, return a list containing only n.
	"""
	factors = []
	# The largest possible prime factor of n is sqrt(n)
	for p in sieve(ceil(sqrt(n))):
		while n % p == 0:
			factors.append(p)
			n /= p
		if n == 1:
			break
	if not factors:
		return [n]
	return factors

################################################################################

# Problem Group 4 - weighty problem
def weight_set(total):
	"""
	Inductive approach to solving The Weight Problem of Bachet de Meziriac.

	Args:
		total (int): Total weight of the initial weight.

	Returns:
		list: Returns list of pieces that can weigh any integer less than or equal to total.
	"""
	pieces = [1]

	while sum(pieces) < total:
		pieces.append(2 * sum(pieces) + 1)

	return pieces