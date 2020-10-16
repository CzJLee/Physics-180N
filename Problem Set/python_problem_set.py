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

## Problem 2.2 - Bisection Search

## Problem 2.3 - Bisection Root Finding

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