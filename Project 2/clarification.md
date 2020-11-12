Helpful Clarification on Eigen Analysis Step #2. 

I had some confusion on the form of solution that #2 is asking for. Thankfully, Pratik was able clarify this step. Here is my best interpretation of it, which may help you. 

Ultimately, we want to compute and plot the wave function of the anharmonic oscillator for the first four energy levels. That is, $$\psi_0(x)$$, $$\psi_1(x)$$, $$\psi_2(x)$$, and $$\psi_3(x)$$.

In order to solve the Hamiltonian, we must choose a value for $$\lambda$$ in $$\hat{H}(\lambda) = \hat{H} + \lambda \hat{x^4}$$, and a size $$N$$ for the Hermitian matrix. I will choose $$\lambda = 1$$, and set $$N$$ to be large, like $$N = 20$$. 

The Hermitain matrix $$\hat{H}(\lambda)$$ can be easily computed. We then diagonalize this with our `hermitian_eigensystem` method. We then get a vector of 20 eigenvalues, and a matrix with 20 columns, the eigenvectors corresponding to the eigenvalues. 

Consider the ground state $$\psi_0(x)$$. From our `hermitian_eigensystem` solution, the first element in the eigenvalue vector is the energy of the ground state, $$E_0(\lambda)$$. The first column of the eigenvector matrix is the eigenvector that corresponds with this eigenvalue for $$E_0(\lambda)$$. Let $$c = [c_0, c_1, \cdots,c_n]$$ be this eigenvector. 

We then say that for the anharmonic oscillator, 

$$\psi_0(x)= c_0 \phi_0(x) + c_1 \phi_1(x) + \dotsc, + c_n \phi_n(x)$$ 

where 

$$\phi _n(x) = (2^n n! \sqrt{\pi})^{-1/2}e^{-x^2/2}H_n(x)$$ 

as given on page 15 of `project2-lecture notes 1.pdf`. Since we chose a matrix size of $$N = 20$$, our $$n$$ ranges from 0 to 19. So, 


$$\psi_0(x)= c_0 \phi_0(x) + c_1 \phi_1(x) + \dotsc, + c_{19} \phi_{19}(x)$$ 

We can then plot the Real part of this wave function as a function of $$x$$. This is the wave function solution for the first energy level. We can then repeat this process for the next three energy levels, for $$\psi_0(x)$$, $$\psi_1(x)$$, $$\psi_2(x)$$, and $$\psi_3(x)$$.