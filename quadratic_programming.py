"""
Quadratic programming problem test test

The standard quadratic programming prolbem is depicted as follows

"""
import numpy as np
from matplotlib import pyplot as plt
from problem_formulation.regression import RegressionProblem

regression_problem = RegressionProblem()

# Generate the samples
A, b = regression_problem.run(m=10, n=5)

# Compute the coefficient
inv_A = np.linalg.inv(np.matmul(A.transpose(), A))
x = np.matmul(np.matmul(inv_A, A.transpose()), b)

result = np.matmul(A, x) - b

print(result)

plt.plot(result)
plt.plot(b)
plt.show()
# print(b)
