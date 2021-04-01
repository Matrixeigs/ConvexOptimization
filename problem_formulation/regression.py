"""
Regression problem formulation using random sample

"""

from numpy import random


class RegressionProblem():
    def run(self, m, n):
        """

        :param m: the number of samples
        :param n: the dimension of features for each sample
        :return:
        """
        # generate a random matrix, with m samples and
        random.seed(0)
        A = random.random((m, n))
        b = random.random((m, 1))

        return A, b




