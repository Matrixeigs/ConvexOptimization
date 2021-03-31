"""
Test cases for linear programming based on CVXOPT

This case is an economic dispatch code
"""

from cvxopt import matrix
from cvxopt import solvers


class LinearProgramming():
    def run(self, PG_MAX, PG_MIN, CG, PD):
        """

        :param PG_MAX: maximal generator output
        :param PG_MIN: minimal generator output
        :param PD: total demand
        :return:
        """
        ## Step 1: Physical model ##
        # 1) PG_MIN <= pg, \forall g
        # 2) pg <= PG_MAX, \forall g
        # 3) sum(pg) = PD
        # obj: sum_{g} CG_{g}*pg_{g}
        ## Step 2: Compact Model ##
        ## min_{x} c^{T}x
        ## s.t. Gx <= h
        ##      Ax = b
        nx = len(PG_MAX)
        ng = len(PG_MIN)
        assert len(PG_MAX) == len(PG_MIN)
        # 1) PG_MIN <= pg, \forall g
        # -pg_{g} <= -PG_MIN_{g}, \forall g
        G = matrix([0.0] * ng * nx, (ng, nx))
        h = matrix([0.0] * ng, (ng, 1))
        for i in range(ng):
            G[i, i] = -1
            h[i, 0] = -PG_MIN[i]
        # 2) pg <= PG_MAX, \forall g
        # pg_{g} <= PG_MAX_{g}, \forall g
        G_temp = matrix([0.0] * ng * nx, (ng, nx))
        h_temp = matrix([0.0] * ng, (ng, 1))
        for i in range(ng):
            G_temp[i, i] = 1
            h_temp[i, 0] = PG_MAX[i]

        G = matrix([G, G_temp])
        h = matrix([h, h_temp])
        # 3) sum(pg) = PD
        A = matrix([1.0] * ng, (1, nx))
        b = matrix([PD])
        # 4) Objective function
        c = matrix(CG, (nx, 1))
        sol = solvers.lp(c=c, G=G, h=h, A=A, b=b)

        pg = sol['x']

        return pg


if __name__ == "__main__":
    linear_programming = LinearProgramming()
    linear_programming.run(PG_MAX=[10., 10., 10.], PG_MIN=[0., 0., 0.], CG=[1., 2., 3.], PD=24.)
