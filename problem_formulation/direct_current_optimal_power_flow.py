"""
The AC OPF power flow problem formulation based on Pypower

"""

from numpy import find_common_type
import numpy as np

from pypower.idx_gen import GEN_BUS, PG, QG, PMAX, PMIN
from pypower.idx_bus import BUS_I, BUS_TYPE, PD, QD, VMAX, VMIN
from pypower.idx_brch import F_BUS, T_BUS, RATE_A, BR_X, BR_R

from cvxopt.modeling import variable, constraint, op, dot
from pypower.loadcase import loadcase  # import the functions
from cvxopt.solvers import qp  # import the quadratic programming solution method
# from cvxopt.solvers import  # import the quadratic programming solution method
from cvxopt import matrix, mul

from pypower.ext2int import ext2int
from numpy import sum, ones, zeros


class RunDCOpf():
    def __init__(self, filename):
        self.mpc = loadcase(filename)

    def main(self):
        mpc = self.mpc
        # mpc=ext2int(mpc)
        bus, gen, branch, gencost = mpc["bus"], mpc["gen"], mpc["branch"], mpc["gencost"]
        ## Step 1: define the variables ##
        nb = bus.shape[0]
        ng = gen.shape[0]
        nl = branch.shape[0]
        # modify the indexes of buses and branches
        bus[:, BUS_I] -= 1
        branch[:, F_BUS] -= 1
        branch[:, T_BUS] -= 1

        # Convert the data format
        bus_matrix = matrix(bus)
        gen_matrix = matrix(gen)
        branch_matrix = matrix(branch)
        gencost_matrix = matrix(gencost)

        pg = variable(size=ng, name='pg')
        theta = variable(size=nb, name='theta')
        pij = variable(size=nl, name='pij')

        ## Step 2: define the objective function##
        obj = 0
        for i in range(ng):
            obj += gencost_matrix[i, 4] * mul(pg[i], pg[i])+ gencost_matrix[i, 5] * pg[i]

        ## Step 3: define the constraints ##
        cons = []
        xij = matrix(branch[:, BR_X])
        # 3.1) KCL equation
        for i in range(nb):  # The power balance constraint
            # Find the generators located at bus i
            gen_id = np.where(gen[:, GEN_BUS] == i)[0].tolist()
            branch_f_id = np.where(branch[:, F_BUS] == i)[0].tolist()
            branch_t_id = np.where(branch[:, T_BUS] == i)[0].tolist()
            if len(gen_id) > 0:
                if len(branch_f_id) > 0:
                    if len(branch_t_id) > 0:
                        cons.append(
                            sum(pg[gen_id]) - sum(pij[branch_f_id]) + sum(pij[branch_t_id]) == bus_matrix[i, PD])
                    else:
                        cons.append(sum(pg[gen_id]) - sum(pij[branch_f_id]) == bus_matrix[i, PD])
                else:
                    if len(branch_t_id) > 0:
                        cons.append(
                            sum(pg[gen_id]) + sum(pij[branch_t_id]) == bus_matrix[i, PD])
                    # the rest case is the nodes are isolated.
            else:
                if len(branch_f_id) > 0:
                    if len(branch_t_id) > 0:
                        cons.append(- sum(pij[branch_f_id]) + sum(pij[branch_t_id]) == bus_matrix[i, PD])
                    else:
                        cons.append(- sum(pij[branch_f_id]) == bus_matrix[i, PD])
                else:
                    if len(branch_t_id) > 0:
                        cons.append(sum(pij[branch_t_id]) == bus_matrix[i, PD])
                        # the rest case is the nodes are isolated.
        # 3.2) relation between power flow and angle
        for i in range(nl):
            cons.append(pij[i] == (theta[int(branch[i, F_BUS])] - theta[int(branch[i, T_BUS])]) / xij[i])
        # 3.3) relation of boundaries
        cons.append(pg >= gen_matrix[:, PMIN])
        cons.append(pg <= gen_matrix[:, PMAX])
        cons.append(theta >= -360)
        cons.append(theta <= 360)
        cons.append(pij >= -branch_matrix[:, RATE_A])
        cons.append(pij <= branch_matrix[:, RATE_A])

        # 4) Formulate the problem
        qp = op(obj, cons)
        # 5) Solve the problem
        qp.solve()
        # 6) Return the values
        pg = pg.value
        pij = pij.value
        theta = theta.value

        return pg, pij, theta


if __name__ == "__main__":
    from pypower.case30 import case30

    run_dc_opf = RunDCOpf(case30())
    run_dc_opf.main()
