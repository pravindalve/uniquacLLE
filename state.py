from itertools import combinations
import numpy as np
from scipy.optimize import newton
from functions import der_bin_uniquac_delG_mix


class State():

    def __init__(self, Texp, xIexp, xIIexp, aguess, q, r):
        self.Texp = Texp
        self.xIexp = xIexp
        self.xIIexp = xIIexp
        self.a = aguess
        self.q = q
        self.r = r
        self.gI = None
        self.gII = None
        self.xIcalc = None
        self.xIIcalc = None
        self.xIiter = self.xIexp
        self.xIIiter = self.xIIexp
        self.LLE = None


    def update_LLE(self):
        combs = combinations(np.arange(0, len(self.xIexp), 1), 2)
        ls = list(combs)
        bool_LLE = np.full(len(self.xIexp), None)
        for comb in ls:
            try:
                rt1 = newton(der_bin_uniquac_delG_mix, 0.01, maxiter = 10000, args = (self.Texp,
                        np.array([[0, self.a[comb[0]][comb[1]]], [self.a[comb[1]][comb[0]], 0]]),
                        np.array([self.q[comb[0]], self.q[comb[1]]]), np.array([self.r[comb[0]], self.r[comb[1]]])))
            except RuntimeError:
                rt1 = None

            try:
                rt2 = newton(der_bin_uniquac_delG_mix, 0.99, maxiter = 10000, args = (self.Texp,
                        [[0, self.a[comb[0]][comb[1]]], [self.a[comb[1]][comb[0]], 0]], [self.q[comb[0]], self.q[comb[1]]],
                                                                                      [self.r[comb[0]], self.r[comb[1]]]))
            except RuntimeError:
                rt2 = None

            if rt1 is not None and rt2 is not None and abs(rt1 - rt2) < 1e-4:
                bool_LLE[ls.index(comb)] = False
            elif rt1 is not None and rt2 is not None and rt1 != rt2:
                bool_LLE[ls.index(comb)] = True
            else:
                bool_LLE[ls.index(comb)] = False
                # because of this condition if there is inconsistency in BIP and root not found by newton method both are shown as miscible.
                # when inconsistency is handled just make change in this condition

        if True in bool_LLE:
            self.LLE = True
        else:
            self.LLE = False

