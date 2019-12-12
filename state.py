from itertools import combinations
import numpy as np
import scipy
from scipy.optimize import newton, fsolve, minimize
from functions import der_bin_uniquac_delG_mix, uniquac_gamma, uniquac_delG_mix, der_uniquac_delG_mix


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
        self.xIcalc = self.xIexp
        self.xIIcalc = self.xIIexp
        self.xIiter = self.xIexp
        self.xIIiter = self.xIIexp
        self.LLE = None
        self.no_of_comp = len(xIexp)


    def update_LLE(self):
        #updates boolean variable LLE by checking existence of LLE 
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


    def equations(self, m):
        #returns list of equations for ternary system 
        xIe = np.array(list(m[:2])+ list([1- sum(m[:2])]))
        xIIe = np.array(list(m[2:4])+ list([1- sum(m[2:4])]))
        bta = m[-1]
        ze = (xIe + xIIe ) / 2
        gmaI = np.empty(self.no_of_comp)
        gmaII = np.empty(self.no_of_comp)

        for i in range(self.no_of_comp):
            gmaI[i] = uniquac_gamma(self.Texp, i, xIe, self.a, self.q, self.r)
            gmaII[i] = uniquac_gamma(self.Texp, i, xIIe, self.a, self.q, self.r)

        ls = xIe * gmaI - xIIe * gmaII
        ls = np.append(ls, bta*xIe[:len(xIe) - 1] + (1 - bta)*xIIe[:len(xIe) - 1] - ze[:len(xIe) - 1])

        if any(np.isnan(ls[i]) for i in range(5)):
            ls = np.array([1, 1, 1, 1, 1])
        return ls

    def update_xcalc(self):
        # updates xcalc for both the phases of ternary system
        xIe = np.full(3,-1)
        xIIe = np.full(3,-1)
        diff1 = np.full(2, -1)
        diff2 = np.full(2, -1)
        xie = np.full(5, -1)
        while any(xie<=0) or any(xie >= 1) or any(abs(xIe - xIIe) < 1e-2) or diff1 > diff2:
            xIe1 = np.clip(np.random.ranf(), max(0.0,self.xIcalc[0] - 0.1), min(self.xIcalc[0], 0.999))
            xIe2 = np.clip(np.random.ranf(), max(0.0,self.xIcalc[1] - 0.1), min(self.xIcalc[1], 0.999))
            xIe3 = 1 - xIe1 - xIe2
            xIIe1 =np.clip(np.random.ranf(), max(0.0,self.xIIcalc[0] - 0.1), min(self.xIIcalc[0], 0.999))
            xIIe2 = np.clip(np.random.ranf(), max(0.0,self.xIIcalc[1] - 0.1), min(self.xIIcalc[1], 0.999))
            xIIe3 = 1 - xIIe1 - xIIe2

            bta = np.random.ranf()
            inp = np.array([xIe1, xIe2, xIIe1, xIIe2, bta])
            # print(inp)
            if all(0 < inp) and all(inp < 1) and 0 < (1 - np.sum(inp[:2])) < 1 and 0 < (1 - np.sum(inp[2:4])) <1:
                try:
                    xie = fsolve(self.equations, (xIe1, xIe2, xIIe1, xIIe2, bta), xtol = 1e-10)
                    # print(xie)
                except scipy.optimize.nonlin.NoConvergence:
                    print('not converged')
                except OverflowError:
                    print('overflow')

            diff1 = np.linalg.norm(np.array(xie[:2]) - np.array(self.xIexp[:2]))
            diff2 = np.linalg.norm(np.array(xie[2:4]) - np.array(self.xIIexp[:2]))
            xIe = np.array(list(xie[:2]) + list([1 - np.sum(xie[:2])]))
            xIIe = np.array(list(xie[2:4]) + list([1 - np.sum(xie[2:4])]))

        self.xIcalc = xIe
        self.xIIcalc = xIIe


    def tan_dist(self, x, xp):
        if (x[0] + x[1]) >= 1:
            d = np.inf
            print("this is not good")
        else:
            d = abs(uniquac_delG_mix(self.Texp, x, self.a, self.q, self.r) - self.com_tan_plane(x, xp))

        return d

    def com_tan_plane(self, x, xp):
        
        com_tan_p = uniquac_delG_mix(self.Texp, xp, self.a, self.q, self.r)\
                         + (der_uniquac_delG_mix(xp, 0, self.Texp, self.a, self.q, self.r) * (x[0] - xp[0]))\
                         + (der_uniquac_delG_mix(xp, 1, self.Texp, self.a, self.q, self.r) * (x[1] - xp[1]))

        return com_tan_p


    def min_tan_dist(self, xie):
        def sumrestrict(x):
            return np.sum(x) - 1

        my_constraint = ({'type':'eq', 'fun':sumrestrict})
    
        distres = scipy.optimize.minimize(self.tan_dist, xie, method='SLSQP', args = (xie), options={'disp': False}, tol = 1e-10, 
                                    bounds = tuple(zip(np.maximum(xie - 0.1, 0.001),np.minimum(xie + 0.1, 0.999))), constraints = my_constraint)
    
        return distres

