import numpy as np
import functions
from state import State
from functions import uniquac_gamma
from scipy.optimize import minimize
from itertools import combinations

class UniquacLLE():

    def __init__(self, Texp, xIexp, xIIexp, a, q, r):
        self.Texp = Texp
        self.xIexp = xIexp
        self.xIIexp = xIIexp
        self.aguess = a
        self.q = q
        self.r = r
        self.no_of_datapts = len(Texp)
        self.no_of_comp = len(xIexp[0])
        self.states = self.initialize_states()
        self.discarded_states = None
        self.a_nm = self.aguess   # nelder-mead a

    def initialize_states(self):
        states = []
        # states = np.empty(self.no_of_datapts, dtype=object)
        for i in range(self.no_of_datapts):
            states.append(State(self.Texp[i], self.xIexp[i], self.xIIexp[i], self.aguess, self.q, self.r))
        return states

    def OF1(self, s):
        s = functions.revive_array(s, self.no_of_comp, 1)
        # print(s)
        s = np.reshape(s, (self.no_of_comp, self.no_of_comp))

        a = np.empty((self.no_of_datapts, self.no_of_comp, self.no_of_comp))
        for i in range(self.no_of_datapts):
            a[i][:][:] = s[:][:]


        gammaI = np.zeros((self.no_of_datapts, self.no_of_comp))
        gammaII = np.zeros((self.no_of_datapts, self.no_of_comp))
        # print(a)
        for i in range(self.no_of_datapts):
            for j in range(self.no_of_comp):
                gammaI[i, j] = uniquac_gamma(self.Texp[i], j, self.states[i].xIiter, a[i], self.states[i].q,
                                             self.states[i].r)
                gammaII[i, j] = uniquac_gamma(self.Texp[i], j, self.states[i].xIIiter, a[i],
                                              self.states[i].q, self.states[i].r)
        OF1 = 0
        for i in range(self.no_of_datapts):
            OF1 = OF1 + np.absolute(
                np.sum((self.states[i].xIiter * gammaI[i, :]) - (self.states[i].xIIiter * gammaII[i, :])))

        self.a_nm = a#checkout if a_nm can be eliminated
        print(OF1)
        return OF1

    def optimize_OF1(self):
        a = functions.flatten_array(self.states[0].a, self.no_of_comp, 1)
        opt = minimize(self.OF1, a, method='nelder-mead', options={'xtol': 1e-4, 'disp': True, 'maxiter':10000})
        for i in range(self.no_of_datapts):
            self.states[i].a = self.a_nm[i]
        return opt

    def OF2(self):
        sum = 0 
        for s in self.states:
            sum = sum + (np.sum(abs(s.xIexp - s.xIcalc)) + np.sum(abs(s.xIiter - s.xIcalc)) + np.sum(abs(s.xIexp - s.xIiter)) +
                    np.sum(abs(s.xIIexp - s.xIIcalc)) + np.sum(abs(s.xIIiter - s.xIIcalc)) + np.sum(abs(s.xIIexp - s.xIIiter)))

        return sum

    def OF(self):
        return np.sqrt(np.square(self.OF1()) + np.square(self.OF2()))
        

    def evaluate_BIP(self):
        print('This function is not yet implemented')
        return self.aguess
