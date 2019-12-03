import sys
sys.path.insert(0, 'C:/Users/FOSSEE/PycharmProjects/uniqacLLE')
import numpy as np
from functions import  der_uniquac_delG_mix, uniquac_delG_mix, uniquac_gamma
from scipy.optimize import fsolve
import scipy
from state import State


T = 298.15
acalc = [[0 , -384.913, 10.07412],[210.9491, 0, -187.338], [628.0464, -165.228, 0]]
q = [1.432, 2.4, 4.396]
r = [1.4311, 3.1878, 5.1742]

def comTanPlane(x, xp):
    T = 298.15
    acalc = [[0 , -384.913, 10.07412],[210.9491, 0, -187.338], [628.0464, -165.228, 0]]
    q = [1.432, 2.4, 4.396]
    r = [1.4311, 3.1878, 5.1742]
    comTanP = uniquac_delG_mix(T, xp, acalc, q, r)\
                     + (der_uniquac_delG_mix(xp, 0, T, acalc, q, r) * (x[0] - xp[0]))\
                     + (der_uniquac_delG_mix(xp, 1, T, acalc, q, r) * (x[1] - xp[1]))

    return comTanP



def Dist(x, xp):
    T = 298.15
    acalc = [[0 , -384.913, 10.07412],[210.9491, 0, -187.338], [628.0464, -165.228, 0]]
    q = [1.432, 2.4, 4.396]
    r = [1.4311, 3.1878, 5.1742]
    if (x[0] + x[1]) >= 1:
        d = np.inf
        print("this is not good")
    else:
        d = abs(uniquac_delG_mix(T, x, acalc, q, r) - comTanPlane(x, xp))
        # print('delgmix = ' + str(uniquac_delG_mix(T, x, acalc, q, r)))
        # print('comPlane = ' + str(comTanPlane(x, xp)))
        # print('x  = ' + str(x))
    return d


def sumrestrict(x):
    return np.sum(x) - 1

my_constraint = ({'type':'eq', 'fun':sumrestrict})

x1e = np.array([0.853166507572777, 0.016252662016868, 0.130580830410355])
x2e = np.array([0.192831909302274, 0.038388700754375, 0.768779389943351])

dist1res = scipy.optimize.minimize(Dist, x1e, method='SLSQP', args = (x1e), options={'disp': False}, 
                                    bounds = tuple(zip(np.maximum(x1e - 0.1, 0.001),np.minimum(x1e + 0.1, 0.999))), constraints = my_constraint)

dist2res = scipy.optimize.minimize(Dist, x2e, method='SLSQP', args = (x2e), options={'disp': False},
                                    bounds = tuple(zip(np.maximum(x2e - 0.1, 0.001),np.minimum(x2e + 0.1, 0.999))), constraints = my_constraint)

print(dist1res)
print(dist2res)
coef1 = der_uniquac_delG_mix(x1e, 0, T, acalc, q, r)
coef2 = der_uniquac_delG_mix(x2e, 0, T, acalc, q, r)
print(abs(coef1 - coef2))