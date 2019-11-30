import numpy as np
import sys
sys.path.insert(0, 'C:/Users/FOSSEE/PycharmProjects/uniqacLLE')
from functions import  der_uniquac_delG_mix, uniquac_delG_mix, uniquac_gamma
from scipy.optimize import fsolve
import scipy

T = 298.15
xIexp = [0.8487, 0.0222, 0.1291]
xIIexp = [0.1866, 0.0328, 0.7806]
xIcalc = [0.876101659470136,0.014410918170012,0.102359852]
xIIcalc = [0.191460745339455,0.039411000363692,0.769128254296853]
acalc = [[0 , -384.913, 10.07412],[210.9491, 0, -187.338], [628.0464, -165.228, 0]]
q = [1.432, 2.4, 4.396]
r = [1.4311, 3.1878, 5.1742]
Nc = len(xIexp)

derGmixI = np.empty(Nc - 1)
derGmixII = np.empty(Nc - 1)

for i in range(Nc - 1):
    derGmixI[i] = der_uniquac_delG_mix(xIcalc, i, T, acalc, q, r)
    derGmixII[i] = der_uniquac_delG_mix(xIIcalc, i, T, acalc, q, r)

print(derGmixI, derGmixII)


def comTanPlane(x, xp):
	T = 298.15
	acalc = [[0 , -384.913, 10.07412],[210.9491, 0, -187.338], [628.0464, -165.228, 0]]
	q = [1.432, 2.4, 4.396]
	r = [1.4311, 3.1878, 5.1742]
	comTanP = uniquac_delG_mix(T, xIcalc, acalc, q, r)\
					 + (der_uniquac_delG_mix(xIcalc, 0, T, acalc, q, r) * (x[0] - xp[0]))\
					 + (der_uniquac_delG_mix(xIcalc, 1, T, acalc, q, r) * (x[1] - xp[1]))

	return comTanP



def minDist(xp):
	T = 298.15
	acalc = [[0 , -384.913, 10.07412],[210.9491, 0, -187.338], [628.0464, -165.228, 0]]
	q = [1.432, 2.4, 4.396]
	r = [1.4311, 3.1878, 5.1742]
	x1 = x2 = np.arange(0.01, 0.99, 0.01)
	d = np.empty((98, 98))
	for i in range(98):
		for j in range(98):
			if (x1[i] + x2[j]) >= 1:
				d[i, j] = np.inf
			else:
				d[i, j] = uniquac_delG_mix(T, [x1[i], x2[j], 1- x1[i] - x2[j]], acalc, q, r) - comTanPlane([x1[i], x2[j], 1- x1[i] - x2[j]], xp)


	dmin = np.amin(d)
	return dmin


s = minDist(xIcalc)
print(s)




# this is copied from test2lle code for solving equations
def equations1(m):
    xIe = np.array(list(m[:2])+ list([1- sum(m[:2])]))
    # print(xIe)
    xIIe = np.array(list(m[2:4])+ list([1- sum(m[2:4])]))
    # print(xIIe)
    bta = m[-1]
    # print((np.sum(xIe + xIIe)))
    ze = (xIe + xIIe ) / 2
    # print(ze)
    gmaI = np.array([1.07525714759282, 0.727282395860607, 6.61872966901385])
    gmaII = np.array([4.7573733448669, 0.307910263657295, 1.12422266744842])
    # 1: x = [0.8487, 0.0222, 0.1291], 2: x = [0.1866, 0.0328, 0.7806]
    amm = [[0 , -384.913, 10.07412],[210.9491, 0, -187.338], [628.0464, -165.228, 0]]

    for i in range(len(xIe)):
        gmaI[i] = uniquac_gamma(T = 298.15,ind=i,x = xIe,
                                a = amm,
                                q = [1.432, 2.4, 4.396],r = [1.4311, 3.1878, 5.1742])
        gmaII[i] = uniquac_gamma(T = 298.15,ind=i,x = xIIe,
                                a = amm,
                                q = [1.432, 2.4, 4.396],r = [1.4311, 3.1878, 5.1742])

    ls = []
    # print(gmaI, gmaII)
    for i in range(len(xIe)):
        ls.append(xIe[i] * gmaI[i] - xIIe[i] * gmaII[i])

    for i in range(len(xIe) - 1):
        ls.append(bta*xIe[i] + (1 - bta)*xIIe[i] - ze[i])
    #
    # ls.append(np.sum(xIe) - 1)
    # ls.append(np.sum(xIIe) - 1)
    if any(np.isnan(ls[i]) for i in range(5)):
        print('tufan')
        ls = np.array([1, 1, 1, 1, 1])
    # print(ls)
    return ls




s = 1

dist = 1
while abs(dist) > 1e-4:
	x1e = np.full(3,-1)
	x2e = np.full(3,-1)
	diff1 = np.full(2, -1)
	diff2 = np.full(2, -1)
	xie = np.full(5, -1)
	while any(xie<=0) or any(xie >= 1) or (abs(x1e[0] - x2e[0]) < 1e-2) or diff1 > diff2:
	    xIe1 = np.clip(np.random.ranf(), max(0.0,0.7), 0.9)
	    xIe2 = np.clip(np.random.ranf(), 0, 0.3)
	    xIe3 = 1 - xIe1 - xIe2
	    # xIe3 = np.random.ranf()
	    # xIe3 = -1
	    xIIe1 =np.clip(np.random.ranf(), 0, 0.2)
	    xIIe2 = np.clip(np.random.ranf(), 0, 0.1)
	    xIIe3 = 1 - xIIe1 - xIIe2
	    # xIIe3 = np.random.ranf()
	    # 1: x = [0.8487, 0.0222, 0.1291], 2: x = [0.1866, 0.0328, 0.7806]
	    # xIe1, xIe2, xIe3 = [0.8, 0.02, 0.1]
	    # xIIe1, xIIe2, xIIe3 = [0.1, 0.01, 0.7]


	    bta = np.random.ranf()
	    inp = np.array([xIe1, xIe2, xIIe1, xIIe2, bta])
	    # print(xIe1, xIe2, xIe3, xIIe1, xIIe2, xIIe3, bta)
	    if all(0 < inp) and all(inp < 1) and 0 < (1 - np.sum(inp[:2])) < 1 and 0 < (1 - np.sum(inp[2:4])) <1:
	        # print(xIe1, xIe2,xIIe1, xIIe2, bta)
	        try:
	            xie = fsolve(equations1, (xIe1, xIe2, xIIe1, xIIe2, bta))
	            print(xie)
	        except scipy.optimize.nonlin.NoConvergence:
	            print('not converged')
	        except OverflowError:
	            print('overflow')
	        except RuntimeWarning:
	            pass

	    diff1 = np.linalg.norm(np.array(xie[:2]) - np.array([0.87, 0.01]))
	    diff2 = np.linalg.norm(np.array(xie[2:4]) - np.array([0.87, 0.01]))
	    x1e = list(xie[:2]) + list([1 - np.sum(xie[:2])])
	    x2e = list(xie[2:4]) + list([1 - np.sum(xie[2:4])])


	    # print(s)
	    s+=1
	    # print(s)
	print(x1e, x2e)
	dist1 = minDist(x1e)
	print(dist1)
s1 = State(T, x1e,x2e, acalc,q,r)
s1.check_CTP()
print(s1.CTP)





