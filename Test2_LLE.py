from state import State
from itertools import combinations
import numpy as np
from scipy.optimize import broyden1, fsolve
from functions import uniquac_gamma
import scipy


T = 298.15
xIexp = [0.8487, 0.0222, 0.1291]
xIIexp = [0.1866, 0.0328, 0.7806]
xIcalc = [0.876101659470136,0.014410918170012,0.102359852]
xIIcalc = [0.191460745339455,0.039411000363692,0.769128254296853]
acalc = [[0 , -384.913, 10.07412],[210.9491, 0, -187.338], [628.0464, -165.228, 0]]

a = [[0, -87.1325, 17.6775],[501.439, 0, -232.51],[631.636, 338.947, 0]]
q = [1.432, 2.4, 4.396]
r = [1.4311, 3.1878, 5.1742]

s = State(T, xIcalc,xIIcalc, acalc,q,r)

gI = np.empty(3)
gII = np.empty(3)
for i in range(3):
    gI[i] = uniquac_gamma(T, i, xIcalc, acalc, q, r)
    gII[i] = uniquac_gamma(T, i, xIIcalc, acalc, q, r)
print(gI, gII)

# def equations(m):
#     xIe = np.array(m[:3])
#     xIIe = np.array(m[3:])
#     bta = 0.5
#     # print((np.sum(xIe + xIIe)))
#     ze = (xIe + xIIe ) / 2
#     # print(ze)
#     gmaI = np.array([1.07525714759282, 0.727282395860607, 6.61872966901385])
#     gmaII = np.array([4.7573733448669, 0.307910263657295, 1.12422266744842])
#     # 1: x = [0.8487, 0.0222, 0.1291], 2: x = [0.1866, 0.0328, 0.7806]
#     amm = [[0 , -384.913, 10.07412],[210.9491, 0, -187.338], [628.0464, -165.228, 0]]
#
#     for i in range(len(xIe)):
#         gmaI[i] = uniquac_gamma(T = 298.15,ind=i,x = xIe,
#                                 a = amm,
#                                 q = [1.432, 2.4, 4.396],r = [1.4311, 3.1878, 5.1742])
#         gmaII[i] = uniquac_gamma(T = 298.15,ind=i,x = xIIe,
#                                 a = amm,
#                                 q = [1.432, 2.4, 4.396],r = [1.4311, 3.1878, 5.1742])
#
#     ls = []
#     # print(gmaI, gmaII)
#     for i in range(len(xIe) -1):
#         ls.append(xIe[i] * gmaI[i] - xIIe[i] * gmaII[i])
#
#     for i in range(1,len(xIe),1):
#         ls.append(bta*xIe[i] + (1 - bta)*xIIe[i] - ze[i])
#
#     ls.append(np.sum(xIe) - 1)
#     ls.append(np.sum(xIIe) - 1)
#
#     if (np.isnan(ls[i]) for i in range(6)):
#         ls = [1, 1, 1, 1, 1,1]
#
#     print(ls)
#     return ls


# xie = np.full(6, -1)
#
# s = 1
# gmI = np.empty(3)
# gmII = np.empty(3)
# amm = [[0 , -384.913, 10.07412], [210.9491, 0, -187.338], [628.0464, -165.228, 0]]
#
# x1e = np.full(3,-1)
# x2e = np.full(3,-1)
# diff1 = np.full(2, -1)
# diff2 = np.full(2, -1)
#
#
# while any(xie<=0) or any(xie >= 1) or (abs(x1e[0] - x2e[0]) < 1e-4) or diff1 > diff2:
#     xIe1 = np.clip(np.random.ranf(), max(0.0,0.7), 0.9)
#     xIe2 = np.clip(np.random.ranf(), 0, 0.3)
#     xIe3 = 1 - xIe1 - xIe2
#     # xIe3 = np.random.ranf()
#     # xIe3 = -1
#     xIIe1 =np.clip(np.random.ranf(), 0, 0.2)
#     xIIe2 = np.clip(np.random.ranf(), 0, 0.1)
#     xIIe3 = 1 - xIIe1 - xIIe2
#     # xIIe3 = np.random.ranf()
#     # 1: x = [0.8487, 0.0222, 0.1291], 2: x = [0.1866, 0.0328, 0.7806]
#     # xIe1, xIe2, xIe3 = [0.8, 0.02, 0.1]
#     # xIIe1, xIIe2, xIIe3 = [0.1, 0.01, 0.7]
#
#
#     bta = np.random.ranf()
#     inp = np.array([xIe1, xIe2, xIe3, xIIe1, xIIe2, xIIe3])
#
#
#     if all(inp>0) and all(inp<1):
#         # print(xIe1, xIe2, xIe3, xIIe1, xIIe2, xIIe3, bta)
#         try:
#             xie = fsolve(equations, (xIe1, xIe2, xIe3, xIIe1, xIIe2, xIIe3))
#             print(xie)
#         except scipy.optimize.nonlin.NoConvergence:
#             print('not converged')
#         except OverflowError:
#             print('overflow')
#         except RuntimeWarning:
#             pass
#
#     diff1 = np.linalg.norm(np.array(xie[:2]) - np.array([0.87, 0.01]))
#     diff2 = np.linalg.norm(np.array(xie[3:5]) - np.array([0.87, 0.01]))
#     x1e = xie[:3]
#     x2e = xie[3:6]
#
#     # print(x1e, x2e)
#     # print(s)
#     s+=1
#     # print(s)
#
# s1 = State(T, xie[:3],xie[3:6], acalc,q,r)
# s1.check_CTP()
# print(s1.CTP)

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
    print(ls)
    return ls


xie = np.full(5, -1)

s = 1
gmI = np.empty(3)
gmII = np.empty(3)
amm = [[0 , -384.913, 10.07412], [210.9491, 0, -187.338], [628.0464, -165.228, 0]]

x1e = np.full(3,-1)
x2e = np.full(3,-1)
diff1 = np.full(2, -1)
diff2 = np.full(2, -1)


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
s1 = State(T, x1e,x2e, acalc,q,r)
s1.check_CTP()
print(s1.CTP)