import numpy as np
from scipy.optimize import minimize
from functions import uniquac_gamma, uniquac_delG_mix
from uniquaclle import UniquacLLE
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D


Texp = [298.15, 298.15, 298.15, 298.15, 298.15, 298.15, 298.15, 298.15]
xIexp = [[0.8487, 0.0222, 0.1291], [0.8319, 0.0295, 0.1386], [0.815, 0.0368, 0.1482], [0.7986, 0.0439, 0.1575],
         [0.775, 0.0541, 0.1709], [0.7665, 0.0578, 0.1757], [0.7506, 0.0646, 0.1848], [0.7285, 0.074, 0.1975]]
xIIexp = [[0.1866, 0.0328, 0.7806], [0.2162, 0.058, 0.7258], [0.2394, 0.0748, 0.6858], [0.2568, 0.0863, 0.6569],
          [0.2817, 0.1004, 0.6179], [0.2927, 0.1059, 0.6014], [0.3106, 0.1137, 0.5757], [0.3382, 0.1231, 0.5387]]

a = np.array([[0, -87.1325, 17.6775],[501.439, 0, -232.51],[631.636, 338.947, 0]])


# a = np.array([[0, 0.2002, 2.980], [-3.003, 0, 9.27], [10.54, -1.001, 0]])
# a.fill(100)
q = [1.432, 2.4, 4.396]

r = [1.4311, 3.1878, 5.1742]
# print(Texp[0], xIexp[0], xIIexp[0])
t = [298.15]
xi = [xIexp[0]]
xii = [xIIexp[0]]
# dataset = UniquacLLE(t, xIexp, xIIexp, a, q, r)
dataset = UniquacLLE(t, xi, xii, [[0,1,1], [1, 0, 1], [1, 1, 0]], q, r)

# value = dataset.optimize_OF1()
# print(dataset.a_nm)

del_g = np.empty(8)
for i in range(8):
    del_g[i] = uniquac_delG_mix(Texp[i], xIexp[i], a, q, r)
print(del_g)

Tb = [253.15]
# ab = [[0, 21.81246], [684.4026, 0]]
qb = [1.432, 4.396]
rb = [1.4311, 5.1742]


dg = uniquac_delG_mix(253.15, [0.9563, 0.0437], [[0, 1],[1, 0]], [1.432, 4.396], [1.4311, 5.1742])
print(dg)
gibthree = uniquac_delG_mix(Texp[0], xIexp[0], a, q, r)

for i in range(8):
    x = np.arange(0.01, 0.99, 0.01)
    y = np.arange(0.01, 0.99, 0.01)

    z = np.empty((98, 98))

    for j in range(98):
        for k in range(98):
            if (1-x[j] - y[k]) < 0:
                z[j, k] = 0
            else:
                z[j, k] = uniquac_delG_mix(Texp[i], [x[j], y[k], 1-x[j] - y[k]], a, q, r)
    # plt.title("Gibbs free energy for datapoint " + str(i))
    #
    # # Plot the points using matplotlib
    # plt.plot(x, y)
    # plt.show()

    x, y = np.meshgrid(x, y)
    z[z == 0.] = np.nan
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.contour3D(x, y, z, 100)
    plt.show()
