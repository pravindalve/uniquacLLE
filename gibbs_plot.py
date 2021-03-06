import numpy as np
from scipy.optimize import minimize
from functions import uniquac_gamma, uniquac_delG_mix, der_uniquac_delG_mix
from uniquaclle import UniquacLLE
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D


Texp = [298.15, 298.15, 298.15, 298.15, 298.15, 298.15, 298.15, 298.15]
xIexp = [[0.8487, 0.0222, 0.1291], [0.8319, 0.0295, 0.1386], [0.815, 0.0368, 0.1482], [0.7986, 0.0439, 0.1575],
         [0.775, 0.0541, 0.1709], [0.7665, 0.0578, 0.1757], [0.7506, 0.0646, 0.1848], [0.7285, 0.074, 0.1975]]
xIIexp = [[0.1866, 0.0328, 0.7806], [0.2162, 0.058, 0.7258], [0.2394, 0.0748, 0.6858], [0.2568, 0.0863, 0.6569],
          [0.2817, 0.1004, 0.6179], [0.2927, 0.1059, 0.6014], [0.3106, 0.1137, 0.5757], [0.3382, 0.1231, 0.5387]]

a = np.array([[0 , -384.913, 10.07412],[210.9491, 0, -187.338], [628.0464, -165.228, 0]])
q = [1.432, 2.4, 4.396]
r = [1.4311, 3.1878, 5.1742]

for i in range(1):
    x = np.arange(0.01, 0.99, 0.01)
    y = np.arange(0.01, 0.99, 0.01)

    z = np.empty((98, 98))
    p = np.empty((98, 98))

    for j in range(98):
        for k in range(98):
            if (1-x[j] - y[k]) < 0:
                z[j, k] = 0
            else:
                z[j, k] = uniquac_delG_mix(Texp[i], [x[j], y[k], 1-x[j] - y[k]], a, q, r)
    x, y = np.meshgrid(x, y)
    z[z == 0.] = np.nan
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.contour3D(x, y, z, 200)
    plt.show()
