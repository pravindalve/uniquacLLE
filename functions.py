import numpy as np


def uniquac_gamma(T, ind, x, a, q, r):
    x = np.array(x)
    a = np.array(a)
    q = np.array(q)
    r = np.array(r)
    tau = np.exp(- a/T)
    z = 10.0
    l = (r - q) * (z/2) - (r - 1)
    theta = (x * q) / np.sum(np.dot(x, q))
    phi = (x * r) / np.sum(np.dot(x, r))
    lngc = np.log(phi[ind] / x[ind]) + (z / 2) * q[ind] * np.log(theta[ind] / phi[ind])  \
           + l[ind] - (phi[ind] / x[ind]) * (np.sum(x * l))

    lngr = q[ind] * (1 - np.log(np.sum(theta * tau[:, ind]))
                      - np.sum((theta * tau[ind, :]) / np.matmul(theta, tau)))
    lng = lngc + lngr
    return np.exp(lng)


def uniquac_delG_mix(T, x, a, q, r):

    x = np.array(x)
    gamma = np.empty(len(x))
    for i in range(len(x)):
        gamma[i] = uniquac_gamma(T, i, x, a, q, r)
    R = 8.314472
    delG_mix = R * T * (np.sum(x * np.log(x)) + np.sum(x * np.log(gamma)))

    if delG_mix is np.nan:
        delG_mix = np.inf

    return delG_mix


def der_bin_uniquac_delG_mix(x, T, a, q, r):
    # note that here x is the mole fraction of first compound only and not an array
    h = 1e-7
    return (uniquac_delG_mix(T,[x+h, 1-x-h], a, q, r) - uniquac_delG_mix(T,[x-h, 1-x+h], a, q, r)) / (2*h)


def der_uniquac_delG_mix(x, i, T, a, q, r):
    # here input is array of x and i is the index along which the slope is required
    h = 1e-7
    xf = np.empty(len(x))
    xb = np.empty(len(x))
    for s in range(len(x)):
        if s == i:
            xf[s] = x[s] + h
            xb[s] = x[s] - h
        else:
            xf[s] = x[s]
            xb[s] = x[s]
    return (uniquac_delG_mix(T, xf, a, q, r) - uniquac_delG_mix(T, xb, a, q, r)) / (2*h)





def flatten_array(a, asize, dsize):
    f = np.empty(dsize * asize * (asize - 1))
    k = 0

    # commented code can be used when each a needs to be modified
    # for m in range(dsize):
    #     for i in range(asize):
    #         for j in range(asize):
    #             if i != j:
    #                 f[k] = a[m][i][j]
    #                 k += 1

    for i in range(asize):
        for j in range(asize):
            if i != j:
                f[k] = a[i][j]
                k += 1

    return f


def revive_array(a, asize, dsize):
    r = np.empty((dsize, asize, asize))
    m = 0
    for k in range(dsize):
        for i in range(asize):
            for j in range(asize):
                if i == j:
                    r[k][i][j] = 0
                else:
                    r[k][i][j] = a[m]
                    m += 1

    return r
