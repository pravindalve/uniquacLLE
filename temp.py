import numpy as np 


x1e = np.array([0.8, 0.00001, 0.99])
xmin = np.maximum(x1e - 0.01, 0.001)
xmax = np.minimum(x1e + 0.01, 0.999)
print(tuple(zip(xmin, xmax)))
# bounds = (np.max(x1e-0.01, np.full(3,0.001)), np.min(x1e + 0.01, np.full(3,0.999)))
# print(bounds)