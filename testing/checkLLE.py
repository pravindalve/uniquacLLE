import sys
sys.path.insert(0, 'C:/Users/FOSSEE/PycharmProjects/uniqacLLE')
from state import State

s = State(298.15, [0.8487, 0.0222, 0.1291], [0.1866, 0.0328, 0.7806],
          [[0, -87.1325, 17.6775],[501.439, 0, -232.51],[631.636, 338.947, 0]],
          [1.432, 2.4, 4.396], [1.4311, 3.1878, 5.1742])

# s = State(298.15, [0.8487, 0.0222, 0.1291], [0.1866, 0.0328, 0.7806],
#           [[0, 1, 1],[1, 0, 1],[1, 1, 0]],
#           [1.432, 2.4, 4.396], [1.4311, 3.1878, 5.1742])

print(s)
print(s.LLE)
s.update_LLE()
print(s.LLE)

print(s.xIcalc, s.xIIcalc)
s.update_xcalc()
print(s.xIcalc, s.xIIcalc)
