from uniquaclle import UniquacLLE
from state import State

Texp = [298.15, 298.15, 298.15, 298.15, 298.15, 298.15, 298.15, 298.15]
xIexp = [[0.8487, 0.0222, 0.1291], [0.8319, 0.0295, 0.1386], [0.815, 0.0368, 0.1482], [0.7986, 0.0439, 0.1575],
         [0.775, 0.0541, 0.1709], [0.7665, 0.0578, 0.1757], [0.7506, 0.0646, 0.1848], [0.7285, 0.074, 0.1975]]
xIIexp = [[0.1866, 0.0328, 0.7806], [0.2162, 0.058, 0.7258], [0.2394, 0.0748, 0.6858], [0.2568, 0.0863, 0.6569],
          [0.2817, 0.1004, 0.6179], [0.2927, 0.1059, 0.6014], [0.3106, 0.1137, 0.5757], [0.3382, 0.1231, 0.5387]]

a = [[0, -87.1325, 17.6775],[501.439, 0, -232.51],[631.636, 338.947, 0]]

q = [1.432, 2.4, 4.396]

r = [1.4311, 3.1878, 5.1742]

dataset = UniquacLLE(Texp, xIexp, xIIexp, a, q, r)

# value = dataset.optimize_OF1()
# print(dataset.a_nm)

# bips = dataset.evaluate_BIP()
# print(bips)

s1 = State(Texp[0], xIexp[0],xIIexp[0], a, q, r)
s = 1
m = 1
dist1 = 1
dist2 = 1
coef1 = 1
coef2 = 2
while dist1 > 1e-6 or dist2 > 1e-6 or abs(coef1 - coef2) > 50:
	s1.update_xcalc()
	dist1res = s1.min_tan_dist(s1.xIcalc)
	dist2res = s1.min_tan_dist(s1.xIIcalc)

	if dist1 > dist1res.fun:
		dist1 = dist1res.fun
	if dist2 > dist2res.fun:
		dist2 = dist2res.fun

	s +=1
	if s > 60:
		break


print(s1.xIexp, s1.xIIexp)
print(s1.xIcalc, s1.xIIcalc)