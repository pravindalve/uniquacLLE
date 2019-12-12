import sys
sys.path.insert(0, 'C:/Users/FOSSEE/PycharmProjects/uniqacLLE')

from state import State

Texp = 298.15
xIexp = [0.8487, 0.0222, 0.1291]
xIIexp = [0.1866, 0.0328, 0.7806]

a = [[0, -87.1325, 17.6775],[501.439, 0, -232.51],[631.636, 338.947, 0]]

q = [1.432, 2.4, 4.396]

r = [1.4311, 3.1878, 5.1742]


s1 = State(Texp, xIexp,xIIexp, a, q, r)
s = 1
dist1 = 1
dist2 = 1
while dist1 > 1e-6 or dist2 > 1e-6:
	s1.update_xcalc()
	dist1res = s1.min_tan_dist(s1.xIcalc)
	dist2res = s1.min_tan_dist(s1.xIIcalc)

	# if dist1 > dist1res.fun:
	dist1 = dist1res.fun
	# if dist2 > dist2res.fun:
	dist2 = dist2res.fun

	s +=1
	if s > 60:
		break
print(s)
print(dist1, dist2)
print(s1.xIcalc, s1.xIIcalc)