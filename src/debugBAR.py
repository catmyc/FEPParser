import BAR
import random
'''
F = []
R = []
for i in range(10000):
    F.append(random.gauss(45,3))
    R.append(-1*random.gauss(45,3))
'''

F = [random.gauss(39, 3) for i in range(10000)]
R = [-random.gauss(41, 3) for i in range(10000)]

bar = BAR.BARestimator(F, R, 300)
print("DeltaF = %8.4f" % bar.BARSC(48, 1e-8))
