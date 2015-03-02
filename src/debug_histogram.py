from Histogram import *
import random
data = [random.gauss(0, 5) for i in range(500000)]
H = Histogram(data)
H.stat(50)
H.printhist()