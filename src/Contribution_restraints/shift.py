#!/usr/bin/env python3
def shiftA(filename):
    l = []
    A = []
    dA = []
    outbuff = ''
    for line in open(filename):
        if line.startswith('#'):
            outbuff += line
            continue
        buff = line.split()
        l.append(float(buff[0]))
        A.append(float(buff[1]))
        dA.append(float(buff[2]))

    zeroind = l.index(0)
    print("Lambda 0 has index %d!" % zeroind)
    baseA = A[zeroind]
    print("A = %.4f, where lambda = 0." % baseA)
    newA = []
    for i in A:
        newA.append(i - baseA)

    l.reverse()
    newA.reverse()
    dA.reverse()
    dA = [-f for f in dA] 
    dA.insert(0, dA.pop()) # move the zero at the end to the beginning

    for i in range(len(l)):
        outbuff += "%-8s  %9.4f  %9.4f\n" % (l[i], newA[i], dA[i])
    open(filename, 'w').write(outbuff)
    print("Output written.")

import sys
filename = sys.argv[1]
shiftA(filename)

