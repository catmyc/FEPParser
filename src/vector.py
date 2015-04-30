__Author__ = '''
            Yuncheng Mao
             '''
__Email__ = '''
            catmyc@gmail.com
            maoyuncheng@mail.nankai.edu.cn
            '''
import math

def vecadd(vec1, vec2):
    if len(vec1) == len(vec2):
        return [a + b for (a, b) in zip(vec1, vec2)]
    else:
        print("ERROR: Two vectors must be of equal length!!!")
        return None

def vecaddsingle(vec, num):
    '''
    add a same value to each element of the vector
    '''
    return [a + num for a in vec]

def vecinv(vec):
    return [-f for f in vec]

def vecsub(vec1, vec2): # return vec1 - vec2
    return [a - b for (a, b) in zip(vec1, vec2)]

def vecsubsingle(vec, num):
    '''
    substract a same value from each element of the vector
    '''
    return [a - num for a in vec]

def vecdot(vec1, vec2):
    if len(vec1) == len(vec2):
        return math.fsum([a * b for (a, b) in zip(vec1, vec2)])
    else:
        print("ERROR: Two vectors must be of equal length!!!")
        return None

def vecmul(vec1, vec2):
    if len(vec1) == len(vec2):
        return [a * b for (a, b) in zip(vec1, vec2)]
    else:
        print("ERROR: Two vectors must be of equal length!!!")
        return None

def vecmulsingle(vec, num):
    return [a * num for a in vec]

def veclen2(vec):
    return vecdot(vec, vec)

def veclen(vec):
    return math.sqrt(vecdot(vec, vec))

def vecnorm(vec):
    return vecmulsingle(vec, 1.0 / veclen(vec))

def vecdist(vec1, vec2):
    return veclen(vecsub(vec1, vec2))

def vecsum(vec):
    return math.fsum(vec)

def vecmean(vec):
    return math.fsum(vec) / len(vec)

def vec_centralized(vec):
    return vecsubsingle(vec, vecmean(vec))

def vecvar2(vec):
    '''
    mean squar variance
    '''
    v = vec_centralized(vec)
    return vecdot(v, v) / (len(v) -1)

def vecvar(vec):
    return math.sqrt(vecvar2(vec))

def vec_standardized(vec):
    return vecmulsingle(vec_centralized(vec), 1.0 / vecvar(vec)) 

