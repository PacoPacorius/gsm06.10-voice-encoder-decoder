import numpy

def offset_compensation(s0: numpy.ndarray):
    s0 = s0.astype(numpy.float64)
    alpha = 32735*2**(-15)
    for k in range (1, len(s0)):
        s0[k] = s0[k] - s0[k-1] + alpha*s0[k-1]
    return s0

def pre_emphasis(s0: numpy.ndarray):
    beta = 28180*2**(-15)
    for k in range (1, len(s0)):
        s0[k] = s0[k] - beta*s0[k-1]
    return s0
