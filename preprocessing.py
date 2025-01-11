import numpy

def offset_compensation(s0: numpy.ndarray):
    sof = numpy.zeros(len(s0.astype(numpy.float64)))
    alpha = 32735*(2**(-15))
    print('alpha = ', alpha)
    print('length of s0', len(s0))
    for k in range (1, len(s0)):
        sof[k] = s0[k] - s0[k-1] + alpha*sof[k-1]
        #print('sof[', k, '] = ', sof[k])
    return sof

def pre_emphasis(sof: numpy.ndarray):
    s = numpy.zeros(len(sof.astype(numpy.float64)))
    beta = 28180*(2**(-15))
    for k in range (1, len(s)):
        s[k] = sof[k] - beta*sof[k-1]
    return sof
