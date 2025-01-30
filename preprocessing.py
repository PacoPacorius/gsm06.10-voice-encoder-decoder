import numpy
from scipy.signal import lfilter
beta=28180*(2**(-15))
alpha =32735*(2**(-15))
def offset_compensation(s0: numpy.ndarray):
    sof = numpy.empty(len(s0.astype(numpy.float64)))

    print('alpha = ', alpha)
    print('length of s0', len(s0))
    b1 = [1, -1]
    a1 = [1, -alpha]


    sof = lfilter(b1, a1, s0)
    #for k in range (1, len(s0)):
        #sof[k] = s0[k] - s0[k-1] + alpha*sof[k-1]
        #print('sof[', k, '] = ', sof[k])
    return sof

def pre_emphasis(sof: numpy.ndarray):
    s = numpy.empty(len(sof.astype(numpy.float64)))

    b2 = [1]
    a2 = [1, -beta]


    s = lfilter(b2, a2, sof)
    #for k in range (1, len(s)):
        #s[k] = sof[k] - beta*sof[k-1]
    return s
