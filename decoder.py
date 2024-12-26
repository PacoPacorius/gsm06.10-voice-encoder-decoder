import numpy

def RPE_frame_st_coder(s0: numpy.ndarray):
    # calculate autocorrelations
    rs = numpy.zeros((9,), 'int32')
    for k in range (0,9):
        for i in range (k, 160):
            # estimate autocorrelation in accordance to section 3.1.4
            rs[k] = s0[i]*s0[i-k]
            print('\ns(i) = ', s0[i], 's(i-k) = ', s0[i-k])
    print('\n\n', rs)
