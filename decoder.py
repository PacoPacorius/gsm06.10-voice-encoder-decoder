import numpy
import hw_utils

def RPE_frame_st_coder(s0: numpy.ndarray):
    # calculate autocorrelations
    rs = numpy.zeros((9,), 'int32')
    for k in range (0,9):
        for i in range (k, 160):
            # estimate autocorrelation in accordance to section 3.1.4
            rs[k] = s0[i]*s0[i-k]
            print('\ns(i) = ', s0[i], 's(i-k) = ', s0[i-k])
    print('\n\n', rs)

    # Create w, R and r, matrices of the normal equations
    r = numpy.array([[rs[1]],
                     [rs[2]],
                     [rs[3]],
                     [rs[4]],
                     [rs[5]],
                     [rs[6]],
                     [rs[7]],
                     [rs[8]]])

    R = numpy.array([[rs[0], rs[1], rs[2], rs[3], rs[4], rs[5], rs[6], rs[7]],
                     [rs[1], rs[0], rs[1], rs[2], rs[3], rs[4], rs[5], rs[6]],
                     [rs[2], rs[1], rs[0], rs[1], rs[2], rs[3], rs[4], rs[5]],
                     [rs[3], rs[2], rs[1], rs[0], rs[1], rs[2], rs[3], rs[4]],
                     [rs[4], rs[3], rs[2], rs[1], rs[0], rs[1], rs[2], rs[3]],
                     [rs[5], rs[4], rs[3], rs[2], rs[1], rs[0], rs[1], rs[2]],
                     [rs[6], rs[5], rs[4], rs[3], rs[2], rs[1], rs[0], rs[1]],
                     [rs[7], rs[6], rs[5], rs[4], rs[3], rs[2], rs[1], rs[0]]])

    a = numpy.zeros(8)
    w = numpy.zeros((0,8))

    # solve normal equation to w
    w = numpy.linalg.solve(R, r)
    for i in range(0,8):
        a[i] = w[i]
    print('w after solving = ', w)
    print('a after solving = ', a)

    # calculate reflection coefficients

    kr = polynomial_coeff_to_reflection_coeff(a)
