import numpy
import hw_utils

def RPE_frame_st_coder(s0: numpy.ndarray):
    # calculate autocorrelations
    rs = numpy.zeros((9,),numpy.float64)
    for k in range (0,9):
        for i in range (k, 160):
            # estimate autocorrelation in accordance to section 3.1.4
            rs[k] += s0[i]*s0[i-k]
            #print('\ns(i) = ', s0[i], 's(i-k) = ', s0[i-k])
    print('\n\nrs = ', rs)

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

    print('r = ', r)
    print('R = ', R)
    print('w = ', w)
    # solve normal equation to w
    w = numpy.linalg.solve(R, r)

    #try:
        #w = numpy.linalg.solve(R, r)
    #except numpy.linalg.LinAlgError as err:
        #if 'Singular matrix' in str(err):
            #print("""\nR is a singular matrix, approximating using the method 
                     #of least squares""")
            #w = numpy.linalg.lstsq(R, r)[0]
        #else:
            #raise

    print('w after solving = ', w)
    for i in range(0,8):
        a[i] = w[i]
    print('a after solving = ', a, ' size of a = ', a.size)

    # calculate reflection coefficients
    kr = hw_utils.polynomial_coeff_to_reflection_coeff(a)
    print('kr after solving = ', kr, ' size of kr = ', kr.size)

    # convert reflection coefficients to Log-Area-Ratios
    LARc = numpy.zeros(8)
    for i in range(0,8):
        abs_kr = numpy.absolute(kr[i-1])
        if abs_kr < 0.675:
            LARc[i] = kr[i-1]
        elif abs_kr >= 0.675 and abs_kr < 0.950:
            LARc[i] = numpy.sign(kr[i-1]) * (2*abs_kr - 0.675)
        elif abs_kr >= 0.950 and abs_kr <= 1:
            LARc[i] = numpy.sign(kr[i-1]) (8*abs_kr - 6.375)

    print('LARc = ', LARc)

