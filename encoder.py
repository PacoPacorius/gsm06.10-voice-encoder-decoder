import numpy
import hw_utils
import scipy.signal

def RPE_frame_st_coder(s: numpy.ndarray):
    # calculate autocorrelations
    rs = numpy.zeros((9,),numpy.float64)
    for k in range (0,9):
        for i in range (k, 160):
            # estimate autocorrelation in accordance to section 3.1.4
            rs[k] += s[i]*s[i-k]
            #print('\ns(i) = ', s[i], 's(i-k) = ', s[i-k])
    print('\n\nrs = ', rs, ' shape of rs = ', rs.shape)

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
    #print('w = ', w)
    # solve normal equation to w
    w = numpy.linalg.solve(R, r)


    print('w after solving = ', w, ' size of w = ', w.size, ' shape of w = ', w.shape)
    for i in range(0,8):
        a[i] = w[i]
    print('a after solving = ', a, ' size of a = ', a.size, ' shape of a = ', a.shape)

    # prepare polynomial_coeff_to_reflection_coeff's input array
    a_mod = numpy.zeros(9)
    for i in range(0,9):
        if i == 0:
            a_mod[i] = 1
        else:
            a_mod[i] = -a[i-1]

    print('Input array to polynomial_coeff_to_reflection_coeff: ', a_mod)
    # calculate reflection coefficients
    kr = hw_utils.polynomial_coeff_to_reflection_coeff(a_mod)
    print('kr after solving = ', kr, ' size of kr = ', kr.size)
    
    # convert reflection coefficients to Log-Area-Ratios
    LAR = numpy.zeros(8)
    LARc = numpy.zeros(8)
    for i in range(0,8):
        abs_kr = numpy.absolute(kr[i])
        if abs_kr < 0.675:
            LAR[i] = kr[i]
        elif abs_kr >= 0.675 and abs_kr < 0.950:
            LAR[i] = numpy.sign(kr[i]) * (2*abs_kr - 0.675)
        elif abs_kr >= 0.950 and abs_kr <= 1:
            LAR[i] = numpy.sign(kr[i]) * (8*abs_kr - 6.375)

    print('LAR = ', LAR, ' size of LAR = ', LAR.size)

    # quantize and encode LAR to LARc
    for i in range(0,8):
        LARc[i] = Nint(A(i)*LAR[i] + B(i))

    print('LARc = ', LARc, ' size of LARc = ', LARc.size)

    LARd = numpy.zeros(8)
    # decode LARc to LARd
    for i in range(0,8):
        LARd[i] = (LARc[i] - B(i)) / A(i)
    print('LARd = ', LARd, ' size of LARd = ', LARd.size)

    # not implementing interpolation

    # LARd to reflection coefficients krd
    krd = numpy.zeros(8)
    for i in range(0,8):
        abs_LARd = numpy.absolute(LARd[i])
        if abs_LARd < 0.675:
            krd[i] = LARd[i]
        elif abs_LARd >= 0.675 and abs_LARd < 1.225:
            krd[i] = numpy.sign(LARd[i]) * ( 0.5 * abs(LARd[i]) + 0.3375 )
        elif abs_LARd >= 1.225 and abs_LARd <= 1.625:
            krd[i] = numpy.sign(LARd[i]) * ( 0.125 * abs(LARd[i]) + 0.796875 )

    print('krd = ', krd, ' size of krd = ', krd.size)

    # get decoded akd from krd
    akd = numpy.zeros(9)
    akd, e_final = hw_utils.reflection_coeff_to_polynomial_coeff(krd)
    print('akd = ', akd, ' size of akd = ', len(akd))
    print('a = ', a, ' size of a = ', a.size, ' shape of a = ', a.shape)

    
    # apply FIR filter and calculate residual
    curr_frame_st_residual = scipy.signal.convolve(s, akd, 'same')
    print('curr_frame_st_residual = ', curr_frame_st_residual, ' size of curr_frame_st_residual = ', curr_frame_st_residual.size)


    return LARc, curr_frame_st_residual

###################################
######## HELPER FUNCTIONS #########
###################################

# round to closest integer value
def Nint(z):
    return int(z + numpy.sign(z)*0.5)

# define LAR quantization and coding coefficients
def A(i):
    if i == 1 or i == 2 or i == 3 or i == 0:
        return 20.
    elif i == 4:
        return 13.637
    elif i == 5:
        return 15.
    elif i == 6:
        return 8.334
    elif i == 7:
        return 8.824
    else:
        return None

def B(i):
    if i == 1 or i == 0:
        return 0.
    elif i == 2: 
        return 4.
    elif i == 3:
        return -5.
    elif i == 4:
        return 0.184
    elif i == 5:
        return -3.5
    elif i == 6:
        return -0.666
    elif i == 7:
        return -2.235
    else:
        return None
