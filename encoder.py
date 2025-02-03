import numpy
import hw_utils
import scipy.signal

def RPE_frame_st_coder(s: numpy.ndarray, prev_frame_st_residual: numpy.ndarray):
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


    #print('r = ', r)
    #print('R = ', R)
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
            a_mod[i] =-a[i-1]

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
        elif (abs_kr >= 0.675) and (abs_kr < 0.950):
            LAR[i] = (numpy.sign(kr[i])) * ((2*abs_kr) - 0.675)
        elif (abs_kr >= 0.950) and (abs_kr <= 1):
            LAR[i] = (numpy.sign(kr[i])) * ((8*abs_kr) - 6.375)

    print('LAR = ', LAR, ' size of LAR = ', LAR.size)

    # quantize and encode LAR to LARc
    for i in range(0,8):
        LARc[i] = Nint((A(i)*LAR[i]) + B(i))

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
        elif (abs_LARd >= 0.675) and (abs_LARd < 1.225):
            krd[i] = numpy.sign(LARd[i]) * ( (0.500 * abs(LARd[i])) + 0.337500 )
        elif (abs_LARd >= 1.225) and (abs_LARd <= 1.625):
            krd[i] = numpy.sign(LARd[i]) * ( (0.125 * abs(LARd[i])) + 0.796875 )

    print('krd = ', krd, ' size of krd = ', krd.size)

    # get decoded akd from krd
    akd = numpy.zeros(9)
    akd, e_final = hw_utils.reflection_coeff_to_polynomial_coeff(krd)

    print('a = ', a, ' size of a = ', a.size, ' shape of a = ', a.shape)



    print('akd = ', akd, ' size of akd = ', len(akd))
    # apply FIR filter and calculate residual
    curr_frame_st_residual =numpy.convolve(s, akd, 'same')
    print('curr_frame_st_residual = ', curr_frame_st_residual, ' size of curr_frame_st_residual = ', curr_frame_st_residual.size)


    #return LARc, curr_frame_st_residual

    #############################
    ######## 2ο Επίπεδο #########
    #############################
    
    print()
    print()
    print("========================")
    print("== 2o epipedo encoder ==")
    print("========================")
    print()
    print()

    j = 0
    d_prev = prev_frame_st_residual 
    d_current = curr_frame_st_residual
    d_reconstruct = numpy.zeros(160)
    N = [0] * 4
    b = [0] * 4
    d_total = numpy.zeros([])
    bc = [0] * 4
    bd = [0] * 4

    ## Estimation ## 
    

    for j in range(0,4):
        ## construct prev_d
        #
        # for the 1st subframe, we should have 120 samples from the 
        # previous frame
        #
        # for the 2nd subframe, we should have 80 samples from the
        # previous frame followed by 40 reconstructed samples from 
        # the current frame
        #
        # for the 3rd subframe, we should have 40 samples from the
        # previous frame followed by 80 reconstructed samples from
        # the current frame
        #
        # for the 4th subframe, we should have 120 reconstructed
        # samples from the current frame
        prev_d = numpy.concatenate((prev_frame_st_residual[range((j+1) * 40, 160)], d_reconstruct[range(0, j*40)]))
        d = curr_frame_st_residual[range(j*40, (j+1)*40)]
        print("prev_d = ", prev_d, "size of prev_d = ", len(prev_d))
        # calculate N and b
        N[j], b[j] = RPE_subframe_slt_lte(d, prev_d)

        # N is already an int, it's already quantized

        # quantize b
        print("b = ", b[j])
        if b[j] <= DLB(0):
            bc[j] = 0
        elif b[j] > DLB(0) and b[j] <= DLB(1):
            bc[j] = 1
        elif b[j] > DLB(1) and b[j] <= DLB(2):
            bc[j] = 2
        elif b[j] > DLB(2):
            bc[j] = 3
        
        ## Prediction ##

        # N is just an int, no need to decode
        # decode b
        bd[j] = QLB(bc[j])

        e = numpy.zeros(160)
        d_predict = numpy.zeros(40)
        for i in range(0, 40):
            # calculate prediction
            d_predict[i] = bd[j] * prev_d[120 + i - N[j]]

            # calculate residual
            e[j*40 + i] = d[i] - d_predict[i]

            ## Synthesis ##
            # calculate reconstructed st residual
            d_reconstruct[j*40 + i] = e[j*40 + i] + d_predict[i]

    print("N = ", N)
    #print("bc = ", bc, "type of bc = ", type(bc[j]), ", number of bc bits = ", bc[j].bit_count())
    #print("bd = ", bd, "type of bd = ", type(bd[j]))
    print("d_predict = ", d_predict, ", size of d_predict = ", len(d_predict))
    print("e = ", e, ", size of e = ", len(e))
    print("d_reconstruct = ", d_reconstruct, ", size of d_reconstruct = ", len(d_reconstruct))


    return LARc, curr_frame_st_residual, N, bc, e





def RPE_subframe_slt_lte(d: numpy.ndarray, prev_d: numpy.ndarray):

    # find N=λ maximizer of auto-correlation Rj(λ)
    R = 0
    max_R = 0
    maximizer_lamda = 40

    for lamda in range(40,121):
        for i in range(0,40):
            R = R + d[i] * prev_d[120 + i - lamda]
        # keep max R and maximizing λ
        if R > max_R:
            max_R = R
            maximizer_lamda = lamda

    N = maximizer_lamda
    #print("N = ", N)
    #print("max R = ", max_R)

    # calculate b
    b_numerator = 0
    b_denominator = 0

    for i in range(0,40):
        b_numerator = b_numerator + d[i] * prev_d[120 + i - lamda]
        b_denominator = b_denominator + prev_d[120 + i - lamda] * prev_d[120 + i - lamda] 

    #print("b_numerator = ", b_numerator)
    #print("b_denominator = ", b_denominator)
    b = b_numerator / b_denominator

    return N, b

###################################
######## HELPER FUNCTIONS #########
###################################


# round to closest integer value
def Nint(z):
    return int(z + (numpy.sign(z)*0.5))

# define LAR quantization and coding coefficients
def A(i):
    if i == 1 or i == 2 or i == 3 or i == 0:
        return 20.000
    elif i == 4:
        return 13.637
    elif i == 5:
        return 15.000
    elif i == 6:
        return 8.334
    elif i == 7:
        return 8.824
    else:
        return None

def B(i):
    if i == 1 or i == 0:
        return 0.000
    elif i == 2: 
        return 4.000
    elif i == 3:
        return -5.000
    elif i == 4:
        return 0.184
    elif i == 5:
        return -3.500
    elif i == 6:
        return -0.666
    elif i == 7:
        return -2.235
    else:
        return None

def QLB(i):
    if i == 0:
        return 0.1
    elif i == 1:
        return 0.35
    elif i == 2:
        return 0.65
    elif i == 3:
        return 1
        
def DLB(i):
    if i == 0:
        return 0.2
    elif i == 1:
        return 0.5
    elif i == 2:
        return 0.8
