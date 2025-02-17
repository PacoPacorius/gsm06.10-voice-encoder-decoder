import numpy
import hw_utils
import scipy.signal
import bitstring
from bitstring import BitStream, BitArray, pack


def RPE_frame_st_coder(s: numpy.ndarray, prev_frame_st_residual: numpy.ndarray):
    # calculate autocorrelations
    rs = numpy.zeros((9,),numpy.float64)
    for k in range (0,9):
        for i in range (k, 160):
            # estimate autocorrelation in accordance to section 3.1.4
            rs[k] += s[i]*s[i-k]

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
        abs_kr = numpy.abs(kr[i])
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

    # convert LARc from float to int while respecting recommendation min-max LARc values
    LARc=LARc.astype(int)
    LARc_ranges = [
        (-32, 31),
        (-32, 31),
        (-16, 15),
        (-16, 15),
        (-8, 7),
        (-8, 7),
        (-4, 3),
        (-4, 3)]
    LARc_clipped = [numpy.clip(LARc[i], LARc_ranges[i][0], LARc_ranges[i][1]) for i in range(len(LARc))]
    LARc_clipped=numpy.asarray(LARc_clipped)
    #print('LARc = ', LARc, ' size of LARc = ', LARc.size)
    #print('LARc_clipped',LARc_clipped)
    bitstream = BitStream()
    bitstream.append(pack('int:6', LARc_clipped[0]))
    bitstream.append(pack('int:6', LARc_clipped[1]))
    bitstream.append(pack('int:5', LARc_clipped[2]))
    bitstream.append(pack('int:5', LARc_clipped[3]))
    bitstream.append(pack('int:4', LARc_clipped[4]))
    bitstream.append(pack('int:4', LARc_clipped[5]))
    bitstream.append(pack('int:3', LARc_clipped[6]))
    bitstream.append(pack('int:3', LARc_clipped[7]))

    #print(bitstream.bin, len(bitstream))


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



    # get decoded akd from krd
    akd = numpy.zeros(9)
    akd, e_final = hw_utils.reflection_coeff_to_polynomial_coeff(krd)

    print('a = ', a, ' size of a = ', a.size, ' shape of a = ', a.shape)

    akd[1:] = -akd[1:]

    # apply FIR filter and calculate residual
    curr_frame_st_residual =numpy.convolve(s, akd, 'same')
    #print('curr_frame_st_residual = ', curr_frame_st_residual, ' size of curr_frame_st_residual = ', curr_frame_st_residual.size)


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

    # declare vars
    j = 0
    prev_d = numpy.zeros(120)
    d_current = curr_frame_st_residual
    d_reconstruct = numpy.zeros(160)
    N  = [0] * 4
    b  = [0] * 4
    bc = [0] * 4
    bd = [0] * 4

    H   = numpy.array([-134,-374,0,2054,5741,8192,5741,2054,0,-374,-134])
    He  = H.astype(float)/(2**13)
    x   = numpy.empty(40)
    x0  = numpy.zeros(13)
    x1  = numpy.zeros(13)
    x2  = numpy.zeros(13)
    x3  = numpy.zeros(13)
    e_full = numpy.zeros(160)
    
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
        #print("prev_d = ", prev_d, "size of prev_d = ", len(prev_d))
        # calculate N and b
        N[j], b[j] = RPE_subframe_slt_lte(d, prev_d)

        # N is already an int, it's already quantized
        Nc = N[j]
        bitstream.append(pack('uint:7', Nc)) # appending Nc for each subframe with 7 bits , as the loop moves through j, each Nc is appended in its correct position

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
        
        bitstream.append(pack('uint:2', bc[j])) #appending b encoded to the bitstream

        ## Prediction ##

        # N is just an int, no need to decode
        # decode b
        bd[j] = QLB(bc[j])

        e = numpy.zeros(160)
        d_predict = numpy.zeros(40)
        for i in range(0,40):
            # calculate prediction
            d_predict[i] = bd[j] * prev_d[120 + i - N[j]]

            # calculate residual
            e[j*40 + i] = d[i] - d_predict[i]


        # the process below is done for each subframe, each subframe has 4 xm sequences each
        for k in range(40):

            # convolving according to the standard
            x[k] = sum(He[n] * e[j*40 + k] for n in range(len(He))) 

        # these xm are for one subframe only
        for d in range(13):
            x0[d] = x[0 + 3 * d]
            x1[d] = x[1 + 3 * d]
            x2[d] = x[2 + 3 * d]
            x3[d] = x[3 + 3 * d]

        # computing the sequence xm with the highest energy
        sequences = [x0,x1,x2,x3]
        E = [numpy.sum(x ** 2) for x in sequences]
        m_max = numpy.argmax(E)                 # m_max is the number of the sequence, m
        bitstream.append(pack('uint:2', m_max)) # appending M, the RPE grid position
        XM = sequences[m_max]                   # XM, the sequence with the highest energy

        x_max = numpy.max(numpy.abs(XM))        # max value of every element in sequence
        xmaxc = x_quant(x_max)

        #print('xmax',x_max)
        #print('xmaxc',xmaxc)


        #print('xmaxc',xmaxc_int)
        bitstream.append(pack('uint:6', xmaxc))

        # appending the sequence with the highest energy of each subframe, XM
        xnorm = XM / (x_dequant(xmaxc))  # using the dequantized xmaxc
        #print('xnorm',xnorm)
        Xmc = xnorm_quant(xnorm)           # quantizing xnormalized according to the standard using the function defined at the end of the code
        #print('Xmc =',Xmc)
        Xmc = Xmc.astype(int)
        bitstream.append(pack('uint:3', Xmc[0]))
        bitstream.append(pack('uint:3', Xmc[1]))
        bitstream.append(pack('uint:3', Xmc[2]))
        bitstream.append(pack('uint:3', Xmc[3]))
        bitstream.append(pack('uint:3', Xmc[4]))
        bitstream.append(pack('uint:3', Xmc[5]))
        bitstream.append(pack('uint:3', Xmc[6]))
        bitstream.append(pack('uint:3', Xmc[7]))
        bitstream.append(pack('uint:3', Xmc[8]))
        bitstream.append(pack('uint:3', Xmc[9]))
        bitstream.append(pack('uint:3', Xmc[10]))
        bitstream.append(pack('uint:3', Xmc[11]))
        bitstream.append(pack('uint:3', Xmc[12]))

        ##- Synthesis -##

        print('Xmc', Xmc)
        # XM approximation, mimicking the decoder
        Xc = xnorm_dequant(Xmc)
        XMapr = Xc * (x_dequant(xmaxc))

        print('X approximation', XMapr)
        # print('original XM',XM)

        e_aprox = numpy.zeros(40)

        # upsampling, keeping only the elements that are in XMapr, the rest are all zero

        for p in range(13):
            e_aprox[m_max + 3*p] = XMapr[p]
        # calculating d_reconstruct with the new e, e_aprox
        for u in range(40):
            d_reconstruct[j*40 + u] = e_aprox[u] + d_predict[u]
            e_full[j*40 + u]=e_aprox[u]


        #print('xnormalized',xnorm)
        #print('e_aprox = ', e_aprox,len(e_aprox))
        #print('x_max value',x_max)
        #print('quantized value  xmaxc',xmaxc)
        #print('xm',x0,x1,x2,x3)
    #print("x = ", x, ' size of x = ', x.size)
    #print("N = ", N)
    #print("bc = ", bc, "type of bc = ", type(bc[j]), ", number of bc bits = ", bc[j].bit_count())
    #print("bd = ", bd, "type of bd = ", type(bd[j]))
    #print("d_predict = ", d_predict, ", size of d_predict = ", len(d_predict))
    #print("e_full = ", e_full, ", size of e = ", len(e))
    #print("d_reconstruct = ", d_reconstruct, ", size of d_reconstruct = ", len(d_reconstruct))



    frame_bit_stream = bitstream
    #print(frame_bit_stream,len(frame_bit_stream))
    return  frame_bit_stream, curr_frame_st_residual





def RPE_subframe_slt_lte(d: numpy.ndarray, prev_d: numpy.ndarray):

    # find N=λ maximizer of auto-correlation Rj(λ)
    R = 0
    max_R = 0
    maximizer_lamda = 40

    for lamda in range(40, 121):
        for i in range(0, 40):
            R = R + (d[i] * prev_d[120 + i - lamda])
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

    for i in range(0, 40):
        # not taking any chances
        b_numerator   = b_numerator + (d[i] * prev_d[120 + i - N])
        b_denominator = b_denominator + (prev_d[120 + i - N] * prev_d[120 + i - N])

    #print("b_numerator = ", b_numerator)
    #print("b_denominator = ", b_denominator)
    b = b_numerator / b_denominator

    return N, b

###################################
######## HELPER FUNCTIONS #########
###################################


# round to closest integer value
def Nint(z):
    return int(z + numpy.sign(z)*0.5)

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

# quantization function for xmax according to the standard ,page 30
dequant_levels = [31,63,95,127,159,191,223,255,287,319,351,383,415,447,479,511,575,639,703,767,831,895,959,1023,1151,1279,1407,1535,1663,1791,1919,2047,2303,2559,2815,3071,3327,3583,3839,4095,4607,5119,5631,6143,6655,7167,7679,8191,9215,10239,11263,12287,13311,14335,15359,16383,18431,20479,22527,24575,26623,28671,30719,32767]
ranges = [0,31,63,95,127,159,191,223,255,287,319,351,383,415,447,479,511,575,639,703,767,831,895,959,1023,1151,1279,1407,1535,1663,1791,1919,2047,2303,2559,2815,3071,3327,3583,3839,4095,4607,5119,5631,6143,6655,7167,7679,8191,9215,10239,11263,12287,13311,14335,15359,16383,18431,20479,22527,24575,26623,28671,30719,32767]
levels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63]

def x_quant(x):
    # quantization according to the ranges in the standard, using the digitize function
    quant_index = numpy.digitize(x, ranges, right=True)

    return levels[quant_index-1]


def x_dequant(x_maxc):
    return dequant_levels[x_maxc]


def xnorm_quant(x):
    xscal = x * 2**15

    q_index = xscal // 8192 #step size of 8192
    Xmc = numpy.clip(q_index + 4, 0, 7) #quantized xnorm ,Xmc

    return Xmc


# function for de quantizing Xmc, the reverse of the above
def xnorm_dequant(x):
    values = ([-28672,-20480,-12288,-4096,4096,12288,20480,28672])
    # since Xmc is an integer from 0 to, de quantization can be done as so

    Xc = numpy.take(values, x)
    Xc = Xc / (2**15)
    return Xc

