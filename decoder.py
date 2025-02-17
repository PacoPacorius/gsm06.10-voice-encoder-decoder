import sys
import wave
import numpy as np
import hw_utils as ut
from scipy.signal import lfilter
import scipy
from scipy.signal import dimpulse
import bitstring
from bitstring import BitStream, BitArray, pack

def RPE_frame_st_decoder(frame_bit_stream,prev_frame_st_residual: np.ndarray):
 frame_bit_stream.pos=0

 (LARc0,LARc1,LARc2,LARc3,LARc4,LARc5,LARc6,LARc7,
  Nj0,bcj0,m_maxj0,xmaxc0,Xmc0j0,Xmc1j0,Xmc2j0,Xmc3j0,Xmc4j0,Xmc5j0,Xmc6j0,Xmc7j0,Xmc8j0,Xmc9j0,Xmc10j0,Xmc11j0,Xmc12j0,
  Nj1,bcj1,m_maxj1,xmaxc1,Xmc0j1,Xmc1j1,Xmc2j1,Xmc3j1,Xmc4j1,Xmc5j1,Xmc6j1,Xmc7j1,Xmc8j1,Xmc9j1,Xmc10j1,Xmc11j1,Xmc12j1,
  Nj2,bcj2,m_maxj2,xmaxc2,Xmc0j2,Xmc1j2,Xmc2j2,Xmc3j2,Xmc4j2,Xmc5j2,Xmc6j2,Xmc7j2,Xmc8j2,Xmc9j2,Xmc10j2,Xmc11j2,Xmc12j2,
  Nj3,bcj3,m_maxj3,xmaxc3,Xmc0j3,Xmc1j3,Xmc2j3,Xmc3j3,Xmc4j3,Xmc5j3,Xmc6j3,Xmc7j3,Xmc8j3,Xmc9j3,Xmc10j3,Xmc11j3,Xmc12j3)=frame_bit_stream.unpack('int:6,int:6,int:5,int:5,int:4,int:4,int:3,int:3,'#LARc
                                                                                                                                          'uint:7,uint:2,uint:2,uint:6,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,'#first subframe
                                                                                                                                          'uint:7,uint:2,uint:2,uint:6,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,'#second subframe
                                                                                                                                          'uint:7,uint:2,uint:2,uint:6,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,'#third
                                                                                                                                          'uint:7,uint:2,uint:2,uint:6,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3,uint:3')#forth

 #after unpacking the bitstream , we seperate in arrays
 LARc=[LARc0,LARc1,LARc2,LARc3,LARc4,LARc5,LARc6,LARc7]
 N=[Nj0,Nj1,Nj2,Nj3]
 bc=[bcj0,bcj1,bcj2,bcj3]
 m_max=[m_maxj0,m_maxj1,m_maxj2,m_maxj3]
 Xmc0=[Xmc0j0,Xmc1j0,Xmc2j0,Xmc3j0,Xmc4j0,Xmc5j0,Xmc6j0,Xmc7j0,Xmc8j0,Xmc9j0,Xmc10j0,Xmc11j0,Xmc12j0]
 Xmc1=[Xmc0j1,Xmc1j1,Xmc2j1,Xmc3j1,Xmc4j1,Xmc5j1,Xmc6j1,Xmc7j1,Xmc8j1,Xmc9j1,Xmc10j1,Xmc11j1,Xmc12j1]
 Xmc2=[Xmc0j2,Xmc1j2,Xmc2j2,Xmc3j2,Xmc4j2,Xmc5j2,Xmc6j2,Xmc7j2,Xmc8j2,Xmc9j2,Xmc10j2,Xmc11j2,Xmc12j2]
 Xmc3=[Xmc0j3,Xmc1j3,Xmc2j3,Xmc3j3,Xmc4j3,Xmc5j3,Xmc6j3,Xmc7j3,Xmc8j3,Xmc9j3,Xmc10j3,Xmc11j3,Xmc12j3]
 xmaxc_allsfr=[xmaxc0,xmaxc1,xmaxc2,xmaxc3]
 Xmc_list=([Xmc0,Xmc1,Xmc2,Xmc3])
 e_full=np.zeros(160)

 #############################
 ######## 2ο Επίπεδο #########
 #############################

 QLB = np.array([0.1, 0.35, 0.65, 1])
 prev_d = np.zeros(120)
 d_reconstruct = np.zeros(160)
 d_predict = np.zeros(40)
 bd = [0] * 4

 for j in range(0,4):
 # N is already decoded
 
 # decode bc
  bd[j] = QLB[bc[j]]


  # prediction can be in subframes of previous frame
  prev_d = np.concatenate((prev_frame_st_residual[range((j+1) * 40, 160)], d_reconstruct[range(0, j*40)]))

  Xc = xnorm_dequant(Xmc_list[j])
  XMapr = Xc * (x_dequant(xmaxc_allsfr[j]))
  print('X approximation decoder', XMapr)
  e_aprox = np.zeros(40)
  for p in range(13):

      e_aprox[m_max[j] + 3*p] = XMapr[p]

  for u in range(40):
      e_full[j*40 + u] = e_aprox[u]

  for i in range(0,40):
   d_predict[i] = bd[j] * prev_d[120 + i - N[j]]
 
   # calculate reconstructed st residual
   d_reconstruct[j*40 + i] = e_full[j*40 + i] + d_predict[i]

 

 print()
 print()
 print("========================")
 print("== 2o epipedo decoder ==")
 print("========================")
 print()
 print()
 #print("decoder bd = ", bd)
 print("decoder d_reconstruct = ", d_reconstruct)
 #print("decoder curr_frame_st_resd = ", curr_frame_st_resd)

 curr_frame_st_resd = d_reconstruct

 #############################
 ######## 1ο Επίπεδο #########
 #############################

 print()
 print()
 print("========================")
 print("== 1o epipedo decoder ==")
 print("========================")
 print()
 print()


 r = np.empty(8) # pre defining the r array
 A = np.array([20.000,20.000,20.000,20.000,13.637,15.000,8.334,8.824])
 B = np.array([0.000,0.000,4.000,-5.000,0.184,-3.500,-0.666,-2.235])
 LAR = np.zeros(8)


 # code to convert larc into lar
 for j in range(len(LARc)):
  LAR[j] = (LARc[j] - B[j]) / A[j]


 # code to convert lar to r according to equation 3.5
 for i in range(len(LAR)):
  if abs(LAR[i]) < 0.675:
    r[i] = LAR[i]
  elif (0.675 <= abs(LAR[i])) and (abs(LAR[i]) < 1.225):
    r[i] = np.sign(LAR[i]) * ((0.500*abs(LAR[i])) + 0.337500)
  elif (1.225 <= abs(LAR[i])) and (abs(LAR[i]) <= 1.625):
    r[i] = np.sign(LAR[i]) * ((0.125*abs(LAR[i])) + 0.796875)

 # converting r to ak
 print('r',r)
 a_k, e_final =ut.reflection_coeff_to_polynomial_coeff(r)

 # constructing the decoding filter H using the ak values
 H = np.empty(9)
 c = 1e-32

 a_k = -a_k[1:]
 print('a_k', a_k)
 for z in range(1, len(H)):
  H[z] = 1 / (1 - sum(a_k[k]*((z+c)**(-k-1)) for k in range(len(a_k))))
 H[0] = 1


 # by applying the H filter to curr_frame_st_resd, we get s after preprocessing
 S = np.convolve(curr_frame_st_resd, H, mode='same')

 # reverting preprocessing
 beta  = 28180  * (2**(-15))
 alpha = 32735 * (2**(-15))

 # coding equations ,for reference
 #Sof(k) = S(k) - S(k-1) + alpha*Sof(k-1)
 #s(k) = Sof(k) - beta*Sof(k-1)


 #Sof[k] = S[k] + beta * Sof[k-1]            # this is the post processing
 #S0[k]  = Sof[k] + S0[k−1] − alpha*Sof[k−1] # reverse offset

 # turning the decoding equations to filters
 b1 = [1]  # filter coefficients
 a1 = [1, -beta]

 # applying lfilter to S
 S_poste = lfilter(b1, a1, S)  # reverting pre emphasis

 # creating filter to apply to S_poste
 b2 = [1, -alpha]
 a2 = [1, -1]

 S0 = lfilter(b2, a2, S_poste)  # reverting offset

 return S0, curr_frame_st_resd



 ######################
 ## HELPER FUNCTIONS ##
 ######################

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

def xnorm_dequant(x):
    values=([-28672,-20480,-12288,-4096,4096,12288,20480,28672])
    #since Xmc is an integer from 0 to  , de quantization can be done as so

    Xc = np.take(values,x)
    Xc = Xc / (2**15)
    return Xc


dequant_levels=[31,63,95,127,159,191,223,255,287,319,351,383,415,447,479,511,575,639,703,767,831,895,959,1023,1151,1279,1407,1535,1663,1791,1919,2047,2303,2559,2815,3071,3327,3583,3839,4095,4607,5119,5631,6143,6655,7167,7679,8191,9215,10239,11263,12287,13311,14335,15359,16383,18431,20479,22527,24575,26623,28671,30719,32767]

def x_dequant(x_maxc):
    return dequant_levels[x_maxc]
