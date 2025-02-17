import sys
import wave
import numpy as np
import hw_utils as ut
from scipy.signal import lfilter
import scipy
from scipy.signal import dimpulse

def RPE_frame_st_decoder(LARc: np.ndarray, 
                         N: list[int], bc: list[int], curr_frame_ex_full: np.ndarray,
                         prev_frame_st_residual: np.ndarray):

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
  
  for i in range(0,40):
   d_predict[i] = bd[j] * prev_d[120 + i - N[j]]
 
   # calculate reconstructed st residual
   d_reconstruct[j*40 + i] = curr_frame_ex_full[j*40 + i] + d_predict[i]

 

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
