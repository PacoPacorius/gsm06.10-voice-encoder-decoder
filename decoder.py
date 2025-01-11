import sys
import wave
import numpy as np
import hw_utils as ut
from scipy.signal import lfilter
import scipy

def RPE_frame_st_decoder(LARc: np.ndarray,curr_frame_st_resd: np.ndarray
 )-> np.ndarray:
 r=np.empty(8) # pre defining the r array
 A=np.array([20.000,20.000,20.000,20.000,13.637,15.000,8.334,8.824])
 B=np.array([0.000,0.000,4.000,-5.000,0.184,-3.500,-0.666,-2.235])
 LAR=np.zeros(8)


 #code to convert larc into lar
 for j in range(len(LARc)):
  LAR[j]=(LARc[j]-B[j])/A[j]


 #code to convert lar to r according to equation 3.5
 for i in range(len(LAR)):
  if abs(LAR[i])<0.675:
    r[i]=LAR[i]
  elif (0.675<=abs(LAR[i])) and (abs(LAR[i])<1.225):
    r[i]=np.sign(LAR[i])*(0.500*abs(LAR[i])+0.337500)
  elif (1.225<=abs(LAR[i])) and (abs(LAR[i])<=1.625):
    r[i]=np.sign(LAR[i])*(0.125*abs(LAR[i])+0.796875)

 #code for converting r to ak

 kr=r
 a,e_final=ut.reflection_coeff_to_polynomial_coeff(kr)

 #constructing the decoding filter H using the ak values
 H=np.zeros(9)
 epsilon = 1e-12

 for z in range(len(H)):
  H[z]=1/(1-sum(a[k]*(z + epsilon)**(-k) for k in range(len(a))))


 #by applying the H filter to curr_frame_st_resd ,aka "d" , we get s after preprocessing

 S = scipy.signal.convolve(curr_frame_st_resd,H,mode='same')


 #reverting preprocessing
 beta=28180*2**(-15)
 alpha = 32735*2**(-15)

 #coding equations ,for reference

 #Sof(k)=S0(k)-S0 (k-1)+alpha*Sof(k-1) , alpha = 32735*2**(-15)
 #s(k)=Sof(k)-beta*Sof(k-1),beta=28180*2**(-15)


 #Sof[k] = S[k] + beta * Sof[k-1] ,this is the post processing
 #
 #S0[k]=Sof[k]+S0[k−1]−alpha*Sof[k−1],reverse offset

 #turning the decoding equations to filters
 #might change that

 b = [1] #filter coefficients
 a = [1, -beta]

 #applying lfilter to S

 Sof = lfilter(b, a, S)

 #creating filter to apply to Sof

 b = [1, -alpha]
 a = [1, -1]

 S0= lfilter(b, a, Sof)



 #it should be the same as before
 return S0