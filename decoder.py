import sys
import wave
import numpy as np
import hw_utils as ut
from scipy.signal import lfilter
import scipy
from scipy.signal import dimpulse

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
 print('kr',kr)
 a_k, e_final =ut.reflection_coeff_to_polynomial_coeff(kr)

 #constructing the decoding filter H using the ak values
 H=np.empty(9)
 c=1e-32

 a_k=a_k[1:]
 for z in range(1,len(H)):
  H[z]=1/(1-sum(a_k[k]*((z+c)**(-k-1)) for k in range(len(a_k))))
 H[0]=1






 #by applying the H filter to curr_frame_st_resd ,aka "d" , we get s after preprocessing

 S =np.convolve(curr_frame_st_resd,H,mode='same')

 #Sof=np.empty(len(S.astype(np.float64)))
 #S0= np.empty(len(S.astype(np.float64)))
 #reverting preprocessing
 beta=28180*(2**(-15))
 alpha =32735*(2**(-15))

 #coding equations ,for reference

 #Sof(k)=S0(k)-S0 (k-1)+alpha*Sof(k-1) , alpha = 32735*2**(-15)
 #s(k)=Sof(k)-beta*Sof(k-1),beta=28180*2**(-15)


 #Sof[k] = S[k] + beta * Sof[k-1] ,this is the post processing
 #
 #S0[k]=Sof[k]+S0[k−1]−alpha*Sof[k−1],reverse offset

 #turning the decoding equations to filters
 #might change that


 b1 = [1] #filter coefficients
 a1 = [1, -beta]

 #applying lfilter to S

 Sof = lfilter(b1, a1, S)

 #creating filter to apply to Sof

 b2 = [1, -alpha]
 a2 = [1, -1]

 S0= lfilter(b2, a2, Sof)



 #it should be the same as before
 return S0
