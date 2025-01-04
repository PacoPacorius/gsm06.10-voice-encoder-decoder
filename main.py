import sys
import wave
import numpy

import audio_wrapper
import encoder
import preprocessing

# read data from wav file 
audio_data = audio_wrapper.scipy_read_data("ena_dio_tria.wav")
# resize audio_data array
s0 = numpy.zeros(160)
for i in range (0,160):
    s0[i] = audio_data[i]
#print('scipy s0 = ', s0, ' s0 length: ', len(s0))

# offset compensation and pre-emphasis
s0 = preprocessing.offset_compensation(s0)
s0 = preprocessing.pre_emphasis(s0)

#print('main after pre-processing s0 = ', s0, ' s0 length: ', len(s0))
# short term analysis
encoder.RPE_frame_st_coder(s0)
