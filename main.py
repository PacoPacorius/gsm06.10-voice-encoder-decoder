import sys
import wave
import numpy

import audio_wrapper
import decoder
import preprocessing

# read data from wav file 
s0 = audio_wrapper.read_data("ena_dio_tria.wav", 160)
#s0 = audio_wrapper.read_data("/home/pacopacorius/test2.wav", 160)
print('main before pre-processing s0 = ', s0, ' s0 length: ', len(s0))

# offset compensation and pre-emphasis
s0 = preprocessing.offset_compensation(s0)
s0 = preprocessing.pre_emphasis(s0)

print('main after pre-processing s0 = ', s0.astype(numpy.int64), ' s0 length: ', len(s0))
# short term analysis
decoder.RPE_frame_st_coder(s0)
