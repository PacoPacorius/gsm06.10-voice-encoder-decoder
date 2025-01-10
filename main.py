import sys
import wave
import numpy

import audio_wrapper
import encoder
import preprocessing


# read data from wav file 
audio_data = audio_wrapper.scipy_read_data("ena_dio_tria.wav")
iterations = len(audio_data) // 160     # // for integer division
for j in range(0,iterations):
    # resize audio_data array
    s0 = numpy.zeros(160)
    offset = j * 160
    for i in range (offset, offset + 160):
        s0[i - offset] = audio_data[i]

    # offset compensation and pre-emphasis
    sof = preprocessing.offset_compensation(s0)
    s = preprocessing.pre_emphasis(sof)
    #print('s0 = ', s0, ' s0 length: ', len(s0))
    #print('sof = ', sof, ' sof length: ', len(sof))
    #print('s = ', s, ' s length: ', len(s))

    #print('main after pre-processing s0 = ', s0, ' s0 length: ', len(s0))
    # short term analysis
    encoder.RPE_frame_st_coder(s)
    print('iteration j = ', j, ', samples [', j * 160, ', ', (j+1) * 160, '] out of ', len(audio_data))
