import sys
import wave
import numpy
import scipy

import audio_wrapper
import encoder
import preprocessing
import decoder


# read data from wav file 
sample_rate,audio_data = audio_wrapper.scipy_read_data("ena_dio_tria.wav")
iterations = len(audio_data) // 160     # // for integer division

all_frames = []
audio_array = numpy.array([], dtype=numpy.int16)

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
    LARd,curr_frame_st_residual=encoder.RPE_frame_st_coder(s)

    #decoder

    S0=decoder.RPE_frame_st_decoder(LARd,curr_frame_st_residual)
    #print('iteration j = ', j, ', samples [', j * 160, ', ', (j+1) * 160, '] out of ', len(audio_data))

    all_frames.append(S0)
    s = numpy.ravel(S0)
    audio_array = numpy.concatenate((audio_array, s))
audio_array = audio_array.astype(numpy.int16)
output_filename = 'reconstructed_audio.wav'
scipy.io.wavfile.write(output_filename, sample_rate, audio_array)

print(audio_array,audio_data)