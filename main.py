import sys
import wave
import numpy
import scipy

import audio_wrapper
import encoder
import preprocessing
import decoder

# declare global variables
all_frames = []
audio_array = numpy.array([], dtype=numpy.float64)
prev_frame_st_residual = numpy.zeros(160)

# read data from wav file 
sample_rate,audio_data_o = audio_wrapper.scipy_read_data("ena_dio_tria.wav")
iterations = len(audio_data_o) // 160     # // for integer division


# pre-processing
audio_data_of = preprocessing.offset_compensation(audio_data_o)
audio_data = preprocessing.pre_emphasis(audio_data_of)

#iterations = 1
for j in range(0,iterations):
    # initialize s0
    s = numpy.zeros(160)
    offset = j * 160
    for i in range (offset, offset + 160):
        s[i - offset] = audio_data[i]

    # encoder
    LARc,curr_frame_st_residual=encoder.RPE_frame_st_coder(s, prev_frame_st_residual)

    # decoder
    S0=decoder.RPE_frame_st_decoder(LARc,curr_frame_st_residual)
    #print('iteration j = ', j, ', samples [', j * 160, ', ', (j+1) * 160, '] out of ', len(audio_data))

    all_frames.append(S0)
    prev_frame_st_residual = curr_frame_st_residual

audio_array = numpy.asarray(all_frames)
audio_array = numpy.ravel(audio_array)
#print('type of audio_array element = ', type(audio_array[0]))
audio_array = audio_array.astype(numpy.int16)       # our wav files will always have 16-bit samples

output_filename = 'reconstructed_audio.wav'
scipy.io.wavfile.write(output_filename, sample_rate, audio_array)

print('audio array = ', audio_array)
print('audio data/starter = ', audio_data)
