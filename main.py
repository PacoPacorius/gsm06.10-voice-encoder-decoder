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
audio_array = numpy.array([], dtype=numpy.int16)
prev_frame_st_resd = numpy.zeros(160)
prev_frame_lt_resd = numpy.zeros(160)

# read data from wav file 
input_filename = 'ena_dio_tria.wav'
sample_rate,audio_data_o = scipy.io.wavfile.read(input_filename)
# calculate how many times we're going to loop, equal to the number of frames 
iterations = len(audio_data_o) // 160     # // for integer division

# preprocessing
audio_data_of = preprocessing.offset_compensation(audio_data_o)
audio_data = preprocessing.pre_emphasis(audio_data_of)

#iterations = 1     # uncomment for easier debugging
for j in range(0,iterations):
    # initialize empty s
    s = numpy.zeros(160)
    # get new frame from audio data, store in it s
    offset = j * 160
    for i in range (offset, offset + 160):
        s[i - offset] = audio_data[i]

    #print('main after pre-processing s0 = ', s0, ' s0 length: ', len(s0))
    # short term analysis encoder
    LARc,curr_frame_st_residual = encoder.RPE_frame_slt_coder(s, prev_frame_st_resd, 
                                                              prev_frame_lt_resd)

    # st decoder
    S0 = decoder.RPE_frame_st_decoder(LARc,curr_frame_st_residual)
    # print debugging information
    print()
    print('main iteration j = ', j, ', samples [', j * 160, ', ', (j+1) * 160, '] out of ', len(audio_data))
    print()

    # χρειάζονται και τα all_frames και το audio_array ή απλά τα αποθηκεύουμε να τα έχουμε;
    all_frames.append(S0)
    s = numpy.ravel(S0)
    audio_array = numpy.concatenate((audio_array, s))
    prev_frame_st_resd = curr_frame_st_residual

#print('type of audio_array element = ', type(audio_array[0]))
audio_array = audio_array.astype(numpy.int16)       # our wav files will always have 16-bit samples! Otherwise encoding doesn't work
output_filename = 'reconstructed_audio.wav'
scipy.io.wavfile.write(output_filename, sample_rate, audio_array)

print('audio array = ', audio_array)
print('audio data = ', audio_data)
print('samples ignored near end-of-file = ', len(audio_data) - (j+1) * 160)     # this should always output less than 160
