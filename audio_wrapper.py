import sys
import wave
import numpy
import scipy




def scipy_read_data(file_wav_audio):
    sample_rate, audio_data = scipy.io.wavfile.read(file_wav_audio)
    #print("***\nScipy method\n***\nSample rate = ", sample_rate, "\naudio_data = ", audio_data, ", length = ", len(audio_data))
    return audio_data



#def read_data(file_wav_audio, no_of_samples, offset):
    ## open wav file in read format
    #audio_data = wave.open(file_wav_audio, 'rb')
#
    ## read 160 frames (samples) from the audio file
    ##print('total audio frames = ', audio_data.getnframes())
    #audio_data.setpos(offset)  # test different positions of the same file
    #samples = audio_data.readframes(no_of_samples)
#
    ## does our sound file have 8-bits or 16-bits per sample?
    #sample_width = audio_data.getsampwidth()
    #audio_data.close()
#
    ## return one-dimensional numpy array with element size dependent on 
    ## sample width (in bytes)
    #if sample_width == 2:
        #return numpy.frombuffer(samples, 'i2')
    #elif sample_width == 1:
        #return numpy.frombuffer(samples, 'i1')
    #else:
        #sys.exit("""
                #Non-standard bits-per-sample WAV file. 
                #Only 8 or 16-bits are allowed, quitting...""")
