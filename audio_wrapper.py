import sys
import wave
import numpy


def read_data(file_wav_audio, no_of_samples):
    # open wav file in read format
    audio_data = wave.open(file_wav_audio, 'rb')

    # read 160 frames (samples) from the audio file
    samples = audio_data.readframes(no_of_samples)

    # does our sound file have 8-bits or 16-bits per sample?
    sample_width = audio_data.getsampwidth()
    audio_data.close()

    # return one-dimensional numpy array with element size dependent on 
    # sample width (in bytes)
    if sample_width == 2:
        return numpy.frombuffer(samples, numpy.int16)
    elif sample_width == 1:
        return numpy.frombuffer(samples, numpy.int8)
    else:
        sys.exit("""
                Non-standard bits-per-sample WAV file. 
                Only 8 or 16-bits are allowed, quitting...""")



# get s0, then make it one-dimensional
#s0 = numpy.array([audio_wrapper_read_data("/home/pacopacorius/test2.wav", 160)])
#s0 = numpy.ravel(s0)
#print("Dimensions of numpy array s0 = ", numpy.ndim(s0), 
      #"\nSize of numpy array s0 = ", numpy.size(s0),
      #"\nShape of numpy array s0 = ", numpy.shape(s0),
      #"\nLength of numpy array s0 = ", len(s0),
      #"\nNumpy array s0 = ", s0,
      #"\ns0[5] = ", s0[5])
