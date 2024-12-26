import sys
import wave
import numpy

import audio_wrapper
import decoder

s0 = audio_wrapper.read_data("ena_dio_tria.wav", 160)
decoder.RPE_frame_st_coder(s0)
