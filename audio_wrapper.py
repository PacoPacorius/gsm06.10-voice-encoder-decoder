import wave

"""
Is it read or write?

First make it a Wave_read object
"""
#audio_data = wave.open("ena_dio_tria.wav", 'rb')
audio_data = wave.open("/home/pacopacorius/test2.wav", 'rb')
print(audio_data.getframerate(), '\n')
bytes_per_sample = audio_data.getsampwidth()
bytes_one_frame = audio_data.readframes(160) # this succesfully reads 160 frames from the audio_data
print(audio_data.getsampwidth(), '\n')
print(len(bytes_one_frame)) # the argument to hex() requires python 3.8 
audio_data.close();
