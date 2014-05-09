import numpy as np
import matplotlib.pyplot as plt
import wave
from scikits.audiolab import *
import scipy
from scipy import signal
import math
from pylab import*


def read_audio(filename):
	spf = wave.open("majorE.wav",'r')
	signal = spf.readframes(-1)
	signal = np.fromstring(signal, 'Int16')
	p = spf.getnchannels()
	f = spf.getframerate()
	sound_info = np.zeros(len(signal),dtype=float)
	signal = signal.astype(np.float)
	sound_info = signal/max(signal)

	#sound_info = sound_info[1:len(sound_info):2]
	if p==2:
		sound_info = scipy.signal.decimate(sound_info,2)

	return p ,f , sound_info


def spectrogram(sound_info,f,nfft,hop):
	Pxx, freqs, bins, im = specgram(sound_info, Fs = f, NFFT = nfft, noverlap=nfft-hop, scale_by_freq=True,sides='default')
	return Pxx, freqs, bins, im


def hz_to_oct(freq):
	#A6 = 27.5*(2^6)
	#A0 = 27.5
	#A4 = 440
	fmin = 27.5
	b = 24
	return np.log2(freq/fmin)*b


def oct_to_hz(oct):
	fmin = 27.5
	b = 24.0
	return fmin*(2**(oct/b))


def generate_filterbank(NFFT,fs,b,z):
	#b is bins per octave
	#z is number of octaves
	#b = 24
	#z = 6
	#fs(downsampled) = 44100/4 = 11025

	octmax = b*z

	octpts = np.arange(-1,octmax+1)
	#print 'octpts',octpts
	#print len(octpts)

	ctrfrq = oct_to_hz(octpts)
	#print "ctrfrq",ctrfrq
	#print len(ctrfrq)

	ctrrep = np.floor((NFFT+2)*ctrfrq/((fs/2)))
	#print "ctrrep",ctrrep
	#print len(ctrrep)

	bank = np.zeros([len(octpts)-2,NFFT/2+1],dtype=float)

	for j in xrange(0,len(octpts)-2):
		y = np.hamming(ctrrep[j+2]-ctrrep[j])
		area = trapz(y, dx=5)
		if area==0:
			area=1
		y2 = (y/area)
		bank[j,ctrrep[j]:ctrrep[j+2]] = y2
		#plot(bank[j,:])

	#show()
	return bank


p, f, sound_info_old = read_audio('majorE.wav')
#print "frequency is",f
#print "channels are",p
sound_info_deci = scipy.signal.decimate(sound_info_old,4)
#165375 is 15 sec 330750 is 30 sec for downsampled signal
sound_info = sound_info_deci[0:165375]
#print "length of audio is ",len(sound_info)
f = f/4

#plot(sound_info)
#show()


"""
Compute spectrogram
"""

Pxx, freqs, bins, im = spectrogram(sound_info, f, 6000, 128)
#print "shape of Spectrogram is",shape(Pxx)
#plot(Pxx,freqs)
#show()

# log filterbank
bank = generate_filterbank(NFFT=6000,fs=f,b=24,z=6)
#print "shape of bank",shape(bank)
#im = imshow(bank,aspect='auto',interpolation='nearest',origin='lower')
#show()


"""
Generate log frequency spectrogram
"""

sal = np.dot(bank,Pxx)
b,col = shape(sal)

#Replace 0 by 1 in matrix before taking log
for i in range(0,len(sal)):
	for j in range(0,len(sal[0])):
		if sal[i][j] == 0:
			sal[i][j]+=1

sal = 10*np.log10(sal)
#Normalize
salm = np.zeros(shape(sal))
salm = (sal - sal.min())/(sal.max()-sal.min())
#salm = sal

"""
for i in range(col):
	salm[:,i] = sal[:,i]/sum(sal[:,i])
"""

#print "shape of CQ Spectrogram is",shape(salm)
subplot(2,1,1)
title('Log Frequency Spectrogram')
#xlim(0,(len(sound_info)/11025.00))
im = imshow(salm,aspect='auto',interpolation='nearest',origin='lower')
xticks(np.arange(0,col,172),[0,2,4,6,8,10,12,14])

#show()


"""
Generate Chroma
"""
row,col = shape(salm)
#print row

#print col
b = 24

chrm = np.zeros(b*col).reshape(b,col)

for i in range(col):
	c = salm[:,i]
	for j in range(b):
		chrm[j,i] = sum(c[j:row:b])

#print chrm

chrm_new = np.zeros(b*col).reshape(b,col)

#median
"""
for i in range(b):
	for j in range(col):
		chrm_new[i,j] = np.median(chrm[i,j-3:j+3])
"""

chrm_new = chrm

chrm_new_norm = np.zeros(b*col).reshape(b,col)

val = [0,0,0,0,0,0,0,0,0,0,0,0]

#normalise
for i in range(col):
	chrm_new_norm[:,i] = chrm_new[:,i]/sum(chrm_new[:,i])
	#chrm_new_norm[i,j] = np.median(chrm_new[i,j-10:j+10])

for i in range(12):
        val[i] = max(sum(chrm_new[2*i-1,56:120]), sum(chrm_new[2*i,52:120]))

#print val

max_fund = 0
chord_sec = 0
chord_trd = 0
chord_sev = 0
for i in range(12):
        if val[i] > max_fund:
                max_fund = val[i]
                max_fund_index = i

#print max_fund_index
val_new = val

if max_fund_index > 0:
        val_new[max_fund_index-1] = 0
if max_fund_index < 11:
        val_new[max_fund_index+1] = 0


for i in range(12):
        if val_new[i] > chord_sec:
                if i != max_fund_index:
                        chord_sec = val[i]
                        chord_sec_index = i

if chord_sec_index > 0:
        val_new[chord_sec_index-1] = 0
if chord_sec_index < 11:
        val_new[chord_sec_index+1] = 0
#print val_new

for i in range(12):
        if val_new[i] > chord_trd:
                if i != max_fund_index and i != chord_sec_index:
                        chord_trd = val[i]
                        chord_trd_index = i
"""
for i in range(12):
        if val_new[i] > chord_sev:
                if i != max_fund_index and i != chord_sec_index and chord_trd_index:
                        chord_sev = val[i]
                        chord_sev_index = i
"""

chord_extract = [0,0,0,0,0,0,0,0,0,0,0,0]
chord_extract[max_fund_index] = 1
chord_extract[chord_sec_index] = 1
chord_extract[chord_trd_index] = 1
# chord_extract[chord_sev_index] = 1

#print chord_extract

"""
major chords
"""

if chord_extract == [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0]:
        print "The best chord match is C major"
if chord_extract == [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1]:
        print "The best chord match is C# major"
if chord_extract == [1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0]:
        print "The best chord match is D major"
if chord_extract == [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0]:
        print "The best chord match is Eb major"
if chord_extract == [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1]:
        print "The best chord match is E major"
if chord_extract == [1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0]:
        print "The best match chord is F major"
if chord_extract == [0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0]:
        print "The best match chord is F# major"
if chord_extract == [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0]:
        print "The best match chord is G major"
if chord_extract == [0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1]:
        print "The best match chord is Ab major"
if chord_extract == [1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0]:
        print "The best match chord is A major"
if chord_extract == [0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0]:
        print "The best match chord is Bb major"
if chord_extract == [0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0]:
        print "The best match chord is B major"

"""
minor chords
"""

if chord_extract == [0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0]:
        print "The best match chord is C minor"
if chord_extract == [0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1]:
        print "The best match chord is C# minor"
if chord_extract == [1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0]:
        print "The best match chord is D minor"
if chord_extract == [0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0]:
        print "The best match chord is Eb minor"
if chord_extract == [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0]:
        print "The best match chord is E minor"
if chord_extract == [0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1]:
        print "The best match chord is F minor"
if chord_extract == [1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0]:
        print "The best match chord is F# minor"
if chord_extract == [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0]:
        print "The best match chord is G minor"
if chord_extract == [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1]:
        print "The best match chord is Ab minor"
if chord_extract == [1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0]:
        print "The best match chord is A minor"
if chord_extract == [0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]:
        print "The best match chord is Bb minor"
if chord_extract == [0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0]:
        print "The best match chord is B minor"


"""
major 7th chords
"""

if chord_extract == [0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0]:
        print "The best match chord is C major 7th"
if chord_extract == [0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1]:
        print "The best match chord is C# major 7th"
if chord_extract == [1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0]:
        print "The best match chord is D major 7th"
if chord_extract == [0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0]:
        print "The best match chord is Eb major 7th"
if chord_extract == [0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1]:
        print "The best match chord is E major 7th"
if chord_extract == [1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0]:
        print "The best match chord is F major 7th"
if chord_extract == [0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0]:
        print "The best match chord is F# major 7th"
if chord_extract == [0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0]:
        print "The best match chord is G major 7th"
if chord_extract == [0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1]:
        print "The best match chord is Ab major 7th"
if chord_extract == [1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1]:
        print "The best match chord is A major 7th"
if chord_extract == [1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0]:
        print "The best match chord is Bb major 7th"
if chord_extract == [0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0]:
        print "The best match chord is B major"



#print chrm_new_norm[9, 76]


subplot(2,1,2)
title('Chroma')

im2 = imshow(chrm_new_norm,aspect='auto',interpolation='nearest',origin='lower')
#xlim(0,(len(sound_info)/11025.00))
yticks([1,3,5,7,9,11,13,15,17,19,21,23], ['A','A#','B','C','C#','D','D#','E','F','F#','G','G#'])
xticks(np.arange(0,col,172),[0,2,4,6,8,10,12,14])

show()

print "finish"

