import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from scipy.optimize import curve_fit


fileList = ['8_8.1_trimmed.csv', '8.1_8.2_trimmed.csv',
            '8.2_8.3_trimmed.csv', '8.3_8.4_trimmed.csv',
            '8.4_8.5_trimmed.csv', '8.5_8.6_trimmed.csv',
            '8.6_8.7_trimmed.csv', '8.7_8.8_trimmed.csv',
            '8.8_8.9_trimmed.csv', '8.9_9.0_trimmed.csv',
            '9.0_9.1_trimmed.csv', '9.1_9.2_trimmed.csv',
            '9.2_9.3_trimmed.csv', '9.3_9.4_trimmed.csv',
            '9.4_9.5_trimmed.csv', '9.5_9.6_trimmed.csv',
            '9.6_9.7_trimmed.csv', '9.7_9.8_trimmed.csv',
            '9.8_9.9_trimmed.csv', '9.9_10.0_trimmed.csv',
            '10.0_10.1_trimmed.csv', '10.1_10.2_trimmed.csv',
            '10.2_10.3_trimmed.csv', '10.3_10.4_trimmed.csv',
            '10.4_10.5_trimmed.csv', '10.5_10.6_trimmed.csv',
            '10.6_10.7_trimmed.csv', '10.7_10.8_trimmed.csv',
            '10.8_10.9_trimmed.csv', '10.9_11.0_trimmed.csv',
            '11.0_11.1_trimmed.csv', '11.1_11.2_trimmed.csv',
            '11.2_11.3_trimmed.csv', '11.3_11.4_trimmed.csv',
            '11.4_11.5_trimmed.csv', '11.5_11.6_trimmed.csv',
            '11.6_11.7_trimmed.csv', '11.7_11.8_trimmed.csv',
            '11.8_11.9_trimmed.csv', '11.9_12.0_trimmed.csv']

### set up histogram
numBins = 250
binMin = -2
binMax = 2
binSize = (binMax - binMin)/(numBins*1.0)
tickMarks = np.linspace(binMin, binMax, numBins)
hist = np.zeros(numBins)

### define a function to add data to an exisiting histogram
def addData(histogram, dataPoint):
   bin_num = int(np.floor((dataPoint-binMin)/binSize))
   if ((bin_num >= 0) & (bin_num < numBins)):
       histogram[bin_num] +=1
   return histogram

#set up pass filter
def PassFilter(frequencyArr, kspaceData, cutoffFreq):
    for i in range(0, len(frequencyArr)):
        if abs(frequencyArr[i]) >= cutoffFreq:
            kspaceData.imag[i] = 0
##        if abs(frequencyArr[i]) >= 2: #for even terms
##            kspaceData.real[i] = 0
    return kspaceData


### calculate error bar and <z>
S_avez = np.zeros(len(fileList))
avez = np.zeros(len(fileList))
phi_tick = np.zeros(len(fileList))

for i in range (0, len(fileList)):
    data = np.loadtxt(fileList[i], delimiter = ',', skiprows = 1)
    z = data[:] #does not take solar height into account
    zarry = list(filter((-999.7).__ne__, z)) ### remove all the -999.7 from the list
    
    N = len(zarry)
    hist = np.zeros(numBins)
    for k in zarry:
        addData(hist, k)
    peak = max(hist)
    nhist = hist/peak

    # compute average z
    kspaceData = np.fft.rfft(nhist)
    freq = np.fft.rfftfreq(tickMarks.shape[-1], d = binSize)
    xspaceData = np.fft.irfft(kspaceData)
    lowpassk = PassFilter(freq, kspaceData, 0.5) # cutoffFreq is for odd terms
    lowpassz = np.fft.irfft(lowpassk) # xspaceData after subtracting waves
    zf = tickMarks*lowpassz
    avez[i] = np.sum(zf)/np.sum(lowpassz)


##    plt.plot(tickMarks, xspaceData, 'b-', markersize = 0.1)
##    plt.plot(tickMarks, lowpassz, 'k.', markersize = 1.4)
##    plt.plot(tickMarks, nhist, 'r--', markersize = 0.1)
##    plt.xlabel('z (kpc)')
##    plt.ylabel('Normalized counts')
##    plt.legend(['Regular fft', 'Lowpass fft', 'Raw data'])
##    #plt.savefig(str(i+1) + '.pdf')
##    plt.title('odd cutoff 0.5__even cutoff 0')
##    plt.show()

##      plt.plot(freq, kspaceData.real)
##      plt.show()

    # compute error bar
    Sn = np.sqrt(hist) #sigma_norm
    multiplier = tickMarks*Sn
    Multiplier = sum(multiplier**2)
    S_avez[i] = np.sqrt((Multiplier/N**2)+((avez[i])**2/N)) #sigma_average z
    
   

    # phi tickMarks
    R_tick = [8.05, 8.15, 8.25, 8.35, 8.45, 8.55, 8.65, 8.75, 8.85, 8.95,
              9.05, 9.15, 9.25, 9.35, 9.45, 9.55, 9.65, 9.75, 9.85, 9.95,
              10.05,10.15,10.25,10.35,10.45,10.55,10.65,10.75,10.85,10.95,
              11.05,11.15,11.25,11.35,11.45,11.55,11.65,11.75,11.85,11.95]


##print (S_avez)
##print (avez)


plt.errorbar(R_tick, avez, S_avez, marker = '.', ls = 'none')
plt.title('Error bars for R[8.0, 12.0] l[225, 245] with G & b cuts')
plt.xlabel('R (kpc)')
plt.ylabel('<z> (kpc)')
##plt.savefig('error bar Fig_R[8.0_12.0]_L[225_245] G smaller than 18.pdf')
plt.show()

    
    

