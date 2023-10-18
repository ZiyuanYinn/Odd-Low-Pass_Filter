import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from scipy.optimize import curve_fit


fileList = [] ### put file lists here. Ex. '8.0_8.1_trimmed.csv'

avez = np.zeros(len(fileList))

for i in range (0, len(fileList)):
    data = np.loadtxt(fileList[i], delimiter = ',', skiprows = 1)
    z = data[:] #does not take solar height into account
    zarry = list(z)
    avez[i] = sum(zarry)/len(zarry)

 ### phi tickMarks
    R_tick = [8.05, 8.15, 8.25, 8.35, 8.45, 8.55, 8.65, 8.75, 8.85, 8.95,
              9.05, 9.15, 9.25, 9.35, 9.45, 9.55, 9.65, 9.75, 9.85, 9.95,
              10.05,10.15,10.25,10.35,10.45,10.55,10.65,10.75,10.85,10.95,
              11.05,11.15,11.25,11.35,11.45,11.55,11.65,11.75,11.85,11.95]

plt.plot(R_tick, avez, '.')
plt.xlabel('R (kpc)')
plt.ylabel('<z> (kpc)')
plt.show()
