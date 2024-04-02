import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

### Insert the raw file from Gaia here
data = np.loadtxt('R8.0_10.0L225_245_2024.csv', delimiter = ',', skiprows = 1)

### Assign values based on your query
l = data[:,0]
b = data[:,1]
w = data[:,2]
G = data[:,3] # apparent magnitude
bp_rp = data[:,4]

d = 1/w
MG = G - 5*np.log10(d/0.01) # absolute magnitude

### GC coordinate
x = d*np.cos(b*(np.pi/180))*np.cos(l*(np.pi/180))-8
y = d*np.cos(b*(np.pi/180))*np.sin(l*(np.pi/180))
R = np.sqrt(x*x+y*y)

### Calculate vertical height and phi
z = d*np.sin((b*(np.pi/180))) # does not take solar height into account

### Select stars within R = 8.0 to 10.0 kpc
nbp_rp = bp_rp[(R<10)*(R>8)]
nMG = MG[(R<10)*(R>8)]
nz = z[(R<10)*(R>8)]
nR = R[(R<10)*(R>8)]


##### Make a CM Diagram to help decide CMD cuts. Please comment out
##### all the sections below before running this CMD

plt.plot(nbp_rp,nMG,',', markersize = 0.3) 
plt.gca().invert_yaxis()
plt.title('CMD Diagram')
plt.xlabel('bp-rp (mag)')
plt.ylabel('Absolute Magnitude (mag)')
plt.show()

### The upperline and lower line here correspond to the CMD cuts;
### in this case we mainly selected stars on the main sequence. 
lowerline = nbp_rp*(3.4251)+(1.517) # CMD cuts by your choice
upperline = nbp_rp*(3.8557)-(0.4312) # CMD cuts by your choice
zf = nz[(nMG<lowerline)*(nMG>upperline)] # z_final values
Rf = nR[(nMG<lowerline)*(nMG>upperline)] # R_final values


print ('Number of stars after the CMD and R range cuts = ', len(Rf))

dataf = {'z': zf, 'R': Rf}
df = pd.DataFrame(dataf, columns = ['z','R'])
df.to_csv(r'C:\Users\scott\Desktop\Reproduce_the_Gaia_Research\R8.0_10.0_initial_cuts.csv',
          index = False, header = True) # change the output directory to your own directory


#############################################################
##### Cut the stars from the initial_cuts.csv file into bins of R = 0.1 kpc
data = np.loadtxt('R8.0_10.0_initial_cuts.csv', delimiter = ',', skiprows = 1)
z = data[:,0]
R = data[:,1]

#### cut into wedges
zarry = np.ones((20,3000000))*(-999.7) ### -999.7 are merely place holders
Rmin = 8
Rstep = 0.1
for i in range (0, 20):
    ztemp = z[(R>(Rmin + i*Rstep))*(R<(Rmin + (i+1)*Rstep))]
    for j in range (0, len(ztemp)):
        zarry[i,j] = ztemp[j]

#### export
# Remember to change the prefix to your own directory
prefix = r'C:\Users\scott\Desktop\Reproduce_the_Gaia_Research'
postfix = r'_trimmed.csv'
list = [r'\8_8.1', r'\8.1_8.2', r'\8.2_8.3', r'\8.3_8.4', r'\8.4_8.5',
        r'\8.5_8.6', r'\8.6_8.7', r'\8.7_8.8', r'\8.8_8.9', r'\8.9_9.0',
        r'\9.0_9.1', r'\9.1_9.2', r'\9.2_9.3', r'\9.3_9.4', r'\9.4_9.5',
        r'\9.5_9.6', r'\9.6_9.7', r'\9.7_9.8', r'\9.8_9.9', r'\9.9_10.0'] # change file names as necessary

for i in range (0, len(list)):
    dataf = {'z': zarry[i]}
    df = pd.DataFrame(dataf, columns = ['z'])
    df.to_csv(prefix + list[i] + postfix, index = False, header = True)





        

