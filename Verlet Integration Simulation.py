import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

### Constants
G = 4.3009e-3 ### gravitational constant in pc (km/s)^2 M_sun^-1
TIMESTEP = 3.156e15 ### time step in seconds
SIZE = 5
M = [100, 200, 150, 300, 500] ### in solar mass
X = [10, 20, 30, 40, 50] ### in pc
Z = [-30, 20, 10, 30, -10] ### initial conditions in z direction (pc)
VZINIT = [0, 0, 0, 0, 0] ### in km/s
iteration = 10
TimeTick = np.linspace(0, iteration*TIMESTEP, iteration+2)
BM = sum(M) * 100 ### background mass in solar mass

### Constants for background potential. May subject to change.
a = 2700 #pc
b = 200 #pc

##def BP (mass, z, x): ### background potential for mass M at position R
##    return -G*mass/sp.sqrt(x**2 + (sp.sqrt(z**2 + b**2) + a**2))



### Find the z partial derivative expression for two masses
x1, x2, z1, z2, MM= sp.symbols('x1 x2 z1 z2 MM', positive = True, real = True)
pot = -(G * MM) / sp.sqrt((x2 - x1) ** 2 + (z2 - z1) ** 2) ### in km^2/s^2
a_dz =  (- sp.diff(pot, z1))/(3.0856 * 10**13) ### partial derivative for z component accel in km/s^2
a_BM_dz = (- sp.diff((-G*BM/sp.sqrt(x1**2 + (sp.sqrt(z1**2 + b**2) + a)**2)),z1))/(3.0856 * 10**13) ### accel due to background mass
print (a_dz)


def accel_z (x1_, x2_, z1_, z2_, m_):
    result = sp.N(a_dz.subs([(x1, x1_), (x2, x2_), (z1, z1_), (z2, z2_), (MM, m_)]))
    return result

def Zaccel (listZ: list[int]):
    accel = [] # acceleration components for the ith mass (e.g. M1 + M2... for M0)
    acceleration = []  # z direction acceleration in km/s^2
    for i in range(0, SIZE):
        for j in range(0, SIZE):
            if i != j:
                accel.append(accel_z(X[i], X[j], listZ[i], listZ[j], M[j]))
        acceleration.append(sum(accel) + sp.N(a_BM_dz.subs([(x1, X[i]), (x2, 0), (z1, listZ[i]), (z2, 0), (MM, 0)])))
        accel = [] 
    return acceleration # the ith term is the acceleration for the ith mass


### Find the positions of the first step for each mass
def firstStep ():
    acceleration = Zaccel(Z)
    ret = []
    for i in range(0, SIZE):
        ret.append(Z[i] + VZINIT[i]/(3.086e13) * TIMESTEP + 0.5 * acceleration[i]/(3.086e13) * TIMESTEP ** 2)
    return ret
init_z = firstStep ()


### Find the position of every step via Verlet Integration
def TimeStep (z_prev, z_curr):
    positions = [Z, init_z]
    currAccel = Zaccel(init_z) ### z direction current acceleration of a planet
    for i in range (iteration): ### number of iterations
        z_next = []
        for j in range (SIZE):
            z_next.append(2 * z_curr[j] - z_prev[j] + currAccel[j]/(3.086e13) * TIMESTEP ** 2)
            
        # update position
        positions.append(z_next)
        ## update variables for next iteration
        currAccel = Zaccel(z_next)
        z_prev = z_curr
        z_curr = z_next
    return positions

PositionArray = TimeStep (Z, init_z)

print (PositionArray)

###Plot every step of Verlet integration
##plt.plot(TimeTick,[item[0] for item in PositionArray])
##plt.show()

