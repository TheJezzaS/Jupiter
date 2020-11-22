import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

mod = np.linalg.norm # Re-lable confusing syntax; this modulises vectors
years = 3000

def vectorise(x1,y1,z1,x2,y2,z2): #takes the position of two masses and outputs the vector difference, and the modulus
    v = np.array([x2-x1,y2-y1,z2-z1])
    return v, mod(v)

def derivatives(y,t):
    G =4*np.pi**2
    Mj = 0.001
    Ms = 1

    #Vik gives the vector pointing from i to k (leading to force on k from i)
    Vjs = vectorise(y[3],y[4],y[5],y[0],y[1],y[2])
    Vsg = vectorise(y[0],y[1],y[2],y[6],y[7],y[8])
    Vjg = vectorise(y[3],y[4],y[5],y[6],y[7],y[8])
    Vst = vectorise(y[0],y[1],y[2],y[9],y[10],y[11])
    Vjt = vectorise(y[3],y[4],y[5],y[9],y[10],y[11])

    return [y[12],y[13],y[14],#first differentials of sun position
    y[15],y[16],y[17],#first differentials of Jupiter position
    y[18],y[19],y[20],#first differentials of Greek position
    y[21],y[22],y[23], #first differentials of Trojan position
    -G*Mj*1/(Vjs[1]**3) *Vjs[0][0], #second differential of y[12] (sun x)
    -G*Mj*1/(Vjs[1]**3) *Vjs[0][1], #second differential of y[13] (sun y)
    -G*Mj*1/(Vjs[1]**3) *Vjs[0][2], #second differential of y[14] (sun z)
    G*Ms*1/(Vjs[1]**3) *Vjs[0][0], #second differential of y[15] (Jupiter x)
    G*Ms*1/(Vjs[1]**3) *Vjs[0][1], #second differential of y[16] (Jupiter y)
    G*Ms*1/(Vjs[1]**3) *Vjs[0][2], #second differential of y[17] (Jupiter z)
    -G*(Ms*1/(Vsg[1]**3) * Vsg[0][0] + Mj*1/(Vjg[1]**3) * Vjg[0][0]), #second differential of y[18] (Greek x)
    -G*(Ms*1/(Vsg[1]**3) * Vsg[0][1] + Mj*1/(Vjg[1]**3) * Vjg[0][1]), #second differential of y[19] (Greek y)
    -G*(Ms*1/(Vsg[1]**3) * Vsg[0][2] + Mj*1/(Vjg[1]**3) * Vjg[0][2]), #second differential of y[20] (Greek z)
    -G*(Ms*1/(Vst[1]**3) * Vst[0][0] + Mj*1/(Vjt[1]**3) * Vjt[0][0]), #second differential of y[21] (Trojan x)
    -G*(Ms*1/(Vst[1]**3) * Vst[0][1] + Mj*1/(Vjt[1]**3) * Vjt[0][1]), #second differential of y[22] (Trojan y)
    -G*(Ms*1/(Vst[1]**3) * Vst[0][2] + Mj*1/(Vjt[1]**3) * Vjt[0][2])] #second differential of y[23] (Trojan z)


# Outputs (no. evaluation times) by (24) matrix, each collum, y[i],
# represents y[i]'s evolution through time
def solver(times,TrojanX, TrojanY, TrojanZ, TrojanVx, TrojanVy, TrojanVz):
    return scipy.integrate.odeint(
        func = derivatives,
        t=times,
        # using maximum distance to sun 5.455 AU, and therefore minimum
        # velocity 2.622 AU/yr (in sun frame) (1AU = 1.496 *10^11 m)
        y0 = (0, 0, 0, #initial values of y[0-2]; sun position
        0, 5.455, 0, #initial values of y[3-5] jupiter position
        -5.455*np.sin(np.pi/3) , 5.455*np.cos(np.pi/3), 0,#initial values of y[6-8]; greek position
        TrojanX, TrojanY, TrojanZ, #initial values of y[9-11]; trojan position
        0,0,0,#initial values of y[12-14]; sun velocity
        2.622, 0, 0,#initial values of y[15-17]; jupiter velocity in AU per year
        2.622*np.cos(np.pi/3), 2.622*np.sin(np.pi/3), 0,#initial values of y[18-20]; greek velocity in AU per year
        TrojanVx, TrojanVy, TrojanVz),#initial values of y[21-23]; trojan velocity in AU per year
    )
times = np.arange(0,years,0.1) # 10 years aproximately = 1 orbit
