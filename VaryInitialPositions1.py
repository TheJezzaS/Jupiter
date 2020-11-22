import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

mod = np.linalg.norm # Re-lable confusing syntax
years = 1000

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
    -G*Mj*1/(Vjs[1]**3) * Vjs[0][0], #second differential of y[12] (sun x)
    -G*Mj*1/(Vjs[1]**3) * Vjs[0][1], #second differential of y[13] (sun y)
    -G*Mj*1/(Vjs[1]**3) * Vjs[0][2], #second differential of y[14] (sun z)
    G*Ms*1/(Vjs[1]**3) * Vjs[0][0], #second differential of y[15] (Jupiter x)
    G*Ms*1/(Vjs[1]**3) * Vjs[0][1], #second differential of y[16] (Jupiter y)
    G*Ms*1/(Vjs[1]**3) * Vjs[0][2], #second differential of y[17] (Jupiter z)
    -G*(Ms*1/(Vsg[1]**3) * Vsg[0][0] + Mj*1/(Vjg[1]**3) * Vjg[0][0]), #second differential of y[18] (Greek x)
    -G*(Ms*1/(Vsg[1]**3) * Vsg[0][1] + Mj*1/(Vjg[1]**3) * Vjg[0][1]), #second differential of y[19] (Greek y)
    -G*(Ms*1/(Vsg[1]**3) * Vsg[0][2] + Mj*1/(Vjg[1]**3) * Vjg[0][2]), #second differential of y[20] (Greek z)
    -G*(Ms*1/(Vst[1]**3) * Vst[0][0] + Mj*1/(Vjt[1]**3) * Vjt[0][0]), #second differential of y[21] (Trojan x)
    -G*(Ms*1/(Vst[1]**3) * Vst[0][1] + Mj*1/(Vjt[1]**3) * Vjt[0][1]), #second differential of y[22] (Trojan y)
    -G*(Ms*1/(Vst[1]**3) * Vst[0][2] + Mj*1/(Vjt[1]**3) * Vjt[0][2])] #second differential of y[23] (Trojan z)

def solver(times): # Outputs (no. evaluation times) by (24) matrix, each collum, y[][i], represents y[][i]'s evolution through time
    return scipy.integrate.odeint(
        func = derivatives,
        t=times,
        # using maximum distance to sun 5.455 AU, and therefore minimum velocity 2.622 AU/yr (in sun frame) (1AU = 1.496 *10^11 m)
        y0 = (0, 0, 0, #initial values of y[0-2]; sun position
        0, 5.455, 0, #initial values of y[3-5] jupiter position
        -5.455*np.sqrt(3/4), 5.455*1/2, 0, #initial values of y[6-8]; greek position
        5.455*np.sqrt(3/4), 5.455*1/2 , 0, #initial values of y[9-11]; trojan position
        0,0,0, #initial values of y[12-14]; sun velocity
        2.622, 0, 0,#initial values of y[15-17]; jupiter velocity in AU per year
        2.622*0.5, 2.622*np.sqrt(3/4), 0, #initial values of y[18-20]; greek velocity in AU per year
        2.622*0.5, -2.622*np.sqrt(3/4), 0), #initial values of y[21-23]; trojan velocity in AU per year
    )
times = np.arange(0,years,0.1) # 10 years ≈ 1 orbit
sol = solver(times)

# Plots orbits of jupiter, and the 2 astroids
plt.plot(sol[:,3]-sol[:,0],sol[:,4]-sol[:,1], 'k', label = 'Jupiter') # Position is given in the Sun's frame, thus sun position is subtracted
plt.plot(sol[:,6]-sol[:,0],sol[:,7]-sol[:,1], label = 'Greeks')
plt.plot(sol[:,9]-sol[:,0],sol[:,10]-sol[:,1], label = 'Trojans')
sun = plt.Circle((sol[0][0],sol[0][1]), .2, color='r') # Creates a red circle in the centre of the orbits; the sun
plt.gcf().gca().add_artist(sun) # Get Current Figure and Axis, add sun

plt.annotate('The Sun',xy=(-0.7,0.3))

plt.title('Orbits of Jupiter, Greek and Trojan around the Sun')
plt.xlabel('X-Distance from the Sun in AU')
plt.ylabel('Y-Distance from the Sun in AU')
plt.legend(loc='upper left', shadow=True)

plt.show()

# Z variation
plt.plot(times,sol[:,5]-sol[:,2],'r', label='Jupiter',linewidth=4)
plt.plot(times,sol[:,8]-sol[:,2],'g', label = 'Greeks',linewidth=3)
plt.plot(times,sol[:,11]-sol[:,2],'k', label = 'Trojans',linewidth=1)
plt.xlim(0,10)
plt.legend()
plt.show()

# Test stability of Lagrange Points

# dot product rule: a.b = |a||b|cos(θ)  ⇒   θ = arccos(a.b/(|a||b|))
def thetaFinder(a,b):
    return np.arccos(np.dot(a,b)/(mod(a)*mod(b)))

thetas = np.zeros(len(times))
distances = np.zeros(len(times))
for i in range(len(times)):
    Pj = np.array([sol[i,3],sol[i,4]]) # position of jupiter
    Pg = np.array([sol[i,6],sol[i,7]]) # position of greek
    Ps = np.array([sol[i,0],sol[i,1]]) # position of sun
    thetas[i] = thetaFinder(Ps-Pj,Ps-Pg)
    distances[i] = mod(Pj-Pg)

# Find Standdard Deviation of angle
print(f'Standard Deviation of the angle between Jupiter and Greek over a period of {years} years: ', np.std(thetas))

plt.figure
plt.subplot(211)
plt.plot(times,thetas *3/np.pi, 'k') # Plots angle between the vectors: sun → jupiter and sun → greek, centered on 1 (normalised by pi/3)
plt.title('Angle Between Jupiter, the Sun and Greeks (Normalised)')
plt.xlabel('Time in Years')
plt.ylabel('Angle/(π/3)')

plt.subplot(212)
plt.plot(times, distances / distances[0], label='Distance') # Plots distance between jupiter and greek normalised reletive to original distance between them
plt.title('Distance Between Jupiter and Greeks (Normalised)')
plt.xlabel('Time in Years')
plt.xlim(0,200)

x = np.arange(0,400*np.pi,0.1)   # start,stop,step
y = 0.097*np.sin(np.pi/11.819*x)**2
plt.plot(x,1-y, 'k:', label='Theoretical Sinusoid')
plt.ylabel('Normalised Distance') # Note; never goes above 1 beacause the program starts with Jupiter at maximum distance from sun; corresponding to maximum distance to Greeks

plt.legend(loc='upper left', shadow=True)
plt.show()
