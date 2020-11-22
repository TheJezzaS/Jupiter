from CoreFunctions import *
from StabilityPlot import thetaFinder
import os
from time import time # Allows run time to be found
start_time = time() # records time at start
# Test stability of Lagrange Points
steps = 101
maxAngleList1 = []
Perterbation = np.linspace(-0.2,0.2,steps)

maxAngleList2 = []
anglePerterbed = np.linspace(-np.pi*5/9,np.pi*5/9,steps)

maxAngleList3 = []
maxAngleList4 = []

for dDistance in Perterbation:
    PerterbedSol = solver( # Perterbed solution
        times,
        5.455*np.sin(np.pi/3)+dDistance, 5.455*np.cos(np.pi/3)+dDistance, 0, # initial values of Trojan position
        2.622*np.cos(np.pi/3), -2.622*np.sin(np.pi/3), 0 # initial values of Trojan Velocity
        ).transpose() # Creates 24 x len(times) array; each row representing one of the Ys
    JSTthetas = np.zeros(len(times))
    for i in range(len(times)):
        Pj = np.array([PerterbedSol[3,i],PerterbedSol[4,i]]) # positions of jupiter through time
        Pt = np.array([PerterbedSol[9,i],PerterbedSol[10,i]]) # positions of trojan through time
        Ps = np.array([PerterbedSol[0,i],PerterbedSol[1,i]]) # positions of sun through time
        JSTthetas[i] = thetaFinder(Ps-Pj,Ps-Pt)
    maxAngleList1.append( [np.max(JSTthetas)*180/np.pi])

for dAngle in anglePerterbed:
    PerterbedSol = solver( # Perterbed solution
        times,
        5.455*np.sin(np.pi/3+dAngle), 5.455*np.cos(np.pi/3+dAngle), 0, #initial values of Trojan position
        2.622*np.cos(np.pi/3), -2.622*np.sin(np.pi/3), 0 #initial values of Trojan Velocity
        ).transpose() # Creates 24 x len(times) array; each row representing one of the Ys

    JSTthetas = np.zeros(len(times))
    for i in range(len(times)):
        Pj = np.array([PerterbedSol[3,i],PerterbedSol[4,i]]) # position of jupiter
        Pt = np.array([PerterbedSol[9,i],PerterbedSol[10,i]]) # position of trojan
        Ps = np.array([PerterbedSol[0,i],PerterbedSol[1,i]]) # position of sun
        JSTthetas[i] = thetaFinder(Ps-Pj,Ps-Pt)
    maxAngleList2.append( np.max(JSTthetas)*180/np.pi )

for dVx in Perterbation:
    PerterbedSol = solver( # Perterbed solution
        times,
        5.455*np.sin(np.pi/3), 5.455*np.cos(np.pi/3), 0, # initial values of Trojan position
        2.622*np.cos(np.pi/3)+dVx, -2.622*np.sin(np.pi/3), 0 # initial values of Trojan Velocity
        ).transpose() # Creates 24 x len(times) array; each row representing one of the Ys
    JSTthetas = np.zeros(len(times))
    for i in range(len(times)):
        Pj = np.array([PerterbedSol[3,i],PerterbedSol[4,i]]) # positions of jupiter through time
        Pt = np.array([PerterbedSol[9,i],PerterbedSol[10,i]]) # positions of trojan through time
        Ps = np.array([PerterbedSol[0,i],PerterbedSol[1,i]]) # positions of sun through time
        JSTthetas[i] = thetaFinder(Ps-Pj,Ps-Pt)
    maxAngleList3.append( [np.max(JSTthetas)*180/np.pi])

for dVy in Perterbation:
    PerterbedSol = solver( # Perterbed solution
        times,
        5.455*np.sin(np.pi/3), 5.455*np.cos(np.pi/3), 0, # initial values of Trojan position
        2.622*np.cos(np.pi/3), -2.622*np.sin(np.pi/3)+dVy, 0 # initial values of Trojan Velocity
        ).transpose() # Creates 24 x len(times) array; each row representing one of the Ys
    JSTthetas = np.zeros(len(times))
    for i in range(len(times)):
        Pj = np.array([PerterbedSol[3,i],PerterbedSol[4,i]]) # positions of jupiter through time
        Pt = np.array([PerterbedSol[9,i],PerterbedSol[10,i]]) # positions of trojan through time
        Ps = np.array([PerterbedSol[0,i],PerterbedSol[1,i]]) # positions of sun through time
        JSTthetas[i] = thetaFinder(Ps-Pj,Ps-Pt)
    maxAngleList4.append( [np.max(JSTthetas)*180/np.pi])


print('This took', time()-start_time, 'to run')

plt.figure()
plt.subplot(211)
plt.plot(np.sqrt(2)*Perterbation, maxAngleList1) # Should stay at 60, then shoot up to 180
plt.ylabel('Max Angle in degrees')
plt.xlabel('Perterbation Distance, Radially out from the Sun, in AU')
plt.title(f'Maximum Angle of Wander Plotted Against Original Perterbation Distance from the Lagrange Point, over {years} Years')

plt.subplot(212)
plt.plot(anglePerterbed*180/np.pi,maxAngleList2) # Should stay at 60, then shoot up to 180
plt.ylabel('Max Angle in degrees')
plt.xlabel('Original Perterbation Angle, Clockwise in Degrees')
plt.title(f'Maximum Angle of Wander Plotted Against Original Perterbation Angle from the Lagrange Point, over {years} Years')

os.system('say "Done"') # Program is slow, so announce when finished
plt.show()



plt.plot(Perterbation, maxAngleList3, label = 'Velocity Perterbation in x Direction',linewidth = 2) # Should stay at 60, then shoot up to 180
plt.plot(Perterbation,maxAngleList4, label = 'Velocity Perterbation in y Direction') # Should stay at 60, then shoot up to 180

plt.ylabel('Max Angle in degrees')
plt.xlabel('Velocity, in AU/year')
plt.title(f'Maximum Angle of Wander Plotted Against Velocity Perterbation, over {years} Years')
plt.legend(loc='lower left')

plt.show()
