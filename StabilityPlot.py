from CoreFunctions import *

# Test stability of Lagrange Points
# Using un-perterbed initial conditions:
sol = solver(
    times,
    5.455*np.sin(np.pi/3), 5.455*np.cos(np.pi/3) , 0, #initial values of Trojan position
    2.622*np.cos(np.pi/3), -2.622*np.sin(np.pi/3), 0 #initial values of Trojan Velocity
    ).transpose() # Creates 24 x len(times) array; each row representing one of the Ys

# dot product rule: a.b = |a||b|cos(θ)  ⇒   θ = arccos(a.b/(|a||b|))
def thetaFinder(a,b):
    return np.arccos(np.dot(a,b)/(mod(a)*mod(b)))

if __name__ == "__main__": # Allowed thetaFinder to be imported to other files without importing the plot bellow:
    thetas = np.zeros(len(times))
    distances = np.zeros(len(times))
    for i in range(len(times)):
        Pj = np.array([sol[3,i],sol[4,i]]) # position of jupiter
        Pg = np.array([sol[6,i],sol[7,i]]) # position of greek
        Ps = np.array([sol[0,i],sol[1,i]]) # position of sun
        thetas[i] = thetaFinder(Ps-Pj,Ps-Pg) # angle in radians between Jupiter, Greek and the Sun
        distances[i] = mod(Pj-Pg)


    # Find mean and Standdard Deviation of angle
    print(f'Mean of the angle between Jupiter and Greek over a period of {years} years: \
        {np.mean(thetas)*180/np.pi} degrees, which is ({np.mean(thetas)/(np.pi/3)*100 - 100})% above 60')
    print(f'Standard Deviation of the angle between Jupiter and Greek over a period of {years} years: '\
        , np.std(thetas*180/np.pi))
    print(f'The maximum angle Trojan wanders to is \
        {np.max(abs(thetas-np.pi/3))*180/np.pi} degrees')

    plt.figure()
    plt.subplot(211)
    # Plots angle between the vectors: sun → jupiter and sun → greek, centered on 1 (normalised by pi/3)
    plt.plot(times,thetas *3/np.pi, 'k')
    plt.title('Angle Between Jupiter, the Sun and Greeks (Normalised) \n')
    plt.xlabel('Time in Years')
    plt.ylabel('Angle/(π/3)')

    plt.subplot(212)
    # Plots distance between jupiter and greek normalised reletive to original distance between them
    plt.plot(times, distances / distances[0], label='Distance')
    plt.title('Distance Between Jupiter and Greeks (Normalised)')
    plt.xlabel('Time in Years')
    plt.xlim(0,200)

    x = np.arange(0,400*np.pi,0.1)   # start,stop,step
    y = 0.097*np.sin(np.pi/11.819*x)**2
    plt.plot(x,1-y, 'k:', label='Theoretical Sinusoid')
    plt.ylabel('Normalised Distance')
    # Note; never goes above 1 beacause the program starts with Jupiter at
    # maximum distance from sun; corresponding to maximum distance to Greeks

    plt.legend(loc='upper left', shadow=True)
    plt.show()

    print(max(thetas)*180/np.pi)
