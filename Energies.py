from CoreFunctions import *

# Using un-perterbed initial conditions:
sol = solver(
    times,
    #initial values of Trojan position:
    5.455*np.sin(np.pi/3), 5.455*np.cos(np.pi/3) , 0,
    #initial values of Trojan Velocity:
    2.622*np.cos(np.pi/3), -2.622*np.sin(np.pi/3), 0
     # Creates 24 x len(times) array; each row representing one of the Ys:
    ).transpose()

# Plot graph of energies to check conservation
potentialE = [-4*np.pi**2*0.001*1 / (vectorise(sol[0,i],sol[1,i],sol[2,i],sol[3,i],sol[4,i],sol[5,i])[1]) for i in range(len(times))]
kineticE = [1/2 * 0.001 * vectorise(sol[12,i],sol[13,i],sol[14,i],sol[15,i],sol[16,i],sol[17,i])[1]**2 for i in range(len(times))]
totalE  = [potentialE[i]+kineticE[i] for i in range(len(times))]

fig, (ax1, ax3) = plt.subplots(2,1)
plt.title('Total Energy over Time')
color = 'tab:red'
ax1.set_xlabel('Time in Years')
ax1.set_ylabel('Potential Energy in M☉ AU² /years²', color=color)
ax1.plot(times, potentialE, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'tab:blue'
ax2.set_ylabel('Kinetic Energy in M☉ AU² /years²', color=color)
ax2.plot(times, kineticE, color=color, linestyle='--')
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.xlim(0,35.58)

ax3.plot(times,totalE) # Should be negative for bound orbit
plt.title('Potential and Kinetic Energies over Time')
ax3.set_xlabel('Time in Years')
ax3.set_ylabel('Total Energy in M☉ AU² /years²')
ax3.set_xlim(0,35.58)

plt.show()
# note; sinusoid stays in a range.

print(f'The maximum total energy is {np.max(totalE)}')
