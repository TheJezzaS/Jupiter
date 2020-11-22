from CoreFunctions import *

# Using un-perterbed initial conditions:
sol = solver(
    times,
    5.455*np.sin(np.pi/3) , 5.455*np.cos(np.pi/3) , 0, #initial values of Trojan position
    2.622*np.cos(np.pi/3), -2.622*np.sin(np.pi/3), 0 #initial values of Trojan Velocity
    ).transpose() # Creates 24 x len(times) array; each row representing one of the Ys

#### Plots orbits of jupiter, and the 2 astroids ####


plt.plot(sol[3]-sol[0],sol[4]-sol[1], 'k', label = 'Jupiter') # Position is given in the Sun's frame, thus sun position is subtracted
plt.plot(sol[6]-sol[0],sol[7]-sol[1], label = 'Greeks')
plt.plot(sol[9]-sol[0],sol[10]-sol[1], label = 'Trojans')
sun = plt.Circle((sol[0][0],sol[1][0]), .2, color='r') # Creates a red circle at a focus of the orbits; the sun
plt.gcf().gca().add_artist(sun) # Get Current Figure and Axis, add sun
plt.annotate('The Sun',xy=(-0.7,0.3))

plt.title('Orbits of Jupiter, Greek and Trojan around the Sun')
plt.xlabel('X-Distance from the Sun in AU')
plt.ylabel('Y-Distance from the Sun in AU')
plt.legend(loc='upper right', shadow=True)
plt.plot([0,-6],[5.46,5.46],'y:')
plt.plot([0,-6],[-4.92,-4.92],'y:')
plt.yticks(list(plt.yticks()[0]) + [5.46, -4.92])

plt.show()

# # Z variation Plotter
# plt.plot(times,sol[5]-sol[2],'r', label='Jupiter',linewidth=4)
# plt.plot(times,sol[8]-sol[2],'g', label = 'Greeks',linewidth=3)
# plt.plot(times,sol[11]-sol[2],'k', label = 'Trojans',linewidth=1)
# # plt.ylim(-0.1,0.1)
# plt.legend()
# plt.title('Z Variation in Time')
# plt.xlabel('Time in Years')
# plt.ylabel('Z-Distance from the Sun in AU')
# plt.show()

print(f'The range of Z variation of Jupiter is {np.min(sol[5]-sol[2])} ≤ z ≤ {np.max(sol[5]-sol[2])}')
print(f'The range of Z variation of Greek is {np.min(sol[8]-sol[2])} ≤ z ≤ {np.max(sol[8]-sol[2])}')
print(f'The range of Z variation of Trojan is {np.min(sol[11]-sol[2])} ≤ z ≤ {np.max(sol[11]-sol[2])}')
