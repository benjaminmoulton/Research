import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# initialize file name
filename = "states_ailkick_t1em1.csv"
# filename = "state_outputs_trim.csv"
# filename = "state_outputs_elevator_kick.csv"
# filename = "state_outputs_flight_test.csv"

# read in file
with open(filename,"r") as f:

    # initialize array
    states = []

    # read header line
    f.readline()

    # run through each line
    for line in f:
        
        # split the line
        states.append([ float(var) for var in line.split()])

    # close file
    f.close()

# turn into numpy array
states = np.array(states)

# determine if the aircraft ever gimbal locks
for i in range(states.shape[0]):
    if (states[i,9]*states[i,11] - states[i,10]*states[i,12] == 0.5):
        print("Gimbal Lock  0.5")
    elif (states[i,9]*states[i,11] - states[i,10]*states[i,12] == -0.5):
        print("Gimbal Lock -0.5")

# take out angle of attack
alpha = np.arctan(states[:,3] / states[:,1])

# take out sideslip angle
beta = np.arctan(states[:,2] / states[:,1])

# take out euler angles
# bank (phi)
phi = np.arctan2(2.*(states[:,9]*states[:,10] + states[:,11]*states[:,12]),\
    (states[:,9]**2. + states[:,12]**2. - states[:,10]**2. - states[:,11]**2.))

# elevation (theta)
theta = np.arcsin(2.*(states[:,9]*states[:,11] - states[:,10]*states[:,12]))

# azimuth or heading (psi)
psi = np.arctan2(2.*(states[:,9]*states[:,12] + states[:,10]*states[:,11]),\
    (states[:,9]**2. + states[:,10]**2. - states[:,11]**2. - states[:,12]**2.))



# # show craft altitude over time
# plt.plot(states[:,0],-states[:,9]-states[0,9])
# plt.xlabel("time (s)")
# plt.ylabel("altitude - {}ft (ft)".format(int(states[0,9])))
# plt.show()

# # show craft xy location # over time
# plt.plot(states[:,7],states[:,8])
# plt.xlabel("x (ft)")
# plt.ylabel("y (ft)")
# plt.show()

# show 3D craft position # over time
if False:
    fig = plt.figure()
    axe = fig.add_subplot(111,projection='3d')
    face = 'z'

    axe.plot(states[:,7],states[:,8],states[:,9],zdir=face)

    # initialize min and max vars
    xmin = np.min(states[:,7]); xmax = np.max(states[:,7])
    ymin = np.min(states[:,8]); ymax = np.max(states[:,8])
    zmin = np.min(states[:,9]); zmax = np.max(states[:,9])

    # solve for center
    xcent = (xmax + xmin) / 2.
    ycent = (ymax + ymin) / 2.
    zcent = (zmax + zmin) / 2.

    # solve for differences
    xdiff = np.abs(xmax - xmin)
    ydiff = np.abs(ymax - ymin)
    zdiff = np.abs(zmax - zmin)

    # solve for max difference
    max_diff = max([xdiff,ydiff,zdiff])

    # define limits
    x_lims = [xcent + 0.5*max_diff,xcent - 0.5*max_diff]
    y_lims = [ycent + 0.5*max_diff,ycent - 0.5*max_diff]
    z_lims = [zcent + 0.5*max_diff,zcent - 0.5*max_diff]

    # set limits
    axe.set_xlim3d(x_lims[1], x_lims[0])
    axe.set_ylim3d(y_lims[0], y_lims[1])
    axe.set_zlim3d(z_lims[1], z_lims[0])

    # set labels and finish plot
    axe.set_xlabel('x (ft)')
    axe.set_ylabel('y (ft)')
    axe.set_zlabel('z (ft)')
    # axe.view_init(200,-60)
    # axe.view_init(20,-120)
    plt.show()

# show short period mode plot
fig, axs = plt.subplots(2)
axs[0].plot(states[:,0],-states[:,9])
axs[0].set(xlabel="time (s)",ylabel="altitude (ft)")
axs[1].plot(states[:,0],np.rad2deg(alpha))
axs[1].set(xlabel="time (s)",ylabel="angle of attack (deg)")
fig.suptitle("Short Period Mode")
plt.show()

# show phugoid mode plot
fig, axs = plt.subplots(2)
axs[0].plot(states[:,0],-states[:,9])
axs[0].set(xlabel="time (s)",ylabel="altitude (ft)")
axs[1].plot(states[:,0],states[:,1])
axs[1].set(xlabel="time (s)",ylabel="body fixed airspeed 'u' (ft/s)")
fig.suptitle("Phugoid Mode")
plt.show()

# show roll mode plot
fig, axs = plt.subplots(2)
axs[0].plot(states[:,0],np.rad2deg(phi))
axs[0].set(xlabel="time (s)",ylabel="bank angle (deg)")
axs[1].plot(states[:,0],states[:,4])
axs[1].set(xlabel="time (s)",ylabel="roll rate (rad/s)")
fig.suptitle("Roll Mode")
plt.show()

# show spiral mode # over time
plt.plot(states[:,7],states[:,8])
plt.xlabel("x (ft)")
plt.ylabel("y (ft)")
plt.title("Spiral Mode")
plt.show()

# show Dutch roll plot over time
fig, (ax1, ax2, ax3) = plt.subplots(1,3)
ax1.plot(np.rad2deg(psi),states[:,0],label="azimuth")
ax1.plot(np.rad2deg(beta),states[:,0],label="sideslip")
ax1.set(xlabel="angle (deg)",ylabel="time (s)")
ax1.legend()
ax2.plot(states[:,8],states[:,0])
ax2.set(xlabel="y (ft)",ylabel="time (s)")
ax3.plot(np.rad2deg(phi),states[:,0],label="bank")
ax3.plot(np.rad2deg(beta),states[:,0],label="sideslip")
ax3.set(xlabel="angle (deg)",ylabel="time (s)")
ax3.legend()
fig.suptitle("Dutch Roll Mode")
plt.show()