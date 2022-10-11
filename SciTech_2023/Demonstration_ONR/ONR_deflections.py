import numpy as np
from matplotlib import pyplot as plt

data = [["06L",	0.9222747738938101, 28, -26],
        ["05L",	0.9441582434772273, 30, -29],
        ["04L",	0.9857830767248648, 30, -20],
        ["03L",	1.0430747438385428, 25, -26],
        ["02L",	1.1104251366970952, 25, -27],
        ["00C",	2.75,               22, -27],
        ["02R",	1.1104251366970952, 24, -28],
        ["03R",	1.0430747438385428, 27, -22],
        ["04R",	0.9857830767248648, 22, -29],
        ["05R",	0.9441582434772273, 23, -28],
        ["06R",	0.9222747738938101, 35, -33]
]

proto = "P18I28S"

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
fig.subplots_adjust(hspace=0.05)  # adjust space between axes
is_first_left = is_first_center = is_first_right = True

for i in range(len(data)):
    is_left = data[i][0][-1] == "L"
    
    if is_left:
        col = "ks"
        if is_first_left:
            lbl = "Left"
            is_first_left = False
        else:
            lbl = ""
    elif data[i][0][-1] == "C":
        col = "k*"
        if is_first_center:
            lbl = "Center"
            is_first_center = False
        else:
            lbl = ""
    else:
        col = "ko"
        if is_first_right:
            lbl = "Right"
            is_first_right = False
        else:
            lbl = ""
    
    ax1.plot(data[i][1],data[i][2],col,label=lbl,fillstyle="none")
    ax1.plot(data[i][1],data[i][3],col,fillstyle="none")
    ax2.plot(data[i][1],data[i][2],col,label=lbl,fillstyle="none")
    ax2.plot(data[i][1],data[i][3],col,fillstyle="none")

    if not(is_left):
        ax1.text(data[i][1]-0.005,-5.0,proto+data[i][0][:-1],rotation=90)
        ax2.text(data[i][1]-0.005,-5.0,proto+data[i][0][:-1],rotation=90)

# zoom in
ax1.set_xlim(.9, 1.12)
ax2.set_xlim(2.64, 2.86)
# hide the spines between ax and ax2
ax1.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax1.yaxis.tick_left()
# ax1.tick_params(labelright='off')
ax2.yaxis.tick_right()

d = .015 # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((1-d,1+d), (-d,+d), **kwargs)
ax1.plot((1-d,1+d),(1-d,1+d), **kwargs)

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d,+d), (1-d,1+d), **kwargs)
ax2.plot((-d,+d), (-d,+d), **kwargs)

# fig.xlabel("Section Average Chord (ft)")
# fig.ylabel("Deflection (deg)")
fig.text(0.5, 0.04, "Section Average Chord (ft)", ha='center')
fig.text(0.04, 0.5, "Deflection (deg)", va='center', rotation='vertical')
plt.legend()
plt.show()


quit()

names = [data[i][0] for i in range(len(data))]
c     = np.array([data[i][1] for i in range(len(data))])
down  = np.array([data[i][2] for i in range(len(data))])
up    = np.array([data[i][3] for i in range(len(data))])

plt.plot(c,down,"ko",label="+ Deflection",fillstyle="none")
plt.plot(c,up,"ks",label="- Deflection",fillstyle="none")

plt.text(c[-1],-5.0,proto+names[-1],rotation=90)
plt.text(c[3],-5.0,proto+names[3],rotation=90)
plt.text(c[1],-5.0,proto+names[1],rotation=90)
plt.text(c[0],-5.0,proto+names[0],rotation=90)
