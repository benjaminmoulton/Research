import numpy as np
from matplotlib import pyplot as plt

# read in input file
with open("input.txt","r") as f:

    # read in headers
    names = f.readline().split()

    # ititialize data array
    data = []

    # read first line
    data.append([float(val) for val in f.readline().split()])
    f.readline()

    # read in data
    for line in f:
        data.append([float(val) for val in line.split()])

    # close file
    f.close()

# turn into numpy array
data = np.array(data)

c = "k"; ls = "-"; lbl = "zyl"
w = "gray"; k = "k"; b = "b"

# plot
# a = [0,0]; c = [80,80]
# plt.plot(a,c,color=w,marker="None",linestyle="-",label="White Color")
# plt.plot(a,c,color=k,marker="None",linestyle="-",label="Black Color")
# plt.plot(a,c,c=b,marker="None",linestyle="-",label="Blue  Color")
# plt.plot(a,c,color=k,marker="None",linestyle="-",label="ZYLtech PLA")
# plt.plot(a,c,color=k,marker="None",linestyle=":",label="PushPlastic PLA")
# plt.plot(a,c,color=k,marker="None",linestyle="-.",label="NinjaFlex TPU")

for i in range(len(names)-1):

    if i == 0: # Zyltech_White_PLA
        c = w; ls = "-";  lbl = "White ZYLtech PLA"

    if i == 1: # NinjaFlex_Black
        c = k; ls = "-."; lbl = "Black NinjaFlex TPU"

    if i == 2: # Zyltech_Black_PLA
        c = k; ls = "-";  lbl = "Black ZYLtech PLA"

    if i == 3: # PushPlastic_White_PLA
        c = w; ls = ":";  lbl = "White PushPlastic PLA"

    if i == 4: # Zyltech_Blue_PLA
        c = b; ls = "-";  lbl = "Blue ZYLtech PLA"

    if i == 5: # NinjaFlex_White
        c = w; ls = "-."; lbl = "White NinjaFlex TPU"

    if i == 6: # Ground
        c = "orange"; ls = "--"; lbl = "Ground"

    if i == 7: # Outside_Temp
        c = "green"; ls = "--"; lbl = "Outdoors"

    if i == 8: # Indoors
        c = "purple"; ls = "--"; lbl = "Indoors"

    if i == 9: # PLA_glass_transition_temp
        c = "red"; ls = "--"; lbl = "PLA glass transition"

    plt.plot(data[:,0],data[:,i+1],c=c,linestyle=ls,label=lbl)


plt.xlabel("Time [min]")
plt.ylabel("Temperature [F]")
plt.legend()
plt.show()