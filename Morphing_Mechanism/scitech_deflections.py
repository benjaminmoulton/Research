import numpy as np
from matplotlib import pyplot as plt

data = [["03",	12.51689693,    21	, -13.5],
        ["04",  11.82939692 , 	22  ,	-17],
        ["04", 	11.82939692, 	33.5,	-28],
        ["05",	11.32989893, 	15	,   -26],
        ["06",	11.06718224, 	12.5,	-16],
        ["03",	12.51689693, 	20  ,	-13],
        ["04",  11.82939692 , 	22  ,	-10],
        ["04",  11.82939692 , 	23  ,	-23],
        ["05",	11.32989893, 	13  ,	-20],
        ["06",	11.06718224, 	13  , -19.5],
        ["06",	11.06718224, 	27  ,	-29],
        ["06",	11.06718224, 	27  ,	-29]]

names = [data[i][0] for i in range(len(data))]
c     = np.array([data[i][1] for i in range(len(data))])
down  = np.array([data[i][2] for i in range(len(data))])
up    = np.array([data[i][3] for i in range(len(data))])

proto = "P18I28S"

plt.plot(c,down,"ko",label="+ Deflection",fillstyle="none")
plt.plot(c,up,"o",c="0.5",label="- Deflection",fillstyle="none")

plt.text(c[-1],-5.0,proto+names[-1],rotation=90)
plt.text(c[3],-5.0,proto+names[3],rotation=90)
plt.text(c[1],-5.0,proto+names[1],rotation=90)
plt.text(c[0],-5.0,proto+names[0],rotation=90)

plt.legend()
plt.xlabel("Section Average Chord (in)")
plt.ylabel("Deflection (deg)")
plt.show()