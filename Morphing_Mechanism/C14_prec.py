import C14_geoPRECISE as c14

# run file
jsonfile = 'input_c14.json'
# vals = c14.read_json(jsonfile)
# from matplotlib import pyplot as plt
# import numpy as np
x2,y2,_ = c14.main(jsonfile)
# vals["C14"]["arc type"] = "90%t"
# x3,y3,_ = c14.main(vals)
# y3 += 0.9
# vals["C14"]["arc type"] = "hinge_point"
# x4,y4,_ = c14.main(vals)
# y4 += 1.8

# plt.axis("equal")
# c2 = "k"
# c3 = "0.3"
# c4 = "0.6"
# for i in range(x2.shape[0]):
#     if i == 0:
#         label2 = "P14I15M2 - Exponential"
#         label3 = "P14I15M3 - 90% t"
#         label4 = "P14I15M4 - hinge point"
#     else:
#         label2 = ""
#         label3 = ""
#         label4 = ""
#     plt.plot(x2[i],y2[i],c=c2,label=label2)
#     plt.plot(x3[i],y3[i],c=c3,label=label3)
#     plt.plot(x4[i],y4[i],c=c4,label=label4)
# plt.xlabel("X location (in)")
# plt.ylabel("Y location (in)")
# plt.legend()
# plt.show()
