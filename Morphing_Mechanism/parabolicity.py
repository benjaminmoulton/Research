import numpy as np
from matplotlib import pyplot as plt

# function that finds a naca 4 digit airfoil...
def get_type(airfoil_name):
    # split off numbers
    str_type = airfoil_name.split(" ")[-1]

    # split up numbers
    atype = [0,0,0]
    atype[0] = float(str_type[0])/100.
    atype[1] = float(str_type[1])/10.
    atype[2] = float(str_type[2] + str_type[3])/100.

    return atype

# create camberline of a set of x points from naca 4digit airfoil definition
def get_camberline(x,naca,c):
    # calculate camberline y values
    yc0 = np.zeros(x.shape)

    for i in range(x.shape[0]):
        if x[i] <= naca[1]:
            yc0[i] = naca[0] * (2.*(x[i]/naca[1]) - (x[i]/naca[1])**2.)
        else:
            yc0[i] = naca[0] * (2.*((c-x[i])/(c-naca[1])) - \
                ((c-x[i])/(c-naca[1]))**2.)

    return yc0
   
# create camberline of a set of x points from naca 4digit airfoil definition
def get_camberline_deriv(x,naca,c):
    # calculate camberline y values
    dyc0 = np.zeros(x.shape)

    if naca[1] != 0.0:
        for i in range(x.shape[0]):
            if x[i] <= naca[1]:
                dyc0[i] = 2.*naca[0]/naca[1]*(1-x[i]/naca[1])
            else:
                dyc0[i] = -2.*naca[0]/(c-naca[1])*(1-(c-x[i])/(c-naca[1]))

    return dyc0

# get thickness
def get_thickness(x,naca,c):
    # calculate thickness
    t = naca[2]*(2.98*(x**0.5)-1.32*x-3.286*x**2.+2.441*x**3.-0.815*x**4.)

    return t

# get radius
def get_radius(x0,yc0,hp):
    # calculate radius
    r = ( (yc0-hp[1])**2. + (x0-hp[0])**2. )**0.5

    return r

# get psi
def get_psi(x0,yc0,hp):
    # calculate psi
    psi = np.arctan((yc0-hp[1])/(x0-hp[0]))

    return psi

# get deflected camberline
def get_deflected_camberline(x0,yc0,df,psi,r,hp):
    # initialize arrays to save to
    xc = np.zeros(x0.shape)
    yc = np.zeros(x0.shape)

    # run through each point and create camberline
    for i in range(x0.shape[0]):
        if x0[i] < hp[0]:
            xc[i] =  x0[i]
            yc[i] = yc0[i]
        else:
            xc[i] = hp[0] + r[i] * np.cos(df-psi[i])
            yc[i] = hp[1] - r[i] * np.sin(df-psi[i])

    return xc,yc

# get deflected camberline derivative
def get_deflected_camberline_deriv(x0,dyc0,df,hp):
    # initialize arrays to save to
    dyc = np.zeros(x0.shape)

    # run through each point and create camberline derivative
    for i in range(x0.shape[0]):
        if x0[i] < hp[0]:
            dyc[i] = dyc0[i]
        else:
            dyc[i] = (dyc0[i]-np.tan(df))/(1+dyc0[i]*np.tan(df))

    return dyc

# get upper and lower surfaces
def get_surfaces(x,y,t,dyc):
    # calculate upper and lower surface x y values
    xu = x - t/2./(1+(dyc**2.))**0.5 * dyc
    xl = x + t/2./(1+(dyc**2.))**0.5 * dyc
    yu = y + t/2./(1+(dyc**2.))**0.5
    yl = y - t/2./(1+(dyc**2.))**0.5

    return xu,yu,xl,yl

# calculate R 
def get_R(dp):
    R = (4*np.tan(dp)**2. +1)**.5 + np.arcsinh(2*np.tan(dp))/2./np.tan(dp)
    return R

# secant method
def Secant(fun,x0,x1, Err = 1e-5, maxiter = 1000):
    
    # initialize a large number so the err is less than
    E = 10000000
    
    # initialize the icounter
    i = 0
    
    # till threshold error is reached, continue
    while(E > Err and i < maxiter):
        
        x2 = x1 - fun(x1)*(x1-x0)/(fun(x1)-fun(x0))
        
        E = abs( (x2 - x1) / x1 )
        
        # start on the next value
        x0 = x1
        x1 = x2
        
        # add counter value
        i += 1
    
    # return the root
    return x2

# find Zp
def get_Zp_Np(Z0,l,R,dp,ZTE):
    # initialize zp array
    Zparr = np.zeros(Z0.shape)

    for i in range(Z0.shape[0]):
        # define function for secant method
        def func(Zp):
            a = Zp/2.*(Zp**2./l**2.*R**2.*np.tan(dp)**2. + 1)**0.5
            b = l/2./R/np.tan(dp)*np.arcsinh(Zp/l*R*np.tan(dp))
            return a+b - Z0[i]
        
        # determine value
        Zparr[i] = Secant(func,ZTE*Z0[i]/l,0.1)
    
    # calculate Np
    Np = - Zparr**2./ZTE * np.tan(dp)

    return Zparr,Np

# calculate xp yp values
def get_xpyp(Zp,Np,hp,phi):
    # calculate values
    xp = hp[0] + Zp*np.cos(phi) - Np*np.sin(phi)
    yp = hp[1] + Zp*np.sin(phi) + Np*np.cos(phi)

    return xp,yp

# get deflected camberline
def get_parabolic_camberline(x0,yc0,xp,yp,Dyc,Zp,ZTE,dp,hp):
    # initialize arrays to save to
    xc = np.zeros(x0.shape)
    yc = np.zeros(x0.shape)

    # run through each point and create camberline
    for i in range(x0.shape[0]):
        if x0[i] < hp[0]:
            xc[i] =  x0[i]
            yc[i] = yc0[i]
        else:
            xc[i] = xp[i] + Dyc[i]*np.sin(np.arctan(2*Zp[i]/ZTE*np.tan(dp)))
            yc[i] = yp[i] + Dyc[i]*np.cos(np.arctan(2*Zp[i]/ZTE*np.tan(dp)))

    return xc,yc

# get deflected camberline derivative
def get_parabolic_camberline_deriv(x0,dyc0,Zp,ZTE,dp,hp):
    # initialize arrays to save to
    dyc = np.zeros(x0.shape)

    # run through each point and create camberline derivative
    for i in range(x0.shape[0]):
        if x0[i] < hp[0]:
            dyc[i] = dyc0[i]
        else:
            dyc[i] = (dyc0[i]-2.*Zp[i]*np.tan(dp/ZTE))/(1+2.*Zp[i]*\
                (np.tan(dp/ZTE))*dyc0[i])

    return dyc

def write(name,xu,yu,xl,yl):

    with open("/home/ben/Desktop/"+name+".dat","w") as f:
        for i in range(-1,-xu.shape[0],-1):
            f.write("{:<18.16f} {:<18.16f}\n".format(xu[i],yu[i]))
        for i in range(xl.shape[0]):
            f.write("{:<18.16f} {:<18.16f}\n".format(xl[i],yl[i]))
        f.close()
    return

def main(airfoil_name,c=1.0,xf_c=0.7,df=15.0,num_points=100):
    # determine airfoil type
    naca = get_type(airfoil_name)
    
    # initialze x0 # cosine cluster the points
    theta = np.linspace(0.0,np.pi,num_points)
    x0 = 0.5*(-np.cos(theta)) + 0.5

    # determine camberline
    chord = 1.0
    yc0 = get_camberline(x0,naca,chord)

    # determine camberline derivative
    dyc0 = get_camberline_deriv(x0,naca,chord)

    # determine thickness
    t = get_thickness(x0,naca,chord)

    # get hinge point
    yhp = get_camberline(np.array([xf_c]),naca,chord)[0]

    # set up tuple of hinge point
    hp = (xf_c,yhp)

    # determine radius function
    r = get_radius(x0,yc0,hp)

    # determine psi
    psi = get_psi(x0,yc0,hp)

    # turn degree df to radians
    df = np.deg2rad(df)

    # determine deflected camberline
    xc,yc = get_deflected_camberline(x0,yc0,df,psi,r,hp)

    # determine deflected camberline derivative
    dyc = get_deflected_camberline_deriv(x0,dyc0,df,hp)

    # calculate l
    l = (hp[1]**2. + (chord-hp[0])**2.)**0.5

    # calculate flap neutral line angle
    phi = - np.arctan2(hp[1],chord-hp[0])

    # calculate R value
    dp = df
    R = get_R(dp)

    # calculate ZTE
    ZTE = 2.*l/R

    # calculate NTE
    NTE = -2.*l/R*np.tan(dp)

    # calculate Z0
    Z0 = (x0-hp[0])/(chord-hp[0])*l

    # calculate Zp
    Zp,Np = get_Zp_Np(Z0,l,R,dp,ZTE)

    # calculate xp and yp
    xp,yp = get_xpyp(Zp,Np,hp,phi)

    # calculate ynl
    ynl = hp[1]*(1-(x0-hp[0])/(chord-hp[0]))

    # calculate Dyc
    Dyc = yc0 - ynl

    # calculate parabolic camberline
    xcp,ycp = get_parabolic_camberline(x0,yc0,xp,yp,Dyc,Zp,ZTE,dp,hp)

    # calculate parabolic camberline derivative
    dycp = get_parabolic_camberline_deriv(x0,dyc0,Zp,ZTE,dp,hp)


    # determine upper and lower surfaces
    xu0,yu0,xl0,yl0 = get_surfaces(x0 ,yc0,t,dyc0)
    # write(airfoil_name,xu0,yu0,xl0,yl0)
    # determine surfaces deflected
    xud,yud,xld,yld = get_surfaces(xc ,yc ,t,dyc )
    # determine surfaces deflected parabolically
    xup,yup,xlp,ylp = get_surfaces(xcp,ycp,t,dycp)

    x_hinge_i = np.min(np.argwhere(x0>=hp[0]))

    # # plot stuffs
    # plt.axis("equal")
    # # plot undeflected
    # color = "k"
    # gray = "0.3"
    # mid  = "0.6"
    # lin = "dashed"
    # g = " Geometry"; c = " Camber"; f = " Flap"
    # t = "Articulated"; p = "Parabolic"; a = airfoil_name
    # # plt.plot(xu0,yu0,color=color,linestyle=lin,label=a+g)
    # # plt.plot(xl0,yl0,color=color,linestyle=lin)
    # # plt.plot(x0,yc0,color=gray,linestyle=lin,label=a+c)
    # # plot deflected
    # lin = "dashdot"
    # plt.plot(xud,yud,color=color,linestyle=lin,label=t+f+g)
    # plt.plot(xld,yld,color=color,linestyle=lin)
    # plt.plot(xc,yc,color=gray,linestyle=lin,label=t+f+c)
    # # plot parabolic
    # lin = "solid"
    # plt.plot(xup,yup,color=color,linestyle=lin,label=p+f+g)
    # plt.plot(xlp,ylp,color=color,linestyle=lin)
    # plt.plot(xcp,ycp,color=gray,linestyle=lin,label=p+f+c)
    # plt.xlabel("x/c")
    # plt.ylabel("y/c")

    # plt.twinx()
    # s = " Slope"
    # lin = "dashdot"
    # plt.plot(x0,dyc0,color=mid,linestyle=lin,label=a+c+s)
    # lin = "dashed"
    # plt.plot(x0[:x_hinge_i],dyc[:x_hinge_i],color=mid,linestyle=lin,\
    #     label=t+f+c+s)
    # plt.plot(x0[x_hinge_i+1:],dyc[x_hinge_i+1:],color=mid,linestyle=lin)
    # lin = "dotted"
    # plt.plot(x0,dycp,color=mid,linestyle=lin,label=p+f+c+s)
    # plt.ylabel("dyc/c")

    # plt.legend()
    # plt.xlim(xf_c-0.1,1.0)
    # # plt.savefig(fname="/home/ben/Desktop/pbc_vs_trad.png")
    # # plt.savefig(fname="/home/ben/Desktop/pbc_vs_trad_noUND.png")
    # plt.savefig(fname="/home/ben/Desktop/cam_slope.png")
    # plt.show()

    # create section flap effectiveness plot
    cf_c = np.linspace(0.00001,1.0,num=40)

    # calculate thetaf
    th_f = np.arccos(2*cf_c-1)

    # calcualte e_ip
    e_ip = ((1+2.*np.cos(th_f)) * (np.pi-th_f) +np.sin(th_f)*(2+np.cos(th_f)))\
        / np.pi / (1+np.cos(th_f))
    
    # calculate e_if
    e_if = 1 - (th_f-np.sin(th_f)) / np.pi

    # # plot
    # c = "k"
    # marker = "o"
    # lin = ""
    # plt.plot(cf_c,e_if,c=c,marker=marker,linestyle=lin,label="Articulated")
    # marker = "s"
    # plt.plot(cf_c,e_ip,c=c,marker=marker,linestyle=lin,label="Parabolic")
    # plt.xlabel("c$_{f}$/c")
    # plt.ylabel("Ideal Flap Effectiveness")
    # plt.legend()
    # plt.savefig(fname="/home/ben/Desktop/flap_eff.png")
    # plt.show()



    return xup,yup,xlp,ylp


import C14_geoPRECISE as c14
import airfoil_db as adb
import deflect
from scipy.interpolate import interp1d as interp

# changes
# df =  -5.0 ; name =  "arcs_m05_rot.txt" #  arcs_m05  -5.0  deg
# df =  -8.5 ; name =  "arcs_m10_rot.txt" #  arcs_m10  -8.5  deg
# df = -11.5 ; name =  "arcs_m15_rot.txt" #  arcs_m15 -11.5  deg
# df =  10.25; name =  "arcs_p05_rot.txt" #  arcs_p05  10.25 deg
# df =  13.75; name =  "arcs_p10_rot.txt" #  arcs_p10  13.75 deg
# df =  17.0 ; name =  "arcs_p15_rot.txt" #  arcs_p15  17.0  deg
# df =  -4.0 ; name = "kincs_m05_rot.txt" # kincs_m05  -4.0  deg
# df = -11.0 ; name = "kincs_m10_rot.txt" # kincs_m10 -11.0  deg
# df = -13.5 ; name = "kincs_m15_rot.txt" # kincs_m15 -13.5  deg
# df =   5.25; name = "kincs_p05_rot.txt" # kincs_p05   5.25 deg
# df =   9.0 ; name = "kincs_p10_rot.txt" # kincs_p10   9.0  deg
# df =  11.5 ; name = "kincs_p15_rot.txt" # kincs_p15  11.5  deg


dfs = [ -5.0, -8.5, -11.5, 10.25, 13.75, 17.0, -4.0, -11.0, -13.5, 5.25, 9.0, 11.5]

names =  ["arcs_m05_rot.txt","arcs_m10_rot.txt","arcs_m15_rot.txt",
          "arcs_p05_rot.txt","arcs_p10_rot.txt","arcs_p15_rot.txt",
          "kincs_m05_rot.txt","kincs_m10_rot.txt","kincs_m15_rot.txt",
          "kincs_p05_rot.txt","kincs_p10_rot.txt","kincs_p15_rot.txt"]

avgs   = np.zeros((len(dfs),2))
stdevs = np.zeros((len(dfs),2))

for i in range(11,12):#len(dfs)):
    df = dfs[i]
    name = names[i]

    if name[:4] == "arcs":
        xf_c = 0.33
        Eppler = True
    elif name[:5] == "kincs":
        xf_c = 0.5
        Eppler = False
    filecompare = "/home/ben/Desktop/parabolicity/" + name
    if Eppler:
        airfoil = "/home/ben/foilshapes/E335_200pts.dat"


        ########################################
        # define parabolic shape for eppler
        add = {}
        add["geometry"] = {}
        add["geometry"]["type"] = "outline_points"
        add["geometry"]["outline_points"] = airfoil
        isNACA = False

        # initialize airfoil database
        foil = adb.Airfoil("foil",add)
        coords = foil.get_outline_points(close_te=not(isNACA),top_first=True,N=400)
        x_af = coords[:,0]; y_af = coords[:,1]

        # initialize interpolation dictionary
        I = {}

        # initialize x and y dictionaries
        x = {}; y = {}

        # create interpolation of outer surface
        O = c14.interpolate(x_af,y_af,make_camberline=True)
        # save airfoil outline for future use
        O["x af"] = x_af; O["y af"] = y_af
        # for key in O:
        #     print(key)

        info = {
            "dty" : "parabolic",
            "f" : xf_c,
            "da" : df,
            "th" : True,
            "t" : 0.8 / 12.0 / 25.4
        }
        xup,yup = deflect.main_single(O,O["tx"],O["ty"],info)
        xlp,ylp = deflect.main_single(O,O["bx"],O["by"],info)
    else:
        airfoil = "NACA 2412"

        xup,yup,xlp,ylp = main(airfoil,num_points=200,xf_c=xf_c,df=df)

    # flip upper line
    xup = np.flip(xup)
    yup = np.flip(yup)

    # combine
    xa = np.append(xup,xlp,axis=0)
    ya = np.append(yup,ylp,axis=0)

    # import comparison
    with open(filecompare,"r") as f:
        xcomp = []
        ycomp = []

        for line in f:
            split = line.split()
            xcomp.append(float(split[0]))
            ycomp.append(float(split[1]))
        
        xc = np.array(xcomp)
        yc = np.array(ycomp)

        f.close()

    # define length function
    def make_xy(x,y):

        # define l array for interpolation
        l = np.linspace(0.0,1.0,x.shape[0])

        # define interpolation functions
        l_to_x = interp(l,x)
        l_to_y = interp(l,y)
        nbeforehalf = int(x.shape[0]/2-10)
        x_to_ltop = interp(x[:nbeforehalf],l[:nbeforehalf])
        nafterhalf = int(x.shape[0]/2+10)
        x_to_lbot = interp(x[nafterhalf:],l[nafterhalf:])

        # define xy function
        def xy(lval):
            return l_to_x(lval),l_to_y(lval)
        
        # define ltop function
        def ltop(xval):
            return x_to_ltop(xval)
        
        # define lbot function
        def lbot(xval):
            return x_to_lbot(xval)

        
        return xy,ltop,lbot


    # define xy function for each
    xya,lta,lba = make_xy(xa,ya)
    xyc,ltc,lbc = make_xy(xc,yc)

    # create l array
    n = 1000
    p30 = int(lta(xf_c) * n)
    p60 = int(lba(xf_c) * n)
    l = np.linspace(0.0,1.0,n)
    # print(xya(l[p30]),xya(l[p60]))

    # spit out x and y for each
    xal,yal = xya(l)
    xcl,ycl = xyc(l)

    xdif = xal - xcl
    ydif = yal - ycl

    dif = (xdif**2. + ydif**2.)**0.5

    # print(np.sum(xdif),np.sum(ydif))
    # print(np.average(xdif),np.average(ydif))
    # print(np.std(xdif),np.std(ydif))
    favg = np.average(dif)
    fstdev = np.std(dif)
    tavg = np.average( np.append( dif[:p30], dif[p60:], axis=0) )
    tstdev = np.std(np.append(dif[:p30],dif[p60:],axis=0))
    print(name,df,"deg")
    print("full : avg : {:<11.8} x/c, 3xstdev : {:<11.8} x/c".format(\
        favg, 3*fstdev))

    print("TE   : avg : {:<11.8} x/c, 3xstdev : {:<11.8} x/c".format(\
        tavg,3*tstdev))
    print()


    plot = True
    if plot:
        # e = -0
        # plt.plot(xa[:e],ya[:e])
        # plt.plot(xc[:e],yc[:e])
        plt.plot(xa,ya,label="analytic")
        plt.plot(xc,yc,label="FDM")
        plt.axis("equal")
        plt.xlabel("x/c")
        plt.ylabel("y/c")
        plt.legend()
        plt.show()
    
    avgs[i,0] = favg
    avgs[i,1] = tavg
    stdevs[i,0] = fstdev
    stdevs[i,1] = tstdev

plotavg = False
if plotavg:
    # plot avg and stdev as function of df
    j = 6

    plt.plot([0,0],[0,0],"k",label="Whole")
    plt.plot([0,0],[0,0],"b",label="Flap")
    plt.plot([0,0],[0,0],"ko",fillstyle="none",label="ARCS")
    plt.plot([0,0],[0,0],"ks",fillstyle="none",label="KINCS")
    plt.plot(dfs[:j],avgs[:j,0],"ko",fillstyle="none")
    plt.plot(dfs[:j],avgs[:j,1],"bo",fillstyle="none")
    plt.plot(dfs[j:],avgs[j:,0],"ks",fillstyle="none")
    plt.plot(dfs[j:],avgs[j:,1],"bs",fillstyle="none")
    plt.legend()
    plt.xlabel("Flap Deflection (deg)")
    plt.ylabel("Average Difference (x/c)")
    plt.tight_layout()
    plt.show()

    plt.plot([0,0],[0,0],"k",label="Whole")
    plt.plot([0,0],[0,0],"b",label="Flap")
    plt.plot([0,0],[0,0],"ko",fillstyle="none",label="ARCS")
    plt.plot([0,0],[0,0],"ks",fillstyle="none",label="KINCS")
    plt.plot(dfs[:j],3*stdevs[:j,0],"ko",fillstyle="none")
    plt.plot(dfs[:j],3*stdevs[:j,1],"bo",fillstyle="none")
    plt.plot(dfs[j:],3*stdevs[j:,0],"ks",fillstyle="none")
    plt.plot(dfs[j:],3*stdevs[j:,1],"bs",fillstyle="none")
    plt.legend()
    plt.xlabel("Flap Deflection (deg)")
    plt.ylabel("Difference 3 x Standard Deviation (x/c)")
    plt.tight_layout()
    plt.show()

