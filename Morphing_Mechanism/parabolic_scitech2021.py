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
    write(airfoil_name,xu0,yu0,xl0,yl0)
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

airfoil = "/home/ben/foilshapes/E335_200pts.dat"






airfoil = "NACA 2412"

xup,yup,xlp,ylp = main(airfoil,num_points=200)