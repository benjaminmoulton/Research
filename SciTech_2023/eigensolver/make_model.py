import numpy as np
import json

if __name__ == "__main__":
    # read in json
    json_string = open("make_model_input.json").read()
    input_dict = json.loads(json_string)["aircraft"]

    # create dictionary for each aircraft
    for airplane in input_dict:
        # input variables
        # AIRCRAFT GEOMETRIC PROPERTIES
        geom = input_dict[airplane]["geometry/mass"]
        S = geom["S"]
        b = geom["b"]
        c = S/b
        W = geom["W"]
        m = geom["m"]
        Ix = geom["Ix"]
        Iy = geom["Iy"]
        Iz = geom["Iz"]
        Ixz = geom["Ixz"]
        # FLIGHT CONDITION
        flight = input_dict[airplane]["flight_condition"]
        h       = flight["h"]
        M       = flight["M"]
        a       = flight["a"]
        p       = flight["p"]
        Vto     = flight["Vto"]
        qdyn    = flight["qdyn"]
        ao      = np.deg2rad(flight["ao[deg]"])
        Uo      = flight["Uo"]
        Wo      = flight["Wo"]
        deo     = np.deg2rad(flight["deo[deg]"])
        yo      = np.deg2rad(flight["yo[deg]"])
        # LONGITUDINAL
        lon = input_dict[airplane]["longitudinal"]
        Xw      = lon["Xw"]
        Xu      = lon["Xu"]
        Xde     = lon["Xde"]
        Zw      = lon["Zw"]
        Zwdot   = lon["Zwdot"]
        Zu      = lon["Zu"]
        Zde     = lon["Zde"]
        Mw      = lon["Mw"]
        Mwdot   = lon["Mwdot"]
        Mq      = lon["Mq"]
        Mu      = lon["Mu"]
        Mde     = lon["Mde"]
        # LATERAL
        lat = input_dict[airplane]["lateral"]
        Yv      = lat["Yv"]
        Ydastar = lat["Ydastar"]
        Ydrstar = lat["Ydrstar"]
        Lbprime = lat["Lbprime"]
        Lpprime = lat["Lpprime"]
        Lrprime = lat["Lrprime"]
        Ldaprime = lat["Ldaprime"]
        Ldrprime = lat["Ldrprime"]
        Nbprime = lat["Nbprime"]
        Npprime = lat["Npprime"]
        Nrprime = lat["Nrprime"]
        Ndaprime = lat["Ndaprime"]
        Ndrprime = lat["Ndrprime"]
        # OTHER PARTICULARS
        other = input_dict[airplane]["other"]
        name = other["name"]
        short_name = other["short_name"]
        file_name = "aircraft_database/" + other["file_name"]
        from_string = other.get("from","NASA CR-96008 -- Teper, Gary L. 1969")
        from_file = other.get("file","https://ntrs.nasa.gov/api/citations/" + \
                "19690022405/downloads/19690022405.pdf")
        notes = other.get("notes","")
        rho = 2.37689261479e-03
        # startoff zeros
        # LONGITUDINAL
        CL_q = 0.0
        CL_adot = 0.0
        CL_mudot = 0.0
        CD_q = 0.0
        CD_adot = 0.0
        CD_mudot = 0.0
        CD_Lq = 0.0
        CD_L2q = 0.0
        Cmo = 0.0
        Cm_mudot = 0.0
        # LATERAL
        Cy_p = 0.0
        Cy_Lp = 0.0
        Cy_r = 0.0
        Cn_Lp = 0.0
        Cl_Lr = 0.0

        # determine stability derivatives, assuming mach derivs to be zero
        # determine body-fixed axis derivs and values
        # LONGITUDINAL
        A = [[-1.,Wo/2./Uo],
            [-2.*Wo/Uo,-1.]]
        B = [[Xu/(p*S*Uo/m)],
            [Xw/(p*S*Uo/2./m)]]
        CX,CX_a  = np.matmul(np.linalg.inv(A),B)[:,0].tolist()
        A = [[-1.,Wo/2./Uo],
            [-2.*Wo/Uo,-1.]]
        B = [[Zu/(p*S*Uo/m)],
            [Zw/(p*S*Uo/2./m)]]
        CN,CN_a  = np.matmul(np.linalg.inv(A),B)[:,0].tolist()
        CN_adot = Zwdot / (-p*S*c/4./m * Uo/Vto)
        A = [[1.,-Wo/2./Uo],
            [2.*Wo/Uo,1.]]
        B = [[Mu/(p*S*c*Uo/2./Iy)],
            [Mw/(p*S*c*Uo/2./Iy)]]
        Cm,Cm_a  = np.matmul(np.linalg.inv(A),B)[:,0].tolist()
        Cm_adot = Mwdot / (p*S*c*c/4./Iy * Uo/Vto)
        Cm_q = Mq / (p*S*c*c*Vto/4./Iy)
        CX_de = Xde / (-p*S*Vto**2./2./m)
        CN_de = Zde / (-p*S*Vto**2./2./m)
        Cm_de = Mde / (p*S*c*Vto**2./2./Iy)
        # LATERAL
        #    Convert out of inertial ratioing
        if Ix != 0.0 and Iz != 0.0:
            G = 1. / (1. - Ixz**2./Ix/Iz)
        else:
            G = 1.0
        if Ix == 0.0: Ix = 1e-16
        if Iz == 0.0: Iz = 1e-16
        A = [[1.,Ixz/Ix],
            [Ixz/Iz,1.]]
        B = [[Lbprime/G],
            [Nbprime/G]]
        Lb,Nb  = np.matmul(np.linalg.inv(A),B)[:,0].tolist()
        A = [[1.,Ixz/Ix],
            [Ixz/Iz,1.]]
        B = [[Lpprime/G],
            [Npprime/G]]
        Lp,Np  = np.matmul(np.linalg.inv(A),B)[:,0].tolist()
        A = [[1.,Ixz/Ix],
            [Ixz/Iz,1.]]
        B = [[Lrprime/G],
            [Nrprime/G]]
        Lr,Nr  = np.matmul(np.linalg.inv(A),B)[:,0].tolist()
        A = [[1.,Ixz/Ix],
            [Ixz/Iz,1.]]
        B = [[Ldaprime/G],
            [Ndaprime/G]]
        Lda,Nda  = np.matmul(np.linalg.inv(A),B)[:,0].tolist()
        A = [[1.,Ixz/Ix],
            [Ixz/Iz,1.]]
        B = [[Ldrprime/G],
            [Ndrprime/G]]
        Ldr,Ndr  = np.matmul(np.linalg.inv(A),B)[:,0].tolist()
        #    Convert to body fixed
        Cy_b__bf = Yv / (p*S*Vto/2./m)
        Cl_b__bf = Lb / (p*S*Vto**2.*b/2./Ix)
        Cl_p__bf = Lp / (p*S*Vto*b**2./4./Ix)
        Cl_r__bf = Lr / (p*S*Vto*b**2./4./Ix)
        Cn_b__bf = Nb / (p*S*Vto**2.*b/2./Iz)
        Cn_p__bf = Np / (p*S*Vto*b**2./4./Iz)
        Cn_r__bf = Nr / (p*S*Vto*b**2./4./Iz)
        Cy_da__bf = Ydastar / (p*S*Vto/2./m)
        Cy_dr__bf = Ydrstar / (p*S*Vto/2./m)
        Cl_da__bf = Lda / (p*S*Vto**2.*b/2./Ix)
        Cl_dr__bf = Ldr / (p*S*Vto**2.*b/2./Ix)
        Cn_da__bf = Nda / (p*S*Vto**2.*b/2./Iz)
        Cn_dr__bf = Ndr / (p*S*Vto**2.*b/2./Iz)
        

        # determine stability-axis derivs
        # LONGITUDINAL
        ca = np.cos(ao); sa = np.sin(ao)
        A = [[ ca, sa],
            [-sa, ca]]
        B = [[CN],
            [CX]]
        CLo,CDo  = np.matmul(np.linalg.inv(A),B)[:,0].tolist()
        CL_ahat = CN_adot / ca
        A = [[ ca, sa],
            [-sa, ca]]
        B = [[CN_a + CLo*sa - CDo*ca],
            [CX_a + CDo*sa + CLo*ca]]
        CL_a,CD_a  = np.matmul(np.linalg.inv(A),B)[:,0].tolist()
        A = [[ ca, sa],
            [-sa, ca]]
        B = [[CN_de],
            [CX_de]]
        CL_de,CD_de  = np.matmul(np.linalg.inv(A),B)[:,0].tolist()
        #   Split CDo and CD_a to CD0,CD1,CD2;  assume CD1 = 0.0
        CD1 = 0.0
        if CLo != 0.0 and CL_a != 0.0:
            CD2 = ( CD_a - CD1*CL_a ) / (2. * CLo * CL_a)
        else:
            CD2 = 0.0
        CD0 = CDo - CD1 * CLo - CD2 * CLo*CLo
        # LATERAL
        #    Convert from body-fixed to stability
        Cy_b = Cy_b__bf * 1.
        A = [[ ca,-sa],
            [ sa, ca]]
        B = [[Cl_b__bf],
            [Cn_b__bf]]
        Cl_b,Cn_b  = np.matmul(np.linalg.inv(A),B)[:,0].tolist()
        Cy_da = Cy_da__bf * 1.
        A = [[ ca,-sa],
            [ sa, ca]]
        B = [[Cl_da__bf],
            [Cn_da__bf]]
        Cl_da,Cn_da  = np.matmul(np.linalg.inv(A),B)[:,0].tolist()
        Cy_dr = Cy_dr__bf * 1.
        A = [[ ca,-sa],
            [ sa, ca]]
        B = [[Cl_dr__bf],
            [Cn_dr__bf]]
        Cl_dr,Cn_dr  = np.matmul(np.linalg.inv(A),B)[:,0].tolist()
        A = [[ ca*ca,-sa*ca,-sa*ca, sa*sa],
            [ sa*ca, ca*ca,-sa*sa,-sa*ca],
            [ sa*ca,-sa*sa, ca*ca,-sa*ca],
            [ sa*sa, sa*ca, sa*ca, ca*ca]]
        B = [[Cl_p__bf],
            [Cl_r__bf],
            [Cn_p__bf],
            [Cn_r__bf]]
        Cl_p,Cl_r,Cn_p,Cn_r  = np.matmul(np.linalg.inv(A),B)[:,0].tolist()


        # CLo,CDo = CN,CX
        # CL_ahat = CN_adot
        # CL_a,CD_a = CN_a,CX_a
        # CL_de,CD_de = CN_de,CX_de
        # Cl_b,Cn_b  = Cl_b__bf,Cn_b__bf
        # Cl_da,Cn_da  = Cl_da__bf,Cn_da__bf
        # Cl_dr,Cn_dr  = Cl_dr__bf,Cn_dr__bf
        # Cl_p,Cl_r,Cn_p,Cn_r  = Cl_p__bf,Cl_r__bf,Cn_p__bf,Cn_r__bf

        # make file 
        craft = {}
        #    aircraft components
        craft["aircraft"] = {}
        craft["aircraft"]["name"] = name
        craft["aircraft"]["short_name"] = short_name
        craft["aircraft"]["wing_area[ft^2]"] = S
        craft["aircraft"]["wing_span[ft]"] = b
        craft["aircraft"]["weight[lbf]"] = W
        craft["aircraft"]["Ixx[slug-ft^2]"] = Ix
        craft["aircraft"]["Iyy[slug-ft^2]"] = Iy
        craft["aircraft"]["Izz[slug-ft^2]"] = Iz
        craft["aircraft"]["Ixz[slug-ft^2]"] = Ixz
        craft["aircraft"]["Ixy[slug-ft^2]"] = 0.0
        craft["aircraft"]["Iyz[slug-ft^2]"] = 0.0
        craft["aircraft"]["from"] = from_string
        craft["aircraft"]["file"] = from_file
        craft["aircraft"]["notes"] = notes
        #   analysis components
        craft["analysis"] = {}
        craft["analysis"]["density[slugs/ft^3]"] = rho
        #   aerodynamics components
        craft["aerodynamics"] = {}
        #       CL components
        craft["aerodynamics"]["CL"] = {}
        craft["aerodynamics"]["CL"]["0"] = CLo
        craft["aerodynamics"]["CL"]["alpha"] = CL_a
        craft["aerodynamics"]["CL"]["qbar"] = CL_q
        craft["aerodynamics"]["CL"]["alpha_hat"] = CL_adot
        craft["aerodynamics"]["CL"]["mu_hat"] = CL_mudot
        craft["aerodynamics"]["CL"]["de"] = CL_de
        #       CY components
        craft["aerodynamics"]["CS"] = {}
        craft["aerodynamics"]["CS"]["beta"] = Cy_b
        craft["aerodynamics"]["CS"]["pbar"] = Cy_p
        craft["aerodynamics"]["CS"]["Lpbar"] = Cy_Lp
        craft["aerodynamics"]["CS"]["rbar"] = Cy_r
        craft["aerodynamics"]["CS"]["da"] = Cy_da
        craft["aerodynamics"]["CS"]["dr"] = Cy_dr
        #       CD components
        craft["aerodynamics"]["CD"] = {}
        craft["aerodynamics"]["CD"]["L0"] = CD0
        craft["aerodynamics"]["CD"]["L"]  = CD1
        craft["aerodynamics"]["CD"]["L2"] = CD2
        craft["aerodynamics"]["CD"]["qbar"]   = CD_q
        craft["aerodynamics"]["CD"]["Lqbar"]  = CD_Lq
        craft["aerodynamics"]["CD"]["L2qbar"] = CD_L2q
        craft["aerodynamics"]["CD"]["alpha_hat"] = CD_adot
        craft["aerodynamics"]["CD"]["mu_hat"] = CD_mudot
        craft["aerodynamics"]["CD"]["de"] = CD_de
        #       Cl components
        craft["aerodynamics"]["Cl"] = {}
        craft["aerodynamics"]["Cl"]["beta"] = Cl_b
        craft["aerodynamics"]["Cl"]["pbar"] = Cl_p
        craft["aerodynamics"]["Cl"]["rbar"] = Cl_r
        craft["aerodynamics"]["Cl"]["Lrbar"] = Cl_Lr
        craft["aerodynamics"]["Cl"]["da"] = Cl_da
        craft["aerodynamics"]["Cl"]["dr"] = Cl_dr
        #       Cm components
        craft["aerodynamics"]["Cm"] = {}
        craft["aerodynamics"]["Cm"]["0"] = Cmo
        craft["aerodynamics"]["Cm"]["alpha"] = Cm_a
        craft["aerodynamics"]["Cm"]["qbar"] = Cm_q
        craft["aerodynamics"]["Cm"]["alpha_hat"] = Cm_adot
        craft["aerodynamics"]["Cm"]["mu_hat"] = Cm_mudot
        craft["aerodynamics"]["Cm"]["de"] = Cm_de
        #       Cn components
        craft["aerodynamics"]["Cn"] = {}
        craft["aerodynamics"]["Cn"]["beta"] = Cn_b
        craft["aerodynamics"]["Cn"]["pbar"] = Cn_p
        craft["aerodynamics"]["Cn"]["Lpbar"] = Cn_Lp
        craft["aerodynamics"]["Cn"]["rbar"] = Cn_r
        craft["aerodynamics"]["Cn"]["da"] = Cn_da
        craft["aerodynamics"]["Cn"]["dr"] = Cn_dr
        craft["handling_qualities"] =input_dict[airplane]["handling_qualities"]

        # save to json file
        with open(file_name, "w") as outfile:
            json.dump(craft, outfile, indent = 4)
    
    # report files created
    print("created files for:",end="")

    for airplane in input_dict:
        print("\"{0}\", ".format(input_dict[airplane]["other"]["file_name"]),\
                end="")
    print()




