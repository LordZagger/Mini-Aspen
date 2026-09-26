# this module contains function I think could be one day useful, and other random stuff like tips and constants
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp
import pandas as pd
import io

'''
Format of a good docstring:
    1. Type contract
    2. Descriptions
    3. Preconditions
    4. Examples (2-4)
'''

sp.init_printing()
x, y, z, t, T, u, v, r, s, p, k = sp.symbols(
    'x y z t T u v r s p k')  # for functions requiring sympy
pi = sp.pi #simpler to write
inf = sp.S.Infinity #simpler to write, positive infinity

def washing_mccabe_thiele(U, O, x_out, x_in, yNplus1=0):
    '''
    Warning: this is for countercurrent flow!
    U - underflow liquid's mass flow rate
    O - overflow liquid's mass flow rate
    x_out - xN, desired purity
    x_in - x0, salt mass fraction of source entrained liquid (underflow)
    y_N+1 - y_in, salt mass fraction of source wash liquid (overflow), typically 0
 
    Prints a table of mass fractions at each stage, using the algebraic method to calculate all mass fractions
    And prints required number of stages to achieve desired purity
    Then presents a McCabe-Thiele stepping graph corresponding to the printed table
    Please use floats for inputs!
    '''
    plot_xs = [x_in]
    plot_ys = []
    print("[i, x_i, y_i+1]")
    x_i = x_in  # x0
    i = 0
    while x_i > x_out: #stepping down
        y_iplus1 = U/O*x_i + yNplus1 - U/O * x_out  # operating line (to get y_i+1)
        plot_ys.append(y_iplus1)
        current_row = [i, x_i, y_iplus1]
        print(current_row)
        i += 1
        x_i = y_iplus1  # equilibrium line (to get x_i of the new (current) i)
        plot_xs.append(x_i)
    plot_ys.append(U/O*x_i + yNplus1 - U/O * x_out) #last y for plotting
    print(f"[{i}, {x_i}, ]")
    print(f"N = {i} stages")
    
    #plot
    op_line_x = np.linspace(0,1,100)
    op_line_y = U/O*op_line_x + yNplus1 - U/O * x_out
    plt.plot(op_line_x, op_line_y, color="blue", label="Operating Line")
    plt.plot(op_line_x, op_line_x, color="orange", label="Equilibrium Line")
    plt.scatter(plot_xs,plot_ys, color="black")
    plt.scatter(x_in, plot_ys[0], color="red", label="(x0,y1)")
    plt.scatter(x_out, yNplus1, color="red", label="(xN, yN+1)")
    j = 0
    while j < i:
        #horizontal line, then vertical line
        plt.plot(plot_xs[j:j+2],[plot_ys[j],plot_ys[j]],color="grey", label="Stepping" if j==0 else "_nolegend_")
        plt.text(sum(plot_xs[j:j+2])/2, sum([plot_ys[j],plot_ys[j]])/2, str(j+1), ha="center", va="bottom") #step number label
        plt.plot([plot_xs[j+1],plot_xs[j+1]],plot_ys[j:j+2],color="grey")
        j += 1
    plt.xlim(0,x_in*1.1)
    plt.ylim(0,plot_ys[0]*1.1)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("McCabe-Thiele Stepping")
    plt.legend(loc="upper left")
    #save graph as png
    img_buffer = io.BytesIO()
    plt.savefig(img_buffer, format="png", bbox_inches="tight")
    plt.close()
    #Return the graph data
    img_buffer.seek(0)
    return img_buffer.getvalue()

def absorption_stripping_single_stage(L_in, V_in, x_in, y_in, P, H, dilute, AbsStrip, linear=False):
    '''
    absorption = solute transfers from gas to liquid
    stripping = solute transfers from liquid to gas
    L - molar flowrate of liquid (typically pure if absorption)
    V - molar flowrate of gas (typically pure if stripping)
    x_in - mole fraction of solute in liquid feed
    y_in - mole fraction of solute in gas feed
    P - operating pressure
    H - Henry's law constant (unit of pressure) for given pair; P and H must have the same units!!!
    dilute - boolean (True or False), True if the liquid solution is dilute, False is concentrated liquid solution
    AbsStrip - boolean (True or False), True if this is for absorption, False for Stripping
    V' - flowrate of carrier gas
    L' - flowrate of pure liquid
    linear - for concentrated solutions, put true if the eq. line given is linear

    For dilute solutions, Henry's Law applies to solute and Raoult's to the solvent, assume L_out=L_in=L and V_out=V_in=V
    For ideal solution, Raoult's Law applies to both components
    can also maybe use Aspen Plus for a problem like this (Jin idea)

    Returns the factor (absorption or stripping), the solute mole fractions of the products streams (x_out, then y_out) and %recovery as a tuple, in that order
    '''
    if AbsStrip == True and dilute == True:
        L = L_in
        V = V_in
        K = H/P  # the smaller K is, the more the solute prefers the liquid; the higher, the more it prefers the gas
        A = (L/V)/K

        # having solved system of equations of op. line and eq. line
        y_out = ((L/V)*x_in + y_in)/(1+A)  # solving for y_out since absorption
        x_out = y_out/K
        percent_recovery = (y_in - y_out)/y_in * 100
        return (A, x_out, y_out, percent_recovery)  # y_out should be low

    elif AbsStrip == True and dilute == False:
        Y_in = y_in/(1-y_in)
        X_in = x_in/(1-x_in)
        K = H/P
        Lp = L_in*(1-x_in)
        Vp = V_in*(1-y_in)

        # solve op.line + eq. line with solveset, save answers in variables X_out and Y_out
        # op. line: Y_out - Y_in = -L'/V'*(X_out-X_in) -> y = -Lp/Vp*(x-X_in) + Y_in
        # eq. line: Y_out = K*X_out/(1+X_out(1-K)) -> y = (K*x)/(1+x*(1-K)); or Y_out=K*X_out if linear, so y=K*x
        # solve by comparison, which gives an equation for solveset (format: op. line - eq. line)
        if linear == False:
            f = -Lp/Vp*(x-X_in)+Y_in - (K*x)/(1+x*(1-K))
        else:
            f = -Lp/Vp*(x-X_in)+Y_in - K*x
        soln_set = sp.solveset(f, x)

        for soln in soln_set:
            # the correct X_out can't be bigger than the x-intercept of the op. line
            soln_set2 = sp.solveset(-Lp/Vp*(t-X_in)+Y_in, t)
            for answer in soln_set2:  # soln_set2 should contain only 1 number, the x-intercept
                x_intercept = answer
            if soln < x_intercept:
                X_out = soln
        if linear == False:
            Y_out = (K*X_out)/(1+X_out*(1-K))
        else:
            Y_out = K*X_out

        # use X_out and Y_out to get x_out, y_out, L_out and V_out and %recovery
        x_out = X_out/(1+X_out)
        y_out = Y_out/(1+Y_out)
        L_out = Lp + Lp*X_out
        V_out = Vp + Vp*Y_out
        percent_recovery = (Y_in-Y_out)/Y_in * 100
        return (x_out, y_out, L_out, V_out, percent_recovery)

    elif AbsStrip == False and dilute == True:
        L = L_in
        V = V_in
        K = H/P
        S = K/(L/V)

        x_out = (x_in + (V/L)*y_in)/(1+S)  # solving for x_out since stripping
        y_out = K*x_out
        percent_recovery = (x_in - x_out)/x_in * 100
        return (S, x_out, y_out, percent_recovery)  # x_out should be low

    elif AbsStrip == False and dilute == False:
        Y_in = y_in/(1-y_in)
        X_in = x_in/(1-x_in)
        K = H/P
        Lp = L_in*(1-x_in)
        Vp = V_in*(1-y_in)

        # solve op.line + eq. line with solveset, save answers in variables X_out and Y_out
        # op. line: Y_out - Y_in = -L'/V'*(X_out-X_in) -> y = -Lp/Vp*(x-X_in) + Y_in
        # eq. line: Y_out = K*X_out/(1+X_out(1-K)) -> y = (K*x)/(1+x*(1-K)); or if linear, Y_out = K*X_out, so y=K*x
        # solve by comparison, which gives an equation for solveset (format: op. line - eq. line)
        if linear == False:
            f = -Lp/Vp*(x-X_in)+Y_in - (K*x)/(1+x*(1-K))
        else:
            f = -Lp/Vp*(x-X_in)+Y_in - K*x
        soln_set = sp.solveset(f, x)

        for soln in soln_set:
            if soln < X_in:
                X_out = soln
        if linear == False:
            Y_out = (K*X_out)/(1+X_out*(1-K))
        else:
            Y_out = K*X_out

        # use X_out and Y_out to get x_out, y_out, L_out and V_out and %recovery
        x_out = X_out/(1+X_out)
        y_out = Y_out/(1+Y_out)
        L_out = Lp + Lp*X_out
        V_out = Vp + Vp*Y_out
        percent_recovery = (X_in-X_out)/X_in * 100
        return (x_out, y_out, L_out, V_out, percent_recovery)


def absorption_stripping_multi_stage(L_in, V_in, x_in, y_in, P, H, dilute, AbsStrip, x_out=None, y_out=None):
    '''
    L - molar flowrate of liquid (typically pure if absorption)
    V - molar flowrate of gas (typically pure if stripping)
    x_in - mole fraction of solute in liquid feed (x0)
    y_in - mole fraction of solute in gas feed (y_N+1)
    P - operating pressure
    H - Henry's law constant (unit of pressure) for given pair; P and H must have the same units!!!
    dilute - boolean (True or False), True if the liquid solution is dilute, False is concentrated liquid solution
    AbsStrip - boolean (True or False), True if this is for absorption, False for Stripping
    V' - flowrate of carrier gas
    L' - flowrate of pure liquid
    x_out - x_N (desired purity - stripping)
    y_out - y1 (desired purity - absorption)

    For dilute solutions, Henry's Law applies to solute and Raoult's to the solvent, assume L_out=L_in=L and V_out=V_in=V
    For ideal solution, Raoult's Law applies to both components
    can also maybe use Aspen Plus for a problem like this (Jin idea)

    Prints a table of mole fractions at each stage, using the algebraic (mccabe-thiele) method to calculate all mole fractions
    And prints required number of stages to achieve desired purity
    Calculates Lmin, L'min, Vmin or V'min using the fact that at minimum, eq. line and op. line intersect, thus there is the same point on both
    lines, thus the equations can now be made equal to each other
    For dilutes, also uses Kremser to confirm answer or give a hint
    '''
    if AbsStrip == True and dilute == True:
        L = L_in
        V = V_in
        K = H/P

        # mccabe-thiele method (aka what was used for washing)
        # operating line: y_n+1 = (L/V)*x_n + y_1 - (L/V)*x_0
        # equilibrium line: y_n = K*x_n
        print("[n, x_n, y_n+1]")
        n = 0
        x_n = x_in
        # start at (x0, y1) on op. line, then keep stepping up until (xN, yN+1) on op. line
        y_nplus1 = y_out
        while y_nplus1 < y_in:
            current_row = [n, x_n, y_nplus1]
            print(current_row)
            n += 1
            # equilibrium line (to get x_n of the new (current) n)
            x_n = y_nplus1/K
            # operating line (to get y_n+1)
            y_nplus1 = L/V*x_n + y_out - L/V*x_in
        current_row = [n, x_n, y_nplus1]
        print(current_row)
        print(f"N = {n} stages")
        
        soln_set = sp.solveset(L*x_in+V*y_in-L*x-V*y_out,x)
        for soln in soln_set:
            x_out = soln
        #Kremser equation: N = ln((y_b-y_b*)/(y_a-y_a*))/ln(A)
        y_a = y_out
        x_a = x_in
        y_as = K*x_a
        y_b = y_in
        x_b = x_out
        y_bs = K*x_b
        A = (L/V)/K
        N = sp.log((y_b-y_bs)/(y_a-y_as))/sp.log(A)
        print(f"Kremser yields {N} stages.")
        
        #realize that for L_min, y_N+1=(L/V)*(x_N-x_0)=K*x_N
        x_out = y_in/K
        #slope_min -> L_min/V
        slope_min = (y_in-y_out)/(x_out-x_in)
        Lp_min = slope_min*V
        print(f"L_min = {Lp_min}")

    elif AbsStrip == True and dilute == False:
        Y_in = y_in/(1-y_in)
        X_in = x_in/(1-x_in)
        Y_out = y_out/(1-y_out)
        K = H/P
        Lp = L_in*(1-x_in)
        Vp = V_in*(1-y_in)

        # operating line: Y_n+1 = (L'/V')*X_n + Y_1 - (L'/V')*X_0
        # equilibrium line: Y_n = (K*X_n)/(1+(1-K)*X_n)
        print("[n, X_n, Y_n+1]")
        n = 0
        X_n = X_in
        # start at (X_0, Y_1) on op. line, then keep stepping up until (X_N, Y_N+1) on op. line
        Y_nplus1 = Y_out
        while Y_nplus1 < Y_in:
            current_row = [n, X_n, Y_nplus1]
            print(current_row)
            n += 1
            # equilibrium line (to get X_n of the new (current) n)
            soln_set = sp.solveset((K*x)/(1+(1-K)*x)-Y_nplus1,x)
            for soln in soln_set:
                X_n = soln
            # operating line (to get Y_n+1)
            Y_nplus1 = (Lp/Vp)*X_n + Y_out - (Lp/Vp)*X_in
        current_row = [n, X_n, Y_nplus1]
        print(current_row)
        print(f"N = {n} stages")
        
        #realize that for L_min, Y_N+1=(K*X_N)/(1+X_N*(1-K))        
        soln_set2 = sp.solveset((K*x)/(1+x*(1-K))-Y_in,x)
        for soln in soln_set2:
            X_out = soln
        #slope_min -> L'min/V'
        slope_min = (Y_in-Y_out)/(X_out-X_in)
        Lp_min = slope_min*Vp
        print(f"L'_min = {Lp_min}")

    elif AbsStrip == False and dilute == True:
        L = L_in
        V = V_in
        K = H/P
        soln_set = sp.solveset(L*x_in+V*y_in-L*x_out-V*y,y)
        for soln in soln_set:
            y_out = soln

        # operating line: y_n+1 = (L/V)*x_n + y_1 - (L/V)*x_0
        # equilibrium line: y_n = K*x_n
        print("[n, x_n, y_n+1]")
        n = 0
        x_n = x_in
        # start at (x0, y1) on op. line, then keep stepping down until (xN, yN+1) on op. line
        y_nplus1 = y_out
        while x_n > x_out:
            current_row = [n, x_n, y_nplus1]
            print(current_row)
            n += 1
            # equilibrium line (to get x_n of the new (current) n)
            x_n = y_nplus1/K
            # operating line (to get y_n+1)
            y_nplus1 = L/V*x_n + y_out - L/V*x_in
        current_row = [n, x_n, y_nplus1]
        print(current_row)
        print(f"N = {n} stages")
        
        #Kremser equation: N = ln((y_b-y_b*)/(y_a-y_a*))/ln(A)
        y_a = y_out
        x_a = x_in
        y_as = K*x_a
        y_b = y_in
        x_b = x_out
        y_bs = K*x_b
        A = (L/V)/K
        N = sp.log((y_b-y_bs)/(y_a-y_as))/sp.log(A)
        print(f"Kremser yields {N} stages.")
        
        #realize that for V_min, y_1=K*x_0       
        y_1 = K*x_in
        slope = (y_1-y_in)/(x_in-x_out) #L/V
        V_min = L/slope
        print(f"V_min = {V_min}")

    elif AbsStrip == False and dilute == False:
        Y_in = y_in/(1-y_in)
        X_in = x_in/(1-x_in)
        X_out = x_out/(1-x_out)
        K = H/P
        Lp = L_in*(1-x_in)
        Vp = V_in*(1-y_in)
        soln_set2 = sp.solveset(Lp*X_in+Vp*Y_in-Lp*X_out-Vp*y,y)
        for soln in soln_set2:
            Y_out = soln

        # operating line: Y_n+1 = (L'/V')*X_n + Y_1 - (L'/V')*X_0
        # equilibrium line: Y_n = (K*X_n)/(1+(1-K)*X_n)
        print("[n, X_n, Y_n+1]")
        n = 0
        X_n = X_in
        # start at (X_0, Y_1) on op. line, then keep stepping down until (X_N, Y_N+1) on op. line
        Y_nplus1 = Y_out
        while X_n > X_out:
            current_row = [n, X_n, Y_nplus1]
            print(current_row)
            n += 1
            # equilibrium line (to get X_n of the new (current) n)
            soln_set = sp.solveset((K*x)/(1+(1-K)*x)-Y_nplus1,x)
            for soln in soln_set:
                X_n = soln
            # operating line (to get Y_n+1)
            Y_nplus1 = (Lp/Vp)*X_n + Y_out - (Lp/Vp)*X_in
        current_row = [n, X_n, Y_nplus1]
        print(current_row)
        print(f"N = {n} stages")
        
        #realize that for V'_min, Y_1=(K*X_0)/(1+X_0*(1-K))        
        Y_1 = (K*X_in)/(1+(1-K)*X_in)
        slope = (Y_1-Y_in)/(X_in-X_out) #L'/V'
        Vp_min = Lp/slope
        print(f"V'_min = {Vp_min}")

def binary_flash_drum_sizing(p_a, p_b, MW_a, MW_b, L, V, x_a, y_a, P, T):
    '''
    Determine the diameter of a flash drum from the inlet and outlet flowrates and their compositions
    a - first component (x_a: liquid mole fraction, y_a: vapour mole fraction)
    b - second component
    L - liquid flowrate
    V - vapour flowrate
    Returns diameter of drum in ft
    
    Densities (p) must be given in g/cm³, molar masses (MW) must be given in g/mol, temperature (T) in °C, pressure (P) in atm
    '''
    x_b = 1 - x_a
    y_b = 1 - y_a
    T += 273.15
    #V, L, compositions, P, T and component molar masses are now known
    
    #outlet molar masses
    MW_l = x_a*MW_a + x_b*MW_b #g/mol
    MW_v = y_a*MW_a + y_b*MW_b #g/mol
    
    #outlet vapour density
    p_v = P*MW_v/(0.082057*T) * 1/1000 #g/L -> g/cm³
    
    #outlet liquid density (use a 1 g basis to start)
    n = 1/MW_l #mol
    V_a = x_a*n*MW_a/p_a
    V_b = x_b*n*MW_b/p_b
    p_l = 1/(V_a+V_b) #g/cm³
    
    Flv = L/V*MW_l/MW_v*np.sqrt(p_v/p_l)
    Y = np.log(Flv)
    
    A=-1.8775
    B=-0.81458
    C=-0.18707
    D=-0.014523
    E=-0.0010149
    K_drum = np.exp(A + B*Y + C*Y**2 + D*Y**3 + E*Y**4) #ft/s
    
    u_lim = K_drum*np.sqrt((p_l-p_v)/p_v) #ft/s
    
    Q = V * 1/3600 * MW_v * 454/1 * 1/p_v * 1/28316.8 #ft³/s
    A = Q/u_lim
    D = np.sqrt(4*A/np.pi)
    return D #ft

def binary_distillation_mccabe_thiele(xD,xB,a,z,full_mode=True,R=None,F=None,q=None,D=None,B=None,L=None,V=None,Lp=None,Vp=None,factor=None):
    '''
    McCabe-Thiele stepping for binary distillation, assuming CMO and feed entering at optimum tray
    Also assumes Total condenser and partial reboiler; steps from the top
    
    When Provided at the minimum xD, xB, alpha and z, the function determines Rmin, Nmin (only in full_mode=True), and N 
    COUNT THE NUMBER OF TRAYS APPROPRIATELY
    full_mode=True is recommended, unless you don't need Nmin or Rmin
    
    Useful equations:
    Top op. line: y_i+1 = L/V*x_i + (1-L/V)*xD = R/(R+1)*x_i + xD/(R+1)
    Bottom op. line: y_i+1 = L'/V'*x_i + (1-L'/V')*xB
    q-line: y = q/(q-1)*x - z/(q-1)
    Equilibrium line (Raoult's): y=a*x/(1+(a-1)*x)
    '''
    #Calculating possibly useful data from given arguments
    if R==None and D!=None and L!=None:
        R = L/D
    if R==None and L!=None and V!=None:
        R = (L/V)/(1-(L/V))
    if D==None and B==None and F!=None:
        #system of equations with total mass balance and component mass balance, x is D and y is B
        soln_set = sp.linsolve([x+y-F,xD*x+xB*y-z*F],(x,y))
        for item in soln_set:
            soln_tup = item
        D = soln_tup[0]
        B = soln_tup[1]
    if R!=None and D!=None and L==None:
        L=R*D
        if V==None:
            V=L+D
    if q==None and Lp!=None and L!=None and F!=None:
        q = (Lp-L)/F
    if q==None and R!=None and Lp!=None and Vp!=None:
        soln_set3 = sp.solveset(R/(R+1)*x+xD/(R+1)-(Lp/Vp*x+(1-Lp/Vp)*xB),x)
        for soln in soln_set3:
            intx2 = soln
        inty2 = R/(R+1)*intx2+xD/(R+1)
        soln_set4 = sp.solveset(s/(s-1)*intx2-z/(s-1)-inty2,s) #s represents q
        for soln in soln_set4:
            q = soln
    if Lp==None and q!=None and F!=None and L!=None:
        Lp=q*F+L
    if Vp==None and B!=None and Lp!=None:
        Vp=Lp-B
    
    if full_mode:
        #From Fenske equation (at total reflux)
        N_min = int(np.ceil(sp.log((xD*(1-xB))/((1-xD)*xB))/sp.log(a)))
        print(f"N_min = {N_min} ({N_min-1} stages + reboiler)")
        
        #min reflux occurs when top op. line intersects the eq. line, so find slope of this theoretical top op line using q-line
        soln_set2 = sp.solveset(a*x/(1+(a-1)*x)-(q/(q-1)*x-z/(q-1)),x) #make q-line and eq. line equal to find intersection x, then y
        for soln in soln_set2:
            if 0 <= soln <= 1:
                intx = soln
        inty = a*intx/(1+(a-1)*intx)
        min_slope = (xD-inty)/(xD-intx)
        R_min = min_slope/(1-min_slope)
        print(f"R_min = {R_min}")
        if R == None and factor != None:
            R = factor*R_min
    
    #real number of stages (stepping starting from the top of the column, with (x0,y1) being (xD,xD))
    print("[i, x_i, y_i+1]")
    i = 0
    x_i = xD
    # start at (X_0, Y_1) on op. line, then keep stepping down until (X_N, Y_N+1) on op. line
    y_iplus1 = xD
    while x_i > xB:
        current_row = [i, x_i, y_iplus1]
        print(current_row)
        i += 1
        # equilibrium line (to get x_i of the new (current) i)
        soln_set = sp.solveset((a*x)/(1+(a-1)*x)-y_iplus1,x)
        for soln in soln_set:
            x_i = soln
        # operating line (to get y_i+1) - step to op. line with the lower value (before intersection, bottom is higher than top, so step on top; after intersection, bottom is lower than top, so step on bottom)
        y_iplus1 = min(R/(R+1)*x_i+xD/(R+1),Lp/Vp*x_i+(1-Lp/Vp)*xB)
    current_row = [i, x_i, y_iplus1]
    print(current_row)
    print(f"N = {i} stages")

def column_diameter(WL, pL, WV, pV, sigma, Q, f, n, spacing):
    '''
    returns the diameter of a distillation column based on provided arguments
    WL - liquid mass flowrate
    WV - vapour mass flowrate
    (WL and WV same units)
    pV - vapour density
    pL - liquid density (should be much higher than pV)
    (pV and pL same units)
    sigma - surface tension (dyne/cm)
    Q - vapour volumetric flowrate (ft³/h)
    f - fraction of flooding velocity
    n - ratio of tray area to cross-sectional area of column
    spacing - tray spacing (in), integer to match key in coefficients dictionary
    
    unlike the flash drum function, please determine pV and pL on your own before inputting
    MAKE SURE UNITS ARE CONSISTENT
    '''
    tray_spacing_coefficients = {6: [1.1977,0.53143,0.18760], 9: [1.1622,0.56014,0.18168], 12: [1.0674,0.55780,0.17919],
                                  18: [1.02262,0.63513,0.20097], 24: [0.94506,0.70234,0.22618], 36: [0.85984,0.73980,0.23725]}
    F_LV = WL/WV*np.sqrt(pV/pL)
    logr = np.log10(F_LV)
    a0, a1, a2 = tray_spacing_coefficients[spacing]
    
    Csb_f = 10**(-a0 - a1*logr - a2*logr**2)
    u_f = Csb_f * (sigma/20)**0.2 * np.sqrt((pL-pV)/pV) #ft/s
    u = f * u_f * 3600 #ft/h
    A = Q/u #ft²
    Ac = A/n
    d = np.sqrt(4*Ac/float(pi))
    return d
    

def FUGK(components, F, q, R=None, factor=None):
    '''
    components: dictionary, where each component is a key (name as string or number as an int), and their data is the value (keep all data about the problem in one spot)
    value should be the list: z_i, fD,i, fB,i, a_i,ref, key (as a string)
    assume class 2 separation (i.e. LNK only at top and HNK only at bottom)
    
    the function uses FUGK to calculate the actual number of stages and optimum feed tray based on provided data
    will take into account number of valid phi values
    works for up to 2 valid phis (2 valid phis with 3 components, one being a DNK, that is)
    
    make F=1 if the value isn't provided
    q is feed quality
    R is either given or is calculated from Rmin*factor
    function returns N (actual number of stages) and Nf (optimum feed tray)
    '''
    if F==None:
        F=1
    
    #1. Fenske
    for key in components:
        if components[key]["key"] == 'LK':
            fD_A = components[key]["fD_i"]
            fB_A = components[key]["fB_i"]
            aA = components[key]["a_i_ref"]
            zLK = components[key]["z_i"]
        if components[key]["key"] == 'HK':
            fD_B = components[key]["fD_i"]
            fB_B = components[key]["fB_i"]
            aB = components[key]["a_i_ref"]
            zHK = components[key]["z_i"]
    aAB = aA/aB #aAB=KA/KB=KA/Kref / KB/Kref = aA/aB
    Nmin = np.log(fD_A*fB_B/((1-fD_A)*(1-fB_B)))/np.log(aAB)
    
    #2. Underwood 1
    valid_phis = []
    LHS = F*(1-q)
    RHS = 0
    for key in components:
        RHS += (components[key]["a_i_ref"]*F*components[key]["z_i"])/(components[key]["a_i_ref"]-p)
    f = RHS - LHS
    soln_set = sp.solveset(f,p)
    for soln in soln_set:
        if aB <= soln <= aA:
            valid_phis.append(soln)
    
    #3. Underwood 2
    Vmin = 0
    if len(valid_phis) == 1:
        D = 0
        for key in components:
            D += components[key]["fD_i"]*components[key]["z_i"]*F
        B = F - D
        
        phi = float(valid_phis[0])
        for key in components:
            Vmin += (components[key]["a_i_ref"]*components[key]["fD_i"]*F*components[key]["z_i"])/(components[key]["a_i_ref"]-phi)
        Rmin = (Vmin-D)/D
    elif len(valid_phis) == 2:
        #x->missing Dxi,D, y->Vmin, f->Vmin with phi1, g->Vmin with phi2
        phi1 = float(valid_phis[0])
        phi2 = float(valid_phis[1])
        f = 0
        g = 0
        D = 0
        for key in components:
            if components[key]["fD_i"] != None and components[key]["fB_i"] != None:
                f += (components[key]["a_i_ref"]*components[key]["fD_i"]*F*components[key]["z_i"])/(components[key]["a_i_ref"]-phi1)
                g += (components[key]["a_i_ref"]*components[key]["fD_i"]*F*components[key]["z_i"])/(components[key]["a_i_ref"]-phi2)
                D += components[key]["fD_i"]*components[key]["z_i"]*F
            else:
                f += (components[key]["a_i_ref"]*x)/(components[key]["a_i_ref"]-phi1)
                g += (components[key]["a_i_ref"]*x)/(components[key]["a_i_ref"]-phi2)
        f -= y
        g -= y
        soln_set = sp.linsolve([f,g],(x,y))
        for soln in soln_set:
            soln_as_list = list(soln)
        Vmin = soln_as_list[1]
        D += soln_as_list[0]
        B = F - D
        Rmin = (Vmin-D)/D
    
    #4. Gilliland
    if factor!=None:
        R = factor*Rmin
    xG = (R-Rmin)/(R+1)
    if 0 <= xG < 0.01:
        RHS = 1 - 18.5715*xG
    elif 0.01 <= xG < 0.9:
        RHS = 0.545827 - 0.591422*xG + 0.002743/xG
    elif 0.9 <= xG < 1:
        RHS = 0.16595-0.16595*xG
    f = (y-Nmin)/(y+1) - RHS
    soln_set = sp.solveset(f,y)
    for soln in soln_set:
        N = soln
    
    #5. Kirkbride
    xLK_B = fB_A*F/B*zLK
    xHK_D = fD_B*F/D*zHK
    f = sp.log(((x-1)/(N-x)),10)-0.206*sp.log((B/D*zHK/zLK*(xLK_B/xHK_D)**2),10)
    soln_set = sp.solveset(f,x)
    for soln in soln_set:
        Nf = soln
    
    return (np.ceil(N),Nf)

# testing (if module is not imported)
if __name__ == '__main__':
    #washing_mccabe_thiele(0.5,1,0.01,0.635)
    #print(absorption_stripping_single_stage(100,25,0,0.015,1,2.5,True,True))
    #absorption_stripping_multi_stage(100,160,0.05,0,1,0.85,True,False,x_out=0.00120671)
    #print(binary_flash_drum_sizing(0.7914,1,32.04,18.01,736,264,0.2,0.579,1,81.7))
    #binary_distillation_mccabe_thiele(0.96,0.02,1.76,0.1,full_mode=False,R=4)
    #print(column_diameter(17500,700,19000,3.5,20,19000/3.5*35.3147,0.8,0.9,24))
    print(FUGK({'P': [0.25,0.995,0.005,13.33,'LK'],'HX': [0.25,0.05,0.95,5.42,'HK'],'HP': [0.25,0,1,2.33,'HNK'],'O': [0.25,0,1,1,'HNK']}, None, 0, R=6)) #1 valid phi (straightforward)
    #print(FUGK({'B': [0.397,0.9992,0.0008,2.25,'LK'],'T': [0.167,None,None,1,'DNK'], 'C': [0.436,0.0001,0.9999,0.21,'HK']}, 1000, 1, R=1.2)) #2 valid phis (system of equations)
