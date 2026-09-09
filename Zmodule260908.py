# this module contains function I think could be one day useful, and other random stuff like tips and constants
import sympy as sp
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp
import pandas as pd

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

# transposes a value within a range to its corresponding value in another range


def mapping(value, fromMin, fromMax, toMin, toMax):
    fromRange = fromMax - fromMin
    toRange = toMax - toMin
    valueScaled = float(value - fromMin) / float(fromRange)
    return toMin + (valueScaled * toRange)

# not amazing, but reverses a list using recursion (worth saving cz it's cool and i'm proud of myself for figuring this out without using the hint! :)
def reverse(seq):
    if len(seq) == 0 or len(seq) == 1:
        return seq

    if len(seq) == 2:
        lily = []
        lily.append(seq[1])
        lily.append(seq[0])
        return lily

    newseq = []
    seq2 = seq[0:]
    i = 0
    j = 1
    while i < int(np.ceil(len(seq)/2)):
        if i == 0:
            seq3 = []
            seq3.append(seq2[0])
            seq3.append(seq2[len(seq2)-1])
            seq3 = reverse(seq3)
            newseq.insert(0, seq3[0])
            newseq.insert(1, seq3[1])
            seq2 = seq2[1:len(seq2)-1]
            i += 1
        if i > 0:
            for number in reverse(seq2):
                newseq.insert(j, number)
                j += 1
            i += 1
    return newseq

# from excel, inspired by doing WEO budgets. package is number per package; required is number of item required
def determine_required_packages(package, required):
    if required // package == 0:
        return required // package
    else:
        return required // package + 1

# temperature conversion function between °C, F, K
def temp_convert(temperature_now, new_temperature, value):
    if temperature_now in ['C', 'F', 'K'] and new_temperature in ['C', 'F', 'K'] and (type(value) == int or type(value) == float):
        if temperature_now == new_temperature:
            return value

        if temperature_now == 'C':
            if new_temperature == 'F':
                return 9/5*value + 32
            else:  # if 'K'
                return value + 273.15
        elif temperature_now == 'F':
            if new_temperature == 'C':
                return (value-32)*5/9
            else:  # if 'K'
                return (value-32)*5/9 + 273.15
        else:  # if 'K'
            if new_temperature == 'C':
                return value - 273.15
            else:  # if 'F'
                return (value - 273.15)*9/5 + 32
    else:
        return None

# taken from CHE120 jupyter notebook notes
def quadratic_root(a, b, c):
    """(float, float, float) -> complex

    Return the maximum root of a quadratic function $f(x) = a x^2 + b x + c$ with real coefficients.
    The maximum root is determined from the absolute value or modulus of the values.


        >>> quadratic_root(2., 3., 1.)
        -1.0
        >>> quadratic_root(1., -2., -3.)
        3.0
        >>> quadratic_root(1., 0., 1.)
        -1j
    """
    # using the quadratic formula, compute the two roots of the
    # polynomial
    root1 = (-b + (b**2 - 4. * a * c)**0.5)/(2. * a)
    root2 = (-b - (b**2 - 4. * a * c)**0.5)/(2. * a)

    # use the built-in function `abs()` to determine which is maximum, note:
    # - `abs()` that this function computes the modulus of a complex number
    # - this works for repeated roots, does not matter which is returned!
    if abs(root1) > abs(root2):
        return (root1, root2)
    else:
        return (root1, root2)

# taken from CHE120 jupyter notebook notes, apparently the fastest of 3 methods
def find_two_smallest(L):
    """ (list of float) -> tuple of (int, int)

    Return a tuple of the indices of the two smallest values in list L.

    >>> find_two_smallest([809, 834, 477, 478, 307, 122, 96, 102, 324, 476])
    (6, 7)
    """

    # Find the index of the minimum and remove that item.
    smallest = min(L)
    min1 = L.index(smallest)
    L.remove(smallest)

    # Find the index of the new minimum.
    next_smallest = min(L)
    min2 = L.index(next_smallest)

    # Put smallest back into L.
    L.insert(min1, smallest)

    # Fix min2 in case it was affected by the re-insertion.
    if min1 <= min2:
        min2 += 1

    return (min1, min2)

# sympy documentation function to make continued fraction
def list_to_frac(l):
    expr = sp.Integer(0)
    for i in reversed(l[1:]):
        expr += i
        expr = 1/expr
    return l[0] + expr


def trapezoid_rule_calculator(f, a, b, n):
    total = 0
    delta_x = (b-a)/n
    list_of_values = []

    for j in range(n+1):
        list_of_values.append(a + j*delta_x)

    for i in list_of_values:
        if i == a or i == b:
            total += f.subs(x, i)/2
        else:
            total += f.subs(x, i)

    f2 = abs(sp.diff(f, x, 2))
    x_values = np.linspace(a, b, 1000)
    y_values = [f2.subs(x, k) for k in x_values]
    M = max(y_values)
    error = abs(M*(b-a)**3/(12*n**2))  # |Tn|

    return (total*delta_x, error)


def simpson_rule_calculator(f, a, b, n, datapoints=False):
    total = 0
    delta_x = (b-a)/n
    list_of_values = []
    for j in range(n+1):
        list_of_values.append(a + j*delta_x)
        
    if not datapoints:
        i = 0
        while i < len(list_of_values):
            if i == 0 or i == len(list_of_values) - 1:
                total += f.subs(x, list_of_values[i])
            elif i != 0 and i != len(list_of_values) - 1 and i % 2 == 0:
                total += 2*f.subs(x, list_of_values[i])
            elif i != 0 and i != len(list_of_values) - 1 and i % 2 == 1:
                total += 4*f.subs(x, list_of_values[i])
            i += 1
    
        f4 = abs(sp.diff(f, x, 4))
        x_values = np.linspace(a, b, 1000)
        y_values = [f4.subs(x, k) for k in x_values]
        M = max(y_values)
        error = abs(M*(b-a)**5/(180*n**4))  # |Sn|
    
        return (total*delta_x/3, error)
    else: #modification to use simpson's based on lab data; but moll already made a program for that
        i = 0
        while i < len(list_of_values):
            pass


def delta_enthalpy_with_shomates(A, B, C, D, E, T1, T2):
    '''only sensible heat for these delta enthalpy hat functions'''
    t1 = (T1+273.15)/1000
    t2 = (T2+273.15)/1000
    shomates = A + B*t + C*t**2 + D*t**3 + E/t**2
    shomates_integrated = sp.integrate(shomates, t)
    CpdT = shomates_integrated.subs(t, t2) - shomates_integrated.subs(t, t1)
    return CpdT  # kJ/mol


def delta_enthalpy_mix(x1, CpdT1, x2, CpdT2):
    return x1*CpdT1 + x2*CpdT2  # kJ/mol

def time_duration(hours: int, minutes: int) -> None:
    """
    Based on the number of hours and minutes provided,
    Returns a float of the number of hours, minutes and seconds that have passed
    Float of number of hours is rounded to 2 decimals
    Ex:
        >>>time_duration(5,30)
        5.5 hours
        330 minutes
        19800 seconds
        
    Constraints: hours and minutes must be nonnegative integers; minutes in particular must be between 0 and 59
    """
    if type(hours) == int and type(minutes) == int and hours >= 0 and 0 <= minutes <= 59:
        num_hours = round(hours + minutes/60, 2)
        num_minutes = hours*60 + minutes
        num_seconds = num_minutes * 60
        print(f"{num_hours} hours")
        print(f"{num_minutes} minutes")
        print(f"{num_seconds} seconds")
    else:
        print("Please review and reenter your inputs.")
        hours = int(input("Hours: "))
        minutes = int(input("Minutes: "))
        time_duration(hours, minutes)
    

def number_triangle():  # from Melissa Vaziri's COMP 1405 Tutorial 3
    while True:
        # Asking The User For An Integer
        value = input("Enter An Integer Between 1 and 9: ")

        # Checking If The Input Is A Digit and if value is between 1 and 9 inclusive
        # until these conditions are met, we keep asking for a valid integer
        if value.isdigit() and 1 <= int(value) <= 9:
            # Basically, I Assigned "I" As A Variable For The Integer Of The Variable Value
            i = int(value)

            # This Is Used For A Loop That Repeats A Specific Amount Of Times
            for row in range(1, i+1):

                # This Is Saying That Col Repeats, And Range Produces Order of Numbers
                for col in range(row):
                    print(str(row), end="")

                # Printing Statement so that after a number row is printed, we go to a new line for the next row
                print()
            break  # postcondition loop using break. we break out of while True after the triangle is successfully printed; we don't wanna stay stuck in the while True forever like the OG code...


def rref_matrix(A, b):
    # A is the coefficient matrix, b is the vector, quick function for solving systems of equations
    # OR USE sp.linsolve([list of equations in the form = 0], (tuple of variable(s) to solve for))
    # sp.linsolve for a system of equations, sp.solveset for a single equation
    A_b = A.row_join(b)
    return A_b.rref()


def antoines_vapour_pressure(A, B, C, T):
    # T in Celsius originally, must turn to Kelvin first
    T += 273.15
    vapour_pressure_at_T = 10**(A-B/(T+C))
    return vapour_pressure_at_T  # in bar if NIST constants


def washing_mccabe_thiele(U, O, x_out, x_in, yNplus1=0):
    '''
    U - underflow liquid's mass flow rate
    O - overflow liquid's mass flow rate
    x_out - xN, desired purity
    x_in - x0, salt mass fraction of source entrained liquid (underflow)
    y_N+1 - y_in, salt mass fraction of source wash liquid (overflow), typically 0

    Prints a table of mass fractions at each stage, using the algebraic method to calculate all mass fractions
    And prints required number of stages to achieve desired purity
    Please use floats for inputs
    '''
    print("[i, x_i, y_i+1]")
    x_i = x_in  # x0
    i = 0
    while x_i > x_out: #stepping down
        y_iplus1 = U/O*x_i + yNplus1 - U/O * x_out  # operating line (to get y_i+1)
        current_row = [i, x_i, y_iplus1]
        print(current_row)
        i += 1
        x_i = y_iplus1  # equilibrium line (to get x_i of the new (current) i)
    print(f"[{i}, {x_i}, ]")
    print(f"N = {i} stages")


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
    linear - for concentrated solutions, in case the eq. line given is linear, give K as the coefficient of x

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
        print(f"L'_min = {Lp_min}")

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
    
    make F=1 if a value isn't provided (if you won't provide a value or any value is known =None)
    q is feed quality
    R is either given or is calculated from Rmin*factor
    function returns N (actual number of stages) and Nf (optimum feed tray)
    '''
    if F==None:
        F=1
    
    #1. Fenske
    for key in components:
        if components[key][4] == 'LK':
            fD_A = components[key][1]
            fB_A = components[key][2]
            aA = components[key][3]
            zLK = components[key][0]
        if components[key][4] == 'HK':
            fD_B = components[key][1]
            fB_B = components[key][2]
            aB = components[key][3]
            zHK = components[key][0]
    aAB = aA/aB #aAB=KA/KB=KA/Kref / KB/Kref = aA/aB
    Nmin = np.log(fD_A*fB_B/((1-fD_A)*(1-fB_B)))/np.log(aAB)
    
    #2. Underwood 1
    valid_phis = []
    LHS = F*(1-q)
    RHS = 0
    for key in components:
        RHS += (components[key][3]*F*components[key][0])/(components[key][3]-p)
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
            D += components[key][1]*components[key][0]*F
        B = F - D
        
        phi = float(valid_phis[0])
        for key in components:
            Vmin += (components[key][3]*components[key][1]*F*components[key][0])/(components[key][3]-phi)
        Rmin = (Vmin-D)/D
    elif len(valid_phis) == 2:
        #x->missing Dxi,D, y->Vmin, f->Vmin with phi1, g->Vmin with phi2
        phi1 = float(valid_phis[0])
        phi2 = float(valid_phis[1])
        f = 0
        g = 0
        D = 0
        for key in components:
            if components[key][1] != None and components[key][2] != None:
                f += (components[key][3]*components[key][1]*F*components[key][0])/(components[key][3]-phi1)
                g += (components[key][3]*components[key][1]*F*components[key][0])/(components[key][3]-phi2)
                D += components[key][1]*components[key][0]*F
            else:
                f += (components[key][3]*x)/(components[key][3]-phi1)
                g += (components[key][3]*x)/(components[key][3]-phi2)
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
    if R == None and factor!=None:
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

def series_convergence():
    pass

def two_point_interpolation(x1,y1,x2,y2,x_value):
    '''
    2-point linear interpolation to determine an unknown value from given values
    '''
    result = y1 + (y2-y1)/(x2-x1)*(x_value-x1)
    return result

def hypothesis_testing(alpha, data_pop1, var_pop1, data_pop2=None, var_pop2=None):
    '''
    CHE220 Hypothesis testing function
    data_pop1 and data_pop2, if provided, are lists
    for a paired t-test, please provided measurements for the same object/person as the same index in each list of data
    for a proportion test, enter data_pop arguments as lists with 2 entries: index 0 is x, index 1 is n (e.g. [x,n])
    use None as input if information is unknown (i.e. for required arguments like var_pop1 but variance is unknown)
    the function will figure out which test to use from the given inputs, assuming the population and sample are normally distributed
    Will print whether or not to reject the null hypothesis, as well as the p-value
    '''
    if data_pop2 != None:
        num_pop = input("To confirm, do you wish to do 2 population hypothesis testing? (yes, no) ")
    if num_pop.lower() == 'no':
        num_pop = 1
    else:
        num_pop = 2
    #see what they want to test
    subject = input("What do you wish to test for? (mean, variance, proportion) ")
    num_sides = input("2-sided, lower or upper? ")
    h0 = float(input("What hypothesis value to test? (please enter a number) "))
    #usual data about sample(s)
    x_bar1 = np.mean(data_pop1)
    s1 = np.std(data_pop1, ddof=1)
    n1 = len(data_pop1)
    if data_pop2 != None and num_pop == 2:
        x_bar2 = np.mean(data_pop2)
        n2 = len(data_pop2)
        s2 = np.std(data_pop2, ddof=1)
    
    #tests for mean
    if subject.lower() == 'mean':
        if num_pop == 1:
            if var_pop1 != None:
                #1 pop, variance known
                z0 = (x_bar1-h0)/np.sqrt(var_pop1/n1)
                if num_sides == '2-sided':
                    if abs(z0) >= scp.stats.norm.ppf(1-alpha/2,loc=0,scale=1):
                        reject = True
                    else:
                        reject = False
                    p_value = 2*(1-scp.stats.norm.cdf(abs(z0)))
                elif num_sides == 'upper':
                    if z0 >= scp.stats.norm.ppf(1-alpha,loc=0,scale=1):
                        reject = True
                    else:
                        reject = False
                    p_value = 1-scp.stats.norm.cdf(z0)
                elif num_sides == 'lower':
                    if z0 <= -scp.stats.norm.ppf(1-alpha,loc=0,scale=1):
                        reject = True
                    else:
                        reject = False
                    p_value = scp.stats.norm.cdf(z0)
                return (reject, p_value)
            else:
                #1 pop, variance unknown
                t0 = (x_bar1-h0)/(s1/np.sqrt(n1))
                if num_sides == '2-sided':
                    if abs(t0) >= scp.stats.t.ppf(1-alpha/2,df=n1-1):
                        reject = True
                    else:
                        reject = False
                    p_value = 2*scp.stats.t.sf(abs(t0),df=n1-1)
                elif num_sides == 'upper':
                    if t0 >= scp.stats.t.ppf(1-alpha,df=n1-1):
                        reject = True
                    else:
                        reject = False
                    p_value = scp.stats.t.sf(t0,df=n1-1)
                elif num_sides == 'lower':
                    if t0 <= scp.stats.t.ppf(alpha,df=n1-1):
                        reject = True
                    else:
                        reject = False
                    p_value = 1-scp.stats.t.sf(t0,df=n1-1)
                return (reject, p_value)
        
        elif num_pop == 2:
            #2 pop, variances known
            if var_pop1 != None and var_pop2 != None:
                z0 = (x_bar1-x_bar2-h0)/np.sqrt(var_pop1/n1 + var_pop2/n2)
                if num_sides == '2-sided':
                    if abs(z0) >= scp.stats.norm.ppf(1-alpha/2,loc=0,scale=1):
                        reject = True
                    else:
                        reject = False
                    p_value = 2*(1-scp.stats.norm.cdf(abs(z0)))
                elif num_sides == 'upper':
                    if z0 >= scp.stats.norm.ppf(1-alpha,loc=0,scale=1):
                        reject = True
                    else:
                        reject = False
                    p_value = 1-scp.stats.norm.cdf(z0)
                elif num_sides == 'lower':
                    if z0 <= -scp.stats.norm.ppf(1-alpha,loc=0,scale=1):
                        reject = True
                    else:
                        reject = False
                    p_value = scp.stats.norm.cdf(z0)
                return (reject, p_value)
            else:
                variance_equality = input("Are the variances equal? (yes, no) ")
                if variance_equality == 'yes':
                    paired_data = input("Are the data points pairs across each list (aka is this for a paired t-test?) (yes, no) ")
                    if paired_data == 'yes':
                        #paired t-test
                        differences = [data_pop2[i]-data_pop1[i] for i in range(len(data_pop1))]
                        d_bar = np.mean(differences)
                        sD = np.std(differences, ddof=1)
                        n = len(differences)
                        t0 = (d_bar-h0)/(sD/np.sqrt(n))
                        if num_sides == '2-sided':
                            if abs(t0) >= scp.stats.t.ppf(1-alpha/2,df=n1-1):
                                reject = True
                            else:
                                reject = False
                            p_value = 2*scp.stats.t.sf(abs(t0),df=n1-1)
                        elif num_sides == 'upper':
                            if t0 >= scp.stats.t.ppf(1-alpha,df=n1-1):
                                reject = True
                            else:
                                reject = False
                            p_value = scp.stats.t.sf(t0,df=n1-1)
                        elif num_sides == 'lower':
                            if t0 <= scp.stats.t.ppf(alpha,df=n1-1):
                                reject = True
                            else:
                                reject = False
                            p_value = 1-scp.stats.t.sf(t0,df=n1-1)
                        return (reject, p_value)
                    else:
                        #2 pop, variances unknown but equal
                        sp2 = ((n1-1)*s1**2+(n2-1)*s2**2)/(n1+n2-2)
                        t0 = (x_bar1-x_bar2-h0)/(np.sqrt(sp2)*np.sqrt(1/n1+1/n2))
                        if num_sides == '2-sided':
                            if abs(t0) >= scp.stats.t.ppf(1-alpha/2,df=n1+n2-2):
                                reject = True
                            else:
                                reject = False
                            p_value = 2*scp.stats.t.sf(abs(t0),df=n1+n2-2)
                        elif num_sides == 'upper':
                            if t0 >= scp.stats.t.ppf(1-alpha,df=n1+n2-2):
                                reject = True
                            else:
                                reject = False
                            p_value = scp.stats.t.sf(t0,df=n1+n2-2)
                        elif num_sides == 'lower':
                            if t0 <= scp.stats.t.ppf(alpha,df=n1+n2-2):
                                reject = True
                            else:
                                reject = False
                            p_value = 1-scp.stats.t.sf(t0,df=n1+n2-2)
                        return (reject, p_value)
                else:
                    #2 pop, variances unknown and unequal
                    v = (s1**2/n1+s2**2/n2)**2/((s1**2/n1)**2/(n1-1)+(s2**2/n2)**2/(n2-1))
                    t0 = (x_bar1-x_bar2-h0)/np.sqrt(s1**2/n1+s2**2/n2)
                    if num_sides == '2-sided':
                        if abs(t0) >= scp.stats.t.ppf(1-alpha/2,df=v):
                            reject = True
                        else:
                            reject = False
                        p_value = 2*scp.stats.t.sf(abs(t0),df=v)
                    elif num_sides == 'upper':
                        if t0 >= scp.stats.t.ppf(1-alpha,df=v):
                            reject = True
                        else:
                            reject = False
                        p_value = scp.stats.t.sf(t0,df=v)
                    elif num_sides == 'lower':
                        if t0 <= scp.stats.t.ppf(alpha,df=v):
                            reject = True
                        else:
                            reject = False
                        p_value = 1-scp.stats.t.sf(t0,df=v)
                    return (reject, p_value)
        
    #tests for variance
    if subject.lower() == 'variance':
        if num_pop == 1:
            #1 pop, variance
            chi0 = (n1-1)*s1**2/h0
            if num_sides == '2-sided':
                if chi0 <= scp.stats.chi2.ppf(alpha/2,n1-1) or chi0 >= scp.stats.chi2.ppf(1-alpha/2,n1-1):
                    reject = True
                else:
                    reject = False
                p_value = 2*scp.stats.chi2.sf(chi0, n1-1)
            elif num_sides == 'upper':
                if chi0 >= scp.stats.chi2.ppf(1-alpha,n1-1):
                    reject = True
                else:
                    reject = False
                p_value = scp.stats.chi2.sf(chi0, n1-1)
            elif num_sides == 'lower':
                if chi0 <= scp.stats.chi2.ppf(alpha,n1-1):
                    reject = True
                else:
                    reject = False
                p_value = 1-scp.stats.chi2.sf(chi0, n1-1)
            return (reject, p_value)
    
        elif num_pop == 2:
            #2 pop, (ratio of) variances
            f0 = (s1**2/s2**2)*(1/h0)
            if num_sides == '2-sided':
                if f0 <= scp.stats.f.ppf(alpha/2,n1-1,n2-1) or chi0 >= scp.stats.chi2.ppf(1-alpha/2,n1-1,n2-1):
                    reject = True
                else:
                    reject = False
                p_value = 2*scp.stats.f.sf(f0, n1-1, n2-1)
            elif num_sides == 'upper':
                if f0 >= scp.stats.f.ppf(1-alpha,n1-1,n2-1):
                    reject = True
                else:
                    reject = False
                p_value = scp.stats.f.sf(f0, n1-1, n2-1)
            elif num_sides == 'lower':
                if f0 <= scp.stats.f.ppf(alpha,n1-1,n2-1):
                    reject = True
                else:
                    reject = False
                p_value = 1-scp.stats.f.sf(f0, n1-1, n2-1)
            return (reject, p_value)
    
    #tests for proportion
    if subject.lower() == 'proportion':
        if num_pop == 1:
            #1 pop, proportion
            x, n = data_pop1[0], data_pop1[1]
            z0 = (x/n-h0)/np.sqrt(h0*(1-h0)/n)
            if num_sides == '2-sided':
                if abs(z0) >= scp.stats.norm.ppf(1-alpha/2,loc=0,scale=1):
                    reject = True
                else:
                    reject = False
                p_value = 2*(1-scp.stats.norm.cdf(abs(z0)))
            elif num_sides == 'upper':
                if z0 >= scp.stats.norm.ppf(1-alpha,loc=0,scale=1):
                    reject = True
                else:
                    reject = False
                p_value = 1-scp.stats.norm.cdf(z0)
            elif num_sides == 'lower':
                if z0 <= -scp.stats.norm.ppf(1-alpha,loc=0,scale=1):
                    reject = True
                else:
                    reject = False
                p_value = scp.stats.norm.cdf(z0)
            return (reject, p_value)
        elif num_pop == 2:
            #2 pop, proportion
            x1, n1, x2, n2 = data_pop1[0], data_pop1[1], data_pop2[0], data_pop2[1]
            p_hat = (x1+x2)/(n1+n2)
            p1_hat = x1/n1
            p2_hat = x2/n2
            z0 = (p1_hat-p2_hat-h0)/np.sqrt(p_hat*(1-p_hat)*(1/n1+1/n2))
            if num_sides == '2-sided':
                if abs(z0) >= scp.stats.norm.ppf(1-alpha/2,loc=0,scale=1):
                    reject = True
                else:
                    reject = False
                p_value = 2*(1-scp.stats.norm.cdf(abs(z0)))
            elif num_sides == 'upper':
                if z0 >= scp.stats.norm.ppf(1-alpha,loc=0,scale=1):
                    reject = True
                else:
                    reject = False
                p_value = 1-scp.stats.norm.cdf(z0)
            elif num_sides == 'lower':
                if z0 <= -scp.stats.norm.ppf(1-alpha,loc=0,scale=1):
                    reject = True
                else:
                    reject = False
                p_value = scp.stats.norm.cdf(z0)
            return (reject, p_value)
        
def least_squares(alpha, k, y_values, x1_values, *args):
    '''
    CHE220 Least Squares calculations and ANOVA for simple and multiple linear regression
    k: number of regressor variables (PLEASE SPECIY AS INT)
    x_values and y_values are lists of the data, should have same lengths
    first data point for each list should correspond, i.e. represent (x11,x12,...,y1)
    provided as many lists of x_values as the number of regressor variables
    the first x_values (required argument x1_values) is [x11,x21,...,xn1]
    if you wish to do polynomial fitting, the extra args should be the appropriate lists (they will not be combined within the function)
    '''
    #finding equation
    n = len(y_values)
    p = k+1
    
    y_array = np.array(y_values).reshape(n,1) #shape nx1
    X_array = np.ones((n,1))
    x_array = np.array(x1_values).reshape(n,1)
    X_array = np.concatenate((X_array,x_array),axis=1) #will eventually be shape nxk
    for arg in args:
        x_array = np.array(arg).reshape(n,1)
        X_array = np.concatenate((X_array,x_array),axis=1)
    B = np.linalg.inv(X_array.T @ X_array) @ X_array.T @ y
    B = B.reshape(1,k)
    print(B)
    
    #ANOVA
    i = 0
    y_hat_list = []
    while i < n:
        y_hat = X_array[i] @ B
        y_hat_list.append(y_hat)
        i += 1
    y_hat_array = np.array(y_hat_list).reshape(n,1)
    residuals_array = y_array - y_hat_array
    residuals_array **= 2
    SSE = residuals_array.sum(axis=0)
    sigma_hat_squared = SSE/(n-p)
    SST = y_array.T@y_array - sum(y_values)**2/n
    SSR = SST - SSE
    MSR = SSR/k
    MSE = sigma_hat_squared
    f0 = MSR/MSE
    p_value = scp.stats.f.sf(f0,k,n-p)
    if p_value <= alpha:
        print("reject the null hypothesis. there is strong evidence of a linear relationship")
    else:
        print("fail to reject the null hypothesis. there is no strong evidence of a linear relationship")
    R_squared = SSR/SST
    adjusted_R_squared = 1 - (SSE/(n-p))/(SST/(n-1))
    print(p_value, R_squared, adjusted_R_squared)
    
"""
def convergence_divergence(series):
    '''
    convergence/divergence of series from MATH118
    the function will ask which tests to use to determine if the series if conditionally converges, absolutely converges or diverges
    test options are: geometric, integral, p-series, comparison, LCT, AST, ratio, root (divergence is automatically done by the function as the first test and as an excuse to calculate the limit)
    please use y to write the series sequence (not n)
    '''
    #start with divergence test
    print("Hello there! Let's start with Divergence Test.")
    lim = sp.limit(series,y,inf)
    if lim != 0 or lim == inf or lim == -inf:
        return 'diverges'
    else:
        print("Divergence test inconclusive (meaning the limit is 0).")
        conclusion_reached = False
        while not conclusion_reached:
            #try different tests until one reaches a conclusion
            next_test = input("Which test to try next? ")
            if next_test == 'geometric':
                pass
                if abs(r) >= 1:
                    return 'diverges'
                    conclusion_reached = True
                else:
                    return 'absolutely converges'
                    conclusion_reached = True
            elif next_test == 'integral':
                #ENSURE SERIES SEQUENCE IS CONTINOUS, POSITIVE AND DECREASING ON [1,INF)
                lily = [series.subs(y,i) for i in range(50)] #if this runs without any problems, probably continuous
                h = 0
                ouch =  False
                while h < len(lily)-1:
                    if lily[h] <= lily[h+1] or lily[h] <= 0 or lily[h+1] <= 0:
                        ouch = True
                        break
                if ouch == True:
                    print("Function not decreasing or not positive")
                else:
                    result = sp.integrate(series,(y,1,inf))
                    if result == inf or result == -inf:
                        return 'diverges'
                        conclusion_reached = True
                    else:
                        return 'converges, run a different test with abs(series) to confirm absolute or conditional convergence'
                        conclusion_reached = True
            elif next_test == 'p-series':
                pass
                if p > 1:
                    return 'absolutely converges'
                    conclusion_reached = True
                else:
                    return 'diverges'
                    conclusion_reached = True
            elif next_test == 'comparison':
                pass
            elif next_test == 'LCT':
                pass
            elif next_test == 'AST':
                pass
                lim = sp.limit(f,y,inf)
                lily = [f.subs(y,i) for i in range(50)]
                h = 0
                ouch =  False
                while h < len(lily)-1:
                    if lily[h] <= lily[h+1]:
                        ouch = True
                        break
                if ouch == True:
                    print("Function not decreasing")
                if lim == 0 and ouch == False:
                    return 'converges, run a different test with abs(series) to confirm absolute or conditional convergence'
                    conclusion_reached = True
            elif next_test == 'ratio':
                pass
            elif next_test == 'root':
                pass
"""

# testing (if module is not imported)
if __name__ == '__main__':
    #print(round(mapping(-42,-95,-42,0,9),2))
    #print(reverse([1,2,3]))
    #print(determine_required_packages(20,30))
    #print(temp_convert('F','C',32))
    #print(quadratic_root(2.,3.,1.))
    #print(find_two_smallest([1,4,5,62,4,6,7,2,8,9]))
    #print(list_to_frac([x,y,z]))
    #print(trapezoid_rule_calculator(1/x**3,1,4,6))
    #print(simpson_rule_calculator(1/x**3,1,4,6))
    #print(delta_enthalpy_with_shomates(20.17,0.4001,0,0,0,25,100))
    #print(delta_enthalpy_mix(0.9989,delta_enthalpy_with_shomates(-203.6060,1523.290,-3196.413,2474.455,3.855326,12.5,134.7),0.011,delta_enthalpy_with_shomates(50.72389,6.72267,-2.517167,10.15934,-0.200675,12.5,134.7))) #seawater module 3-4 assignment
    #number_triangle()
    #print(rref_matrix(sp.Matrix([[1,1],[2,1],[1,-1]]),sp.Matrix([[-1],[2],[7]])))
    #print(antoines_vapour_pressure(5.20409,1581.341,-33.5,37))
    #washing_mccabe_thiele(0.5,1,0.01,0.635)
    #print(absorption_stripping_single_stage(100,25,0,0.015,1,2.5,True,True))
    #absorption_stripping_multi_stage(100,160,0.05,0,1,0.85,True,False,x_out=0.00120671)
    #print(binary_flash_drum_sizing(0.7914,1,32.04,18.01,736,264,0.2,0.579,1,81.7))
    #binary_distillation_mccabe_thiele(0.96,0.02,1.76,0.1,full_mode=False,R=4)
    #print(column_diameter(17500,700,19000,3.5,20,19000/3.5*35.3147,0.8,0.9,24))
    #print(FUGK({'P': [0.25,0.995,0.005,13.33,'LK'],'HX': [0.25,0.05,0.95,5.42,'HK'],'HP': [0.25,0,1,2.33,'HNK'],'O': [0.25,0,1,1,'HNK']}, None, 0, R=6)) #1 valid phi (straightforward)
    #print(FUGK({'B': [0.397,0.9992,0.0008,2.25,'LK'],'T': [0.167,None,None,1,'DNK'], 'C': [0.436,0.0001,0.9999,0.21,'HK']}, 1000, 1, R=1.2)) #2 valid phis (system of equations)
    #print(convergence_divergence((y-1)/(3*y-1)))
    time_duration(5,-1)
