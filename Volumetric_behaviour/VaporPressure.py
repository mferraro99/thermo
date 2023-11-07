# Code for the automatic estimation of vapor pressure of pure components
# Single component, real gas/liquid
# Given the conditions of the system (T or P) and the properties of the
# compound (critical T, critical P, normal boiling T @ 1 bar
# and the Pitzer acentric factor omega),
# the codes estimates the value of P0 @ Tsystem or T0 @ Psystem
# Author: Marcello Ferraro
# Last update: 07/11/2023

import sys
import math
from scipy.optimize import fsolve

# Define constants
R = 8.31446261815324 # [J/mol/K]
pi = 3.14159265358979 # [-]

# Define the model used for the calculation:
# vdW = van der Waals
# RK = Redlich-Kwong
# RKS = Redlich-Kwong-Soave
# PR = Peng-Robinson

model = "PR"

### Define data

# System properties --> "unknown" is the quantity we want to calculate
P = 2.115e5 # System pressure [Pa]
T = "unknown" # System temperature [K]

# Compound properties
Pc = 24.97*1e5 # Critical pressure [Pa]
Tc = 569.40 # Critical temperature [K]
omega = 0.398 # Pitzer acentric factor [-]
TebN = 398.9 # Normal boiling temperature [K]
PebN = 1e5 # Normal pressure (change only if != 1 bar) [Pa]

# The Poynting correction is neglected in a first approximation

if model != "vdW" and model != "RK" and model != "RKS" and \
    model != "PR" and model != "Virial":
    print("ERROR: define a valid model!")
    sys.exit()

### Cubic root finder
def Zroot_calculation(a,b,c):
    p = b-1/3*a**2
    q = 2/27*a**3 - 1/3*a*b + c
    D = 1/4*q**2+1/27*p**3

    if D > 0:
        Z = (-q/2+D**(1/2))**(1/3) + (-q/2-D**(1/2))**(1/3) - a/3
        return [Z,0,0]
    
    elif D == 0:
        Z1 = -2*(q/2)**(1/3) - a/3
        Z2 = (q/2)**(1/3) - a/3
        if Z1 >= Z2: return [Z1,Z2,Z2]
        else: return [Z2,Z2,Z1]

    else:
        r = (-1/27 * p**3)**(1/2)
        teta = math.acos(-q/2/r)

        Z1 = 2 * r**(1/3)*math.cos(teta/3) - a/3
        Z2 = 2 * r**(1/3)*math.cos((teta+2*pi)/3) - a/3
        Z3 = 2 * r**(1/3)*math.cos((teta+4*pi)/3) - a/3
        Z = [Z1,Z2,Z3]
        for i in range(len(Z)):
            # Flag to optimize the sorting process
            swapped = False

            for j in range(0, len(Z)-i-1):
                # Compare adjacent elements
                if Z[j] > Z[j+1]:
                    # Swap the elements
                    Z[j], Z[j+1] = Z[j+1], Z[j]
                    swapped = True

            # If no two elements were swapped in the inner loop, 
            # the array is already sorted
            if not swapped:
                break
        Z[2],Z[0] = Z[0],Z[2]
        return Z
    
# Define cubic coefficients
def alfa_beta_gamma(A,B,u,w):
    alfa = -1-B+u*B
    beta = A+w*B**2-u*B-u*B**2
    gamma = -A*B-w*B**2-w*B**3
    return alfa,beta,gamma

# Calculate ln(phi) for V and L phase
def calculate_phi(Z,A,B,model):
    if model == "vdW":
        lnphi = Z - 1 - A/Z - math.log(Z-B)
    elif model == "RK":
        lnphi = Z - 1 - A/B * math.log(1+B/Z) - math.log(Z-B)
    elif model == "RKS":
        lnphi = Z - 1 - A/B * math.log(1+B/Z) - math.log(Z-B)
    elif model == "PR":
        lnphi = Z - 1 - math.log(Z-B) - A/B/2/(2**(1/2)) *\
            math.log((Z+B*(1+2**(1/2)))/(Z+B*(1-2**(1/2))))
    
    return lnphi

# Determine a first guess value for P0 or T0, 
# based on Clausius-Clapeyron equation

dHev = R * math.log(Pc/PebN) / (1/TebN-1/Tc)
C1 = math.log(Pc) + dHev/R/Tc
if P == "unknown":
    check = "P"
    P = math.exp(C1-dHev/R/T)
elif T == "unknown":
    check = "T"
    T = dHev/R/(C1-math.log(P))

# van der Waals model
if model == "vdW":
    def equation_solver(fg,x1,check,Tc,Pc):

        if check == "P": P,T = fg,x1
        elif check == "T": T,P = fg,x1

        a = 0.421875*(R*Tc)**2/Pc
        b = 0.125*R*Tc/Pc
        A = a*P/(R*T)**2
        B = b*P/R/T

        u = 0
        w = 0
        alfa,beta,gamma = alfa_beta_gamma(A,B,u,w)
        Zvector = Zroot_calculation(alfa,beta,gamma)
        ZV = max(Zvector)
        ZL = min(Zvector)

        lnphiV = calculate_phi(ZV,A,B,model)
        lnphiL = calculate_phi(ZL,A,B,model)

        return math.exp(lnphiV) / math.exp(lnphiL) - 1
    
    # Solve the 1 eq/1 unk
    if check == "P": 
        fg = P
        parameters = (T,"P",Tc,Pc)
    elif check == "T": 
        fg = T
        parameters = (P,"T",Tc,Pc)

    solution = fsolve(equation_solver, [fg], args = parameters, xtol=1e-10)

elif model == "RK":
    def equation_solver(fg,x1,check,Tc,Pc):

        if check == "P": P,T = fg,x1
        elif check == "T": T,P = fg,x1

        Tr = T/Tc
        k = 1/((Tr)**(1/2))

        a = 0.42748*(R*Tc)**2 * k/Pc
        b = 0.08664*R*Tc/Pc
        A = a*P/(R*T)**2
        B = b*P/R/T

        u = 1
        w = 0
        alfa,beta,gamma = alfa_beta_gamma(A,B,u,w)
        Zvector = Zroot_calculation(alfa,beta,gamma)
        ZV = max(Zvector)
        ZL = min(Zvector)

        lnphiV = calculate_phi(ZV,A,B,model)
        lnphiL = calculate_phi(ZL,A,B,model)

        return math.exp(lnphiV) / math.exp(lnphiL) - 1
    
    # Solve the 1 eq/1 unk
    if check == "P": 
        fg = P
        parameters = (T,"P",Tc,Pc)
    elif check == "T": 
        fg = T
        parameters = (P,"T",Tc,Pc)
    solution = fsolve(equation_solver, [fg], args = parameters, xtol=1e-10)

elif model == "RKS":
    def equation_solver(fg,x1,check,Tc,Pc,omega):

        if check == "P": P,T = fg,x1
        elif check == "T": T,P = fg,x1

        Tr = T/Tc
        S_omega = 0.48 + 1.574*omega - 0.176*omega**2
        k = (1+S_omega*(1-(Tr)**(1/2)))**2

        a = 0.42748*(R*Tc)**2 * k/Pc
        b = 0.08664*R*Tc/Pc
        A = a*P/(R*T)**2
        B = b*P/R/T

        u = 1
        w = 0
        alfa,beta,gamma = alfa_beta_gamma(A,B,u,w)
        Zvector = Zroot_calculation(alfa,beta,gamma)
        ZV = max(Zvector)
        ZL = min(Zvector)

        lnphiV = calculate_phi(ZV,A,B,model)
        lnphiL = calculate_phi(ZL,A,B,model)

        return math.exp(lnphiV) / math.exp(lnphiL) - 1
    
    # Solve the 1 eq/1 unk
    if check == "P": 
        fg = P
        parameters = (T,"P",Tc,Pc,omega)
    elif check == "T": 
        fg = T
        parameters = (P,"T",Tc,Pc,omega)

    solution = fsolve(equation_solver, [fg], args = parameters, xtol=1e-10)

elif model == "PR":
    def equation_solver(fg,x1,check,Tc,Pc,omega):

        if check == "P": P,T = fg,x1
        elif check == "T": T,P = fg,x1

        Tr = T/Tc
        S_omega = 0.37464 + 1.54226*omega - 0.26992*omega**2
        k = (1+S_omega*(1-(Tr)**(1/2)))**2

        a = 0.45724*(R*Tc)**2 * k/Pc
        b = 0.0778*R*Tc/Pc
        A = a*P/(R*T)**2
        B = b*P/R/T

        u = 2
        w = -1
        alfa,beta,gamma = alfa_beta_gamma(A,B,u,w)
        Zvector = Zroot_calculation(alfa,beta,gamma)
        ZV = max(Zvector)
        ZL = min(Zvector)

        lnphiV = calculate_phi(ZV,A,B,model)
        lnphiL = calculate_phi(ZL,A,B,model)

        return math.exp(lnphiV) / math.exp(lnphiL) - 1
    
    # Solve the 1 eq/1 unk
    if check == "P": 
        fg = P
        parameters = (T,"P",Tc,Pc,omega)
    elif check == "T": 
        fg = T
        parameters = (P,"T",Tc,Pc,omega)

    solution = fsolve(equation_solver, [fg], args = parameters, xtol=1e-10)

if check == "P":
        print("Estimated vapor pressure with the "+model+" model: "\
            +str(round(solution[0]/1e5,8))+" [bar]")
elif check == "T": 
        print("Estimated vapor temperature with the "+model+" model: "\
            +str(round(solution[0],8))+" [K]")
