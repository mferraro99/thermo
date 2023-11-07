# Code for the automatic estimation of volumetric properties (T,P,v)
# Single component, real gas
# Given two state properties (T,P or T,v or P,v) the code applies
# different models, selected by the user, for the estimation of the
# third property
# Author: Marcello Ferraro
# Last update: 07/11/2023

import sys
import math
from scipy.optimize import fsolve

# Define constants
R = 8.31446261815324 # [J/mol/K]
pi = 3.14159265358979 # [-]

# Define the model used for the calculation:
# IG = Ideal Gas
# vdW = van der Waals
# RK = Redlich-Kwong
# RKS = Redlich-Kwong-Soave
# PR = Peng-Robinson
# Virial = Virial EoS (Onnes)

model = "RKS"

### Define data

# System properties --> "unknown" is the quantity we want to calculate
P = 2.15*1e5 # System pressure [Pa]
T = 427.9 # System temperature [K]
v = "unknown" # System molar volume [m3/mol]

# Compound properties
Pc = 24.97*1e5 # Critical pressure [Pa]
Tc = 569.40 # Critical temperature [K]
omega = 0.398 # Pitzer acentric factor [-]

### Define a first guess value for our quantity (Ideal gas law)

if P == "unknown": 
    P = R*T/v
    check = "P"
elif T == "unknown": 
    T = P*v/R
    check = "T"
elif v == "unknown": 
    v = R*T/P
    check = "v"

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
    
def alfa_beta_gamma(A,B,u,w):
    alfa = -1-B+u*B
    beta = A+w*B**2-u*B-u*B**2
    gamma = -A*B-w*B**2-w*B**3
    return alfa,beta,gamma

### Ideal Gas model
if model == "IG":
    if P == "unknown": 
        print("Estimated pressure with the ideal gas law: "\
              +str(round(P/1e5,5))+" [bar]")
    elif T == "unknown": 
        print("Estimated temperature with the ideal gas law: "\
              +str(round(T,5))+" [K]")
    elif v == "unknown": 
        print("Estimated molar volume with the ideal gas law: "\
              +str(round(v,5))+" [m3/mol]")
    
    sys.exit()

### van der Waals model
if model == "vdW":
    def equation_solver(fg,x1,x2,check,Tc,Pc):

        if check == "P": P,T,v = fg,x1,x2
        elif check == "T": T,P,v = fg,x1,x2
        elif check == "v": v,T,P = fg,x1,x2

        a = 0.421875*(R*Tc)**2/Pc
        b = 0.125*R*Tc/Pc
        A = a*P/(R*T)**2
        B = b*P/R/T

        u = 0
        w = 0
        alfa,beta,gamma = alfa_beta_gamma(A,B,u,w)
        Z = max(Zroot_calculation(alfa,beta,gamma))

        if check == "P": eq = P-Z*R*T/v
        elif check == "T": eq = T-P*v/R/Z
        elif check == "v": eq = v-R*T*Z/P

        return eq
    
    # Solve the 1 eq/1 unk
    if check == "P": 
        fg = P
        parameters = (T,v,"P",Tc,Pc)
    elif check == "T": 
        fg = T
        parameters = (P,v,"T",Tc,Pc)
    elif check == "v": 
        fg = v
        parameters = (T,P,"v",Tc,Pc)

    solution = fsolve(equation_solver, [fg], args = parameters, xtol=1e-10)

        
### Redlich-Kwong model
if model == "RK":
    def equation_solver(fg,x1,x2,check,Tc,Pc):

        if check == "P": P,T,v = fg,x1,x2
        elif check == "T": T,P,v = fg,x1,x2
        elif check == "v": v,T,P = fg,x1,x2

        Tr = T/Tc
        k = 1/((Tr)**(1/2))

        a = 0.42748*(R*Tc)**2 * k/Pc
        b = 0.08664*R*Tc/Pc
        A = a*P/(R*T)**2
        B = b*P/R/T

        u = 1
        w = 0
        alfa,beta,gamma = alfa_beta_gamma(A,B,u,w)
        Z = max(Zroot_calculation(alfa,beta,gamma))

        if check == "P": eq = P-Z*R*T/v
        elif check == "T": eq = T-P*v/R/Z
        elif check == "v": eq = v-R*T*Z/P

        return eq
    
    # Solve the 1 eq/1 unk
    if check == "P": 
        fg = P
        parameters = (T,v,"P",Tc,Pc)
    elif check == "T": 
        fg = T
        parameters = (P,v,"T",Tc,Pc)
    elif check == "v": 
        fg = v
        parameters = (T,P,"v",Tc,Pc)

    solution = fsolve(equation_solver, [fg], args = parameters, xtol=1e-10)

### Redlich-Kwong-Soave model
if model == "RKS":
    def equation_solver(fg,x1,x2,check,Tc,Pc,omega):

        if check == "P": P,T,v = fg,x1,x2
        elif check == "T": T,P,v = fg,x1,x2
        elif check == "v": v,T,P = fg,x1,x2

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
        Z = max(Zroot_calculation(alfa,beta,gamma))

        if check == "P": eq = P-Z*R*T/v
        elif check == "T": eq = T-P*v/R/Z
        elif check == "v": eq = v-R*T*Z/P

        return eq
    
    # Solve the 1 eq/1 unk
    if check == "P": 
        fg = P
        parameters = (T,v,"P",Tc,Pc,omega)
    elif check == "T": 
        fg = T
        parameters = (P,v,"T",Tc,Pc,omega)
    elif check == "v": 
        fg = v
        parameters = (T,P,"v",Tc,Pc,omega)

    solution = fsolve(equation_solver, [fg], args = parameters, xtol=1e-10)

### Peng-Robinson model
if model == "PR":
    def equation_solver(fg,x1,x2,check,Tc,Pc,omega):

        if check == "P": P,T,v = fg,x1,x2
        elif check == "T": T,P,v = fg,x1,x2
        elif check == "v": v,T,P = fg,x1,x2

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
        Z = max(Zroot_calculation(alfa,beta,gamma))

        if check == "P": eq = P-Z*R*T/v
        elif check == "T": eq = T-P*v/R/Z
        elif check == "v": eq = v-R*T*Z/P

        return eq
    
    # Solve the 1 eq/1 unk
    if check == "P": 
        fg = P
        parameters = (T,v,"P",Tc,Pc,omega)
    elif check == "T": 
        fg = T
        parameters = (P,v,"T",Tc,Pc,omega)
    elif check == "v": 
        fg = v
        parameters = (T,P,"v",Tc,Pc,omega)

    solution = fsolve(equation_solver, [fg], args = parameters, xtol=1e-10)

### Virial equation model
if model == "Virial":
    def equation_solver(fg,x1,x2,check,Tc,Pc,omega):

        if check == "P": P,T,v = fg,x1,x2
        elif check == "T": T,P,v = fg,x1,x2
        elif check == "v": v,T,P = fg,x1,x2

        Tr = T/Tc
        Pr = P/Pc
        B0 = 0.083 - 0.422/(Tr**(1.6))
        B1 = 0.139 - 0.172/(Tr**(4.2))
        BPc_RTc = B0 + B1 * omega
        Z = 1 + BPc_RTc *Pr/Tr

        if check == "P": eq = P-Z*R*T/v
        elif check == "T": eq = T-P*v/R/Z
        elif check == "v": eq = v-R*T*Z/P

        return eq
    
    # Solve the 1 eq/1 unk
    if check == "P": 
        fg = P
        parameters = (T,v,"P",Tc,Pc,omega)
    elif check == "T": 
        fg = T
        parameters = (P,v,"T",Tc,Pc,omega)
    elif check == "v": 
        fg = v
        parameters = (T,P,"v",Tc,Pc,omega)

    solution = fsolve(equation_solver, [fg], args = parameters, xtol=1e-10)


if check == "P":
    try:
        print("Estimated pressure with the "+model+" model: "\
            +str(round(solution[0]/1e5,8))+" [bar]")
    except:
        print("ERROR: check model definition!")
elif check == "T": 
    try:
        print("Estimated temperature with the "+model+" model: "\
            +str(round(solution[0],8))+" [K]")
    except:
        print("ERROR: check model definition!")
    
elif check == "v": 
    try:
        print("Estimated molar volume with the "+model+" model: "\
            +str(round(solution[0],8))+" [m3/mol]")
    except:
        print("ERROR: check model definition!")
