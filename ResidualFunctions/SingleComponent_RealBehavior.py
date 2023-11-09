# Code for the automatic estimation of residual gRes, hRes, sRes
# Single component, real gas/liquid
# Given the conditions of the system (T,P) and the properties of the
# compound (critical T, critical P and the Pitzer acentric factor omega),
# the codes estimates the residual value of Gibbs free energy, enthalpy
# and entropy, in [J/mol]
# Author: Marcello Ferraro
# Last update: 07/11/2023

import sys
import math

# Define constants
R = 8.31446261815324 # [J/mol/K]
pi = 3.14159265358979 # [-]

# Define the model used for the calculation:
# vdW = van der Waals
# RK = Redlich-Kwong
# RKS = Redlich-Kwong-Soave
# PR = Peng-Robinson
# Virial = Virial EoS (Onnes)

model = "RKS"
phase = "L" # Define phase of the component (G=Gas , L=liquid)

### Define data

# Compound properties
Pc = 24.97*1e5 # Critical pressure [Pa]
Tc = 569.40 # Critical temperature [K]
omega = 0.398 # Pitzer acentric factor [-]

# System properties 
P = 2.15*1e5 # System pressure [Pa]
T = 427.9 # System temperature [K]

if model != "vdW" and model != "RK" and model != "RKS" and \
    model != "PR" and model != "Virial":
    print("ERROR: define a valid model!")
    sys.exit()

if phase != "G" and phase != "L":
    print("ERROR: define a valid phase (G/L)!")
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

def alfa_beta_gamma(A,B,u,w):
    alfa = -1-B+u*B
    beta = A+w*B**2-u*B-u*B**2
    gamma = -A*B-w*B**2-w*B**3
    return alfa,beta,gamma

# Calculate the value of Z @ system conditions and the Res. functions
if model == "vdW":
    a = 0.421875*(R*Tc)**2/Pc
    b = 0.125*R*Tc/Pc
    A = a*P/(R*T)**2
    B = b*P/R/T

    u = 0
    w = 0
    alfa,beta,gamma = alfa_beta_gamma(A,B,u,w)
    if phase == "G":
        Z = max(Zroot_calculation(alfa,beta,gamma))
    else:
        Z = min(Zroot_calculation(alfa,beta,gamma))
    
    hRT = Z - 1 - A/Z
    hRes = hRT * R * T

    sR = math.log(Z-B)
    sRes = sR * R

    gRT = Z - 1 - A/Z - math.log(Z-B)
    gRes = gRT * R * T

elif model == "RK":
    Tr = T/Tc
    k = 1/((Tr)**(1/2))

    a = 0.42748*(R*Tc)**2 * k/Pc
    b = 0.08664*R*Tc/Pc
    A = a*P/(R*T)**2
    B = b*P/R/T

    u = 1
    w = 0
    alfa,beta,gamma = alfa_beta_gamma(A,B,u,w)
    if phase == "G":
        Z = max(Zroot_calculation(alfa,beta,gamma))
    else:
        Z = min(Zroot_calculation(alfa,beta,gamma))
    
    hRT = Z - 1 - 3*A/2/B * math.log(1+B/Z)
    hRes = hRT * R * T

    sR = math.log(Z-B) - A/2/B*math.log(1+B/Z)
    sRes = sR * R

    gRT = Z - 1 - A/B*math.log(1+B/Z) - math.log(Z-B)
    gRes = gRT * R * T

elif model == "RKS":
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
    if phase == "G":
        Z = max(Zroot_calculation(alfa,beta,gamma))
    else:
        Z = min(Zroot_calculation(alfa,beta,gamma))
    
    hRT = Z - 1 - A/B * math.log(1+B/Z) * \
        (1 + S_omega*(Tr)**(1/2) / (1+S_omega*(1-(Tr)**(1/2))))
    hRes = hRT * R * T

    sR = math.log(Z-B) - A/B * math.log(1+B/Z) * \
        (S_omega*(Tr)**(1/2) / (1+S_omega*(1-(Tr)**(1/2))))
    sRes = sR * R

    gRT = Z - 1 - A/B*math.log(1+B/Z) - math.log(Z-B)
    gRes = gRT * R * T

elif model == "PR":
    Tr = T/Tc
    S_omega = 0.37464 + 1.54226*omega - 0.26992*omega**2
    k = (1+S_omega*(1-(Tr)**(1/2)))**2

    a = 0.45724*(R*Tc)**2 * k/Pc
    b = 0.0778*R*Tc/Pc
    A = a*P/(R*T)**2
    B = b*P/R/T
    print(A,B,S_omega)
    u = 2
    w = -1
    alfa,beta,gamma = alfa_beta_gamma(A,B,u,w)
    if phase == "G":
        Z = max(Zroot_calculation(alfa,beta,gamma))
    else:
        Z = min(Zroot_calculation(alfa,beta,gamma))
    
    hRT = Z - 1 - A/B/2/(2**(1/2)) * \
        math.log((Z+B*(1+2**(1/2)))/(Z+B*(1-2**(1/2))))*\
        (1 + S_omega*(Tr)**(1/2) / (1+S_omega*(1-(Tr)**(1/2))))
    hRes = hRT * R * T

    sR = math.log(Z-B) - A/B/2/(2**(1/2)) * \
        math.log((Z+B*(1+2**(1/2)))/(Z+B*(1-2**(1/2))))*\
        (S_omega*(Tr)**(1/2) / (1+S_omega*(1-(Tr)**(1/2))))
    sRes = sR * R

    gRT = Z - 1 - math.log(Z-B) - A/B/2/(2**(1/2)) *\
        math.log((Z+B*(1+2**(1/2)))/(Z+B*(1-2**(1/2))))
    gRes = gRT * R * T

elif model == "Virial":

    if phase == "L":
        print("ERROR: the Virial model can be used only for gaseous phase!")
        sys.exit()

    Tr = T/Tc
    Pr = P/Pc
    B0 = 0.083 - 0.422/(Tr**(1.6))
    B1 = 0.139 - 0.172/(Tr**(4.2))
    BPc_RTc = B0 + B1 * omega
    B = BPc_RTc * R * Tc / Pc
    dBdT = R/Pc * (0.675/(Tr**(2.6)) + omega*0.722/(Tr**(5.2)))

    hRT = P/R/T*(B-T*dBdT)
    hRes = hRT * R * T

    sR = -P/R*dBdT
    sRes = sR * R

    gRT = B*P/R/T
    gRes = gRT * R * T

print("Model: "+model)
print("Residual enthalpy hRes ("+phase+" phase): \
"+str(round(hRes,8))+" [J/mol]")
print("Residual entropy sRes ("+phase+" phase): \
"+str(round(sRes,8))+" [J/mol]")
print("Residual Gibbs free energy gRes ("+phase+" phase): \
"+str(round(gRes,8))+" [J/mol]")
