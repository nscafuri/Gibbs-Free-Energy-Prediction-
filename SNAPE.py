import os
import math
import readline
import sys
from scipy import constants
from scipy.constants import pi,c,h,k,Avogadro,R

sp = " "

# convert c from m/s  to cm/s
c_cm = c*100  

# Declare lists
nI, v, ThetaRot, ThetaVib = [], [], [], []


# Check if input files exist
if(os.path.exists("snape.inp") == False):
  print("<< ERROR >>: snape.inp file does not exist.")
  quit()

if(os.path.exists("vibrations") == False):
  print("<< ERROR >>: vibrations file does not exist.")
  quit()


# Open snape.inp  and read parameters
with open("snape.inp","r") as s:

     for line in s:
         values = line.split()

         if line.startswith("T"):
            T = float(values[1])
            RT = R*T

         if line.startswith("P"):
            P = float(values[1])

         if line.startswith("sigma"):
            sigma = float(values[1])

         if line.startswith("P"):
            P = float(values[1])

         if line.startswith("M"):
            M = float(values[1])
     
         if line.startswith("Escf"):
            Escf = float(values[1])

         if line.startswith("nI"):
            nI = float(values[1])

         if line.startswith("I1"):
            I1 = nI.append(float(values[1]))
     
         if line.startswith("I2"):
            I2 = nI.append(float(values[1]))

         if line.startswith("I3"):
            I3 = nI.append(float(values[1]))

     nI.sort()

# Open vibrations file and read values
with open("vibrations","r") as vib:
     for line in vib:
         if(len(line)>1):
            v.append(float(line))
     v.sort()

def Rotational_Temperatures(T,nI):
    global ThetaRot
    """
    Calculation of the rotational temperatures 

    Args = [T, nI]: Temperature, Moments of Inertia
    Returns = [ThetaRot]: Rotational Temperatures (1/K)  
    """
    if len(nI) > 0:
       for i in range (len(nI)):
           ThetaRot.append((h**2)/(8*pi*pi*nI[i]*k))
    else:
        ThetaRot = []
    return ThetaRot       


def Vibrational_Temperatures(T,v,c_cm):
    global ThetaVib
    """
    Calculation of the vibrational temperatures 

    Args = [T, v, c_cm]: Temperature, vibrations, speed of light
    Returns = [ThetaVib]: Vibrational Temperatures  (1/K) 
    """

    for i in range(len(v)):
        ThetaVib.append((h*c_cm*v[i])/k)
    
    return ThetaVib
    

def Zero_Point_Energy(v,c_cm):
    global ZPE
    """
    Calculation of the Zero Point Energy 

    Args = [T, v, c_cm, Avogadro]: vibrations, speed of light, Avogadr's number 
    Returns = [ZPE]: Zero Point Energy (J/mol)
    """
    ZPE = 0.
    for i in range(len(v)):
        ZPE = ZPE + (0.5*h*c_cm*v[i])
    
    ZPE = ZPE*Avogadro
    #ZPE = ZPE/1000
    
    return ZPE

def Vibrational_Energy(T,ThetaVib):
    global Uvib
    """
    Calculation of Vibrational Energy (J/(mol))

    Args = [T, v,ThetaVib]: Temperature, vibrations,speed of light, Vibrational temperatures 
    Returns = [Uvib]: Vibrational Energy (J/mol) 
    """
    a = 0.
    for i in range(len(v)):
        a = a + ThetaVib[i]*((0.5 + (1/(math.exp(ThetaVib[i]/T) -1))))
        Uvib = a*R

    return Uvib


def Vibrational_Entropy(T,ThetaVib):
    global Svib
    """
    Calculation of Vibrational Entropy 

    Args = [T, v,ThetaVib]: Temperature, vibrations, Vibrational temperatures 
    Returns = [Svib]: Vibrational Entropy (J/mol*K) 
    """
    b = 0.
    for i in range(len(v)):
        b = b + (ThetaVib[i]/T)/((math.exp(ThetaVib[i]/T - 1) - math.log(1 - math.exp(-ThetaVib[i]/T))))
    Svib = b*R
    
    return Svib 

def Translational_Partition_Function(T,M,P):
    global qtr, Utr, Str
    """
    Calculation of the: Translational Parition Function (if gas phase molecule present)
                        Translational Energy
                        Translational Entropy
    Args = [T, M, P]: Temperature, Mass, Pression
    Returns = [qtr, Utr, Str]: Translational Parition Function (adimensional), Trans. energy (J/mol), Trans. entropy (J/mol K))
    """
    if len(nI) > 0:
       M = M*(1.661e-27)
       qtr = ((2*pi*M*k*T)/(h*h))
       qtr = qtr**1.5
       qtr = qtr*((k*T)/P)
   
       Str = R*((math.log(qtr) + 1 + 1.5))
       Utr = 1.5*(R*T)
    else:
        qtr = 0
        Utr = 0
        Str = 0
    return qtr, Utr, Str


def Rotational_Partition_Function(T, ThetaRot, sigma):
    global qrot, Urot, Srot
    """
    Calculation of the: Rotational Parition Function (if gas phase molecule present)
                        Rotational Energy
                        Rotational Entropy
    Args = [T, M, P]: Temperature, Mass, Pression
    Returns = [qrot, Urot, Str]: Translational Parition Function (adimensional), Rot. energy (J/mol), Rot. entropy (J/mol K))
    """
    if len(nI) > 0:
     if len(ThetaRot) == 1:  # Linear molecule 
        qrot = (1/sigma)/(T/ThetaRot[0])
        Srot = R*(math.log(qrot) + 1)
        Urot = 1.5*(R*T)

     elif len(ThetaRot) > 1:
          qrot = ((pi**0.5)/sigma)*((T**1.5)/(ThetaRot[0]*ThetaRot[1]*ThetaRot[2])**0.5)
          Srot = R*(math.log(qrot) + 1.5)
          Urot = 1.5*(R*T)
    else:
     qrot = 0
     Srot = 0
     Urot = 0
    return qrot, Urot, Srot

def write_ouput(Escf,T,sigma,M,nI,v,P,ThetaRot, ThetaVib, qrot, qtr, ZPE, Uvib, Utr, Urot, Str, Svib, Srot,sp):
    """
    Write the output file (snape.out)
    """
    #emptyline = print()
    with open("snape.out","w") as o:
         o.write("********SNAPE OUTPUT******************\t")
         o.write("\n")
         o.write("INPUT PARAMETERS\n")
         o.write("T= {} {}".format(T, "K",) + "\n")
         o.write("P= {} {}".format(P, "Pascal",) + "\n")
         o.write("M= {} {}".format(M, "a.u.",) + "\n")
         o.write("sigma= {}".format(sigma) + "\n")
         o.write("Escf= {} {}".format(Escf, "Hartree") + "\n")
         o.write("\n")
         
         if len(nI) > 0:
             o.write("MOMENTS OF INERTIA {} = {}".format("(Kg/m**2)", len(nI)) + "\n")
             for i in range (len(nI)):
                 o.write("{} {}".format(i, nI[i]) + "\n")
             o.write("\n")
             
             o.write("TRANLATIONAL ENERGY = {} {}".format(Utr, "J/mol") + "\n")
             o.write("TRANLATIONAL ENTROPY = {} {}".format(Str, "J/mol*K") + "\n")
             

             
         else:
             o.write("NO MOMENT OF INERTIA" + "\n") 
            
         o.write("\n")
         o.write("RT = {} {}".format(RT, "J/mol") + "\n") 
         o.write("\n")
         if len(nI) > 0:
             o.write("ROTATIONAL TEMPERATURES(1/K) = {}".format(len(ThetaRot)) + "\n")
             for i in range (len(ThetaRot)):
                 o.write("{} = {}".format(i, ThetaRot[i]) + "\n") 
             o.write("\n")
                 
             o.write("ROTATIONAL ENERGY = {} {}".format(Urot, "J/mol") + "\n") 
             o.write("ROTATIONAL ENTROPY = {} {}".format(Srot, "J/mol K") + "\n")
         else:
             o.write("NO ROTATIONAL TEMPERATURE" + "\n") 

         o.write("\n")
         o.write("VIBRATIONAL TEMPERATURES(1/K) = {}".format(len(ThetaVib)) + "\n") 
         for i in range (len(ThetaVib)):
             o.write("{} = {}".format(i, ThetaVib[i]) + "\n")
         o.write("\n")
         o.write("ZERO POINT ENERGY (ZPE) = {} {}".format(ZPE, "J/mol") + "\n") 
         o.write("VIBRATIONAL ENERGY = {} {}".format(Uvib, "J/mol") + "\n") 
         o.write("VIBRATIONAL ENTROPY = {} {}".format(Svib, "J/mol*K") + "\n") 
         
         Escf = Escf*(2623*1000) 
         Stot = Str + Srot + Svib
         Utot = Escf + Utr + Urot + Uvib
         H = Utot + (R*T)
         G = H - (T*Stot)

         o.write("\n")
         o.write("Stot = {} {}".format(Stot, "J/mol*K") + "\n") 
         o.write("Utot = {} {}".format(Utot, "J/mol*K") + "\n") 
         o.write("H = {} {}".format(H, "J/mol") + "\n") 
         o.write("G = {} {}".format(G, "J/mol") + "\n") 
         
Rotational_Temperatures(T,nI)
Rotational_Partition_Function(T, ThetaRot, sigma)
Translational_Partition_Function(T,M,P)

Zero_Point_Energy(v,c_cm)
Vibrational_Temperatures(T,v,c_cm)
Vibrational_Energy(T,ThetaVib)
Vibrational_Entropy(T,ThetaVib)

write_ouput(Escf,T,sigma,M,nI,v,P,ThetaRot, ThetaVib, qrot, qtr, ZPE, Uvib, Utr, Urot, Str, Svib, Srot,sp)

print("NORMAL COMPLETATION: data written in snape.out file.")






