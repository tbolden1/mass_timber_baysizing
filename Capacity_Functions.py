# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 14:50:59 2024

@author: tbolden
"""

# Wood Design Section Functions

# Beam Flexural, Shear Capacity Function
# -----------------------------------------------------------------------------
def flexural_capacity_beam_ASD(beam, L, s, CD=1.0, CL=1.0, Cfu=1.0, Cc=1.0, CI=1.0):  # Can make this a method of the beam class to have functionality in the class
    
  #Need to incorporate Beam Stability Factor
    # if d>b:
    #     CL = 1.0
    # elif compression_bracing = 'Yes':
    #     CL = 1.0
    
    #Determine the Volume Factor
    if beam.species == "SPF":
        x = 20
    else:
        x = 10
        
    CV = (21/(L/12))**(1/x)*(12/beam.d)**(1/x)*(5.125/beam.b)**(1/x)
    
    #Determine Adjusted Design Value F'b
    F_b = beam.Fb * CD * beam.CM * beam.Ct * min(CL,CV) * Cfu * Cc * CI
    
    #Allowable Flexural Moment
    Ma = F_b * beam.Sx
    
    return Ma

def shear_capacity_beam_ASD(beam, CD=1.0):
    
    #Shear Reduction Factor Cvr per AWC 5.3.10
    if beam.conn_type == 'seat':
        Cvr = 1.0
        de = beam.d
    if beam.conn_type == 'shear':
        Cvr = 0.72
        de = 0.8*beam.d     
    
    #Shear Parrallel to Grain
    F_v = beam.Fv*beam.CM*beam.Ct*Cvr
    print("Value of 'F_v':", F_v)
    
    #Allowable Shear
    Va = (2/3*F_v)*beam.b*de*(de/beam.d)**2
    
    return Va

# Plank Flexural, Shear Capacity Function
# -----------------------------------------------------------------------------

def flexural_capacity_plank_ASD(CLT,CD=1.0):
    Ma = CLT.fbseff_0*CLT.CM*CLT.Ct*CD
    
    return Ma

def shear_capacity_plank_ASD(CLT,CD=1.0):
    Va = CLT.vs_0*CLT.CM*CLT.Ct*CD
    
    return Va
    