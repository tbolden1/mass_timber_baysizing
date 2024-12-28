# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:08:28 2024

@author: tbolden
"""

# Define Member Classes that can be referenced, for function. Contains innate properties of the components of a bay system

#Beam Member Class
class Beam:
    def __init__(self, Fb, Fv, E, CM, Ct, species, conn_type, max_depth = 120):
        self.Fb = Fb                # Allowable bending stress (in psi)
        self.Fv = Fv                # Allowable shear stress (in psi)
        self.E = E                  # Modulus of elasticity (in psi)
        #self.b = b                  # Width of the beam (in inches)
        #self.d = d                  # Depth of the beam (in inches)
        self.CM = CM                # Wet Service Factor  *(Can these be input in a global properties)*
        self.Ct = Ct                # Temperature Factor  *(Can these be input in a global properties)*
        self.species = species      # Wood Species
        self.conn_type = conn_type  # Beam connection type: seat, shear (could this be a separate class that calculates the connection?)
        #self.Sx = b * (d ** 2) / 6  # Section modulus (in^3)
        self.max_depth = max_depth
    
#Column Member Class
class Column:
    def __init__(self, Fbx, Fby, Fvx, Fvy, Fc, E, Emin, CM, Ct, species):  
        self.Fbx = Fbx              # Allowable bending stress (in psi)
        self.Fby = Fby              # Allowable bending stress (in psi)
        self.Fvx = Fvx              # Allowable shear stress (in psi)
        self.Fvy = Fvy              # Allowable shear stress (in psi)
        self.Fc = Fc                # Allowable compression parrallel to grain (in psi)
        self.E = E                  # Modulus of elasticity for deflection (in psi)
        self.Emin = Emin            # Modulus of elasticity for stability (in psi)
        self.CM = CM                # Wet Service Factor  *(Can these be input in a global properties)*
        self.Ct = Ct                # Temperature Factor  *(Can these be input in a global properties)*
        self.species = species      # Wood Species

#CLT Panel Member Class
class CLTSize:
    def __init__(self, row, CM, Ct):
        self.clt_size = row['CLT Size']
        self.layers = row['Layers']
        self.depth = row['Depth']
        self.tlayer = row['tlayer']
        self.wt = row['Wt']
        self.fbseff_0 = row['FbSeff_0']
        self.eieff_0 = row['Eieff_0']
        self.gaeff_0 = row['Gaeff_0']
        self.vs_0 = row['Vs_0']
        self.fbseff_90 = row['FbSeff_90']
        self.eieff_90 = row['Eieff_90']
        self.gaeff_90 = row['Gaeff_90']
        self.vs_90 = row['Vs_90']
        self.EIapp_0 = row['EIapp_0']
        self.CM = CM                        # Wet Service Factor  *(Can these be input in a global properties)*
        self.Ct = Ct                        # Temperature Factor  *(Can these be input in a global properties)*
    
    def __repr__(self):
        return (f"CLTSize(clt_size={self.clt_size}, layers={self.layers}, depth={self.depth}, tlayer={self.tlayer}, "
                f"wt={self.wt}, fbseff_0={self.fbseff_0}, eieff_0={self.eieff_0}, gaeff_0={self.gaeff_0}, vs_0={self.vs_0}, "
                f"fbseff_90={self.fbseff_90}, eieff_90={self.eieff_90}, gaeff_90={self.gaeff_90}, vs_90={self.vs_90})")

