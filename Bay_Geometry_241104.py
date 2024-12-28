# -*- coding: utf-8 -*-
"""
Created on Wed May 22 12:53:46 2024

@author: tbolden
"""
import numpy as np
import pandas as pd
import openpyxl
from Member_Classes import Beam, Column, CLTSize
from Capacity_Functions import flexural_capacity_plank_ASD, shear_capacity_plank_ASD
from Support_Functions import deflection_due_to_point_load, calculate_total_deflection
# ----------------------------------------------------------------------------
#Read Databases for Member Design
# ----------------------------------------------------------------------------
#Define Global Properties
CM = 1.0        # Wet Service Factor  
Ct = 1.0        # Temperature Factor
gamma_GLB = 36  # [pcf] Unit Weight of Glulam

DBL_Grider = 1  # Consider DBL Girder? 1 = Yes , 0 = No

# ----------------------------------------------------------------------------
#Geometry of Bays
# ----------------------------------------------------------------------------

# Load your matrix data

input_bay_df = pd.read_excel('Input_Bay_Matrix.xlsx')
results_xlsx = 'Results_Matrix_Girder.xlsx'

# input_bay_df = pd.read_excel('Input_Bay_Matrix_Test.xlsx')
# results_xlsx = 'Results_Bay_Matrix_Test.xlsx'

# ----------------------------------------------------------------------------
# Importing Member Section Properties
# ----------------------------------------------------------------------------

# Glulam Properties ----------------------------------------------------------
glb_sections = pd.read_csv('NDS_GLB_STDSizes.csv')
DBL_glb_sections = pd.read_csv('NDS_DBL_GLB_STDSizes.csv')

# CLT Plank Properties -------------------------------------------------------
plank_DF = pd.read_csv('SmartLam_Alt_CLT.csv') # Read database of CLT Plank Design Properties (Currently Loading SmartLam DF Alt Properties)  ** Update to PRG-320 **

#Define Beam/Girder Properties -----------------------------------------------
max_BM_depth = 36
max_GR_depth = 36
beam_24F_V4 = Beam(2400, 265, 1800000, CM, Ct, 'DF', 'seat', max_BM_depth)
girder_24F_V4 = Beam(2400, 265, 1800000, CM, Ct, 'DF', 'seat', max_GR_depth)

#Define Column Properties ----------------------------------------------------
col_lay3 = Column(2000, 2100, 265, 230, 2300, 2000000, 1000000, CM, Ct, 'DF')

# ----------------------------------------------------------------------------
# Dead Load, Live Load , Deflection Criteria 
# ----------------------------------------------------------------------------
t_topping = 3       #Topping Slab Thickness
lambda_conc = 115   #Lightweight Concrete Topping
super_DL = 12       #Superimposed Deadload
# ----------------------------------------------------------------------------
LL = 80 #Liveload
# ----------------------------------------------------------------------------    
Def_LL_lim = 360
Def_LT_lim = 240


# Placeholder list to store the results
results = []

def min_beam_size(Vu_beam, Mu_beam, beam, glb_sections, w_beam_DL, w_beam_LL, CD=1.0, CL=1.0, Cfu=1.0, Cc=1.0, CI=1.0):
    
    gamma_GLB = 36 # [pcf] Unit Weight of Glulam
    
    #Determine the Volume Factor
    if beam.species == "SPF":
        x = 20
    else:
        x = 10
    beam_GLB_sections = glb_sections
    beam_GLB_sections['CV'] = (21/(W))**(1/x)*(12/beam_GLB_sections['d'])**(1/x)*(5.125/beam_GLB_sections['b'])**(1/x)
    
    #Shear Reduction Factor Cvr per AWC 5.3.10
    if beam.conn_type == 'seat':
        Cvr = 1.0
        beam_GLB_sections['de'] = beam_GLB_sections.d 
    if beam.conn_type == 'shear':
        Cvr = 0.72
        beam_GLB_sections['de'] = 0.8*beam_GLB_sections.d 
    
    beam_GLB_sections['F_b'] = beam.Fb * CD * beam.CM * beam.Ct * np.minimum(CL,beam_GLB_sections['CV']) * Cfu * Cc * CI
    F_v = Cvr * beam.CM * beam.Ct * beam.Fv
    beam_GLB_sections['M_total'] = ( beam_GLB_sections.A / 144 * gamma_GLB ) * W ** 2 / 8 + Mu_beam
    beam_GLB_sections['V_total'] = ( beam_GLB_sections.A / 144 * gamma_GLB ) * W / 2 + Vu_beam
    beam_GLB_sections['f_v'] = 3 * beam_GLB_sections.V_total / ( 2 * beam_GLB_sections.b * beam_GLB_sections.de)
    beam_GLB_sections['Sx_min'] = beam_GLB_sections['M_total'] * 12 / beam_GLB_sections['F_b']
    beam_GLB_sections['Def_LL'] = 5 * w_beam_LL / 12 * (W * 12) ** 4 / (384 * beam.E * beam_GLB_sections.Ix)
    beam_GLB_sections['Def_LT'] = 5 * (0.5 * (w_beam_DL + ( beam_GLB_sections.A / 144 * gamma_GLB )) + w_beam_LL) / 12 * (W * 12) ** 4 / (384 * beam.E * beam_GLB_sections.Ix)
    
    condition_1 = beam_GLB_sections['Sx_min'] < beam_GLB_sections['Sx'] #Flexural Capacity Check
    condition_2 = beam_GLB_sections['f_v'] < F_v                        #Shear Capacity Check
    condition_3 = beam_GLB_sections['d'] < beam.max_depth               #Max Depth Check
    condition_4 = (W * 12) / beam_GLB_sections['Def_LL'] > Def_LL_lim   #Live Load Deflection Check
    condition_5 = (W * 12) / beam_GLB_sections['Def_LT'] > Def_LT_lim   #Long-Term Deflection Check
    
    
    filtered_df = beam_GLB_sections[condition_1 & condition_2 & condition_3 & condition_4 & condition_5]
    
    if filtered_df.empty:
        # Return a row of zeros if no valid girder size is found
        zero_row = pd.Series(0, index=beam_GLB_sections.columns)
        return zero_row, zero_row
    
    
    min_vol_index = filtered_df['A'].idxmin()
    min_d_index = filtered_df['d'].idxmin()
    min_vol_row = filtered_df.loc[min_vol_index]
    min_d_row = filtered_df.loc[min_d_index]

    return min_vol_row, min_d_row

# ----------------------------------------------------------------------------
# Function to determine minimum grider size - Considers the Length Equal to L (length)

def min_girder_size(Vu_beam, Mu_beam, beam, glb_sections, DLB_glb_sections, point_loads_LT, point_loads_LL, DBL_Grider=0, CD=1.0, CL=1.0, Cfu=1.0, Cc=1.0, CI=1.0):
    
    gamma_GLB = 36 # [pcf] Unit Weight of Glulam
    
    #Determine the Volume Factor
    if beam.species == "SPF":
        x = 20
    else:
        x = 10
    beam_GLB_sections = glb_sections
    DBL_bm_GLB_sections = DBL_glb_sections
    beam_GLB_sections['CV'] = (21/(L))**(1/x)*(12/beam_GLB_sections['d'])**(1/x)*(5.125/beam_GLB_sections['b'])**(1/x)
    DBL_bm_GLB_sections['CV'] = (21/(L))**(1/x)*(12/beam_GLB_sections['d'])**(1/x)*(5.125/beam_GLB_sections['b'])**(1/x)
    
    #Shear Reduction Factor Cvr per AWC 5.3.10
    if beam.conn_type == 'seat':
        Cvr = 1.0
        beam_GLB_sections['de'] = beam_GLB_sections.d
        DBL_bm_GLB_sections['de'] = DBL_bm_GLB_sections.d
    if beam.conn_type == 'shear':
        Cvr = 0.72
        beam_GLB_sections['de'] = 0.8*beam_GLB_sections.d
        DBL_bm_GLB_sections['de'] = 0.8*DBL_bm_GLB_sections.d
    
    beam_GLB_sections['F_b'] = beam.Fb * CD * beam.CM * beam.Ct * np.minimum(CL,beam_GLB_sections['CV']) * Cfu * Cc * CI
    F_v = Cvr * beam.CM * beam.Ct * beam.Fv
    beam_GLB_sections['M_total'] = ( beam_GLB_sections.A / 144 * gamma_GLB ) * L ** 2 / 8 + Mu_beam
    beam_GLB_sections['V_total'] = ( beam_GLB_sections.A / 144 * gamma_GLB ) * L / 2 + Vu_beam
    beam_GLB_sections['f_v'] = 3 * beam_GLB_sections.V_total / ( 2 * beam_GLB_sections.b * beam_GLB_sections.de)
    beam_GLB_sections['Sx_min'] = beam_GLB_sections['M_total'] * 12 / beam_GLB_sections['F_b']
    beam_GLB_sections['n_section'] = 1
    
    condition_1 = beam_GLB_sections['Sx_min'] < beam_GLB_sections['Sx'] #Flexural Capacity Check
    condition_2 = beam_GLB_sections['f_v'] < F_v                        #Shear Capacity Check
    condition_3 = beam_GLB_sections['d'] < beam.max_depth               #Max Depth Check
    
    filtered_df = beam_GLB_sections[condition_1 & condition_2 & condition_3]
    
    # Check deflection criteria for each filtered girder size
    deflection_criteria = []
    
    for i, row in filtered_df.iterrows():
    # Calculate the maximum deflection due to long term loads           
        max_deflection_long_term = calculate_total_deflection(L * 12, beam.E, row['Ix'], point_loads_LT)
        max_deflection_live_load = calculate_total_deflection(L * 12, beam.E, row['Ix'], point_loads_LT)
                    
        # Check both long term and live load deflections
        deflection_check = (max_deflection_long_term < (L * 12) / Def_LT_lim) and (max_deflection_live_load < (L * 12) / Def_LL_lim)
        deflection_criteria.append(deflection_check)
    
    # Filter further by deflection criteria (both load cases)
    filtered_df = filtered_df[deflection_criteria]
    
    if filtered_df.empty:
        
        if DBL_Grider == 1:
            DBL_bm_GLB_sections['F_b'] = beam.Fb * CD * beam.CM * beam.Ct * np.minimum(CL,beam_GLB_sections['CV']) * Cfu * Cc * CI
            F_v = Cvr * beam.CM * beam.Ct * beam.Fv
            DBL_bm_GLB_sections['M_total'] = ( DBL_bm_GLB_sections.A / 144 * gamma_GLB ) * L ** 2 / 8 + Mu_beam
            DBL_bm_GLB_sections['V_total'] = ( DBL_bm_GLB_sections.A / 144 * gamma_GLB ) * L / 2 + Vu_beam
            DBL_bm_GLB_sections['f_v'] = 3 * DBL_bm_GLB_sections.V_total / ( 2 * DBL_bm_GLB_sections.b * DBL_bm_GLB_sections.de)
            DBL_bm_GLB_sections['Sx_min'] = DBL_bm_GLB_sections['M_total'] * 12 / DBL_bm_GLB_sections['F_b']
            DBL_bm_GLB_sections['n_section'] = 2
            
            condition_1 = DBL_bm_GLB_sections['Sx_min'] < DBL_bm_GLB_sections['Sx'] #Flexural Capacity Check
            condition_2 = DBL_bm_GLB_sections['f_v'] < F_v                        #Shear Capacity Check
            condition_3 = DBL_bm_GLB_sections['d'] < beam.max_depth               #Max Depth Check
            
            filtered_df = DBL_bm_GLB_sections[condition_1 & condition_2 & condition_3]
            for i, row in filtered_df.iterrows():
            # Calculate the maximum deflection due to long term loads           
                max_deflection_long_term = calculate_total_deflection(L * 12, beam.E, row['Ix'], point_loads_LT)
                max_deflection_live_load = calculate_total_deflection(L * 12, beam.E, row['Ix'], point_loads_LT)
                            
                # Check both long term and live load deflections
                deflection_check = (max_deflection_long_term < (L * 12) / Def_LT_lim) and (max_deflection_live_load < (L * 12) / Def_LL_lim)
                deflection_criteria.append(deflection_check)
                # Filter further by deflection criteria (both load cases)
            filtered_df = filtered_df[deflection_criteria]
            
        else:    
            # Return a row of zeros if no valid girder size is found
            zero_row = pd.Series(0, index=beam_GLB_sections.columns)
            return zero_row, zero_row
    
    
    min_vol_index = filtered_df['A'].idxmin()
    min_d_index = filtered_df['d'].idxmin()
    min_vol_row = filtered_df.loc[min_vol_index]
    min_d_row = filtered_df.loc[min_d_index]

    return min_vol_row, min_d_row

def min_column_size(Pu_col, le, column, glb_sections, beam_b=0, girder_b=0, CD=1.0):
    
    col_sections = glb_sections        
    Fc_star = column.Fc * column.CM * column.Ct * CD
    c = 0.9
    col_sections['F_cE'] = 0.822 * column.Emin / (le / col_sections.d) ** 2
    col_sections['C_P'] = (1 + (col_sections.F_cE/Fc_star))/(2 * c) - (((1 + (col_sections.F_cE/Fc_star))/ (2 * c)) ** 2 - (col_sections.F_cE/Fc_star)/c) ** 0.5
    col_sections['F_c'] = column.Fc * column.CM * column.Ct * CD * col_sections.C_P
    col_sections['f_c'] = Pu_col / col_sections.A       
    condition_1 = col_sections['f_c'] < col_sections['F_c'] #Flexural Capacity Check
    condition_2 = col_sections.d < 1.5 * col_sections.b
    condition_3 = col_sections.d >= max(beam_b, girder_b)
    condition_4 = col_sections.b >= min(beam_b, girder_b)
    
    filtered_df = col_sections[condition_1 & condition_3 & condition_2 & condition_4]
                
    min_index = filtered_df['A'].idxmin()
    min_row = filtered_df.loc[min_index]
    
    return min_row

# Iterate over each row in the matrix
for _, row in input_bay_df.iterrows():
    # Extract the parameters for this row
    W = row['W']                    #Bay Width [ft]
    L = row['L']                    #Bay Length [ft]
    H = row['H']                    #Bay Height [ft]
    n_beam = row['n_beam']          #Number beam bays in span
    s_beam = L/n_beam               #Spacing of beams
    bay_type = row['bay_type']      #Bay Layout Type:
                                        #1-Beam spans the Width, Girder Spans the Length, Columns Only Located at the corners
                                        #2-Beam spans the Width, drop Girder, Columns located at every bay
    
    # ----------------------------------------------------------------------------
    # CLT Plank Design - Determine minimum required plank depth
    # ----------------------------------------------------------------------------
        
    # Define the panel thickness that is to be used
    CLT_Type = ['3 PLY', '5 PLY', '7 PLY', '9 PLY' ]  #Suggestion keeping the same name but adding a suffix i.e s or _list
    
    Ks = 4.8    # Fixed-Pinned Shear Shear Adjustment Factor
    plank_DF['EIapp_0'] = plank_DF['Eieff_0'] / (1 + Ks * plank_DF.Eieff_0 / plank_DF.Gaeff_0 / (s_beam * 12) ** 2)   # *** Dot notation can break if there is a space in the variable, suggestion to use the [] to directly access the attribute. Upper case vs lower case. Upper case - global variables 
    Kcr = 2     # Long-term creep factor for CLT
    
    # Iterate through CLT Types to find suitable panel
    suitable_plank = None
    for clt_type in CLT_Type:
        # Filter the DataFrame for the current CLT type
        plank_data = plank_DF[plank_DF['CLT Size'] == clt_type].iloc[0]
        
        # Create an instance of CLTSize for the current plank type
        typ_plank = CLTSize(plank_data, CM, Ct)
        DL = typ_plank.wt + t_topping/12*lambda_conc + super_DL  # Self Weight of Plank + Topping Slab + Super Imposed Deadload 
        w_u = DL + LL     
        
        # Panel Demands - Assumes a minimum of two bay span
        Mu_panel = w_u * s_beam ** 2 / 8
        Vu_panel = 5/8 * w_u * s_beam
        Def_panel_Live = (LL)/12 * (s_beam * 12) ** 4 / (185 * typ_plank.EIapp_0)
        Def_panel_LT = (DL + LL)/12 * (s_beam * 12) ** 4 / (185 * typ_plank.EIapp_0)
        
        # Calculate the plank's moment and shear capacities
        plank_moment = flexural_capacity_plank_ASD(typ_plank)
        plank_shear = shear_capacity_plank_ASD(typ_plank)  # !!!!!!!!!!!!!!!!!!!!!!!! Does the 2/3 factor need to be applied for the shear capacity?
        
        # Check if the current plank meets both moment and shear demands
        if plank_moment > Mu_panel and plank_shear > Vu_panel and (s_beam * 12) / Def_panel_Live > Def_LL_lim and (s_beam * 12) / Def_panel_LT > Def_LT_lim:
            suitable_plank = typ_plank
            break
    
    if suitable_plank is not None:
        print(f"Suitable CLT type found: {clt_type}")
    else:
        print("No suitable CLT type found that meets both the moment and shear demands.")
        
    # ----------------------------------------------------------------------------
    # Equations for minimumm required member sizes 
    # ----------------------------------------------------------------------------
    # Function to determine minimum beam size - Considers the Length Equal to W (width)
    
    
        
    # ---------Member Loading Based on Bay Layout Type-------------------------------------------------------------------
    # Bay Type 1 - Beam spans the Width, Girder Spans the Length, Columns Only Located at the corners
    #   Panel is determined based on the beam spacing
    # -------------------------------------------------------------------------------------------------------------------
    
    if bay_type == 1:
        # Beam Demands
        
        #Beam Live Load Reduction
        AT_beam = s_beam * W
        KLL = 2
        LLr_beam = LL * max(min((0.25 + 15/(KLL*AT_beam)**0.5), 1),0.5)
        
        #Beam Design Load
        w_u = DL + LLr_beam
        w_beam_DL = DL * s_beam
        w_beam_LL = LLr_beam * s_beam
        
        
        Mu_beam = (1.25 * s_beam * w_u) * W ** 2 / 8
        Vu_beam = (1.25 * s_beam * w_u) * W / 2
        
        # Use functions to size members based on demands
        beam_minV, beam_mind = min_beam_size(Vu_beam, Mu_beam, beam_24F_V4, glb_sections, w_beam_DL, w_beam_LL, CD=1.0, CL=1.0, Cfu=1.0, Cc=1.0, CI=1.0)
        
        #Girder Demands
        # Assumes two beam frames directly into the column, remaining beams frame into the girder. 
        # Assumes this is an interior span bay, typical
        
        AT_girder = L * W
        KLL = 2
        LLr_girder = LL * max(min((0.25 + 15/(KLL*AT_girder)**0.5), 1),0.5)
        
        #Girder Design Load
        w_u = (DL + beam_minV.A * gamma_GLB /144 / s_beam) + LLr_girder
        w_DL = (DL + beam_minV.A * gamma_GLB /144 / s_beam)
        w_LL = LLr_girder
        
        Vu_point = (s_beam * w_u) * W / 2
        V_DL_point = (s_beam * w_DL) * W / 2
        V_LL_point = (s_beam * w_LL) * W / 2
        
        # Define the beam length and point loads
        beam_length = L  # Beam length is the Girder Span
        ## Set Point Loads up to be iterative depending on the number of beam spacing
        load_magnitude = 2 * Vu_point
        load_LT = 2 * (0.5 * V_DL_point + V_LL_point)
        load_LL = 2 * V_LL_point
        point_loads = [{"position": s_beam * (i + 1), "magnitude": load_magnitude} for i in range(n_beam-1)]
        point_loads_LT = [{"position": s_beam * (i + 1), "magnitude": load_LT} for i in range(n_beam-1)]
        point_loads_LL = [{"position": s_beam * (i + 1), "magnitude": load_LL} for i in range(n_beam-1)]
        
        # Function to calculate reactions at supports for a simply supported beam
        def calculate_support_reactions(beam_length, point_loads):
            total_load = sum(load["magnitude"] for load in point_loads)
            moment_about_A = sum(load["magnitude"] * load["position"] for load in point_loads)
            reaction_B = moment_about_A / beam_length
            reaction_A = total_load - reaction_B
            return reaction_A, reaction_B
        
        # Function to calculate shear force at a given position
        def calculate_shear_force(position, point_loads, reaction_A):
            shear_force = reaction_A
            for load in point_loads:
                if load["position"] <= position:
                    shear_force -= load["magnitude"]
            return shear_force
        
        # Function to calculate bending moment at a given position
        def calculate_bending_moment(position, point_loads, reaction_A):
            bending_moment = reaction_A * position
            for load in point_loads:
                if load["position"] <= position:
                    bending_moment -= load["magnitude"] * (position - load["position"])
            return bending_moment
        
        # Define positions along the beam to evaluate shear and moment
        num_positions = 100
        positions = np.linspace(0, beam_length, num_positions)
        
        # Calculate shear force and bending moment at each position
        reaction_A, reaction_B = calculate_support_reactions(beam_length, point_loads)
        shear_forces = [calculate_shear_force(pos, point_loads, reaction_A) for pos in positions]
        bending_moments = [calculate_bending_moment(pos, point_loads, reaction_A) for pos in positions]
        
        # Find the maximum shear force and bending moment
        Vu_girder = max(shear_forces)
        Mu_girder = max(bending_moments)
        
        # Column Live Load Reduction
        AT_col = W * L
        KLL = 4
        LLr_col = LL * max(min((0.25 + 15/(KLL*AT_beam)**0.5), 1),0.5)
        
        # Column Design Load
        w_u = DL + LLr_col
        
        # Maximum Column Demands
        Pu_col =  w_u * W * L
                  
        girder_minV, girder_mind = min_girder_size(Vu_girder, Mu_girder, girder_24F_V4, glb_sections, DBL_glb_sections, point_loads_LT, point_loads_LL, DBL_Grider, CD=1.0, CL=1.0, Cfu=1.0, Cc=1.0, CI=1.0)
        
        # Check if either beam or girder results are zeros
        if (beam_minV.eq(0).all() or beam_mind.eq(0).all() or girder_minV.eq(0).all() or girder_mind.eq(0).all()):
            # Set all results to zero
            bay_volume_minV = 0
            bay_volume_mind = 0
            beam_volume_minV = 0
            beam_volume_mind = 0
            girder_volume_minV = 0
            girder_volume_mind = 0
            plank_volume = 0
            column_volume = 0
            n_girder = 0
            n_column = 0
            
        else:
        
            le = H * 12 - min(beam_mind.d , girder_mind.d)           #Effective Column Length
            
            column_size = min_column_size(Pu_col, le, col_lay3, glb_sections, beam_minV.b, girder_minV.b)
            
            #Bay Results
            bay_volume_minV = (W - girder_minV.b/12) * beam_minV.A / 144 * (n_beam + 1) + (L - column_size.b / 12) * girder_minV.A / 144 * 2 + suitable_plank.depth / 12 * W * L + column_size.A / 144 * (H -suitable_plank.depth / 12) * 4
            bay_volume_mind = (W - girder_mind.b/12) * beam_mind.A / 144 * (n_beam + 1) + (L - column_size.b / 12) * girder_mind.A / 144 * 2 + suitable_plank.depth / 12 * W * L + column_size.A / 144 * (H -suitable_plank.depth / 12) * 4
            beam_volume_minV = (W - girder_minV.b/12) * beam_minV.A / 144 * (n_beam + 1)
            beam_volume_mind = (W - girder_mind.b/12) * beam_mind.A / 144 * (n_beam + 1)
            girder_volume_minV = (L - column_size.b / 12) * girder_minV.A / 144 * 2
            girder_volume_mind = (L - column_size.b / 12) * girder_mind.A / 144 * 2
            plank_volume = suitable_plank.depth / 12 * W * L
            column_volume = column_size.A / 144 * H * 4
            n_girder = 2
            n_column = 4

    # ---------Member Loading Based on Bay Layout Type-------------------------------------------------------------------
    # Bay Type 2 - Beam spans the Width, Drop the Girder, Columns Only Located at the corners
    #   Panel is determined based on the beam spacing
    # -------------------------------------------------------------------------------------------------------------------

    if bay_type == 2:
        
        #Beam Live Load Reduction
        AT_beam = s_beam * W
        KLL = 2
        LLr_beam = LL * max(min((0.25 + 15/(KLL*AT_beam)**0.5), 1),0.5)
         
        #Beam Design Load
        w_beam_DL = DL * s_beam
        w_beam_LL = LLr_beam * s_beam
        w_u = DL + LLr_beam
        
        # Beam Demands
        Mu_beam = (1.25 * s_beam * w_u) * W ** 2 / 8
        Vu_beam = (1.25 * s_beam * w_u) * W / 2
        
        # Column Live Load Reduction
        AT_col = s_beam * W
        KLL = 4
        LLr_col = LL * max(min((0.25 + 15/(KLL*AT_beam)**0.5), 1),0.5)
        
        # Column Design Load
        w_u = DL + LLr_col
        
        # Maximum Column Demands
        Pu_col = (s_beam * w_u) * W
        
        # Use functions to size members based on demands
        beam_minV, beam_mind = min_beam_size(Vu_beam, Mu_beam, beam_24F_V4, glb_sections, w_beam_DL, w_beam_LL, CD=1.0, CL=1.0, Cfu=1.0, Cc=1.0, CI=1.0)
        
        # Create a girder_size with all zeros
        girder_minV = pd.Series(0, index=glb_sections.columns)
        girder_mind = pd.Series(0, index=glb_sections.columns)
        
        le = H * 12                                             #Effective Column Length
        
        column_size = min_column_size(Pu_col, le, col_lay3, glb_sections, beam_minV.b)
        
        #Bay Results
        bay_volume_minV = (W - column_size.b/12) * beam_minV.A / 144 * (n_beam + 1) + suitable_plank.depth / 12 * W * L + column_size.A / 144 * (H -suitable_plank.depth / 12) * 4
        bay_volume_mind = (W - column_size.b/12) * beam_mind.A / 144 * (n_beam + 1) + suitable_plank.depth / 12 * W * L + column_size.A / 144 * (H -suitable_plank.depth / 12) * 4
        
        beam_volume_minV = W * beam_minV.A / 144 * (n_beam + 1)
        girder_volume_minV = 0
        plank_volume_minV = suitable_plank.depth / 12 * W * L
        column_volume_minV = column_size.A / 144 * H * (n_beam + 1) * 2
        n_girder = 0
        n_column = ((n_beam + 1) * 2)         

    # ----------------------------------------------------------------------------
    # Summary or bay parameter results
    # ----------------------------------------------------------------------------
    # print(f"Bay width: {W} ft")
    # print(f"Bay length: {L} ft")
    # print(f"Beam Bays: {n_beam}")
    # print(f"Total Dead Load: {DL} psf")
    # print(f"Total Live Load: {LL} psf")
    # print(f"Beam size: {beam_size.b} x {beam_size.d}")
    # print(f"Beam connection demand: {beam_size.V_total} lbf")
    # print(f"Girder size: {girder_size.b} x {girder_size.d}")
    # print(f"Girder connection demand: {girder_size.V_total} lbf")
        
    volume_p_SF_minV = bay_volume_minV / (W * L)
    volume_p_SF_mind = bay_volume_mind / (W * L)
    # print(f"Bay volume / square foot: {volume_p_SF} ft3 / ft2")
    
    # Collect the results for this row
    result = {
        'W': W,
        'L': L,
        'Beam Bays': n_beam,
        's_beam': s_beam,
        'CLT Type': clt_type,
        'Bay Type': bay_type,
        'Total Dead Load': DL,
        'Total Live Load': LL,
        # Min Volume Results
        'Beam_Min_V': f"{beam_minV.b} x {beam_minV.d}",
        'Beam Connection Demand_MinV': beam_minV.V_total,
        'Girder_Min_V': f"({girder_minV.n_section}) - {girder_minV.b/girder_minV.n_section} x {girder_minV.d}",
        'Girder Connection Demand_minV': girder_minV.V_total,
        'volume_per_sqft_MinV': volume_p_SF_minV,
        'Bay Volume_Min V': bay_volume_minV,
        'Beam Volume_Min V': beam_volume_minV,
        'Girder Volume_Min V': girder_volume_minV,
        # Min Depth Results
        'Beam_Min_D': f"{beam_mind.b} x {beam_mind.d}",
        'Beam Connection Demand_Mind': beam_mind.V_total,
        'Girder_Min_D': f"({girder_mind.n_section}) - {girder_mind.b/girder_mind.n_section} x {girder_mind.d}",
        'Girder Connection Demand': girder_mind.V_total,
        'volume_per_sqft_Mind': volume_p_SF_mind,
        'Bay Volume_Min d': bay_volume_mind,
        'Beam Volume_Min d': beam_volume_mind,
        'Girder Volume_Min d': girder_volume_mind,
        # Uniform Variables
        'Column_Size': f"{column_size.b} x {column_size.d}",
        'Panel Volume': plank_volume,
        'Column Volume': column_volume,
        '# Beam Connections': (n_beam+1) * 2,
        '# Girder Connections': (n_girder) * 2,
        '# Column Connections': n_column * 2
    }   
    
    results.append(result)
    
# Convert results to DataFrame and save to Excel
results_df = pd.DataFrame(results)
results_df.to_excel(results_xlsx, index=False)

