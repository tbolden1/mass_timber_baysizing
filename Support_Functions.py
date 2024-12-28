# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 13:29:13 2024

@author: tbolden
"""

import numpy as np

def deflection_due_to_point_load(L, E, I, P, a, x):
    """
    Calculate deflection at a point 'x' due to a point load 'P' at distance 'a' from the left support.
    """
    b = L - a
    if x <= a:
        delta = (P * b * x) / (6 * L * E * I) * (L**2 - b**2 - x**2)
    else:
        delta = (P * a * (L - x)) / (6 * L * E * I) * (L**2 - a**2 - (L - x)**2)
    return delta

def calculate_total_deflection(L, E, I, point_loads, num_points=100):
    """
    Calculate the deflections across the length of the girder due to multiple point loads.
    """
    x_positions = np.linspace(0, L, num_points)
    total_deflections = np.zeros(num_points)
    
    # Iterate over each point load and calculate deflections due to it
    for load in point_loads:
        P = load['magnitude']
        a = load['position']
        deflections = [deflection_due_to_point_load(L, E, I, P, a, x) for x in x_positions]
        total_deflections += deflections  # Sum the deflections from each point load
    
    # Max deflection
    max_deflection = np.max(total_deflections)
    return max_deflection