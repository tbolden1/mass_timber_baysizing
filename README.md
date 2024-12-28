# mass_timber_baysizing
A python/excel based tool to explore various mass timber bay designs



This code sizes CLT Panels, beams, girders, and columns for various bay sizes.  There are two bay types that can be considered with in the code. Bay Type 1 is a system where beams span the ‘W’ length to girders which span the ‘L’ length to columns. Bay Type 2 is a system where the girders are dropped and the beams span the ‘W’ length directly to columns. The CLT panels span the ‘L’ length between beams for both bay types.  

Steps to utilize code: (_Bay_Geometry_241104.py_)
1.	The user fills out an array of bays to be sized in an excel spreadsheet that is then read by the code. _(Input_Bay_Matrix.xlsx)_
2.	The user fills out the material properties, loading criteria and deflection limits to be considered. _(Code lines 46-65)_
3.	Run the code to produce a matrix of results parameters for each bay size input by the user. _(Results_Matrix_Girder.xlsx)_
4.	This matrix can then be utilized to visualize the results of all of bay parameters input by the user. _(Plotting_Function.py)_

