# Hydro_Seper
This repository contains the Matlab program for the Characteristic Point Method (CPM) and some other codes related to the topic of hydrograph separation and properties. CPM is a method for automatic separation of hydrograph and extraction of flood events from long-term time series reported in Mei & Anagnostou (2015, A Hydrograph Separation Method Based on Information from Rainfall and Runoff Records). A user manual of the program and demonstration of implementation are also provided.

Details of CPM:
 1. The "Demo_data.mat" contains inputs of the CPM program and the "Results.mat" contains results produced by the program.
 2. The "Demo.m" file is an example of how to run the CPM program.
 3. The "RCK", "CPM.m", "CPM_FE.m", "CPM_peak.m", and "CPM_RE.m" are components of the CPM program.
 4. The "CPM User Manual_MMMYYYY.pdf" is the user manual of the program.
 5. Please see the user manual for instructions on how to set up the CPM program. The user manual also contains the instructions for
    running the Demo code.
 6. The program may subject to update ocassionally. Please check for the new version on our website.

Other related Matlab codes:
 1. The "hydro_pro.m" file is a code for calculation of the rainfall-runoff event properties.
 2. The "SMM.m" file is a code to implement the United Kindom Institute of Hydrology Smooth Minima Method (SMM) baseflow separation method by Gustard
    et al. (1992, Low Flow Estimation in the United Kingdom).
 3. The "RDF.m" file is a code to implement the recursive digital filter (RDF) baseflow separation method by Eckhardt (2005, How to
    construct recursive digital filters for baseflow separation).
 4. The "BRM.m" file is a code to implement the Bump and Rise Method (BRM) baseflow separation method by Steward (2015, Promising new baseflow 
    separation and recession analysis methods applied to streamflow at Glendhu Catchment, New Zealand).
