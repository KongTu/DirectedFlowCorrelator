# DirectedFlowCorrelator
The latest development of the code is in Branch "CodeWithOnlyOneAngle"

- The source code can be found in plugin, parameters are defined in python folder. 
- The CMSSW run config and crab config are under test folder.
- The outputs are stored in the format of histograms. 
- Macros folder has all needed macros to analyze the histograms, e.g., obtaining the correlation, mass, etc.

#Inportant macros:

- macros/massfitvn_combine.C
- macros/massfitvn_combine_3fit.C

These two can simultaneously fit mass and Vn distribution, and extract the signal Vn. The input files are generated from macros/D0_v1VSmass_odd.C (or with the mixed extension).
