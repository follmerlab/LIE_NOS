# LIE_NOS
Analysis scripts used for "Application of the Linear Interaction Energy Method to Nitric Oxide Synthase Structure Based Inhibitor Design"

The 'runGBLIE.ipynb' is used to perform an MMGBSA calculation on a user supplied trajectory stripped of explicit waters. The functions are imported from 'gbLIE.py' and use the 'mmpbsa.in' template file when it calls the 'mm_pbsa.pl' script. It can be run here with the example folder for ligand JI4.

Application_of_LIE_NOS.ipynb notebook is for fitting alpha,beta,and gamma reading EEL and VDW values from the 'mmgbsa_vals_4fit.csv' and 'SI_NOS_InhibitorTable_Final.xlsx' spreadsheets.


