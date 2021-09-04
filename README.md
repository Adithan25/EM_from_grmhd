# EM_from_grmhd
Calculate afterglow emission using structure of blastwaves from GRMHD simulations as input

This code calculates the afterglow emission from a forward shock for a structured blastwave in a 2D approximation for a jet (jet_AG) and Kilonova (KN_AG). 
The input requires three arrays of the polar angle (theta) values of the blast wave, the energy per solid angle and Lorentz factor of the blastwave at each value of theta.
The input parameters for the blastwave can be chosen as described in the preamble, the notation and descrption for these parameters follows the work of Sari et. al (1999) and the synchrotron calcualtion mainly follows this work as well, with some modifications and improvements to generalize the blastwave, see Kathirgamaraju et al. (2019).
The output is a text file which either outputs the light curve or spectrum. The light curve consists of the emitted flux (in micro Jy by default) and time. The spectrum data will consist of time and frequency of emission. 
Please email adithankathirgamaraju@gmail.com for questions or more info.
