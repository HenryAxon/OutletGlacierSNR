# OutletGlacierSNR
Analysis scripts used to compute signal to noise and amplification via UO HPC resources.

# snr_new.m
This file contains the script used to compute signal to noise ratio for a large aleatoric ensemble using a linear regression to find the change in the ensemble mean of a parameter.

# snr_means.m 
This script computes SNR similarly to snr_new.m, but instead of a linear regression simply takes the differnece between the ensemble mean at time t and the ensemble mean at time t+1 to get the noise component. This method is NOT as effective but was used as a possible alternative, despite a bit of hope that it might be a simpler, faster way to compute the same things are the linear regression based method.

# gates_fixedX.m
This script computes velocity, flux and thickness in time at locations given to it.

# SSA_SNR_Control.m
This script reads in the outputs of the gate function, and then computes the SNR via the standard method described in snr_new.m, returns ice dynamical SNR values and its components. 

# SSA_SRN_Means_Control.m 
The exact same as the last script, but for the alternative method of SNR computation.

# PlottingSNRBasic.m
The plotting script for SNR and its components. Amplification factor plots were created in a seperate file, but to make those plots and do that analysis is simply to divide by SNR of the forcing. 
