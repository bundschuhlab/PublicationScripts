#!/usr/bin/env python

#This program calculates the partition function and probability 
#to be in an open state for the composite AuNP-DNA nano-hinge 
#described in the manuscript 
#"A quantitative model for a nanoscale switch accurately predicts thermal actuation behavior" 
#by Kyle Crocker, Joshua Johnson, Wolfgang Pfeifer, Carlos Castro, and Ralf Bundschuh.


#Copyright (C) <2020>  <The Ohio State University>       

#This program is free software: you can redistribute it and/or modify                              
#it under the terms of the GNU General Public License as published by 
#the Free Software Foundation, either version 3 of the License, or    
#(at your option) any later version.                                                                                       
#This program is distributed in the hope that it will be useful, 
#but WITHOUT ANY WARRANTY; without even the implied warranty of           
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
#GNU General Public License for more details.                                                                             

#You should have received a copy of the GNU General Public License 
#along with this program.  If not, see <https://www.gnu.org/licenses/>.

import sys
import numpy as np
import scipy.special as spec
from ast import literal_eval

##define constants##
#Boltzmann constant
k = 1.38e-23 #J/K                                                                       
k = k/1000.0 #kJ/K                                              
k = k*6.02e23 #kJ/mol K

#number of bases on ssDNAs coating the AuNP                                        
T_bases = 23
num_Ts = T_bases - 1

#SantaLucia terminal base pairing terms, see  https://doi.org/10.1073/pnas.95.4.1460
H_term = 9.6
S_term = 0.0172

####
##partition function calculation as described in Crocker et al.##
####

def Z_closed_states(num_types, N_c_max, N_c, n, Z_S, i, n_sum_old, len_T):  
    Z_closed_states_sum = np.zeros(len_T)
    while n[i] <= N_c_max-n_sum_old:
        if i == num_types - 1:
            Z_closed_states_sum = Z_closed_states_sum + spec.binom(N_c[i], n[i])*(Z_S[i,:]**n[i])
            n[i] = n[i] + 1
        elif i < num_types - 1:
            Z_closed_states_sum = Z_closed_states_sum + spec.binom(N_c[i], n[i])*(Z_S[i,:]**n[i])*Z_closed_states(num_types, N_c_max, N_c, n, Z_S, i + 1, n_sum_old + n[i], len_T)
            n[i] = n[i] + 1
    #reset n[i]
    n[i] = 0
    return Z_closed_states_sum

def ideal_Z_calc(N_stacks, T, H_c, H_s, S_s, S_c, N_c, N_c_max, num_types):
    norm_G_c = (H_c - T * S_c) / (k*T)
    norm_G_s = np.zeros((num_types,np.max(np.asarray(N_stacks)), len(T)))
    for l in range(num_types):
        for i in range(N_stacks[l]):
            for j in range(len(T)):
                m = i + 1
                norm_G_s[l][i][j] = (m * (H_s - T[j] * S_s)+ 2*(H_term - T[j]*S_term)) / (k*T[j])
    Z_S = np.zeros((num_types, len(T)))
    for j in range(num_types):
        for i in range(N_stacks[j]):
            m = i + 1
            Z_S[j,:] = (N_stacks[j] - m + 1)*(num_Ts - m + 1)*np.exp(-norm_G_s[j,i,:]) + Z_S[j,:]
    n = np.zeros(num_types)
    Z_c_new = Z_closed_states(num_types, N_c_max, N_c, n, Z_S, 0, 0, len(T))
    Z = Z_c_new*np.exp(-norm_G_c) + 1 #prefactors are multiplicities, 1 is by definition the partition function of the reference open state
    return Z

def main():
    print('usage is python Z_calc.py H_cl S_cl H_bp S_bp num_types [N_stacks_1,N_stacks_2,...] [N_c_1,N_c_2,...] N_c_max T_min T_max DeltaT')
    print('')
    print('')

    print('parameters are: hinge closing enthalpy; hinge closing entropy; base pairing enthalpy; base pairing entropy;  number of types of overhangs; number of stacks (N_bp - 1) in each overhang type as an array i.e. [number of stacks of type 1, number of stacks of type 2, ...]; number of overhangs for each overhang type as an array i.e. [number of overhangs of type 1, number of overhangs of type 2, ...]; maximum number of strands allowed to bind simultaneously; minimum temperature at which to calculate probability; maximum temperature at which to calculate probability; spacing of temperature values at which to calculate the probability')
    print('')
    print('')
    
    print('for example, to generate a data file that predicts the actuation curve of A_{6,6,8} from 20 C to 55 C, spaced by 0.1 C, using the best fit parameters, type: \n python Z_calc.py 0 -0.09664 -33.08 -0.09974 2 [5,7] [2,1] 2 20 55 0.1')
    print('')
    print('')

    H_c = float(sys.argv[1])
    S_c = float(sys.argv[2])
    H_s = float(sys.argv[3])
    S_s = float(sys.argv[4])
    num_types = int(sys.argv[5])
    N_stacks = np.asarray(literal_eval(sys.argv[6]))
    N_c = np.asarray(literal_eval(sys.argv[7]))
    N_c_max = float(sys.argv[8])
    T = np.arange(float(sys.argv[9]), float(sys.argv[10]), float(sys.argv[11]))
    T_K = T + 273.15 # convert to Kelvin for the sake of the calculation
    if num_types != len(N_stacks) or num_types != len(N_c):
        raise NameError('num_types should be the same as the length of N_stacks and N_c')
    
    Z = ideal_Z_calc(N_stacks, T_K, H_c, H_s, S_s, S_c, N_c, N_c_max, num_types)

    #save T, Z, and the corresponding probability that the hinge is open
    filename = ('output_H_cl='+sys.argv[1]+'_S_cl='+sys.argv[2]+'_H_bp='+sys.argv[3]+'_S_bp='+sys.argv[4]+'_num_types='+sys.argv[5]+'_N_stacks='+sys.argv[6]+'_N_c='+sys.argv[7]+'_N_c_max='+sys.argv[8]+'_T_min='+sys.argv[9]+'C_T_max='+sys.argv[10]+'C_DeltaT='+sys.argv[11]+'C.txt')
    f = open(filename, "w")
    f.write("#T          Z               p_open \n")
    for i in range(len(T)):
        f.write(str(T[i]) + '      ' + str(Z[i])+'     ' + str(1/Z[i])+'\n')
    f.close()



if __name__ == "__main__":
    main()
