#Copyright (C) <2019>  <The Ohio State University>

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
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

outfilename = sys.argv[1]
outfile = open(outfilename, 'w')

infilename_1 = sys.argv[2] #file of H2B up genes
infile_1 = open(infilename_1, 'r')                                                                   

infilename_2 = sys.argv[3] #file of H2B down genes                                                                                                                                      
infile_2 = open(infilename_2, 'r')

infilename_3 = sys.argv[4] #file of H2B even genes                                                                                                                           
infile_3 = open(infilename_3, 'r')

infilename_4 = sys.argv[5] #file of H4 even genes                                                                                                                 
infile_4 = open(infilename_4, 'r')

upgenes = list()
downgenes = list()
evengenes = list()
H4evengenes = list()
chr_list = list()

U_a_list = [0]*147
U_t_list = [0]*147
U_c_list = [0]*147
U_g_list = [0]*147

D_a_list = [0]*147
D_t_list = [0]*147
D_c_list = [0]*147
D_g_list = [0]*147

E_a_list = [0]*147
E_t_list = [0]*147
E_c_list = [0]*147
E_g_list = [0]*147

#read gene data into memory
for x in range(16) :          
    
    infilename = #path to each genome
    infile = open(infilename, 'r')
    chr_list.append(infile.readline())
    infile.close

for line1 in infile_1 :
    
    upgenes.append(line1.strip())

for line2 in infile_2 :

    downgenes.append(line2.strip())

for line3 in infile_3 :

    evengenes.append(line3.strip())

for line4 in infile_4 :

    H4evengenes.append(line4.strip())

for x in range(1,4) :

    infile_5 = #path to +1, +2, and +3 nucleosome files
    infile_5.readline()

    for line in infile_5 :

        fields = line.split(',')
        chrom  = fields[0]
        chrom = int(chrom[3:])
        gene = fields[1]
        nuc = gene + ' ' + str(x)
        dyad = int(fields[3])
        genome = chr_list[chrom-1]
        temp = genome[dyad-73:dyad+74]
        tempnew = ''
        if(gene[-1] == 'C' or gene[-3] == 'C') :
            temp = temp[::-1]
            for z in range(len(temp)) :
                if(temp[z] == 'A') :
                    tempnew += 'T'
                elif(temp[z] == 'T') :
                    tempnew += 'A'
                elif(temp[z] =='G') :
                    tempnew += 'C'
                elif(temp[z] =='C') :
                    tempnew += 'G'
            temp = tempnew

        if(nuc in upgenes and nuc in H4evengenes) :

            for p in range(len(temp)) :
                if(temp[p] == 'A') :
                    U_a_list[p] += 1
                if(temp[p] == 'T') :
                    U_t_list[p] += 1
                if(temp[p] == 'C') :
                    U_c_list[p] += 1
                if(temp[p] == 'G') :
                    U_g_list[p] += 1

        if(nuc in downgenes and nuc in H4evengenes) :

            for p in range(len(temp)) :
                if(temp[p] == 'A') :
                    D_a_list[p] += 1
                if(temp[p] == 'T') :
                    D_t_list[p] += 1
                if(temp[p] == 'C') :
                    D_c_list[p] += 1
                if(temp[p] == 'G') :
                    D_g_list[p] += 1

        if(nuc in evengenes and nuc in H4evengenes) :

            for p in range(len(temp)) :
                if(temp[p] == 'A') :
                    E_a_list[p] += 1
                if(temp[p] == 'T') :
                    E_t_list[p] += 1
                if(temp[p] == 'C') :
                    E_c_list[p] += 1
                if(temp[p] == 'G') :
                    E_g_list[p] += 1

    infile_5.close()

plot_list = [U_a_list,U_t_list,U_c_list,U_g_list,D_a_list,D_t_list,D_c_list,D_g_list,E_a_list,E_t_list,E_c_list,E_g_list]                                                           

totalU = U_a_list[0]+U_t_list[0]+U_c_list[0]+U_g_list[0]
totalD = D_a_list[0]+D_t_list[0]+D_c_list[0]+D_g_list[0]
totalE = E_a_list[0]+E_t_list[0]+E_c_list[0]+E_g_list[0]

for y in range(0,4) :
    plot_list[y] = [float(x)/totalU for x in plot_list[y]]

for y in range(4,8) :
    plot_list[y] = [float(x)/totalD for x in plot_list[y]]

for y in range(8,12) :
    plot_list[y] = [float(x)/totalE for x in plot_list[y]]

for i in range(0,4,1) :                                                                                                                                                             
    for j in range(3,len(plot_list[i])-9,10) :                                                                                                                                               
        ttest_UD = stats.ttest_1samp(np.divide(plot_list[i][j:j+10],plot_list[i+4][j:j+10]),1)                                                                                             
        outfile.write(str(ttest_UD[1]) + '\n')                                                                                                                                                
    outfile.write('\n')                                                                                                                                                                   

for i in range(0,4,1) :                                                                                                                                                             
    for j in range(3,len(plot_list[i])-9,10) :                                                                                                                                               
        ttest_UE = stats.ttest_1samp(np.divide(plot_list[i][j:j+10],plot_list[i+8][j:j+10]),1)                                                                                             
        outfile.write(str(ttest_UE[1]) + '\n')                                                                                                                                            
    outfile.write('\n')                                                                                                                                                                    

for i in range(4,8,1) :                                                                                                                                                             
    for j in range(3,len(plot_list[i])-9,10) :                                                                                                                                               
        ttest_DE = stats.ttest_1samp(np.divide(plot_list[i][j:j+10],plot_list[i+4][j:j+10]),1)                                                                                         
        outfile.write(str(ttest_DE[1]) + '\n')                                                                                                             

outfile.close
