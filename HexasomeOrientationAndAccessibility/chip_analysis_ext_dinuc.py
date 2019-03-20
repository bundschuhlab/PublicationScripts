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

legend_properties = {'weight':'bold', 'size':8}

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

U_at_list = [0]*146
U_cg_list = [0]*146
U_ag_list = [0]*146
U_tg_list = [0]*146
U_ct_list = [0]*146
U_ac_list = [0]*146
U_aa_list = [0]*146
U_tt_list = [0]*146
U_cc_list = [0]*146
U_gg_list = [0]*146
U_ta_list = [0]*146
U_gc_list = [0]*146
U_ga_list = [0]*146
U_gt_list = [0]*146
U_tc_list = [0]*146
U_ca_list = [0]*146

D_at_list = [0]*146
D_cg_list = [0]*146
D_ag_list = [0]*146
D_tg_list = [0]*146
D_ct_list = [0]*146
D_ac_list = [0]*146
D_aa_list = [0]*146
D_tt_list = [0]*146
D_cc_list = [0]*146
D_gg_list = [0]*146
D_ta_list = [0]*146
D_gc_list = [0]*146
D_ga_list = [0]*146
D_gt_list = [0]*146
D_tc_list = [0]*146
D_ca_list = [0]*146

E_at_list = [0]*146
E_cg_list = [0]*146
E_ag_list = [0]*146
E_tg_list = [0]*146
E_ct_list = [0]*146
E_ac_list = [0]*146
E_aa_list = [0]*146
E_tt_list = [0]*146
E_cc_list = [0]*146
E_gg_list = [0]*146
E_ta_list = [0]*146
E_gc_list = [0]*146
E_ga_list = [0]*146
E_gt_list = [0]*146
E_tc_list = [0]*146
E_ca_list = [0]*146

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
                
            for p in range(0,len(temp)-1,1) :

                if(temp[p] == 'A' and temp[p+1] == 'T') :
                    U_at_list[p] += 1
                if(temp[p] == 'C' and temp[p+1] == 'G') :
                    U_cg_list[p] += 1
                if(temp[p] == 'A' and temp[p+1] == 'G') :
                    U_ag_list[p] += 1
                if(temp[p] == 'T' and temp[p+1] == 'G') :
                    U_tg_list[p] += 1
                if(temp[p] == 'C' and temp[p+1] == 'T') :
                    U_ct_list[p] += 1
                if(temp[p] == 'A' and temp[p+1] == 'C') :
                    U_ac_list[p] += 1
                if(temp[p] == 'A' and temp[p+1] == 'A') :
                    U_aa_list[p] += 1
                if(temp[p] == 'T' and temp[p+1] == 'T') :
                    U_tt_list[p] += 1
                if(temp[p] == 'C' and temp[p+1] == 'C') :
                    U_cc_list[p] += 1
                if(temp[p] == 'G' and temp[p+1] == 'G') :
                    U_gg_list[p] += 1
                if(temp[p] == 'T' and temp[p+1] == 'A') :
                    U_ta_list[p] += 1
                if(temp[p] == 'G' and temp[p+1] == 'C') :
                    U_gc_list[p] += 1
                if(temp[p] == 'G' and temp[p+1] == 'A') :
                    U_ga_list[p] += 1
                if(temp[p] == 'G' and temp[p+1] == 'T') :
                    U_gt_list[p] += 1
                if(temp[p] == 'T' and temp[p+1] == 'C') :
                    U_tc_list[p] += 1
                if(temp[p] == 'C' and temp[p+1] == 'A') :
                    U_ca_list[p] += 1

        if(nuc in downgenes and nuc in H4evengenes) :
            
            for p in range(0,len(temp)-1,1) :

                if(temp[p] == 'A' and temp[p+1] == 'T') :
                    D_at_list[p] += 1
                if(temp[p] == 'C' and temp[p+1] == 'G') :
                    D_cg_list[p] += 1
                if(temp[p] == 'A' and temp[p+1] == 'G') :
                    D_ag_list[p] += 1
                if(temp[p] == 'T' and temp[p+1] == 'G') :
                    D_tg_list[p] += 1
                if(temp[p] == 'C' and temp[p+1] == 'T') :
                    D_ct_list[p] += 1
                if(temp[p] == 'A' and temp[p+1] == 'C') :
                    D_ac_list[p] += 1
                if(temp[p] == 'A' and temp[p+1] == 'A') :
                    D_aa_list[p] += 1
                if(temp[p] == 'T' and temp[p+1] == 'T') :
                    D_tt_list[p] += 1
                if(temp[p] == 'C' and temp[p+1] == 'C') :
                    D_cc_list[p] += 1
                if(temp[p] == 'G' and temp[p+1] == 'G') :
                    D_gg_list[p] += 1
                if(temp[p] == 'T' and temp[p+1] == 'A') :
                    D_ta_list[p] += 1
                if(temp[p] == 'G' and temp[p+1] == 'C') :
                    D_gc_list[p] += 1
                if(temp[p] == 'G' and temp[p+1] == 'A') :
                    D_ga_list[p] += 1
                if(temp[p] == 'G' and temp[p+1] == 'T') :
                    D_gt_list[p] += 1
                if(temp[p] == 'T' and temp[p+1] == 'C') :
                    D_tc_list[p] += 1
                if(temp[p] == 'C' and temp[p+1] == 'A') :
                    D_ca_list[p] += 1

        if(nuc in evengenes and nuc in H4evengenes) :

            for p in range(0,len(temp)-1,1) :
                
                if(temp[p] == 'A' and temp[p+1] == 'T') :
                    E_at_list[p] += 1
                if(temp[p] == 'C' and temp[p+1] == 'G') :
                    E_cg_list[p] += 1
                if(temp[p] == 'A' and temp[p+1] == 'G') :
                    E_ag_list[p] += 1
                if(temp[p] == 'T' and temp[p+1] == 'G') :
                    E_tg_list[p] += 1
                if(temp[p] == 'C' and temp[p+1] == 'T') :
                    E_ct_list[p] += 1
                if(temp[p] == 'A' and temp[p+1] == 'C') :
                    E_ac_list[p] += 1
                if(temp[p] == 'A' and temp[p+1] == 'A') :
                    E_aa_list[p] += 1
                if(temp[p] == 'T' and temp[p+1] == 'T') :
                    E_tt_list[p] += 1
                if(temp[p] == 'C' and temp[p+1] == 'C') :
                    E_cc_list[p] += 1
                if(temp[p] == 'G' and temp[p+1] == 'G') :
                    E_gg_list[p] += 1
                if(temp[p] == 'T' and temp[p+1] == 'A') :
                    E_ta_list[p] += 1
                if(temp[p] == 'G' and temp[p+1] == 'C') :
                    E_gc_list[p] += 1
                if(temp[p] == 'G' and temp[p+1] == 'A') :
                    E_ga_list[p] += 1
                if(temp[p] == 'G' and temp[p+1] == 'T') :
                    E_gt_list[p] += 1
                if(temp[p] == 'T' and temp[p+1] == 'C') :
                    E_tc_list[p] += 1
                if(temp[p] == 'C' and temp[p+1] == 'A') :
                    E_ca_list[p] += 1

    infile_5.close()

U_plot_list = [U_at_list,U_cg_list,U_ag_list,U_tg_list,U_ct_list,U_ac_list,U_aa_list,U_tt_list,U_cc_list,U_gg_list,U_ta_list,U_gc_list,U_ga_list,U_gt_list,U_tc_list,U_ca_list]
D_plot_list = [D_at_list,D_cg_list,D_ag_list,D_tg_list,D_ct_list,D_ac_list,D_aa_list,D_tt_list,D_cc_list,D_gg_list,D_ta_list,D_gc_list,D_ga_list,D_gt_list,D_tc_list,D_ca_list]
E_plot_list = [E_at_list,E_cg_list,E_ag_list,E_tg_list,E_ct_list,E_ac_list,E_aa_list,E_tt_list,E_cc_list,E_gg_list,E_ta_list,E_gc_list,E_ga_list,E_gt_list,E_tc_list,E_ca_list] 

totalU = U_at_list[0]+U_cg_list[0]+U_ag_list[0]+U_tg_list[0]+U_ct_list[0]+U_ac_list[0]+U_aa_list[0]+U_tt_list[0]+U_cc_list[0]+U_gg_list[0]+U_ta_list[0]+U_gc_list[0]+U_ga_list[0]+U_gt_list[0]+U_tc_list[0]+U_ca_list[0]
totalD = D_at_list[0]+D_cg_list[0]+D_ag_list[0]+D_tg_list[0]+D_ct_list[0]+D_ac_list[0]+D_aa_list[0]+D_tt_list[0]+D_cc_list[0]+D_gg_list[0]+D_ta_list[0]+D_gc_list[0]+D_ga_list[0]+D_gt_list[0]+D_tc_list[0]+D_ca_list[0]
totalE = E_at_list[0]+E_cg_list[0]+E_ag_list[0]+E_tg_list[0]+E_ct_list[0]+E_ac_list[0]+E_aa_list[0]+E_tt_list[0]+E_cc_list[0]+E_gg_list[0]+E_ta_list[0]+E_gc_list[0]+E_ga_list[0]+E_gt_list[0]+E_tc_list[0]+E_ca_list[0]

for y in range(0,16) :
    U_plot_list[y] = [float(x)/totalU for x in U_plot_list[y]]

for y in range(0,16) :
    D_plot_list[y] = [float(x)/totalD for x in D_plot_list[y]]

for y in range(0,16) :
    E_plot_list[y] = [float(x)/totalE for x in E_plot_list[y]]

for i in range(len(U_plot_list)) :
    for j in range(0,len(U_plot_list[i])-9,10) :
        ttest_UD = stats.ttest_1samp(np.divide(U_plot_list[i][j:j+10],D_plot_list[i][j:j+10]),1)
        outfile.write(str(ttest_UD[1]) + '\n')
    outfile.write('\n')

for i in range(len(U_plot_list)) :
    for j in range(0,len(U_plot_list[i])-9,10) :
        ttest_UE = stats.ttest_1samp(np.divide(U_plot_list[i][j:j+10],E_plot_list[i][j:j+10]),1)
        outfile.write(str(ttest_UE[1]) + '\n')
    outfile.write('\n')

for i in range(len(U_plot_list)) :
    for j in range(0,len(U_plot_list[i])-9,10) :
        ttest_DE = stats.ttest_1samp(np.divide(D_plot_list[i][j:j+10],E_plot_list[i][j:j+10]),1)
        outfile.write(str(ttest_DE[1]) + '\n')
    outfile.write('\n')

for line in range(len(U_plot_list)) :
    newlist = []
    for element in range(0,len(U_plot_list[line])-9,10) :
        newlist.append(sum(U_plot_list[line][element:element+10])/10)
    U_plot_list[line] = newlist[:]

for line in range(len(D_plot_list)) :
    newlist = []
    for element in range(0,len(D_plot_list[line])-9,10) :
        newlist.append(sum(D_plot_list[line][element:element+10])/10)
    D_plot_list[line] = newlist[:]

for line in range(len(E_plot_list)) :
    newlist = []
    for element in range(0,len(E_plot_list[line])-9,10) :
        newlist.append(sum(E_plot_list[line][element:element+10])/10)
    E_plot_list[line] = newlist[:]

outfile.close
