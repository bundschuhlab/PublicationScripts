#!/usr/bin/env python

#Copyright (C) <2021>  <The Ohio State University>       

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

#This code was written by Jackson Killian

import argparse
import numpy as np
import matplotlib.pyplot as plt

def parse_degs(open_file):

	d = {}
	
	#consume header
	open_file.readline()
	
	for line in open_file:
		line = line.strip().split('\t')
		gene = line[0]
		fchange = 0
		try:
			fchange = float(line[2])
		except:
			fchange = 0
		qval = 0
		try:
			qval = float(line[6])
		except:
			qval = 1
			
		if gene not in d and fchange != 0:
			d[gene] = {'dir':fchange < 0,'q':qval}
		
	return d
	
def parse_dmrs(open_file):
	
	d = {}
	
	#read header
	head = open_file.readline().split('\t')
	
	gcol = -1
	qcol = -1
	acol = -1
	bcol = -1
	for i,item in enumerate(head):
		if item == "genesym": gcol = i
		if item == "qval": qcol = i
		if item == "avgGrpA": acol = i
		if item == "avgGrpB": bcol = i
	if gcol == -1: raise ValueError("Could not find genesym column.")
	if qcol == -1: raise ValueError("Could not find qval column.")
	if acol == -1: raise ValueError("Could not find avgGrpA column.")
	if bcol == -1: raise ValueError("Could not find avgGrpB column.")
	
	
	for line in open_file:
		line = line.strip().split('\t')
		genestr = line[gcol]
		qval = float(line[qcol])
		a = float(line[acol])
		b = float(line[bcol])
		
		genes = genestr.split(';')
		
		for gene in genes:
			if gene not in d:
				d[gene] = {'dir': b > a,'q':qval}
	
	return d	
		
	

parser = argparse.ArgumentParser(description = "")
parser.add_argument("dmrs", type = str, help = "MethMAGE output")
parser.add_argument("degs", type = str, help = "DEGs")
parser.add_argument("-c", "--cutoff", type = float, default = 0.05, help = "cutoff")
parser.add_argument("-t", "--title", type = str, default = "DMRs vs DEGs", help = "cutoff")
parser.add_argument("-o", "--output", type = str, default = "dmrs_v_degs.png", help = "cutoff")
parser.add_argument("-sfn", "--significant_same_filename", type = str, default = "matched_dmrs_degs_same.txt", help = "filename for list of dmrs and degs that are significant and have the same direction")
parser.add_argument("-ofn", "--significant_opposite_filename", type = str, default = "matched_dmrs_degs_opposite.txt", help = "filename for list of dmrs and degs that are significant and have the opposite direction")
parser.add_argument("-plots", "--make_plots", action = "store_true", help = "make the plots")
args = parser.parse_args()

dmrs = parse_dmrs(open(args.dmrs))
degs = parse_degs(open(args.degs))

print args.title
print len(dmrs)
print len(degs)

print len(set(dmrs.keys()) & set(degs.keys()))

sf = open(args.significant_same_filename,'w')
of = open(args.significant_opposite_filename,'w')


# for DMRS
# [sig-same, sig-different, nonsig-same, nonsig-different]
dmr_bars = [0 for i in range(4)]
for gene in dmrs:
	m_dir = dmrs[gene]['dir']
	m_qval = dmrs[gene]['q']
	if m_qval <= args.cutoff and gene in degs:
		g_dir = degs[gene]['dir']
		g_qval = degs[gene]['q']
		
		#sig-same
		if   g_qval <= args.cutoff and m_dir == g_dir:
			dmr_bars[0]+=1
			print >> sf, gene
		#sig-different
		elif g_qval <= args.cutoff and m_dir != g_dir:
			dmr_bars[1]+=1
			print >> of, gene
		#nonsig-same
		elif g_qval > args.cutoff and m_dir == g_dir: dmr_bars[2]+=1
		#nonsig-different
		elif g_qval > args.cutoff and m_dir != g_dir: dmr_bars[3]+=1

sf.close()
of.close()

#for DEGs
# [sig-same, sig-different, nonsig-same, nonsig-different]
deg_bars = [0 for i in range(4)]
for gene in degs:
	g_dir = degs[gene]['dir']
	g_qval = degs[gene]['q']
	if g_qval <= args.cutoff and gene in dmrs:
		m_dir = dmrs[gene]['dir']
		m_qval = dmrs[gene]['q']
		
		#sig-same
		if   m_qval <= args.cutoff and m_dir == g_dir: deg_bars[0]+=1
		#sig-different
		elif m_qval <= args.cutoff and m_dir != g_dir: deg_bars[1]+=1
		#nonsig-same
		elif m_qval > args.cutoff and m_dir == g_dir: deg_bars[2]+=1
		#nonsig-different
		elif m_qval > args.cutoff and m_dir != g_dir: deg_bars[3]+=1

print dmr_bars
print deg_bars

if args.make_plots:
	dmr_sigs = [dmr_bars[0],dmr_bars[1],0]
	dmr_nonsigs = [dmr_bars[2],dmr_bars[3],0]
	deg_sigs = [deg_bars[0],deg_bars[1],0]
	deg_nonsigs = [deg_bars[2],deg_bars[3],0]

	index = np.arange(3)
	bar_width = 0.1

	opacity = 0.6
	lpacity = opacity - 0.2
	error_config = {'ecolor': '0.3'}

	rects1 = plt.bar(index, dmr_sigs, bar_width,
		             alpha=opacity,
		             color='b',
		             label='DMR/DEG')
	rects2 = plt.bar(index, dmr_nonsigs, bar_width,
		             alpha=lpacity,
		             color='b',
		             bottom=dmr_sigs,
		             label='DMR/nonDEG')
		             

	rects3 = plt.bar(index + 2*bar_width, deg_sigs, bar_width,
		             alpha=opacity,
		             color='r',
		             label='DEG/DMR')
	rects4 = plt.bar(index + 2*bar_width, deg_nonsigs, bar_width,
		             alpha=lpacity,
		             color='r',
		             bottom=deg_sigs,
		             label='DEG/nonDMRs')
	plt.ylim([0,575])
	plt.xlabel('Group')
	plt.ylabel('Count')
	plt.title(args.title)
	plt.xticks(index + bar_width, ( 'Same','Opposite'))
	plt.legend()

	plt.tight_layout()
	plt.savefig(args.output)
	plt.show()
