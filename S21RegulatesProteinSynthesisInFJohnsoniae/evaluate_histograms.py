#!/usr/bin/env python

#This program calculates cumulative histograms of SD free energies for
#multiple bacterial genomes as well as a histogram of the number of
#genes with unsually strong SD sequences per organism

#Copyright (C) <2022>  <The Ohio State University>       

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


#
# import necessary libraries
#
import csv
import glob
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as mticker
import scipy.optimize as optimization
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import font_manager
#
# select font for figures
#
font_manager.findSystemFonts(fontpaths=None, fontext="ttf")
font_manager.findfont("Arial") # Test with "Special Elite" too
csfont = {'fontname':'arial'}
axis_label_size = 16 # was 24 for Bryan
title_label_size = 34 # was 34 for Bryan
tick_label_size = 14 # was 24 for Bryan
#
# linear function for fitting
#
def linearfunction(x,m,b):
    return m*x+b
#
# set parameters
#
minbin=-6
maxbin=-1
offsetbins=10
maxexcess=6
plot_histograms=True
#
# initialize the overall histogram
#
excess_genes = np.zeros((offsetbins,maxexcess+1))
#
# find all files
#
files=glob.glob('freescan_summaries/*-sd_min_tirs.csv')
#
# loop over all files
#
for file in files:
    filename=file[19:]
    firstdash=filename.find('-')
    accession=filename[:firstdash]
    #
    # read SD TIR energies from csv file
    #
    tirenergies = []
    with open('freescan_summaries/' + accession + '-sd_min_tirs.csv') as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader)
        for row in csvreader:
            tirenergies.append(float(row[1]))
    #
    # read MSD TIR energies from csv file
    #
    msdtirenergies = []
    with open('freescan_summaries/' + accession + '-msd_min_tirs.csv') as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader)
        for row in csvreader:
            msdtirenergies.append(float(row[1]))
    #
    # make a histogram of TIR energies
    #
    histogram, edges = np.histogram(tirenergies, bins=25, range=(-25,0))
    cumhist = np.cumsum(histogram)
    msdhistogram, msdedges = np.histogram(msdtirenergies, bins=25, range=(-25,0))
    msdcumhist = np.cumsum(msdhistogram)
    #
    # fit to it
    #
    if cumhist[minbin]==0:
        fitparams = optimization.curve_fit(linearfunction, edges[(minbin+1):maxbin]-0.5,list(map(math.log10,cumhist[(minbin+1):maxbin])))[0]
        print('Shortened fit range by one bin for ' + accession)
    else:
        fitparams = optimization.curve_fit(linearfunction, edges[minbin:maxbin]-0.5,list(map(math.log10,cumhist[minbin:maxbin])))[0]
    fit_scale=fitparams[0]
    fit_prefactor=10**fitparams[1]
    #
    # plot the histogram
    #
    if plot_histograms:
        fig = plt.figure()
        ax=fig.add_subplot(111)
        plt.rc('xtick',labelsize=tick_label_size)
        plt.tick_params(axis='both',which='major',width=2,length=10,labelsize=tick_label_size)
        plt.tick_params(axis='both',which='minor',width=2,length=3,labelsize=tick_label_size)
        plt.xlabel('TIR free energy [kcal/mol]',**csfont,fontsize=axis_label_size)
        plt.ylabel('number of genes', **csfont,fontsize=axis_label_size)
        plt.yscale('log')
        ax.spines['top'].set_linewidth(2)
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.spines['right'].set_linewidth(2)
        ax.xaxis.set_minor_locator(mticker.AutoMinorLocator(10))
        plt.yticks([1,10,100,1000,10000])
        ax.yaxis.set_major_formatter(mticker.LogFormatter(labelOnlyBase=True))
        ax.yaxis.set_minor_locator(mticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks=50))
        ax.plot(edges[1:]-0.5,cumhist,'-ro')
        ax.plot(edges[1:]-0.5,msdcumhist,'-bo')
        ax.plot(edges[1:]-0.5,fit_prefactor*10**(fit_scale*(edges[1:]-0.5)),'g--',linewidth=2)
        ax.plot(edges[1:]-0.5,0*edges[1:]+1,'g',linewidth=2)
        plt.ylim(bottom=0.2)
        plt.tight_layout()
        plt.savefig('sd_tir_histograms/' + accession + '-sd_min_tirs_fitted.png', bbox_inches='tight')
        plt.close()
    #
    # extract energy where fit line his 1 
    #
    energy_for_one=-math.log10(fit_prefactor)/fit_scale 
    #
    # calculate all TIR energies relative to this one
    #
    offset_energies=[energy for energy in tirenergies-energy_for_one if energy<0]
    #
    # take their histogram
    #
    offset_histogram, offset_edges = np.histogram(offset_energies, bins=offsetbins, range=(-offsetbins,0))
    offset_cumhist = np.cumsum(offset_histogram)
    #
    # plot offset histogram
    #
    #plt.bar(offset_edges[1:]-0.5,offset_cumhist)
    #plt.show()
    #
    # store excess gene counts
    #
    for i in range(len(offset_cumhist)):
        excess_genes[i, min(offset_cumhist[i], maxexcess)]+=1

#
# plot 2D bar plot
#
def update_ticks(x, pos):
    if x == maxexcess:
        return '>' + str(maxexcess-1)
    else:
        return str(int(x))
#
cmap=plt.get_cmap('viridis_r')
_x = [n-0.5 for n in range(maxexcess+1)]
_y = range(offsetbins-1,-1,-1)
_xx, _yy = np.meshgrid(_x,_y)
x,y = _xx.ravel(), _yy.ravel()
_l = np.linspace(0.0,1.0,maxexcess+1)
_lx, _ly = np.meshgrid(_l, _y)
c = cmap(_lx.ravel())
fig = plt.figure()
plt.rc('xtick',labelsize=12)
ax = fig.add_subplot(projection='3d')
ax.bar3d(x, y, np.zeros_like(excess_genes.ravel()), 1, 1, excess_genes.ravel(), shade=True, color=c) 
ax.set_xlabel('number of genes',**csfont,fontsize=14)
ax.set_ylabel('excess TIR ΔG [kcal/mol]',**csfont,fontsize=14)
ax.set_zlabel('number of species',**csfont,fontsize=14)
ax.xaxis.set_ticks(range(maxexcess+1))
ax.xaxis.set_major_formatter(mticker.FuncFormatter(update_ticks))
plt.savefig('ExcessTIRFreeEnergyHistograms.png')
plt.show()

