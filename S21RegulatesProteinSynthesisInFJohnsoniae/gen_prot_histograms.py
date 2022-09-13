# Copyright (C) <2022>  <The Ohio State University>       

# This program is free software: you can redistribute it and/or modify                              
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or    
# (at your option) any later version.                                                                                       
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of           
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
# GNU General Public License for more details.                                                                             

# You should have received a copy of the GNU General Public License 
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import csv
import os

mpl.rcParams['axes.linewidth'] = 2 #set the value globally
font_manager.findSystemFonts(fontpaths=None, fontext="ttf")
font_manager.findfont("Arial") # Test with "Special Elite" too
csfont = {'fontname':'arial'}


def get_list_assembly_accs(pipeline_output_loc, tree):
    """
    """
    assembly_accs_to_grab = []
    with open(pipeline_output_loc, "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            taxonomy = row[1].lower()
            assembly_acc = row[0]

            parent_nodes = tree.ascend(taxonomy)
            is_bact = False
            for node in parent_nodes:
                node_taxonomy = node.taxid
                if node_taxonomy == "c__bacteroidia":
                    is_bact = True
            
            if is_bact == True:
                assembly_accs_to_grab.append(assembly_acc)
    f.close()

    return assembly_accs_to_grab


def get_min_delG_from_freescan_output(freescan_output_file):
    """
    """
    min_e = 0.0
    with open(freescan_output_file, "r") as f:
        for row in f:
            if row[0] != "#" and row.strip() != "":
                signal = float(row.split()[0])

                if signal < min_e:
                    min_e = signal
    f.close()

    return min_e


def mine_out_freescan_output(prot_dict, prot_loc):
    """
    """
    for filename in os.listdir(prot_loc):
        prot_name = filename.split("-")[0]
        min_E = get_min_delG_from_freescan_output(prot_loc + filename)

        if prot_name not in prot_dict:
            prot_dict[prot_name] = []
        prot_dict[prot_name].append(min_E)

    return prot_dict


def plot_species_s_l_tirs(png_out, png_title, tir_list, plt_color):
    """
    """
    # write out histogram
    plt.figure()
    plt.hist(tir_list, bins = np.linspace(-21, 0, num = 0 - -21 + 1), color = plt_color, density = True)

    plt.ylim([0, 1])
    plt.xlim([-25, 0])
    plt.title(png_title, **csfont, fontsize = 34)
    plt.xlabel("ΔG [kcal/mol]", **csfont, fontsize = 24)
    plt.ylabel("Frequency", **csfont, fontsize = 24)
    plt.tick_params(axis='both', which='major', width=2, length=10, labelsize=24)
    plt.savefig(png_out + png_title.replace(" ","-") + ".png", bbox_inches='tight')
    plt.close()

    return


def get_s_l_prot_histograms(pipeline_output_loc, bio_data_loc, histogram_plot_loc, \
                                tree):
    """
    """
    assembly_accs_to_grab = get_list_assembly_accs(pipeline_output_loc, tree)

    s_prot_sd_dict = {}
    s_prot_msd_dict = {}
    l_prot_sd_dict = {}
    l_prot_msd_dict = {}

    for assembly_acc in assembly_accs_to_grab:
        # grab SD, MSD for all available S and L proteins
        assembly_bio_data_loc = bio_data_loc + assembly_acc + "/"

        s_prot_loc = assembly_bio_data_loc + "s_proteins/"
        s_prot_sd_loc = s_prot_loc + "tir_bind/"
        s_prot_msd_loc = s_prot_loc + "tir_bind_MSD/"

        l_prot_loc = assembly_bio_data_loc + "l_proteins/"
        l_prot_sd_loc = l_prot_loc + "tir_bind/"
        l_prot_msd_loc = l_prot_loc + "tir_bind_MSD/"

        # get TIRs
        s_prot_sd_dict = mine_out_freescan_output(s_prot_sd_dict, s_prot_sd_loc)
        s_prot_msd_dict = mine_out_freescan_output(s_prot_msd_dict, s_prot_msd_loc)
        l_prot_sd_dict = mine_out_freescan_output(l_prot_sd_dict, l_prot_sd_loc)
        print(l_prot_sd_loc)
        l_prot_msd_dict = mine_out_freescan_output(l_prot_msd_dict, l_prot_msd_loc)
    
    # write out each TIR's for histograms
    histogram_plot_loc_specific = histogram_plot_loc + "s_prot_SD/"
    os.system("rm -rf " + histogram_plot_loc_specific)
    os.system("mkdir " + histogram_plot_loc_specific)
    for prot_name in s_prot_sd_dict:
        tir_list = s_prot_sd_dict[prot_name]
        plot_species_s_l_tirs(histogram_plot_loc_specific, prot_name + " " + "SD", tir_list, "blue")

    histogram_plot_loc_specific = histogram_plot_loc + "s_prot_MSD/"
    os.system("rm -rf " + histogram_plot_loc_specific)
    os.system("mkdir " + histogram_plot_loc_specific)
    for prot_name in s_prot_msd_dict:
        tir_list = s_prot_msd_dict[prot_name]
        plot_species_s_l_tirs(histogram_plot_loc_specific, prot_name + " " + "MSD", tir_list, "grey")

    histogram_plot_loc_specific = histogram_plot_loc + "l_prot_SD/"
    os.system("rm -rf " + histogram_plot_loc_specific)
    os.system("mkdir " + histogram_plot_loc_specific)
    for prot_name in l_prot_sd_dict:
        tir_list = l_prot_sd_dict[prot_name]
        plot_species_s_l_tirs(histogram_plot_loc_specific, prot_name + " " + "SD", tir_list, "blue")

    histogram_plot_loc_specific = histogram_plot_loc + "l_prot_MSD/"
    os.system("rm -rf " + histogram_plot_loc_specific)
    os.system("mkdir " + histogram_plot_loc_specific)
    for prot_name in l_prot_msd_dict:
        tir_list = l_prot_msd_dict[prot_name]
        plot_species_s_l_tirs(histogram_plot_loc_specific, prot_name + " " + "MSD", tir_list, "grey")

    return