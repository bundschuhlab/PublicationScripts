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

import csv
import matplotlib.pyplot as plt


def get_list_annotation_genomes(pipeline_output_loc):
    """
    """
    assembly_acc_energy_dict = {}
    assembly_acc_taxonomy_name_dict = {}
    with open(pipeline_output_loc, "r") as f:
        reader = csv.reader(f)

        counter = 0        
        for row in reader:
            counter += 1
            if counter == 1:
                # header
                s_l_prot_names = row[4:]
            else:
                assembly_acc, taxonomy_name, s21_energy, s21_status_call = row[:4]
                s_l_prot_tirs = row[4:]

                assembly_acc_energy_dict[assembly_acc] = s_l_prot_tirs
                assembly_acc_taxonomy_name_dict[assembly_acc] = taxonomy_name
    f.close()

    return assembly_acc_energy_dict, assembly_acc_taxonomy_name_dict, s_l_prot_names


def gen_species_tir_plots(pipeline_output_loc, prot_histogram_energy_loc):
    """
    """
    assembly_acc_energy_dict, assembly_acc_taxonomy_name_dict, s_l_prot_names = get_list_annotation_genomes(pipeline_output_loc)

    for assembly_acc in assembly_acc_energy_dict:
        taxonomy_name = assembly_acc_taxonomy_name_dict[assembly_acc]
        s_l_prot_tirs = assembly_acc_energy_dict[assembly_acc]

        s21_energy = []
        s18_energy = []
        other_energy = []

        for i in range(len(s_l_prot_names)):
            s_l_prot_name = s_l_prot_names[i]
            energy = s_l_prot_tirs[i]

            if energy != "":
                energy = float(energy)

                if s_l_prot_name == "S21":
                    s21_energy.append(energy)
                elif s_l_prot_name == "S18":
                    s18_energy.append(energy)
                else:
                    other_energy.append(energy)
        
        # write out histogram
        plt.figure()
        plt.hist(other_energy, bins = 50, range = (-25,0), label="Other S and L", color = "blue")
        plt.hist(s18_energy, bins = 50, range = (-25,0), label="S18", color = "green")
        plt.hist(s21_energy, bins = 50, range = (-25,0), label="S21", color = "red")

        plt.xlim([-25, 0])
        plt.title(taxonomy_name)
        plt.xlabel("Binding Energy [kcal/mol]")
        plt.ylabel("Frequency")
        plt.legend(loc = "upper left")
        plt.savefig(prot_histogram_energy_loc + taxonomy_name + ".png", bbox_inches='tight')
        plt.close()


    return
