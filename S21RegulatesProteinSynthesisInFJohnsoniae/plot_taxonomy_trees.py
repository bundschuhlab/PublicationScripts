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

import os
import csv


def mine_out_species(pipeline_output_loc, min_delg_for_s21_winner):
    """
     each dict {taxonomy_name : S21 energy}
    """
    s21_winner_dict = {}
    s21_loser_dict = {}
    s21_missing_dict = {}

    with open(pipeline_output_loc, "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        counter = 0        
        for row in reader:
            assembly_acc, taxonomy_name, s21_energy, s21_status_call = row[:4]
            taxonomy_name = taxonomy_name.lower()

            if s21_energy == "NA":
                s21_missing_dict[taxonomy_name] = ""
            else:
                if s21_status_call == "S21_Winner" and float(s21_energy) <= min_delg_for_s21_winner:
                    s21_winner_dict[taxonomy_name] = str(round(float(s21_energy), 2))
                else:
                    s21_loser_dict[taxonomy_name] = str(round(float(s21_energy), 2))
    f.close()

    return s21_winner_dict, s21_loser_dict, s21_missing_dict


def check_for_bact(taxonomy_list_in, tree):
    """
    """
    taxonomy_out = []
    for taxonomy in taxonomy_list_in:
        taxonomy = taxonomy.lower()
        parent_nodes = tree.ascend(taxonomy)

        is_bact = False
        for node in parent_nodes:
            node_taxonomy = node.taxid
            if node_taxonomy == "p__bacteroidota":
                is_bact = True
        
        if is_bact == True:
            taxonomy_out.append(taxonomy)

    return taxonomy_out


def gen_tree_plots(tree_plot_loc, min_delg_for_s21_winner, tree, \
                        pipeline_output_loc):
    """
    """
    # get dicts of s21_winner and s21_loser
    # use deltaG condition
    s21_winner_dict, s21_loser_dict, s21_missing_dict = \
        mine_out_species(pipeline_output_loc, min_delg_for_s21_winner)

    # downselect to just bacteroidetes
    bact_s21_winners = check_for_bact(s21_winner_dict, tree)
    bact_s21_losers = check_for_bact(s21_loser_dict, tree)

    from colour import Color
    red = Color("#facecb") 
    colors = list(red.range_to(Color("#F91607"),11))
    print(colors)

    txt_lineages = []
    with open(tree_plot_loc + "tree-labels.txt", "w") as f:
        f.write("TREE_COLORS" + "\n")
        f.write("SEPARATOR TAB" + "\n")
        f.write("DATA" + "\n")

        # write out each order taxonomy color
        f.write("o__chlorobiales" + "\t" + "range" + "\t" + "#eeeeee" + "\t" + "Chlorobiales" + "\n")
        f.write("o__chitinophagales" + "\t" + "range" + "\t" + "#ddffff" + "\t" + "Chitinophagales" + "\n")
        f.write("o__sphingobacteriales" + "\t" + "range" + "\t" + "#ffeeee" + "\t" + "Sphingobacteriales" + "\n")
        f.write("o__bacteroidales" + "\t" + "range" + "\t" + "#ffddff" + "\t" + "Bacteroidales" + "\n")
        f.write("o__cytophagales" + "\t" + "range" + "\t" + "#ffffdd" + "\t" + "Cytophagales" + "\n")
        f.write("o__flavobacteriales" + "\t" + "range" + "\t" + "#eeeeff" + "\t" + "Flavobacteriales" + "\n")

        for taxonomy in bact_s21_winners:
            parents = tree.ascend(taxonomy)
            lineage = []
            for node in parents:
                node_taxonomy = node.taxid

                if node_taxonomy not in lineage:
                    lineage.append(node_taxonomy)
            lineage.reverse()

            itol_str = ""
            for l_taxonomy in lineage:
                """
                if l_taxonomy in s21_winner_dict:
                    energy = s21_winner_dict[l_taxonomy]
                    l_taxonomy = l_taxonomy + " [" + str(energy) + "]"
                """
                itol_str = itol_str + l_taxonomy + ";"
            itol_str = itol_str

            #n_taxonomy = taxonomy + " [" + str(s21_winner_dict[taxonomy]) + "]"
            f.write(taxonomy + "\tlabel\t#F91607" + "\n")
            f.write(taxonomy + "\tbranch\t#F91607\tnormal\t4\n")
            txt_lineages.append(itol_str)

        for taxonomy in bact_s21_losers:
            parents = tree.ascend(taxonomy)
            lineage = []
            for node in parents:
                node_taxonomy = node.taxid

                if node_taxonomy not in lineage:
                    lineage.append(node_taxonomy)
            lineage.reverse()

            itol_str = ""
            for l_taxonomy in lineage:
                """
                if l_taxonomy in s21_loser_dict:
                    energy = s21_loser_dict[l_taxonomy]
                    l_taxonomy = l_taxonomy + " [" + str(energy) + "]"
                """
                itol_str = itol_str + l_taxonomy + ";"
            itol_str = itol_str

            #n_taxonomy = taxonomy + " [" + str(s21_loser_dict[taxonomy]) + "]"
            f.write(taxonomy + "\tlabel\t#0719F9" + "\n")
            f.write(taxonomy + "\tbranch\t#0719F9\tnormal\t4\n")
            txt_lineages.append(itol_str)
    f.close()

    with open(tree_plot_loc + "tree-leaf-labels.txt", "w") as f:
        f.write("LABELS" + "\n")
        f.write("SEPARATOR TAB" + "\n")
        f.write("DATA" + "\n")

        for taxonomy in bact_s21_winners:
            n_taxonomy = taxonomy + " [" + str(s21_winner_dict[taxonomy]) + "]"
            f.write(taxonomy + "\t" + n_taxonomy + "\n")

        for taxonomy in bact_s21_losers:
            n_taxonomy = taxonomy + " [" + str(s21_loser_dict[taxonomy]) + "]"
            f.write(taxonomy + "\t" + n_taxonomy + "\n")
    f.close()

    command = 'echo "'
    for line in txt_lineages:
        command = command + line + "\n"
    command = command[:-1]

    command = command + '" | java Biostar52895'
    command = command + " > " + tree_plot_loc + "itol_tree.txt"
    os.system(command)

    return