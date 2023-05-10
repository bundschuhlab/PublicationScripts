import csv
from ete3 import Tree


def generate_tree_of_results(species_genome_dict, genome_species_dict, \
                                genome_wide_results, \
                                tree, sp_clusters, gtdb_full_tree_loc):
    """
    """
    # get downselected tree
    accs_to_grab = set()
    with open("../codon_freqs/results/bacteria-prfb_asd_freescan_updated.csv", "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            acc = row[0]
            num_uniq_asd = int(row[8])

            if num_uniq_asd != 1:
                continue

            accs_to_grab.add(acc)
    f.close()
    #accs_to_grab = set()


    #tree_file_fldr = "results/tree_diffs_full/"
    tree_file_fldr = "results/tree_fractions_prfbfs-compare/"

    import os
    if os.path.isdir(tree_file_fldr) == False:
        os.system("mkdir " + tree_file_fldr)
    
    acc_signals = {}
    with open(genome_wide_results, "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            acc, species, genus, family, \
                            order, classs, phylum, domain, \
                            average_diff, fract_below_thresh, \
                            fract_sd, fract_msd = row

            average_diff = float(average_diff)
            fract_below_thresh = float(fract_below_thresh)

            if len(accs_to_grab) > 0 and acc not in accs_to_grab:
                continue

            if acc not in acc_signals:
                acc_signals[acc] = []
            
            #acc_signals[acc].append(average_diff)
            acc_signals[acc].append(fract_below_thresh)
    f.close()

    acc_signal_dict = {}
    all_signals = []
    for acc in acc_signals:
        signals = acc_signals[acc]

        avg_signal = sum(signals) / len(signals)

        acc_signal_dict[acc] = avg_signal
        all_signals.append(avg_signal)

    t = Tree(gtdb_full_tree_loc, format=1, quoted_node_names = True)

    # extract the region of interest (p__Bacteroidota)
    clade_tree_search = t.search_nodes(name = "d__Bacteria")[0] # d__Bacteria 99.0:c__Gammaproteobacteria 45.0:p__Firmicutes 100.0:p__Bacteroidota
    downselected_node = clade_tree_search.detach()
    
    acc_fullnode_dict = {}
    with open(sp_clusters, "r") as f:
        next(f)
        for row in f:
            data = row.split("\t")[0]
            
            if "GB_" in data:
                acc = data.split("GB_")[1]
            elif "RS_" in data:
                acc = data.split("RS_")[1]
            
            acc_fullnode_dict[acc] = data
    f.close()

    # now prune the tree to only those nodes of interest
    desired_nodes = set()
    acc_node_map = {}
    for acc in acc_signal_dict:
        node_name = acc_fullnode_dict[acc]

        acc_node_map[acc] = node_name
        desired_nodes.add(node_name)

    # note- if we want to preserve internal nodes (we do I think) - need to get a list of them
    # BUT - we only want internal nodes that are a parent of a leaf we're keeping
    internal_labels = set()
    clean_external_nodes = set()
    for node in downselected_node.search_nodes():
        if node.is_leaf() == False:
            if "__" in node.name:
                # check if leaves in set
                node_leaves = node.get_leaf_names()
                for n in node_leaves:
                    if n in desired_nodes:
                        internal_labels.add(node.name)
        else:
            if node.name in desired_nodes:
                clean_external_nodes.add(node.name)

    downselected_node.prune(list(clean_external_nodes.union(internal_labels)), preserve_branch_length = True)
    downselected_node.write(format = 1, outfile = tree_file_fldr + "itol.tree")

    # create a color gradient from min to max observed z-score
    min_signal, max_signal = min(all_signals), max(all_signals)
    from colour import Color
    light_red = Color("#ffcccb") # white: #FFFFFF "#ffcccb"
    dark_red = Color("#F91607")
    colors = list(light_red.range_to(dark_red,101))

    white = Color("#5f96fe") # "#ADD8E6" light blue b/c white --> blue includes green which is weird
    blue = Color("#0000FF")
    above_0_colors = list(white.range_to(blue, 101))

    # get lookup of nodes for desired color shading
    prfb_compar_desired_taxa_to_shade = [
        "p__Proteobacteria", \
        "p__Firmicutes", \
        "p__Bacteroidota", \
        "p__Cyanobacteria", \
        "p__Firmicutes_A", \
        "p__Desulfobacterota", \
        "p__Planctomycetota", \
        "p__Spirochaetota", \
        "p__Firmicutes_B", \
        "p__Myxococcota", \
        "p__Fusobacteriota", \
        "p__Desulfobacterota_I", \
        "p__Deinococcota", \
        "p__Actinobacteriota", \
        #"p__Campylobacterota"
    ]
    #prfb_compar_desired_taxa_to_shade = []
    taxa_node_map = {}
    for taxa in prfb_compar_desired_taxa_to_shade:
        found_node = False

        for node in downselected_node.search_nodes():
            node_name = node.name
            if taxa in node_name:
                # control for underscore after taxa
                if len(node_name.split(taxa)[1]) == 0 or node_name.split(taxa)[1][0] != "_":
                    if taxa in taxa_node_map:
                        print("WARNING - duplicate taxa nodes??", "t:", taxa, "n1:", node_name, "n2:", taxa_node_map[taxa])
                    taxa_node_map[taxa] = node_name
                    found_node = True
        
        if found_node == False:
            print("WARNING - unable to find node for taxa:", taxa)

    with open(tree_file_fldr + "tree-labels.txt", "w") as f:
        f.write("TREE_COLORS" + "\n")
        f.write("SEPARATOR TAB" + "\n")
        f.write("DATA" + "\n")

        # clade shading
        for taxa in taxa_node_map:
            node_name = taxa_node_map[taxa]
            f.write(node_name.replace(":","_").replace(";","_") + "\trange\t#ffffff\t" + taxa + "\n")

        for acc in acc_signal_dict:
            node_name = acc_node_map[acc]
            acc_signal = acc_signal_dict[acc]

            #color = colors[int(100 * (acc_signal - min_signal) / (max_signal - min_signal))]
            if acc_signal > 0.0:
                #color = above_0_colors[int(100 * (acc_signal - 0.0) / (max_signal - 0.0))]
                if "tree_diffs" in tree_file_fldr:
                    color = above_0_colors[int(100 * (acc_signal - 0.0) / (max_signal - 0.0))]
                else:
                    color = colors[int(100 * (acc_signal - 0.0) / (max_signal - 0.0))]
            else:
                #color = colors[100 - int(100 * (acc_signal - min_signal) / (0.0 - min_signal))]
                if "tree_diffs" in tree_file_fldr:
                    color = colors[100 - int(100 * (acc_signal - min_signal) / (0.0 - min_signal))]
                else:
                    color = above_0_colors[100 - int(100 * (acc_signal - min_signal) / (0.0 - min_signal))]

            f.write(node_name + "\tlabel\t" + str(color) + "\n")
            f.write(node_name + "\tbranch\t" + str(color) + "\tnormal\t4\n")
    f.close()

    with open(tree_file_fldr + "tree-names.txt", "w") as f:
        f.write("LABELS" + "\n")
        f.write("SEPARATOR TAB" + "\n")
        f.write("DATA" + "\n")

        for acc in acc_signal_dict:
            taxonomy = genome_species_dict[acc]
            node_name = acc_node_map[acc]
            acc_signal = acc_signal_dict[acc]

            f.write(node_name + "\t" + taxonomy + " [" + str(round(acc_signal, 2)) + "]" + "\n")
    f.close()

    return