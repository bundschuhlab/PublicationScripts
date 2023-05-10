import os
import csv
from Bio.Seq import Seq


def fasta_to_dict(fasta_in):
    """
    """
    seqid_dict = {}
    with open(fasta_in, "r") as f:
        subj, seq = "", ""

        for row in f:
            if row.startswith(">"):
                if len(subj) != 0:
                    seqid_dict[subj] = seq
                subj = row[1:].replace("\n","").split()[0]
                seq = ""
            else:
                seq = seq + row.replace("\n","")
        
        if len(subj) != 0:
            seqid_dict[subj] = seq
    f.close()

    return seqid_dict


def get_16s_positions(barrnap_annotation_loc):
    """
    """
    barrnap_16s_predictions = []
    with open(barrnap_annotation_loc, "r") as f:
        for row in f:
            if row.startswith(">"):
                seqid = row[1:].replace("\n","")
                if "16S_rRNA" in seqid and "partial" not in seqid:
                    genome_seqid = seqid.split("::")[1].split(":")[0]
                    genome_positions = seqid.split("::")[1].split(":")[1].split("(")[0].split("-")
                    direction = seqid.split("(")[1].split(")")[0]

                    barrnap_16s_predictions.append([genome_seqid, genome_positions, direction])
    f.close()

    return barrnap_16s_predictions


def write_out_16s_asdg_area(barrnap_seqment_loc, barrnap_16s_predictions, \
                                    barrnap_num_downstream, \
                                    barrnap_num_upstream_from_end, \
                                    genome_seq_dict):
    """
    """
    with open(barrnap_seqment_loc, "w") as f:
        for genome_seqid, genome_positions, direction in barrnap_16s_predictions:
            genome_seq = genome_seq_dict[genome_seqid]
            position_start, position_end = genome_positions
            position_start = int(position_start)
            position_end = int(position_end)

            # grab snippet of DNA for outfile writing. Note we set min/maxes to 
            # not grab over genome edges
            if direction == "+":
                snip_start = max(1, position_end - barrnap_num_upstream_from_end)
                snip_end = min(len(genome_seq), position_end + barrnap_num_downstream)
                seq_snip = genome_seq[snip_start:snip_end]
            elif direction == "-":
                snip_start = max(1, position_start - barrnap_num_downstream)
                snip_end = min(position_start + barrnap_num_upstream_from_end, len(genome_seq))
                seq_snip = genome_seq[snip_start:snip_end]
                dna = Seq(seq_snip)
                seq_snip = str(dna.reverse_complement())
            else:
                print("issue with direction: ", genome_positions, direction)
                break

            f.write(">" + str(genome_positions) + "|" + direction + "\n")
            f.write(seq_snip + "\n")
    f.close()

    return 


def downselect_cutadapt_trimmed(cutadapt_out_loc, asd_prediction_loc):
    """
    """
    cutadapt_seqid_dict = fasta_to_dict(cutadapt_out_loc)
    
    with open(asd_prediction_loc, "w") as f:
        for seqid in cutadapt_seqid_dict:
            cutadapt_seq = cutadapt_seqid_dict[seqid]

            snip_of_interest = cutadapt_seq[17:17+15]

            f.write(">" + seqid + "-asdg-predict" + "\n")
            f.write(snip_of_interest + "\n")
    f.close()

    return


def identify_asd_seq(asd_prediction_loc):
    """
    """
    asd_seqid_dict = fasta_to_dict(asd_prediction_loc)

    uniq_seqs = set()
    for seqid in asd_seqid_dict:
        seq = asd_seqid_dict[seqid]
        uniq_seqs.add(seq)

    return uniq_seqs


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


def extract_sequence_snip(locus_tag, genome_fasta_dict, positions, grab_extra_nt_at_end = False):
    """
    """
    genome_seq = genome_fasta_dict[locus_tag]

    position_start, position_end = positions.split("..")

    if grab_extra_nt_at_end == True:
        position_end = int(position_end) + 1

    genome_snip = genome_seq[int(position_start) - 1 : int(position_end)]

    return genome_snip


def get_tax_rank_of_interest(taxonomy, prefix_of_interest, tree):
    """
    """
    parent_nodes = tree.ascend(taxonomy.lower())
    for node in parent_nodes:
        node_t = node.taxid
        if prefix_of_interest in node_t:
            return node_t

    return ""


def get_delg_values(arfa_summary_loc, assembly_data_summary, \
                        assembly_dl_root, arfa_output_dir, \
                        tree):
    """
    """
    barrnap_exe = "barrnap"
    cutadapt_exe = "$HOME/.local/bin/cutadapt"
    freescan_exe = "perl ../tools_resources/free2bind/free_scan.pl"

    arfa_call_genomeacc = {}
    arfa_call_genomeacc["True"] = []
    arfa_call_genomeacc["False"] = []

    acc_taxonomy_lookup = {}
    with open(assembly_data_summary, "r") as f:
        for row in f:
            row = row.replace("\n","").split("\t")
            taxonomy, acc = row

            acc_taxonomy_lookup[acc] = taxonomy
    f.close()

    # get list of accessions with true/false prfB call from ARFA
    with open(arfa_summary_loc, "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            genome_acc, contains_prfb_annotation, prfb_has_orf0, prfb_orf0_evalue, prfb_status_call = row

            if prfb_status_call == "True":
                arfa_call_genomeacc["True"].append(genome_acc)
            elif prfb_status_call == "False":
                arfa_call_genomeacc["False"].append(genome_acc)
    f.close()

    with open("results/bacteria-prfb_asd_freescan_updated.csv", "w") as f:
        out = csv.writer(f)
        out.writerow(["Accession", "Taxonomy", "Family Tax", "Order Tax", "Class Tax", "Phylum Tax", \
                        "prfB FS", "ASD", "# of Uniq ASDs", \
                        "SD Region", "MSD Region", \
                        "SD Freescan MinE", "MSD Freescan MinE"])

        for status in ["True"]:
            nucl_prfb_seqdict = fasta_to_dict("results/prfB_" + status.upper() + "_nuc.fa")
            genome_accs = arfa_call_genomeacc[status]

            for genome_acc in genome_accs:
                taxonomy = acc_taxonomy_lookup[genome_acc]

                # get higher rank info
                family_tax = get_tax_rank_of_interest(taxonomy, "f__", tree)
                order_tax = get_tax_rank_of_interest(taxonomy, "o__", tree)
                class_tax = get_tax_rank_of_interest(taxonomy, "c__", tree)
                phylum_tax = get_tax_rank_of_interest(taxonomy, "p__", tree)


                # run barrnap on the genome
                fasta_loc = assembly_dl_root + genome_acc + "_genomic.fna"

                barrnap_out = "resources/barrnap.fasta"
                command = barrnap_exe + \
                        " " + fasta_loc + \
                        " --outseq " + barrnap_out + \
                        " --quiet"
                os.system(command)

                barrnap_16s_predictions = get_16s_positions(barrnap_out)

                genome_seq_dict = fasta_to_dict(fasta_loc)

                barrnap_seqment_loc = "resources/barrnap-16S-segs.fa"
                write_out_16s_asdg_area(barrnap_seqment_loc, barrnap_16s_predictions, \
                                        10, \
                                        100, \
                                        genome_seq_dict)

                # run cutadapt
                cutadapt_out_loc = "resources/cutadapt.fa"
                command = cutadapt_exe + \
                                " -g AAGTCGTAACAAGGTAGCCGT" + \
                                " -e 0.25 -O 21" + \
                                " --discard-untrimmed" + \
                                " -o " + cutadapt_out_loc + \
                                " " + barrnap_seqment_loc + \
                                " -j 0 --quiet"
                os.system(command)

                asd_prediction_loc = "resources/cutadapt-snips.fa"
                downselect_cutadapt_trimmed(cutadapt_out_loc, asd_prediction_loc)

                uniq_asd_seqs = identify_asd_seq(asd_prediction_loc)
                num_asd_seqs = len(uniq_asd_seqs)
                #if num_asd_seqs == 1 and genome_acc in nucl_prfb_seqdict:
                if num_asd_seqs == 1:
                    for s in uniq_asd_seqs:
                        asd_seq = s
                else:
                    asd_seq = ""

                    out.writerow([genome_acc, taxonomy, family_tax, order_tax, class_tax, phylum_tax, \
                                    status, asd_seq, num_asd_seqs, "", "", "", ""])
                    continue

                # if asd_seq is longer than 0, 
                #genome_prfb_sequence = nucl_prfb_seqdict[genome_acc]
                # extract out region of prfB sequence of interest

                arfa_output_loc = arfa_output_dir + genome_acc + "_arfa.txt"

                with open(arfa_output_loc, "r") as r:
                    prfb_annotation = False
                    translated_region = False
                    prfb_sequence = ""
                    for row in r:
                        if row.startswith('                     /gene="prfB"'):
                            prfb_annotation = True
                        elif prfb_annotation == True and row.startswith("     CDS             "):
                            positions = row.split("CDS")[1].replace("\n","").strip()
                        elif prfb_annotation == True and row.startswith("                     /locus_tag="):
                            locus_tag = row.split('/locus_tag="')[1].split("|")[0]
                        elif prfb_annotation == True and row.startswith('                     /translation='):
                            prfb_sequence = prfb_sequence + row.split("/translation=")[1]
                            translated_region = True
                        elif row.startswith("     gene"):
                            prfb_annotation = False
                            translated_region = False
                        elif prfb_annotation == True and translated_region == True:
                            prfb_sequence = prfb_sequence + row.strip()
                r.close()

                # get prfb nucleotide sequence
                first_set_pos = positions.split("join(")[1].split(",")[0]
                second_set_pos = positions.split(",")[1].split(")")[0]

                if "complement" in positions:
                    fs_start_pos = int(second_set_pos.split("..")[1]) - int(second_set_pos.split("..")[0])
                else:
                    fs_start_pos = int(first_set_pos.split("..")[1]) - int(first_set_pos.split("..")[0])

                # get prfb nucleotide sequence
                first_set_pos = positions.split("join(")[1].split(",")[0]
                second_set_pos = positions.split(",")[1].split(")")[0]

                if "complement" in positions:
                    num_bp_in_fs = int(second_set_pos.split("..")[1]) - int(second_set_pos.split("..")[0])
                else:
                    num_bp_in_fs = int(first_set_pos.split("..")[1]) - int(first_set_pos.split("..")[0])

                genome_fasta_dict = fasta_to_dict("/home/fung/Documents/codon_freqs/resources/ncbi-data/" + genome_acc + "_genomic.fna")
                nuc_seq_portion = extract_sequence_snip(locus_tag, genome_fasta_dict, first_set_pos, True) + \
                                extract_sequence_snip(locus_tag, genome_fasta_dict, second_set_pos, False)

                if "complement" in positions:
                    dna = Seq(nuc_seq_portion)
                    nuc_seq_portion = str(dna.reverse_complement())

                # get sd/msd
                #sd_seq = genome_prfb_sequence[fs_start_pos - 25 : fs_start_pos + 1]
                #msd_seq = genome_prfb_sequence[fs_start_pos - 50 : fs_start_pos - 25 + 1]
                sd_seq = nuc_seq_portion[fs_start_pos - 20 : fs_start_pos + 4]
                msd_seq = nuc_seq_portion[fs_start_pos + 50 : fs_start_pos + 50 + 20 + 4]

                with open("resources/sd_seq.fa", "w") as f:
                    f.write(">SD" + "\n" + sd_seq)
                f.close()
                with open("resources/msd_seq.fa", "w") as f:
                    f.write(">MSD" + "\n" + msd_seq)
                f.close()

                sd_out_loc = "resources/sd_out.txt"
                msd_out_loc = "resources/msd_out.txt"
                # run freescan for each
                command = freescan_exe + " -q " + asd_seq[::-1] + \
                                " " + "resources/sd_seq.fa" + \
                                " > " + sd_out_loc
                os.system(command)
                command = freescan_exe + " -q " + asd_seq[::-1] + \
                                " " + "resources/msd_seq.fa" + \
                                " > " + msd_out_loc
                os.system(command)

                sd_min_e = get_min_delG_from_freescan_output(sd_out_loc)
                msd_min_e = get_min_delG_from_freescan_output(msd_out_loc)

                out.writerow([genome_acc, taxonomy, family_tax, order_tax, class_tax, phylum_tax, \
                                status, asd_seq, num_asd_seqs, sd_seq, msd_seq, sd_min_e, msd_min_e])
    f.close()

    return

def gen_itol_tree_2(sp_clusters, tree):
    """
    """
    from ete3 import Tree
    acc_taxid_dict = {}
    passing_taxonomies, failing_taxonomies = [], []
    middle_taxonomies = []
    extra_passing_taxonomies = []
    acc_signal_dict = {}
    with open("results/bacteria-prfb_asd_freescan_updated.csv", "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            acc = row[0] # row[1].lower()
            taxonomy = row[1]

            acc_taxid_dict[acc] = taxonomy
            num_uniq_asd = int(row[8])

            if num_uniq_asd != 1:
                continue

            phylum_tax = row[5]
            if phylum_tax != "p__bacteroidota":
                continue

            sd_e = float(row[11])
            acc_signal_dict[acc] = sd_e

            if sd_e <= -10.0:
                extra_passing_taxonomies.append(acc)
            elif sd_e <= -7.0:
                passing_taxonomies.append(acc)
            elif sd_e <= -4.0:
                middle_taxonomies.append(acc)
            else:
                failing_taxonomies.append(acc)
    f.close()

    gtdb_full_tree_loc = "resources/bac120.tree"
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
    for acc in passing_taxonomies + failing_taxonomies + middle_taxonomies + extra_passing_taxonomies:
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
    downselected_node.write(format = 1, outfile = "results/itol/itol.tree")

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
    prfb_compar_desired_taxa_to_shade = []
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


    with open("results/itol/" + "tree-labels.txt", "w") as f:
        f.write("TREE_COLORS" + "\n")
        f.write("SEPARATOR TAB" + "\n")
        f.write("DATA" + "\n")

        # clade shading
        for taxa in taxa_node_map:
            node_name = taxa_node_map[taxa]
            f.write(node_name.replace(":","_").replace(";","_") + "\trange\t#ffffff\t" + taxa + "\n")

        for acc in failing_taxonomies:
            node_name = acc_node_map[acc]
            f.write(node_name + "\tlabel\t#0000FF" + "\n") # blue
            f.write(node_name + "\tbranch\t" + "#0000FF" + "\tnormal\t4\n")
        for acc in middle_taxonomies:
            node_name = acc_node_map[acc]
            f.write(node_name + "\tlabel\t#FFC0CB" + "\n") # pink
            f.write(node_name + "\tbranch\t" + "#FFC0CB" + "\tnormal\t4\n")
        for acc in passing_taxonomies:
            node_name = acc_node_map[acc]
            f.write(node_name + "\tlabel\t#FF0000" + "\n") # red
            f.write(node_name + "\tbranch\t" + "#FF0000" + "\tnormal\t4\n")
        for acc in extra_passing_taxonomies:
            node_name = acc_node_map[acc]
            f.write(node_name + "\tlabel\t#800020" + "\n") # burgundy
            f.write(node_name + "\tbranch\t" + "#800020" + "\tnormal\t4\n")
    f.close()

    with open("results/itol/" + "tree-names.txt", "w") as f:
        f.write("LABELS" + "\n")
        f.write("SEPARATOR TAB" + "\n")
        f.write("DATA" + "\n")

        for acc in passing_taxonomies + failing_taxonomies + middle_taxonomies + extra_passing_taxonomies:
            taxonomy = acc_taxid_dict[acc]
            node_name = acc_node_map[acc]
            acc_signal = acc_signal_dict[acc]

            f.write(node_name + "\t" + taxonomy + " [" + str(round(acc_signal, 2)) + "]" + "\n")
    f.close()

    return


def gen_itol_tree_input(tree):
    """
    """
    passing_taxonomies, failing_taxonomies = [], []
    with open("results/prfb_asd_freescan.csv", "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            taxonomy = row[1].lower()
            num_uniq_asd = int(row[4])

            if num_uniq_asd != 1:
                continue

            sd_e = float(row[7])

            if sd_e <= -5.0:
                passing_taxonomies.append(taxonomy)
            else:
                failing_taxonomies.append(taxonomy)
    f.close()
    
    txt_lineages = []
    with open("results/" + "tree-labels.txt", "w") as f:
        f.write("TREE_COLORS" + "\n")
        f.write("SEPARATOR TAB" + "\n")
        f.write("DATA" + "\n")

        for taxonomy in passing_taxonomies:
            parents = tree.ascend(taxonomy)
            lineage = []
            for node in parents:
                node_taxonomy = node.taxid

                if node_taxonomy not in lineage:
                    lineage.append(node_taxonomy)
            lineage.reverse()

            itol_str = ""
            for l_taxonomy in lineage:
                itol_str = itol_str + l_taxonomy + ";"
            itol_str = itol_str

            f.write(taxonomy + "\tlabel\t#F91607" + "\n")
            #f.write(taxonomy + "\tbranch\t#F91607\tnormal\t4\n")
            txt_lineages.append(itol_str)

        for taxonomy in failing_taxonomies:
            parents = tree.ascend(taxonomy)
            lineage = []
            for node in parents:
                node_taxonomy = node.taxid

                if node_taxonomy not in lineage:
                    lineage.append(node_taxonomy)
            lineage.reverse()

            itol_str = ""
            for l_taxonomy in lineage:
                itol_str = itol_str + l_taxonomy + ";"
            itol_str = itol_str

            f.write(taxonomy + "\tlabel\t#0719F9" + "\n")
            #f.write(taxonomy + "\tbranch\t#0719F9\tnormal\t4\n")
            txt_lineages.append(itol_str)
    f.close()

    command = 'echo "'
    for line in txt_lineages:
        command = command + line + "\n"
    command = command[:-1]

    command = command + '" | java Biostar52895'
    command = command + " > " + "results/" + "itol_tree.txt"
    os.system(command)

    return


def id_interesting_clades(tree):
    #list_level_of_interest = ["family", "order", "class", "phylum"]
    list_level_of_interest = ["phylum"]

    acc_tax_dict = {}
    with open("results/bacteria-ncbi-download-summary.txt", "r") as f:
        for row in f:
            row = row.replace("\n","").split("\t")
            tax, acc = row
            acc_tax_dict[acc] = tax.lower()
    f.close()

    for level_of_interest in list_level_of_interest:
        tax_dict_fracts = {}
        with open("results/bacteria-prfb_asd_freescan_updated.csv", "r") as f:
            reader = csv.reader(f)
            next(reader, None)

            for row in reader:
                num_uniq_asd = int(row[8])

                if num_uniq_asd != 1:
                    continue

                if level_of_interest == "family":
                    tax = row[2]
                elif level_of_interest == "order":
                    tax = row[3]
                elif level_of_interest == "class":
                    tax = row[4]
                elif level_of_interest == "phylum":
                    tax = row[5]

                if tax not in tax_dict_fracts:
                    tax_dict_fracts[tax] = [0, 0]
                tax_dict_fracts[tax][0] += 1

                sd_e = float(row[11])
                if sd_e > -5.0:
                    tax_dict_fracts[tax][1] += 1
        f.close()

        tax_tot_num_counts = {}
        with open("results/bacteria-arfa-summary.csv", "r") as f:
            reader = csv.reader(f)
            next(reader, None)

            for row in reader:
                acc = row[0]
                tax = acc_tax_dict[acc]

                p_lineage = tree.ascend(tax)
                phylum_tax = "NA"
                for node in p_lineage:
                    node_tax = node.taxid

                    if "p__" in node_tax:
                        phylum_tax = node_tax

                if phylum_tax not in tax_tot_num_counts:
                    tax_tot_num_counts[phylum_tax] = 0
                tax_tot_num_counts[phylum_tax] += 1
        f.close()
            
        out_loc = "results/tax-levels/" + level_of_interest + "-splices.csv"
        with open(out_loc, "w") as f:
            out = csv.writer(f)
            out.writerow(["Taxonomy", "Tot # of Species with prfB FS", \
                            "Num Species w/ dG > -5.0", \
                            "Fraction"])
            
            for tax in tax_dict_fracts:
                tot_num, num_blue = tax_dict_fracts[tax]
                tot_num_screened = tax_tot_num_counts[tax]

                out.writerow([tax, tot_num_screened, tot_num, num_blue, num_blue / tot_num])
        f.close()


    return


def scan_sd_seqs(arfa_summary_loc, assembly_dl_root, arfa_output_dir, assembly_data_summary, tree, sp_clusters):
    """
    """
    #get_delg_values(arfa_summary_loc, assembly_data_summary, \
    #                    assembly_dl_root, arfa_output_dir, \
    #                    tree)

    #gen_itol_tree_input(tree)
    gen_itol_tree_2(sp_clusters, tree)

    #id_interesting_clades(tree)

    return
