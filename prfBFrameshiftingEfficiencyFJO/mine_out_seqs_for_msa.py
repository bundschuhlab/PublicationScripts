import os
import csv
from Bio.Seq import Seq
import pandas as pd
import logomaker
import matplotlib.pyplot as plt

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


def extract_sequence_snip(locus_tag, genome_fasta_dict, positions, grab_extra_nt_at_end = False):
    """
    """
    genome_seq = genome_fasta_dict[locus_tag]

    position_start, position_end = positions.split("..")

    if grab_extra_nt_at_end == True:
        position_end = int(position_end) + 1

    genome_snip = genome_seq[int(position_start) - 1 : int(position_end)]

    return genome_snip


def prot_msa(arfa_summary_loc, arfa_output_dir, desired_prefix, tree):
    """
    """
    acc_taxid_dict = {}
    """
    with open("results/prfb_asd_freescan.csv", "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            acc = row[0] # row[1].lower()
            taxonomy = row[1].lower()

            acc_taxid_dict[acc] = taxonomy
    f.close()
    """
    with open("results/ncbi-download-summary.txt", "r") as f:
        for row in f:
            taxonomy, acc = row.replace("\n","").split("\t")
            acc_taxid_dict[acc] = taxonomy.lower()
    f.close

    arfa_call_genomeacc = {}
    arfa_call_genomeacc["True"] = []
    arfa_call_genomeacc["False"] = []

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

    # mine out seqs from prfB FS = FALSE annotations
    false_nuc_fs_neighborhood_dict = {}
    for genome_acc in arfa_call_genomeacc["False"] + arfa_call_genomeacc["True"]:
        taxonomy = acc_taxid_dict[genome_acc]
        msa_fasta = "results/msa-prot-aligns/" + genome_acc + ".fa"
        with open(msa_fasta, "r") as f:
            seq = ""
            for row in f:
                if row.startswith(">"):
                    subj = row
                else:
                    seq = seq + row.replace("\n","")
        f.close()

        # extract portion of sequence depending on relative FS site
        fs_pos = 58
        extracted_seq_por = seq[fs_pos - 10 : fs_pos + 4]
        false_nuc_fs_neighborhood_dict[taxonomy] = extracted_seq_por

    # mine out protein sequences from prfB ARFA annotation
    nuc_seqid_dict = {}
    nuc_fs_neighborhood_dict = {}
    with open("results/prfB_TRUE_protf.fa", "w") as f:
        for genome_acc in arfa_call_genomeacc["True"]:
            arfa_output_loc = arfa_output_dir + genome_acc + "_arfa.txt"

            genome_fasta_dict = fasta_to_dict("/home/fung/Documents/codon_freqs/resources/ncbi-data/" + genome_acc + "_genomic.fna")

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

            prfb_sequence = prfb_sequence.replace("\n","").replace('"','').replace(" ","").replace("nohit!","")

            # get prfb nucleotide sequence
            first_set_pos = positions.split("join(")[1].split(",")[0]
            second_set_pos = positions.split(",")[1].split(")")[0]

            if "complement" in positions:
                num_bp_in_fs = int(second_set_pos.split("..")[1]) - int(second_set_pos.split("..")[0])
            else:
                num_bp_in_fs = int(first_set_pos.split("..")[1]) - int(first_set_pos.split("..")[0])

            nuc_seq_portion = extract_sequence_snip(locus_tag, genome_fasta_dict, first_set_pos, True) + \
                            extract_sequence_snip(locus_tag, genome_fasta_dict, second_set_pos, False)

            if "complement" in positions:
                dna = Seq(nuc_seq_portion)
                nuc_seq_portion = str(dna.reverse_complement())

            # get area around FS site for tree building/seq logos
            seq_por = nuc_seq_portion[num_bp_in_fs - 30 : num_bp_in_fs + 10]
            #nuc_fs_neighborhood_dict[genome_acc] = seq_por

            # get prot seq por
            seq_por = prfb_sequence[int(num_bp_in_fs / 3) - 10 : int(num_bp_in_fs / 3) + 4]
            nuc_fs_neighborhood_dict[genome_acc] = seq_por

            nuc_seqid_dict[genome_acc] = nuc_seq_portion
            

            f.write(">" + genome_acc + "\n")
            f.write(prfb_sequence + "\n")
    f.close()

    tax_seq_dict = {}
    with open("results/prfb_fs_neighorhood.fa", "w") as f:
        for acc in nuc_fs_neighborhood_dict:
            seq = nuc_fs_neighborhood_dict[acc]

            taxonomy = acc_taxid_dict[acc]

            tax_seq_dict[taxonomy] = seq

            f.write(">" + taxonomy + "\n")
            f.write(seq + "\n")
    f.close()

    # get ASD seqs
    all_asd_seqs = []
    tax_asd_dict = {}
    with open("results/prfb_asd_freescan.csv", "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            genome_acc, taxonomy, status, asd_seq, num_asd_seqs, sd_seq, msd_seq, sd_min_e, msd_min_e = row
            if len(asd_seq) > 0:
                all_asd_seqs.append(asd_seq)
                tax_asd_dict[taxonomy.lower()] = asd_seq
    f.close()

    # get all found families
    found_families = set()
    for taxonomy in tax_seq_dict:
        p_lineage = tree.ascend(taxonomy)

        for node in p_lineage:
            node_t = node.taxid

            if desired_prefix in node_t:
                found_families.add(node_t)
    
    # now write out FASTA's for each family neighborhood
    os.system("mkdir results/" + desired_prefix + "/")
    for family_t in found_families:
        family_asd_seqs = []
        with open("results/" + desired_prefix + "/" + family_t + "-fs_n.fa", "w") as f:
            for taxonomy in tax_seq_dict:
                seq = tax_seq_dict[taxonomy]

                p_lineage = tree.ascend(taxonomy)
                f_in_lineage = False
                for node in p_lineage:
                    node_t = node.taxid

                    if node_t == family_t:
                        f_in_lineage = True

                if f_in_lineage == True:
                    f.write(">" + taxonomy + "\n")
                    f.write(seq + "\n")

                    if taxonomy in tax_asd_dict:
                        family_asd_seqs.append(tax_asd_dict[taxonomy])
        f.close()

        sites_df = pd.read_csv("results/" + desired_prefix + "/" + family_t + "-fs_n.fa", comment='>', names=['site'])
        sites_list = sites_df['site'].values

        # write out SD region wordlogo
        fig = plt.figure()

        counts_df = logomaker.alignment_to_matrix(sequences=sites_list, to_type='counts')
        logomaker.Logo(counts_df, color_scheme='NajafabadiEtAl2017')

        fig.tight_layout()
        plt.title(family_t)
        plt.savefig("results/" + desired_prefix + "/" + family_t + "-seq_logo.png")
        plt.close()

        # write out corresponding ASD wordlogo
        fig = plt.figure()

        if len(family_asd_seqs) > 0:
            counts_df = logomaker.alignment_to_matrix(sequences=family_asd_seqs, to_type='counts')
            logomaker.Logo(counts_df, color_scheme='NajafabadiEtAl2017')

            fig.tight_layout()
            plt.title(family_t + " ASDs")
            plt.savefig("results/" + desired_prefix + "/" + family_t + "_ASD-seq_logo.png")
            plt.close()
        
        # get FS = False sequences per the MSA method
        num_in_lin = 0
        with open("results/" + desired_prefix + "/" + family_t + "-MSA_FS_FALSE.fa", "w") as f:
            for genome_acc in arfa_call_genomeacc["False"]:
                taxonomy = acc_taxid_dict[genome_acc]
                if taxonomy in false_nuc_fs_neighborhood_dict:
                    seq = false_nuc_fs_neighborhood_dict[taxonomy]

                    p_lineage = tree.ascend(taxonomy)
                    f_in_lineage = False
                    for node in p_lineage:
                        node_t = node.taxid

                        if node_t == family_t:
                            f_in_lineage = True

                    if f_in_lineage == True:
                        num_in_lin += 1
                        f.write(">" + taxonomy + "\n")
                        f.write(seq + "\n")
        f.close()

        if num_in_lin > 0:
            sites_df = pd.read_csv("results/" + desired_prefix + "/" + family_t + "-MSA_FS_FALSE.fa", comment='>', names=['site'])
            sites_list = sites_df['site'].values

            # write out SD region wordlogo
            fig = plt.figure()

            counts_df = logomaker.alignment_to_matrix(sequences=sites_list, to_type='counts')
            logomaker.Logo(counts_df, color_scheme='NajafabadiEtAl2017')

            fig.tight_layout()
            plt.title(family_t + " from MSA FS=FALSE")
            plt.savefig("results/" + desired_prefix + "/" + family_t + "-seq_logo-MSA_FS_FALSE.png")
            plt.close()

        # get FS = True sequences per the MSA method
        num_in_lin = 0
        with open("results/" + desired_prefix + "/" + family_t + "-MSA_FS_TRUE.fa", "w") as f:
            for genome_acc in arfa_call_genomeacc["True"]:
                taxonomy = acc_taxid_dict[genome_acc]
                if taxonomy in false_nuc_fs_neighborhood_dict:
                    seq = false_nuc_fs_neighborhood_dict[taxonomy]

                    p_lineage = tree.ascend(taxonomy)
                    f_in_lineage = False
                    for node in p_lineage:
                        node_t = node.taxid

                        if node_t == family_t:
                            f_in_lineage = True

                    if f_in_lineage == True:
                        num_in_lin += 1
                        f.write(">" + taxonomy + "\n")
                        f.write(seq + "\n")
        f.close()
        
        print("here!", num_in_lin)
        if num_in_lin > 0:
            sites_df = pd.read_csv("results/" + desired_prefix + "/" + family_t + "-MSA_FS_TRUE.fa", comment='>', names=['site'])
            sites_list = sites_df['site'].values

            # write out SD region wordlogo
            fig = plt.figure()

            counts_df = logomaker.alignment_to_matrix(sequences=sites_list, to_type='counts')
            logomaker.Logo(counts_df, color_scheme='NajafabadiEtAl2017')

            fig.tight_layout()
            plt.title(family_t + " from MSA FS=TRUE")
            plt.savefig("results/" + desired_prefix + "/" + family_t + "-seq_logo-MSA_FS_TRUE.png")
            plt.close()


    fig = plt.figure()
    counts_df = logomaker.alignment_to_matrix(sequences=all_asd_seqs, to_type='counts')
    logomaker.Logo(counts_df, color_scheme='NajafabadiEtAl2017')

    fig.tight_layout()
    plt.title("ASD Seqs - ALL")
    plt.savefig("results/" + "ASD-ALL-seg_logo.png")
    plt.close()

    return


def nuc_msa(arfa_summary_loc, arfa_output_dir, desired_prefix, tree):
    """
    """
    acc_taxid_dict = {}
    """
    with open("results/prfb_asd_freescan.csv", "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            acc = row[0] # row[1].lower()
            taxonomy = row[1].lower()

            acc_taxid_dict[acc] = taxonomy
    f.close()
    """
    with open("results/ncbi-download-summary.txt", "r") as f:
        for row in f:
            taxonomy, acc = row.replace("\n","").split("\t")
            acc_taxid_dict[acc] = taxonomy.lower()
    f.close

    arfa_call_genomeacc = {}
    arfa_call_genomeacc["True"] = []
    arfa_call_genomeacc["False"] = []

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

    # mine out seqs from prfB FS = FALSE annotations
    false_nuc_fs_neighborhood_dict = {}
    for genome_acc in arfa_call_genomeacc["False"] + arfa_call_genomeacc["True"]:
        taxonomy = acc_taxid_dict[genome_acc]
        msa_fasta = "results/msa-nuc-aligns/" + genome_acc + ".fa"
        if os.path.isfile(msa_fasta) == True:
            with open(msa_fasta, "r") as f:
                seq = ""
                for row in f:
                    if row.startswith(">"):
                        subj = row
                    else:
                        seq = seq + row.replace("\n","")
            f.close()

            # extract portion of sequence depending on relative FS site
            fs_pos = 177
            extracted_seq_por = seq[fs_pos - 30 : fs_pos + 10]
            false_nuc_fs_neighborhood_dict[taxonomy] = extracted_seq_por

    # mine out protein sequences from prfB ARFA annotation
    nuc_seqid_dict = {}
    nuc_fs_neighborhood_dict = {}
    with open("results/prfB_TRUE_protf.fa", "w") as f:
        for genome_acc in arfa_call_genomeacc["True"]:
            arfa_output_loc = arfa_output_dir + genome_acc + "_arfa.txt"

            genome_fasta_dict = fasta_to_dict("/home/fung/Documents/codon_freqs/resources/ncbi-data/" + genome_acc + "_genomic.fna")

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

            prfb_sequence = prfb_sequence.replace("\n","").replace('"','').replace(" ","").replace("nohit!","")

            # get prfb nucleotide sequence
            first_set_pos = positions.split("join(")[1].split(",")[0]
            second_set_pos = positions.split(",")[1].split(")")[0]

            if "complement" in positions:
                num_bp_in_fs = int(second_set_pos.split("..")[1]) - int(second_set_pos.split("..")[0])
            else:
                num_bp_in_fs = int(first_set_pos.split("..")[1]) - int(first_set_pos.split("..")[0])

            nuc_seq_portion = extract_sequence_snip(locus_tag, genome_fasta_dict, first_set_pos, True) + \
                            extract_sequence_snip(locus_tag, genome_fasta_dict, second_set_pos, False)

            if "complement" in positions:
                dna = Seq(nuc_seq_portion)
                nuc_seq_portion = str(dna.reverse_complement())

            # get area around FS site for tree building/seq logos
            seq_por = nuc_seq_portion[num_bp_in_fs - 30 : num_bp_in_fs + 10]
            nuc_fs_neighborhood_dict[genome_acc] = seq_por

            # get prot seq por
            seq_por = prfb_sequence[int(num_bp_in_fs / 3) - 10 : int(num_bp_in_fs / 3) + 4]
            #nuc_fs_neighborhood_dict[genome_acc] = seq_por

            nuc_seqid_dict[genome_acc] = nuc_seq_portion
            

            f.write(">" + genome_acc + "\n")
            f.write(prfb_sequence + "\n")
    f.close()

    tax_seq_dict = {}
    with open("results/prfb_fs_neighorhood.fa", "w") as f:
        for acc in nuc_fs_neighborhood_dict:
            seq = nuc_fs_neighborhood_dict[acc]

            taxonomy = acc_taxid_dict[acc]

            tax_seq_dict[taxonomy] = seq

            f.write(">" + taxonomy + "\n")
            f.write(seq + "\n")
    f.close()

    # get ASD seqs
    all_asd_seqs = []
    tax_asd_dict = {}
    with open("results/prfb_asd_freescan.csv", "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            genome_acc, taxonomy, status, asd_seq, num_asd_seqs, sd_seq, msd_seq, sd_min_e, msd_min_e = row
            if len(asd_seq) > 0:
                all_asd_seqs.append(asd_seq)
                tax_asd_dict[taxonomy.lower()] = asd_seq
    f.close()

    # get all found families
    found_families = set()
    for taxonomy in tax_seq_dict:
        p_lineage = tree.ascend(taxonomy)

        for node in p_lineage:
            node_t = node.taxid

            if desired_prefix in node_t:
                found_families.add(node_t)
    
    # now write out FASTA's for each family neighborhood
    os.system("mkdir results/" + desired_prefix + "/")
    for family_t in found_families:
        family_asd_seqs = []
        with open("results/" + desired_prefix + "/" + family_t + "-fs_n.fa", "w") as f:
            for taxonomy in tax_seq_dict:
                seq = tax_seq_dict[taxonomy]

                p_lineage = tree.ascend(taxonomy)
                f_in_lineage = False
                for node in p_lineage:
                    node_t = node.taxid

                    if node_t == family_t:
                        f_in_lineage = True

                if f_in_lineage == True:
                    f.write(">" + taxonomy + "\n")
                    f.write(seq + "\n")

                    if taxonomy in tax_asd_dict:
                        family_asd_seqs.append(tax_asd_dict[taxonomy])
        f.close()

        sites_df = pd.read_csv("results/" + desired_prefix + "/" + family_t + "-fs_n.fa", comment='>', names=['site'])
        sites_list = sites_df['site'].values

        # write out SD region wordlogo
        fig = plt.figure()

        counts_df = logomaker.alignment_to_matrix(sequences=sites_list, to_type='counts')
        logomaker.Logo(counts_df, color_scheme='NajafabadiEtAl2017')

        fig.tight_layout()
        plt.title(family_t)
        plt.savefig("results/" + desired_prefix + "/" + family_t + "-seq_logo.png")
        plt.close()

        # write out corresponding ASD wordlogo
        fig = plt.figure()

        if len(family_asd_seqs) > 0:
            counts_df = logomaker.alignment_to_matrix(sequences=family_asd_seqs, to_type='counts')
            logomaker.Logo(counts_df, color_scheme='NajafabadiEtAl2017')

            fig.tight_layout()
            plt.title(family_t + " ASDs")
            plt.savefig("results/" + desired_prefix + "/" + family_t + "_ASD-seq_logo.png")
            plt.close()
        
        # get FS = False sequences per the MSA method
        num_in_lin = 0
        gutcheck = 0
        with open("results/" + desired_prefix + "/" + family_t + "-MSA_FS_FALSE.fa", "w") as f:
            for genome_acc in arfa_call_genomeacc["False"]:
                taxonomy = acc_taxid_dict[genome_acc]

                p_lineage = tree.ascend(taxonomy)
                f_in_lineage = False
                for node in p_lineage:
                    node_t = node.taxid

                    if node_t == family_t:
                        f_in_lineage = True

                if f_in_lineage == True:
                    gutcheck += 1

                if taxonomy in false_nuc_fs_neighborhood_dict:
                    seq = false_nuc_fs_neighborhood_dict[taxonomy]

                    p_lineage = tree.ascend(taxonomy)
                    f_in_lineage = False
                    for node in p_lineage:
                        node_t = node.taxid

                        if node_t == family_t:
                            f_in_lineage = True

                    if f_in_lineage == True:
                        num_in_lin += 1
                        f.write(">" + taxonomy + "\n")
                        f.write(seq + "\n")
        f.close()

        print(family_t, gutcheck)

        if num_in_lin > 0:
            sites_df = pd.read_csv("results/" + desired_prefix + "/" + family_t + "-MSA_FS_FALSE.fa", comment='>', names=['site'])
            sites_list = sites_df['site'].values

            # write out SD region wordlogo
            fig = plt.figure()

            counts_df = logomaker.alignment_to_matrix(sequences=sites_list, to_type='counts')
            logomaker.Logo(counts_df, color_scheme='NajafabadiEtAl2017')

            fig.tight_layout()
            plt.title(family_t + " from MSA FS=FALSE")
            plt.savefig("results/" + desired_prefix + "/" + family_t + "-seq_logo-MSA_FS_FALSE.png")
            plt.close()

        # get FS = True sequences per the MSA method
        num_in_lin = 0
        with open("results/" + desired_prefix + "/" + family_t + "-MSA_FS_TRUE.fa", "w") as f:
            for genome_acc in arfa_call_genomeacc["True"]:
                taxonomy = acc_taxid_dict[genome_acc]
                if taxonomy in false_nuc_fs_neighborhood_dict:
                    seq = false_nuc_fs_neighborhood_dict[taxonomy]

                    p_lineage = tree.ascend(taxonomy)
                    f_in_lineage = False
                    for node in p_lineage:
                        node_t = node.taxid

                        if node_t == family_t:
                            f_in_lineage = True

                    if f_in_lineage == True:
                        num_in_lin += 1
                        f.write(">" + taxonomy + "\n")
                        f.write(seq + "\n")
        f.close()
        
        print("here!", num_in_lin)
        if num_in_lin > 0:
            sites_df = pd.read_csv("results/" + desired_prefix + "/" + family_t + "-MSA_FS_TRUE.fa", comment='>', names=['site'])
            sites_list = sites_df['site'].values

            # write out SD region wordlogo
            fig = plt.figure()

            counts_df = logomaker.alignment_to_matrix(sequences=sites_list, to_type='counts')
            logomaker.Logo(counts_df, color_scheme='NajafabadiEtAl2017')

            fig.tight_layout()
            plt.title(family_t + " from MSA FS=TRUE")
            plt.savefig("results/" + desired_prefix + "/" + family_t + "-seq_logo-MSA_FS_TRUE.png")
            plt.close()


    fig = plt.figure()
    counts_df = logomaker.alignment_to_matrix(sequences=all_asd_seqs, to_type='counts')
    logomaker.Logo(counts_df, color_scheme='NajafabadiEtAl2017')

    fig.tight_layout()
    plt.title("ASD Seqs - ALL")
    plt.savefig("results/" + "ASD-ALL-seg_logo.png")
    plt.close()

    return




def mine_seqs_for_msa(arfa_summary_loc, arfa_output_dir, tree):
    """
    """
    desired_prefix = "o__"

    #prot_msa(arfa_summary_loc, arfa_output_dir, desired_prefix, tree)
    nuc_msa(arfa_summary_loc, arfa_output_dir, desired_prefix, tree)

    raise Exception

    with open("results/prfB_TRUE_nuc.fa", "w") as f:
        for genome_acc in nuc_seqid_dict:
            nuc_seq_portion = nuc_seqid_dict[genome_acc]

            f.write(">" + genome_acc + "\n")
            f.write(nuc_seq_portion + "\n")
    f.close()

    nuc_seqid_dict = {}
    with open("results/prfB_FALSE_prot.fa", "w") as f:
        for genome_acc in arfa_call_genomeacc["False"]:
            arfa_output_loc = arfa_output_dir + genome_acc + "_arfa.txt"

            genome_fasta_dict = fasta_to_dict("/home/fung/Documents/codon_freqs/resources/ncbi-data/" + genome_acc + "_genomic.fna")

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

            if "join" not in positions:
                # get prfb nucleotide sequence
                if "complement" in positions:
                    position_temp = positions.split("(")[1].split(")")[0]
                else:
                    position_temp = positions

                nuc_seq_portion = extract_sequence_snip(locus_tag, genome_fasta_dict, position_temp)

                if "complement" in positions:
                    dna = Seq(nuc_seq_portion)
                    nuc_seq_portion = str(dna.reverse_complement())

                nuc_seqid_dict[genome_acc] = nuc_seq_portion
            prfb_sequence = prfb_sequence.replace("\n","").replace('"','').replace(" ","").replace("nohit!","")

            f.write(">" + genome_acc + "\n")
            f.write(prfb_sequence + "\n")
    f.close()

    with open("results/prfB_FALSE_nuc.fa", "w") as f:
        for genome_acc in nuc_seqid_dict:
            nuc_seq_portion = nuc_seqid_dict[genome_acc]

            f.write(">" + genome_acc + "\n")
            f.write(nuc_seq_portion + "\n")
    f.close()

    return