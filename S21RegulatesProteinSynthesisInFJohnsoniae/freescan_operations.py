import csv
from itertools import product
import os
from utils import fasta_to_dict


def get_list_annotation_genomes(s_l_annotation_loc, min_number_s_l_annotations):
    """
    """
    assembly_accs_w_annotations = []
    with open(s_l_annotation_loc, "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            assembly_acc, number_uniq_s_prot_annotations, \
                            number_uniq_l_prot_annotations = row
            
            tot_num = int(number_uniq_s_prot_annotations) + int(number_uniq_l_prot_annotations)
            if tot_num >= min_number_s_l_annotations:
                assembly_accs_w_annotations.append(assembly_acc)
    f.close()

    return assembly_accs_w_annotations


def get_asd_seqs(asd_identification_loc):
    """
    """
    assembly_acc_asd_dict = {}
    with open(asd_identification_loc, "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            assembly_acc, status_call, num_asd_seqs, asd_seq = row
            assembly_acc_asd_dict[assembly_acc] = asd_seq
    f.close()

    return assembly_acc_asd_dict


def split_into_bind_msd_bind(prot_fasta_locs, prot_tir_bind, prot_tir_bind_msd):
    """
    """
    # write the first # bps to TIR_Bind
    with open(prot_tir_bind, "w") as f:
        for filename in os.listdir(prot_fasta_locs):
            seqid_dict = fasta_to_dict(prot_fasta_locs + filename)

            for seqid in seqid_dict:
                seq = seqid_dict[seqid]

                f.write(">" + seqid + "||" + "TIR_BIND" + "\n")
                f.write(seq[25:] + "\n")
    f.close()

    # write the first # bps to TIR_Bind_MSD
    with open(prot_tir_bind_msd, "w") as f:
        for filename in os.listdir(prot_fasta_locs):
            seqid_dict = fasta_to_dict(prot_fasta_locs + filename)

            for seqid in seqid_dict:
                seq = seqid_dict[seqid]

                f.write(">" + seqid + "||" + "TIR_BIND_MSD" + "\n")
                f.write(seq[:29] + "\n")
    f.close()    

    return


def run_freescan(freescan_exe, tir_fasta, asd_seq, prot_loc, sd_type):
    """
    """
    temp_fa_loc = "resources/temp/temp.fa"

    output_fldr = prot_loc + sd_type + "/"
    os.system("rm -rf " + output_fldr)
    os.system("mkdir " + output_fldr)

    seqid_dict = fasta_to_dict(tir_fasta)
    product_num_dict = {}
    for seqid in seqid_dict:
        _, _, product_num, _ = seqid.split("||")
        seq = seqid_dict[seqid]

        if product_num not in product_num_dict:
            product_num_dict[product_num] = 0
        product_num_dict[product_num] += 1

        with open(temp_fa_loc, "w") as f:
            f.write(">" + seqid + "\n")
            f.write(seq + "\n")
        f.close()

        command = freescan_exe + " -q " + asd_seq[::-1] + \
                                " " + temp_fa_loc + \
                                " > " + output_fldr + str(product_num) + "-" + str(product_num_dict[product_num]) + ".out"
        os.system(command)

    return


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


def get_min_tir_array(s_prot_loc, l_prot_loc):
    """
    """
    prot_mine_dict = {}
    # parse over s/l protein freescan output
    for filename in os.listdir(s_prot_loc + "tir_bind/"):
        filename_label = filename.split("-")[0]
        freescan_min_e = get_min_delG_from_freescan_output(s_prot_loc + "tir_bind/" + filename)

        if filename_label not in prot_mine_dict:
            prot_mine_dict[filename_label] = 0.0
        if freescan_min_e < prot_mine_dict[filename_label]:
            prot_mine_dict[filename_label] = freescan_min_e

    # note for L proteins, combine L7/L12
    for filename in os.listdir(l_prot_loc + "tir_bind/"):
        filename_label = filename.split("-")[0]
        if filename_label == "L7" or filename_label == "L12":
            filename_label = "L7_12"

        freescan_min_e = get_min_delG_from_freescan_output(l_prot_loc + "tir_bind/" + filename)

        if filename_label not in prot_mine_dict:
            prot_mine_dict[filename_label] = 0.0
        if freescan_min_e < prot_mine_dict[filename_label]:
            prot_mine_dict[filename_label] = freescan_min_e

    # determine if S21 is min E
    if "S21" not in prot_mine_dict:
        s21_energy = "NA"
        s21_status_call = "NA"
    else:
        s21_energy = prot_mine_dict["S21"]
        global_min_e = min(prot_mine_dict.values())

        if s21_energy == global_min_e:
            s21_status_call = "S21_Winner"
        else:
            s21_status_call = "S21_Loser"

    return prot_mine_dict, s21_energy, s21_status_call


def get_prot_labels():
    """
    """
    prot_labels = []
    # S proteins
    for x in range(1, 22, 1):
        label = "S" + str(x)
        prot_labels.append(label)
    # L proteins
    for x in range(1, 37, 1):
        # note - L7/L12 are equivalent 
        if x == 12:
            continue
        elif x == 7:
            label = "L7_12"
        else:
            label = "L" + str(x)
        prot_labels.append(label)

    return prot_labels


def run_freescan_pipeline(s_l_annotation_loc, pipeline_output_loc, \
                            asd_identification_loc, \
                            min_number_s_l_annotations, \
                            bio_data_loc, freescan_exe, \
                            assembly_acc_name_dict):
    """
    """
    assembly_accs_w_annotations = get_list_annotation_genomes(s_l_annotation_loc, min_number_s_l_annotations)
    assembly_acc_asd_dict = get_asd_seqs(asd_identification_loc)
    prot_labels = get_prot_labels()

    with open(pipeline_output_loc, "w") as f:
        out = csv.writer(f)
        header_row = ["Assembly Acc", "Taxonomy Name", "S21 Energy", "S21 Winner?"]
        for prot_label in prot_labels:
            header_row.append(prot_label)
        out.writerow(header_row)

        for assembly_acc in assembly_accs_w_annotations:
            print(assembly_acc, assembly_accs_w_annotations.index(assembly_acc), "/", len(assembly_accs_w_annotations))
            taxonomy_name = assembly_acc_name_dict[assembly_acc]
            assembly_bio_data_loc = bio_data_loc + assembly_acc + "/"

            s_prot_loc = assembly_bio_data_loc + "s_proteins/"
            s_prot_fasta_locs = s_prot_loc + "fastas/"
            s_prot_tir_bind = s_prot_loc + "TIR_Bind.fa"
            s_prot_tir_bind_msd = s_prot_loc + "TIR_Bind_MSD.fa"
            split_into_bind_msd_bind(s_prot_fasta_locs, s_prot_tir_bind, s_prot_tir_bind_msd)

            l_prot_loc = assembly_bio_data_loc + "l_proteins/"
            l_prot_fasta_locs = l_prot_loc + "fastas/"
            l_prot_tir_bind = l_prot_loc + "TIR_Bind.fa"
            l_prot_tir_bind_msd = l_prot_loc + "TIR_Bind_MSD.fa"
            split_into_bind_msd_bind(l_prot_fasta_locs, l_prot_tir_bind, l_prot_tir_bind_msd)

            asd_seq = assembly_acc_asd_dict[assembly_acc]

            run_freescan(freescan_exe, s_prot_tir_bind, asd_seq, s_prot_loc, "tir_bind")
            run_freescan(freescan_exe, s_prot_tir_bind_msd, asd_seq, s_prot_loc, "tir_bind_MSD")
            run_freescan(freescan_exe, l_prot_tir_bind, asd_seq, l_prot_loc, "tir_bind")
            run_freescan(freescan_exe, l_prot_tir_bind_msd, asd_seq, l_prot_loc, "tir_bind_MSD")

            prot_mine_dict, s21_energy, s21_status_call = get_min_tir_array(s_prot_loc, l_prot_loc)
            out_row = [assembly_acc, taxonomy_name, s21_energy, s21_status_call]
            for prot_label in prot_labels:
                if prot_label in prot_mine_dict:
                    energy = prot_mine_dict[prot_label]
                else:
                    energy = ""
                out_row.append(energy)
            out.writerow(out_row)
    f.close()


    return