import os
import csv
from Bio.Seq import Seq
from utils import fasta_to_dict


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


def get_assembly_accs_of_interest(pipeline_output_loc, tree):
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
            is_flavo = False
            for node in parent_nodes:
                node_taxonomy = node.taxid
                if node_taxonomy == "o__flavobacteriales":
                    is_flavo = True
            
            if is_flavo == True:
                assembly_accs_to_grab.append(assembly_acc)
    f.close()

    return assembly_accs_to_grab


def get_gene_priors_to_fasta(genome_loc, gff_loc, assembly_gene_prior_fa_loc, illegal_letters):
    """
    """
    genome_seqid_dict = fasta_to_dict(genome_loc)

    counter = 0
    counter_feature_dict = {}
    with open(assembly_gene_prior_fa_loc, "w") as f:
        with open(gff_loc, "r") as r:
            for row in r:
                if row.startswith("#"):
                    continue
                # feature row
                row = row.replace("\n","").split("\t")

                annotation_type = row[2]
                if annotation_type == "CDS":
                    genome_seqid = row[0]
                    genome_seq = genome_seqid_dict[genome_seqid]
                    position_start = int(row[3])
                    position_end = int(row[4])
                    strand_direction = row[6]

                    feature_ids = row[8]
                    feature_id_substring = feature_ids.split("ID=")[1].split(";")[0]

                    # positive strand handling
                    if strand_direction == "+":
                        seq_por = genome_seq[position_start - 55 : position_start - 1].upper()
                        # if fails feature conditions, move on to loop
                        if position_start < 56 or any((nt in illegal_letters) for nt in seq_por):
                            continue

                    # negative strand handling
                    if strand_direction == "-":
                        seq_por = genome_seq[position_end : position_end + 54].upper()
                        # if fails feature conditions, move on to loop
                        if (position_end > len(genome_seq) - 56) or any((nt in illegal_letters) for nt in seq_por):
                            continue

                        dna = Seq(seq_por)
                        seq_por = str(dna.reverse_complement())

                    # write out sequences
                    counter += 1
                    f.write(">" + str(counter) + "\n")
                    f.write(seq_por + "\n")
                    counter_feature_dict[counter] = feature_ids
        r.close()
    f.close()

    return counter_feature_dict


def split_into_sd_msd(assembly_all_gene_sd_loc, assembly_all_gene_msd_loc, assembly_gene_prior_fa_loc):
    """
    """
    seqid_dict = fasta_to_dict(assembly_gene_prior_fa_loc)

    with open(assembly_all_gene_sd_loc, "w") as f:
        for seqid in seqid_dict:
            seq = seqid_dict[seqid]

            f.write(">" + seqid + "||" + "TIR_BIND" + "\n")
            f.write(seq[25:] + "\n")
    f.close()

    with open(assembly_all_gene_msd_loc, "w") as f:
        for seqid in seqid_dict:
            seq = seqid_dict[seqid]

            f.write(">" + seqid + "||" + "TIR_BIND_MSD" + "\n")
            f.write(seq[:29] + "\n")
    f.close()

    return


def run_freescan(freescan_exe, tir_fasta, asd_seq, output_fldr):
    """
    """
    temp_fa_loc = "resources/temp/temp.fa"

    os.system("rm -rf " + output_fldr)
    os.system("mkdir " + output_fldr)

    seqid_dict = fasta_to_dict(tir_fasta)
    product_num_dict = {}
    for seqid in seqid_dict:
        counter, _ = seqid.split("||")
        seq = seqid_dict[seqid]

        with open(temp_fa_loc, "w") as f:
            f.write(">" + seqid + "\n")
            f.write(seq + "\n")
        f.close()

        command = freescan_exe + " -q " + asd_seq[::-1] + \
                                " " + temp_fa_loc + \
                                " > " + output_fldr + str(counter) + ".out"
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


def run_freescan_on_all_genes(flavo_gene_freescan_loc, pipeline_output_loc, \
                                illegal_letters, bio_data_loc, \
                                asd_identification_loc, tree, \
                                freescan_exe):
    """
    """
    assembly_accs_to_grab = get_assembly_accs_of_interest(pipeline_output_loc, tree)
    assembly_acc_asd_dict = get_asd_seqs(asd_identification_loc)

    # get list of all gene annotations
    for assembly_acc in assembly_accs_to_grab:
        print(assembly_acc, assembly_accs_to_grab.index(assembly_acc), "/", len(assembly_accs_to_grab))
        assembly_bio_data_loc = bio_data_loc + assembly_acc + "/"

        genome_loc = assembly_bio_data_loc + assembly_acc + "_genomic.fna"
        gff_loc = assembly_bio_data_loc + assembly_acc + "_genomic.gff"

        assembly_gene_prior_fa_loc = assembly_bio_data_loc + "all_genes_priors.fa"
        counter_feature_dict = get_gene_priors_to_fasta(genome_loc, gff_loc, assembly_gene_prior_fa_loc, illegal_letters)

        assembly_all_gene_loc = assembly_bio_data_loc + "all_genes/"
        os.system("mkdir " + assembly_all_gene_loc)
        assembly_all_gene_sd_loc = assembly_all_gene_loc + "all_gene_SD_tirs.fa"
        assembly_all_gene_msd_loc = assembly_all_gene_loc + "all_gene_MSD_tirs.fa"
        split_into_sd_msd(assembly_all_gene_sd_loc, assembly_all_gene_msd_loc, assembly_gene_prior_fa_loc)

        asd_seq = assembly_acc_asd_dict[assembly_acc]

        assembly_all_gene_sd_folder = assembly_all_gene_loc + "all_gene_SDs/"
        assembly_all_gene_msd_folder = assembly_all_gene_loc + "all_gene_MSDs/"

        run_freescan(freescan_exe, assembly_all_gene_sd_loc, asd_seq, assembly_all_gene_sd_folder)
        run_freescan(freescan_exe, assembly_all_gene_msd_loc, asd_seq, assembly_all_gene_msd_folder)

        with open(flavo_gene_freescan_loc + assembly_acc + ".csv", "w") as f:
            out = csv.writer(f)
            out.writerow(["CDS Counter", "Feature String", "ASD Seq", \
                                "SD TIR", "MSD TIR"])

            for counter in counter_feature_dict:
                feature_string = counter_feature_dict[counter]

                # get MSD and SD freescan minimum TIRs
                sd_loc = assembly_all_gene_sd_folder + str(counter) + ".out"
                msd_loc = assembly_all_gene_msd_folder + str(counter) + ".out"

                sd_min_E = get_min_delG_from_freescan_output(sd_loc)
                msd_min_E = get_min_delG_from_freescan_output(msd_loc)

                out.writerow([counter, feature_string, asd_seq, sd_min_E, msd_min_E])
        f.close()

    return
