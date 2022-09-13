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
from Bio.Seq import Seq
from utils import fasta_to_dict


def extract_s_proteins(genome_seqid_dict, gff_loc, \
                            illegal_letters, \
                            assembly_acc, \
                            s_prot_fasta_locs):
    """
    """
    uniq_prot_numbers = set()
    # parse over all S protein #s
    for x in range(1, 22, 1):
        with open(s_prot_fasta_locs + assembly_acc + "-S" + str(x) + ".fa", "w") as f:
            with open(gff_loc, "r") as r:
                for row in r:
                    if row.startswith("#"):
                        continue
                    # feature row
                    row = row.replace("\n","").split("\t")

                    genome_seqid = row[0]
                    genome_seq = genome_seqid_dict[genome_seqid]
                    position_start = int(row[3])
                    position_end = int(row[4])
                    strand_direction = row[6]

                    feature_ids = row[8]
                    feature_id_substring = feature_ids.split("ID=")[1].split(";")[0]
 
                    if('product=ribosomal protein S' + str(x) + ';' in feature_ids or 'product=30S ribosomal protein S' + str(x) + ';' in feature_ids or 'product=small subunit ribosomal protein S' + str(x) + ';' in feature_ids or 'product=SSU ribosomal protein S' + str(x) + 'P;' in feature_ids or 'product=Ribosomal protein S' + str(x) + ';' in feature_ids or 'product=SSU ribosomal protein S' + str(x) + 'p;' in feature_ids or 'product=S' + str(x) + 'p: ribosomal protein S' + str(x) + ';' in feature_ids or 'product=30S ribosomal subunit S' + str(x) + ';' in feature_ids or 'product=SSU ribosomal protein S' + str(x) + 'p RspU;' in feature_ids or 'product=putative ribosomal protein S' + str(x) + ';' in feature_ids or 'product=SSU ribosomal protein s' + str(x) + 'p;' in feature_ids or 'product=30S ribosomal protein S' + str(x) + ' {ECO:0000255|HAMAP-Rule:MF_00358};' in feature_ids or 'product=putative 30S ribosomal protein S' + str(x) + ';' in feature_ids or 'product=S' + str(x) + ': ribosomal protein S' + str(x) + ';' in feature_ids or 'product=30S ribosomal subunit protein S' + str(x) + ';' in feature_ids or 'product=SSU ribosomal protein S' + str(x) + 'p / SSU ribosomal protein S' + str(x) + 'p%2C zinc-independent;' in feature_ids or 'product=30S ribosomal protein S' + str(x) + ' {ECO:0000255|HAMAP-Rule:MF_00270};' in feature_ids or 'product=ribosomal protein S' + str(x) + ' (BS' + str(x) + ')' in feature_ids):
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
                        f.write(">" + feature_id_substring + "||" + strand_direction + "||" + "S" + str(x) + "\n")
                        f.write(seq_por + "\n")
                        uniq_prot_numbers.add(x)
            r.close()
        f.close()

    return len(uniq_prot_numbers)


def extract_l_proteins(genome_seqid_dict, gff_loc, \
                            illegal_letters, \
                            assembly_acc, \
                            l_prot_fasta_locs):
    """
    """
    uniq_prot_numbers = set()
    # parse over all S protein #s
    for x in range(1, 37, 1):
        with open(l_prot_fasta_locs + assembly_acc + "-L" + str(x) + ".fa", "w") as f:
            with open(gff_loc, "r") as r:
                for row in r:
                    if row.startswith("#"):
                        continue
                    # feature row
                    row = row.replace("\n","").split("\t")

                    genome_seqid = row[0]
                    genome_seq = genome_seqid_dict[genome_seqid]
                    position_start = int(row[3])
                    position_end = int(row[4])
                    strand_direction = row[6]

                    feature_ids = row[8]
                    feature_id_substring = feature_ids.split("ID=")[1].split(";")[0]
 
                    if('product=ribosomal protein L' + str(x) + ';' in feature_ids or 'product=50S ribosomal protein L' + str(x) + ';' in feature_ids or 'product=large subunit ribosomal protein L' + str(x) + ';' in feature_ids or 'product=LSU ribosomal protein L' + str(x) + 'P;' in feature_ids or 'product=Ribosomal protein L' + str(x) + ';' in feature_ids or 'product=LSU ribosomal protein L' + str(x) + 'p;' in feature_ids or 'product=L' + str(x) + 'p: ribosomal protein L' + str(x) + ';' in feature_ids or 'product=50S ribosomal subunit L' + str(x) + ';' in feature_ids or 'product=LSU ribosomal protein L' + str(x) + 'p RspU;' in feature_ids or 'product=putative ribosomal protein L' + str(x) + ';' in feature_ids or 'product=LSU ribosomal protein l' + str(x) + 'p;' in feature_ids or 'product=50S ribosomal protein L' + str(x) + ' {ECO:0000255|HAMAP-Rule:MF_00358};' in feature_ids or 'product=putative 50S ribosomal protein L' + str(x) + ';' in feature_ids or 'product=L' + str(x) + ': ribosomal protein L' + str(x) + ';' in feature_ids or 'product=50S ribosomal subunit protein L' + str(x) + ';' in feature_ids or 'product=LSU ribosomal protein L' + str(x) + 'p / LSU ribosomal protein L' + str(x) + 'p%2C zinc-independent;' in feature_ids or 'product=50S ribosomal protein L' + str(x) + ' {ECO:0000255|HAMAP-Rule:MF_00270};' in feature_ids or 'product=ribosomal protein L' + str(x) + ' (BL' + str(x) + ')' in feature_ids):
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
                        f.write(">" + feature_id_substring + "||" + strand_direction + "||" + "L" + str(x) + "\n")
                        f.write(seq_por + "\n")
                        uniq_prot_numbers.add(x)
            r.close()
        f.close()

    return len(uniq_prot_numbers)


def get_list_asd_genomes(asd_identification_loc):
    """
    """
    assembly_accs_w_asd = []
    with open(asd_identification_loc, "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            assembly_acc, status_call, num_asd_seqs, asd_seq = row

            if status_call == "True" and int(num_asd_seqs) == 1:
                assembly_accs_w_asd.append(assembly_acc)
    f.close()

    return assembly_accs_w_asd


def extract_s_l_annotations(asd_identification_loc, bio_data_loc, \
                                s_l_annotation_loc, illegal_letters):
    """
    """
    assembly_accs_w_asd = get_list_asd_genomes(asd_identification_loc)

    with open(s_l_annotation_loc, "w") as f:
        out = csv.writer(f)
        out.writerow(["Assembly Acc", "Number S Protein Annotations", \
                        "Number L Protein Annotations"])

        for assembly_acc in assembly_accs_w_asd:
            print(assembly_acc, assembly_accs_w_asd.index(assembly_acc), "/", len(assembly_accs_w_asd))
            assembly_bio_data_loc = bio_data_loc + assembly_acc + "/"
            genome_loc = assembly_bio_data_loc + assembly_acc + "_genomic.fna"

            ###### TEMP fix gbff --> gff
            if os.path.isfile(assembly_bio_data_loc + assembly_acc + "_genomic.gbff") == True:
                os.system("mv " + assembly_bio_data_loc + assembly_acc + "_genomic.gbff" + \
                                " " + assembly_bio_data_loc + assembly_acc + "_genomic.gff")

            gff_loc = assembly_bio_data_loc + assembly_acc + "_genomic.gff"

            genome_seqid_dict = fasta_to_dict(genome_loc)

            # S proteins
            s_prot_loc = assembly_bio_data_loc + "s_proteins/"
            os.system("rm -rf " + s_prot_loc)
            os.system("mkdir " + s_prot_loc)
            s_prot_fasta_locs = s_prot_loc + "fastas/"
            os.system("rm -rf " + s_prot_fasta_locs)
            os.system("mkdir " + s_prot_fasta_locs)
            number_uniq_s_prot_annotations = extract_s_proteins(genome_seqid_dict, gff_loc, \
                                                                    illegal_letters, \
                                                                    assembly_acc, \
                                                                    s_prot_fasta_locs)

            # L proteins
            l_prot_loc = assembly_bio_data_loc + "l_proteins/"
            os.system("rm -rf " + l_prot_loc)
            os.system("mkdir " + l_prot_loc)
            l_prot_fasta_locs = l_prot_loc + "fastas/"
            os.system("rm -rf " + l_prot_fasta_locs)
            os.system("mkdir " + l_prot_fasta_locs)
            number_uniq_l_prot_annotations = extract_l_proteins(genome_seqid_dict, gff_loc, \
                                                                    illegal_letters, \
                                                                    assembly_acc, \
                                                                    l_prot_fasta_locs)

            out.writerow([assembly_acc, number_uniq_s_prot_annotations, \
                            number_uniq_l_prot_annotations])
    f.close()

    return