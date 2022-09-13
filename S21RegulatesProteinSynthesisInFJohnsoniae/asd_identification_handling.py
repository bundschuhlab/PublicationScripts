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
from gtdb_genome_handling import download_from_assembly_ftp
from utils import fasta_to_dict


def get_list_downloaded_assembly(genome_assembly_ncbi_data_loc):
    """
    """
    downloaded_assembly_accs = []
    with open(genome_assembly_ncbi_data_loc, "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            assembly_acc, taxonomy_name, assembly_type, ftp_address, success = row

            if success == "True" and assembly_type == "refseq":
                downloaded_assembly_accs.append(assembly_acc)
    f.close()

    return downloaded_assembly_accs


def get_16s_positions(barrnap_annotation_loc):
    """
    """
    barrnap_16s_predictions = []
    with open(barrnap_annotation_loc, "r") as f:
        for row in f:
            if row.startswith(">"):
                seqid = row[1:].replace("\n","")
                if "16S_rRNA" in seqid:
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


def run_cutadapt(cutadapt_exe, barrnap_seqment_loc, cutadapt_out_loc):
    """
    """
    command = cutadapt_exe + \
                            " -g AAGTCGTAACAAGGTAGCCGT" + \
                            " -e 0.25 -O 21" + \
                            " --discard-untrimmed" + \
                            " -o " + cutadapt_out_loc + \
                            " " + barrnap_seqment_loc + \
                            " -j 0 --quiet"
    os.system(command)

    return 


def run_barrnap_pipeline(barrnap_exe, cutadapt_exe, \
                            assembly_bio_data_loc, assembly_acc, \
                            barrnap_num_downstream, \
                            barrnap_num_upstream_from_end):
    """
    """
    genome_loc = assembly_bio_data_loc + assembly_acc + "_genomic.fna"
    barrnap_seq_loc = assembly_bio_data_loc + assembly_acc + "_barrnap.fa"
    barrnap_seqment_loc = assembly_bio_data_loc + assembly_acc + "_16S-segments.fa"
    cutadapt_out_loc = assembly_bio_data_loc + assembly_acc + "_cutadapt.fa"
    # TEMP!!
    if os.path.isfile(assembly_bio_data_loc + assembly_acc + "asd-predictions.fa") == True:
        os.system("rm " + assembly_bio_data_loc + assembly_acc + "asd-predictions.fa")
    asd_prediction_loc = assembly_bio_data_loc + assembly_acc + "_asd-predictions.fa"

    try:
        # run barrnap
        command = barrnap_exe + \
                    " " + genome_loc + \
                    " --outseq " + barrnap_seq_loc + \
                    " --quiet"
        os.system(command)

        # extract 16S positions
        barrnap_16s_predictions = get_16s_positions(barrnap_seq_loc)

        # get full genome dictionary
        genome_seq_dict = fasta_to_dict(genome_loc)

        # extract snips from around 16S end
        write_out_16s_asdg_area(barrnap_seqment_loc, barrnap_16s_predictions, \
                                    barrnap_num_downstream, \
                                    barrnap_num_upstream_from_end, \
                                    genome_seq_dict)

        # run cutadapt over snips
        run_cutadapt(cutadapt_exe, barrnap_seqment_loc, cutadapt_out_loc)

        # extract snip areas of interest to grab ASD
        downselect_cutadapt_trimmed(cutadapt_out_loc, asd_prediction_loc)

        return True, asd_prediction_loc

    except:
        return False, asd_prediction_loc


def identify_asd_seq(asd_prediction_loc):
    """
    """
    asd_seqid_dict = fasta_to_dict(asd_prediction_loc)

    uniq_seqs = set()
    for seqid in asd_seqid_dict:
        seq = asd_seqid_dict[seqid]
        uniq_seqs.add(seq)

    return uniq_seqs


def extract_asd_seqs(genome_assembly_ncbi_data_loc, asd_identification_loc, \
                        bio_data_loc, barrnap_exe, cutadapt_exe, \
                            barrnap_num_downstream, \
                            barrnap_num_upstream_from_end):
    """
    """
    downloaded_assembly_accs = get_list_downloaded_assembly(genome_assembly_ncbi_data_loc)

    with open(asd_identification_loc, "w") as f:
        out = csv.writer(f)
        out.writerow(["Assembly Acc", "Barrnap/Cutadapt Successful?", \
                        "Number of Uniq ASD Seqs", "ASD Seq"])

        for assembly_acc in downloaded_assembly_accs:
            assembly_bio_data_loc = bio_data_loc + assembly_acc + "/"
            
            print(assembly_acc, downloaded_assembly_accs.index(assembly_acc), "/", len(downloaded_assembly_accs))
            status_call, asd_prediction_loc = \
                    run_barrnap_pipeline(barrnap_exe, cutadapt_exe, \
                                            assembly_bio_data_loc, assembly_acc, \
                                            barrnap_num_downstream, \
                                            barrnap_num_upstream_from_end)

            if status_call == False:
                asd_seq = ""
                num_asd_seqs = "NA"

            elif status_call == True:
                uniq_asd_seqs = identify_asd_seq(asd_prediction_loc)
                num_asd_seqs = len(uniq_asd_seqs)

                if num_asd_seqs == 1:
                    for s in uniq_asd_seqs:
                        asd_seq = s
                else:
                    asd_seq = ""
            
            out.writerow([assembly_acc, status_call, num_asd_seqs, asd_seq])
    f.close()

    return