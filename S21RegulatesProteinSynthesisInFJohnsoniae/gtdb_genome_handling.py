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
import time


def find_gtdb_species_genomes(gtdb_clusters, gtdb_bact_tsv):
    """
    """
    genomes_that_are_complete = []
    with open(gtdb_bact_tsv, "r") as f:
        next(f)
        for row in f:
            row = row.replace("\n","").split("\t")
            gtdb_assembly_acc = row[0]
            ncbi_genome_annotation = row[45]

            if ncbi_genome_annotation == "Complete Genome":
                genomes_that_are_complete.append(gtdb_assembly_acc)
    f.close()

    assembly_acc_name_dict = {}
    with open(gtdb_clusters, "r") as f:
        next(f)
        for row in f:
            row = row.replace("\n","").split("\t")
            gtdb_assembly_acc, name = row[:2]
            assembly_acc = gtdb_assembly_acc.split("_")[1] + "_" + gtdb_assembly_acc.split("_")[2]
            #if gtdb_assembly_acc in genomes_that_are_complete and assembly_acc not in already_done_accs:
            if gtdb_assembly_acc in genomes_that_are_complete:
                assembly_acc_name_dict[assembly_acc] = name
    f.close()

    return assembly_acc_name_dict


def get_assembly_dict(assembly_report_loc):
    """
    """
    assembly_ftp_dict = {}
    with open(assembly_report_loc, "r") as f:
        for row in f:
            if row[0] != "#":
                row = row.replace("\n","").split("\t")
                assembly_acc = row[0]
                asm_acc = row[17]
                ftp_path = row[19]

                # add both assembly_acc and ftp_path to dictionary lookup
                assembly_ftp_dict[assembly_acc] = ftp_path
                assembly_ftp_dict[asm_acc] = ftp_path
    f.close()

    return assembly_ftp_dict


def download_from_assembly_ftp(assembly_bio_data_loc, assembly_acc, ftp_address, ending):
    """
    """
    time.sleep(2) # pause 2 seconds - try not to piss off NCBI...
    #t_stop = time.time() + 60*2 # allow to run for 2 minutes then kill

    web_assembly_report_loc = ftp_address + "/" + ftp_address.split("/")[-1] + ending
    
    assembly_report_dl = assembly_bio_data_loc + assembly_acc + ending
    if ".gbff" in ending:
        assembly_report_dl = assembly_report_dl.replace(".gbff", ".gff")

    command = "timeout 120 wget " + web_assembly_report_loc + " -c -o /dev/null"
    dl_success = os.system(command)
    if dl_success != 0:
        return False

    os.system("mv " + ftp_address.split("/")[-1] + ending + " " + assembly_report_dl)
    
    if ".gz" in ending:
        os.system("gunzip " + assembly_report_dl)		   

    return True


def identify_download_assembly_reports(assembly_acc_name_dict, genbank_assembly_summary_loc, \
                                            refseq_assembly_summary_loc, bio_data_loc, \
                                                genome_assembly_ncbi_data_loc):
    """
    """
    # reset the biological data grab
    os.system("rm -rf " + bio_data_loc + "*")

    genbank_assembly_ftp_dict = get_assembly_dict(genbank_assembly_summary_loc)
    refseq_assembly_ftp_dict = get_assembly_dict(refseq_assembly_summary_loc)

    with open(genome_assembly_ncbi_data_loc, "w") as f:
        out = csv.writer(f)
        out.writerow(["Assembly Acc", "GTDB Taxonomy Name", \
                        "Assembly Type", "FTP address", "Download Success?"])

        for assembly_acc in assembly_acc_name_dict:
            taxonomy_name = assembly_acc_name_dict[assembly_acc]

            # check for refseq first, then genbank if not
            if assembly_acc in refseq_assembly_ftp_dict:
                ftp_address = refseq_assembly_ftp_dict[assembly_acc]
                assembly_type = "refseq"
            elif assembly_acc in genbank_assembly_ftp_dict:
                ftp_address = genbank_assembly_ftp_dict[assembly_acc]
                assembly_type = "genbank"
            else:
                ftp_address = ""
                assembly_type = "not_found"
            
            success = False
            if assembly_type == "refseq" or assembly_type == "genbank":
                assembly_bio_data_loc = bio_data_loc + assembly_acc + "/"
                os.system("mkdir " + assembly_bio_data_loc)

                # download assembly report
                assem_download_status = download_from_assembly_ftp(assembly_bio_data_loc, assembly_acc, ftp_address, "_assembly_report.txt")
                
                # download FASTA + GFF files
                fasta_download_status = download_from_assembly_ftp(assembly_bio_data_loc, assembly_acc, ftp_address, "_genomic.fna.gz")
                if assembly_type == "refseq":
                    gff_download_status = download_from_assembly_ftp(assembly_bio_data_loc, assembly_acc, ftp_address, "_genomic.gff.gz")
                elif assembly_type == "genbank":
                    gff_download_status = download_from_assembly_ftp(assembly_bio_data_loc, assembly_acc, ftp_address, "_genomic.gbff.gz")
        
                if fasta_download_status == True and gff_download_status == True:
                    success = True

            # tracking
            out.writerow([assembly_acc, taxonomy_name, assembly_type, ftp_address, success])

        f.close()

    return