from prfB_FS_config import *
from download_data import genome_download
from parse_arfa_output import mine_arfa_output_files
from mine_out_seqs_for_msa import mine_seqs_for_msa
from find_scan_sds import scan_sd_seqs

def main():
    """
    """
    if "download_data" in tasks:
        genome_download(refseq_assembly_summary_loc, \
                        clade_of_interest, rep_genome_only, \
                        sp_clusters, bact_metadata, 
                        assembly_dl_root, assembly_data_summary, \
                        tree)

    if "parse_arfa_output" in tasks:
        mine_arfa_output_files(arfa_output_dir, arfa_summary_loc, assembly_data_summary)

    if "get_seqs_for_msa" in tasks:
        mine_seqs_for_msa(arfa_summary_loc, arfa_output_dir, tree)

    if "sd_scanning" in tasks:
        scan_sd_seqs(arfa_summary_loc, assembly_dl_root, arfa_output_dir, assembly_data_summary, tree, sp_clusters)

    return


main()
