from config import *
from download_ncbi_ftp import load_ncbi_data
from gtdb_genome_handling import find_gtdb_species_genomes, identify_download_assembly_reports
from asd_identification_handling import extract_asd_seqs
from extract_s_l_seqs import extract_s_l_annotations
from freescan_operations import run_freescan_pipeline
from plot_species_tirs import gen_species_tir_plots
from plot_taxonomy_trees import gen_tree_plots
from gen_prot_histograms import get_s_l_prot_histograms
from run_all_genes import run_freescan_on_all_genes


def main():
    """
    """
    print("load data/taxtree as needed\n\n")
    tree = load_ncbi_data(temp_data_root, taxdmp_loc, genbank_assembly_summary_loc, \
                        refseq_assembly_summary_loc, force_redownload)

    if "download_data" in tasks or "run_freescan" in tasks:
        print("load GTDB genomes/accessions\n\n")
        assembly_acc_name_dict = find_gtdb_species_genomes(gtdb_clusters, gtdb_bact_tsv)

    if "download_data" in tasks:
        print("identify NCBI address for GTDB genomes/accessions\n\n")
        identify_download_assembly_reports(assembly_acc_name_dict, genbank_assembly_summary_loc, \
                                                refseq_assembly_summary_loc, bio_data_loc, \
                                                    genome_assembly_ncbi_data_loc)

    if "extract_asd" in tasks:
        print("run barrnap to identify 16S sequences, extract ASD from cutadapt motif\n\n")
        extract_asd_seqs(genome_assembly_ncbi_data_loc, asd_identification_loc, \
                            bio_data_loc, barrnap_exe, cutadapt_exe, \
                            barrnap_num_downstream, \
                            barrnap_num_upstream_from_end)

    if "extract_s_l_seqs" in tasks:
        print("extract S and L sequences from GFF file\n\n")
        extract_s_l_annotations(asd_identification_loc, bio_data_loc, \
                                    s_l_annotation_loc, illegal_letters)


    if "run_freescan" in tasks:
        print("freescan activities for S/L protein annotations\n\n")
        run_freescan_pipeline(s_l_annotation_loc, pipeline_output_loc, \
                                asd_identification_loc, \
                                min_number_s_l_annotations, \
                                bio_data_loc, freescan_exe, \
                                assembly_acc_name_dict)

    if "gen_species_plots" in tasks:
        print("generating per-species S/L protein TIR plots")
        gen_species_tir_plots(pipeline_output_loc, prot_histogram_energy_loc)

    if "gen_tree_plots" in tasks:
        print("generating itol tree")
        gen_tree_plots(tree_plot_loc, min_delg_for_s21_winner, tree, \
                        pipeline_output_loc)

    if "gen_s_l_prot_histos" in tasks:
        print("getting S and L protein MSD and SD TIRs")
        get_s_l_prot_histograms(pipeline_output_loc, bio_data_loc, histogram_plot_loc, \
                                    tree)

    if "get_flavo_genes" in tasks:
        print("getting flavo genes MSD and SD TIRs")
        run_freescan_on_all_genes(flavo_gene_freescan_loc, pipeline_output_loc, \
                                    illegal_letters, bio_data_loc, \
                                    asd_identification_loc, tree, \
                                    freescan_exe)


    return


main()