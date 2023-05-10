from genome_sdusage_utils import check_mk_fldrs, load_in_species_accs, \
                    load_genome_asds
from run_freescan import run_all_genomes_sd
import gtdb_tree
from examine_genome_sd_results import examine_genome_wide_results
from plot_taxonomies import generate_tree_of_results
from config import *


def main():
    """
    perlbrew use perl-5.8.9
    """
    # load gtdb taxonomy tree
    print("loading in taxtree")
    tree = gtdb_tree.TaxonomyTree(taxdmp_loc)
    print("taxtree loaded\n\n\n")	

    # get list of species/accs that have successful data download
    species_genome_dict, genome_species_dict = load_in_species_accs(assembly_data_summary)

    # get ASD of each genome
    genome_asds = load_genome_asds(asd_loc)

    # run free_scan on each genome's annotations' SD and MSD regions
    run_all_genomes_sd(per_genome_sd_msd_out, illegal_letters, \
                            genome_asds, assembly_dl_root, \
                            freescan_exe, threads)

    # concat free_scan results using both methods (averages vs threshold)
    examine_genome_wide_results(sd_dg_threshold, per_genome_sd_msd_out, \
                                    genome_wide_results, \
                                    genome_species_dict, tree)

    # generate tree figures
    generate_tree_of_results(species_genome_dict, genome_species_dict, \
                                genome_wide_results, \
                                tree, sp_clusters, gtdb_full_tree_loc)

    return


main()
