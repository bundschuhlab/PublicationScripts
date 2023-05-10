import gtdb_tree
taxdmp_loc =  "../leader_trailer/resources/temp/taxdmp/"

print("loading in taxtree")
tree = gtdb_tree.TaxonomyTree(taxdmp_loc)
print("taxtree loaded\n\n\n")	

tasks = [
    #"download_data", \
    #"get_stop_codon_freqs", \
    #"parse_arfa_output", \
    #"gen_scatterplots", \
    #"get_seqs_for_msa", \
    #"s21_arfa_study", \
    "sd_scanning"
]



clade_of_interest = "d__bacteria" 
rep_genome_only = True

refseq_assembly_summary_loc = "../leader_trailer/resources/temp/assembly_summary_refseq.txt"
sp_clusters = "../leader_trailer/resources/temp/sp_clusters.tsv"
bact_metadata = "../leader_trailer/resources/temp/bac120_metadata_r207.tsv"

assembly_dl_root = "resources/ncbi-data/"

assembly_data_summary = "results/bacteria-ncbi-download-summary.txt"
last_codon_freq_summary = "results/bacteria-last-codon-summary.csv"

arfa_output_dir = "bacteria-arfa-output/"
arfa_summary_loc = "results/bacteria-arfa-summary.csv"