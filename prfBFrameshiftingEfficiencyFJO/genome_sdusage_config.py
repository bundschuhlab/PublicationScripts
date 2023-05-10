

threads = 10
sd_dg_threshold = -4

genome_wide_results = "results/genome_sd_msd_usage.csv"

assembly_data_summary = "resources/" + "taxonomy-genome-summary.txt"
asd_loc = "../s21_study/resources/assembly-asd.csv"
assembly_dl_root = "../leader_trailer/resources/assembly_dls/"

per_genome_sd_msd_out = "per_genome_sd_msd/"
taxdmp_loc = "../leader_trailer/resources/temp/" + "taxdmp/"
sp_clusters = "../leader_trailer/resources/temp/sp_clusters.tsv"
gtdb_full_tree_loc = "../codon_freqs/resources/bac120.tree"

illegal_letters = ["N", "Y", "R", "S", "W", "K", "M", "B", "V", "D", "H"]
freescan_exe = "perl /home/fung/Documents/tools_resources/free2bind/free_scan.pl"