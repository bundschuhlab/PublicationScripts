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

tasks = [
    "download_data", \
    "extract_asd", \
    "extract_s_l_seqs", \
    "run_freescan", \
    "gen_species_plots", \
    "gen_tree_plots", \
    "gen_s_l_prot_histos", \
    "get_flavo_genes"
]
force_redownload = False

barrnap_num_downstream = 10
barrnap_num_upstream_from_end = 100
illegal_letters = ["N","Y","R","S","W","K","M","B","V","D","H"]
min_number_s_l_annotations = 25
min_delg_for_s21_winner = -13.0

genome_assembly_ncbi_data_loc = "resources/assembly-download-report.csv"
asd_identification_loc = "resources/assembly-asd.csv"
s_l_annotation_loc = "resources/assembly-s_l-annotations.csv"

pipeline_output_loc = "results/pipeline-summary.csv"
prot_histogram_energy_loc = "results/gene_tir_histograms/"
tree_plot_loc = "results/trees/"
histogram_plot_loc = "results/s_l_prot_histograms/"
flavo_gene_freescan_loc = "results/all_genes_tirs/"

temp_data_root = "resources/temp/"
genbank_assembly_summary_loc = temp_data_root + "assembly_summary_genbank.txt"
refseq_assembly_summary_loc = temp_data_root + "assembly_summary_refseq.txt"
taxdmp_loc = temp_data_root + "taxdmp/"

bio_data_loc = "resources/bio_data/"

gtdb_clusters = temp_data_root + "sp_clusters.tsv"
gtdb_bact_tsv = temp_data_root + "bac120_metadata_r207.tsv"

cutadapt_exe = "$HOME/.local/bin/cutadapt"
barrnap_exe = "barrnap"
freescan_exe = "perl ../tools_resources/free2bind/free_scan.pl"