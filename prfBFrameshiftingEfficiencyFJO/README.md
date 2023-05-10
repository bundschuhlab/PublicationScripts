# prfBFrameshiftingEfficiencyFJO
This directory contains the custom scripts used for the data analysis in the manuscript "Occlusion of the anti-Shine-Dalgarno on the 30S platform reduces the efficiency of prfB programmed frameshifting in Flavobacterium johnsoniae" by Fawwaz Naeem, Bryan T. Gemler, Zakkary A. McNutt, Ralf Bundschuh, and Kurt Fredrick.

# Requirements
The following packages are required to run this pipeline:
  -Python3
  -Barrnap
  -Cutadapt
  -Freescan (note: requires old version of perl, tested using perl-5.8.9 which was spun up using perlbrew)
Paths can be set in respective ...config.py, i.e., "genome_sdusage_config.py" and "prfB_FS_config.py"

# Quick user guide
All code is driven by each ...main.py, which gathers  configuration parameters, paths, and modules to run from its respective ...config.py. The pipeline is set up to be modular and run in series, set modules to be run in the tasks list. Run the pipeline after setting desired configurations by opening a command line and running "python3 ...main.py".

* prfB_FS_main.py: this script parses over ARFA output to identify genomes with a prfB FS site and extracts the surrounding sequence region. Freescan is run using the detected ASD to determine free energy dG (kcal/mol). Phylogenetic trees of SD strength for iTOL input are generated. Nucleotide and protein MSA files are processed to generate taxonomy-specific word logos.
* genome_sdusage_main.py: this script parses over the entire bacteria GTDB tree and runs freescan against the SD and MSD region of every coding sequence. Genome-wide signals of SD - MSD are compiled and phylogenetic trees for iTOL input are generated.
