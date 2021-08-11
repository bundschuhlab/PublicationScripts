# WoundEdgeHypermethylation 
This directory contains the custom scripts used to analyze RNA-Seq and DNA methylation data for the manuscript "Genome-Wide Hypermethylation in the Wound-Edge of Chronic Wound Patients Opposes Closure by Impairing Epithelial to Mesenchymal Transition" by Kanhaiya Singh, Yashika Rustagi, Ahmed Abouhashem, Saba Tabasum, Priyanka Verma, Edward Hernandez, Durba Pal, Dolly K Khona, Sujit K Mohanty, Manishekhar Kumar, Rajneesh Srivastava, Poornachander R Guda, Sanskruti Mahajan, Subhadip Ghatak, Shomita Steiner, Kristen E Wanczyk, Sheng Liu, Jun Wan, Pearlly Yan, Ralf Bundschuh, Savita Khanna, Gayle M Gordillo, Michael P Murphy, Sashwati Roy, and Chandan K Sen.

## System requirements

The code consists of two python scripts, which should run under any version of python 2.7. Specifically, they have been tested under python 2.7.5 on Linux version 3.10.0. The dmr_v_deg.py script depends on matplotlib while the sorting.py script does not have any dependencies outside of standard python libraries.

## Installation

No installation of the software is necessary. Assuming python 2.7 is installed on the computer the scripts can be run immediately.

## Sorting

The script sorting.py takes a list of differentially expressed genes as reported by DESeq2, filters for those with a baseMean of at least 10 and sorts the list by baseMean.

The demo is run via

    python sorting.py sorting_in.csv sorting_out.txt

It runs in fractions of a second and should create the output file sorting_out_expected.txt. Running it on any other data file is done by replacing sorting_in.dat by the name of the input file to be processed.

