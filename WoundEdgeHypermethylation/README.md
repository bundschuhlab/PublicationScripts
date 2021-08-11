# WoundEdgeHypermethylation 
This directory contains the custom scripts used to analyze RNA-Seq and DNA methylation data for the manuscript "Genome-Wide Hypermethylation in the Wound-Edge of Chronic Wound Patients Opposes Closure by Impairing Epithelial to Mesenchymal Transition" by Kanhaiya Singh, Yashika Rustagi, Ahmed Abouhashem, Saba Tabasum, Priyanka Verma, Edward Hernandez, Durba Pal, Dolly K Khona, Sujit K Mohanty, Manishekhar Kumar, Rajneesh Srivastava, Poornachander R Guda, Sanskruti Mahajan, Subhadip Ghatak, Shomita Steiner, Kristen E Wanczyk, Sheng Liu, Jun Wan, Pearlly Yan, Ralf Bundschuh, Savita Khanna, Gayle M Gordillo, Michael P Murphy, Sashwati Roy, and Chandan K Sen.

## System requirements

The code consists of two python scripts, which should run under any version of python 2.7. Specifically, they have been tested under python 2.7.5 on Linux version 3.10.0. The `comparing.py` script depends on the commonly used packages `numpy` and `matplotlib` while the `sorting.py` script does not have any dependencies outside of standard python libraries.

## Installation

No installation of the software is necessary. Assuming python 2.7 with `numpy` and `matplotlib` is installed on the computer the scripts can be run immediately.

## Sorting

The script `sorting.py` takes a list of differentially expressed genes as reported by DESeq2, filters for those with a baseMean of at least 10 and sorts the list by baseMean.

The demo is run via

    python sorting.py degs.tsv sorting_out.txt

It runs in fractions of a second and creates the output file `sorting_out.txt`, which should agree with the file `sorting_out_expected.txt`. Running it on any other data file is done by replacing `sorting_in.tsv` by the name of the input file to be processed.

## Comparing

The script `comparing.py` takes a list of differentially expressed genes as reported by DESeq2 and a list of differentially methylated genes as reported from PrEMeR-CG and extracts genes that are both differentially expressed and methylated.

The demo is run via

    python comparing.py dmrs.tsv degs.tsv --make_plots

It runs in fractions of a second and should create output files `matched_dmrs_degs_same.txt`, `matched_dmrs_degs_opposite.txt`, and `dmrs_v_degs.png`. These should agree with the provided files  `matched_dmrs_degs_same_expected.txt`, `matched_dmrs_degs_opposite_expected.txt`, and `dmrs_v_degs_expected.png`, respectively (the image file might differ in binary format but should show the same image).  Running it on any other data file is done by replacing `dmrs.tsv` and `degs.tsv` by the names of the input files to be processed.
