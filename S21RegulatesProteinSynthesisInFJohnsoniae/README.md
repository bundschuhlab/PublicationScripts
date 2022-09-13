# S21RegulatesProteinSynthesisInFJohnsoniae
This directory contains the custom scripts used for the data analysis in the manuscript "Ribosomes lacking bS21 gain function to regulate protein synthesis in Flavobacterium johnsoniae" by Zakkary A. McNutt, Bappaditya Roy, Bryan T. Gemler, Elan A. Shatoff, Kyung-Mee Moon, Leonard J. Foster, Ralf Bundschuh, and Kurt Fredrick.

# Requirements
The following packages are required to run this pipeline:
  -Python3
  -Barrnap
  -Cutadapt
  -Freescan (note: requires old version of perl, tested using perl-5.8.9 which was spun up using perlbrew)
Paths can be set in config.py

# Quick user guide
All code is driven by main.py, which gathers configuration parameters, paths, and modules to run from config.py. The pipeline is set up to be modular and run in series, set modules to be run in the tasks list. Run the pipeline after setting desired configurations by opening a command line and running "python3 main.py".

The pipeline requires an external internet conection to download data from NCBI and GTDB's FTP sites. Note that one can download the data separately and import if desired, see code called by "download_data" task for web paths and download locations.

We observed that the pipeline will run in approximately a day to complete all default/standard steps.

Notes for future optimizations/nice to haves:
-We don't multithread well, if at all. Could optimize pipeline for 1 organism at a time, and kick off Python processes for each organism on all available threads
-The code was tested in a unix environment, where barrnap was installed using miniconda. It was observed that the Ubuntu barrnap installation did not predict 23S rRNA genes, so conda installation method was used instead. But the conda environment does not behave well with perlbrew, which is requried to roll perl version back to run freescan. Getting all tools to run on a single environment would be a nice to have.
-Currently, taxtree.py have hard-coded taxonomy file locations. Could set these as inputs to TaxonomyTree function that is called to generate tree utility.
