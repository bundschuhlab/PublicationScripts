import csv
import os
import numpy as np


def measure_genome_signals(genome_sd_msd_loc, sd_dg_threshold):
    """
    """
    all_signals = []
    all_sd_signals, all_msd_signals = [], []
    with open(genome_sd_msd_loc, "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            genome_id, locus_tag, sd_min_e, msd_min_e, asd, sd_seq, msd_seq = row

            signal = float(sd_min_e) - float(msd_min_e)
            all_signals.append(signal)
            all_sd_signals.append(float(sd_min_e))
            all_msd_signals.append(float(msd_min_e))
    f.close()

    average_diff = np.average(all_signals)
    
    num_sd_below_thresh, num_msd_below_thresh = 0, 0
    for signal in all_sd_signals:
        if signal <= sd_dg_threshold:
            num_sd_below_thresh += 1
    for signal in all_msd_signals:
        if signal <= sd_dg_threshold:
            num_msd_below_thresh += 1

    fract_below_thresh = (num_sd_below_thresh - num_msd_below_thresh) / len(all_signals)

    return average_diff, fract_below_thresh, num_sd_below_thresh / len(all_signals), num_msd_below_thresh / len(all_signals)


def find_lvls_from_species(species, tree):
    """
    """
    genus, family, order, classs, phylum, domain = "", "", "", "", "", ""

    parent_lineage = tree.ascend(species)
    for node in parent_lineage:
        taxid = node.taxid

        if "g__" in taxid:
            genus = taxid
        elif "f__" in taxid:
            family = taxid
        elif "o__" in taxid:
            order = taxid
        elif "c__" in taxid:
            classs = taxid
        elif "p__" in taxid:
            phylum = taxid
        elif "d__" in taxid:
            domain = taxid

    return genus, family, order, classs, phylum, domain


def examine_genome_wide_results(sd_dg_threshold, per_genome_sd_msd_out, \
                                    genome_wide_results, \
                                    genome_species_dict, tree):
    """
    """
    genome_signal_dict = {}
    for filename in os.listdir(per_genome_sd_msd_out):
        acc = filename.split("-")[0]

        average_diff, fract_below_thresh, sd_below, msd_below = measure_genome_signals(per_genome_sd_msd_out + filename, sd_dg_threshold)
        genome_signal_dict[acc] = [average_diff, fract_below_thresh, sd_below, msd_below]

    with open(genome_wide_results, "w") as f:
        out = csv.writer(f)
        out.writerow(["Acc", "Species", "Genus", "Family", \
                        "Order", "Class", "Phylum", "Domain", \
                        "Average SD - MSD Across Genome", \
                        "Fraction of SD - MSD Annotations Below Threshold", \
                        "Fraction of SDs Below Threshold", \
                        "Fraction of MSDs Below Threshold"])

        for acc in genome_signal_dict:
            average_diff, fract_below_thresh, sd_below, msd_below = genome_signal_dict[acc]
            species = genome_species_dict[acc].lower()

            genus, family, order, classs, phylum, domain = find_lvls_from_species(species, tree)

            out.writerow([acc, species, genus, family, \
                            order, classs, phylum, domain, \
                            average_diff, fract_below_thresh, \
                            sd_below, msd_below])
    f.close()


    return