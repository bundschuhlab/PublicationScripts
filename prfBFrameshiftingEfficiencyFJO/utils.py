import os
import csv


def check_mk_fldrs(fldrs_in):
    """
    """
    for fldr in fldrs_in:
        if os.path.isdir(fldr) == False:
            os.makedirs(fldr)

    return


def load_in_species_accs(assembly_data_summary):
    """
    """
    species_genome_dict = {}
    genome_species_dict = {}
    with open(assembly_data_summary, "r") as f:
        for row in f:
            row = row.replace("\n","").split("\t")
            species, accs = row

            accs = accs.split(";")
            species_genome_dict[species] = accs

            for acc in accs:
                genome_species_dict[acc] = species
    f.close()

    return species_genome_dict, genome_species_dict


def load_genome_asds(asd_loc):
    """
    """
    genome_asds = {}
    with open(asd_loc, "r") as f:
        reader = csv.reader(f)

        next(reader, None)
        for row in reader:
            assembly_acc, success, number, asd = row

            if asd == "":
                continue
            
            genome_asds[assembly_acc] = asd
    f.close()

    return genome_asds
