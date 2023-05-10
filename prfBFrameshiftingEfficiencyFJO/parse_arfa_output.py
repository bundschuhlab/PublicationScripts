import os
import csv


def mine_arfa_output_files(arfa_output_dir, arfa_summary_loc, assembly_data_summary):
    """
    """
    desired_genome_accs = set()
    with open(assembly_data_summary, "r") as f:
        for row in f:
            row = row.replace("\n","").split("\t")
            taxonomy, acc = row
            desired_genome_accs.add(acc)
    f.close()

    genome_arfa_result_dict = {}
    for filename in os.listdir(arfa_output_dir):
        genome_acc = filename.split("_arfa.txt")[0]

        if genome_acc not in desired_genome_accs:
            continue

        contains_prfb_annotation = False
        prfb_has_orf0 = False
        prfb_orf0_evalue = ""
        with open(arfa_output_dir + filename, "r") as f:
            grab_prfb_line = False
            for row in f:
                if row.startswith('                     /gene="prfB"'):
                    grab_prfb_line = True
                elif grab_prfb_line == True and "/note" in row:
                    note = row.split('"')[1]
                    grab_prfb_line = False

                    contains_prfb_annotation = True
                    if "ORF0" in note:
                        prfb_has_orf0 = True
                        e_value = note.split("ORF0: ")[1].split()[0].strip()
                        prfb_orf0_evalue = float(e_value)
        f.close()

        if contains_prfb_annotation == True and prfb_has_orf0 == True and prfb_orf0_evalue <= 0.1:
            prfb_status_call = True
        elif contains_prfb_annotation == True:
            prfb_status_call = False
        else:
            prfb_status_call = "NA"
        
        genome_arfa_result_dict[genome_acc] = [contains_prfb_annotation, prfb_has_orf0, prfb_orf0_evalue, prfb_status_call]

    with open(arfa_summary_loc, "w") as f:
        out = csv.writer(f)
        out.writerow(["Genome Acc", "Contains prfB annotation", \
                        "prfB has RF2 ORF0", "prfB RF2 ORF0 e-value", \
                        "Genome acc prfB call (True/False)"])

        for genome_acc in genome_arfa_result_dict:
            contains_prfb_annotation, prfb_has_orf0, prfb_orf0_evalue, prfb_status_call = genome_arfa_result_dict[genome_acc]

            out.writerow([genome_acc, contains_prfb_annotation, prfb_has_orf0, prfb_orf0_evalue, prfb_status_call])
    f.close()

    return