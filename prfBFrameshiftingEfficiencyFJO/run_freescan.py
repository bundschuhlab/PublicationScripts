import os
import csv
from Bio.Seq import Seq
from multiprocessing import Process, Queue, Pool, Manager


def load_fasta_to_dict(fasta):
    """
    """
    seqid_dict = {}
    with open(fasta, "r") as f:
        subj, seq = "", ""
        for row in f:
            if row.startswith(">"):
                if len(subj) != "":
                    seqid_dict[subj] = seq
                subj = row[1:].split()[0].replace("\n","")
                seq = ""
            else:
                seq = seq + row.replace("\n","").upper()

        if len(subj) != 0:
            seqid_dict[subj] = seq
    f.close()

    return seqid_dict


def run_freescan(freescan_exe, fasta_loc, asd, acc, region_type):
    """
    """
    freescan_out = "temp/freescan_out_" + acc + "-" + region_type + ".txt"

    command = freescan_exe + " -q " + asd[::-1] + \
                            " " + fasta_loc + \
                            " > " + freescan_out
    os.system(command)

    min_e = 0.0
    with open(freescan_out, "r") as f:
        for row in f:
            if row[0] != "#" and row.strip() != "":
                signal = float(row.split()[0])

                if signal < min_e:
                    min_e = signal
    f.close()

    return min_e


def load_annotations(annotation_loc):
    """
    """
    genome_annotations = []
    with open(annotation_loc, "r") as f:
        for row in f:
            if row.startswith("#"):
                continue

            row = row.replace("\n","").split("\t")
            genome_id = row[0]

            annot_type = row[2]
            if annot_type != "CDS":
                continue

            pos_start = int(row[3])
            pos_end = int(row[4])

            annot_strand = row[6]
            if annot_strand == "+" or annot_strand == "-":
                strand = row[6]
            elif annot_strand == "1":
                strand = "+"
            elif annot_strand == "-1":
                strand = "-"
            else:
                print("NOTE: Unmapped strand here - ", annot_strand)

            id_string = row[8]

            genome_annotations.append([genome_id, id_string, pos_start, pos_end, strand])
    f.close()

    return genome_annotations


def run_genome_get_sd_msd(acc, genome_loc, annot_loc, \
                            sd_msd_loc, illegal_letters, \
                            freescan_exe, asd):
    """
    """
    # load full genome seq
    genome_seq_dict = load_fasta_to_dict(genome_loc)    

    # load annotations
    genome_annotations = load_annotations(annot_loc)

    with open(sd_msd_loc, "w") as f:
        out = csv.writer(f)
        
        out.writerow(["Genome ID", "Locus Tag", \
                        "SD Min E", "MSD Min E", \
                        "ASD", "SD", "MSD"])
        
        for genome_id, locus_tag, pos_start, pos_end, strand in genome_annotations:
            genome_seq = genome_seq_dict[genome_id]

            # get sd/msd regions
            sd_region, msd_region = "temp/" + acc + "sd.fa", "temp/" + acc + "msd.fa"

            if strand == "+":
                seq_por = genome_seq[pos_start - 55 : pos_start - 1].upper()

                if pos_start < 56 or any((nt in illegal_letters) for nt in seq_por):
                    continue

            elif strand == "-":
                seq_por = genome_seq[pos_end : pos_end + 54].upper()
                if (pos_end > len(genome_seq) - 56) or any((nt in illegal_letters) for nt in seq_por):
                    continue

                dna = Seq(seq_por)
                seq_por = str(dna.reverse_complement())


            sd_seq = seq_por[25:]
            msd_seq = seq_por[:29]
            with open(sd_region, "w") as f:
                f.write(">" + "SD" + "\n")
                f.write(sd_seq + "\n")
            f.close()

            with open(msd_region, "w") as f:
                f.write(">" + "MSD" + "\n")
                f.write(msd_seq + "\n")
            f.close()

            # run freescan
            sd_min_e = run_freescan(freescan_exe, sd_region, asd, acc, "SD")
            msd_min_e = run_freescan(freescan_exe, msd_region, asd, acc, "MSD")

            out.writerow([genome_id, locus_tag, sd_min_e, msd_min_e, asd, sd_seq, msd_seq])
    f.close()

    # reset memory (why does Python make me do this...)
    genome_seq_dict = {}
    genome_annotations = []

    return


def run_all_genomes_sd(per_genome_sd_msd_out, illegal_letters, \
                        genome_asds, assembly_dl_root, \
                        freescan_exe, threads):
    """
    """
    m = Manager()
    q = m.Queue()
    p = Pool(threads)

    already_done_accs = set()
    for filename in os.listdir(per_genome_sd_msd_out):
        acc = filename.split("-")[0]
        already_done_accs.add(acc)

    tasks = []
    for acc in genome_asds:
        if acc in already_done_accs:
            continue
        print(acc)
        asd = genome_asds[acc]

        genome_loc = assembly_dl_root + acc + "_genomic.fna"
        annot_loc = assembly_dl_root + acc + "_genomic.gff"

        if os.path.isfile(genome_loc) == False or os.path.isfile(annot_loc) == False:
            continue

        sd_msd_loc = per_genome_sd_msd_out + acc + "-annot_sd_msd.csv"

        tasks.append(p.apply_async(
            run_genome_get_sd_msd, (acc, genome_loc, annot_loc, \
                                sd_msd_loc, illegal_letters, \
                                freescan_exe, asd, )))

    print("# of genomes to run:", len(tasks))
    print("# of genomes already done:", len(already_done_accs))
    [r.get() for r in tasks]

    return