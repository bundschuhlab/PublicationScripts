

def fasta_to_dict(fasta_in):
    """
    """
    seqid_dict = {}
    with open(fasta_in, "r") as f:
        subj, seq = "", ""

        for row in f:
            if row.startswith(">"):
                if len(subj) != 0:
                    seqid_dict[subj] = seq
                subj = row[1:].replace("\n","").split()[0]
                seq = ""
            else:
                seq = seq + row.replace("\n","")
        
        if len(subj) != 0:
            seqid_dict[subj] = seq
    f.close()

    return seqid_dict

