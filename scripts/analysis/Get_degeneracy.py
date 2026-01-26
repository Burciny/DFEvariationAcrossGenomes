import gzip
import pandas as pd
import re, sys
from Bio import SeqIO

input_ref = sys.argv[1]  #
output = sys.argv[2]  #

# codon degeneracy information
degeneracy_table = {
    "GCT": ["0fold", "0fold", "4fold"],
    "GCC": ["0fold", "0fold", "4fold"],
    "GCA": ["0fold", "0fold", "4fold"],
    "GCG": ["0fold", "0fold", "4fold"],
    "CGT": ["0fold", "0fold", "4fold"],
    "CGC": ["0fold", "0fold", "4fold"],
    "CGA": ["2fold", "0fold", "4fold"],
    "CGG": ["2fold", "0fold", "4fold"],
    "AGA": ["2fold", "0fold", "2fold"],
    "AGG": ["2fold", "0fold", "2fold"],
    "GGT": ["0fold", "0fold", "4fold"],
    "GGC": ["0fold", "0fold", "4fold"],
    "GGA": ["0fold", "0fold", "4fold"],
    "GGG": ["0fold", "0fold", "4fold"],
    "TTA": ["2fold", "0fold", "2fold"],
    "TTG": ["2fold", "0fold", "2fold"],
    "CTT": ["0fold", "0fold", "4fold"],
    "CTC": ["0fold", "0fold", "4fold"],
    "CTA": ["2fold", "0fold", "4fold"],
    "CTG": ["2fold", "0fold", "4fold"],
    "CCT": ["0fold", "0fold", "4fold"],
    "CCC": ["0fold", "0fold", "4fold"],
    "CCA": ["0fold", "0fold", "4fold"],
    "CCG": ["0fold", "0fold", "4fold"],
    "TCT": ["0fold", "0fold", "4fold"],
    "TCC": ["0fold", "0fold", "4fold"],
    "TCA": ["0fold", "0fold", "4fold"],
    "TCG": ["0fold", "0fold", "4fold"],
    "AGT": ["0fold", "0fold", "2fold"],
    "AGC": ["0fold", "0fold", "2fold"],
    "ACT": ["0fold", "0fold", "4fold"],
    "ACC": ["0fold", "0fold", "4fold"],
    "ACA": ["0fold", "0fold", "4fold"],
    "ACG": ["0fold", "0fold", "4fold"],
    "GTT": ["0fold", "0fold", "4fold"],
    "GTC": ["0fold", "0fold", "4fold"],
    "GTA": ["0fold", "0fold", "4fold"],
    "GTG": ["0fold", "0fold", "4fold"],
    "TTT": ["0fold", "0fold", "2fold"],
    "TTC": ["0fold", "0fold", "2fold"],
    "TAT": ["0fold", "0fold", "2fold"],
    "TAC": ["0fold", "0fold", "2fold"],
    "CAT": ["0fold", "0fold", "2fold"],
    "CAC": ["0fold", "0fold", "2fold"],
    "CAA": ["0fold", "0fold", "2fold"],
    "CAG": ["0fold", "0fold", "2fold"],
    "AAT": ["0fold", "0fold", "2fold"],
    "AAC": ["0fold", "0fold", "2fold"],
    "AAA": ["0fold", "0fold", "2fold"],
    "AAG": ["0fold", "0fold", "2fold"],
    "GAT": ["0fold", "0fold", "2fold"],
    "GAC": ["0fold", "0fold", "2fold"],
    "GAA": ["0fold", "0fold", "2fold"],
    "GAG": ["0fold", "0fold", "2fold"],
    "TGT": ["0fold", "0fold", "2fold"],
    "TGC": ["0fold", "0fold", "2fold"],
    "ATT": ["0fold", "0fold", "3fold"],
    "ATC": ["0fold", "0fold", "3fold"],
    "ATA": ["0fold", "0fold", "3fold"],
    "ATG": ["0fold", "0fold", "0fold"],
    "TGG": ["0fold", "0fold", "0fold"],
    "TAA": ["stop", "stop", "stop"],
    "TGA": ["stop", "stop", "stop"],
    "TAG": ["stop", "stop", "stop"]
}

codon_aa_dict = {
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
    "TTA": "Leu", "TTG": "Leu", "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser", "AGT": "Ser", "AGC": "Ser",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TTT": "Phe", "TTC": "Phe",
    "TAT": "Tyr", "TAC": "Tyr",
    "CAT": "His", "CAC": "His",
    "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn",
    "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp",
    "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", # Can be classified as twofold also?
    "ATG": "Met",
    "TGG": "Trp",
    "TAA": "Stop", "TGA": "Stop", "TAG": "Stop"
}


# This is specific to Drosophila reference fasta header: Change for others
def extract_gene_name(header):
    match = re.search(r'parent=([^;]+)', header)
    return match.group(1) if match else None

## Mouse CDS reference fasta header
# def extract_gene_name(header):
#     base= header.split("(")[0]
#     return base

## Yeast CDS reference fasta header
# def extract_gene_name(header):
#    match = re.search(r"gene:([^\s]+)", header)
#    return match.group(1) if match else None


def annotate_fasta(fasta_file, output_file):
    """Annotate each position in the reference CDS fasta file with degeneracy"""
    data = []
    if fasta_file.endswith(".gz"):
        handle = gzip.open(fasta_file, "rt")
    else:
        handle = open(fasta_file, "r")

    with handle:
        for record in SeqIO.parse(handle, "fasta"):
            gene_name = extract_gene_name(record.description)
            if not gene_name:
                continue

            seq = str(record.seq).upper()
            if len(seq) % 3 != 0:
                continue  # Skip sequences not divisible by 3

            position = 1
            for i in range(0, len(seq), 3):
                codon = seq[i:i + 3]
                if codon in degeneracy_table:
                    for j in range(3):
                        data.append([
                            gene_name,
                            position + j,
                            codon[j],
                            degeneracy_table[codon][j],
                            codon,
                            codon_aa_dict[codon],
                            j+1
                        ])
                position += 3

    # Convert to pandas DataFrame
    df = pd.DataFrame(data, columns=["Gene_name", "Position", "Ref_base", "Degeneracy", "codon", "aa", "codon_pos"])

    # Write to a tab-separated file
    df.to_csv(output_file, sep="\t", index=False)


annotate_fasta(input_ref, output)

