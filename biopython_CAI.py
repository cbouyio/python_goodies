#!/usr/bin/env python3
"""A python tool module that exends with some functionalities the CodonUsage module of BioPython and provides interface to calculate CAI (Codon Adaptation Index) from a fasta file.

Runs as a command line application, or can be imported as a module.

Author: Costas Bouyioukos, @UMR7216, 2019
"""

# TODO develope it as a biopython extension to calculate and update the CAI indexices from direct connection to ENSEMBL.! An interesting addition.

import argparse
import json
import datetime
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import CodonUsage

CurrentCodonIndex = {}

def generate_codon_dict_fasta(fasta_file):
    """Calculate a codon counts dictionary from a fasta file.

    Return the dict.
    """
    codon_counts = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0,
                    'CTT': 0, 'CTC': 0, 'CTA': 0, 'CTG': 0,
                    'ATT': 0, 'ATC': 0, 'ATA': 0, 'ATG': 0,
                    'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0,
                    'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
                    'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0,
                    'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0,
                    'GAT': 0, 'GAC': 0, 'GAA': 0, 'GAG': 0,
                    'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0,
                    'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
                    'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0,
                    'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0,
                    'TGT': 0, 'TGC': 0, 'TGA': 0, 'TGG': 0,
                    'CGT': 0, 'CGC': 0, 'CGA': 0, 'CGG': 0,
                    'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
                    'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}
    with fasta_file as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            if len(rec.seq) < 60:
                continue;
            # make all UPPERcase.
            if str(rec.seq).islower():
                dna_seq = str(rec.seq).upper()
            else:
                dna_seq = str(rec.seq)
            for i in range(0, len(dna_seq), 3):
                codon = dna_seq[i:i+3]
                if codon in codon_counts:
                    codon_counts[codon] += 1
                elif len(codon) < 3:
                    pass
                elif "N" in codon:
                    pass
                else:
                    raise TypeError("Illegal codon {} in gene: {}".format(codon, rec.id))
    with open("codon_counts_" + datetime.datetime.now().strftime('%Y%m%d_%H%M%S') + ".json", "w") as js:
        json.dump(codon_counts, js)
    return codon_counts


def generate_index_dict(codon_dict):
    """Generate a codon usage index from a codon usage dictionary.

    Calculate the codon index as the CurrentCodonIndex module constant.

    RCSU values.
    """
    # TODO ABSOLUTELY INCORPORATE ALL THIS IN BIOPYTHON.
    #To calculate the index we first need to sum the number of times synonymous codons were used all together.
    for aa in CodonUsage.SynonymousCodons:
        total = 0.0
        # RCSU values are CodonCount/((1/num of synonymous codons) * sum of all synonymous codons)
        codons = CodonUsage.SynonymousCodons[aa]
        for codon in codons:
            total += codon_dict[codon]
        # calculate the RSCU value for each of the codons
        rcsu = []
        for codon in codons:
            denominator = float(total) / len(codons)
            rcsu.append(codon_dict[codon] / denominator)
        # now generate the index W=RCSUi/RCSUmax:
        rcsu_max = max(rcsu)
        for codon_index, codon in enumerate(codons):
            CurrentCodonIndex[codon] = rcsu[codon_index] / rcsu_max


def calculate_CAI(file, cci=CurrentCodonIndex):
    """Calculate the Codon Adaptation Index from a FASTA file.
    """
    caidf = pd.DataFrame(columns=["CAI"])
    SeqCai = CodonUsage.CodonAdaptationIndex()
    SeqCai.set_cai_index(cci)
    for seq in SeqIO.parse(file, "fasta"):
        cai = SeqCai.cai_for_gene(str(seq.seq))
        idt = seq.id
        caidf.loc[idt] = float(cai)
    return caidf


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='biopython_CAI.py', description='Fetch transcript features (sequence and data) from ENSEMBL, for each gene ID in a list and output a multi FASTA file.', epilog="Author: Costas Bouyioukos, 2019, Paris UMR7216.")
    parser.add_argument('infile', nargs='?', default='-', type=argparse.FileType('r'), metavar="input_file", help='Path to the input FASTA file. (or STDIN).')
    parser.add_argument("outfile", nargs='?', default='-', type=argparse.FileType('w'), metavar='output_file', help="Path to output table file. (or STDOUT).")
    codon_counts = parser.add_mutually_exclusive_group(required=True)
    codon_counts.add_argument('-o', '--orfs', nargs = "?", help="Path to a FASTA file contianing all ORFs that we need to calculate CAI for a new organism. (Default=None).", type=argparse.FileType('r'), default=None, dest="orfs", metavar="ORFsFile")
    codon_counts.add_argument('-j', '--jsonFile', nargs = "?", help="Path to a JSON file containing the dictionary of codon counts for a known organism. (Default=None).", type=argparse.FileType('r'), default=None, dest="jscd", metavar="CountsJSON")

    # Parse the command line arguments.
    optArgs = parser.parse_args()

    # Load or generate the codon counts.
    if optArgs.orfs:
        codon_counts = generate_codon_dict_fasta(optArgs.orfs)
    elif optArgs.jscd:
        codon_counts = json.load(optArgs.jscd)

    # Generate the CAIndex dictionary.
    generate_index_dict(codon_counts)

    print(CurrentCodonIndex)

    import sys
    sys.exit()
    # Calculate the CAI pre gene.
    pdCai = calculate_CAI(optArgs.infile)

    # Print it to outfile.
    pdCai.to_csv(optArgs.outfile, sep = "\t")
