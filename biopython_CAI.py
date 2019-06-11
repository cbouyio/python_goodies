#!/usr/env ptynon
"""A python module to calculate CAI (Codon Adaptation Index) from a fasta file.

Author: Costas Bouyioukos, @UMR7216, 2019
"""

import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import CodonUsage

# Taken from https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932
CodonsDictYeast = {'TTT': 170666, 'TTC': 120510, 'TTA': 170884, 'TTG': 177573,
                   'CTT': 80076, 'CTC': 35545, 'CTA': 87619, 'CTG': 68494,
                   'ATT': 196893, 'ATC': 112176, 'ATA': 116254, 'ATG': 136805,
                   'GTT': 144243, 'GTC': 76947, 'GTA': 76927, 'GTG': 70337,
                   'TAT': 122728, 'TAC': 96596, 'TAA': 6913, 'TAG': 3312,
                   'CAT': 89007, 'CAC': 50785, 'CAA': 178251, 'CAG': 79121,
                   'AAT': 233124, 'AAC': 162199, 'AAA': 273618, 'AAG': 201361,
                   'GAT': 245641, 'GAC': 132048, 'GAA': 297944, 'GAG': 125717,
                   'TCT': 153557, 'TCC': 92923, 'TCA': 122028, 'TCG': 55951,
                   'CCT': 88263, 'CCC': 44309, 'CCA': 119641, 'CCG': 34597,
                   'ACT': 132522, 'ACC': 83207, 'ACA': 116084, 'ACG': 52045,
                   'GCT': 138358, 'GCC': 82357, 'GCA': 105910, 'GCG': 40358,
                   'TGT': 52903, 'TGC': 31095, 'TGA': 4447, 'TGG': 67789,
                   'CGT': 41791, 'CGC': 16993, 'CGA': 19562, 'CGG': 11351,
                   'AGT': 92466, 'AGC': 63726, 'AGA': 139081, 'AGG': 60289,
                   'GGT': 156109, 'GGC': 63903, 'GGA': 71216, 'GGG': 39359}

CurrentCodonIndex = {}


def generate_index(codon_dict):
    """Generate a codon usage index from a FASTA file of CDS sequences.

    Take a codon usage dictionary and calculate the codon index in the CurrentCodonIndex module constant.

    RCSU values.
    """
    # TODO ABSOLUTELY INCORPORATE ALL THIS IN BIOPYTHON.
    #To calculate the index we first need to sum the number of times synonymous codons were used all together.
    for aa in SynonymousCodons:
        total = 0.0
        # RCSU values are CodonCount/((1/num of synonymous codons) * sum of all synonymous codons)
        codons = SynonymousCodons[aa]
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


def calculate_CAI(file, cdi):
    """Calculate the Codon Adaptation Index from a FASTA file.
    """
    caidf = pd.DataFrame(columns=["CAI"])
    SeqCai = CodonUsage.CodonAdaptationIndex()
    SeqCai.set_cai_index(cdi)
    for seq in SeqIO.parse(file, "fasta"):
        cai = SeqCai.cai_for_gene(str(seq.seq))
        idt = seq.id
        caidf.loc[idt] = float(cai)
    return caidf


if name = __main():
    parser = argparse.ArgumentParser(prog='biopython_CAI.py', description='Fetch transcript features (sequence and data) from ENSEMBL, for each gene ID in a list and output a multi FASTA file.', epilog="Author: Costas Bouyioukos, 2019, Paris UMR7216.")
    parser.add_argument('infile', nargs='?', default='-', type=argparse.FileType('r'), metavar="input_file", help='Path to the input FASTA file. (or STDIN).')
    parser.add_argument("outfile", nargs='?', default='-', type=argparse.FileType('w'), metavar='output_file', help="Path to output table file. (or STDOUT).")
    # TODO include an option with the codon usage table.

    # Parse the command line arguments.
    optArgs = parser.parse_args()

    # Generate the CA Index dictionary.
    generate_index(CodonsDictYeast)

    # Calculate the CAI pre gene.
    pdCai = calculate_CAI(optArgs.infile, CurrentCodonIndex)

    # Print it to outfile.
    pd.write.csv(pdCai, file = optArgs.outfile)
