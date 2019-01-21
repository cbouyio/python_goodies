#!/usr/bin/env python

"""A python tool to calculate local score of DNA/RNA sequences."""


SCORING_DNA = {"A" : -1.0,
               "C" : 2.0,
               "G" : 1.0,
               "T" : -1.0}
# TODO perhaps have the dictionary also as an external file for example for proteins.


def local_score(sequence, scoring = SCORING_DNA):
    """Main recursive function to calculate local score.

    sequence: A sequence of symbols (i.e string)
    scoring : A scoring scheme (i.e. a dictionary of symbols -> scores)
    return: A list of the local score
    """
    local_scores = [0.0]  # Initial local score is zero by definition.
    for symbol in sequence:
        local_scores.append(max(0, local_scores[-1] + scoring[symbol]))
    return(local_scores[1:])


if __name__ == '__main__':
    # If we launch it as a script load the biopython and argparse stuff.
    import argparse
    from Bio import SeqIO

    # Set up the command line arguments.
    parser = argparse.ArgumentParser(prog='localScore', description='Calculate and prints the local score of biosequences.', epilog="Author: Costas Bouyioukos, 2018")
    parser.add_argument('infile', nargs='?', default='-', type=argparse.FileType('r'), metavar="input_file", help='Path to a fasta DNA/RNA/protein file. (or STDIN)')
    parser.add_argument("outfile", nargs='?', default='-', type=argparse.FileType('w'), metavar='output_file', help="Path for the results file. (or STDOUT)")
    parser.add_argument('-s', '--scoring', type=str, help='A scoring system. Default: a hardcoded SCORING_DNA dictionary.')
    optArgs = parser.parse_args()

    # Iterate over all fasta sequences.
    for seq in SeqIO.parse(optArgs.infile, "fasta"):
        ls = local_score(seq.seq)
        optArgs.outfile.write("ID: {}\nSeq: {}\nScore: {}".format(seq.id, seq.seq, ' '.join(map(str, ls))))
    optArgs.outfile.write("\n")
