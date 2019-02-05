#!/usr/bin/env python

"""A python tool to calculate local score of DNA/RNA sequences."""


def local_score(sequence, scoring):
    """Main recursive function to calculate local score.

    sequence: A sequence of symbols (i.e string)
    scoring : A scoring scheme (i.e. a dictionary of symbols -> scores)
    return: A list of the local score
    """
    local_scores = [0.0]  # Initial local score is zero by definition.
    for symbol in sequence:
        try:
            local_scores.append(max(0, local_scores[-1] + scoring[symbol]))
        except KeyError:
            raise KeyError('Symbol {} not found in scoring dictionary. Please check your scoring function and/or your input file.'.format(symbol))
    return(local_scores[1:])


if __name__ == '__main__':
    # If we launch it as a script load the biopython and argparse stuff.
    import argparse
    from Bio import SeqIO

    # Set up the command line arguments.
    parser = argparse.ArgumentParser(prog='localScore', description='Calculate and prints the local score of biosequences.', epilog="Author: Costas Bouyioukos, 2018")
    parser.add_argument('infile', nargs='?', default='-', type=argparse.FileType('r'), metavar="input_file", help='Path to a fasta DNA/RNA/protein file. (or STDIN)')
    parser.add_argument("outfile", nargs='?', default='-', type=argparse.FileType('w'), metavar='output_file', help="Path for the results file. (or STDOUT)")
    parser.add_argument('-s', '--scoring', type=str, help='A scoring system file (a whitespace separated file with nucl./aa. first column and score as seconf)). Default: a hardcoded SCORING_DNA dictionary for GC richness.', dest='scr')
    optArgs = parser.parse_args()

    # Choose the scoring functionself.
    SCORING = {}
    if optArgs.scr :
        with open(optArgs.scr) as f:
            for line in f:
                (key, val) = line.split()
                SCORING[key] = float(val)
    else:
        # By default the local score finds GC reach regions.
        SCORING = {"A" : -1.0,
                   "C" : 1.0,
                   "G" : 1.0,
                   "T" : -1.0}
                   
    # Iterate over all fasta sequences.
    for seq in SeqIO.parse(optArgs.infile, "fasta"):
        ls = local_score(seq.seq, SCORING)
        optArgs.outfile.write("ID: {}  LS_max: {}\nSeq: {}\nScore: {}".format(seq.id, max(ls), seq.seq, ' '.join(map(str, ls))))
    optArgs.outfile.write("\n")
