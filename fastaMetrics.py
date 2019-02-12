#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A collection of useful functions to report statistics of for fasta files.

Works as a script if it is called by its name.

@author: Costas Bouyioukos
@organization: The Sainsbury Laboratory
@since: January 2011
@copyright: The program is coming as it is. You have the right to redistribute,
transform and change the source code presuming the appropriate reference and
the license is kept free.
@license: GNU GPL3 or newer.
@contact: U{Costas Bouyioukos<mailto:cbouyio@gmail.com>}
@version: 0.2.1"""


import sys
import re

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_fasta(fastaFile, typ = "fasta") :
  """Low level function to return an iterator of sequence record objects.

  This is the primitive function that all the rest of the functions in that
  module should call."""
  if typ not in ("fasta", "fastq") :
    raise StandardError, 'Only filetypes of fasta or fastq are supported.'
  seqList = []
  for seqRec in SeqIO.parse(fastaFile, typ) :
    seqList.append(seqRec)
  return seqList

def count_fasta(fastaFile) :
  """Return the number of sequences a fasta file contains."""
  return len(parse_fasta(fastaFile))

def count_lengths(fastaFile, sort = False) :
  """Return a list with the lengths of the sequences."""
  dist = []
  for seqRec in parse_fasta(fastaFile) :
    dist.append(len(seqRec))
  if sort :
    dist.sort(reverse = True)
  return dist

def print_lengths(fastaFile, sort = False) :
  """Print a list with the lengths of the sequences. Elsewhere, pandas.cut() is a convenient way to bin values into arbitrary intervals. Let’s say you have some data on ages of individuals and want to bucket them sensibly:"""
  # TODO I do not know if we need this function.
  print 'length'
  for seqLen in count_lengths(fastaFile, sort) :
    print str(seqLen)

def count_nucleotides(fastaFile) :
  """Return the number of bases the file contains."""
  return sum(count_lengths(fastaFile))

def length_histogram(fastaFile, bins) :
  """Return a tuple of frequencies and bin edges."""
  lengthsList = []
  for seqRec in parse_fasta(fastaFile) :
    lengthsList.append(len(seqRec))
  return np.histogram(lengthsList, bins)

def print_histogram(fastaFile, bin = 'auto') :
  """Return a "pretty" print out of the histogram."""
  #Compute the bin centres
  hist, bins = length_histogram(fastaFile, bin)
  # binRanges = []
  # for i in xrange(bins) :
  #   binRange = [int(round(hist[1] + i * hist[2], 0)), int(round(hist[1] + (i + 1) * hist[2], 0))]
  #   if binRange[0] < 0.0 :
  #     binRange[0] = 0
  #   binRanges.append(binRange)
  # freqs = hist[0]
  maxDigits1 = 2 * len(comma_me(str(bins[-1]))) + 1
  maxDigits2 = len(comma_me(str(hist[-1]))) + 1
  print 'Length'.center(maxDigits1), ':', 'Counts (freq)'.center(maxDigits2), ':', 'Hist'.center(maxDigits1 + maxDigits2)
  for i in xrange(len(hist)) :
    binStr = "{}-{}".format(comma_me(str(int(bins[i]))), comma_me(str(int(bins[i+1]) - 1)))
    freqStr = comma_me(str(hist[i])) + " ({0:.3f})".format(hist[i]/float(count_fasta(fastaFile)))
    histStr = '{}'.format('+' * hist[i])
    print binStr.rjust(maxDigits1), ':', freqStr.rjust(maxDigits2), ':', histStr.ljust(maxDigits1 + maxDigits2)

def min_max_lengths(fastaFile) :
  """Return the min and max of the sequence lengths."""
  lengths = count_lengths(fastaFile, sort = True)
  return (lengths[-1], lengths[0])

def calculate_N50(fastaFile) :
  """Return the N50 of a fasta file."""
  #TODO Include a reference file as an option for the calculation of the N50.
  halfLenTotal = count_nucleotides(fastaFile) / 2.0
  currentLen = 0
  n50 = 0
  for length in count_lengths(fastaFile, sort = True) :
    currentLen = currentLen + length
    if currentLen >= halfLenTotal :
      n50 = length
      break
  return n50

#TODO Replace the use of stalib by numpy!!!
# def mean_sd_median_lengths(fastaFile) :
#   """Return some descriptive statistics of the sequence lengths.
#
#   Mean, SD, Median, Coefficient of variation.
#   """
#   lengths = count_lengths(fastaFile)
#   mean = stats.lmean(lengths)
#   sd = stats.stdev(lengths)
#   median = stats.medianscore(lengths)
#   coefVar = stats.variation(lengths)
#   return (mean, sd, median, coefVar)

# def mode(fastaFile) :
#   """Return the mode (the most frequent value of the data.)"""
#   lengths = count_lengths(fastaFile)
#   md = stats.lmode(lengths)[1][-1]
#   return md

def comma_me(amount):
  """Printout recipe.
  Taken from: http://code.activestate.com/recipes/146461-commafying-an-integer/"""
  orig = amount
  new = re.sub("^(-?\d+)(\d{3})", '\g<1>,\g<2>', amount)
  if orig == new:
    return new
  else:
    return comma_me(new)


if __name__ == "__main__" :
  import argparse

  parser = argparse.ArgumentParser(description = 'Python script (and module) to calculate a series of statistics related to sequences contained in a Fasta/q file.')
  parser.add_argument('fasta', nargs = '?', default = '-', type = argparse.FileType('r'), metavar = 'path_fasta/q', help = 'A fasta/q input file path or STDIN. (Default: STDIN)')
  parser.add_argument('out', nargs = '?', default = '-', type = argparse.FileType('w'), metavar = 'path_output', help = 'The output file path or STDOUT. (Default: STDOUT)')
  parser.add_argument('-b', '--bins', type = int, metavar = 'no_of_bins', help = "The number of bins for the calculation of the length distribution. (Default: 'auto')", default = 10, dest = 'bins')
  parser.add_argument('-t', '--type', type = str, metavar = 'type_of_file', help = "Designate the type of the Fasta/q file. (Default: fasta)", default = 'fasta', dest = 'tp')

  ns = parser.parse_args()
  bins = ns.bins
  tp = ns.tp
  fastaFile = ns.fasta

  print 'File "' + fastaFile.name + '" has the following statistics:'
  print 'Number of sequences          : %s' % comma_me(str(count_fasta(fastaFile)))
  print 'Number of nucleotides        : %s' % comma_me(str(count_nucleotides(fastaFile)))
  mn, mx = min_max_lengths(fastaFile)
  print 'Longest Sequence             : %s' % comma_me(str(mx))
  print 'Shortest Sequence            : %s' % comma_me(str(mn))
  mean, sd, median, cv = mean_sd_median_lengths(fastaFile)
  print 'Mean sequence length         : %.2f' % mean
  print 'STDV sequence length         : %.2f' % sd
  print 'Variation of sequence length : %.2f' % cv
  print 'Median sequence length       : %s' % comma_me(str(median))
  print 'Mode of sequence length      : %s' % comma_me(str(mode(fastaFile)))
  print 'N50 of the current file      : %s' % comma_me(str(calculate_N50(fastaFile)))
  print '-Sequence length histogram:'
  print_histogram(fastaFile, bins)
