#!/usr/bin/env python


"""
A collection of usefull functions for fasta/q file.

@author: Costas Bouyioukos
@organization: The Sainsbury Laboratory
@since: January 2011
@copyright: The program is coming as it is. You have the right to redistribute,
transform and change the source code presuming the apropriate reference and
the lisence is kept free.
@license: GNU GPL3 or newer.
@contact: U{Costas Bouyioukos<mailto:k.bouyioukos@uea.ac.uk>}
@version: 0.0.1
"""


from statlib import stats
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_fastx(fastaFile, typ = "fasta") :
  """Low level function to return a list of sequence record objects.

  This is the primitive function that all the rest of the functions in that
  module should call.
  """
  if typ not in ("fasta", "fastq") :
    raise StandardError, 'Only filetypes of fasta or fastq are supported.'
  seqList = []
  for seqRec in SeqIO.parse(fastaFile, typ) :
    seqList.append(seqRec)
  return seqList



def count_fasta(fastaFile) :
  """Return the number of sequences a fasta file contains.

  """
  return len(parse_fastx(fastaFile))



def length_distribution(fastaFile, sort = False) :
  """Return a list with the lengths of the sequences.

  """
  dist = []
  for seqRec in parse_fastx(fastaFile) :
    dist.append(len(seqRec))
  if sort :
    dist.sort(reverse = True)
  return dist


def print_length_distribution(fastaFile, sort = False) :
  """Return a list with the lengths of the sequences.

  """



def count_nucleotides(fastaFile) :
  """Return the number of bases the file contains.

  """
  nucs = 0
  for length in length_distribution(fastaFile) :
    nucs = nucs + length
  return nucs



def length_histogram(fastaFile, bins = 10) :
  """Return the length frequencies in 10 bins.

  """
  lengthsList = []
  for seqRec in parse_fastx(fastaFile) :
    lengthsList.append(len(seqRec))
  return stats.lhistogram(lengthsList, bins)



def print_histogram(fastaFile, bins = 10) :
  """Return a "pretty" print out of the histogram.

  """
  #Compute the bin centres
  hist = get_length_histogram(fastaFile, bins)
  binCentres = []
  for i in xrange(bins) :
    binCentres.append((i + 1) * hist[2] + hist[1])
  freqs = hist[0]
  print "Length\tFreq"
  for i in xrange(len(freqs)) :
    print "%i\t%i\n" % (int(binCentres[i]), freqs[i])



def calculate_N50(fastaFile) :
  """Return the N50 of a fasta file.

  """
  halfLenTotal = count_nucleotides(fastaFile) / 2.0
  currentLen = 0
  for length in length_distribution(fastaFile, sort = True) :
    currentLen = currentLen + length
    if currentLen >= halfLenTotal :
      n50 = length
      break
  return n50



if __name__ == "__main__" :

  #TODO Make the module a callable script for generic Fasta/q statistics.
  # Run the module as a script.
  import argparse
  import sys
  import shlex
  parser = argparse.ArgumentParser(description='Calculate some Fasta/q statistics.')
  parser.add_argument('-c', '--count',  )

