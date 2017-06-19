import sys
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna
import argparse

parser=argparse.ArgumentParser(description='Grab the first sequence from a multifasta to use as a reference for baitstools.')

parser.add_argument('fastafile', metavar='fasta infile', help="Path to fasta file")

args=parser.parse_args()

fastapath = str(args.fastafile)

seq_store = []
names = []
for record in SeqIO.parse(fastapath, "fasta"):
    names.append(record.id)
    seq_store.append(record.seq)

ref_seq = SeqRecord(seq_store[0], id = names[0], description="")

SeqIO.write(ref_seq, "reference.fasta", "fasta")
