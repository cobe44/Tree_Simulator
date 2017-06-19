
import os
import sys
import math
import numpy as np
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna
import argparse

parser=argparse.ArgumentParser(description='Chop up your FASTA based on provided parameters and sequences outlined in input BED file.')

parser.add_argument('bedfile', metavar='bed,', help="Path to bed file")
parser.add_argument('len', metavar='Fragment length,', type=int, help="Length of desired bait fragment, must be less than or equal to 120")
parser.add_argument('outfile', metavar='Outfile name', help="Name of output file")

args=parser.parse_args()

bedfile=str(args.bedfile)

lenny = int(args.len)

outfile = str(args.outfile)




with open(bedfile, 'r') as infile:
    bed_data = (line.split('\t') for line in infile)
    bed_data = np.asarray(list(bed_data))
    infile.close()

for i in range(bed_data.shape[0]):
    bed_data[i][2] = bed_data[i][2].rstrip('\n')

bed_data = np.delete(bed_data, 0, 1).astype(int)

bed_data = np.subtract(bed_data, 1)

seq_store = []
names = []
for record in SeqIO.parse("final.fasta", "fasta"):
    tmp_seq = ""
    names.append(record.id)
    for i in range(bed_data.shape[0]):
        tmp_seq += record.seq[bed_data[i][0]+(120-lenny)/2:bed_data[i][1]-(120-lenny)/2]
    seq_store.append(tmp_seq)
    print(len(tmp_seq))

proper_seqs = []

for i in range(len(names)):
    proper_seqs.append(SeqRecord(seq_store[i], id = names[i], description=""))




SeqIO.write(proper_seqs, outfile, "fasta")
