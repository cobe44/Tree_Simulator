
import os
import sys
import math
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna
import argparse

parser=argparse.ArgumentParser(description='Chop up your FASTA based on provided parameters and sequences outlined in input BED file.')

parser.add_argument('bedfile', metavar='bed,', help="Path to bed file")
parser.add_argument('len', metavar='Fragment length,', type=int, help="Length of desired bait fragment, must be less than or equal to 120")
parser.add_argument('outfile', metavar='Outfile name', help="Name of output file")
parser.add_argument('fastafile', metavar='fasta infile', help="Path to fasta file")
parser.add_argument('exclude', nargs='?', type=int, default=0, help="Whether to exclude baits with more than a given number of SNPs (Optional)")


args=parser.parse_args()

bedfile=str(args.bedfile)

lenny = int(args.len)

outfile = str(args.outfile)

fastapath = str(args.fastafile)

exclude = int(args.exclude)




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

if(exclude==0):
    print("Normal mode achieved")
    sys.exit()
    for record in SeqIO.parse(fastapath, "fasta"):
        tmp_seq = ""
        names.append(record.id)
        for i in range(bed_data.shape[0]):
            if(lenny%2==0):
                tmp_seq += record.seq[int(bed_data[i][0]+(120-lenny)/2):int(bed_data[i][1]-(120-lenny)/2)+1]
            else:
                tmp_seq += record.seq[int(math.floor(bed_data[i][0]+(120-lenny)/2)):int(math.ceil(bed_data[i][1]-(120-lenny)/2))]
        seq_store.append(tmp_seq)

else:
    ref_list = []
    for record in SeqIO.parse(fastapath, "fasta"):
        ex_count = 0
        tmp_seq = ""
        names.append(record.id)
        if(record.id == "Common_Rosefinch"):
            for i in range(bed_data.shape[0]):
                if(lenny%2==0):
                    ref_list.append(record.seq[int(bed_data[i][0]+(120-lenny)/2):int(bed_data[i][1]-(120-lenny)/2)+1])
                else:
                    ref_list.append(record.seq[int(math.floor(bed_data[i][0]+(120-lenny)/2)):int(math.ceil(bed_data[i][1]-(120-lenny)/2))])
        for i in range(bed_data.shape[0]):
            if(lenny%2==0):
                cur_bait = record.seq[int(bed_data[i][0]+(120-lenny)/2):int(bed_data[i][1]-(120-lenny)/2)+1]

                count = sum(1 for a, b in zip(ref_list[i], cur_bait) if a != b)
                if(count<=exclude):
                    tmp_seq += cur_bait
                else:
                    ex_count += 1
            else:
                cur_bait = record.seq[int(math.floor(bed_data[i][0]+(120-lenny)/2)):int(math.ceil(bed_data[i][1]-(120-lenny)/2))]
                count = sum(1 for a, b in zip(ref_list[i], cur_bait) if a != b)
                if(count<=exclude):
                    tmp_seq += cur_bait
                else:
                    ex_count += 1
        seq_store.append(tmp_seq)
        print ex_count
print bed_data.shape[0]


proper_seqs = []


for i in range(len(names)):
    proper_seqs.append(SeqRecord(seq_store[i], id = names[i], description=""))




SeqIO.write(proper_seqs, outfile, "fasta")
