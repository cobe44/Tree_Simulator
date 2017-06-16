#!/bin/bash

DIRNAME=$1

SIZE=$2

ALPHA=$3

python Gammas.py amakihi.tree $ALPHA $SIZE

mkdir $DIRNAME

mv alpha.fa site_rates* $DIRNAME

python3 /home/genetics/jacob/PoMo/FastaToVCF.py $DIRNAME/alpha.fa variants001.vcf

mv variants001.vcf $DIRNAME

./iqtree-omp -s $DIRNAME/alpha.fa -m HKY+G4+I -bb 1000 -nt AUTO
