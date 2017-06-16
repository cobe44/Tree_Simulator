import pyvolve
import argparse
import sys


parser=argparse.ArgumentParser(description='Execute instances of pyvolve with various parameters.')

parser.add_argument('treefile', metavar='tree,', help="Path to treefile")
parser.add_argument('alpha', metavar='alpha,', type=float, help="Alpha parameter for gamma distribution")
parser.add_argument('size', metavar='partition_size', type=int, help="Size of partition to be created")

args=parser.parse_args()

treefile=str(args.treefile)

param_alpha = float(args.alpha)

partition_size = int(args.size)



#Read in phylogeny along which to simulate
mammal_tree = pyvolve.read_tree(file="amakihi.tree")

#Define nucleotide frequencies
freqs = [0.25, 0.25, 0.25, 0.25]

#Define evolutionary models
nuc_model = pyvolve.Model("nucleotide", {"kappa":1.75}, rate_factors=[0.01, 0.02, 0.03, 0.14, 0], alpha=param_alpha, num_categories=4, pinv=0.9)

#Define partitions; change this if you want to model heterogeneity within your tree
my_partition = pyvolve.Partition(models = nuc_model, size = partition_size)

# Evolve partitions with the callable Evolver class
my_evolver = pyvolve.Evolver(tree = mammal_tree, partitions = my_partition)

my_evolver(seqfile="alpha.fa")
