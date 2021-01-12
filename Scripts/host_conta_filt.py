#!/usr/bin/env python
from Bio import SeqIO
from itertools import compress
import argparse

def get_arguments():
    parser = argparse.ArgumentParser(description='Filter a genome based on a set of given clusters',
				     epilog="Developed by Simón Villanueva Corrales @sivico26",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input', type=str,
                        help='Path to genome input')
    parser.add_argument('-g', '--gff', type=str,
                        help='Path to gff')
    parser.add_argument('-o', '--output', type=str, default="host_genome.fasta" ,
                        help='Path to put output clean genome [%(default)s]')
    parser.add_argument('-c', '--contamination', type=str, default="contamination.fasta" ,
                        help='Path to put output contamination sequences [%(default)s]')
    parser.add_argument('-f', '--format', type=str, default="fasta", metavar="STR",
                        help='Input file format [%(default)s]')
    args = parser.parse_args()
    return args


def filt_clust(genoma, conta_list):
    index_conta = [True if req.name in conta_list else False for req in genoma]
    index_host = [not i for i in index_conta]
    genome_conta = compress(genoma, index_conta)
    genome_host = compress(genoma, index_host)
    return genome_host, genome_conta

def main():
    args = get_arguments()
    genome = list(SeqIO.parse(args.input, args.format))
    with open(args.gff,"r") as file:
        conta_list = [line.split()[0] for line in file.readlines()[1:]]
    genome_host, genome_conta = filt_clust(genome, conta_list)
    SeqIO.write(genome_host, args.output, args.format)
    SeqIO.write(genome_conta, args.contamination, args.format)
    
if __name__ == '__main__':
    main()
