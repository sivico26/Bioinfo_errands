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
    parser.add_argument('-d', '--dat', type=str,
                        help='Path to clusters file(.dat)')
    parser.add_argument('-c', '--clusters', type=int, nargs='+',
                        help='Spaced list of clusters to be included')
    parser.add_argument('-o', '--output', type=str, default="genome_filt.fasta" ,
                        help='Path to put output genome [%(default)s]')
    parser.add_argument('-f', '--format', type=str, default="fasta", metavar="STR",
                        help='Input file format [%(default)s]')
    args = parser.parse_args()
    return args


def filt_clust(genoma, clust_list, clusters):
    index_list = [True if int(line.split()[0]) in clusters else False for line in clust_list]
    return compress(genoma, index_list)

def main():
    args = get_arguments()
    genome = list(SeqIO.parse(args.input, args.format))
    with open(args.dat,"r") as file:
        clust_list = file.readlines()
    genome_filt = filt_clust(genome, clust_list,args.clusters)
    SeqIO.write(genome_filt, args.output, args.format)

if __name__ == '__main__':
    main()
