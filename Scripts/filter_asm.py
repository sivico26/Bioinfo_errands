#!/usr/bin/env python
from Bio import SeqIO
from itertools import accumulate
import argparse
import sys

def get_arguments():
    parser = argparse.ArgumentParser(description='Filter shortest sequences of an assembly',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=str,
                        help='Paths to assembly')
    parser.add_argument('-o', '--output', type=str, default="./asm.fasta",
                        help='Path to filtered assembly (%(default)s)')
    
    criteria = parser.add_mutually_exclusive_group(required = True)
    criteria.add_argument('-g', '--genome_size', type=int, default=None,
                        help='Filter the shortest sequences after genome size is reached')
    criteria.add_argument('-s', '--size', type=int, default=None,
                       help='Filter sequences smaller than size')
    criteria.add_argument('-n', '--number', type=int, default=None,
                       help='Keep the N longest sequences')
    
    args = parser.parse_args()
    return args

def filter_asm(asm, output, genome_size = None, seq_size = None, seq_num = None):
    reqs = sorted(list(SeqIO.parse(asm,"fasta")), key=lambda x: len(x), reverse = True)
    if genome_size:
        higher = lambda x: x > genome_size
        idx = tuple(map(higher, accumulate(map(len, reqs)))).index(True) + 1 
        SeqIO.write(reqs[:idx], output, "fasta")
    elif seq_size:
        SeqIO.write(filter(lambda x: len(x) > seq_size, reqs), output, "fasta")
    elif seq_num:
        SeqIO.write(reqs[:seq_num], output, "fasta")
    else:
        print("No filtering?")
        sys.exit()
        
def main():
    args = get_arguments()
    filter_asm(args.input, args.output, args.genome_size, args.size, args.number)
    #NGx(args.output, *args.input, genome_size = args.genome_size, ylog= args.ylog)

if __name__ == "__main__":
    main()
