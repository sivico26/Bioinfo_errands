#!/usr/bin/env python
# reformat_headers.py
# Changes the headers of sequences files to format "scaffold<#>_<length>"
# where <#> is the number of the sequence in the file(by given order) and <length> is the length of the sequence

import sys, argparse
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser(description='Sequence headers reformatter.')
    parser.add_argument("-i", "--input", required=True, help='Input sequence file.')
    parser.add_argument("-o", "--output", help='Output file with kmer frequencies.')
    parser.add_argument("-f", "--format", default='fasta', dest="format", help='Input and output sequence format.')
    args = parser.parse_args()
    return args

def change_name(reqs):
    for i,req in enumerate(reqs,1):
        req.id = f"scaffold{i}_{len(req)}"
        req.description = ""
        yield req

def main():
    args = get_args()
    SeqIO.write(change_name(SeqIO.parse(args.input, args.format)), args.output, args.format)

if __name__ == "__main__":
    sys.exit(main())
