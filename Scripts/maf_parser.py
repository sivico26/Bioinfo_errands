#! /usr/bin/env python

import sys
import argparse

def get_arguments():
    parser = argparse.ArgumentParser(description='Parse MAF files to include species name in alingments',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', type=str,
                        help='Path to input MAF file')
    parser.add_argument('output', type=str,
                        help='Path to put new MAF file')
    parser.add_argument('-p', '--patterns', type=str,
                        help='Comma separated list of patterns (Align with species)')
    parser.add_argument('-s', '--species', type=str, 
                        help='Comma separated list of species (Align with patterns)')
    args = parser.parse_args()
    return args

def inc_sp(maf, pattern, sp):
    new_maf = []
    for line in maf:
        aux = line.split(" ")
        if line.startswith("s"):
            if aux[1].startswith(pattern):
                aux[1] = sp + aux[1]
        new_maf.append(aux)

    return new_maf

def main():
    args = get_arguments()
    with open(args.input,"r") as inp, open(args.output,"w") as out:
        new_maf = inc_sp(inp, args.patterns, args.sp)
        for line in new_maf:
            out.write(line)

if __name__ == '__main__':
    main()
