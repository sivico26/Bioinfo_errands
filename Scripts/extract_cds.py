#!/usr/bin/env python
from Bio import SeqIO
import argparse
import pathlib


def get_arguments():
    parser = argparse.ArgumentParser(description='Extract CDS from a genbank to output a fasta',
			formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input', type=str,
			help='Path to input genbank file')
    parser.add_argument('output', type=str,help='Path to put file/folder output')

    parser.add_argument('-i', '--ignore', type=str, metavar = 'KEY', default=None, help="if 'key' matches a CDS name it won't be included in the output")
    parser.add_argument('-m', '--multi', action ='store_true', help = "Specify if the input file is a multigenbank, in which case the CDS of each entry would be extracted in a different fasta file in an output directory in the specified output path")

    args = parser.parse_args()

    return args

def get_features(record, key):
    cds = {}
    if key == None:
        for i,ft in enumerate(record.features):
            if ft.type == "CDS":
                if "gene" in ft.qualifiers.keys():
                    gene = ft.qualifiers["gene"][0]
                    cds[gene] = ft.extract(record)
    else:
        for i,ft in enumerate(record.features):
            if ft.type == "CDS":
                if "gene" in ft.qualifiers.keys():
                    if key not in ft.qualifiers["gene"][0]:
                        gene = ft.qualifiers["gene"][0]
                        cds[gene] = ft.extract(record)
    return cds

def reformat(cds):
    for gene, record in cds.items():
        record.id = gene
        record.description = ""
    return cds

def main():
    args = get_arguments()
    #if args.ignore == None:
    #    args.ignore == ""

    if args.multi is True:
        recs = SeqIO.parse(args.input,"gb")

        taxa = {}
        for rec in recs:
            specie = rec.annotations["organism"].replace(" ","_")
            taxa[specie] = reformat(get_features(rec, args.ignore))

        ## Create directory
        pathlib.Path(args.output.rstrip("/")+'/extract_cds_output').mkdir(parents=True, exist_ok=True)

        ## Write fastas
        for specie, genes in taxa.items():
            filepath = args.output.rstrip("/")+'/extract_cds_output'+"/"+specie+".fasta"
            SeqIO.write(genes.values(),filepath,"fasta")

    else:
        rec = SeqIO.read(args.input, "gb")
        aux = get_features(rec, args.ignore)
        cds = reformat(aux)

        ## Write filenames
        filename = args.output.strip("/")
#        filename = args.output.strip("/") + "/" + rec.annotations["organism"].replace(" ","_") + ".fasta"
        SeqIO.write(cds.values(), filename, "fasta")


if __name__ == '__main__':
    main()
