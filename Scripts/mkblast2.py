#!/usr/bin/env python
from Bio import SeqIO
from Bio.Blast import NCBIWWW,NCBIXML
import pathlib
import argparse
#import asyncio

def get_arguments():
    parser = argparse.ArgumentParser(description='Make blasts for sequences from a file and summarise them',
				     epilog="Developed by Simón Villanueva Corrales @sivico26",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input', type=str,
                        help='Path to file input')
    parser.add_argument('output', type=str, default=".",
                        help='Path to put output')
    parser.add_argument('-f', '--format', type=str, default="fasta", metavar="STR",
                        help='Input file format [%(default)s]')
    parser.add_argument('-l', '--longest', type=int, default=None, metavar="INT",
                        help='Whether to take just the N longest sequences in file [%(default)s]')
    parser.add_argument('-u', '--utility', type=str, default="blastn", metavar="STR",
                        help='Blast utility to use [%(default)s]')
    parser.add_argument('-p', '--parse', action='store_true',
                        help='Do not make the blast searches, only generate summary [%(default)s]')
    parser.add_argument('-d', '--database', type=str, default="nt", metavar="STR",
                        help='Database to search in [%(default)s]')
    parser.add_argument('-m','--matches', type=int, default=5, metavar="INT",
                        help='Number of hits to report in summary per sequence [%(default)s]')
    args = parser.parse_args()
    return args

def read_file(file,form):
    return list(SeqIO.parse(file,form))

def filt_longest(reqs, ret = 10):
    lengths = [len(req.seq) for req in reqs]
    return [req for _,req in sorted(zip(lengths,reqs),reverse=True)][:ret]

def make_blast(req, output_folder, db ="nt", utility = "blastn"):
    output = output_folder.rstrip("/")
    pathlib.Path(output+'/blast_results').mkdir(parents=True, exist_ok=True)
    result_handle = NCBIWWW.qblast(utility, db, req.format("fasta"))
    with open(output+'/blast_results/'+req.name, "w") as out_handle:
        out_handle.write(result_handle.read())

def infoBLAST (blast_record, n=5):
    with open(blast_record,"r") as file:
        blast_record = NCBIXML.read(file)
    info_tot = []
    for alignment in blast_record.alignments[:n]:
        info = []
        for hsp in alignment.hsps:
            info.append(blast_record.query)
            info.append(alignment.hit_def) #name
            info.append(alignment.hit_id.split('|')[-2]) #accesion
            info.append(str(hsp.expect)) #evalue
            info.append(str((hsp.positives/hsp.align_length)*100)) #identity
            info.append(str(hsp.align_length/min(len(hsp.query),len(hsp.sbjct))))
        info_tot.append(info)
    return info_tot
## Orden de las sublistas:
## Query_#,Query_name, hit_#, hit_name, accsion, evalue,identity, coverage

def main():
    args = get_arguments()
    if args.parse:
        with open(args.output.strip("/")+"/blast_summary.tsv","w") as summary:
            summary.write("Query_#\thit_#\tQuery_name\thit_name\taccesion\tevalue\tidentity[%]\tcoverage[%]\n")
            for seq_n,xml in enumerate(pathlib.Path(args.input).iterdir()):
                for hit_n,hit in enumerate(infoBLAST(xml, n=args.matches)):
                    summary.write(f'{seq_n+1}\t{hit_n+1}\t'+"\t".join(hit)+"\n")
    else:
        files = read_file(args.input,args.format)
        if args.longest:
            files = filt_longest(files, ret=args.longest)
        for req in files:
            make_blast(req, args.output, db=args.database, utility=args.utility)
        with open(args.output.strip("/")+"/blast_summary.tsv","w") as summary:
            summary.write("Query_#\thit_#\tQuery_name\thit_name\taccesion\tevalue\tidentity[%]\tcoverage[%]\n")
            for seq_n,xml in enumerate(pathlib.Path(args.output+"/blast_results").iterdir()):
                for hit_n,hit in enumerate(infoBLAST(xml, n=args.matches)):
                    summary.write(f'{seq_n+1}\t{hit_n+1}\t'+"\t".join(hit)+"\n")

if __name__ == '__main__':
    main()
